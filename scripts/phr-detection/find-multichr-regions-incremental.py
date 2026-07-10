#!/usr/bin/env python3
"""
Incremental inter-chromosomal region finder with early stopping.

Integrates impg query execution with early stopping to avoid unnecessary queries.
Processes sequences in parallel but windows sequentially within each sequence.
"""
import argparse
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Thread-safe print lock
print_lock = threading.Lock()

def parse_args():
    parser = argparse.ArgumentParser(
        description='Find multi-chromosomal signal regions with incremental impg queries'
    )
    parser.add_argument('--paf-list', required=True, help='File containing list of PAF files')
    parser.add_argument('--impg-index', '-i', required=True, help='Path to impg index')
    parser.add_argument('--window', type=int, default=10000, help='Window size (default: 10000)')
    parser.add_argument('--step', type=int, default=10000, help='Step size (default: 10000)')
    parser.add_argument('--max', type=int, default=500000, help='Maximum distance from telomere (default: 500000)')
    parser.add_argument('--min-identity', type=float, default=0.98, help='Minimum alignment identity (default: 0.98)')
    parser.add_argument('--min-output-length', type=int, default=5000, help='Minimum output length (default: 5000)')
    parser.add_argument('--min-diff-chrs', type=int, default=2,
                        help='Minimum different chromosomes with signal (default: 2)')
    parser.add_argument('--min-count', type=int, default=1,
                        help='Minimum count to consider a chromosome present (default: 1)')
    parser.add_argument('--min-consecutive-stop', type=int, default=1,
                        help='Consecutive windows failing criteria to stop (default: 1)')
    parser.add_argument('--threads', type=int, default=1,
                        help='Parallel sequences to process (default: 1)')
    parser.add_argument('--per-window-output', default=None,
                        help='Optional: output per-window chromosome counts TSV')
    parser.add_argument('--fai', default=None,
                        help='FAI index of input sequences. When provided, per-sequence '
                             'max is capped to actual sequence length (prevents out-of-bounds '
                             'coordinates when sequences are shorter than --max, e.g. after '
                             'telomere trimming).')
    parser.add_argument('--verbose', action='store_true',
                        help='Print progress and query counts to stderr')
    return parser.parse_args()

def _strip_pansn(name):
    """Strip PanSN prefix (e.g., 'B6#1#chr1_parm:...' → 'chr1_parm:...').
    Also handles non-PanSN names (e.g., 'chr1_parm:...' stays unchanged)."""
    if '#' in name:
        return name.split('#')[-1]
    return name

def get_self_chr(seq_name):
    """Extract chromosome from sequence name.
    Flank names end with _chrN_parm or _chrN_qarm; extract chrN from the suffix
    (not the contig name, which may itself contain 'chr' in WashU/CHM13 assemblies)."""
    parts = seq_name.rsplit('_', 2)
    if len(parts) >= 3 and parts[-1] in ('parm', 'qarm') and parts[-2].startswith('chr'):
        return parts[-2]
    # Fallback for non-standard names
    stripped = _strip_pansn(seq_name)
    for p in stripped.split('_'):
        if p.startswith('chr'):
            return p
    return None

def get_arm(seq_name):
    """Extract arm (parm or qarm) from sequence name"""
    if '_parm' in seq_name:
        return 'parm'
    elif '_qarm' in seq_name:
        return 'qarm'
    return None

def get_sequences(paf_list, impg_index):
    """Get list of sequences from impg stats"""
    cmd = ['impg', 'stats', '--alignment-list', paf_list, '-i', impg_index, '--list-sequences', '-v', '0']
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"Error running impg stats: {result.stderr}")

    sequences = []
    for line in result.stdout.strip().split('\n')[1:]:  # Skip header
        if line:
            seq_name = line.split('\t')[0]
            sequences.append(seq_name)
    return sequences

def count_chrs_with_signal(chr_counts, min_count):
    """Return set of chromosomes with count >= min_count"""
    return {chr_name for chr_name, count in chr_counts.items() if count >= min_count}

def get_target_chrarm(target_name):
    """Extract chromosome+arm (e.g., 'chr10q') from target name.
    Flank names end with _chrN_parm or _chrN_qarm; extract from the suffix
    (not the contig name, which may itself contain 'chr' in WashU/CHM13 assemblies)."""
    parts = target_name.rsplit('_', 2)
    if len(parts) >= 3:
        arm_part = parts[-1].split(':')[0]  # Handle 'parm:coords' → 'parm'
        chr_part = parts[-2]
        if arm_part in ('parm', 'qarm') and chr_part.startswith('chr'):
            arm = 'p' if arm_part == 'parm' else 'q'
            return f"{chr_part}{arm}"
    # Fallback for non-standard names
    stripped = _strip_pansn(target_name)
    chr_name = None
    arm = None
    for p in stripped.split('_'):
        if p.startswith('chr'):
            chr_name = p
        p_bare = p.split(':')[0]
        if p_bare == 'parm':
            arm = 'p'
        elif p_bare == 'qarm':
            arm = 'q'
    if chr_name and arm:
        return f"{chr_name}{arm}"
    return chr_name


def query_impg_window(paf_list, impg_index, seq_name, start, end, min_identity, min_output_length):
    """Query impg for a single window, return (chr_counts, chrarm_counts) dicts."""
    region = f"{seq_name}:{start}-{end}"
    cmd = [
        'impg', 'query',
        '--alignment-list', paf_list,
        '-i', impg_index,
        '-r', region,
        '--min-identity', str(min_identity),
        '--min-output-length', str(min_output_length),
        '--threads', '1',
        '-v', '0'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return defaultdict(int), defaultdict(int)

    chr_counts = defaultdict(int)
    chrarm_counts = defaultdict(int)
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        if len(fields) < 1:
            continue
        target_name = fields[0]
        # Chromosome-level counts (for early stopping logic)
        # Extract chr from _chrN_parm/_qarm suffix (not contig name)
        t_parts = target_name.rsplit('_', 2)
        chr_found = False
        if len(t_parts) >= 3 and t_parts[-1].split(':')[0] in ('parm', 'qarm') and t_parts[-2].startswith('chr'):
            chr_counts[t_parts[-2]] += 1
            chr_found = True
        if not chr_found:
            # Fallback for non-standard names
            stripped = _strip_pansn(target_name)
            for p in stripped.split('_'):
                if p.startswith('chr'):
                    chr_counts[p] += 1
                    break
        # Chromosome+arm counts (for output)
        chrarm = get_target_chrarm(target_name)
        if chrarm:
            chrarm_counts[chrarm] += 1

    return chr_counts, chrarm_counts


def process_sequence(seq_name, paf_list, args, seq_lengths=None):
    """
    Process a single sequence with early stopping.
    Returns (seq_name, arm, self_chr, region_start, region_end, chrs_involved, arms_involved, query_count, window_details)
    """
    self_chr = get_self_chr(seq_name)
    arm = get_arm(seq_name)

    if not self_chr or not arm:
        return None

    # Per-sequence max: use actual sequence length if FAI was provided
    seq_max = args.max
    if seq_lengths and seq_name in seq_lengths:
        actual_len = seq_lengths[seq_name]
        if actual_len < seq_max:
            if args.verbose:
                print(f"  {seq_name}: capping max from {seq_max} to {actual_len} (actual length)", file=sys.stderr)
            seq_max = actual_len

    # Generate windows
    if seq_max < args.window:
        if args.verbose:
            print(f"  {seq_name}: sequence too short ({seq_max} < window {args.window}), skipping", file=sys.stderr)
        return None
    num_windows = (seq_max - args.window) // args.step + 1
    all_windows = [(i * args.step, i * args.step + args.window) for i in range(num_windows)]

    # Sort windows: parm from start (0, 10000, ...), qarm from end (..., 10000, 0)
    if arm == 'qarm':
        all_windows = all_windows[::-1]

    region_windows = []
    all_chrs_involved = set()
    all_arms_involved = set()
    consecutive_fail = 0
    pending_windows = []
    query_count = 0
    window_details = []

    for start, end in all_windows:
        chr_counts, chrarm_counts = query_impg_window(
            paf_list, args.impg_index, seq_name, start, end,
            args.min_identity, args.min_output_length
        )
        query_count += 1

        # Early stopping uses chromosome-level counts (unchanged behavior)
        chrs_present = count_chrs_with_signal(chr_counts, args.min_count)
        arms_present = count_chrs_with_signal(chrarm_counts, args.min_count)
        num_diff_chrs = len(chrs_present)
        passes = num_diff_chrs >= args.min_diff_chrs

        # Collect per-window detail
        if args.per_window_output:
            # Distance from telomere: parm telomere at 0, qarm telomere at max
            if arm == 'parm':
                dist = start
            else:
                dist = seq_max - end
            window_details.append({
                'start': start, 'end': end,
                'dist_from_telomere': dist,
                'n_chrs': len(chrs_present),
                'n_arms': len(arms_present),
                'chrs': ','.join(sorted(chrs_present, key=lambda x: (len(x), x))) if chrs_present else '.',
                'arms': ','.join(sorted(arms_present, key=lambda x: (len(x), x))) if arms_present else '.',
                'passes': passes
            })

        if passes:
            # Passes criteria - add pending windows and this one
            region_windows.extend(pending_windows)
            for pw in pending_windows:
                all_chrs_involved.update(pw['chrs'])
                all_arms_involved.update(pw.get('arms', set()))
            pending_windows = []
            consecutive_fail = 0
            region_windows.append({'start': start, 'end': end, 'chrs': chrs_present, 'arms': arms_present})
            all_chrs_involved.update(chrs_present)
            all_arms_involved.update(arms_present)
        else:
            consecutive_fail += 1
            pending_windows.append({'start': start, 'end': end, 'chrs': chrs_present, 'arms': arms_present})
            if consecutive_fail >= args.min_consecutive_stop:
                break

    # Compute final region
    if region_windows:
        if arm == 'parm':
            region_start = region_windows[0]['start']
            region_end = region_windows[-1]['end']
        else:
            region_start = region_windows[-1]['start']
            region_end = region_windows[0]['end']
        chrs_str = ','.join(sorted(all_chrs_involved, key=lambda x: (len(x), x)))
        arms_str = ','.join(sorted(all_arms_involved, key=lambda x: (len(x), x)))
    else:
        region_start = '.'
        region_end = '.'
        chrs_str = '.'
        arms_str = '.'

    return (seq_name, arm, self_chr, region_start, region_end, chrs_str, arms_str, query_count, window_details)

def main():
    args = parse_args()

    # Load sequence lengths from FAI if provided
    seq_lengths = None
    if args.fai:
        seq_lengths = {}
        with open(args.fai) as f:
            for line in f:
                parts = line.strip().split('\t')
                seq_lengths[parts[0]] = int(parts[1])
        if args.verbose:
            print(f"Loaded {len(seq_lengths)} sequence lengths from {args.fai}", file=sys.stderr)
            n_short = sum(1 for v in seq_lengths.values() if v < args.max)
            if n_short > 0:
                print(f"  {n_short} sequences shorter than --max ({args.max})", file=sys.stderr)

    # Get sequences
    if args.verbose:
        print("Getting sequence list from impg...", file=sys.stderr)
    sequences = get_sequences(args.paf_list, args.impg_index)

    if args.verbose:
        print(f"Found {len(sequences)} sequences to process", file=sys.stderr)

    # Print header
    print("seq\tarm\tself_chr\tregion_start\tregion_end\tchrs_involved\tarms_involved")

    total_queries = 0
    results = []

    if args.threads == 1:
        # Sequential processing
        for i, seq_name in enumerate(sequences):
            result = process_sequence(seq_name, args.paf_list, args, seq_lengths)
            if result:
                results.append(result)
                total_queries += result[7]
            if args.verbose and (i + 1) % 10 == 0:
                print(f"Processed {i + 1}/{len(sequences)} sequences...", file=sys.stderr)
    else:
        # Parallel processing
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures = {
                executor.submit(process_sequence, seq_name, args.paf_list, args, seq_lengths): seq_name
                for seq_name in sequences
            }

            completed = 0
            for future in as_completed(futures):
                result = future.result()
                if result:
                    results.append(result)
                    total_queries += result[7]
                completed += 1
                if args.verbose and completed % 10 == 0:
                    print(f"Processed {completed}/{len(sequences)} sequences...", file=sys.stderr)

    # Sort and print results
    results.sort(key=lambda x: x[0])  # Sort by sequence name
    for result in results:
        seq_name, arm, self_chr, region_start, region_end, chrs_str, arms_str, _, _ = result
        print(f"{seq_name}\t{arm}\t{self_chr}\t{region_start}\t{region_end}\t{chrs_str}\t{arms_str}")

    # Write per-window output if requested
    if args.per_window_output:
        with open(args.per_window_output, 'w') as f:
            f.write("seq\tarm\tself_chr\twindow_start\twindow_end\tdist_from_telomere\tn_chrs\tn_arms\tchrs_involved\tarms_involved\tpasses_threshold\n")
            for result in results:
                seq_name, arm, self_chr, _, _, _, _, _, window_details = result
                for w in window_details:
                    f.write(f"{seq_name}\t{arm}\t{self_chr}\t{w['start']}\t{w['end']}\t{w['dist_from_telomere']}\t{w['n_chrs']}\t{w['n_arms']}\t{w['chrs']}\t{w['arms']}\t{w['passes']}\n")
        if args.verbose:
            total_windows = sum(len(r[8]) for r in results)
            print(f"Per-window output: {total_windows} windows written to {args.per_window_output}", file=sys.stderr)

    if args.verbose:
        max_queries = len(sequences) * ((args.max - args.window) // args.step + 1)
        savings = (1 - total_queries / max_queries) * 100 if max_queries > 0 else 0
        print(f"\nTotal impg queries: {total_queries}", file=sys.stderr)
        print(f"Maximum possible queries: {max_queries}", file=sys.stderr)
        print(f"Queries saved: {max_queries - total_queries} ({savings:.1f}%)", file=sys.stderr)

if __name__ == '__main__':
    main()
