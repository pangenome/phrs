#!/usr/bin/env python3
"""
Contig Classification Script: p-arms, q-arms, and pq-arms with Telomere-based Filtering

Classifies contigs from PAF alignment files based on their alignment to
p-arms and q-arms of CHM13 chromosomes. Only considers alignments where
the query contig's chromosome of origin matches the target CHM13 chromosome
(using a chromosome alias file) and alignments meet minimum length requirements.

Classifications:
- p-contig: Aligns only to p-arm regions of its corresponding chromosome
- q-contig: Aligns only to q-arm regions of its corresponding chromosome
- pq-contig: Aligns to both p-arm and q-arm regions of its corresponding chromosome
- unclassified: No alignments to p-arms or q-arms

Filtering Requirements:
1. Minimum contig length: All contigs must be >= min_contig_length (default 1Mb)
2. Minimum alignment length: Only count alignments >= min_alignment_length (default 0, no filtering)
   Note: PAF files typically have small alignment blocks (~50kb average), so default is 0
3. Actual telomere presence (from .telo.tsv files):
   - p-contigs: Must have p-telomere (or q-telomere if reverse strand)
   - q-contigs: Must have q-telomere (or p-telomere if reverse strand)
   - pq-contigs: Must have at least 1 telomere (p or q)
4. No multiple telomeres of same type (assembly errors)
5. No duplicate telomeres per SAMPLE#HAPLOTYPE + chromosome
"""

import argparse
import sys
import glob
from collections import defaultdict
from typing import Dict, Set, Tuple, List


def load_chrom_alias(chrom_alias_file: str) -> Dict[str, str]:
    """
    Load chromosome alias mapping file.

    Args:
        chrom_alias_file: Path to chromAlias file

    Returns:
        Dictionary: {full_contig_name: canonical_chr}
    """
    chrom_alias = {}

    try:
        with open(chrom_alias_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                full_name = parts[0]
                canonical_chr = parts[1]
                chrom_alias[full_name] = canonical_chr
    except FileNotFoundError:
        print(f"Error: chromAlias file not found: {chrom_alias_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading chromAlias file: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded chromosome aliases for {len(chrom_alias):,} contigs", file=sys.stderr)
    return chrom_alias


def load_arm_coordinates(p_arms_bed: str, q_arms_bed: str) -> Dict[str, Dict[str, int]]:
    """
    Load p-arm and q-arm coordinates from BED files.

    Args:
        p_arms_bed: Path to p-arms BED file
        q_arms_bed: Path to q-arms BED file

    Returns:
        Dictionary: {chr: {'p_end': int, 'q_start': int, 'q_end': int}}
    """
    arm_coords = {}

    # Load p-arm coordinates
    try:
        with open(p_arms_bed, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                chr_name = parts[0]
                start = int(parts[1])
                end = int(parts[2])

                if chr_name not in arm_coords:
                    arm_coords[chr_name] = {}
                arm_coords[chr_name]['p_start'] = start
                arm_coords[chr_name]['p_end'] = end
    except FileNotFoundError:
        print(f"Error: p-arms file not found: {p_arms_bed}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading p-arms file: {e}", file=sys.stderr)
        sys.exit(1)

    # Load q-arm coordinates
    try:
        with open(q_arms_bed, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                chr_name = parts[0]
                start = int(parts[1])
                end = int(parts[2])

                if chr_name not in arm_coords:
                    arm_coords[chr_name] = {}
                arm_coords[chr_name]['q_start'] = start
                arm_coords[chr_name]['q_end'] = end
    except FileNotFoundError:
        print(f"Error: q-arms file not found: {q_arms_bed}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading q-arms file: {e}", file=sys.stderr)
        sys.exit(1)

    # Verify all chromosomes have both p and q arms
    for chr_name, coords in arm_coords.items():
        if 'p_end' not in coords or 'q_start' not in coords:
            print(f"Warning: Chromosome {chr_name} missing p-arm or q-arm coordinates",
                  file=sys.stderr)

    print(f"Loaded arm coordinates for {len(arm_coords)} chromosomes", file=sys.stderr)
    return arm_coords


def load_t2t_assignments(t2t_dir: str) -> Set[str]:
    """
    Load T2T (telomere-to-telomere) contig assignments.

    Args:
        t2t_dir: Directory containing T2T assignment TSV files

    Returns:
        Set of contig names that are T2T
    """
    t2t_files = glob.glob(f"{t2t_dir}/*.tsv")
    if not t2t_files:
        print(f"Warning: No T2T assignment files found in {t2t_dir}", file=sys.stderr)
        return set()

    t2t_contigs = set()

    for t2t_file in t2t_files:
        try:
            with open(t2t_file, 'r') as f:
                for line in f:
                    if not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) > 0:
                        contig_name = parts[0]
                        t2t_contigs.add(contig_name)
        except Exception as e:
            print(f"Warning: Error reading T2T file {t2t_file}: {e}", file=sys.stderr)
            continue

    print(f"Loaded {len(t2t_contigs):,} T2T contig assignments from {len(t2t_files)} files", file=sys.stderr)
    return t2t_contigs


def load_telomere_calls(telomere_dir: str) -> Dict[str, List[Tuple[int, int, int]]]:
    """
    Load telomere calls from all .telo.tsv files.

    Args:
        telomere_dir: Directory containing .telo.tsv files

    Returns:
        Dict mapping contig_name -> [(start, end, contig_length), ...]
    """
    telomere_files = glob.glob(f"{telomere_dir}/*.telo.tsv")
    if not telomere_files:
        print(f"Warning: No telomere files found in {telomere_dir}", file=sys.stderr)
        return {}

    telomere_calls = defaultdict(list)
    total_entries = 0

    for telo_file in telomere_files:
        try:
            with open(telo_file, 'r') as f:
                for line in f:
                    if not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 4:
                        continue

                    contig_name = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    contig_length = int(parts[3])

                    telomere_calls[contig_name].append((start, end, contig_length))
                    total_entries += 1

        except Exception as e:
            print(f"Warning: Error reading telomere file {telo_file}: {e}", file=sys.stderr)
            continue

    print(f"Loaded {total_entries:,} telomere entries from {len(telomere_files)} files", file=sys.stderr)
    return dict(telomere_calls)


def determine_contig_strand(contig_name: str, contig_data: Dict) -> Tuple[str, float]:
    """
    Determine predominant strand orientation for a contig.

    Uses weighted majority vote based on alignment block lengths.

    Args:
        contig_name: Name of the contig
        contig_data: Dictionary containing 'strand_blocks' key with list of (strand, block_length)

    Returns:
        (strand, confidence) where:
        - strand: '+', '-', or 'mixed'
        - confidence: 0.0-1.0 (fraction of aligned bases in majority strand)
    """
    strand_blocks = contig_data.get('strand_blocks', [])
    if not strand_blocks:
        return ('mixed', 0.0)

    forward_bases = sum(length for strand, length in strand_blocks if strand == '+')
    reverse_bases = sum(length for strand, length in strand_blocks if strand == '-')
    total_bases = forward_bases + reverse_bases

    if total_bases == 0:
        return ('mixed', 0.0)

    if forward_bases >= reverse_bases:
        confidence = forward_bases / total_bases
        strand = '+' if confidence >= 0.75 else 'mixed'
    else:
        confidence = reverse_bases / total_bases
        strand = '-' if confidence >= 0.75 else 'mixed'

    return (strand, confidence)


def classify_telomere_position(telomere_start: int, contig_length: int, threshold: int = 100000) -> str:
    """
    Classify telomere as 'p-telomere', 'q-telomere', or 'internal'.

    Args:
        telomere_start: Start position of telomere
        contig_length: Total contig length
        threshold: Distance from ends to consider telomeric (default 100kb)

    Returns:
        'p-telomere' if near start, 'q-telomere' if near end, 'internal' otherwise
    """
    if telomere_start <= threshold:
        return 'p-telomere'
    elif (contig_length - telomere_start) <= threshold:
        return 'q-telomere'
    else:
        return 'internal'


def validate_telomere_placement(
    classification: str,
    telomeres: List[Tuple[int, int, int]],
    strand: str,
    confidence: float
) -> Dict:
    """
    Validate that telomere positions match expected positions for classification.

    Args:
        classification: Contig classification ('p-contig', 'q-contig', 'pq-contig')
        telomeres: List of (start, end, contig_length) tuples
        strand: Predominant strand ('+', '-', 'mixed')
        confidence: Strand confidence (0.0-1.0)

    Returns:
        Dictionary with validation results
    """
    # Classify each telomere position
    p_telomeres = []
    q_telomeres = []

    for start, end, contig_length in telomeres:
        telo_type = classify_telomere_position(start, contig_length)
        if telo_type == 'p-telomere':
            p_telomeres.append((start, end, contig_length))
        elif telo_type == 'q-telomere':
            q_telomeres.append((start, end, contig_length))

    num_p = len(p_telomeres)
    num_q = len(q_telomeres)

    # FAIL if multiple telomeres of same type (assembly problem)
    if num_p > 1:
        return {
            'has_p_telomere': True,
            'has_q_telomere': num_q > 0,
            'num_p_telomeres': num_p,
            'num_q_telomeres': num_q,
            'telomere_positions': [t[0] for t in telomeres],
            'validation_status': 'fail',
            'validation_message': 'multiple_p_telomeres'
        }

    if num_q > 1:
        return {
            'has_p_telomere': num_p > 0,
            'has_q_telomere': True,
            'num_p_telomeres': num_p,
            'num_q_telomeres': num_q,
            'telomere_positions': [t[0] for t in telomeres],
            'validation_status': 'fail',
            'validation_message': 'multiple_q_telomeres'
        }

    # Skip validation if strand is too ambiguous
    if confidence < 0.75:
        return {
            'has_p_telomere': num_p > 0,
            'has_q_telomere': num_q > 0,
            'num_p_telomeres': num_p,
            'num_q_telomeres': num_q,
            'telomere_positions': [t[0] for t in telomeres],
            'validation_status': 'skip',
            'validation_message': 'mixed_strand'
        }

    # Validate based on classification and strand
    status = 'pass'
    message = 'valid'

    if classification == 'p-contig':
        if strand == '+':
            expected = num_p > 0
            status = 'pass' if expected else 'fail'
            message = 'has_p_telomere' if expected else 'missing_p_telomere'
        else:  # strand == '-'
            expected = num_q > 0
            status = 'pass' if expected else 'fail'
            message = 'has_q_telomere_reversed' if expected else 'missing_q_telomere_reversed'

    elif classification == 'q-contig':
        if strand == '+':
            expected = num_q > 0
            status = 'pass' if expected else 'fail'
            message = 'has_q_telomere' if expected else 'missing_q_telomere'
        else:  # strand == '-'
            expected = num_p > 0
            status = 'pass' if expected else 'fail'
            message = 'has_p_telomere_reversed' if expected else 'missing_p_telomere_reversed'

    elif classification == 'pq-contig':
        # pq-contigs should have at least one telomere
        has_any_telomere = (num_p > 0) or (num_q > 0)
        status = 'pass' if has_any_telomere else 'warning'
        if num_p > 0 and num_q > 0:
            message = 'has_both_telomeres'
        elif num_p > 0:
            message = 'has_p_telomere_only'
        elif num_q > 0:
            message = 'has_q_telomere_only'
        else:
            message = 'no_telomeres'

    return {
        'has_p_telomere': num_p > 0,
        'has_q_telomere': num_q > 0,
        'num_p_telomeres': num_p,
        'num_q_telomeres': num_q,
        'telomere_positions': [t[0] for t in telomeres],
        'validation_status': status,
        'validation_message': message
    }


def check_and_fail_duplicate_telomeres(
    classifications: Dict,
    output_file: str,
    t2t_contigs: Set[str]
) -> int:
    """
    Check for duplicate telomeres per SAMPLE#HAPLOTYPE + chromosome.
    Use smart resolution: prioritize T2T, then proximity, then both telomeres, then length.

    Args:
        classifications: Dictionary of contig classifications
        output_file: Base output file path
        t2t_contigs: Set of T2T contig names

    Returns:
        Number of contigs failed due to duplicates
    """
    # Group by (sample, haplotype, chromosome)
    groups = defaultdict(lambda: {'p_telomeres': [], 'q_telomeres': []})

    for contig, info in classifications.items():
        if not info.get('telomeric_filter'):
            continue

        # Parse contig name: SAMPLE#HAPLOTYPE#ACCESSION
        parts = contig.split('#')
        if len(parts) < 3:
            continue

        sample_hap = f"{parts[0]}#{parts[1]}"
        chrom = info.get('target_chr')

        if not chrom:
            continue

        key = (sample_hap, chrom)

        # Track telomere types
        if info.get('has_p_telomere'):
            groups[key]['p_telomeres'].append(contig)
        if info.get('has_q_telomere'):
            groups[key]['q_telomeres'].append(contig)

    # Smart duplicate resolution
    duplicates_file = output_file.replace('.tsv', '_duplicates.tsv')
    failed_count = 0
    rescued_count = 0

    with open(duplicates_file, 'w') as f:
        f.write("sample_haplotype\tchromosome\ttelomere_type\t"
               "kept_contig\tfailed_contigs\tresolution_method\n")

        for key, data in sorted(groups.items()):
            sample_hap, chrom = key

            # Handle p-telomere duplicates
            if len(data['p_telomeres']) > 1:
                contigs = data['p_telomeres']
                kept, failed, method = resolve_duplicates(contigs, classifications, t2t_contigs, telomere_type='p')

                # Fail the duplicates
                for contig in failed:
                    classifications[contig]['telomeric_filter'] = False
                    classifications[contig]['filter_reason'] = 'duplicate_p_telomere'
                    failed_count += 1

                if kept:
                    rescued_count += 1

                f.write(f"{sample_hap}\t{chrom}\tp-telomere\t"
                       f"{kept if kept else 'NONE'}\t{','.join(sorted(failed))}\t{method}\n")

            # Handle q-telomere duplicates
            if len(data['q_telomeres']) > 1:
                contigs = data['q_telomeres']
                kept, failed, method = resolve_duplicates(contigs, classifications, t2t_contigs, telomere_type='q')

                # Fail the duplicates
                for contig in failed:
                    classifications[contig]['telomeric_filter'] = False
                    classifications[contig]['filter_reason'] = 'duplicate_q_telomere'
                    failed_count += 1

                if kept:
                    rescued_count += 1

                f.write(f"{sample_hap}\t{chrom}\tq-telomere\t"
                       f"{kept if kept else 'NONE'}\t{','.join(sorted(failed))}\t{method}\n")

    if failed_count > 0:
        print(f"\nDuplicate telomere resolution:",
              file=sys.stderr)
        print(f"  Rescued (kept best): {rescued_count}",
              file=sys.stderr)
        print(f"  Failed (duplicates): {failed_count}",
              file=sys.stderr)
        print(f"  Details written to: {duplicates_file}", file=sys.stderr)

    return failed_count


def resolve_duplicates(contigs: List[str], classifications: Dict, t2t_contigs: Set[str], telomere_type: str = 'p') -> Tuple[str, List[str], str]:
    """
    Smart duplicate resolution strategy.

    Priority:
    1. Keep T2T contig (if only one is T2T)
    2. Keep contig that aligns closest to chromosome end (most telomeric)
    3. Keep contig with both telomeres (if only one has both)
    4. Among contigs with both telomeres, keep longest
    5. If no contig has both, keep longest
    6. If top candidates are similar (within 10% proximity AND length), fail all

    Args:
        contigs: List of duplicate contig names
        classifications: Dictionary of classifications
        t2t_contigs: Set of T2T contig names
        telomere_type: 'p' for p-telomere duplicates, 'q' for q-telomere duplicates

    Returns:
        (kept_contig, failed_contigs, resolution_method)
    """
    # Get info for all contigs
    contig_info = []
    for contig in contigs:
        info = classifications[contig]
        is_t2t = contig in t2t_contigs
        has_both = info.get('has_p_telomere', False) and info.get('has_q_telomere', False)
        length = info.get('contig_length', 0)

        # Get alignment proximity to chromosome end
        # For p-telomeres: lower min_target_start is better (closer to chr start)
        # For q-telomeres: higher max_target_end is better (closer to chr end)
        if telomere_type == 'p':
            proximity = info.get('min_target_start', float('inf'))
        else:  # q-telomere
            # For q, we want the smallest distance from end, which means largest max_target_end
            # We'll negate it so smaller is still better for sorting
            proximity = -info.get('max_target_end', 0)

        contig_info.append((contig, is_t2t, has_both, length, proximity))

    # Strategy 1: Prioritize T2T contigs
    t2t_contigs_list = [c for c in contig_info if c[1]]  # c[1] is is_t2t
    non_t2t_contigs = [c for c in contig_info if not c[1]]

    if t2t_contigs_list:
        # If only one is T2T, keep it
        if len(t2t_contigs_list) == 1:
            kept = t2t_contigs_list[0][0]
            failed = [c[0] for c in contig_info if c[0] != kept]
            return (kept, failed, 't2t_priority')

        # Multiple T2T - resolve among them by proximity
        contig_info = t2t_contigs_list

    # Strategy 2: Prioritize contig that reaches closest to chromosome end
    sorted_by_proximity = sorted(contig_info, key=lambda x: x[4])  # x[4] is proximity
    closest = sorted_by_proximity[0]
    second_closest = sorted_by_proximity[1] if len(sorted_by_proximity) > 1 else None

    # If one contig is clearly closer to the end (>100kb threshold difference)
    if second_closest:
        proximity_diff = abs(closest[4] - second_closest[4])
        if proximity_diff > 100000:  # >100kb difference means clear winner
            kept = closest[0]
            failed = [c[0] for c in contig_info if c[0] != kept]
            if non_t2t_contigs:
                failed.extend([c[0] for c in non_t2t_contigs])
            return (kept, failed, 'closest_to_telomere')

    # Strategy 3: Prioritize contigs with both telomeres
    with_both = [c for c in contig_info if c[2]]  # c[2] is has_both
    without_both = [c for c in contig_info if not c[2]]

    if with_both:
        # If only one has both telomeres, keep it
        if len(with_both) == 1:
            kept = with_both[0][0]
            failed = [c[0] for c in without_both]
            if non_t2t_contigs:
                failed.extend([c[0] for c in non_t2t_contigs])
            return (kept, failed, 'has_both_telomeres')

        # Multiple have both - keep longest
        with_both_sorted = sorted(with_both, key=lambda x: x[3], reverse=True)  # x[3] is length
        longest = with_both_sorted[0]
        second_longest = with_both_sorted[1] if len(with_both_sorted) > 1 else None

        # Check if genuinely ambiguous (within 10% length)
        if second_longest and second_longest[3] >= longest[3] * 0.9:
            # Too similar - fail all
            failed = [c[0] for c in contig_info]
            if non_t2t_contigs:
                failed.extend([c[0] for c in non_t2t_contigs])
            return (None, failed, 'ambiguous_similar_lengths')

        kept = longest[0]
        failed = [c[0] for c in contig_info if c[0] != kept]
        if non_t2t_contigs:
            failed.extend([c[0] for c in non_t2t_contigs])
        return (kept, failed, 'longest_with_both_telomeres')

    # Strategy 4: No contig has both telomeres - keep longest
    sorted_contigs = sorted(contig_info, key=lambda x: x[3], reverse=True)  # x[3] is length
    longest = sorted_contigs[0]
    second_longest = sorted_contigs[1] if len(sorted_contigs) > 1 else None

    # Check if genuinely ambiguous (within 10% length)
    if second_longest and second_longest[3] >= longest[3] * 0.9:  # x[3] is length
        # Too similar - fail all
        failed = [c[0] for c in contig_info]
        if non_t2t_contigs:
            failed.extend([c[0] for c in non_t2t_contigs])
        return (None, failed, 'ambiguous_similar_lengths')

    kept = longest[0]
    failed = [c[0] for c in sorted_contigs if c[0] != kept]
    if non_t2t_contigs:
        failed.extend([c[0] for c in non_t2t_contigs])
    return (kept, failed, 'longest_contig')


def add_chm13_contigs(classifications: Dict, arm_coords: Dict) -> Dict:
    """
    Add CHM13 reference contigs to classifications as valid T2T contigs.
    CHM13 is the T2T reference genome, so all its chromosomes (except chrY)
    are complete telomere-to-telomere assemblies with both telomeres.

    Args:
        classifications: Dictionary of contig classifications
        arm_coords: Dictionary of chromosome arm coordinates (to get chr lengths)

    Returns:
        Updated classifications dictionary with CHM13 contigs added
    """
    # CHM13 chromosomes to add (all except chrY)
    chm13_chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

    for chrom in chm13_chromosomes:
        chrom_full = f"CHM13#0#{chrom}"

        # Get chromosome length from arm coordinates
        p_arm = arm_coords.get(chrom_full, {})
        q_arm_key = chrom_full
        # Calculate total length from p-arm start to q-arm end
        p_start = p_arm.get('start', 0)
        q_end = arm_coords.get(chrom_full, {}).get('end', 0)

        # For CHM13, we can estimate length from the q-arm end coordinate
        # since the reference goes from 0 to the end of q-arm
        contig_length = q_end if q_end > 0 else 0

        classifications[chrom_full] = {
            'classification': 'pq-contig',
            'p_chrs': {chrom_full},
            'q_chrs': {chrom_full},
            'contig_length': contig_length,
            'telomeric_filter': True,
            'filter_reason': 'pass',
            'has_p_telomere': True,
            'has_q_telomere': True,
            'num_p_telomeres': 1,
            'num_q_telomeres': 1,
            'strand': '+',
            'strand_confidence': 1.0,
            'target_chr': chrom_full,
            'min_target_start': 0,
            'max_target_end': contig_length,
            'validation_status': 'pass',
            'validation_message': 'chm13_reference_t2t'
        }

    print(f"Added {len(chm13_chromosomes)} CHM13 reference contigs as valid T2T", file=sys.stderr)
    return classifications


def process_paf_files(paf_pattern: str, arm_coords: Dict[str, Dict[str, int]], chrom_alias: Dict[str, str], min_alignment_length: int) -> Dict[str, Dict]:
    """
    Process all PAF files and track which contigs align to p-arms and q-arms.
    Only considers alignments where the query contig's chromosome matches the target chromosome.
    Also tracks telomeric proximity metrics for filtering.

    Args:
        paf_pattern: Glob pattern for PAF files
        arm_coords: Dictionary of chromosome arm coordinates
        chrom_alias: Mapping of full contig names to canonical chromosomes
        min_alignment_length: Minimum alignment length to count toward classification

    Returns:
        Dictionary: {contig: {'p_arms': set(chrs), 'q_arms': set(chrs),
                              'contig_length': int, 'min_target_start': int,
                              'max_target_end': int, 'target_chr': str}}
    """
    # Find all PAF files
    paf_files = glob.glob(paf_pattern)
    if not paf_files:
        print(f"Error: No PAF files found matching pattern: {paf_pattern}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(paf_files)} PAF files to process", file=sys.stderr)

    # Track contig alignments with telomeric metrics
    contig_arms = defaultdict(lambda: {
        'p_arms': set(),
        'q_arms': set(),
        'contig_length': 0,
        'min_target_start': float('inf'),
        'max_target_end': 0,
        'target_chr': None,
        'strand_blocks': []  # List of (strand, block_length) tuples
    })

    total_alignments = 0
    processed_files = 0
    skipped_lines = 0
    filtered_alignments = 0
    filtered_by_length = 0

    for paf_file in paf_files:
        try:
            with open(paf_file, 'r') as f:
                for line in f:
                    if not line.strip():
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        skipped_lines += 1
                        continue

                    # PAF columns: 0=query, 1=query_len, 4=strand, 5=target, 7=target_start, 8=target_end
                    contig_name = parts[0]
                    query_len = int(parts[1])
                    strand = parts[4]
                    target_chr = parts[5]
                    target_start = int(parts[7])
                    target_end = int(parts[8])

                    total_alignments += 1

                    # Calculate alignment length
                    alignment_length = target_end - target_start

                    # Skip alignments shorter than threshold
                    if alignment_length < min_alignment_length:
                        filtered_by_length += 1
                        continue

                    # Check if target chromosome has arm coordinates
                    if target_chr not in arm_coords:
                        continue

                    # Only consider alignments where query and target chromosome match
                    query_canonical_chr = chrom_alias.get(contig_name)
                    target_canonical_chr = chrom_alias.get(target_chr)

                    if query_canonical_chr is None or target_canonical_chr is None:
                        filtered_alignments += 1
                        continue

                    if query_canonical_chr != target_canonical_chr:
                        filtered_alignments += 1
                        continue

                    coords = arm_coords[target_chr]

                    # Update contig metrics
                    contig_data = contig_arms[contig_name]
                    contig_data['contig_length'] = query_len
                    contig_data['target_chr'] = target_chr
                    contig_data['min_target_start'] = min(contig_data['min_target_start'], target_start)
                    contig_data['max_target_end'] = max(contig_data['max_target_end'], target_end)

                    # Track strand information
                    block_length = target_end - target_start
                    contig_data['strand_blocks'].append((strand, block_length))

                    # Check if alignment overlaps with p-arm
                    # Alignment overlaps p-arm if alignment_end > p_start and alignment_start < p_end
                    if 'p_end' in coords:
                        p_start = coords.get('p_start', 0)
                        p_end = coords['p_end']
                        if target_start < p_end and target_end > p_start:
                            contig_data['p_arms'].add(target_chr)

                    # Check if alignment overlaps with q-arm
                    # Alignment overlaps q-arm if alignment_end > q_start and alignment_start < q_end
                    if 'q_start' in coords and 'q_end' in coords:
                        q_start = coords['q_start']
                        q_end = coords['q_end']
                        if target_start < q_end and target_end > q_start:
                            contig_data['q_arms'].add(target_chr)

            processed_files += 1
            if processed_files % 50 == 0:
                print(f"Processed {processed_files}/{len(paf_files)} PAF files...", file=sys.stderr)

        except Exception as e:
            print(f"Warning: Error processing {paf_file}: {e}", file=sys.stderr)
            continue

    print(f"\nProcessing complete:", file=sys.stderr)
    print(f"  Total alignments: {total_alignments:,}", file=sys.stderr)
    print(f"  Filtered (too short): {filtered_by_length:,}", file=sys.stderr)
    print(f"  Filtered (chr mismatch): {filtered_alignments:,}", file=sys.stderr)
    print(f"  Processed files: {processed_files}", file=sys.stderr)
    print(f"  Skipped lines: {skipped_lines}", file=sys.stderr)
    print(f"  Unique contigs: {len(contig_arms):,}", file=sys.stderr)

    return dict(contig_arms)


def classify_contigs(contig_arms: Dict[str, Dict[str, Set[str]]]) -> Dict[str, Dict]:
    """
    Classify contigs based on their p-arm and q-arm alignments.

    Args:
        contig_arms: Dictionary of contig alignments to arms

    Returns:
        Dictionary: {contig: {'classification': str, 'p_chrs': set, 'q_chrs': set}}
    """
    classifications = {}

    for contig, arms in contig_arms.items():
        p_arms = arms['p_arms']
        q_arms = arms['q_arms']

        # Determine classification
        has_p = len(p_arms) > 0
        has_q = len(q_arms) > 0

        if has_p and has_q:
            classification = 'pq-contig'
        elif has_p:
            classification = 'p-contig'
        elif has_q:
            classification = 'q-contig'
        else:
            classification = 'unclassified'

        classifications[contig] = {
            'classification': classification,
            'p_chrs': p_arms,
            'q_chrs': q_arms
        }

    return classifications


def apply_telomeric_filter(
    classifications: Dict[str, Dict],
    contig_data: Dict[str, Dict],
    telomere_calls: Dict[str, List[Tuple[int, int, int]]],
    min_contig_length: int
) -> Dict[str, Dict]:
    """
    Apply filtering based on:
    1. Minimum contig length (all contigs)
    2. Actual telomere presence (required)
    3. No multiple telomeres of same type

    Args:
        classifications: Dictionary of contig classifications
        contig_data: Dictionary with contig metrics from PAF processing
        telomere_calls: Dictionary of telomere calls per contig
        min_contig_length: Minimum contig length (bp)

    Returns:
        Enhanced classifications with telomeric_filter and filter_reason fields
    """
    for contig, info in classifications.items():
        classification = info['classification']

        # Get contig data
        data = contig_data.get(contig, {})
        contig_length = data.get('contig_length', 0)

        # FILTER 1: Minimum length (ALL contigs)
        if contig_length < min_contig_length:
            info['telomeric_filter'] = False
            info['filter_reason'] = 'too_short'
            info['contig_length'] = contig_length
            info['target_chr'] = data.get('target_chr', 'NA')
            info['min_target_start'] = data.get('min_target_start', float('inf'))
            info['max_target_end'] = data.get('max_target_end', 0)
            continue

        # Skip unclassified
        if classification == 'unclassified':
            info['telomeric_filter'] = False
            info['filter_reason'] = 'unclassified'
            info['contig_length'] = contig_length
            info['target_chr'] = data.get('target_chr', 'NA')
            info['min_target_start'] = data.get('min_target_start', float('inf'))
            info['max_target_end'] = data.get('max_target_end', 0)
            continue

        # Get telomere validation
        telomeres = telomere_calls.get(contig, [])
        strand, confidence = determine_contig_strand(contig, data)
        validation = validate_telomere_placement(
            classification,
            telomeres,
            strand,
            confidence
        )

        # Add validation info
        info.update(validation)
        info['strand'] = strand
        info['strand_confidence'] = confidence
        info['target_chr'] = data.get('target_chr', 'NA')
        info['contig_length'] = contig_length
        info['min_target_start'] = data.get('min_target_start', float('inf'))
        info['max_target_end'] = data.get('max_target_end', 0)

        # FILTER 2: Multiple telomeres (assembly error)
        if validation.get('num_p_telomeres', 0) > 1 or validation.get('num_q_telomeres', 0) > 1:
            info['telomeric_filter'] = False
            info['filter_reason'] = 'multiple_telomeres'
            continue

        # FILTER 3: Require actual telomere calls
        has_p = validation.get('has_p_telomere', False)
        has_q = validation.get('has_q_telomere', False)

        if classification == 'p-contig':
            # Forward strand: need p-telomere, Reverse: need q-telomere
            if strand == '+':
                passes = has_p
                reason = 'no_p_telomere' if not passes else 'pass'
            elif strand == '-':
                passes = has_q
                reason = 'no_q_telomere_reversed' if not passes else 'pass'
            else:  # mixed strand
                passes = (has_p or has_q)
                reason = 'no_telomere_mixed_strand' if not passes else 'pass'

        elif classification == 'q-contig':
            # Forward strand: need q-telomere, Reverse: need p-telomere
            if strand == '+':
                passes = has_q
                reason = 'no_q_telomere' if not passes else 'pass'
            elif strand == '-':
                passes = has_p
                reason = 'no_p_telomere_reversed' if not passes else 'pass'
            else:  # mixed strand
                passes = (has_p or has_q)
                reason = 'no_telomere_mixed_strand' if not passes else 'pass'

        elif classification == 'pq-contig':
            # Need at least one telomere
            passes = (has_p or has_q)
            reason = 'no_telomeres' if not passes else 'pass'
        else:
            passes = False
            reason = 'unknown_classification'

        info['telomeric_filter'] = passes
        info['filter_reason'] = reason

    return classifications


def write_results(classifications: Dict[str, Dict], output_file: str):
    """
    Write classification results to TSV file.

    Main output: Only valid contigs (passing with at least 1 telomere)
    Debug output: All failed/problematic contigs with detailed information

    Args:
        classifications: Dictionary of contig classifications with telomeric data
        output_file: Path to output TSV file
    """
    try:
        # Write main output file - ONLY passing contigs with at least 1 telomere
        valid_count = 0
        with open(output_file, 'w') as f:
            # Write header
            f.write("contig_name\tclassification\tp_arm_chrs\tq_arm_chrs\tcontig_length\t"
                   "has_p_telomere\thas_q_telomere\t"
                   "num_p_telomeres\tnum_q_telomeres\t"
                   "strand\tstrand_confidence\t"
                   "validation_status\tvalidation_message\n")

            # Write only passing contigs with at least 1 telomere
            for contig in sorted(classifications.keys()):
                info = classifications[contig]

                # Only write if passing filter AND has at least 1 telomere
                passes_filter = info.get('telomeric_filter', False)
                has_telomere = info.get('has_p_telomere', False) or info.get('has_q_telomere', False)

                if passes_filter and has_telomere:
                    classification = info['classification']

                    # Format chromosome sets as comma-separated lists (use NA for empty)
                    p_chrs = ','.join(sorted(info['p_chrs'])) if info['p_chrs'] else 'NA'
                    q_chrs = ','.join(sorted(info['q_chrs'])) if info['q_chrs'] else 'NA'

                    # Get data
                    contig_length = info.get('contig_length', 0)
                    has_p_telo = str(info.get('has_p_telomere', False))
                    has_q_telo = str(info.get('has_q_telomere', False))
                    num_p_telo = info.get('num_p_telomeres', 0)
                    num_q_telo = info.get('num_q_telomeres', 0)
                    strand = info.get('strand', 'NA')
                    strand_conf = info.get('strand_confidence', 0.0)
                    val_status = info.get('validation_status', 'NA')
                    val_message = info.get('validation_message', 'NA')

                    f.write(f"{contig}\t{classification}\t{p_chrs}\t{q_chrs}\t{contig_length}\t"
                           f"{has_p_telo}\t{has_q_telo}\t"
                           f"{num_p_telo}\t{num_q_telo}\t"
                           f"{strand}\t{strand_conf:.3f}\t"
                           f"{val_status}\t{val_message}\n")
                    valid_count += 1

        print(f"\nValid contigs written to: {output_file} ({valid_count} contigs)", file=sys.stderr)

        # Write debug output file - ALL failed/problematic contigs
        debug_file = output_file.replace('.tsv', '_debug.tsv')
        debug_count = 0
        with open(debug_file, 'w') as f:
            f.write("contig_name\tclassification\tp_arm_chrs\tq_arm_chrs\tcontig_length\t"
                   "telomeric_filter\tfilter_reason\t"
                   "has_p_telomere\thas_q_telomere\t"
                   "num_p_telomeres\tnum_q_telomeres\t"
                   "strand\tstrand_confidence\t"
                   "validation_status\tvalidation_message\n")

            # Write all non-passing contigs
            for contig in sorted(classifications.keys()):
                info = classifications[contig]

                # Write if NOT passing OR no telomeres
                passes_filter = info.get('telomeric_filter', False)
                has_telomere = info.get('has_p_telomere', False) or info.get('has_q_telomere', False)

                if not passes_filter or not has_telomere:
                    classification = info['classification']

                    # Format chromosome sets
                    p_chrs = ','.join(sorted(info['p_chrs'])) if info['p_chrs'] else 'NA'
                    q_chrs = ','.join(sorted(info['q_chrs'])) if info['q_chrs'] else 'NA'

                    # Get all data
                    contig_length = info.get('contig_length', 0)
                    telomeric_filter = 'pass' if info.get('telomeric_filter', False) else 'fail'
                    filter_reason = info.get('filter_reason', 'NA')
                    has_p_telo = str(info.get('has_p_telomere', False))
                    has_q_telo = str(info.get('has_q_telomere', False))
                    num_p_telo = info.get('num_p_telomeres', 0)
                    num_q_telo = info.get('num_q_telomeres', 0)
                    strand = info.get('strand', 'NA')
                    strand_conf = info.get('strand_confidence', 0.0)
                    val_status = info.get('validation_status', 'NA')
                    val_message = info.get('validation_message', 'NA')

                    f.write(f"{contig}\t{classification}\t{p_chrs}\t{q_chrs}\t{contig_length}\t"
                           f"{telomeric_filter}\t{filter_reason}\t"
                           f"{has_p_telo}\t{has_q_telo}\t"
                           f"{num_p_telo}\t{num_q_telo}\t"
                           f"{strand}\t{strand_conf:.3f}\t"
                           f"{val_status}\t{val_message}\n")
                    debug_count += 1

        print(f"Failed/problematic contigs written to: {debug_file} ({debug_count} contigs)", file=sys.stderr)

    except Exception as e:
        print(f"Error writing results: {e}", file=sys.stderr)
        sys.exit(1)


def print_statistics(classifications: Dict[str, Dict]):
    """
    Print summary statistics of classifications with telomeric filter breakdown.

    Args:
        classifications: Dictionary of contig classifications with telomeric data
    """
    # Count each classification type and telomeric filter status
    counts = defaultdict(int)
    telomeric_pass = defaultdict(int)
    telomeric_fail = defaultdict(int)
    valid_with_telomere = defaultdict(int)

    for info in classifications.values():
        classification = info['classification']
        counts[classification] += 1

        passes_filter = info.get('telomeric_filter', False)
        has_telomere = info.get('has_p_telomere', False) or info.get('has_q_telomere', False)

        if passes_filter:
            telomeric_pass[classification] += 1
            if has_telomere:
                valid_with_telomere[classification] += 1
        else:
            telomeric_fail[classification] += 1

    total = len(classifications)
    total_telomeric_pass = sum(telomeric_pass.values())
    total_valid = sum(valid_with_telomere.values())

    print("\n" + "="*60, file=sys.stderr)
    print("Classification Summary:", file=sys.stderr)
    print("="*60, file=sys.stderr)

    # Print in specific order with telomeric breakdown
    for cat in ['p-contig', 'q-contig', 'pq-contig', 'unclassified']:
        count = counts[cat]
        percentage = (count / total * 100) if total > 0 else 0
        print(f"{cat:15s} {count:8,d} ({percentage:5.1f}%)", file=sys.stderr)

        # Show telomeric pass/fail breakdown
        if count > 0:
            pass_count = telomeric_pass[cat]
            fail_count = telomeric_fail[cat]
            valid_count = valid_with_telomere[cat]
            pass_pct = (pass_count / count * 100) if count > 0 else 0
            fail_pct = (fail_count / count * 100) if count > 0 else 0
            valid_pct = (valid_count / count * 100) if count > 0 else 0
            print(f"  passed filter:       {pass_count:6,d} ({pass_pct:5.1f}%)", file=sys.stderr)
            print(f"  valid (w/ telomere): {valid_count:6,d} ({valid_pct:5.1f}%)", file=sys.stderr)
            print(f"  failed:              {fail_count:6,d} ({fail_pct:5.1f}%)", file=sys.stderr)

    print("-"*60, file=sys.stderr)
    print(f"{'Total':15s} {total:8,d} (100.0%)", file=sys.stderr)
    telomeric_pass_pct = (total_telomeric_pass / total * 100) if total > 0 else 0
    valid_pct = (total_valid / total * 100) if total > 0 else 0
    print(f"{'Passed Filter':15s} {total_telomeric_pass:8,d} ({telomeric_pass_pct:5.1f}%)", file=sys.stderr)
    print(f"{'Valid (w/ telo)':15s} {total_valid:8,d} ({valid_pct:5.1f}%)", file=sys.stderr)
    print("="*60, file=sys.stderr)


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description='Classify contigs based on alignment to CHM13 p-arms and q-arms with telomeric filtering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --p-arms p_arms.bed --q-arms q_arms.bed \\
           --chrom-alias chromAlias.txt \\
           --paf-pattern '/path/to/*.paf' --output results.tsv \\
           --min-contig-length 1000000 \\
           --telomere-dir /path/to/telomeres

  # Optional: Filter out tiny spurious alignments
  %(prog)s --min-alignment-length 10000 ...
        """
    )

    parser.add_argument(
        '--p-arms',
        required=True,
        help='Path to p-arms BED file'
    )

    parser.add_argument(
        '--q-arms',
        required=True,
        help='Path to q-arms BED file'
    )

    parser.add_argument(
        '--chrom-alias',
        required=True,
        help='Path to chromosome alias file (maps contigs to canonical chromosomes)'
    )

    parser.add_argument(
        '--paf-pattern',
        required=True,
        help='Glob pattern for PAF files (e.g., "/path/to/*.paf")'
    )

    parser.add_argument(
        '--output',
        required=True,
        help='Path to output TSV file'
    )

    parser.add_argument(
        '--min-contig-length',
        type=int,
        default=1000000,
        help='Minimum contig length in bp (default: 1000000)'
    )

    parser.add_argument(
        '--min-alignment-length',
        type=int,
        default=0,
        help='Minimum alignment length to count toward classification in bp (default: 0, no filtering)'
    )

    parser.add_argument(
        '--telomere-dir',
        default='/moosefs/pangenomes/HPRCv2/telomeres',
        help='Directory containing telomere call TSV files (default: /moosefs/pangenomes/HPRCv2/telomeres)'
    )

    parser.add_argument(
        '--t2t-dir',
        default='/moosefs/guarracino/HPRCv2/coverage/t2t-assignments',
        help='Directory containing T2T assignment TSV files (default: /moosefs/guarracino/HPRCv2/coverage/t2t-assignments)'
    )

    args = parser.parse_args()

    # Step 1: Load chromosome aliases
    print("\n[1/7] Loading chromosome aliases...", file=sys.stderr)
    chrom_alias = load_chrom_alias(args.chrom_alias)

    # Step 2: Load arm coordinates
    print("\n[2/7] Loading arm coordinates...", file=sys.stderr)
    arm_coords = load_arm_coordinates(args.p_arms, args.q_arms)

    # Step 3: Load T2T assignments
    print("\n[3/7] Loading T2T assignments...", file=sys.stderr)
    t2t_contigs = load_t2t_assignments(args.t2t_dir)

    # Step 4: Load telomere calls
    print("\n[4/7] Loading telomere calls...", file=sys.stderr)
    telomere_calls = load_telomere_calls(args.telomere_dir)
    print(f"Loaded telomere calls for {len(telomere_calls):,} contigs", file=sys.stderr)

    # Step 5: Process PAF files
    print("\n[5/7] Processing PAF files...", file=sys.stderr)
    contig_data = process_paf_files(args.paf_pattern, arm_coords, chrom_alias, args.min_alignment_length)

    # Step 6: Classify contigs and apply filters
    print("\n[6/7] Classifying contigs and applying filters...", file=sys.stderr)
    classifications = classify_contigs(contig_data)
    print(f"Classified {len(classifications):,} contigs", file=sys.stderr)

    print("Applying telomere-based filtering...", file=sys.stderr)
    classifications = apply_telomeric_filter(
        classifications,
        contig_data,
        telomere_calls,
        args.min_contig_length
    )

    # Check and fail duplicate telomeres
    print("Checking for duplicate telomeres...", file=sys.stderr)
    failed_count = check_and_fail_duplicate_telomeres(classifications, args.output, t2t_contigs)

    # Add CHM13 reference contigs as valid T2T
    print("Adding CHM13 reference contigs...", file=sys.stderr)
    classifications = add_chm13_contigs(classifications, arm_coords)

    # Step 7: Write results and print statistics
    print("\n[7/7] Writing results...", file=sys.stderr)
    write_results(classifications, args.output)
    print_statistics(classifications)

    print("\nDone!", file=sys.stderr)


if __name__ == '__main__':
    main()
