#!/usr/bin/env python3
import os
import re
import subprocess
import sys

from common import COMPARISONS, PACKAGE_DIR, read_tsv, write_tsv


HAP_LABEL = {
    "#1": "maternal",
    "#2": "paternal",
}


def read_fai_names(fasta):
    fai = fasta + ".fai"
    names = []
    with open(fai) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0]:
                names.append(fields[0])
    return names


def stream_fasta_names(fasta):
    names = []
    with subprocess.Popen(["gzip", "-cd", fasta], stdout=subprocess.PIPE, text=True) as proc:
        for line in proc.stdout:
            if line.startswith(">"):
                names.append(line[1:].strip().split()[0])
        ret = proc.wait()
        if ret != 0:
            raise subprocess.CalledProcessError(ret, ["gzip", "-cd", fasta])
    return names


def names_for_source(fasta, cache):
    if fasta in cache:
        return cache[fasta]
    try:
        names = read_fai_names(fasta)
        source = fasta + ".fai"
    except OSError:
        names = stream_fasta_names(fasta)
        source = "streamed full FASTA headers after .fai read failed"
    cache[fasta] = (names, source)
    return cache[fasta]


def select_haplotype_names(fasta, prefix, cache):
    names, source = names_for_source(fasta, cache)
    selected = [name for name in names if name.startswith(prefix + "#chr")]
    if not selected:
        raise SystemExit("No full-chromosome records matched %s in %s.fai" % (prefix, fasta))
    return selected, source


def write_lines(path, lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as out:
        for line in lines:
            out.write(line + "\n")


def faidx_extract(fasta, names_path, fasta_out):
    with open(fasta_out, "w") as out:
        subprocess.check_call(["samtools", "faidx", "-r", names_path, fasta], stdout=out)
    subprocess.check_call(["samtools", "faidx", fasta_out])


def stream_extract(fasta, wanted_names, fasta_out):
    wanted = set(wanted_names)
    keep = False
    with subprocess.Popen(["gzip", "-cd", fasta], stdout=subprocess.PIPE, text=True) as proc, open(fasta_out, "w") as out:
        for line in proc.stdout:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                keep = name in wanted
            if keep:
                out.write(line)
        ret = proc.wait()
        if ret != 0:
            raise subprocess.CalledProcessError(ret, ["gzip", "-cd", fasta])
    subprocess.check_call(["samtools", "faidx", fasta_out])


def extract_records(fasta, names_path, names, fasta_out, name_source):
    if name_source.endswith(".fai"):
        try:
            faidx_extract(fasta, names_path, fasta_out)
            return "samtools faidx using source .fai"
        except (OSError, subprocess.CalledProcessError):
            pass
    stream_extract(fasta, names, fasta_out)
    return "streamed full source FASTA with gzip -cd; selected records by exact header"


def collapse_target_group_for_joint_parent(fasta_path, target_prefixes, target_sample):
    joint_prefix = target_sample + "#joint"
    tmp_path = fasta_path + ".tmp"
    prefix_re = re.compile(r"^(%s)#(.+)$" % "|".join(re.escape(p) for p in target_prefixes))
    with open(fasta_path) as inp, open(tmp_path, "w") as out:
        for line in inp:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                match = prefix_re.match(header)
                if match:
                    original_prefix = match.group(1)
                    hap_token = original_prefix.rsplit("#", 1)[-1]
                    header = "%s#h%s_%s" % (joint_prefix, hap_token, match.group(2))
                out.write(">" + header + "\n")
            else:
                out.write(line)
    os.replace(tmp_path, fasta_path)
    subprocess.check_call(["samtools", "faidx", fasta_path])


def hap_label(prefix):
    for token, label in HAP_LABEL.items():
        if token in prefix:
            return label
    return "unknown"


def main():
    rows = []
    name_cache = {}
    for row in read_tsv(COMPARISONS):
        cid = row["comparison_id"]
        query_names, query_name_source = select_haplotype_names(row["query_source_fasta"], row["query_haplotype"], name_cache)
        target_prefixes = [x for x in row["target_haplotypes"].split("+") if x]
        target_names = []
        target_counts = []
        target_name_sources = []
        for prefix in target_prefixes:
            names, target_name_source = select_haplotype_names(row["target_source_fasta"], prefix, name_cache)
            target_names.extend(names)
            target_counts.append("%s:%d" % (prefix, len(names)))
            target_name_sources.append("%s:%s" % (prefix, target_name_source))

        query_names_path = os.path.join(PACKAGE_DIR, "inputs", cid + ".query.names.txt")
        target_names_path = os.path.join(PACKAGE_DIR, "inputs", cid + ".target.names.txt")
        query_fa = os.path.join(PACKAGE_DIR, "inputs", cid + ".query.fa")
        target_fa = os.path.join(PACKAGE_DIR, "inputs", cid + ".target.fa")
        write_lines(query_names_path, query_names)
        write_lines(target_names_path, target_names)

        query_extraction = "existing nonempty FASTA reused"
        target_extraction = "existing nonempty FASTA reused"
        if not os.path.exists(query_fa) or os.path.getsize(query_fa) == 0:
            query_extraction = extract_records(row["query_source_fasta"], query_names_path, query_names, query_fa, query_name_source)
        if not os.path.exists(target_fa) or os.path.getsize(target_fa) == 0:
            target_extraction = extract_records(row["target_source_fasta"], target_names_path, target_names, target_fa, target_name_sources[0].split(":", 1)[1])
            collapse_target_group_for_joint_parent(target_fa, target_prefixes, row["target_sample"])

        rows.append({
            "comparison_id": cid,
            "query_source_fasta": row["query_source_fasta"],
            "query_haplotype_prefix": row["query_haplotype"],
            "query_haplotype_label": hap_label(row["query_haplotype"]),
            "query_sequence_count": str(len(query_names)),
            "query_first_sequence": query_names[0],
            "query_last_sequence": query_names[-1],
            "query_fasta": query_fa,
            "query_name_source": query_name_source,
            "query_extraction": query_extraction,
            "target_source_fasta": row["target_source_fasta"],
            "target_haplotype_prefixes": row["target_haplotypes"],
            "target_haplotype_counts": ";".join(target_counts),
            "target_name_sources": ";".join(target_name_sources),
            "target_extraction": target_extraction,
            "target_joint_header_prefix": row["target_sample"] + "#joint",
            "target_first_sequence": target_names[0],
            "target_last_sequence": target_names[-1],
            "target_fasta": target_fa,
            "scope": "full whole-genome haplotype FASTA records from source assembly .fai; no telomeric-window FASTA",
        })

    write_tsv(
        os.path.join(PACKAGE_DIR, "summaries", "input_manifest.tsv"),
        rows,
        [
            "comparison_id",
            "query_source_fasta",
            "query_haplotype_prefix",
            "query_haplotype_label",
            "query_sequence_count",
            "query_first_sequence",
            "query_last_sequence",
            "query_fasta",
            "query_name_source",
            "query_extraction",
            "target_source_fasta",
            "target_haplotype_prefixes",
            "target_haplotype_counts",
            "target_name_sources",
            "target_extraction",
            "target_joint_header_prefix",
            "target_first_sequence",
            "target_last_sequence",
            "target_fasta",
            "scope",
        ],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
