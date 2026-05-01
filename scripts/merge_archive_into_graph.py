#!/usr/bin/env python3
"""Merge archive.jsonl back into graph.jsonl.

Reads each line of .wg/archive.jsonl, drops archive-only bookkeeping keys,
and appends the cleaned record (compact one-line JSON) to .wg/graph.jsonl.
"""
import json
import os
import sys
from pathlib import Path

WG = Path("/moosefs/erikg/phrs/.wg")
ARCHIVE = WG / "archive.jsonl"
GRAPH = WG / "graph.jsonl"

ARCHIVE_ONLY_KEYS = {
    "unplaced",
    "last_resurrected_at",
    "resurrection_count",
    "session_id",
    "triage_count",
}


def main() -> int:
    if not ARCHIVE.exists():
        print(f"ERROR: {ARCHIVE} does not exist", file=sys.stderr)
        return 1
    if not GRAPH.exists():
        print(f"ERROR: {GRAPH} does not exist", file=sys.stderr)
        return 1

    # Existing graph IDs (sanity check duplicates).
    existing_ids: set[str] = set()
    with GRAPH.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            obj = json.loads(line)
            existing_ids.add(obj["id"])

    archive_records = []
    with ARCHIVE.open() as f:
        for i, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError as e:
                print(f"ERROR: archive.jsonl line {i}: {e}", file=sys.stderr)
                return 2
            archive_records.append(obj)

    # Drop archive-only bookkeeping keys.
    cleaned = []
    skipped = []
    for obj in archive_records:
        for k in ARCHIVE_ONLY_KEYS:
            obj.pop(k, None)
        if obj["id"] in existing_ids:
            skipped.append(obj["id"])
            continue
        cleaned.append(obj)

    if skipped:
        print(f"WARNING: {len(skipped)} IDs already exist in graph.jsonl, skipping:")
        for sid in skipped:
            print(f"  - {sid}")

    # Append cleaned records (compact JSON, one per line, terminated by \n).
    appended = 0
    with GRAPH.open("a") as f:
        for obj in cleaned:
            f.write(json.dumps(obj, separators=(",", ":")))
            f.write("\n")
            appended += 1

    print(f"Appended {appended} records to {GRAPH}")
    print(f"Skipped {len(skipped)} duplicates")
    print(f"Archive total: {len(archive_records)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
