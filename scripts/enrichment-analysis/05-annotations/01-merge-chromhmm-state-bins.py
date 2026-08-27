#!/usr/bin/env python3

import argparse
from collections import Counter, defaultdict
from pathlib import Path


def emit(handle, pending, chrom, start, end, label):
    if start >= end:
        return pending
    if pending is not None and pending[:2] == (chrom, label) and pending[3] == start:
        return chrom, label, pending[2], end
    if pending is not None:
        handle.write(f"{pending[0]}\t{pending[2]}\t{pending[3]}\t{pending[1]}\n")
    return chrom, label, start, end


def resolve(chrom, intervals, output, ambiguous, output_pending, ambiguous_pending):
    events = defaultdict(list)
    for start, end, state in intervals:
        events[start].append((state, 1))
        events[end].append((state, -1))
    active = Counter()
    previous = None
    for position in sorted(events):
        if previous is not None and previous < position:
            states = sorted(state for state, count in active.items() if count > 0)
            if len(states) == 1:
                output_pending = emit(output, output_pending, chrom, previous, position, states[0])
            elif len(states) > 1:
                ambiguous_pending = emit(
                    ambiguous, ambiguous_pending, chrom, previous, position, ",".join(states)
                )
        for state, delta in events[position]:
            active[state] += delta
            if active[state] == 0:
                del active[state]
        previous = position
    return output_pending, ambiguous_pending


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ambiguous", type=Path, required=True)
    args = parser.parse_args()

    current_chrom = None
    intervals = []
    output_pending = ambiguous_pending = None
    with open(args.input) as source, open(args.output, "w") as output, open(args.ambiguous, "w") as ambiguous:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start_text, end_text, state, *_ = line.rstrip("\n").split("\t")
            start, end = int(start_text), int(end_text)
            if current_chrom is not None and chrom != current_chrom:
                output_pending, ambiguous_pending = resolve(
                    current_chrom, intervals, output, ambiguous, output_pending, ambiguous_pending
                )
                intervals = []
            current_chrom = chrom
            intervals.append((start, end, state))
        if current_chrom is not None:
            output_pending, ambiguous_pending = resolve(
                current_chrom, intervals, output, ambiguous, output_pending, ambiguous_pending
            )
        if output_pending is not None:
            output.write(
                f"{output_pending[0]}\t{output_pending[2]}\t{output_pending[3]}\t{output_pending[1]}\n"
            )
        if ambiguous_pending is not None:
            ambiguous.write(
                f"{ambiguous_pending[0]}\t{ambiguous_pending[2]}\t"
                f"{ambiguous_pending[3]}\t{ambiguous_pending[1]}\n"
            )


if __name__ == "__main__":
    main()
