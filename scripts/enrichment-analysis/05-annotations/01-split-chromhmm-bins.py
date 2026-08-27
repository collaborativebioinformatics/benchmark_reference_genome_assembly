#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bin-size", type=int, default=200)
    args = parser.parse_args()

    identifier = 0
    with gzip.open(args.input, "rt") as source, open(args.output, "w") as output:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start_text, end_text, state, *_ = line.rstrip("\n").split("\t")
            start, end = int(start_text), int(end_text)
            for chunk_start in range(start, end, args.bin_size):
                chunk_end = min(chunk_start + args.bin_size, end)
                identifier += 1
                output.write(
                    f"{chrom}\t{chunk_start}\t{chunk_end}\t{state}\t"
                    f"chromhmm_{identifier}\t{chunk_end - chunk_start}\n"
                )


if __name__ == "__main__":
    main()

