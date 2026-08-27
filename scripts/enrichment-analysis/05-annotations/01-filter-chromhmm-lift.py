#!/usr/bin/env python3

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--accepted", type=Path, required=True)
    parser.add_argument("--rejected", type=Path, required=True)
    parser.add_argument("--minimum-ratio", type=float, default=0.80)
    parser.add_argument("--maximum-ratio", type=float, default=1.20)
    args = parser.parse_args()

    with open(args.input) as source, open(args.accepted, "w") as accepted, open(args.rejected, "w") as rejected:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                rejected.write(line.rstrip("\n") + "\tmalformed_lifted_bed\n")
                continue
            target_length = int(fields[2]) - int(fields[1])
            source_length = int(fields[5])
            ratio = target_length / source_length if source_length else 0
            if args.minimum_ratio <= ratio <= args.maximum_ratio:
                accepted.write(line)
            else:
                rejected.write(line.rstrip("\n") + f"\tlength_ratio={ratio:.6f}\n")


if __name__ == "__main__":
    main()

