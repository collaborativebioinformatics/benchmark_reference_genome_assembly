#!/usr/bin/env python3

import argparse
import math
from collections import OrderedDict
from pathlib import Path

import pysam


def scalar(value):
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def event_key(record):
    return (
        record.contig,
        record.pos,
        record.stop,
        str(scalar(record.info.get("SVTYPE", "MISSING"))),
        str(scalar(record.info.get("SVLEN", "MISSING"))),
        record.ref,
        tuple(record.alts or ()),
    )


def genotype_signature(record):
    return tuple(
        (sample, tuple(record.samples[sample].get("GT") or ()), record.samples[sample].phased)
        for sample in record.samples
    )


def numeric(value, default):
    value = scalar(value)
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def score(record):
    return (
        numeric(record.qual, float("-inf")),
        numeric(record.info.get("SUPPORT"), float("-inf")),
    )


def display(value):
    return "." if value is None else str(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument(
        "--on-genotype-conflict",
        choices=("error", "exclude"),
        default="error",
    )
    args = parser.parse_args()

    source = pysam.VariantFile(str(args.input))
    groups = OrderedDict()
    for record in source:
        groups.setdefault(event_key(record), []).append(record.copy())

    conflicts = []
    output = pysam.VariantFile(str(args.output), "wz", header=source.header.copy())
    with open(args.audit, "w") as audit:
        audit.write(
            "CHROM\tPOS\tEND\tSVTYPE\tSVLEN\tREF\tALT\tretained_ID\tremoved_ID\t"
            "retained_QUAL\tremoved_QUAL\tretained_SUPPORT\tremoved_SUPPORT\tstatus\n"
        )
        for key, records in groups.items():
            signatures = {genotype_signature(record) for record in records}
            if len(signatures) > 1:
                if args.on_genotype_conflict == "error":
                    conflicts.append(key)
                for record in records:
                    audit.write(
                        "\t".join(
                            [
                                key[0], str(key[1]), str(key[2]), key[3], key[4], key[5],
                                ",".join(key[6]), ".", display(record.id), ".",
                                display(record.qual), ".",
                                display(scalar(record.info.get("SUPPORT"))),
                                "genotype_conflict" if args.on_genotype_conflict == "error"
                                else "excluded_ambiguous_liftover_collision",
                            ]
                        ) + "\n"
                    )
                continue
            retained = max(records, key=score)
            output.write(retained)
            for record in records:
                if record is retained:
                    continue
                audit.write(
                    "\t".join(
                        [
                            key[0], str(key[1]), str(key[2]), key[3], key[4], key[5],
                            ",".join(key[6]), display(retained.id), display(record.id),
                            display(retained.qual), display(record.qual),
                            display(scalar(retained.info.get("SUPPORT"))),
                            display(scalar(record.info.get("SUPPORT"))),
                            "removed_exact_sv_duplicate",
                        ]
                    ) + "\n"
                )
    output.close()
    source.close()

    if conflicts:
        args.output.unlink(missing_ok=True)
        raise SystemExit(
            f"Refusing to deduplicate {len(conflicts)} SV event groups with conflicting genotypes; "
            f"see {args.audit}"
        )


if __name__ == "__main__":
    main()
