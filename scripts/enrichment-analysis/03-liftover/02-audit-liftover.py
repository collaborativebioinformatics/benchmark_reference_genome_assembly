#!/usr/bin/env python3

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

import pysam


BND_PATTERN = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def scalar(value):
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def count_unmapped(path):
    reasons = Counter()
    if not path.exists():
        return 0, reasons
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line == "native_hs1":
                continue
            fields = line.split("\t")
            reason = fields[-1]
            if "not lifted" in reason:
                reason = "required_position_not_lifted"
            elif "different chromosome" in reason:
                reason = "different_chromosomes"
            elif "reverse order" in reason or " > lifted_END " in f" {reason} ":
                reason = "reversed_order"
            elif "changes significantly" in reason:
                reason = "length_change_gt_5pct"
            reasons[reason] += 1
    return sum(reasons.values()), reasons


def count_rejected(path):
    if not path.exists():
        return 0
    with open(path) as handle:
        return max(0, sum(1 for _ in handle) - 1)


def audit_vcf(vcf_path, fasta_path):
    fasta = pysam.FastaFile(str(fasta_path))
    contig_lengths = dict(zip(fasta.references, fasta.lengths))
    counts = Counter()
    svtypes = Counter()
    seen = set()
    previous = None

    with pysam.VariantFile(str(vcf_path)) as vcf:
        header_order = {name: index for index, name in enumerate(vcf.header.contigs)}
        for record in vcf:
            counts["accepted"] += 1
            svtype = str(scalar(record.info.get("SVTYPE", "MISSING")))
            svtypes[svtype] += 1
            key = (
                record.contig, record.pos, record.stop, record.ref,
                tuple(record.alts or ()), svtype,
            )
            if key in seen:
                counts["exact_duplicates"] += 1
            seen.add(key)

            if record.contig not in contig_lengths:
                counts["unknown_contig"] += 1
                continue
            if record.pos < 1 or record.pos > contig_lengths[record.contig]:
                counts["invalid_pos"] += 1
            if record.stop < record.pos or record.stop > contig_lengths[record.contig]:
                counts["invalid_end"] += 1

            expected_ref = fasta.fetch(
                record.contig, record.start, record.start + len(record.ref)
            ).upper()
            if record.ref.upper() != expected_ref:
                counts["ref_mismatch"] += 1

            current = (header_order.get(record.contig, 10**9), record.pos)
            if previous is not None and current < previous:
                counts["sort_errors"] += 1
            previous = current

            svlen = scalar(record.info.get("SVLEN"))
            try:
                svlen = abs(int(svlen))
            except (TypeError, ValueError):
                svlen = None
            if svtype == "INS" and record.stop != record.pos:
                counts["insertion_end_not_pos"] += 1
            if svtype in {"DEL", "DUP", "INV"}:
                span = record.stop - record.pos
                if span < 1:
                    counts["nonpositive_span"] += 1
                if svlen is None:
                    counts["missing_svlen"] += 1
                elif span != svlen:
                    counts["span_svlen_disagreement"] += 1
            if svtype == "BND":
                alt = (record.alts or ("",))[0]
                mate = BND_PATTERN.search(alt)
                if mate is None:
                    counts["invalid_bnd_alt"] += 1
                else:
                    mate_contig, mate_pos_text = mate.groups()
                    mate_pos = int(mate_pos_text)
                    if mate_contig not in contig_lengths:
                        counts["bnd_unknown_mate_contig"] += 1
                    elif mate_pos < 1 or mate_pos > contig_lengths[mate_contig]:
                        counts["bnd_invalid_mate_pos"] += 1

    fasta.close()
    return counts, svtypes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--callsets", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    fasta = (args.root / "../references/hs1/hs1.fa").resolve()
    results = []
    with open(args.callsets, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            arm, build = row["arm"], row["build"]
            prefix = args.root / "results/lifted" / arm / build / f"{arm}.{build}.to-hs1"
            vcf = Path(f"{prefix}.vcf.gz")
            if not vcf.exists():
                continue
            counts, svtypes = audit_vcf(vcf, fasta)
            liftover_unmapped, reasons = count_unmapped(Path(f"{prefix}.unmapped.tsv"))
            canonical_rejected = count_rejected(Path(f"{prefix}.canonicalization-rejected.tsv"))
            excluded = liftover_unmapped + canonical_rejected
            source_count_path = (
                args.root / "results/preprocessed" / arm / build / f"{arm}.{build}.pass.vcf.gz.records.txt"
            )
            source_count = int(source_count_path.read_text().strip())
            accounting_delta = source_count - counts["accepted"] - excluded
            results.append({
                "task_id": int(row["task_id"]),
                "arm": arm,
                "build": build,
                "source_pass": source_count,
                "accepted": counts["accepted"],
                "excluded": excluded,
                "liftover_unmapped": liftover_unmapped,
                "canonical_rejected": canonical_rejected,
                "accounting_delta": accounting_delta,
                "acceptance_fraction": counts["accepted"] / source_count if source_count else 0,
                "counts": dict(counts),
                "svtypes": dict(svtypes),
                "unmapped_reasons": dict(reasons),
                "vcf": str(vcf),
            })

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_prefix.with_suffix(".json"), "w") as handle:
        json.dump(results, handle, indent=2, sort_keys=True)
        handle.write("\n")

    count_fields = [
        "unknown_contig", "invalid_pos", "invalid_end", "ref_mismatch", "sort_errors",
        "exact_duplicates", "insertion_end_not_pos", "nonpositive_span", "missing_svlen",
        "span_svlen_disagreement", "invalid_bnd_alt", "bnd_unknown_mate_contig",
        "bnd_invalid_mate_pos",
    ]
    fields = [
        "task_id", "arm", "build", "source_pass", "accepted", "excluded",
        "liftover_unmapped", "canonical_rejected",
        "accounting_delta", "acceptance_fraction", *count_fields, "svtypes",
        "unmapped_reasons", "vcf",
    ]
    with open(args.output_prefix.with_suffix(".tsv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for result in results:
            writer.writerow({
                **{key: result[key] for key in fields[:10]},
                **{key: result["counts"].get(key, 0) for key in count_fields},
                "svtypes": json.dumps(result["svtypes"], sort_keys=True),
                "unmapped_reasons": json.dumps(result["unmapped_reasons"], sort_keys=True),
                "vcf": result["vcf"],
            })

    critical_fields = [
        "unknown_contig", "invalid_pos", "invalid_end", "ref_mismatch", "sort_errors",
        "insertion_end_not_pos", "nonpositive_span", "missing_svlen",
        "span_svlen_disagreement", "invalid_bnd_alt", "bnd_unknown_mate_contig",
        "bnd_invalid_mate_pos",
    ]
    failures = [
        f"{row['arm']}/{row['build']}:{field}={row['counts'].get(field, 0)}"
        for row in results
        for field in critical_fields
        if row["counts"].get(field, 0)
    ]
    failures.extend(
        f"{row['arm']}/{row['build']}:accounting_delta={row['accounting_delta']}"
        for row in results
        if row["accounting_delta"] != 0
    )
    if failures:
        raise SystemExit("Liftover audit failed: " + "; ".join(failures))


if __name__ == "__main__":
    main()
