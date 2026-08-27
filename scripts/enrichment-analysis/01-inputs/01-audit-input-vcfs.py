#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_info(value):
    result = {}
    if value == ".":
        return result
    for field in value.split(";"):
        if "=" in field:
            key, item = field.split("=", 1)
            result[key] = item
        else:
            result[field] = True
    return result


def integer(value):
    if value in (None, "."):
        return None
    try:
        return int(str(value).split(",", 1)[0])
    except ValueError:
        return None


def audit(path):
    contigs = {}
    samples = []
    sources = []
    reference = ""
    counts = Counter()
    svtypes = Counter()
    filters = Counter()
    genotypes = Counter()
    duplicate_keys = set()
    seen_keys = set()
    last = None

    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("##contig=<"):
                body = raw.rstrip()[10:-1]
                fields = dict(
                    item.split("=", 1) for item in body.split(",") if "=" in item
                )
                if "ID" in fields:
                    contigs[fields["ID"]] = integer(fields.get("length"))
                continue
            if raw.startswith("##source="):
                sources.append(raw.rstrip().split("=", 1)[1])
                continue
            if raw.startswith("##reference="):
                reference = raw.rstrip().split("=", 1)[1]
                continue
            if raw.startswith("#CHROM"):
                columns = raw.rstrip().split("\t")
                samples = columns[9:]
                continue
            if raw.startswith("#"):
                continue

            fields = raw.rstrip("\n").split("\t")
            counts["records"] += 1
            if len(fields) < 8:
                counts["malformed"] += 1
                continue
            chrom, pos_text, ident, ref, alt, _, filt, info_text = fields[:8]
            pos = integer(pos_text)
            info = parse_info(info_text)
            svtype = str(info.get("SVTYPE", "MISSING")).split(",", 1)[0]
            end = integer(info.get("END"))
            svlen = integer(info.get("SVLEN"))
            svtypes[svtype] += 1
            filters[filt] += 1

            if ident in ("", "."):
                counts["missing_id"] += 1
            if "SVTYPE" not in info:
                counts["missing_svtype"] += 1
            if "END" not in info:
                counts["missing_end"] += 1
            if "SVLEN" not in info:
                counts["missing_svlen"] += 1
            if alt.startswith("<") and alt.endswith(">"):
                counts["symbolic"] += 1
            elif "[" in alt or "]" in alt:
                counts["breakend"] += 1
            else:
                counts["sequence_resolved"] += 1
            if chrom not in contigs:
                counts["unknown_contig"] += 1
            if pos is None or pos < 1:
                counts["invalid_pos"] += 1
            elif contigs.get(chrom) is not None and pos > contigs[chrom]:
                counts["pos_past_contig"] += 1
            if end is not None:
                if pos is not None and end < pos:
                    counts["end_before_pos"] += 1
                if contigs.get(chrom) is not None and end > contigs[chrom]:
                    counts["end_past_contig"] += 1
            if svlen == 0:
                counts["zero_svlen"] += 1

            key = (chrom, pos, ref, alt, end, svtype)
            if key in seen_keys:
                duplicate_keys.add(key)
            else:
                seen_keys.add(key)

            if last is not None:
                last_chrom, last_pos = last
                if chrom == last_chrom and pos is not None and last_pos is not None and pos < last_pos:
                    counts["within_contig_sort_errors"] += 1
            last = (chrom, pos)

            if len(fields) >= 10:
                formats = fields[8].split(":")
                if "GT" in formats:
                    gt_index = formats.index("GT")
                    for sample_field in fields[9:]:
                        sample_values = sample_field.split(":")
                        gt = sample_values[gt_index] if gt_index < len(sample_values) else "."
                        genotypes[gt] += 1
                else:
                    counts["records_missing_gt_format"] += 1
            else:
                counts["records_without_samples"] += 1

    counts["exact_duplicate_records"] = len(duplicate_keys)
    return {
        "path": str(path),
        "source": ",".join(sources),
        "reference": reference,
        "samples": samples,
        "contig_count": len(contigs),
        "counts": dict(counts),
        "svtypes": dict(sorted(svtypes.items())),
        "filters": dict(sorted(filters.items())),
        "genotypes": dict(sorted(genotypes.items())),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    with open(args.manifest, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["role"] not in {"callset", "truth_vcf"}:
                continue
            path = args.root / row["local_path"]
            result = audit(path)
            result.update({key: row[key] for key in ("role", "arm", "build", "dx_file_id")})
            rows.append(result)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_prefix.with_suffix(".json"), "w") as handle:
        json.dump(rows, handle, indent=2, sort_keys=True)
        handle.write("\n")

    fields = [
        "role", "arm", "build", "dx_file_id", "source", "reference", "samples",
        "contig_count", "records", "PASS", "symbolic", "sequence_resolved", "breakend",
        "missing_svtype", "missing_end", "missing_svlen", "records_missing_gt_format",
        "unknown_contig", "invalid_pos", "end_before_pos", "pos_past_contig",
        "end_past_contig", "exact_duplicate_records", "within_contig_sort_errors",
        "svtypes", "filters", "genotypes", "path",
    ]
    with open(args.output_prefix.with_suffix(".tsv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for result in rows:
            counts = result["counts"]
            filters = result["filters"]
            writer.writerow({
                "role": result["role"],
                "arm": result["arm"],
                "build": result["build"],
                "dx_file_id": result["dx_file_id"],
                "source": result["source"],
                "reference": result["reference"],
                "samples": ",".join(result["samples"]),
                "contig_count": result["contig_count"],
                "records": counts.get("records", 0),
                "PASS": filters.get("PASS", 0),
                "symbolic": counts.get("symbolic", 0),
                "sequence_resolved": counts.get("sequence_resolved", 0),
                "breakend": counts.get("breakend", 0),
                "missing_svtype": counts.get("missing_svtype", 0),
                "missing_end": counts.get("missing_end", 0),
                "missing_svlen": counts.get("missing_svlen", 0),
                "records_missing_gt_format": counts.get("records_missing_gt_format", 0),
                "unknown_contig": counts.get("unknown_contig", 0),
                "invalid_pos": counts.get("invalid_pos", 0),
                "end_before_pos": counts.get("end_before_pos", 0),
                "pos_past_contig": counts.get("pos_past_contig", 0),
                "end_past_contig": counts.get("end_past_contig", 0),
                "exact_duplicate_records": counts.get("exact_duplicate_records", 0),
                "within_contig_sort_errors": counts.get("within_contig_sort_errors", 0),
                "svtypes": json.dumps(result["svtypes"], sort_keys=True),
                "filters": json.dumps(filters, sort_keys=True),
                "genotypes": json.dumps(result["genotypes"], sort_keys=True),
                "path": result["path"],
            })


if __name__ == "__main__":
    main()
