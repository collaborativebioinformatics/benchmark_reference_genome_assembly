#!/usr/bin/env python3

import argparse
from collections import defaultdict
from pathlib import Path


def attributes(text):
    values = {}
    for item in text.rstrip("\n").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            values[key] = value
    return values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--repeatmasker", type=Path, required=True)
    parser.add_argument("--segdup", type=Path, required=True)
    parser.add_argument("--fai", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    contig_lengths = {}
    with open(args.fai) as handle:
        for line in handle:
            chrom, length, *_ = line.rstrip("\n").split("\t")
            contig_lengths[chrom] = int(length)

    exons = defaultdict(list)
    cds_bounds = {}
    with (
        open(args.gff) as source,
        open(args.output_dir / "genes.raw.bed", "w") as genes,
        open(args.output_dir / "features.raw.bed", "w") as features,
    ):
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, _, feature, start_text, end_text, _, strand, _, attr_text = line.split("\t", 8)
            if chrom not in contig_lengths:
                continue
            start, end = int(start_text) - 1, int(end_text)
            attrs = attributes(attr_text)
            gene_name = attrs.get("gene_name", attrs.get("Name", "."))
            gene_id = attrs.get("gene_id", attrs.get("source_gene", "."))
            transcript_id = attrs.get("transcript_id", attrs.get("Parent", "."))
            if feature == "gene":
                biotype = attrs.get("gene_biotype", ".")
                genes.write(
                    f"{chrom}\t{start}\t{end}\t{gene_name}\t{gene_id}\t{strand}\t{biotype}\n"
                )
                tss = start if strand == "+" else end
                if strand == "+":
                    promoter_start, promoter_end = max(0, tss - 2000), min(contig_lengths[chrom], tss + 200)
                else:
                    promoter_start, promoter_end = max(0, tss - 200), min(contig_lengths[chrom], tss + 2000)
                features.write(
                    f"{chrom}\t{promoter_start}\t{promoter_end}\tpromoter\t{gene_name}\t{gene_id}\n"
                )
            elif feature in {"exon", "intron"}:
                category = feature
                features.write(
                    f"{chrom}\t{start}\t{end}\t{category}\t{gene_name}\t{gene_id}\n"
                )
                if feature == "exon":
                    exons[transcript_id].append((chrom, start, end, gene_name, gene_id))
            elif feature == "CDS":
                if transcript_id not in cds_bounds:
                    cds_bounds[transcript_id] = [start, end]
                else:
                    cds_bounds[transcript_id][0] = min(cds_bounds[transcript_id][0], start)
                    cds_bounds[transcript_id][1] = max(cds_bounds[transcript_id][1], end)

        for transcript_id, transcript_exons in exons.items():
            if transcript_id not in cds_bounds:
                continue
            cds_start, cds_end = cds_bounds[transcript_id]
            for chrom, start, end, gene_name, gene_id in transcript_exons:
                if start < cds_start:
                    utr_end = min(end, cds_start)
                    if start < utr_end:
                        features.write(
                            f"{chrom}\t{start}\t{utr_end}\tUTR\t{gene_name}\t{gene_id}\n"
                        )
                if end > cds_end:
                    utr_start = max(start, cds_end)
                    if utr_start < end:
                        features.write(
                            f"{chrom}\t{utr_start}\t{end}\tUTR\t{gene_name}\t{gene_id}\n"
                        )

    with open(args.repeatmasker) as source, open(args.output_dir / "repeatmasker.raw.bed", "w") as output:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or fields[0] not in contig_lengths:
                continue
            output.write(
                f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[6]}\t{fields[7]}\t{fields[3]}\n"
            )

    with open(args.segdup) as source, open(args.output_dir / "segdup.raw.bed", "w") as output:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] not in contig_lengths:
                continue
            output.write("\t".join(fields[:4]) + "\n")


if __name__ == "__main__":
    main()
