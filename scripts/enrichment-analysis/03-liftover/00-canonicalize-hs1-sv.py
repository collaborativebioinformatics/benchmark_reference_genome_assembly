#!/usr/bin/env python3

import argparse
from pathlib import Path

import pysam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--rejected", type=Path, required=True)
    args = parser.parse_args()

    fasta = pysam.FastaFile(str(args.reference))
    source = pysam.VariantFile(str(args.input))
    header = source.header.copy()
    header.add_line(
        "##hs1_canonicalization=DEL,DUP,INV converted to symbolic alleles; "
        "REF set from hs1; span SVLEN recomputed"
    )
    output = pysam.VariantFile(str(args.output), "wz", header=header)
    rejected = open(args.rejected, "w")
    rejected.write("CHROM\tPOS\tID\tSVTYPE\treason\n")

    for record in source:
        svtype_value = record.info.get("SVTYPE")
        if isinstance(svtype_value, tuple):
            svtype_value = svtype_value[0] if svtype_value else None
        svtype = str(svtype_value) if svtype_value is not None else ""
        if record.contig not in fasta.references:
            rejected.write(f"{record.contig}\t{record.pos}\t{record.id}\t{svtype}\tunknown_contig\n")
            continue
        contig_length = fasta.get_reference_length(record.contig)
        if record.pos < 1 or record.pos > contig_length or record.stop < record.pos or record.stop > contig_length:
            rejected.write(f"{record.contig}\t{record.pos}\t{record.id}\t{svtype}\tinvalid_span\n")
            continue
        if svtype in {"DEL", "DUP", "INV"} and record.stop == record.pos:
            rejected.write(f"{record.contig}\t{record.pos}\t{record.id}\t{svtype}\tzero_length_span\n")
            continue
        if svtype in {"DEL", "DUP", "INV"}:
            old_stop = record.stop
            anchor = fasta.fetch(record.contig, record.start, record.start + 1).upper()
            record.ref = anchor
            record.alts = (f"<{svtype}>",)
            record.stop = old_stop
            span = old_stop - record.pos
            record.info["SVLEN"] = -span if svtype == "DEL" else span
        elif svtype == "INS":
            old_ref = record.ref
            old_stop = record.stop
            anchor = fasta.fetch(record.contig, record.start, record.start + 1).upper()
            old_alt = (record.alts or ("",))[0]
            if old_alt.startswith("<") and old_alt.endswith(">"):
                record.ref = anchor
            elif old_alt.upper().startswith(old_ref.upper()):
                record.ref = anchor
                record.alts = (anchor + old_alt[len(old_ref):],)
            else:
                rejected.write(
                    f"{record.contig}\t{record.pos}\t{record.id}\t{svtype}\t"
                    "insertion_alt_does_not_start_with_ref\n"
                )
                continue
            record.stop = old_stop
        output.write(record)

    output.close()
    source.close()
    fasta.close()
    rejected.close()


if __name__ == "__main__":
    main()
