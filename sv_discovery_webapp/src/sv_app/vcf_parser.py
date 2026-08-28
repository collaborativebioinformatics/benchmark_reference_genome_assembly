"""
Lightweight VCF parser for structural variant files produced by Sniffles2.
Avoids external VCF libraries — just regex over the INFO column, which is
all we need for the fields this app displays.
"""

import re
import pandas as pd

# Regex to pull key=value pairs out of the INFO column (e.g. "SVTYPE=DEL;SVLEN=-401;...")
_INFO_KV = re.compile(r"([A-Za-z_]+)=([^;]+)")


def _parse_info(info_str):
    """Turn a VCF INFO string into a dict of its key=value pairs."""
    return dict(_INFO_KV.findall(info_str))


def parse_vcf(vcf_path):
    """
    Read a Sniffles2 VCF file and return a tidy DataFrame with one row per
    structural variant, containing the fields most useful for display:
    CHROM, POS, SVTYPE, SVLEN, END, SUPPORT, VAF, FILTER, QUAL, GENOTYPE.
    """
    rows = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue  # skip header lines
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, _id, _ref, alt, qual, filt, info = fields[:8]
            info_dict = _parse_info(info)

            genotype = None
            if len(fields) >= 10:
                # FORMAT column tells us the order of subfields in SAMPLE column
                format_keys = fields[8].split(":")
                sample_vals = fields[9].split(":")
                sample_dict = dict(zip(format_keys, sample_vals))
                genotype = sample_dict.get("GT")

            rows.append({
                "CHROM": chrom,
                "POS": int(pos),
                "SVTYPE": info_dict.get("SVTYPE", "UNKNOWN"),
                "SVLEN": int(info_dict["SVLEN"]) if "SVLEN" in info_dict else None,
                "END": int(info_dict["END"]) if "END" in info_dict else None,
                "SUPPORT": int(info_dict["SUPPORT"]) if "SUPPORT" in info_dict else None,
                "VAF": float(info_dict["VAF"]) if "VAF" in info_dict else None,
                "QUAL": float(qual) if qual not in (".", "") else None,
                "FILTER": filt,
                "GENOTYPE": genotype,
            })

    return pd.DataFrame(rows)
