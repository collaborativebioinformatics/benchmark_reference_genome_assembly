# Supplementary Table S10. Round-trip coordinate-projection control

A native hs1 callset was projected from hs1 through each legacy reference and returned to hs1 before comparison with unchanged hs1 truth. Outward survival is the proportion of eligible native calls retained after the first hs1-to-legacy projection; return survival is the proportion of outward survivors retained after the legacy-to-hs1 projection.

The verified projection implementation retained PASS INS, DEL, and DUP calls with representable spans greater than 0 bp and at most 1,000,000 bp, excluded BND and INV calls, and applied a 10% maximum span-change filter to DEL and DUP records. The historical implementation records aggregate input and retained counts only; it does not distinguish unmapped, fragmented/ambiguous, and span-filtered outcomes.

| Aligner | Route | Eligible input calls | Calls after outward lift | Outward survival | Return survival | Native hs1 F1 | Round-trip F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| Minimap2 | hs1 → hg38 → hs1 | 22,597 | 11,269 | 49.87% | 98.19% | 0.770 | 0.396 |
| Minimap2 | hs1 → hg19 → hs1 | 22,652 | 10,932 | 48.26% | 98.43% | 0.770 | 0.394 |
| Minimap2 | hs1 → hg18 → hs1 | 22,685 | 10,830 | 47.74% | 98.57% | 0.770 | 0.393 |
| Minimap2 | hs1 → hg17 → hs1 | 22,704 | 10,807 | 47.60% | 98.66% | 0.770 | 0.393 |
| Minimap2 | hs1 → hg16 → hs1 | 22,678 | 10,804 | 47.64% | 98.53% | 0.770 | 0.391 |
| Minimap2 | hs1 → hg15 → hs1 | 22,577 | 10,769 | 47.70% | 98.10% | 0.770 | 0.390 |
| Winnowmap | hs1 → hg38 → hs1 | 23,841 | 11,339 | 47.56% | 96.72% | 0.734 | 0.381 |
| Winnowmap | hs1 → hg19 → hs1 | 24,082 | 10,813 | 44.90% | 97.69% | 0.734 | 0.379 |
| Winnowmap | hs1 → hg18 → hs1 | 24,128 | 10,684 | 44.28% | 97.89% | 0.734 | 0.378 |
| Winnowmap | hs1 → hg17 → hs1 | 24,160 | 10,652 | 44.09% | 98.01% | 0.734 | 0.378 |
| Winnowmap | hs1 → hg16 → hs1 | 24,093 | 10,719 | 44.49% | 97.74% | 0.734 | 0.377 |
| Winnowmap | hs1 → hg15 → hs1 | 23,865 | 10,725 | 44.94% | 96.81% | 0.734 | 0.375 |

Coordinate projection alone removed approximately half of eligible calls on the outward transfer, whereas most calls surviving that first step successfully returned to hs1. Projected benchmarking metrics should therefore be interpreted as applying to a selected survivor subset rather than to the original native callset.
