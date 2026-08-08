# OPTIMuS-PRIME-Del

> **Note:** This repository is under active construction and will be updated with additional code, data, and documentation as the manuscripts are finalized. Scripts provided here are shared for transparency and do not constitute streamlined pipelines.

Code and data for: **Multiplex programming of mappable, molecularly tagged genomic deletions**

## Repository layout

```
OPTIMuS-PRIME-Del/
├── T7_mapping/              # Bulk T7 IVT mapping (v1 18K and v2 top-100 libraries)
│   ├── v1_18K_K562/         # Alignment and deletion count parsing for the 18K library in K562
│   └── v2_top100_K562/      # Alignment and barcode/deletion count parsing for the top-100 library in K562
├── amplicon_seq/            # Amplicon-seq validation of individual programmed deletions (top-100, K562)
├── scT7_mapping/            # Single-cell RNA-seq + T7-capture pipeline for single-cell deletion genotyping (top-100, K562)
├── figures/                 # Notebooks that generate the manuscript's main and supplemental figure panels
└── plasmid_library/         # TODO: plasmid library reference/QC files (v1 and v2)
```

## What's in each directory

- **`T7_mapping/`** – Alignment and pysam-based parsing of T7 IVT reads to map and count deletions, for both deletion libraries/versions.
- **`amplicon_seq/`** – Amplicon sequencing of individually validated deletions: alignment to per-design deletion/WT reference sequences, read/UMI counting per allele, and per-position variant (indel) frequency profiling used to assess editing precision.
- **`scT7_mapping/`** – Cell Ranger processing of the single-cell RNA transcriptome and T7-capture libraries, barcode/UMI-corrected T7 capture-to-cell-barcode association calling, removal of residual T7-IST reads from the transcriptome library, a second Cell Ranger pass on the filtered data, Seurat-based QC/normalization, and empirical FDR-based differential expression testing of genes near each targeted cCRE.
- **`figures/`** – Notebooks that take the processed outputs above and produce the manuscript's figure panels.

## Notes

Filenames referenced in scripts/notebooks use internal cluster paths (see comments within individual scripts for the specific internal-to-GEO naming mapping where naming conventions differ).

## Software

Analyses were run using (versions confirmed at the time of analysis):

- bwa 0.7.17-r1188
- samtools 1.19
- cutadapt 4.9
- bedtools 2.31.1
- Cell Ranger 8.0.0
- R 4.3.2 (Seurat, GenomicRanges, Signac, rtracklayer, Biostrings, ggplot2, and related packages)
- Python 3 (pysam, Biopython, pandas, numpy)
