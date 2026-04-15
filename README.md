# Non-Cancer-Related Pathogenic Germline Variants and Expression Consequences in Ten-Thousand Cancer Genomes

[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs13073--021--00964--1-blue)](https://doi.org/10.1186/s13073-021-00964-1)
[![Language](https://img.shields.io/badge/Language-R-276DC3)](https://www.r-project.org/)

## Overview

This repository contains the analysis code accompanying our publication in *Genome Medicine*:

> Wang Z, *et al.* Non-cancer-related pathogenic germline variants and expression consequences in ten-thousand cancer genomes. *Genome Medicine* **13**, 147 (2021). [https://doi.org/10.1186/s13073-021-00964-1](https://doi.org/10.1186/s13073-021-00964-1)

We systematically characterized non-cancer-related pathogenic/likely pathogenic (NC P/LP) germline variants across ~10,000 TCGA samples, investigating their distribution across ancestries, impact on gene expression, and allele-specific expression (ASE) patterns.

---

## Repository Structure

```
.
├── AlleleSpecificExpression/         # ASE analysis and gene enrichment
├── GeneralFunctions/                 # Shared utility functions
├── VariantDistributionAcrossAncestries/  # Cross-ancestry variant analyses
├── VariantImpactOnExp/               # Variant–expression association analyses
└── gnomad/                           # gnomAD comparison and validation
```

---

## Module Descriptions

### `AlleleSpecificExpression/`

Analysis of allele-specific expression associated with rare NC P/LP variants.

| Script | Description |
|--------|-------------|
| `ASE_basic.R` | Identify rare NC P/LPs associated with significant allele-specific expression |
| `ASE_Gene_Enrichment.R` | Gene-level enrichment analysis for NC P/LPs exhibiting significant ASE |
| `ASE_MultiVars.R` | Distribution of NC P/LPs with distinct ASE enrichment status across predicted variant functional classes |
| `GeneVariant_Distribution.R` | Format NC P/LP annotations for lollipop plot visualization via [ProteinPaint](https://proteinpaint.stjude.org/) |

### `GeneralFunctions/`

Shared utility and helper functions used throughout the analyses, including statistical routines, column averaging, and global aesthetic settings for plots.

### `VariantDistributionAcrossAncestries/`

Cross-ancestry characterization of NC P/LP carrier frequencies and variant counts.

| Script | Description |
|--------|-------------|
| `AR.R` | Carrier frequency and variant counts for autosomal recessive (AR) and autosomal dominant (AD) genes across ancestries |
| `Distribution_for_ACMG_classification.R` | NC P/LP carrier frequencies and variant counts stratified by ACMG classification across ancestries |
| `Distribution_for_genes.R` | Per-ancestry carrier frequencies for ACMG 59 genes and the top 10% of genes ranked by aggregate ancestry frequency (excluding Mixed and Other) |

### `VariantImpactOnExp/`

Analyses assessing the impact of NC P/LP variants on gene expression.

| Script | Description |
|--------|-------------|
| `VariantImpactOnExp.R` | Identify genes whose expression is significantly affected by NC P/LPs |
| `VariantImpactOnExp_Plot.R` | Volcano plots for variant–expression associations |
| `PlotPercentileExp.R` | Percentile expression distributions for NC P/LP carriers in specific cancer types; node color encodes variant function, border color encodes ASE enrichment status |
| `PlotPercentileExp_DiffExpSplitCount.R` | Proportion of sample-variants across expression quantile bins stratified by predicted variant function and ASE status |
| `PlotPercentileExp_DiffExpSplitCount_GeneInfo.R` | Detailed sample-variant–level annotations for genes with significant or suggestive expression effects |

### `gnomad/`

Comparison of TCGA variant frequencies with the Genome Aggregation Database (gnomAD).

| Directory / Script | Description |
|--------------------|-------------|
| `bcftools/` | Shell commands for extracting variants of interest from gnomAD VCFs using `bcftools` |
| `bcftools_process/process.R` | Preprocessing of extracted gnomAD variant annotations |
| `bcftools_process/plot_ancestry.R` | Frequency and count distributions of ACMG-classified variants across gnomAD ancestries |
| `ancestry_variant_distribution/count.R` | Ancestry-matched variant counts between TCGA and gnomAD (gnomAD European = FIN ∪ NFE); population-specific NC P/LPs are indicated by triangles, with top variants labeled |
| `ancestry_variant_distribution/statistic.R` | Correlation analyses of variant frequencies between matched TCGA and gnomAD ancestries |

---

## Requirements

- **R** (≥ 3.6)
- Standard Bioconductor and CRAN packages for genomic analysis and visualization

---

## Citation

If this code is useful for your research, please cite:

```bibtex
@article{wang2021noncancer,
  title={Non-cancer-related pathogenic germline variants and expression consequences in ten-thousand cancer genomes},
  author={Wang, Zishan and others},
  journal={Genome Medicine},
  volume={13},
  number={1},
  pages={147},
  year={2021},
  doi={10.1186/s13073-021-00964-1}
}
```

## License

Please refer to the repository license for terms of use.
