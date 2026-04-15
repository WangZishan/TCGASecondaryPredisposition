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
| `ASE_basic.R` | Identification of rare NC P/LPs associated with allele specific expression (ASE). |
| `ASE_Gene_Enrichment.R` | Gene enrichment analysis for NC P/LPs showing significant ASE. |
| `ASE_MultiVars.R` | Plots/Files for distribution of NC P/LPs with distinct ASE enrichment status across predicted variant function classes. |
| `GeneVariant_Distribution.R` | Extract information of NC P/LPs for genes of interest to generate the format required by online lolliplots software [ProteinPaint](https://proteinpaint.stjude.org/). |

### `GeneralFunctions/`

Functions used in this analysis.

### `VariantDistributionAcrossAncestries/`

Cross-ancestry characterization of NC P/LP carrier frequencies and variant counts.

| Script | Description |
|--------|-------------|
| `AR.R` | Plots for the carrier frequency and NC P/LPs count of autosomal recessive (AR) and autosomal dominant (AD) genes across ancestries. |
| `Distribution_for_ACMG_classification.R` | Plots for the frequency of NC P/LP carriers and count of NC P/LPs across ancestries. |
| `Distribution_for_genes.R` | Plots for the frequency/count of NC P/LP carriers in each ancestry among the ACMG 59 genes and the top 10% genes (ranked by sums of all defined ancestry frequencies, excluding Mix and Other). |

### `VariantImpactOnExp/`

Analyses assessing the impact of NC P/LP variants on gene expression.

| Script | Description |
|--------|-------------|
| `VariantImpactOnExp.R` | Identification of genes whose expression is affected by related NC P/LPs. |
| `VariantImpactOnExp_Plot.R` | Volcano plots for genes whose expression is affected by related NC P/LPs. |
| `PlotPercentileExp.R` | Distribution of percentile expression in a specific cancer at NC P/LP carriers of genes whose expression is significantly/suggestively impacted by NC P/LPs or enriched with significant ASE variants. Color of node represents variant type. Color of node edge represent ASE enrichment status. |
| `PlotPercentileExp_DiffExpSplitCount.R` | Count/Proportion of sample-variants across expression splits vs predicted variant function/ASE status for genes whose expression is significantly/suggestively impacted by NC P/LPs or enriched with significant ASE variants. |
| `PlotPercentileExp_DiffExpSplitCount_GeneInfo.R` | Detailed information for sample-variants of for genes whose expression is significantly/suggestively impacted by NC P/LPs or enriched with significant ASE variants. |

### `gnomad/`

Comparison of TCGA variant frequencies with the Genome Aggregation Database (gnomAD).

| Directory / Script | Description |
|--------------------|-------------|
| `bcftools/` | Command of bcftools to extract information of variant of interests from gnomad dataset. |
| `bcftools_process/process.R` | Preprocess the information of variant of interests from gnomad dataset. |
| `bcftools_process/plot_ancestry.R` | Plot for the frequency/count of ACMG status/genes across different ancestries. |
| `ancestry_variant_distribution/count.R` | Variant count of predisposing variants in the matched gnomAD ancestry (European of gnomAD is the union of FIN and NFE populations). TCGA population-specific NC P/LPs, exclusively found in a specific TCGA ancestry, are shown as a triangle. Top NC P/LP or top TCGA ancestry-specific NC P/LP, ranked by allele counts in TCGA or gnomAD, was labelled. |
| `ancestry_variant_distribution/statistic.R` | (Significance of) Correlations of variant frequencies in the matched ancestries between TCGA and gnomAD. |

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
