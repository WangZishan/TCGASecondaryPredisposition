# **Non-cancer-related pathogenic germline variants and expression consequences in ten-thousand cancer genomes**

[Citation of our manuscript after published]

Code for the analysis in this project.<br /><br />

## **./AlleleSpecificExpression**

**ASE_basic.R**

Identification of rare NC P/LPs associated with allele specific expression (ASE).

**ASE_Gene_Enrichment.R**

Gene enrichment analysis for NC P/LPs showing significant ASE

**ASE_MultiVars.R**

Plots/Files for distribution of NC P/LPs with distinct ASE enrichment status across predicted variant function classes.

**GeneVariant_Distribution.R**

Extract information of NC P/LPs for genes of interest to generate the format required by online lolliplots software ProteinPaint at https://proteinpaint.stjude.org/.

<br /><br /><br />


## **./GeneralFunctions**

Functions used in this analysis.

<br /><br /><br />



## **./VariantDistributionAcrossAncestries**

**AR.R**

Plots for the carrier frequency and NC P/LPs count of autosomal recessive (AR) and autosomal dominant (AD) genes across ancestries.

**Distribution_for_ACMG_classification.R**

Plots for the frequency of NC P/LP carriers and count of NC P/LPs across ancestries.

**Distribution_for_genes.R**

Plots for the frequency/count of NC P/LP carriers in each ancestry among the ACMG 59 genes and the top 10% genes (ranked by sums of all defined ancestry frequencies, excluding Mix and Other).

<br /><br /><br />




## **./VariantImpactOnExp**

**VariantImpactOnExp.R**

Identification of genes whose expression is affected by related NC P/LPs.

**VariantImpactOnExp_Plot.R**

Volcano plots for genes whose expression is affected by related NC P/LPs.

**PlotPercentileExp.R**

Distribution of percentile expression in a specific cancer at NC P/LP carriers of genes whose expression is significantly/suggestively impacted by NC P/LPs or enriched with significant ASE variants. Color of node represents variant type. Color of node edge represent ASE enrichment status.

**PlotPercentileExp_DiffExpSplitCount.R**

Count/Proportion of sample-variants across expression splits vs predicted variant function/ASE status for genes whose expression is significantly/suggestively impacted by NC P/LPs or enriched with significant ASE variants.

**PlotPercentileExp_DiffExpSplitCount_GeneInfo.R**

Detailed information for sample-variants of for genes whose expression is significantly/suggestively impacted by NC P/LPs or enriched with significant ASE variants.

<br /><br /><br />




## **./gnomad/bcftools**

Command of bcftools to extract information of variant of interests from gnomad dataset.

## **./gnomad/ancestry_variant_distribution**

**count.R**

Variant count of predisposing variants in the matched gnomAD ancestry (European of gnomAD is the union of FIN and NFE populations). TCGA population-specific NC P/LPs, exclusively found in a specific TCGA ancestry, are shown as a triangle. Top NC P/LP or top TCGA ancestry-specific NC P/LP, ranked by allele counts in TCGA or gnomAD, was labelled.

**statistic.R**

(Significance of) Correlations of variant frequencies in the matched ancestries between TCGA and gnomAD.

## **./gnomad/bcftools_process**

**process.R**

Preprocess the information of variant of interests from gnomad dataset.

**plot_ancestry.R**

Plot for the frequency/count of ACMG status/genes across different ancestries.
