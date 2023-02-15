---
categories: ["Dimensionality reduction"]
tags: ["dimensionality reduction", "bio"]
title: "Single-cell genomics"
linkTitle: "Single-cell genomics"
weight: 1
date: 2023-02-10
description: >
  Single-cell genomics.
---

{{< alert title="Classic paper" >}}
https://www.nature.com/articles/ncb2709
{{< /alert >}}


{{< alert title="Reference" >}}
A Wagner, A Regev & N Yosef. *Revealing the vectors of cellular identity with Single-cell genomics.* Nature Biotechnology **34**:1145 (2016).

<https://www.nature.com/articles/nbt.3711>
{{< /alert >}}


![Biological and technical factors combine to determine the measured
genomic profiles of single
cells.](single_cell_review_nbt_2016_fig2.png)

## Sources of variation in single-cell RNA-seq

### Technical variation

Due to technical factors in the data generation process.

Main axes of variation in scRNA-seq often dominated by technical factors.

False negatives (expressed but undetected) and false positives (overamplification) are due to minute quantities of starting RNA material and corresponding need for massive amplification.

###  Allele-intrinsic variation

Stochastic factors intrinsic to the molecular mechanisms that control gene expression.

Does not correlate between two alleles of the same gene.

### Allele-extrinsic variation

Factors extrinsic to the process of transcription.

Contribute to establishing differences between cell types or states, either stably or transiently.

Correlated between two alleles of the same gene.

**Most studies aim to understand allele-extrinsic variation and its function, while technical and allele-intrinsic variations are major confounders.**

![Technical confounders of single-cell RNA-seq and computational methods to handle them.](single_cell_review_nbt_2016_fig3.png)

## Addressing overamplification

False-positive detections due to amplification occur at early PCR amplification cycle.

Can be addressed using random molecular tags (RMTs), also called unique molecular identifiers (UMIs):

- Short barcode attached to 3' or 5' end of cDNA molecules before amplification.
- After amplification, the number of unique barcodes associated with reads aligned to a gene/transcript, rather than number of aligned reads serves as the gene/transcript abundance.

Beware:

- Sequencing errors will cause RMT sequences that are close to but not identical, and need to be collapsed.
- "Barcode collision" -- two or more copies of a transcript may receive the same barcode, making them appear as a single copy.

## Addressing false negatives

Limited efficiency of RNA capture and conversion to cDNA leads to dropout events: transcripts that are expressed in the cell but undetected in its mRNA profile.

### Zero-inflated models:

Gene expression is modelled as a mixture of two distributions:

- one in which it is successfully amplified and detected at a level that correlates with its true abundance,
- one in which it is undetected because of technical effects.

(Hidden) component labels can be inferred using EM.

Only cells belonging to the "success" component are used for downstream analysis, e.g. differential expression analysis between    subpopulations of cells.


## Assignment


{{< alert title="Assignment" >}}
Kobak & Berens. *The art of using t-SNE for single-cell transcriptomics*. Nat. Comm. **10**:5416 (2019) <https://doi.org/10.1038/s41467-019-13056-x>

See also: <https://github.com/berenslab/rna-seq-tsne>
{{< /alert >}}