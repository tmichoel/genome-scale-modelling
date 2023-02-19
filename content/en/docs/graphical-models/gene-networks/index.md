---
categories: ["Graphical models"]
tags: ["graphical models", "bio"]
title: "Gene regulatory networks"
linkTitle: "Gene regulatory networks"
weight: 1
date: 2023-02-10
description: >
  Gene regulatory networks.
---




{{< alert title="Classic paper" >}}
J Faith *et al.* [*Large-Scale Mapping and Validation of E. coli Transcriptional Regulation from a Compendium of Expression Profiles.*](https://doi.org/10.1371/journal.pbio.0050008) PLOS Biol 5:e8 (2007).


https://www.nature.com/articles/d41586-022-02826-1
{{< /alert >}}

 

## Genes are organized in hierarchical, multi-tissue causal networks

![Civelek and Lusis, Nat Rev Genet (2014)](civelek_nrg_fig3a)

![Albert and Kruglyak, Nat Rev Genet (2016)](albert2016-fig3A)

Variation in expression of one gene has downstream consequences on expression of other genes.

Example: Introduction of just 4 TFs ("Yamanaka factors") converts adult cells into stem cells.

Hundreds to thousands of genes are differentially expressed between cellular states (e.g. healthy vs. disease).

Gene expression in one tissue can affect gene expression in other tissues.

Phenotype variation also causes gene expression variation ("reverse causation").

Reconstruction of causal pathways and gene regulatory networks (GRNs) is essential to understand how the genotype determines the phenotype.

## GRN reconstruction from transcriptomics data

![Gardner and Faith, Phys. Life Rev. (2005)](fig_gardner2005)
![Gardner and Faith, Phys. Life Rev. (2005)](fig_gardner2005-01)

## GRN reconstruction builds on coexpression analysis

![Mackay et al, Nat. Rev. Genet (2009)](nrg2612-f3)

## From correlation to causation

- The aim of network inference is to reconstruct context-specific **causal** gene regulatory networks.

- To infer causality from correlation, additional data is required:
  - Prior information (e.g. gene annotation).
  - Temporal data
  - Instrumental variables (e.g. genetic markers).

![image](fig-causal.pdf)

![image](fig-causal-prior)

- J Faith *et al.* *Large-Scale Mapping and Validation of E. coli Transcriptional Regulation from a Compendium of Expression Profiles.* PLOS Biol 5:e8 (2007). <https://doi.org/10.1371/journal.pbio.0050008>

- V Huynh-Thu *et al.* *Inferring Regulatory Networks from Expression Data Using Tree-Based Methods.* PLOS One 5:e12776 (2010). <https://doi.org/10.1371/journal.pone.0012776>

- N Friedman. *Inferring Cellular Networks Using Probabilistic Graphical Models.* Science 303:799 (2004).<https://doi.org/10.1126/science.1094068>

![image](fig-causal-snp)

- J Zhu *et al.* *An integrative genomics approach to the reconstruction of gene networks in segregating populations.* Cytogenet Genome Res 105:363--374 (2004). <https://doi.org/10.1159/000078209>

- L Wang & T Michoel. *Efficient and accurate causal inference with hidden confounders from genome-transcriptome variation data.* PLOS Comp Biol 13:e1005703 (2017).   <https://doi.org/10.1371/journal.pcbi.1005703>

## Assignment


{{< alert title="Assignment" >}}


{{< /alert >}}