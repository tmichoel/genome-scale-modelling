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




{{< alert title="Classic papers" >}}

J Zhu et al. [*An integrative genomics approach to the reconstruction of gene networks in segregating populations.*](https://doi.org/10.1159/000078209) Cytogenet Genome Res 105:363–374 (2004)

J Faith *et al.* [*Large-Scale Mapping and Validation of E. coli Transcriptional Regulation from a Compendium of Expression Profiles.*](https://doi.org/10.1371/journal.pbio.0050008) PLOS Biol 5:e8 (2007).


V Huynh-Thu et al. [*Inferring Regulatory Networks from Expression Data Using Tree-Based Methods.*](https://doi.org/10.1371/journal.pone.0012776) PLOS One 5:e12776 (2010).



See also a recent perspective on network inference in Nature: [*Smart software untangles gene regulation in cells*](https://www.nature.com/articles/d41586-022-02826-1)
{{< /alert >}}

## An integrative genomics approach to the reconstruction of gene networks in segregating populations

### Figure 1

{{< alert title="Hierarchical clustering of the data set in the gene expression and eQTL dimensions" >}}
![Hierarchical clustering of the data set in the gene expression and eQTL dimensions](zhu2004-fig1.jpeg)



Figure obtained from [full text on EuropePMC](https://europepmc.org/article/med/15237224).
{{< /alert >}}

### Figure 3

{{< alert title="Sub-networks associated with Hsd11b1" >}}
![Sub-networks associated with Hsd11b1](zhu2004-fig3.jpeg)



Figure obtained from [full text on EuropePMC](https://europepmc.org/article/med/15237224).
{{< /alert >}}

## Large-Scale Mapping and Validation of E. coli Transcriptional Regulation from a Compendium of Expression Profiles

### Software

- [Original (Matlab)](http://m3d.mssm.edu/network\_inference.html)
- [R/Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/minet.html)

### Figure 1

{{< alert title="Overview of the approach" >}}
![Overview of the approach](faith2007-fig1.png)



Figure obtained from [full text on EuropePMC](https://europepmc.org/article/med/17214507).
{{< /alert >}}

### Figure 2

{{< alert title="The CLR algorithm" >}}
![The CLR algorithm](faith2007-fig2.png)



Figure obtained from [full text on EuropePMC](https://europepmc.org/article/med/17214507).
{{< /alert >}}


### Figure 5

{{< alert title="Experimental Validation of Inferred Regulatory Interactions" >}}
![The Transcriptional Regulatory Map Inferred by CLR with an Estimated 60% Precision](faith2007-fig5.png)



Figure obtained from [full text on EuropePMC](https://europepmc.org/article/med/17214507).
{{< /alert >}}

## Inferring Regulatory Networks from Expression Data Using Tree-Based Methods

### Software

[Genie3 (Python, Matlab, R)](https://github.com/vahuynh/GENIE3)

### Figure 1

{{< alert title="GENIE3 procedure" >}}
![GENIE3 procedure](HuynhThu2010-fig1.png)



Figure obtained from [full text on EuropePMC](https://europepmc.org/article/med/20927193).
{{< /alert >}}


### Figure 4

{{< alert title="Precision-Recall curves for the E. coli network" >}}
![Precision-Recall curves for the E. coli network](HuynhThu2010-fig4.png)



Figure obtained from [full text on EuropePMC](https://europepmc.org/article/med/20927193).
{{< /alert >}}
