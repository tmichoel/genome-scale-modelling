
---
title: "Probabilistic and causal modelling in computational biology"
linkTitle: "Lecture notes"
weight: 20
menu:
  main:
    weight: 20
---

Machine learning plays an important role in computational biology. See the [Machine Learning in Computational and Systems Biology community](http://cosi.iscb.org/wiki/MLCSB:Home) or the [Machine Learning in Computational Biology conference series](http://mlcb.org).

These lecture notes focus on **probabilistic** machine learning methods for computational biology, where the experimental data are viewed as random samples from an underlying data-generating process. 

The **"probabilistic modelling"** in the title refers to the use of abstract data-generating processes, not based on any specific biological mechanisms, and derived from generic models and methods. A typical example will be clustering using [Gaussian mixture models](https://en.wikipedia.org/wiki/Mixture_model#Gaussian_mixture_model).

To speak of **"causal modelling"** will require something more, namely that the data-generating process is based on some qualitative prior knowledge or understanding of the true underlying biological process. A typical example will be  [path analysis](https://en.wikipedia.org/wiki/Path_analysis_(statistics)).

The notes are divided in chapters, each focusing on a specific class of methods:

- Clustering
- Regularized regression
- Dimensionality reduction
- Causal inference
- Graphical models
- Spatio-Temporal models

Each chapter follows the same structure:

- A "classic" biological or biomedical research paper is studied where the algorithm (or class of algorithms) of interest was first used. A more recent follow-up or related paper is given as a reading assignment.
- The method used in the classic paper is presented in detail, with additional more modern methods to solve the same type of problem. The methods are put in practice in a programming assignment. Where possible, original data from the papers studied in the first part is used. 

Four appendices contain the minimum required background knowledge on gene regulation, probability theory, linear algebra, and optimization.

The lecture notes are part of the [pcmcb-lectures repository](https://github.com/tmichoel/pcmcb-lectures), which contains linked scripts and notebooks for downloading test data and applying the methods studied in each chapter.

{{< alert title="A note on figures and copyright" >}}
One of the objectives of the course is to learn to read and understand scientific papers. Figures from papers selected for discussion are reproduced in these notes. Attribution to the original authors is always given. Where possible [open access](https://en.wikipedia.org/wiki/Open_access) papers are used, but some "classic" papers date from before the birth of open access. If the full text version of the paper is available on [EuropePMC](https://europepmc.org/), through [Unpaywall](http://blog.europepmc.org/2018/04/unlocking-open-europe-pmc-integrates.html) or otherwise, I've reused figures without seeking any further reprint permission. If anyone feels their copyright is violated, please contact me to take down your content.
{{< /alert >}}
