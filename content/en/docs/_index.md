
---
title: "Probabilistic and causal modelling of genome-scale data"
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
- The method used in the classic paper is presented in detail, along with additional methods to solve the same type of problem. The methods are put in practice in a programming assignment. Where possible, original data from the papers studied in the first part is used. 

Four appendices contain the minimum required background knowledge on gene regulation, probability theory, linear algebra, and optimization.

The theoretical sections contain the basic information to understand a method. For more background, try the following textbooks (with free pdfs!), all used in preparation of this course:

- Trevor Hastie, Robert Tibshirani, and Jerome Friedman. [*The Elements of Statistical Learning (second edition)*](https://hastie.su.domains/ElemStatLearn/) (2009).
- Christopher Bishop. [*Pattern Recognition and Machine Learning*](https://www.microsoft.com/en-us/research/publication/pattern-recognition-machine-learning/) (2006).
- Marc Peter Deisenroth, A. Aldo Faisal, and Cheng Soon Ong. [*Mathematics for Machine Learning*](https://mml-book.github.io/) (2020).

The use of classic or path-breaking papers is motivated by [*Back to the future: education for systems-level biologists*](https://pubmed.ncbi.nlm.nih.gov/16990789/). Since the field of genome-scale data analysis is still relatively young, the choice of papers for study is still a bit open and likely to evolve as the course matures.

These lecture are taught as part of the [master program in bioinformatics at UiB](https://www.uib.no/en/studies/MAMN-INF/BI), making up about half of the [BINF301 Genome-scale Algorithms](https://www.uib.no/en/course/BINF301) course. As such, good background knowledge on basic bioinformatics and omics data is assumed. 

{{< alert title="A note on figures and copyright" >}}
One of the objectives of the course is to learn to read and understand scientific papers. Figures from papers selected for discussion are reproduced in these notes. Attribution to the original authors is always given. Where possible [open access](https://en.wikipedia.org/wiki/Open_access) papers are used, but some "classic" papers date from before the birth of open access. If the full text version of the paper is available on [EuropePMC](https://europepmc.org/), through [Unpaywall](http://blog.europepmc.org/2018/04/unlocking-open-europe-pmc-integrates.html) or otherwise, I've reused figures without seeking any further reprint permission. If anyone feels their copyright is violated, please let me know.
{{< /alert >}}
