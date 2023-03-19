---
categories: ["Cluster analysis"]
tags: ["cluster analysis", "algorithms"]
title: "Combinatorial clustering tutorial"
linkTitle: "Tutorial"
date: 2023-02-10
weight: 4
draft: false
description: >
  Illustrating data exploration and K-means clustering using TCGA breast cancer data.
---


```julia
using DrWatson
@quickactivate "Genome-scale ML"

using DataFrames
using CSV
using HypothesisTests
using Statistics
using StatsBase
using StatsPlots, LaTeXStrings
using MultivariateStats
using LinearAlgebra
using Random
using Clustering
using Distances
```


## Import data



See data-preprocessing.jl

```julia
fname = datadir("exp_pro","TCGA-BRCA","TCGA-BRCA-ER-exp-348-med.csv");
df = DataFrame(CSV.File(fname));
ERpos = df.ER .== 1;
```



create a dataframe view that only includes gene expression values

```julia
dfg = @view df[:,3:end];
```



## Data exploration 



### Select most variable genes



Filter variable genes

```julia
sd = map(x -> std(x), eachcol(dfg));
sdcut = 0.7;
dfg = @view dfg[:,sd.>sdcut]
ng = size(dfg,2);
```



Show histogram of standard deviations and cutoff line

```julia
histogram(sd, xlabel="Standard deviation", label="")
vline!([sdcut],linewidth=2,label="")
```

{{< figure src="cluster-analysis-tcga_5_1.png"  >}}

Have a version of the data with standardized genes

```julia
dfgz = mapcols(zscore, dfg);
```



### PCA 



Run PCA to plot basic 2D visualization of the data, and overlay with ER status

```julia
X = Matrix(dfg);
Xz = Matrix(dfgz);
pc = fit(PCA,Xz; maxoutdim=2);

scatter(pc.proj[ERpos,1],pc.proj[ERpos,2], label="ER+")
scatter!(pc.proj[.!ERpos,1],pc.proj[.!ERpos,2], label="ER-")
xlabel!("PC 1")
ylabel!("PC 2")
```

{{< figure src="cluster-analysis-tcga_7_1.png"  >}}

## Data clustering



The PCA figure suggest that samples segregate well according to ER status, and that the signal is detected with (linear PCA), hence K-means is a good choice for clustering.

Run kmeans with two clusters, to analyze overlap with ER

```julia
km = kmeans(X',2);
kmz = kmeans(Xz',2);
```



Compute cluster sizes and overlaps with ER for both clusterings

```julia
counts = [sum(km.assignments.==1 .&& ERpos) sum(km.assignments.==1 .&& .!ERpos)
          sum(km.assignments.==2 .&& ERpos) sum(km.assignments.==2 .&& .!ERpos)
          sum(kmz.assignments.==1 .&& ERpos) sum(kmz.assignments.==1 .&& .!ERpos)
          sum(kmz.assignments.==2 .&& ERpos) sum(kmz.assignments.==2 .&& .!ERpos)
        ]   

scatter(pc.proj[km.assignments.==1,1],pc.proj[km.assignments.==1,2], label="Cluster 1")
scatter!(pc.proj[km.assignments.==2,1],pc.proj[km.assignments.==2,2], label="Cluster 2")
xlabel!("PC 1")
ylabel!("PC 2")
```

{{< figure src="cluster-analysis-tcga_8_1.png"  >}}