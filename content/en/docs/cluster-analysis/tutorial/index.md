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


## Process TCGA BRCA data



- BRCA expression data available at https://gdc.cancer.gov/about-data/publications/brca_2012, file "BRCA.exp.348.med.txt".
- ER status and cancer stage available in Supplementary Table 1 of the paper at https://www.nature.com/articles/nature11412



I manually converted the clinical data file to XLSX format and have renamed the files as follows:

```julia
fexpr = datadir("exp_raw","TCGA-BRCA","TCGA-BRCA-exp-348-med.txt");
fclin = datadir("exp_raw","TCGA-BRCA","TCGA-BRCA-Supplementary-Tables-1-4.xlsx");
```



### Expression data



Read the data first into a dataframe and then extract the matrix of values, and names of samples and genes. Fix the following:



- Sample IDs in Supp Table only use 12 chars
- Remove genes with missing data 

```julia
df = DataFrame(CSV.File(fexpr));
dropmissing!(df); # drop genes with missing values
df = permutedims(df,1); # turn samples into rows, genes as columns
```



rename samples

```julia
for i=1:lastindex(df[:,1])
    df[i,1] = df[i,1][1:12]
end
```



### Clinical data

```julia
xf = XLSX.readxlsx(fclin); # Read xlsx file
m = xf["SuppTable1"][:][2:end,:]; # Get first sheet, all data, skipping 1st empty row
m_samples = Vector{String}(m[2:end,1]); # the sample names
m_er = Vector{String}(m[2:end,4]); # 4th column is ER status
```



convert ER status to number. Keep only samples with positive or negative ER status.

```julia
pn = unique(m_er)[1:2];
ER = indexin(m_er,pn);
tf = .!isnothing.(ER);
ER = 2*(ER[tf].-1).-1 # Set ER to -1 or +1
df2 = DataFrame("NAME" => m_samples[tf], "ER" => ER);
```



### Merge expression and clinical data

```julia
dat = innerjoin(df2,df, on=:NAME);
```

Save data

```julia
fout = datadir("exp_pro","TCGA-BRCA","TCGA-BRCA-ER-exp-348-med.csv");
CSV.write(fout,dat);
```



## Import data



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