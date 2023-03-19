---
categories: ["Statistical significance"]
tags: ["statistical significance", "algo"]
title: "False discovery rate estimation tutorial"
linkTitle: "Tutorial"
weight: 4
draft: true
date: 2023-02-10
description: >
  Illustrating Storey's FDR and q-value estimation method using TCGA breast cancer data.
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
using LinearAlgebra
using Random
using SmoothingSplines
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



Filter variable genes

```julia
sd = map(x -> std(x), eachcol(dfg));
sdcut = 0.7;
dfg = @view dfg[:,sd.>sdcut]
ng = size(dfg,2);
```



## Differential expression analysis



### Theoretical p-values

```julia
pₜ = map(x -> pvalue(UnequalVarianceTTest(x[ERpos], x[.!ERpos])), eachcol(dfg));
```



### Empirical p-values



Compute the t-statistics

```julia
t = map(x -> UnequalVarianceTTest(x[ERpos], x[.!ERpos]).t, eachcol(dfg));
```



Compute empirical p-values

```julia
B = 100; # number of random permutations
trand = zeros(ng,B);
for b = 1:B
    tf = shuffle(ERpos);
    trand[:,b] = map(x -> abs(UnequalVarianceTTest(x[tf], x[.!tf]).t), eachcol(dfg));
end

pₑ = map(x -> sum(trand .>= x)./(ng*B), abs.(t));
```



### Compare empirical and theoretical p-values



Scatter plot

```julia
scatter(pₜ,pₑ,label="",xlabel="Theoretical p-values",ylabel="Empirical p-values")
```

{{< figure src="statistical-significance_8_1.png"  >}}

### Estimate π₀



Because the theoretical and empirical p-values agree, we will work with the theoretical ones, as they don't have any values identically zero.





Get estimate of π₀ for a range of λ values

```julia
λ = 0.01:0.01:0.95;
π₀_vec = map(x -> sum(pₜ.>x)/(ng*(1-x)),λ);
```



Cubic spline smoothing

```julia
spl = fit(SmoothingSpline, λ, π₀_vec, 0.05);
```



Set final estimate

```julia
π₀ = predict(spl,1.0);
```



Reproduce Fig. 1. The inset zooms in on the histogram of p-values greater than 0.5, which mostly consist of null p-values. The line shows the estimate of the proportion of null p-values. It is unusually low for this type of analysis. Due to the very strong difference in global gene expression between ER positive and ER negative breast tumours, we estimate that the majority of genes (>80%) are truely differentially expressed!

```julia
histogram(pₜ, label="", xlabel="Theoretical p-values", normalize=:pdf)
histogram!(
    pₜ,
    xlims=[0.5, 1],
    ylims = [0,1.01],
    label="",
    normalize=:pdf,
    inset = (1, bbox(0.05, 0.05, 0.5, 0.5, :top, :right)),
    subplot = 2,
    bg_inside = nothing
)
hline!([π₀], label="", linewidth=2, subplot=2)
```

{{< figure src="statistical-significance_12_1.png"  >}}

Reproduce Fig. 3

```julia
scatter(λ,π₀_vec,label="",xlabel=L"\lambda", ylabel=L"\pi_0(\lambda)")
plot!(λ,predict(spl,λ),linewidth=3, label="Cubic smoothing spline")
```

{{< figure src="statistical-significance_13_1.png"  >}}

## Estimate q-values



Create ordered p-values

```julia
per = sortperm(pₜ);
pₒ = pₜ[per];
```



Compute q-values

```julia
qₒ = zeros(size(pₒ));
qₒ[end] = π₀ * pₒ[end];
for k=ng-1:-1:1
    qₒ[k] = min(π₀*ng*pₒ[k]/k, qₒ[k+1])
end
```



number of significant genes and expected false positives at every possible cutoff

```julia
nsig = map(x -> sum(qₒ.<=x),qₒ);
nfp = qₒ .* nsig;
```



Reproduce Fig. 2

```julia
l = @layout [a b;c d];
p1 = scatter(
    t[per],qₒ, 
    label="", xlabel="t-statistics", ylabel="q-values",
    markersize=.5,  markercolor=:blue,  markerstrokewidth=0)

p2 = scatter(
    pₒ,qₒ, 
    label="", xlabel="p-values", ylabel="q-values",
    markersize=.5,  markercolor=:blue,  markerstrokewidth=0)

p3 = scatter(
    qₒ, nsig, 
    label="", xlabel="q-values", ylabel="number of significant genes",
    markersize=.5,  markercolor=:blue,  markerstrokewidth=0)

p4 = scatter(
    nsig, nfp, 
    label="", xlabel="number of significant genes", ylabel="number of exp. FP",
    markersize=.5,  markercolor=:blue,  markerstrokewidth=0)

plot(p1, p2, p3, p4, layout=l)
```

{{< figure src="statistical-significance_16_1.png"  >}}
