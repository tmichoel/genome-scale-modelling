### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 3da358be-6ab5-429f-a14c-44fff2b39bb3
using DrWatson

# ╔═╡ 332265be-40ed-4699-8aab-adfa846f0041
# ╠═╡ show_logs = false
quickactivate(@__DIR__)

# ╔═╡ a5e35cb1-7ae0-4381-8c0b-2c424fce2e4f
begin
	using DataFrames
	using CSV
	using HypothesisTests
	using Statistics
	using StatsBase
	using StatsPlots, LaTeXStrings
	using LinearAlgebra
	using Random
	using SmoothingSplines
	using PlutoUI
	using Printf
end

# ╔═╡ 6f0288d0-dd2c-11ed-30b4-05a5fdd0c0c5
md"# Statistical significance for genome-wide studies
## Setup the environment
"

# ╔═╡ 1a0e649d-89bb-48e4-b7ba-2886588c47d5
md" ## Import data
See the notebook data-preprocessing-TCGA-BRCA.jl
"

# ╔═╡ dee6dc57-1020-4293-9356-538a7f0063f7
fname = datadir("exp_pro","TCGA-BRCA","TCGA-BRCA-ER-exp-348-med.csv");

# ╔═╡ 8c4ce7e1-9dfb-441b-996a-4feca8ab5aee
df = DataFrame(CSV.File(fname));

# ╔═╡ 427ef119-90be-4a41-b9df-11a618913a66
ERpos = df.ER .== 1;

# ╔═╡ b9088e82-6a40-4950-b4ce-4835d151fd95
md"Create a data frame with only gene expression values"

# ╔═╡ 1a92fc7e-6e6f-4fa6-9618-c540b7bcf7ad
dfg = @view df[:,3:end];

# ╔═╡ 55e48b8d-dadd-4c66-b890-2cfbe0ed1089
ns,ng = size(dfg);

# ╔═╡ b8badfe2-8df5-4978-ba78-c38ecd21abe5
md"## Differential expression analysis"

# ╔═╡ 5ba6b68d-5696-4af7-bdb6-deeb5f75cf6b
md"There are $(ng) genes and $(ns) samples in the data. There are $(sum(ERpos)) ER positive and $(sum(.!ERpos)) ER negative samples. We wish to identify genes that are differentially expressed between ER positive and negative samples."

# ╔═╡ 3a7a668c-6c1d-4658-8601-9da1eb47cd3c
md"### Theoretical p-values
We compute theoretical differential expression p-values using a t-test for each gene.
"

# ╔═╡ 157afc07-5e21-41fb-83e9-371e7a9828dd
pₜ = map(x -> pvalue(UnequalVarianceTTest(x[ERpos], x[.!ERpos])), eachcol(dfg));

# ╔═╡ 4cd49e1f-1e16-4081-940a-6dfdd447cfde
md"### Empirical p-values
We compute empirical p-values from the t-test statistics using permuted data.

First compute the t-statistics:
"

# ╔═╡ f8815ef9-40f6-4095-9b8b-2e689ffb3656
t = map(x -> UnequalVarianceTTest(x[ERpos], x[.!ERpos]).t, eachcol(dfg));

# ╔═╡ 581c1761-1b98-45d2-bd56-561f5f15beb0
md"Now compute the empirical p-values:"

# ╔═╡ ae2252af-23cf-444b-ac04-eb1ea1b4bbc9
@bind B Slider(20:20:200)

# ╔═╡ 1131a0f7-9776-4b30-adcb-a66d8b09b1ed
md"Move the slider to select the number of random permutations. B = $(B) selected"

# ╔═╡ 6baefc7c-9a2c-4205-a989-6b42c4c6e525
begin
	trand = zeros(ng,B);
	for b = 1:B
	    tf = shuffle(ERpos);
	    trand[:,b] = map(x -> abs(UnequalVarianceTTest(x[tf], x[.!tf]).t), eachcol(dfg));
	end
end

# ╔═╡ b70d0bf7-629d-42e2-b6a5-6606685eb928
pₑ = map(x -> sum(trand .>= x)./(ng*B), abs.(t));

# ╔═╡ 589be755-6e77-4cd3-9d80-7c2652917ff5
md"### Compare empirical and theoretical p-values"

# ╔═╡ de22ae5e-4a08-43bc-89a2-f5fdeb2bfa29
scatter(pₜ,pₑ,label="",xlabel="Theoretical p-values",ylabel="Empirical p-values")

# ╔═╡ ec841bd7-38f7-43e3-be33-5a86f24d55bb
md"### Estimate $\pi_0$
Because the theoretical and empirical p-values agree, we will work with the theoretical ones, as they don't have any values identically zero.

We implement the method from Storey & Tibshirani.

First get an estimate of $\pi_0$ for a range of $\lambda$ values.
"

# ╔═╡ 92fa03c3-8e99-4f3e-ae3c-4c20cc90c5d4
λ = 0.01:0.01:0.95;

# ╔═╡ 7252cef1-19eb-49f2-ac79-64ea65d61608
π₀_vec = map(x -> sum(pₜ.>x)/(ng*(1-x)),λ);

# ╔═╡ 37e48fb8-0db9-4fe3-885a-71dfbed06be2
md"Cubic spline smoothing"

# ╔═╡ 95bc389b-462a-4cb4-9b52-f5016ea8b4cb
@bind sm Slider(0:0.005:0.02)

# ╔═╡ 6fd6d426-2641-4b9f-acc0-aa1b19ea09a3
spl = fit(SmoothingSpline, λ, π₀_vec, sm);

# ╔═╡ 155e9fea-e277-4711-8df7-d6ee677b5523
md"Extrapolate to get the final estimate:"

# ╔═╡ 24975cba-b06f-4b12-8807-feab36ff53c0
π₀ = predict(spl,1.0);

# ╔═╡ 24b7f603-2683-427f-9b6b-977212d2fd2a
s0 = @sprintf("%1.2f", π₀);

# ╔═╡ 5d64ea0d-d62d-4d87-99c8-1d64912794de
md"Reproduce Fig. 1. The inset zooms in on the histogram of p-values greater than 0.5, which mostly consist of null p-values. The line shows the estimate of the proportion of null p-values ($\pi_0$=$(s0)). It is unusually low for this type of analysis. Due to the very strong difference in global gene expression between ER positive and ER negative breast tumours, we estimate that the majority of genes are truely differentially expressed!"

# ╔═╡ 1321f1e7-2f53-48ce-81a0-98ab5ece6235
begin
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
end

# ╔═╡ 6045e0ad-6c97-4458-bb17-0c095ce031e6
md"Reproduce Fig. 3."

# ╔═╡ a57cc78f-8cdb-4c9d-8b15-7a63ef6fd066
begin
	scatter(λ,π₀_vec,label="",xlabel=L"\lambda", ylabel=L"\pi_0(\lambda)")
	plot!(λ,predict(spl,λ),linewidth=3, label="Cubic smoothing spline")
end

# ╔═╡ 7ff37ea1-8729-4b3a-a1e4-87296af7afd2
md"The spline interpolation curve (red line) depends on a smoothing parameter which determines how closely the line follows the data points being fitted. Go back up and move the slider to see the effect of changing this parameter."

# ╔═╡ 4756c955-c72e-4402-a4f2-51612517e701
md"### Estimate q-values
First order the p-values:
"

# ╔═╡ 39b70507-771e-4c54-9358-4c4566010ff1
begin
	per = sortperm(pₜ);
	pₒ = pₜ[per];
end

# ╔═╡ 8d56bf64-ea33-470a-b1b6-e38fb0fa160c
md"Now compute the q-values using the formulae from Storey & Tibshirani:"

# ╔═╡ 4357e307-2836-4ca7-a4c7-35a1e3aaf787
begin
	qₒ = zeros(size(pₒ));
	qₒ[end] = π₀ * pₒ[end];
	for k=ng-1:-1:1
	    qₒ[k] = min(π₀*ng*pₒ[k]/k, qₒ[k+1])
	end
end

# ╔═╡ 34e5c904-df1d-4609-8cee-1c58225a998c
md"The number of significant genes and expected false positives at every possible cutoff is calcuated as follows:"

# ╔═╡ dd56d34e-163b-46f4-af24-6fbc5388eafb
nsig = map(x -> sum(qₒ.<=x),qₒ);

# ╔═╡ 59d70b1f-914b-47da-9d95-b2e7551c27f3
nfp = qₒ .* nsig;

# ╔═╡ a66a3bdc-6e42-455e-82b5-c03dc43cb194
md"We can now reproduce Fig. 1."

# ╔═╡ 8cb20a58-c329-439f-8769-85c588d1c4ff
begin
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
end

# ╔═╡ Cell order:
# ╟─6f0288d0-dd2c-11ed-30b4-05a5fdd0c0c5
# ╠═3da358be-6ab5-429f-a14c-44fff2b39bb3
# ╠═332265be-40ed-4699-8aab-adfa846f0041
# ╠═a5e35cb1-7ae0-4381-8c0b-2c424fce2e4f
# ╟─1a0e649d-89bb-48e4-b7ba-2886588c47d5
# ╠═dee6dc57-1020-4293-9356-538a7f0063f7
# ╠═8c4ce7e1-9dfb-441b-996a-4feca8ab5aee
# ╠═427ef119-90be-4a41-b9df-11a618913a66
# ╟─b9088e82-6a40-4950-b4ce-4835d151fd95
# ╠═1a92fc7e-6e6f-4fa6-9618-c540b7bcf7ad
# ╠═55e48b8d-dadd-4c66-b890-2cfbe0ed1089
# ╟─b8badfe2-8df5-4978-ba78-c38ecd21abe5
# ╟─5ba6b68d-5696-4af7-bdb6-deeb5f75cf6b
# ╟─3a7a668c-6c1d-4658-8601-9da1eb47cd3c
# ╠═157afc07-5e21-41fb-83e9-371e7a9828dd
# ╟─4cd49e1f-1e16-4081-940a-6dfdd447cfde
# ╠═f8815ef9-40f6-4095-9b8b-2e689ffb3656
# ╟─581c1761-1b98-45d2-bd56-561f5f15beb0
# ╟─ae2252af-23cf-444b-ac04-eb1ea1b4bbc9
# ╟─1131a0f7-9776-4b30-adcb-a66d8b09b1ed
# ╠═6baefc7c-9a2c-4205-a989-6b42c4c6e525
# ╠═b70d0bf7-629d-42e2-b6a5-6606685eb928
# ╟─589be755-6e77-4cd3-9d80-7c2652917ff5
# ╠═de22ae5e-4a08-43bc-89a2-f5fdeb2bfa29
# ╟─ec841bd7-38f7-43e3-be33-5a86f24d55bb
# ╠═92fa03c3-8e99-4f3e-ae3c-4c20cc90c5d4
# ╠═7252cef1-19eb-49f2-ac79-64ea65d61608
# ╟─37e48fb8-0db9-4fe3-885a-71dfbed06be2
# ╠═95bc389b-462a-4cb4-9b52-f5016ea8b4cb
# ╠═6fd6d426-2641-4b9f-acc0-aa1b19ea09a3
# ╟─155e9fea-e277-4711-8df7-d6ee677b5523
# ╠═24975cba-b06f-4b12-8807-feab36ff53c0
# ╟─24b7f603-2683-427f-9b6b-977212d2fd2a
# ╟─5d64ea0d-d62d-4d87-99c8-1d64912794de
# ╠═1321f1e7-2f53-48ce-81a0-98ab5ece6235
# ╟─6045e0ad-6c97-4458-bb17-0c095ce031e6
# ╠═a57cc78f-8cdb-4c9d-8b15-7a63ef6fd066
# ╟─7ff37ea1-8729-4b3a-a1e4-87296af7afd2
# ╟─4756c955-c72e-4402-a4f2-51612517e701
# ╠═39b70507-771e-4c54-9358-4c4566010ff1
# ╟─8d56bf64-ea33-470a-b1b6-e38fb0fa160c
# ╠═4357e307-2836-4ca7-a4c7-35a1e3aaf787
# ╟─34e5c904-df1d-4609-8cee-1c58225a998c
# ╠═dd56d34e-163b-46f4-af24-6fbc5388eafb
# ╠═59d70b1f-914b-47da-9d95-b2e7551c27f3
# ╟─a66a3bdc-6e42-455e-82b5-c03dc43cb194
# ╠═8cb20a58-c329-439f-8769-85c588d1c4ff
