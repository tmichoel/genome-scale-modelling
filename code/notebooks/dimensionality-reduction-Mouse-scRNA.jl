### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ 1a16bf2d-2539-42c7-9437-bd7fcaec1708
using DrWatson

# ╔═╡ c349a9cb-a823-4637-ae33-795ecb9e2fa2
# ╠═╡ show_logs = false
quickactivate(@__DIR__)

# ╔═╡ 2fff1e7b-9755-4c0d-9726-3ab0f68ad570
begin
	using DataFrames
	using DataFramesMeta
	using StatsBase
	using Statistics
	using LaTeXStrings
	using Colors
	using TSne
	using LinearAlgebra
	using Distances
	using CSV
	using Arrow
end

# ╔═╡ 27815b60-df9a-11ed-0190-838687447c20
md"# Dimensionality reduction
## Setup the environment
"

# ╔═╡ 19fdc84b-7877-457d-b2ad-76b25c8911f9
md"
## Import data

See the notebook data-preprocessing-Mouse-scRNA.jl for details. Since most of the intial analysis will be at the level of genes, we permute the dimensions of the data file to have genes as columns. 
"

# ╔═╡ cfd650a4-ede6-4dff-b68c-661dc453c591
fexpr = datadir("exp_pro","mouse-brain-single-cell","mouse_ALM_VISp_gene_expression.arrow");

# ╔═╡ f5936ab1-6fa2-439a-8d90-31122e70562e
expr = DataFrame(Arrow.Table(fexpr))

# ╔═╡ 00d3eacc-219d-48e1-95c9-8f4747d09fa4
md"Split off the gene ids:"

# ╔═╡ 34ea9bbc-a9d7-4420-a38d-c763a1c1a3be
genes = expr.Column1

# ╔═╡ 87e8a805-20ff-4fa4-a709-354ff2c7f8d9
select!(expr,Not(:Column1))

# ╔═╡ a4a553d9-ec75-45ce-aa4b-0d739ab37bbf
md"Read the cluster labels:"

# ╔═╡ bcffb5a9-9421-4c8b-999d-5c4888d65a3b
fclust = datadir("exp_raw","mouse-brain-single-cell","tasic-sample_heatmap_plot_data.csv");

# ╔═╡ b372440c-6887-4998-8dd8-2a9dc8f79662
clust = DataFrame(CSV.File(fclust))

# ╔═╡ 702ccafa-7f8d-47ed-abbb-985d63c99aa6
md"## Kobak & Berens pipeline

### Feature selection

Genes that have non-zero expression (that is, expression greater than a threshold $t=32$) in less than $n_{\min}=10$ cells are discarded.
"

# ╔═╡ 8768e307-c97d-4a5e-9bff-c2268dbc813e
t = 32

# ╔═╡ 3a2ed2f2-abd7-4dff-93c9-dcb71608516c
nmin = 10

# ╔═╡ 846a3f22-d216-4678-bd1d-8648dd501756
tf = sum(eachcol(expr .>= t)) .>= nmin

# ╔═╡ 5e407b3f-0881-4425-b845-d7da7d83e601
expr2 = expr[tf,:]

# ╔═╡ f418525b-3063-485d-91ec-606d7094e8f0
sum(ss .> nmin)

# ╔═╡ Cell order:
# ╟─27815b60-df9a-11ed-0190-838687447c20
# ╠═1a16bf2d-2539-42c7-9437-bd7fcaec1708
# ╠═c349a9cb-a823-4637-ae33-795ecb9e2fa2
# ╠═2fff1e7b-9755-4c0d-9726-3ab0f68ad570
# ╟─19fdc84b-7877-457d-b2ad-76b25c8911f9
# ╠═cfd650a4-ede6-4dff-b68c-661dc453c591
# ╠═f5936ab1-6fa2-439a-8d90-31122e70562e
# ╟─00d3eacc-219d-48e1-95c9-8f4747d09fa4
# ╠═34ea9bbc-a9d7-4420-a38d-c763a1c1a3be
# ╠═87e8a805-20ff-4fa4-a709-354ff2c7f8d9
# ╟─a4a553d9-ec75-45ce-aa4b-0d739ab37bbf
# ╠═bcffb5a9-9421-4c8b-999d-5c4888d65a3b
# ╠═b372440c-6887-4998-8dd8-2a9dc8f79662
# ╟─702ccafa-7f8d-47ed-abbb-985d63c99aa6
# ╠═8768e307-c97d-4a5e-9bff-c2268dbc813e
# ╠═3a2ed2f2-abd7-4dff-93c9-dcb71608516c
# ╠═846a3f22-d216-4678-bd1d-8648dd501756
# ╠═5e407b3f-0881-4425-b845-d7da7d83e601
# ╠═f418525b-3063-485d-91ec-606d7094e8f0
