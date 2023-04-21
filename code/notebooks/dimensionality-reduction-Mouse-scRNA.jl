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
	using CSV
	using StatsBase
	using Statistics
	using LaTeXStrings
	using Colors
	using TSne
	using LinearAlgebra
	using Distances
end

# ╔═╡ 27815b60-df9a-11ed-0190-838687447c20
md"# Dimensionality reduction
## Setup the environment
"

# ╔═╡ 19fdc84b-7877-457d-b2ad-76b25c8911f9
md"
## Import data

See the notebook data-preprocessing-Mouse-scRNA.jl for details.
"

# ╔═╡ cfd650a4-ede6-4dff-b68c-661dc453c591
fexpr = datadir("exp_pro","mouse-brain-single-cell","mouse_ALM_VISp_gene_expression.csv");

# ╔═╡ f5936ab1-6fa2-439a-8d90-31122e70562e
expr = DataFrame(CSV.File(fexpr))

# ╔═╡ 90858cc5-9f3b-4b44-a8c2-e4404c2174d5
X = Matrix(expr)

# ╔═╡ bcffb5a9-9421-4c8b-999d-5c4888d65a3b
fclust = datadir("exp_raw","mouse-brain-single-cell","tasic-sample_heatmap_plot_data.csv");

# ╔═╡ b372440c-6887-4998-8dd8-2a9dc8f79662
clust = DataFrame(CSV.File(fclust))

# ╔═╡ 702ccafa-7f8d-47ed-abbb-985d63c99aa6
md"## Kobak & Berens pipeline

### Remove cell and genes with (almost) all zero counts
"

# ╔═╡ 8768e307-c97d-4a5e-9bff-c2268dbc813e
t = 32

# ╔═╡ 3a2ed2f2-abd7-4dff-93c9-dcb71608516c
nmin = 10

# ╔═╡ ccff2da8-e086-4a5f-901a-11807a8e247f


# ╔═╡ Cell order:
# ╟─27815b60-df9a-11ed-0190-838687447c20
# ╠═1a16bf2d-2539-42c7-9437-bd7fcaec1708
# ╠═c349a9cb-a823-4637-ae33-795ecb9e2fa2
# ╠═2fff1e7b-9755-4c0d-9726-3ab0f68ad570
# ╟─19fdc84b-7877-457d-b2ad-76b25c8911f9
# ╠═cfd650a4-ede6-4dff-b68c-661dc453c591
# ╠═f5936ab1-6fa2-439a-8d90-31122e70562e
# ╠═90858cc5-9f3b-4b44-a8c2-e4404c2174d5
# ╠═bcffb5a9-9421-4c8b-999d-5c4888d65a3b
# ╠═b372440c-6887-4998-8dd8-2a9dc8f79662
# ╟─702ccafa-7f8d-47ed-abbb-985d63c99aa6
# ╠═8768e307-c97d-4a5e-9bff-c2268dbc813e
# ╠═3a2ed2f2-abd7-4dff-93c9-dcb71608516c
# ╠═ccff2da8-e086-4a5f-901a-11807a8e247f
