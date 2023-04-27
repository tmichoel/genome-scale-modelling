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
	using StatsPlots
	using MultivariateStats
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

See the notebook data-preprocessing-Mouse-scRNA.jl for details.
"

# ╔═╡ cfd650a4-ede6-4dff-b68c-661dc453c591
fexpr = datadir("exp_pro","mouse-brain-single-cell","mouse_ALM_VISp_gene_expression.arrow");

# ╔═╡ f5936ab1-6fa2-439a-8d90-31122e70562e
expr = DataFrame(Arrow.Table(fexpr));

# ╔═╡ 00d3eacc-219d-48e1-95c9-8f4747d09fa4
md"Split off the gene ids:"

# ╔═╡ 34ea9bbc-a9d7-4420-a38d-c763a1c1a3be
genes = expr.Column1;

# ╔═╡ 87e8a805-20ff-4fa4-a709-354ff2c7f8d9
select!(expr,Not(:Column1));

# ╔═╡ a4a553d9-ec75-45ce-aa4b-0d739ab37bbf
md"Read the cluster labels:"

# ╔═╡ bcffb5a9-9421-4c8b-999d-5c4888d65a3b
fclust = datadir("exp_raw","mouse-brain-single-cell","tasic-sample_heatmap_plot_data.csv");

# ╔═╡ b372440c-6887-4998-8dd8-2a9dc8f79662
clust = DataFrame(CSV.File(fclust));

# ╔═╡ 702ccafa-7f8d-47ed-abbb-985d63c99aa6
md"## Kobak & Berens pipeline

### Feature selection

Genes that have non-zero expression (that is, expression greater than a threshold $t=32$) in less than $n_{\min}=10$ cells are discarded. We start by computing the library depth per million for each cell as it will be needed later.
"

# ╔═╡ 7c7cd501-118c-4cc6-8184-5eab8c950744
libraryDepth = map( x -> sum(x), eachcol(expr)) / 1e6

# ╔═╡ 8768e307-c97d-4a5e-9bff-c2268dbc813e
t = 32;

# ╔═╡ 3a2ed2f2-abd7-4dff-93c9-dcb71608516c
nmin = 10;

# ╔═╡ 846a3f22-d216-4678-bd1d-8648dd501756
tf = sum(eachcol(expr .>= t)) .>= nmin;

# ╔═╡ 75f7d655-a5a8-407a-94da-14a2e79922a9
md"After filtering the genes, the dataframe is of sufficiently reduced size that swapping rows and columns goes relatively quickly. Having genes as columns will speed up the remaining analyses."

# ╔═╡ 5e407b3f-0881-4425-b845-d7da7d83e601
expr2 = permutedims(expr[tf,:]);

# ╔═╡ 4bd5c51f-9d46-41dc-825d-19e1fe47242e
md"Compute the fraction of near-zero counts and the mean log2 expression over the non-zero cells for each gene. Following Kobak & Berens, this is done on the *raw* counts before sequencing depth normalization."

# ╔═╡ 717fc8d5-742c-4d26-9f7c-285229bfd18e
ncell = nrow(expr2);

# ╔═╡ 88300ee7-fdd4-4c64-91d4-1e58acc09ad2
d = map( x -> sum(x.<=t), eachcol(expr2) ) / ncell;

# ╔═╡ f8a64c96-ccba-41c9-afdc-797f13f475c5
m = map( x -> mean(log2.(x[x.>t])), eachcol(expr2) );

# ╔═╡ aadeb190-5c16-4d90-bdcb-5a741d9e663d
md"Use the same formula as Kobak & Berens to select genes, with parameters from their Supp Fig 4."

# ╔═╡ 218f2608-eddb-4e42-a461-ef988b4a03c9
begin
	a = 1.5;
	b = 6.56;
	featureSelect = vec(d .> exp.(-a.*(m.-b)) .+ 0.02);
end

# ╔═╡ 66676aa3-eb8c-4a29-beae-4aa6e06360c4
md"Reproduce Kobak & Berens Supp Fig 4:"

# ╔═╡ c31d8827-a1ef-4ae0-b70d-50ad3ebe9899
begin
	x = range(minimum(m),maximum(m),length=100);
	scatter(m,d,
		label = "",
		xlabel = "Mean log2 nonzero expression",
		ylabel = "Frequency of nonzero expression"
	)
	plot!(x, exp.(-a.*(x.-b)) .+ 0.02, color=:red, linewidth=3, label="")
	ylims!(0, 1.0)
	annotate!(12.5,0.2,L"y=\exp(-1.5(x-6.56))+0.02", :red)
end

# ╔═╡ 8ba73802-e226-4fe5-a664-d983a035003d
md"### Non-linear transformation

Divide the read counts of each cell by the cell's sequencing depth per million and transform to logarithmic scale. The sequencing depth is counted across all genes (see above), but only the previously selected genes are retained in the final matrix."

# ╔═╡ 533b0237-4c52-4b36-bd49-160d7769f4f8
Y = log2.( select(expr2, featureSelect) ./ libraryDepth .+ 1);

# ╔═╡ 8e4f3869-0279-4da8-a140-184a8afc55b9
md"## Dimension reduction
### t-SNE

Run t-SNE on a reduced dataset:
"

# ╔═╡ 7e4dbad9-48c3-4eef-babd-e90f9d21f5c4
subc = 1:5:nrow(Y);

# ╔═╡ 1790684d-2715-4050-bc42-7c519bc7768b
X = Matrix(Y)[subc,:];

# ╔═╡ 1ec5357a-5266-42cf-a745-b272a6f99c48
# ╠═╡ show_logs = false
Z = tsne(X, 2, 50, 1000, 20.0);

# ╔═╡ ca3c4d86-c3d6-4556-8b2b-c8cdf496f615
md"### PCA"

# ╔═╡ 04223165-df9f-4a7c-a823-592f0525bd6e
P = fit(PCA, X; maxoutdim=2);

# ╔═╡ 5f678b8a-fbab-47c1-a536-4b962bbec068
md"## Figures"

# ╔═╡ 65c3a6da-66af-4598-94fa-01181ab724bc
begin
	cellID = names(expr)[subc];
	ii = indexin(cellID,clust.sample_name);
	clid = clust.cluster_id[ii];
	clcol = clust.cluster_color[ii];
	uu = unique(clcol);
end

# ╔═╡ df4a4c1a-3988-4677-a6bd-8f7d3e3adf05
begin
	f = scatter(Z[:,1],Z[:,2],
		label="",
		xlabel="t-SNE 1",
		ylabel="t-SNE 2",
		markerstrokewidth = 0.5
	)
	for col in uu
	    sel = isequal.(clcol,col)
	    scatter!(Z[sel,1],Z[sel,2],
			color=parse(Colorant, col),
			label="",
			markerstrokewidth = 0.5
		)
	end
	f
end

# ╔═╡ 3e3a1cbe-4a09-41b5-a361-f557c6fdd8f0
begin
	f2 = scatter(P.proj[:,1],P.proj[:,2],
		label="",
		xlabel="PCA 1",
		ylabel="PCA 2",
		markerstrokewidth = 0.5
	)
	for col in uu
	    sel = isequal.(clcol,col)
	    scatter!(P.proj[sel,1],P.proj[sel,2],
			color=parse(Colorant, col),
			label="",
			markerstrokewidth = 0
		)
	end
	f2
end

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
# ╠═7c7cd501-118c-4cc6-8184-5eab8c950744
# ╠═8768e307-c97d-4a5e-9bff-c2268dbc813e
# ╠═3a2ed2f2-abd7-4dff-93c9-dcb71608516c
# ╠═846a3f22-d216-4678-bd1d-8648dd501756
# ╟─75f7d655-a5a8-407a-94da-14a2e79922a9
# ╠═5e407b3f-0881-4425-b845-d7da7d83e601
# ╟─4bd5c51f-9d46-41dc-825d-19e1fe47242e
# ╠═717fc8d5-742c-4d26-9f7c-285229bfd18e
# ╠═88300ee7-fdd4-4c64-91d4-1e58acc09ad2
# ╠═f8a64c96-ccba-41c9-afdc-797f13f475c5
# ╟─aadeb190-5c16-4d90-bdcb-5a741d9e663d
# ╠═218f2608-eddb-4e42-a461-ef988b4a03c9
# ╟─66676aa3-eb8c-4a29-beae-4aa6e06360c4
# ╠═c31d8827-a1ef-4ae0-b70d-50ad3ebe9899
# ╟─8ba73802-e226-4fe5-a664-d983a035003d
# ╠═533b0237-4c52-4b36-bd49-160d7769f4f8
# ╟─8e4f3869-0279-4da8-a140-184a8afc55b9
# ╠═7e4dbad9-48c3-4eef-babd-e90f9d21f5c4
# ╠═1790684d-2715-4050-bc42-7c519bc7768b
# ╠═1ec5357a-5266-42cf-a745-b272a6f99c48
# ╟─ca3c4d86-c3d6-4556-8b2b-c8cdf496f615
# ╠═04223165-df9f-4a7c-a823-592f0525bd6e
# ╟─5f678b8a-fbab-47c1-a536-4b962bbec068
# ╠═65c3a6da-66af-4598-94fa-01181ab724bc
# ╠═df4a4c1a-3988-4677-a6bd-8f7d3e3adf05
# ╠═3e3a1cbe-4a09-41b5-a361-f557c6fdd8f0
