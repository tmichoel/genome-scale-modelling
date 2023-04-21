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

# ╔═╡ 2601aee4-e4ed-4162-a5ab-9130ba2703cd
using DrWatson

# ╔═╡ 03ec1c38-bc72-4fac-b6b2-dc9f593e73e4
# ╠═╡ show_logs = false
quickactivate(@__DIR__)

# ╔═╡ 0f27ea66-102c-4a58-ba4c-77d651d47cc1
begin
	using DataFrames
	using CSV
	using Statistics
	using StatsBase
	using StatsPlots, LaTeXStrings
	using MultivariateStats
	using LinearAlgebra
	using Random
	using Clustering
	using Distances
	using PlutoUI
	using Printf
end

# ╔═╡ 49552810-dd16-11ed-0e3c-efcdc0bbac23
md"# Cluster analysis
## Setup the environment
"

# ╔═╡ 50bca4fb-a282-4b80-afd9-b2187d36fc24
md"
## Import data

See the notebook data-preprocessing-TCGA-BRCA.jl for details.
"

# ╔═╡ dff97de6-1a10-417b-a58c-12bf024d4433
begin
	fname = datadir("exp_pro","TCGA-BRCA","TCGA-BRCA-ER-exp-348-med.csv");
	df = DataFrame(CSV.File(fname));
	ERpos = df.ER .== 1;
end

# ╔═╡ 83e0f4f6-e6b9-4115-aacd-61e3ca3ce702
md"Create a view that only includes gene expression values"

# ╔═╡ 28407ba6-d258-4670-8956-4d3660124010
dfg = @view df[:,3:end];

# ╔═╡ 8f97df78-ad41-45b7-a169-c1ada2bdb8c9
md"## Filter variable genes
In cluster analysis we often standardize features (genes in our case) to give all of them equal weight in the Euclidean distance measure between observations (breast tumours in our case), see Elements of Statistical Learning section ...

In the figure below, we see a histogram of standard deviations of all genes. Move the slider to decide on a standard deviation cutoff where only genes with standard deviation above the cutoff are retained.
"

# ╔═╡ d1133895-2f79-4e40-8310-d88277b0e643
@bind sdcut Slider(0:0.1:2)

# ╔═╡ cfc519e3-f6e4-4e69-b46f-a3140a293262
begin
	sd = map(x -> std(x), eachcol(dfg));
	ng = sum(sd .> sdcut)
	histogram(sd, xlabel="Standard deviation", label="")
	vline!([sdcut],linewidth=2,label="")
	annotate!(2.5,1700,@sprintf("Standard deviation cutoff: %1.1f", sdcut))
	annotate!(2.5,1400,@sprintf("Number of selected genes: %d", ng))
end

# ╔═╡ 9d08f7b3-2536-496c-9bbc-93500e78cffd
dfgz = mapcols(zscore, @view dfg[:,sd.>sdcut]);

# ╔═╡ 6fee0058-c354-4813-99f9-226f099e58c4
ns,ngv = size(dfgz);

# ╔═╡ 6c36c0f0-4eb6-419d-859a-bea6528c3522
md"Make a version of the data with filtered and standardized genes. The filtered data consists of expression values for $(ngv) genes in $(ns) samples."

# ╔═╡ a2f7dfaa-26e4-41b9-8ab4-6f2123f8dcea
md"## PCA
Let's start with some basic data exploration by running PCA.
"

# ╔═╡ 8510cc2f-2cc1-4329-8cc2-d9a12441a231
Xz = Matrix(dfgz);

# ╔═╡ 8887b407-9604-42a7-bdf6-96deb7b89b0e
pc = fit(PCA,Xz; maxoutdim=2);

# ╔═╡ 20619594-5a0f-4f29-9ef3-517641ee7ee9
md"A scatter plot of the first two principal components with samples colored by ER status suggests that expression profiles of breast tumours segregate well according to ER status:"

# ╔═╡ d78b38fe-6adf-4da7-a8d4-64bb2bffd3fd
begin
	scatter(pc.proj[ERpos,1],pc.proj[ERpos,2], label="ER+")
	scatter!(pc.proj[.!ERpos,1],pc.proj[.!ERpos,2], label="ER-")
	xlabel!("PC 1")
	ylabel!("PC 2")
end

# ╔═╡ f9cc2e13-a2ca-468c-9e17-a2335972ab82
md"## K-means clustering
We use K-means clustering with $K=2$ to analyze the overlap between expression clusters and ER status. 
"

# ╔═╡ 0b2b9265-648f-47f3-ab74-200b568156e0
kmz = kmeans(Xz',2);

# ╔═╡ 73ad232d-f45f-416f-8ce4-e3b817a8493a
md"Compute the overlap count for each cluster with ER positive and negative samples."

# ╔═╡ a8b14f6b-eade-4510-9e42-c80fb83aa413
begin
	counts = [
	          sum(kmz.assignments.==1 .&& ERpos) sum(kmz.assignments.==1 .&& .!ERpos)
	          sum(kmz.assignments.==2 .&& ERpos) sum(kmz.assignments.==2 .&& .!ERpos)
	];
	dfc = DataFrame("ER positive" => counts[:,1], "ER negative" => counts[:,2])
end

# ╔═╡ 07e190ab-6094-42be-8570-1b585e643bcc
md"Repeating the scatter plot of the first two principal components, this time coloured by cluster membership, confirms that cluster 1 consists largely of the ER negative samples and cluster 2 of the ER positive samples."

# ╔═╡ 0ed03e10-9897-4fdd-a0d1-18d75eab408b
begin
	scatter(pc.proj[kmz.assignments.==1,1],pc.proj[kmz.assignments.==1,2], label="Cluster 1")
	scatter!(pc.proj[kmz.assignments.==2,1],pc.proj[kmz.assignments.==2,2], label="Cluster 2")
	xlabel!("PC 1")
	ylabel!("PC 2")
end

# ╔═╡ Cell order:
# ╠═49552810-dd16-11ed-0e3c-efcdc0bbac23
# ╠═2601aee4-e4ed-4162-a5ab-9130ba2703cd
# ╠═03ec1c38-bc72-4fac-b6b2-dc9f593e73e4
# ╠═0f27ea66-102c-4a58-ba4c-77d651d47cc1
# ╠═50bca4fb-a282-4b80-afd9-b2187d36fc24
# ╠═dff97de6-1a10-417b-a58c-12bf024d4433
# ╟─83e0f4f6-e6b9-4115-aacd-61e3ca3ce702
# ╠═28407ba6-d258-4670-8956-4d3660124010
# ╟─8f97df78-ad41-45b7-a169-c1ada2bdb8c9
# ╟─d1133895-2f79-4e40-8310-d88277b0e643
# ╠═cfc519e3-f6e4-4e69-b46f-a3140a293262
# ╠═9d08f7b3-2536-496c-9bbc-93500e78cffd
# ╠═6fee0058-c354-4813-99f9-226f099e58c4
# ╟─6c36c0f0-4eb6-419d-859a-bea6528c3522
# ╟─a2f7dfaa-26e4-41b9-8ab4-6f2123f8dcea
# ╠═8510cc2f-2cc1-4329-8cc2-d9a12441a231
# ╠═8887b407-9604-42a7-bdf6-96deb7b89b0e
# ╟─20619594-5a0f-4f29-9ef3-517641ee7ee9
# ╠═d78b38fe-6adf-4da7-a8d4-64bb2bffd3fd
# ╟─f9cc2e13-a2ca-468c-9e17-a2335972ab82
# ╠═0b2b9265-648f-47f3-ab74-200b568156e0
# ╟─73ad232d-f45f-416f-8ce4-e3b817a8493a
# ╠═a8b14f6b-eade-4510-9e42-c80fb83aa413
# ╟─07e190ab-6094-42be-8570-1b585e643bcc
# ╠═0ed03e10-9897-4fdd-a0d1-18d75eab408b
