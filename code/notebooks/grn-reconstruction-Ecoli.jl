### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ d189778d-0919-49b8-8e7f-2920184c9936
using DrWatson

# ╔═╡ f6d998a6-5b81-4d29-8f8c-1f2b64f5612e
# ╠═╡ show_logs = false
quickactivate(@__DIR__)

# ╔═╡ ae8732af-889e-4881-922b-e00eb8703d26
begin
	using DataFrames
	using CSV
	using Arrow
	using SparseArrays
	using Statistics
	using StatsBase
	using LinearAlgebra
	using StatsPlots
end

# ╔═╡ eff4b1f0-edc3-11ed-3721-29ddcd0941db
md"# Reconstruction of the E. coli GRN
## Setup the environment
"

# ╔═╡ a398c501-6de7-4bee-8f38-e07301db6408
md"## Read the data
### Expression data

- URL: http://m3d.mssm.edu/norm/
- Download: E_coli_v4_Build_6.tar.gz
- Extract: avg\_E\_coli\_v4\_Build\_6\_exps466probes4297.tab
- Truncate gene names at first \"_\"
"

# ╔═╡ b6548f9b-4fc4-4f14-af11-5c3fb51882a7
fexpr = datadir("exp_raw","ecoli-grn","E_coli_v4_Build_6", "avg_E_coli_v4_Build_6_exps466probes4297.tab");

# ╔═╡ a809638e-ee6b-4a9c-af2f-c2b28df9b289
begin
	dfexpr = permutedims(DataFrame(CSV.File(fexpr)),1);
	rename!(x -> lowercase(split(x,"_")[1]), dfexpr);
	select!(dfexpr, Not(1));
end

# ╔═╡ c83237c8-b2c7-4e27-8014-3b4d74dd6b74
md"### RegulongDB network of known TF - gene interactions
- URL: https://regulondb.ccg.unam.mx/menu/download/datasets/index.jsp
- Download: Regulatory Network Interactions > TF - gene interactions
- Extract file: network\_tf\_gene.txt
- Remove all comment lines (start with #) and rename file to regulonDB\_network\_tf\_gene.txt
- Extract 2nd (TFs) and 4th (target genes) columns, change TF names to start with lowercase
"

# ╔═╡ f4725470-b451-4c4e-aaff-73d9ac762ed9
fnet = datadir("exp_raw","ecoli-grn","regulonDB_network_tf_gene.txt");

# ╔═╡ 27c85097-d80c-46b4-b2b6-09271fdaf887
begin
	dfnet = DataFrame(CSV.File(fnet, header=false))
	select!(dfnet, :Column2 => (x -> lowercase.(x)), :Column4 => (x -> lowercase.(x)))
	rename!(dfnet,[1 => :TF, 2 => :target])
	unique!(dfnet)
end

# ╔═╡ 16ef3de3-bf98-43a0-ad23-09044b4eb109
md"## Align the datasets
Find genes common in both datasets
"

# ╔═╡ 978fc217-cbf2-473b-b784-c42d6121b13b
genes_net = intersect(names(dfexpr), union(dfnet.:TF,dfnet.:target))

# ╔═╡ 8fc472a9-f7ed-451d-9c5f-b7f0ad99090d
df = select(dfexpr, genes_net)

# ╔═╡ d4145fc0-f529-49fb-bab4-117ccc202c09
begin
	iTF = indexin(dfnet[:,1],genes_net)
	itgt = indexin(dfnet[:,2],genes_net)
	rws = (.!isnothing.(iTF) .& .!isnothing.(itgt)) .& (iTF != itgt)
end

# ╔═╡ b45389fa-b346-4b5f-833c-28749907cf2d
refnet = sparse(iTF[rws],itgt[rws],ones(sum(rws)),length(genes_net),length(genes_net))

# ╔═╡ d9640052-9618-4eee-8665-f09854461c44
md"## Predict networks
### Pairwise absolute correlation and Z-score networks

- Compute pairwise correlations
- Set diagonal to zero
- Compute Z-scores
- Filter TF rows
"

# ╔═╡ 42e7aa16-1bbf-4352-af57-eebe47c18ceb
begin
	C = abs.(cor(Matrix(df)))
	C -= diagm(diag(C))
end

# ╔═╡ 0019fd80-d041-4f50-b151-d25799dc75c3
begin
	dt1 = fit(ZScoreTransform, C, dims=1)
	Z1 = StatsBase.transform(dt1, C).^2
	dt2 = fit(ZScoreTransform, C, dims=2)
	Z2 = StatsBase.transform(dt2, C).^2
	Z = sqrt.(Z1+Z2)
	Z -= diagm(diag(Z))
end

# ╔═╡ 0b6ef026-ecec-4a86-add6-51edad70a960
begin
	tf = vec(sum(refnet,dims=2).>0)
	Cnet = vec(C[tf,:])
	Znet = vec(Z[tf,:])
	ref = vec(refnet[tf,:])
end

# ╔═╡ f85c08fb-9551-4542-b8aa-1580c97fc7fb
md"### Random forest networks
- Standardize expression data for all genes
- Train random forest regressor for each gene using TFs as predictive features
- Extract feature importances
"

# ╔═╡ cab9d2fc-0555-4cb1-aee1-ae892a1a2bc0
md"## Measure performance
- Calculate recall and precision at a number of thresholds
- Plot precision vs recall curve
"

# ╔═╡ c32e778d-c18e-4780-948c-ed3756086a17
begin
	nTPtot = sum(ref)
	nperf = 500
	rngC = range(0.0,maximum(Cnet),nperf)
	perfC = zeros(length(rngC),2)
	rngZ = range(0.0,maximum(Znet),nperf)
	perfZ = zeros(length(rngZ),2)
	for i = 1:nperf
		tpC = sum(ref[Cnet.>=rngC[i]])
    	perfC[i,1] = tpC/nTPtot # recall
    	perfC[i,2] = tpC/sum(Cnet.>=rngC[i])
    
    	tpZ = sum(ref[Znet.>=rngZ[i]])
    	perfZ[i,1] = tpZ/nTPtot # recall
    	perfZ[i,2] = tpZ/sum(Znet.>=rngZ[i])
	end
end

# ╔═╡ 387d941e-b33f-4742-aa91-ff2dfc2cdc7f
begin
	plot(perfC[:,1],perfC[:,2], linewidth=2, label="Correlation")
	plot!(perfZ[:,1],perfZ[:,2], linewidth=2, label="Z-score")
	xlims!(0, 0.2)
	xlabel!("Recall")
	ylabel!("Precision")
end

# ╔═╡ Cell order:
# ╟─eff4b1f0-edc3-11ed-3721-29ddcd0941db
# ╠═d189778d-0919-49b8-8e7f-2920184c9936
# ╠═f6d998a6-5b81-4d29-8f8c-1f2b64f5612e
# ╠═ae8732af-889e-4881-922b-e00eb8703d26
# ╟─a398c501-6de7-4bee-8f38-e07301db6408
# ╠═b6548f9b-4fc4-4f14-af11-5c3fb51882a7
# ╠═a809638e-ee6b-4a9c-af2f-c2b28df9b289
# ╟─c83237c8-b2c7-4e27-8014-3b4d74dd6b74
# ╠═f4725470-b451-4c4e-aaff-73d9ac762ed9
# ╠═27c85097-d80c-46b4-b2b6-09271fdaf887
# ╟─16ef3de3-bf98-43a0-ad23-09044b4eb109
# ╠═978fc217-cbf2-473b-b784-c42d6121b13b
# ╠═8fc472a9-f7ed-451d-9c5f-b7f0ad99090d
# ╠═d4145fc0-f529-49fb-bab4-117ccc202c09
# ╠═b45389fa-b346-4b5f-833c-28749907cf2d
# ╟─d9640052-9618-4eee-8665-f09854461c44
# ╠═42e7aa16-1bbf-4352-af57-eebe47c18ceb
# ╠═0019fd80-d041-4f50-b151-d25799dc75c3
# ╠═0b6ef026-ecec-4a86-add6-51edad70a960
# ╟─f85c08fb-9551-4542-b8aa-1580c97fc7fb
# ╟─cab9d2fc-0555-4cb1-aee1-ae892a1a2bc0
# ╠═c32e778d-c18e-4780-948c-ed3756086a17
# ╠═387d941e-b33f-4742-aa91-ff2dfc2cdc7f
