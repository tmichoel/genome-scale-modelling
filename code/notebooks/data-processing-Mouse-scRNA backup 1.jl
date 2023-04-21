### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ 4e7999ab-48e4-480a-8af8-511f7c1d4bd2
using DrWatson

# ╔═╡ df94308a-5bfb-4247-a2ea-746372200e4b
# ╠═╡ show_logs = false
quickactivate(@__DIR__)

# ╔═╡ 3525f741-008a-4b10-aaba-9f2374b2bc9b
begin
	using DataFrames
	using CSV
	using XLSX
end

# ╔═╡ 5fb9f4f2-df55-11ed-067c-09eeaf9c6bf9
md"# Process Mouse brain scRNA data
## Setup the environment
"

# ╔═╡ e0c0c65d-618a-4e1d-a024-ed5674832a0e
md"## Read the data
- 
"

# ╔═╡ bcaa65e3-2cbf-40c8-9971-df0677985412
### scRNA expression data

# ╔═╡ f4f22162-72cf-4f19-8828-3eace1bb7c2d
fALM = datadir("exp_raw","mouse-brain-single-cell","mouse_ALM_gene_expression_matrices_2018-06-14", "mouse_ALM_2018-06-14_exon-matrix.csv");

# ╔═╡ 8d3bd2f3-7146-4a0b-a58c-f906b9dbe78f
dfALM = DataFrame(CSV.File(fALM))

# ╔═╡ 34b3fc77-5523-4008-aed8-03ba7fcfe115
fVIS = datadir("exp_raw","mouse-brain-single-cell","mouse_VISp_gene_expression_matrices_2018-06-14", "mouse_VISp_2018-06-14_exon-matrix.csv");

# ╔═╡ 7256041f-74a4-495c-a0e5-9f73bf10cfb1
dfVIS = DataFrame(CSV.File(fVIS))

# ╔═╡ d34af095-5821-4812-b786-2cb920e6ac60
md"Merge data:"

# ╔═╡ bf140594-caae-4f97-85ae-f959a634cc13
df = innerjoin(dfALM,dfVIS,on=:Column1)

# ╔═╡ 5ecc65c3-8b59-4bae-a87e-b2e7164f4b0e
### Cluster label data

# ╔═╡ adae9508-8cb2-44d1-babd-f17c728a7af1
fclust = datadir("exp-raw","mouse-brain-single-cell","tasic-sample_heatmap_plot_data.csv")

# ╔═╡ a570b65f-c918-47e7-a4df-a7cf64fd7c41
dfClust = DataFrame(CSV.File(fClust))

# ╔═╡ d744d3a1-ea53-4606-bf70-808872fc1234
md"Find cell that have a cluster label:"

# ╔═╡ d07682b7-f317-48fd-af15-9c86a8b3b4d1
tf = .!isnothing.(indexin(cellID,dfClust.sample_name))

# ╔═╡ de329d2b-b611-470d-a564-9d8238d67de4
select!(df, findall(tf))

# ╔═╡ 7f0e001d-9f08-4f85-a4a5-9b9140435a13
md"### Save data
"

# ╔═╡ aaf6ee66-e5b1-4531-987d-53f9b318cd3d
fexpr_out = datadir("exp_pro","mouse-brain-single-cell","mouse_ALM_VISp_gene_expression.csv");

# ╔═╡ Cell order:
# ╟─5fb9f4f2-df55-11ed-067c-09eeaf9c6bf9
# ╠═4e7999ab-48e4-480a-8af8-511f7c1d4bd2
# ╠═df94308a-5bfb-4247-a2ea-746372200e4b
# ╠═3525f741-008a-4b10-aaba-9f2374b2bc9b
# ╟─e0c0c65d-618a-4e1d-a024-ed5674832a0e
# ╠═bcaa65e3-2cbf-40c8-9971-df0677985412
# ╠═f4f22162-72cf-4f19-8828-3eace1bb7c2d
# ╠═8d3bd2f3-7146-4a0b-a58c-f906b9dbe78f
# ╠═34b3fc77-5523-4008-aed8-03ba7fcfe115
# ╠═7256041f-74a4-495c-a0e5-9f73bf10cfb1
# ╠═d34af095-5821-4812-b786-2cb920e6ac60
# ╠═bf140594-caae-4f97-85ae-f959a634cc13
# ╠═5ecc65c3-8b59-4bae-a87e-b2e7164f4b0e
# ╠═adae9508-8cb2-44d1-babd-f17c728a7af1
# ╠═a570b65f-c918-47e7-a4df-a7cf64fd7c41
# ╠═d744d3a1-ea53-4606-bf70-808872fc1234
# ╠═d07682b7-f317-48fd-af15-9c86a8b3b4d1
# ╠═de329d2b-b611-470d-a564-9d8238d67de4
# ╟─7f0e001d-9f08-4f85-a4a5-9b9140435a13
# ╠═aaf6ee66-e5b1-4531-987d-53f9b318cd3d
# ╠═9a7760d7-1669-4e8d-9f28-64ac97a0408f
