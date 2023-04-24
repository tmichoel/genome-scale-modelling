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
	using Arrow
end

# ╔═╡ 5fb9f4f2-df55-11ed-067c-09eeaf9c6bf9
md"# Process Mouse brain scRNA data
## Setup the environment
"

# ╔═╡ e0c0c65d-618a-4e1d-a024-ed5674832a0e
md"## Convert scRNA csv files to arrow format
"

# ╔═╡ f4f22162-72cf-4f19-8828-3eace1bb7c2d
fALM_in = datadir("exp_raw","mouse-brain-single-cell","mouse_ALM_gene_expression_matrices_2018-06-14", "mouse_ALM_2018-06-14_exon-matrix.csv");

# ╔═╡ 89e512a4-fb2d-46d7-9804-df73805da107
fALM_out = datadir("exp_pro","mouse-brain-single-cell", "mouse_ALM_2018-06-14_exon-matrix.arrow");

# ╔═╡ cedf47e7-96ae-4e91-92bd-b2f55bf4d8e0
if !isfile(fALM_out)
	Arrow.write(fALM_out, CSV.File(fALM_in));
end

# ╔═╡ 34b3fc77-5523-4008-aed8-03ba7fcfe115
fVIS_in = datadir("exp_raw","mouse-brain-single-cell","mouse_VISp_gene_expression_matrices_2018-06-14", "mouse_VISp_2018-06-14_exon-matrix.csv");

# ╔═╡ 7826aef0-fda3-4a9f-b556-d4fac9e0ed03
fVIS_out = datadir("exp_pro","mouse-brain-single-cell", "mouse_VISp_2018-06-14_exon-matrix.arrow");

# ╔═╡ 4280ad62-82d5-4551-b502-9bb06af4bf96
if !isfile(fVIS_out)
	Arrow.write(fVIS_out, CSV.File(fVIS_in));
end

# ╔═╡ 568f73e1-f459-4d7b-820a-6731c12691ec
md"## Read data

### scRNA data"

# ╔═╡ 5c38ec32-0265-4bc5-be16-48b3d899b875
fALM = datadir("exp_pro","mouse-brain-single-cell", "mouse_ALM_2018-06-14_exon-matrix.arrow");

# ╔═╡ 330e3bc0-7fea-4fc6-8397-c0ab7304d661
dfALM = DataFrame(Arrow.Table(fALM))

# ╔═╡ 702db3b6-79c3-4fb6-851c-2cf9047ff081
fVIS = datadir("exp_pro","mouse-brain-single-cell", "mouse_VISp_2018-06-14_exon-matrix.arrow");

# ╔═╡ 3aa88c4d-76cd-459e-94ce-5e448c69bfb9
dfVIS = DataFrame(Arrow.Table(fVIS))

# ╔═╡ d34af095-5821-4812-b786-2cb920e6ac60
md"Merge data:"

# ╔═╡ bf140594-caae-4f97-85ae-f959a634cc13
df = innerjoin(dfALM,dfVIS,on=:Column1)

# ╔═╡ 7256041f-74a4-495c-a0e5-9f73bf10cfb1
# ╠═╡ disabled = true
#=╠═╡
dfVIS = DataFrame(CSV.File(fVIS))
  ╠═╡ =#

# ╔═╡ 5ecc65c3-8b59-4bae-a87e-b2e7164f4b0e
md"### Cluster label data"

# ╔═╡ adae9508-8cb2-44d1-babd-f17c728a7af1
fclust = datadir("exp_raw","mouse-brain-single-cell","tasic-sample_heatmap_plot_data.csv");

# ╔═╡ a570b65f-c918-47e7-a4df-a7cf64fd7c41
dfClust = DataFrame(CSV.File(fclust))

# ╔═╡ d744d3a1-ea53-4606-bf70-808872fc1234
md"Find cell that have a cluster label:"

# ╔═╡ d07682b7-f317-48fd-af15-9c86a8b3b4d1
begin
	tf = .!isnothing.(indexin(names(df),dfClust.sample_name))
	tf[1] = true; # keep the gene IDs
end

# ╔═╡ de329d2b-b611-470d-a564-9d8238d67de4
select!(df, findall(tf))

# ╔═╡ ab5076b6-fcc8-46dc-b941-0313f59a2035
df.Column1 = string.(df.Column1)

# ╔═╡ 224ae815-e34e-48bd-891a-d7b5c8c48f17
# ╠═╡ disabled = true
#=╠═╡
df2 = permutedims(df,1)
  ╠═╡ =#

# ╔═╡ 1a0aace0-3e04-4293-b05a-354154c302f1
# ╠═╡ disabled = true
#=╠═╡
size(df2)
  ╠═╡ =#

# ╔═╡ 7f0e001d-9f08-4f85-a4a5-9b9140435a13
md"### Save data
"

# ╔═╡ aaf6ee66-e5b1-4531-987d-53f9b318cd3d
fexpr = datadir("exp_pro","mouse-brain-single-cell","mouse_ALM_VISp_gene_expression_by_genes.arrow");

# ╔═╡ 7530a342-e562-4f03-a525-eb9dc825bccf
# ╠═╡ disabled = true
#=╠═╡
Arrow.write(fexpr,df2)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─5fb9f4f2-df55-11ed-067c-09eeaf9c6bf9
# ╠═4e7999ab-48e4-480a-8af8-511f7c1d4bd2
# ╠═df94308a-5bfb-4247-a2ea-746372200e4b
# ╠═3525f741-008a-4b10-aaba-9f2374b2bc9b
# ╟─e0c0c65d-618a-4e1d-a024-ed5674832a0e
# ╠═f4f22162-72cf-4f19-8828-3eace1bb7c2d
# ╠═89e512a4-fb2d-46d7-9804-df73805da107
# ╠═cedf47e7-96ae-4e91-92bd-b2f55bf4d8e0
# ╠═34b3fc77-5523-4008-aed8-03ba7fcfe115
# ╠═7826aef0-fda3-4a9f-b556-d4fac9e0ed03
# ╠═4280ad62-82d5-4551-b502-9bb06af4bf96
# ╠═568f73e1-f459-4d7b-820a-6731c12691ec
# ╠═5c38ec32-0265-4bc5-be16-48b3d899b875
# ╠═330e3bc0-7fea-4fc6-8397-c0ab7304d661
# ╠═702db3b6-79c3-4fb6-851c-2cf9047ff081
# ╠═3aa88c4d-76cd-459e-94ce-5e448c69bfb9
# ╟─d34af095-5821-4812-b786-2cb920e6ac60
# ╠═bf140594-caae-4f97-85ae-f959a634cc13
# ╟─7256041f-74a4-495c-a0e5-9f73bf10cfb1
# ╟─5ecc65c3-8b59-4bae-a87e-b2e7164f4b0e
# ╠═adae9508-8cb2-44d1-babd-f17c728a7af1
# ╠═a570b65f-c918-47e7-a4df-a7cf64fd7c41
# ╟─d744d3a1-ea53-4606-bf70-808872fc1234
# ╠═d07682b7-f317-48fd-af15-9c86a8b3b4d1
# ╠═de329d2b-b611-470d-a564-9d8238d67de4
# ╠═ab5076b6-fcc8-46dc-b941-0313f59a2035
# ╠═224ae815-e34e-48bd-891a-d7b5c8c48f17
# ╠═1a0aace0-3e04-4293-b05a-354154c302f1
# ╠═7f0e001d-9f08-4f85-a4a5-9b9140435a13
# ╠═aaf6ee66-e5b1-4531-987d-53f9b318cd3d
# ╠═7530a342-e562-4f03-a525-eb9dc825bccf
