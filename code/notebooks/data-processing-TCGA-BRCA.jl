### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ 97ad8296-6caf-432f-b599-6cb8ba60199c
using DrWatson

# ╔═╡ fe71668a-d75a-4223-b050-85f5586f0604
# ╠═╡ show_logs = false
quickactivate(@__DIR__)

# ╔═╡ 4bed867b-c741-4c62-bdd3-67a423401b88
begin
	using DataFrames
	using CSV
	using XLSX
end

# ╔═╡ eb339710-dd12-11ed-2124-cf73816253ef
md"# Process TCGA BRCA data
## Setup the environment
"

# ╔═╡ 7cdb12c2-02cd-4f95-8f0d-61f3a192dfa8
md"## Read the data
- BRCA expression data available at https://gdc.cancer.gov/about-data/publications/brca_2012, file BRCA.exp.348.med.txt.
- ER status and cancer stage available in Supplementary Table 1 of the paper at https://www.nature.com/articles/nature11412
- I manually converted the clinical data file to XLSX format and have renamed the files as follows:
"

# ╔═╡ e4914491-77eb-4214-9b0e-b3dfe86c956d
fexpr = datadir("exp_raw","TCGA-BRCA","TCGA-BRCA-exp-348-med.txt");

# ╔═╡ 15136b36-cb9b-466a-a38e-a23243e3771e
fclin = datadir("exp_raw","TCGA-BRCA","TCGA-BRCA-Supplementary-Tables-1-4.xlsx");

# ╔═╡ ede73679-a01c-42f4-b8e6-181ffc22a7fb
md"### Expression data

Read the data first into a dataframe and then extract the matrix of values, and names of samples and genes. Fix the following:

- Sample IDs in Supp Table only use 12 chars
- Remove genes with missing data 
"

# ╔═╡ 7481a093-97b6-4b9c-9f47-b82719f83e4d
begin
	df = DataFrame(CSV.File(fexpr));
	dropmissing!(df); # drop genes with missing values
	df = permutedims(df,1);
end

# ╔═╡ 3ec81279-479a-4aff-bcc4-7ab64fe04c91
md"Rename samples:"

# ╔═╡ 21d6df82-2f72-4c08-b07a-fb9b0246a4b6
for i=1:lastindex(df[:,1])
    df[i,1] = df[i,1][1:12]
end

# ╔═╡ 9d164f96-1afa-4da2-b3fa-4649be8752b9
md"### Clinical data"

# ╔═╡ f369ef3f-87c0-40f7-bb63-285896abc925
begin
	xf = XLSX.readxlsx(fclin); # Read xlsx file
	m = xf["SuppTable1"][:][2:end,:]; # Get first sheet, all data, skipping 1st empty row
	m_samples = Vector{String}(m[2:end,1]); # the sample names
	m_er = Vector{String}(m[2:end,4]); # 4th column is ER status
end

# ╔═╡ 67927dbe-dd0d-4230-9612-fb152aeeeee8
md"Convert ER status to number. Keep only samples with positive or negative ER status.
"

# ╔═╡ 40076de8-2d0f-4dbc-b10b-fb17811c5a85
begin
	pn = unique(m_er)[1:2];
	ER = indexin(m_er,pn);
	tf = .!isnothing.(ER);
	ER = 2*(ER[tf].-1).-1 # Set ER to -1 or +1
	df2 = DataFrame("NAME" => m_samples[tf], "ER" => ER);
end

# ╔═╡ dddae96e-40f5-43d0-bc8b-cc4432974891
md"### Merge expression and clinical data"

# ╔═╡ e0ae25b7-4d90-4e9b-a583-e876f1b38024
dat = innerjoin(df2,df, on=:NAME);

# ╔═╡ c993ab7b-e13f-4808-9422-55348f77c9c7
md"## Save data"

# ╔═╡ 86952336-0a80-40e7-b5b5-52c9a6b05a38
fout = datadir("exp_pro","TCGA-BRCA","TCGA-BRCA-ER-exp-348-med.csv");

# ╔═╡ 062baf3b-29b1-4bb2-9cb3-cc3e643cb294
CSV.write(fout,dat);

# ╔═╡ Cell order:
# ╠═eb339710-dd12-11ed-2124-cf73816253ef
# ╠═97ad8296-6caf-432f-b599-6cb8ba60199c
# ╠═fe71668a-d75a-4223-b050-85f5586f0604
# ╠═4bed867b-c741-4c62-bdd3-67a423401b88
# ╠═7cdb12c2-02cd-4f95-8f0d-61f3a192dfa8
# ╠═e4914491-77eb-4214-9b0e-b3dfe86c956d
# ╠═15136b36-cb9b-466a-a38e-a23243e3771e
# ╟─ede73679-a01c-42f4-b8e6-181ffc22a7fb
# ╠═7481a093-97b6-4b9c-9f47-b82719f83e4d
# ╟─3ec81279-479a-4aff-bcc4-7ab64fe04c91
# ╠═21d6df82-2f72-4c08-b07a-fb9b0246a4b6
# ╟─9d164f96-1afa-4da2-b3fa-4649be8752b9
# ╠═f369ef3f-87c0-40f7-bb63-285896abc925
# ╟─67927dbe-dd0d-4230-9612-fb152aeeeee8
# ╠═40076de8-2d0f-4dbc-b10b-fb17811c5a85
# ╠═dddae96e-40f5-43d0-bc8b-cc4432974891
# ╠═e0ae25b7-4d90-4e9b-a583-e876f1b38024
# ╟─c993ab7b-e13f-4808-9422-55348f77c9c7
# ╠═86952336-0a80-40e7-b5b5-52c9a6b05a38
# ╠═062baf3b-29b1-4bb2-9cb3-cc3e643cb294
