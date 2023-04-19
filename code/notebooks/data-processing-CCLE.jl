### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ 2df5510e-02c0-49ca-a47b-b27d2df861e8
using DrWatson

# ╔═╡ 6a03d77b-8c06-4594-9195-099dea1ae792
# ╠═╡ show_logs = false
quickactivate(@__DIR__)

# ╔═╡ 56ae56ff-5835-4d78-899e-be19c359ff7d
begin
	using DataFrames
	using CSV
	using XLSX
end

# ╔═╡ 965744d0-ddee-11ed-2d82-75e91559b905
md"# Process CCLE data
## Setup the environment
"

# ╔═╡ 0c51af6e-b1f5-4c90-95e6-2dc2f06d6d95
md"## Read the data
- CCLE expression data from [GSE36139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36139); use the Series Matrix File GSE36139-GPL15308_series_matrix.txt.
- Drug sensitivities: Supplementary Table 11 from the [original publication](https://www.nature.com/articles/nature11003); use the activitiy area (actarea) as a response variable.
"

# ╔═╡ e9fc18e5-fa1b-4d0d-9d5b-6cb839aceb7e
md"### Expression data
A series matrix file consists of metadata (all lines beginning with \"!\") and data (all lines between \"!series\_matrix\_table\_begin\" and \"!series\_matrix\_table\_end\").

From the metadata, we need the mapping from GEO sample names (GSM...) to cell line names. This information is in the lines \"!Sample\_title\" and \"!Sample\_geo\_accession\". I manually removed the metadata and made the sample\_title line the header line, and saved to a new file:
"

# ╔═╡ 2f4a7a29-e626-4c35-9245-ffaf3449be26
fexpr = datadir("exp_raw","CCLE","GSE36139-GPL15308_series_matrix_nometa.txt");

# ╔═╡ 4659b47a-4089-40a1-86a8-0c2e7ccc092f
md"Read the data into a dataframe. Drop genes with missing data. Extract gene IDs. Transpose to have genes as columns, samples as rows."

# ╔═╡ 60bd64a1-b5ef-439f-8d70-5805315a871e
begin
	df = DataFrame(CSV.File(fexpr));
	dropmissing!(df);
	df = permutedims(df,1);
	rename!(df,map(x -> x[1: findfirst('_',x)-1], names(df)))
end

# ╔═╡ cd36ee87-0c92-4ce1-bd66-8e497afc7d53
md"### Drug sensitivity data
I manually converted the drug sensitivity file to XLSX format:
"

# ╔═╡ a2489ee9-30f4-4ee2-92bb-85e0cbc2fce8
fsens = datadir("exp_raw","CCLE","41586_2012_BFnature11003_MOESM90_ESM.xlsx"); 

# ╔═╡ 3e7c94fb-a2ba-4114-bd61-91c1cb4da274
md"Read the cell lines column into a vector of strings:"

# ╔═╡ 04665b4d-08f9-4350-ae85-5ba45e4dfb70
cl = vec(map(x -> string(x), XLSX.readdata(fsens, "Table S11!B4:B11673")));

# ╔═╡ 10d19a01-539c-4508-adfa-ee410b34fa74
md"Read the compounds column into a vector of strings:"

# ╔═╡ a2b346d0-18cc-4007-8443-0c4e6c71c530
cp = vec(map(x -> string(x), XLSX.readdata(fsens, "Table S11!C4:C11673")));

# ╔═╡ 18bf54c6-642b-47c4-a849-6669a865b888
md"Read the activity area column into a vector of reals:"

# ╔═╡ 26cb5388-a697-489c-8eec-a9298b6067d3
aa = vec(map(x -> float(x), XLSX.readdata(fsens, "Table S11!M4:M11673")));

# ╔═╡ 3714d7dc-b156-41a2-954e-61e96b574f23
md"Merge into dataframe:"

# ╔═╡ ba697dd6-352d-4690-83ed-d2f058724a8a
df2 = DataFrame("Cell_line" => cl, "Compound" => cp, "ActArea" => aa);

# ╔═╡ fa21f44e-1de6-4a8b-b9d5-20ec9c181172
md"### Merge expression and clinical data
Find common cell lines
"

# ╔═╡ f1ca06f5-230f-4a5c-a868-205bdfbbcafe
cl_common = intersect(df.Sample,unique(df2.Cell_line));

# ╔═╡ 525d598e-52d9-4fa8-ac02-c2818cf7dc70
md"Keep only cell lines common in both datasets:"

# ╔═╡ fcfcdb5e-ad33-4895-b404-2b25b0b38d15
begin
	expr_in_common = map(x -> in(x,cl_common), df.Sample);
	expr = df[expr_in_common,:];
	rename!(expr, "Sample" => "Cell_line");
end

# ╔═╡ c83465d7-090c-4439-9e68-3de8abde3169
begin
	sens_in_common = map(x -> in(x,cl_common), df2.Cell_line);
	sens = df2[sens_in_common,:];
end

# ╔═╡ f156151a-62b6-4cf7-acf9-5162358ec3be
md"Create a new dataframe for the activitiy area values that has cell lines as rows and compounds as columns, by joining the entries in the grouped dataframe one-by-one. Use \"leftjoin\", which will put missing values for cell lines not present for a given compound."

# ╔═╡ ba90046b-3283-451f-a21d-30d86be1d1f7
begin
	gf = groupby(sens,:Compound);
	actarea = DataFrame("Cell_line" => expr[:,1]);
	for k=eachindex(gf)
	    leftjoin!(actarea, gf[k][:,[1,3]], on=:Cell_line)
	    rename!(actarea, "ActArea" => gf[k].Compound[1])
	end
end

# ╔═╡ 13fe29e9-475e-42d6-ae01-58f2286d4b90
md"### Save data
"

# ╔═╡ d7b4ee64-7207-40c5-b8ae-ed6c72e441fb
fexpr_out = datadir("exp_pro","CCLE","CCLE-expr.csv");

# ╔═╡ ae83fe97-d170-44c3-bc05-140d88a80e43
CSV.write(fexpr_out,expr);

# ╔═╡ 1c0e0e6b-cd37-4cf2-84de-7d599ba03fbc
fsens_out = datadir("exp_pro","CCLE","CCLE-ActArea.csv");

# ╔═╡ 4c477b3d-97ba-418c-8bc9-d1cfbdc04d35
CSV.write(fsens_out,actarea);

# ╔═╡ Cell order:
# ╠═965744d0-ddee-11ed-2d82-75e91559b905
# ╠═2df5510e-02c0-49ca-a47b-b27d2df861e8
# ╠═6a03d77b-8c06-4594-9195-099dea1ae792
# ╠═56ae56ff-5835-4d78-899e-be19c359ff7d
# ╟─0c51af6e-b1f5-4c90-95e6-2dc2f06d6d95
# ╟─e9fc18e5-fa1b-4d0d-9d5b-6cb839aceb7e
# ╠═2f4a7a29-e626-4c35-9245-ffaf3449be26
# ╟─4659b47a-4089-40a1-86a8-0c2e7ccc092f
# ╠═60bd64a1-b5ef-439f-8d70-5805315a871e
# ╟─cd36ee87-0c92-4ce1-bd66-8e497afc7d53
# ╠═a2489ee9-30f4-4ee2-92bb-85e0cbc2fce8
# ╟─3e7c94fb-a2ba-4114-bd61-91c1cb4da274
# ╠═04665b4d-08f9-4350-ae85-5ba45e4dfb70
# ╟─10d19a01-539c-4508-adfa-ee410b34fa74
# ╠═a2b346d0-18cc-4007-8443-0c4e6c71c530
# ╟─18bf54c6-642b-47c4-a849-6669a865b888
# ╠═26cb5388-a697-489c-8eec-a9298b6067d3
# ╟─3714d7dc-b156-41a2-954e-61e96b574f23
# ╠═ba697dd6-352d-4690-83ed-d2f058724a8a
# ╟─fa21f44e-1de6-4a8b-b9d5-20ec9c181172
# ╠═f1ca06f5-230f-4a5c-a868-205bdfbbcafe
# ╟─525d598e-52d9-4fa8-ac02-c2818cf7dc70
# ╠═fcfcdb5e-ad33-4895-b404-2b25b0b38d15
# ╠═c83465d7-090c-4439-9e68-3de8abde3169
# ╟─f156151a-62b6-4cf7-acf9-5162358ec3be
# ╠═ba90046b-3283-451f-a21d-30d86be1d1f7
# ╟─13fe29e9-475e-42d6-ae01-58f2286d4b90
# ╠═d7b4ee64-7207-40c5-b8ae-ed6c72e441fb
# ╠═ae83fe97-d170-44c3-bc05-140d88a80e43
# ╠═1c0e0e6b-cd37-4cf2-84de-7d599ba03fbc
# ╠═4c477b3d-97ba-418c-8bc9-d1cfbdc04d35
