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

# ╔═╡ 44b0fb80-c13b-490a-a28c-f45fe69ffe97
using DrWatson

# ╔═╡ 0bd0d1d8-0bdb-4755-b3db-b1885de15ae4
# ╠═╡ show_logs = false
quickactivate(@__DIR__)

# ╔═╡ c6664290-bc06-48ae-b586-9ea2f8af5d36
# ╠═╡ show_logs = false
begin
	using DataFrames
	using CSV
	using Statistics
	using StatsBase
	using StatsPlots, LaTeXStrings
	using MLJ
	using MLJScikitLearnInterface
	using MLJModels
	using Bootstrap
	using PlutoUI
end

# ╔═╡ 0b0c3970-de94-11ed-1198-03fb54d95a60
md"# Drug sensitivity prediction using regularized regression
## Setup the environment
"

# ╔═╡ cfbbbabf-dd48-4b70-a6ea-27e5db433f7f
md"## Import data
See the notebook data-preprocessing-CCLE.jl
"

# ╔═╡ e307c340-444e-41d8-9957-84171a383936
fexpr = datadir("exp_pro","CCLE","CCLE-expr.csv");

# ╔═╡ 448a7611-8f82-4560-91a3-401e675d95bb
expr = DataFrame(CSV.File(fexpr));

# ╔═╡ 260cdc0a-e7ee-4c43-816f-05a258b76aa8
fsens =  datadir("exp_pro","CCLE","CCLE-ActArea.csv");

# ╔═╡ 27b21c99-92a9-4a0e-9825-31daa7db8253
actarea = DataFrame(CSV.File(fsens));

# ╔═╡ 5ebf223f-1dac-4250-a613-3ff56853d52a
md"Double-check sample alignment and then drop cell line names"

# ╔═╡ 724eb427-58aa-4f5f-86c9-ad6e6e499cf6
sum(expr.Cell_line .== actarea.Cell_line)

# ╔═╡ cf9a7a42-4a8c-4389-87a2-74634c115508
begin	
	select!(expr, Not(:Cell_line))
	select!(actarea, Not(:Cell_line))
end

# ╔═╡ 9c77c287-1b9a-408f-8b5e-c0c08114f12f
md"## Elastic net regression
Start by choosing the drug sensitivity we want to predict and :
"

# ╔═╡ 26ea5f2f-bf4e-4fdb-8dab-7717e3ee8936
y = actarea."PD-0325901";

# ╔═╡ 38c26212-c25b-4d1f-8412-80af034f58a2
X = Matrix(expr);

# ╔═╡ 1fb45595-087e-412c-822c-4e95762eeca6
md"We will use more conventional ML approaches for hyperparameter tuning and model validation than what was used in the CCLE paper.

Start by splitting into training and testing samples:
"

# ╔═╡ 8ce82b6f-b9a1-46a6-8885-285c7cd15111
(Xtrain, Xtest), (ytrain, ytest) = partition((expr,y), 0.8; shuffle=true, multi=true);

# ╔═╡ 88521c4f-cd8d-484f-a4b6-4e735702fb9d
md"Prefilter genes based on correlation with the outcome:"

# ╔═╡ 3ed5819f-b5f7-413e-a50d-71faa4a9bc0e
cc = map(x -> abs(cor(x,ytrain)), eachcol(Xtrain));

# ╔═╡ 0c1bc2c8-cbc5-4af2-94ba-004dd7272e57
@bind cccut Slider(0:0.05:0.5)

# ╔═╡ ae6af865-a992-4aee-a1d4-c4db8127378a
begin
	histogram(cc, xlabel="Correlation with drug sensitivity", label="")
	vline!([cccut],linewidth=2,label="")
end

# ╔═╡ f2b1a375-236e-429a-b3f9-64dea50f9a2f
md"Reduce predictor data and keep indices of selected genes"

# ╔═╡ 39e157a7-432e-485e-8f42-4152077c52a2
begin
	gene_select = cc.>cccut; #sd.>sdcut .&& cc.>cccut;
	Xtrain_select = Xtrain[:,gene_select]; 
	Xtest_select = Xtest[:,gene_select]; 
end

# ╔═╡ 9d2246c6-8d68-449c-b258-58eb914e46f8
md"Create a standardizer for the training data which we will later apply to the test data."

# ╔═╡ 39ab3efe-46e7-4b40-aba9-2dbac702cf16
# ╠═╡ show_logs = false
begin
	Standardizer = @load Standardizer pkg=MLJModels;
	std_mach = machine(Standardizer(), Xtrain_select);
	fit!(std_mach);
end

# ╔═╡ a63b2bb5-1c40-4dd2-a0c0-2630f65c370a
md"Get the transformed training and testing data. Using MLJ recommended code is extremely slow - why? Instead of using this:
```julia
Xtrain_select_std = MLJ.transform(std_mach,Xtrain_select);
Xtest_select_std = MLJ.transform(std_mach,Xtest_select);
```
manual standardization is much faster:
"

# ╔═╡ 89ed9b99-7125-4ff9-8101-335852282e71
Xtrain_select_std = mapcols(zscore, Xtrain_select);

# ╔═╡ f40d36cb-fe85-4cca-927f-ab38227e8aa4
md"To transform the testing data, we need to use the means and standard deviations from the training data. First make a copy of Xtest\_select with means and standard deviations from the training data added as new rows. Then transform and remove the added rows.
"

# ╔═╡ 92048503-c94e-440f-b046-4fde669ee82f
begin
	Xtest_select_std = Xtest_select;
	prepend!(Xtest_select_std,  mapcols(std, Xtrain_select));
	prepend!(Xtest_select_std,  mapcols(mean, Xtrain_select));
	mapcols!(x -> (x .- x[1])./ x[2], Xtest_select_std);
	deleteat!(Xtest_select_std,1:2);
end

# ╔═╡ e071ecf4-2441-4bf5-8830-9546428e27fa
md"Define an elastic net cross-validation model. After a bit of playing around I didn't see much difference in performance when using multiple values for the L1 ratio parameter, or using the 250 values for the regularization strength paramater as in the CCLE paper. Hence, for simplicity, just use the default L1 ratio of 0.5, with the default 100  automatically determined regularization strength paramater values, and the default 5-fold cross-validation.
"

# ╔═╡ ccdadd48-7776-4937-ba1a-264865d8cf24
# ╠═╡ show_logs = false
begin
	ElNetCV = @load ElasticNetCVRegressor;
	elnetcv_mach = machine(ElNetCV(),Xtrain_select_std,ytrain);
	fit!(elnetcv_mach);
end

# ╔═╡ a2095a59-b776-45d0-8279-2872ba12f318
elnetcv_param = fitted_params(elnetcv_mach);

# ╔═╡ c2529615-388f-4ee5-969a-c172ac46d577
md"Show MSE mean and standard deviation across the folds"

# ╔═╡ fcb84ae1-8b63-4f8d-9e9e-0eb0089169d8
errorline(elnetcv_param.alphas, elnetcv_param.mse_path,
  label = "",
  xlabel = "Regularization strength",
  ylabel = "MSE (mean and std)"
)

# ╔═╡ 03aff45f-206d-44ad-a587-1a8113056f85


# ╔═╡ Cell order:
# ╠═0b0c3970-de94-11ed-1198-03fb54d95a60
# ╠═44b0fb80-c13b-490a-a28c-f45fe69ffe97
# ╠═0bd0d1d8-0bdb-4755-b3db-b1885de15ae4
# ╠═c6664290-bc06-48ae-b586-9ea2f8af5d36
# ╠═cfbbbabf-dd48-4b70-a6ea-27e5db433f7f
# ╠═e307c340-444e-41d8-9957-84171a383936
# ╠═448a7611-8f82-4560-91a3-401e675d95bb
# ╠═260cdc0a-e7ee-4c43-816f-05a258b76aa8
# ╠═27b21c99-92a9-4a0e-9825-31daa7db8253
# ╟─5ebf223f-1dac-4250-a613-3ff56853d52a
# ╠═724eb427-58aa-4f5f-86c9-ad6e6e499cf6
# ╠═cf9a7a42-4a8c-4389-87a2-74634c115508
# ╟─9c77c287-1b9a-408f-8b5e-c0c08114f12f
# ╠═26ea5f2f-bf4e-4fdb-8dab-7717e3ee8936
# ╠═38c26212-c25b-4d1f-8412-80af034f58a2
# ╟─1fb45595-087e-412c-822c-4e95762eeca6
# ╠═8ce82b6f-b9a1-46a6-8885-285c7cd15111
# ╟─88521c4f-cd8d-484f-a4b6-4e735702fb9d
# ╠═3ed5819f-b5f7-413e-a50d-71faa4a9bc0e
# ╟─0c1bc2c8-cbc5-4af2-94ba-004dd7272e57
# ╠═ae6af865-a992-4aee-a1d4-c4db8127378a
# ╟─f2b1a375-236e-429a-b3f9-64dea50f9a2f
# ╠═39e157a7-432e-485e-8f42-4152077c52a2
# ╟─9d2246c6-8d68-449c-b258-58eb914e46f8
# ╠═39ab3efe-46e7-4b40-aba9-2dbac702cf16
# ╟─a63b2bb5-1c40-4dd2-a0c0-2630f65c370a
# ╠═89ed9b99-7125-4ff9-8101-335852282e71
# ╟─f40d36cb-fe85-4cca-927f-ab38227e8aa4
# ╠═92048503-c94e-440f-b046-4fde669ee82f
# ╟─e071ecf4-2441-4bf5-8830-9546428e27fa
# ╠═ccdadd48-7776-4937-ba1a-264865d8cf24
# ╠═a2095a59-b776-45d0-8279-2872ba12f318
# ╟─c2529615-388f-4ee5-969a-c172ac46d577
# ╠═fcb84ae1-8b63-4f8d-9e9e-0eb0089169d8
# ╠═03aff45f-206d-44ad-a587-1a8113056f85
