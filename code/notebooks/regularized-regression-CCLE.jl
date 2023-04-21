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
	using GLMNet
	using MLJ
	using Printf
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

# ╔═╡ cf9a7a42-4a8c-4389-87a2-74634c115508
begin	
	sum(expr.Cell_line .== actarea.Cell_line)
	select!(expr, Not(:Cell_line))
	select!(actarea, Not(:Cell_line))
end

# ╔═╡ 9c77c287-1b9a-408f-8b5e-c0c08114f12f
md"## Elastic net regression
Start by choosing the drug sensitivity we want to predict and :
"

# ╔═╡ 26ea5f2f-bf4e-4fdb-8dab-7717e3ee8936
y = actarea."PD-0325901";

# ╔═╡ 1fb45595-087e-412c-822c-4e95762eeca6
md"We will use more conventional ML approaches for hyperparameter tuning and model validation than what was used in the CCLE paper.

Start by splitting into training and testing samples:
"

# ╔═╡ 8ce82b6f-b9a1-46a6-8885-285c7cd15111
(Xtrain, Xtest), (ytrain, ytest) = partition((expr,y), 0.8; shuffle=true, multi=true);

# ╔═╡ 88521c4f-cd8d-484f-a4b6-4e735702fb9d
md"Prefilter genes based on correlation with the outcome. Move the slider to set a correlation coefficient cutoff (default 0.1)"

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
	Xtrain_select = Matrix(Xtrain[:,gene_select]); 
	Xtest_select = Matrix(Xtest[:,gene_select]); 
end

# ╔═╡ 9d2246c6-8d68-449c-b258-58eb914e46f8
md"Create a standardizer for the training data and apply it to both train and test data."

# ╔═╡ 826b0d23-bcf7-433c-82f7-428a39a044e0
begin
	dt = StatsBase.fit(ZScoreTransform, Xtrain_select; dims=1)
	Xtrain_select_std = StatsBase.transform(dt,Xtrain_select)
	Xtest_select_std = StatsBase.transform(dt,Xtest_select)
end

# ╔═╡ e071ecf4-2441-4bf5-8830-9546428e27fa
md"Define an elastic net cross-validation model. After a bit of playing around I didn't see much difference in performance when using multiple values for the L1 ratio parameter ($\alpha$), or using the 250 values for the regularization strength paramater ($\lambda$) as in the CCLE paper. Hence, for simplicity, use a fixed value for the L1 ratio (default 0.5), and the default 100  automatically determined regularization strength paramater values, and the default 10-fold cross-validation. Move the slider to change the value of the L1 ratio from 0 (ridge regression) to 1 (lasso regression).
"

# ╔═╡ 03aff45f-206d-44ad-a587-1a8113056f85
@bind α Slider(0:0.1:1) 

# ╔═╡ 6519f119-9f91-4ca3-b9ec-d9116149cae3
cv = glmnetcv(Xtrain_select_std, ytrain; 
  alpha=α,
  standardize=false,
  )

# ╔═╡ caabc7ea-2f25-4a88-9621-133e384a08d0
md"Show the mean and standard deviation of the loss across the 10 folds, with the regularization strength on log-scale:"

# ╔═╡ 2fe8a865-00c4-431b-a070-916d0d31a070
begin
	plot(cv.lambda,cv.meanloss, 
	  xscale=:log10, 
	  legend=false, 
	  yerror=cv.stdloss,
	  xlabel=L"\lambda",
	  ylabel="loss")
	vline!([lambdamin(cv)])
end

# ╔═╡ a6a903a9-9b06-436b-a656-36c71a282117
md"Now train a new elastic net model on the training data with this regularization strength to find the final coefficient estimates:"

# ╔═╡ 35dbc788-430b-4491-a792-a369fb81b0a5
mdl = glmnet(Xtrain_select_std, ytrain; 
  alpha=α, 
  lambda=[lambdamin(cv)],
  standardize=false)

# ╔═╡ f7092f32-b332-4a61-80c4-46c92e33f7c7
md"Test the trained model on the test data."

# ╔═╡ 469000fe-d028-4e2f-93c2-fff727d373bd
ypred = GLMNet.predict(mdl, Xtest_select_std);

# ╔═╡ fff9698d-06a5-4e75-9ada-734806b66370
begin
	scatter(ytest,ypred,
	  label = "",
	  xlabel = "True drug sensitivities",
	  ylabel = "Predicted drug sensitivities"
	)
	annotate!(0.5,4.5,Printf.@sprintf("r = %1.2f", cor(ytest,ypred)[1]))
end

# ╔═╡ 33258279-671c-4b29-99eb-a0e7b375b8df


# ╔═╡ b284165c-e351-4a7a-980a-aa66b676de8b
md"## Bootstrapping
In the CCLE paper, bootstrapping is performed using the optimal regularization strength to determine the robustness of selected features. The procedure samples training samples with replacement, relearns the model, and computes how frequently a feature (gene) is chosen as predictor in the bootstrapped models."

# ╔═╡ 8aba9a19-8733-408f-bc4f-e58dd415401c
begin
	B = 200;
	ns = length(ytrain);
	betas_bootstrap = zeros(length(mdl.betas),B);
	for b = 1:B
	  # Sample training samples with replacement
	  bootstrap_samples = sample(1:ns,ns,replace=true);
	  betas_bootstrap[:,b] = glmnet(Xtrain_select_std[bootstrap_samples,:], 
	    ytrain[bootstrap_samples]; 
	    alpha=0.5, 
	    lambda=[lambdamin(cv)],
	    standardize=false).betas;
	end
end

# ╔═╡ 3d8ddba6-9f91-4b9f-8743-ae9072d77556
md"Compute the bootstrap frequencies"

# ╔═╡ 637c83a8-419b-40b4-8327-788cc6a8b51b
freq_bootstrap = vec(sum(betas_bootstrap .!= 0, dims=2)) / B;

# ╔═╡ ec3fd049-7037-4d42-94fd-a5e10a442a73
md"Visualize the bootstrap frequencies and choose an appropriate cutoff by moving the slider."

# ╔═╡ a7398194-f661-42ac-b6d2-0bf21e9c11d1
@bind fcut Slider(0:0.1:1)

# ╔═╡ 574d5df0-4ca4-4ed3-b57a-4d0e868556c3
begin
	histogram(freq_bootstrap, xlabel="Bootstrap frequency", label="")
	vline!([fcut],linewidth=2,label="")
end

# ╔═╡ 6e96c038-e3c7-44eb-81ef-bda3741fab55
md"Print gene names and frequencies above cutoff"

# ╔═╡ a9461142-e2bb-49ed-9c21-f405bb25ee3f
names_select = names(Xtrain[:,gene_select]);

# ╔═╡ 85f4cc21-5e09-4fd0-93a9-e3d436ab8879
begin
	genes_bootstrap = DataFrame("Gene" => names_select[freq_bootstrap .>= fcut], "Frequency" => freq_bootstrap[freq_bootstrap .>= fcut]);
	sort!(genes_bootstrap, :Frequency, rev=true)
end

# ╔═╡ Cell order:
# ╟─0b0c3970-de94-11ed-1198-03fb54d95a60
# ╠═44b0fb80-c13b-490a-a28c-f45fe69ffe97
# ╠═0bd0d1d8-0bdb-4755-b3db-b1885de15ae4
# ╠═c6664290-bc06-48ae-b586-9ea2f8af5d36
# ╟─cfbbbabf-dd48-4b70-a6ea-27e5db433f7f
# ╠═e307c340-444e-41d8-9957-84171a383936
# ╠═448a7611-8f82-4560-91a3-401e675d95bb
# ╠═260cdc0a-e7ee-4c43-816f-05a258b76aa8
# ╠═27b21c99-92a9-4a0e-9825-31daa7db8253
# ╟─5ebf223f-1dac-4250-a613-3ff56853d52a
# ╠═cf9a7a42-4a8c-4389-87a2-74634c115508
# ╟─9c77c287-1b9a-408f-8b5e-c0c08114f12f
# ╠═26ea5f2f-bf4e-4fdb-8dab-7717e3ee8936
# ╟─1fb45595-087e-412c-822c-4e95762eeca6
# ╠═8ce82b6f-b9a1-46a6-8885-285c7cd15111
# ╟─88521c4f-cd8d-484f-a4b6-4e735702fb9d
# ╠═3ed5819f-b5f7-413e-a50d-71faa4a9bc0e
# ╟─0c1bc2c8-cbc5-4af2-94ba-004dd7272e57
# ╠═ae6af865-a992-4aee-a1d4-c4db8127378a
# ╟─f2b1a375-236e-429a-b3f9-64dea50f9a2f
# ╠═39e157a7-432e-485e-8f42-4152077c52a2
# ╟─9d2246c6-8d68-449c-b258-58eb914e46f8
# ╠═826b0d23-bcf7-433c-82f7-428a39a044e0
# ╟─e071ecf4-2441-4bf5-8830-9546428e27fa
# ╟─03aff45f-206d-44ad-a587-1a8113056f85
# ╠═6519f119-9f91-4ca3-b9ec-d9116149cae3
# ╟─caabc7ea-2f25-4a88-9621-133e384a08d0
# ╠═2fe8a865-00c4-431b-a070-916d0d31a070
# ╟─a6a903a9-9b06-436b-a656-36c71a282117
# ╠═35dbc788-430b-4491-a792-a369fb81b0a5
# ╠═f7092f32-b332-4a61-80c4-46c92e33f7c7
# ╠═469000fe-d028-4e2f-93c2-fff727d373bd
# ╠═fff9698d-06a5-4e75-9ada-734806b66370
# ╠═33258279-671c-4b29-99eb-a0e7b375b8df
# ╟─b284165c-e351-4a7a-980a-aa66b676de8b
# ╠═8aba9a19-8733-408f-bc4f-e58dd415401c
# ╟─3d8ddba6-9f91-4b9f-8743-ae9072d77556
# ╠═637c83a8-419b-40b4-8327-788cc6a8b51b
# ╟─ec3fd049-7037-4d42-94fd-a5e10a442a73
# ╟─a7398194-f661-42ac-b6d2-0bf21e9c11d1
# ╠═574d5df0-4ca4-4ed3-b57a-4d0e868556c3
# ╟─6e96c038-e3c7-44eb-81ef-bda3741fab55
# ╠═a9461142-e2bb-49ed-9c21-f405bb25ee3f
# ╠═85f4cc21-5e09-4fd0-93a9-e3d436ab8879
