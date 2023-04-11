---
categories: ["Regularized regression"]
tags: ["regularized regression", "algorithms"]
title: "Elastic net regression tutorial"
linkTitle: "Tutorial"
weight: 3
draft: false
description: >
  Elastic net regression tutorial using CCLE data.
---


## Import packages

```julia
using DrWatson
@quickactivate "Genome-scale ML"

using DataFrames
using CSV
using XLSX
using Statistics
using StatsBase
using StatsPlots, LaTeXStrings
using MLJ
using MLJScikitLearnInterface
using MLJModels
using Bootstrap
```



## Process CCLE data



- CCLE expression data from [GSE36139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36139); use the Series Matrix File GSE36139-GPL15308_series_matrix.txt.
- Drug sensitivities: Supplementary Table 11 from the [original publication](https://www.nature.com/articles/nature11003); use the "activitiy area" (actarea) as a response variable.



### Expression data processing

A series matrix file consists of metadata (all lines beginning with "!") and data (all lines between "!series_matrix_table_begin" and "!series_matrix_table_end"). 

From the metadata, we need the mapping from GEO sample names (GSM...) to cell line names. This information is in the lines "!Sample_title" and "!Sample_geo_accession". I manually removed the metadata and made the sample_title line the header line, and saved to a new file:

```julia
fexpr = datadir("exp_raw","CCLE","GSE36139-GPL15308_series_matrix_nometa.txt");
```



### Expression data



Read the data into a dataframe. Drop genes with missing data. Transpose to have genes as columns, samples as rows. 

```julia
df = DataFrame(CSV.File(fexpr));
dropmissing!(df);
df = permutedims(df,1);
```



### Drug sensitivity data processing

I manually converted the drug sensitivity file to XLSX format:

```julia
fsens = datadir("exp_raw","CCLE","41586_2012_BFnature11003_MOESM90_ESM.xlsx");
```



Read the cell lines column into a vector of strings:

```julia
cl = vec(map(x -> string(x), XLSX.readdata(fsens, "Table S11!B4:B11673")));
```



Read the compounds column into a vector of strings:

```julia
cp = vec(map(x -> string(x), XLSX.readdata(fsens, "Table S11!C4:C11673")));
```



Read the activity area column into a vector of reals:

```julia
aa = vec(map(x -> float(x), XLSX.readdata(fsens, "Table S11!M4:M11673")));
```



Merge into dataframe:

```julia
df2 = DataFrame("Cell_line" => cl, "Compound" => cp, "ActArea" => aa);
```



Create a grouped dataframe, grouped by compound:

```julia
gf = groupby(df2,:Compound);
```



### Merge expression and clinical data



Find common cell lines

```julia
cl_common = intersect(df.Sample_title,unique(df2.Cell_line));
```



Keep only cell lines common in both datasets:

```julia
expr_in_common = map(x -> in(x,cl_common), df.Sample_title);
expr = df[expr_in_common,:];
rename!(expr, "Sample_title" => "Cell_line");

sens_in_common = map(x -> in(x,cl_common), df2.Cell_line);
sens = df2[sens_in_common,:];
```



Create a new dataframe for the activitiy area values that has cell lines as rows and compounds as columns, by joining the entries in the grouped dataframe one-by-one. Use "leftjoin", which will put missing values for cell lines not present for a given compound.

```julia
gf = groupby(sens,:Compound);
actarea = DataFrame("Cell_line" => expr[:,1]);
for k=eachindex(gf)
    leftjoin!(actarea, gf[k][:,[1,3]], on=:Cell_line)
    rename!(actarea, "ActArea" => gf[k].Compound[1])
end
```



### Save data

```julia
fexpr_out = datadir("exp_pro","CCLE","CCLE-expr.csv");
CSV.write(fexpr_out,expr);

fsens_out = datadir("exp_pro","CCLE","CCLE-ActArea.csv");
CSV.write(fsens_out,actarea);
```




## Import data



See data-preprocessing-ccle.jl

```julia
fexpr = datadir("exp_pro","CCLE","CCLE-expr.csv");
expr = DataFrame(CSV.File(fexpr));

fsens =  datadir("exp_pro","CCLE","CCLE-ActArea.csv");
actarea = DataFrame(CSV.File(fsens));
```



Double-check sample alignment and then drop cell line names

```julia
sum(expr.Cell_line .== actarea.Cell_line);
expr = expr[:, Not(:Cell_line)];
actarea = actarea[:, Not(:Cell_line)];
```



## Elastic net regression

Start by choosing the drug sensitivity we want to predict:

```julia
y = actarea."PD-0325901";
X = Matrix(expr);
```



We will use more conventional ML approaches for hyperparameter tuning and model validation than what was used in the CCLE paper.

Start by splitting into training and testing samples

```julia
(Xtrain, Xtest), (ytrain, ytest) = partition((expr,y), 0.8; shuffle=true, multi=true);
```



Prefilter and standardize genes based on training data alone

First find variable genes

```julia
sd = map(x -> std(x), eachcol(Xtrain));
sdcut = 0.7;
```



Show histogram of standard deviations and cutoff line

```julia
histogram(sd, xlabel="Standard deviation", label="")
vline!([sdcut],linewidth=2,label="")
```

{{< figure src="regularized-regression-ccle_7_1.png"  >}}

Now find genes correlated with response

```julia
cc = map(x -> abs(cor(x,ytrain)), eachcol(Xtrain));
cccut = 0.1;
```



Show histogram of correlations and cutoff line

```julia
histogram(cc, xlabel="Correlation with drug sensitivity", label="")
vline!([cccut],linewidth=2,label="")
```

{{< figure src="regularized-regression-ccle_9_1.png"  >}}

Reduce predictor data and keep indices of selected genes

```julia
gene_select = cc.>cccut; #sd.>sdcut .&& cc.>cccut;
Xtrain_select = Xtrain[:,gene_select]; 
Xtest_select = Xtest[:,gene_select];
```



Create a standardizer for the training data which we will later apply to the test data.

```julia
Standardizer = @load Standardizer pkg=MLJModels;
std_mach = machine(Standardizer(), Xtrain_select);
fit!(std_mach);
```

```
import MLJModels ✔
```




Get the transformed training and testing data. This is extremely slow - why?

```julia
Xtrain_select_std = MLJ.transform(std_mach,Xtrain_select);
Xtest_select_std = MLJ.transform(std_mach,Xtest_select);
```


Manual standardization is much faster

```julia
Xtrain_select_std = mapcols(zscore, Xtrain_select);
```



To transform the testing data, we need to use the means and standard deviations from the training data. 

First make a copy of Xtest_select with means and standard deviations from the training data added as new rows. Then transform and remove the added rows.

```julia
Xtest_select_std = Xtest_select;
prepend!(Xtest_select_std,  mapcols(std, Xtrain_select));
prepend!(Xtest_select_std,  mapcols(mean, Xtrain_select));
mapcols!(x -> (x .- x[1])./ x[2], Xtest_select_std);
deleteat!(Xtest_select_std,1:2);
```



Define an elastic net cross-validation model. After a bit of playing around I didn't see much difference in performance when using multiple values for the L1 ratio parameter, or using the 250 values for the regularization strength paramater as in the CCLE paper. Hence, for simplicity, just use the default L1 ratio of 0.5, with the default 100  automatically determined regularization strength paramater values, and the default 5-fold cross-validation.

```julia
ElNetCV = @load ElasticNetCVRegressor pkg=ScikitLearn;
elnetcv_mach = machine(ElNetCV(),Xtrain_select_std,ytrain);
fit!(elnetcv_mach);

elnetcv_param = fitted_params(elnetcv_mach);
```

```
import MLJScikitLearnInterface ✔
```




Show MSE mean and standard deviation across the folds

```julia
errorline(elnetcv_param.alphas, elnetcv_param.mse_path,
  label = "",
  xlabel = "Regularization strength",
  ylabel = "MSE (mean and std)"
)
```

{{< figure src="regularized-regression-ccle_15_1.png"  >}}

Find regularization strength with lowest mean MSE. Note that what is called `alphas`` in the elastic net regressor is what we called λ in the lecture notes (what we called α is called `l1_ratio`).

```julia
λ = elnetcv_param.alphas[argmin(mean(elnetcv_param.mse_path,dims=2))]
```

```
0.02998899770303372
```




Now train a new elastic net model on the training data with this regularization strength to find the final coefficient estimates:

```julia
ElNet = @load ElasticNetRegressor pkg=ScikitLearn;
elnet_mach = machine(ElNet(alpha=λ),Xtrain_select_std,ytrain);
fit!(elnet_mach);
elnet_param = fitted_params(elnet_mach);
```

```
import MLJScikitLearnInterface ✔
```




Test the trained model on the test data.

```julia
ypred = MLJ.predict(elnet_mach,Xtest_select_std);

scatter(ytest,ypred,
  label = "",
  xlabel = "True drug sensitivities",
  ylabel = "Predicted drug sensitivities"
)
```

{{< figure src="regularized-regression-ccle_18_1.png"  >}}

In the CCLE paper, bootstrapping is performed using the optimal regularization strength to determine the robustness of selected features

```julia
B = 200;
ns = length(ytrain);
coef_bootstrap = zeros(length(elnet_param.coef),B);
for b = 1:B
  # Sample training samples with replacement
  bootstrap_samples = sample(1:ns,ns,replace=true);
  coef_bootstrap[:,b] = fitted_params(
    fit!(machine(ElNet(alpha=λ),Xtrain_select_std[bootstrap_samples,:],ytrain[bootstrap_samples]))
  ).coef;
end
```



Compute the bootstrap frequencies

```julia
freq_bootstrap = vec(sum(coef_bootstrap .!= 0, dims=2)) / B;
```


Visualize the bootstrap frequencies

```julia
fcut = 0.6;
histogram(freq_bootstrap, xlabel="Bootstrap frequency", label="")
vline!([fcut],linewidth=2,label="")
```

{{< figure src="regularized-regression-ccle_21_1.png"  >}}

Print gene names and frequencies above cutoff

```julia
genes_bootstrap = DataFrame("Gene" => names(Xtrain_select)[freq_bootstrap .>= fcut], "Frequency" => freq_bootstrap[freq_bootstrap .>= fcut]);
sort!(genes_bootstrap, :Frequency, rev=true)
```


<div><div style = "float: left;"><span>25×2 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "header"><th class = "rowNumber" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">Gene</th><th style = "text-align: left;">Frequency</th></tr><tr class = "subheader headerLastRow"><th class = "rowNumber" style = "font-weight: bold; text-align: right;"></th><th title = "String" style = "text-align: left;">String</th><th title = "Float64" style = "text-align: left;">Float64</th></tr></thead><tbody><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">646600_at</td><td style = "text-align: right;">0.97</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">10253_at</td><td style = "text-align: right;">0.925</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">1594_at</td><td style = "text-align: right;">0.87</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">5728_at</td><td style = "text-align: right;">0.85</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">5</td><td style = "text-align: left;">653483_at</td><td style = "text-align: right;">0.815</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">6</td><td style = "text-align: left;">83897_at</td><td style = "text-align: right;">0.815</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">7</td><td style = "text-align: left;">26031_at</td><td style = "text-align: right;">0.78</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">8</td><td style = "text-align: left;">89801_at</td><td style = "text-align: right;">0.75</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">9</td><td style = "text-align: left;">4023_at</td><td style = "text-align: right;">0.74</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">10</td><td style = "text-align: left;">140825_at</td><td style = "text-align: right;">0.73</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">11</td><td style = "text-align: left;">26330_at</td><td style = "text-align: right;">0.72</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">12</td><td style = "text-align: left;">9315_at</td><td style = "text-align: right;">0.72</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">13</td><td style = "text-align: left;">341405_at</td><td style = "text-align: right;">0.71</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">14</td><td style = "text-align: left;">728819_at</td><td style = "text-align: right;">0.7</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">15</td><td style = "text-align: left;">54823_at</td><td style = "text-align: right;">0.675</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">16</td><td style = "text-align: left;">112616_at</td><td style = "text-align: right;">0.67</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">17</td><td style = "text-align: left;">10326_at</td><td style = "text-align: right;">0.66</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">18</td><td style = "text-align: left;">28659_at</td><td style = "text-align: right;">0.66</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">19</td><td style = "text-align: left;">6847_at</td><td style = "text-align: right;">0.635</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">20</td><td style = "text-align: left;">340602_at</td><td style = "text-align: right;">0.63</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">21</td><td style = "text-align: left;">3950_at</td><td style = "text-align: right;">0.625</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">22</td><td style = "text-align: left;">729296_at</td><td style = "text-align: right;">0.625</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">23</td><td style = "text-align: left;">83890_at</td><td style = "text-align: right;">0.615</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">24</td><td style = "text-align: left;">3781_at</td><td style = "text-align: right;">0.61</td></tr><tr><td class = "rowNumber" style = "font-weight: bold; text-align: right;">25</td><td style = "text-align: left;">84034_at</td><td style = "text-align: right;">0.6</td></tr></tbody></table></div>
