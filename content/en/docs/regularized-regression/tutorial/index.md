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

```julia
using DrWatson
@quickactivate "Genome-scale ML"
```


## Import packages

```julia
using DataFrames
using CSV
using XLSX
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