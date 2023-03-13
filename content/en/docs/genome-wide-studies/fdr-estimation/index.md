---
categories: ["Statistical significance"]
tags: ["statistical significance", "algorithms"]
title: "False discovery rate estimation"
linkTitle: "False discovery rate estimation"
weight: 2
date: 2023-02-10
description: >
  False discovery rate estimation.
---




{{< alert title="Reference" >}}
[ESL] Trevor Hastie, Robert Tibshirani, and Jerome
Friedman. *The Elements of Statistical Learning (second edition)* (2009).

<https://hastie.su.domains/ElemStatLearn/>\
<https://link.springer.com/book/10.1007%2F978-0-387-84858-7>

Section 18.7
{{< /alert >}}

## Statistical significance

When we analyze genome-wide data, we are usually interested in identifying **statistically significant** features or patterns, that is, patterns we would not expect to see purely by chance. For example:

- Genes differentially expressed between healthy and disease samples, treated and untreated samples, etc.
- Genes coexpressed across tumour samples
- ...

We will use differential expression between two conditions, treated and untreated, as a running example. Of course this applies equally to other types of comparisons, e.g. ER positive vs ER negative breast tumours. 

Statistical significance is expressed in terms of comparing a **null hypothesis** $H_0$ to an **alternative hypothesis** $H_1$. For instance

$$
\begin{aligned}
  H_0 &= \text{treatment has no effect on gene expression}\\\\
  H_1 &= \text{treatment does have an effect on gene expression}
\end{aligned}
$$

To compare gene expression between the two groups, we can for instance do a [two-sample $t$-test](https://en.wikipedia.org/wiki/Student%27s_t-test): compute for each gene

$$
t = \frac{\overline{x_1} - \overline{x_2}}{\text{se}}
$$

where $\overline{x_1}$ and $\overline{x_2}$ are the average expression levels of the gene in each group, and $\text{se}$ is the pooled within-group standard error:

$$
\text{se} = \sqrt{\frac{s_1^2}{n_1}+\frac{s_2^2}{n_1}}
$$

where $s_i^2$ and $n_i$ are the variance and number of samples in group $i$, respectively

$$
\begin{aligned}
  \overline{x_i} &= \frac{1}{n_i}\sum_{j\in G_i} x_{j} \\\\
  s_i^2 &= \frac{1}{(n_i-1)}\sum_{j\in G_i} (x_i - \overline{x_i})^2
\end{aligned}
$$

The $t$-statistic measures the difference in average expression between the two groups, relative to the natural variation within each group. The greater $|t|$, the more likely there is a true difference between the groups. Nevertheless, even if there is no true difference, some variation between the groups will be observed in a finite number of samples due to random sampling noise. The distribution of $t$-values observed between two groups when in fact there is no true difference between them is called the **null distribution**. Statistical significance is usually expressed by the $p$-value, which, for a given value of $t$, expresses the probability to observe a $t$-statistic greater than $t$ (in absolute value) under the null hypothesis:

$$
\begin{aligned}
  p = \text{Pr}(|T|\geq |t| \mid H_0)
\end{aligned}
$$

For some test statistics, the true null distribution to compute the $p$-value is known. Otherwise random permutations can be used: randomly permute the group labels of the samples a large number of times, compute the test statistic for each permutation, and use these values to construct an **empirical null distribution**.

An important property of $p$-values is the following: under the null hypothesis, $p$-values are uniformly distributed between 0 and 1:

$$
\text{Pr}(P\leq p \mid H_0) = p
$$

Conventionally, a $p$-value threshold of 0.05 is often used. This means that the probability that we declare a gene as being affected by the treatment, even though in reality it is not affected, is less than 5%.

In genome-wide studies, we typically measure around 20,000 genes (in humans). From the uniformity of the $p$-value distribution under the null hypothesis, it follows that even if none of the genes are affected by treatment, we would still expect 5% or 1,000 genes to have a $p$-value less than 0.05. Clearly it would be wrong to call any of these genes differentially expressed!

## False discovery rate

The $p$-value is often wrongly interpreted as saying that the probability that our gene is truly affected by treatment is greater than 95% if $p<0.05$. The latter is the probability that the alternative hypothesis is true given that the $p$-value is below a certain threshold, or $\text{Pr}(H_1 \mid P\leq p)$. More commonly, we express this through the **false discrovery rate (FDR)**, that is, the probability that the null hypothesis is true given that the $p$-value is below a certain threshold, $\text{Pr}(H_0 \mid P\leq p)=1-\text{Pr}(H_1 \mid P\leq p)$.

To get a feel for what such probabilities mean, consider the following table:

|                 | Not Significant | Significant | Total |
|-----------------|-----------------|-------------|-------|
| **$H_0$ true**  | $U$             | $V$         | $M_0$ |
| **$H_0$ false** | $T$             | $S$         | $M_1$ |
|**Total**        | $M-R$           | $R$         | $M$   |

Then the false discovery rate is

$$
\text{FDR} = \mathbb{E}\left(\frac{V}{R}\right)
$$

Obviously, we never know the true value of the numbers in the table, and hence FDR must be estimated. There exist several approaches for doing this. We consider two popular ones: the plug-in estimator and the Bayesian approach from Storey's paper.

### Plug-in estimator of the FDR

Consider again the differential expression problem. Compute $t$-statistics for each gene and consider a range of thresholds $C$, either on the absolute $t$-values or on the corresponding $p$-values. For each value of $C$, we know $R_{obs}(C)$, the observed number of significant genes at threshold $C$.

Now generate $K$ random permutations of the group labels, and compute for each permutation $k$, $V_{obs}(C,k)$, the number of significant genes at threshold $C$ in permutation $k$, which by virtue of the randomization, must all correspond to cases where $H_0$ is true. Compute the average over all permutations:

$$
V_{obs}(C) = \frac{1}{K}\sum_{k=1}^K V_{obs}(C,k)
$$

The plug-in estimate of the FDR at threshold $C$ is then defined as

$$
\widehat{\text{FDR}}(C) = \frac{V_{obs}(C)}{R_{obs}(C)}
$$

We can now vary $C$ until we reached a desired FDR value, e.g. 10%, such that we expect no more than 10% of the genes we call differentially expressed to be false positives.

Note that in practice, the number of permutations $K$ need not be very large to reach a stable value for the average. This is because in practice

### Bayesian estimate of the FDR

When we perform a genome-wide test for differential expression, we obtain on the order of $10^4$ test statistic values, with corresponding $p$-values under the null hypothesis. The distribution of these $p$-values, visualized using a histogram, typically looks like this:

*Figure*

See also: [How to interpret a $p$-value histogram](http://varianceexplained.org/statistics/interpreting-pvalue-histogram/)

We know that under the null hypothesis, $p$-values are uniformly distributed. We recognize the uniform distribution in the flat profile of the histogram for unsignificant ($p\to 1$) $p$-values, which with very high probability come from genes unaffected by the treatment. At very small $p$-values, we see a deviation from the uniform distribution, which indicates the presence of truly affected genes, which obsviously don't follow the null distribution.

Hence it seems like we can model the observed $p$-value distribution as a **mixture distribution** with one component corresponding to the null distribution and one component corresponding to the alternative distribution:

$$
f(p) = \pi_0 f_0(p) + (1-\pi_0) f_1(p)
$$

where $f(p)$ is the probability density function (p.d.f.) of the observed $p$-value distribution, $f_0$ and $f_1$ are the p.d.f. of the $p$-value distribution under the null and alternative distribution, respectively, and $\pi_0 = \text{Pr}(H_0)$, the prior probability of a gene being unaffected by treatment (prior as in, before observing any data), . 

As before, we can imagine that our $p$-values are generated by a process that first samples a hidden variable $Z=0$ with probability $\pi_0$ and $Z=1$ with probability $1-\pi_0$, and then samples a $p$-value (or $t$-statistic) from the null or alternative distribution depending on the value of $Z$. We know the null distribution is the uniform distribution. If we knew a parametric form for the alternative distribution, we could apply EM to estimate $\pi_0$ and the parameters of the alternative distribution, and then compute for each gene the *recognition distribution* or *local false discovery rate* (lfdr)

$$
\begin{aligned}
  \text{lfdr}(p)  &= \text{Pr}(H_0 \mid p)
  = \text{Pr}(Z=0 \mid p)
  = \frac{\pi_0 f_0(p)}{f(p)}
\end{aligned}
$$

The word *"local"* refers to the fact that the above expression gives the expected rate of false positive among all genes with $p$-value equal to $p$, as opposed to the *"tail"* FDR estimated in the previous section giving the he expected rate of false positive among all genes with $p$-value less than or equal to $p$. 

In reality, EM is rarely used in this context, because we rarely have a good enough idea about a parametric form for the alternative distribution. Instead a non-parametric approach is used that exploits the knowledge that:

1. $p$-values are uniformly distributed under the null hypothesis, that is

   $$
   f_0(p) = 1
   $$

2. $p$-values close to 1 almost surely come from the null distribution, that is,

   $$
  \begin{aligned}
    \lim_{p\to 1} f_1(p) &= 0 \\\\
    \lim_{p\to 1} f(p) &= \pi_0 \lim_{p\to 1} f_0(p) = \pi_0
  \end{aligned}
  $$

Hence in the observed $p$-value histogram, $\pi_0$ can be estimated from the height of the "flat" region near $p\approx 1$.

We obtain

$$
\text{lfdr}(p) = \frac{\hat{\pi}_0}{f(p)}
$$

## Assignment


{{< alert title="Assignment" >}}


{{< /alert >}}