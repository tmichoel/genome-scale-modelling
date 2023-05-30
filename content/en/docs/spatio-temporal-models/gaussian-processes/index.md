---
categories: ["Spatio-temporal models"]
tags: ["spatio-temporal models", "algorithms"]
title: "Variance component models and Gaussian processes"
linkTitle: "Gaussian processes"
weight: 2
date: 2023-02-10
description: >
  Variance component models and Gaussian processes.
---


{{< alert title="Reference" >}}
Christopher Bishop. [*Pattern Recognition and Machine Learning*](https://www.microsoft.com/en-us/research/publication/pattern-recognition-machine-learning/) (2006). Chapter 6
{{< /alert >}}


## Linear regression, from fixed to random coefficients

In a standard linear regression model for an outcome $Y$ on one or more predictors $X_1,\dots,X_p$

$$
\begin{aligned}
  Y &= \sum_i\beta_i X_i + \epsilon = X^T\beta + \epsilon, \\qquad
  \epsilon \sim \mathcal{N}(0,\sigma_e^2)
\end{aligned}
$$

the effects $\beta_i$ are treated as **fixed**, such that the distribution of $Y$ is:

$$
p(Y\mid X) = \mathcal{N}(X^T\beta,\sigma_e^2),
$$

that is, in a **fixed effect model**, the predictors affect the average or expected value of $Y$. In a fixed effect model, we always assume that the data to fit the parameters are obtained from **independent samples** of an underlying model.

If our data consists of **non-independent samples** (e.g. due to temporal, spatial or familial relations among the observations), we can use a **random effect model** to explain these correlations as a function of covariates $X$. As [before](../../regularized-regression/ridge-lasso-elnet/), we collect the data for the inputs and ouptput in a $N\times p$ matrix $\mathbf{X}$ and $N$-vector $\mathbf{y}$, respectively,
$$\begin{aligned}
    \mathbf{X}= (x_1,x_2,\dots,x_p) = \begin{pmatrix}
      x_{11} & x_{12} & \dots & x_{1p}\\\
      x_{21} & x_{22} & \dots & x_{2p}\\\
      \vdots & \vdots & & \vdots\\\
      x_{N1} & x_{N2} & \dots & x_{Np}
    \end{pmatrix} &&
    \mathbf{y}= \begin{pmatrix}
      y_1\\\
      y_2\\\
      \vdots\\\
      y_N
    \end{pmatrix}
  \end{aligned}$$

We model the data as

$$
\begin{aligned}
  \mathbf{y} &= \mathbf{X} \alpha + \epsilon\\
\end{aligned}
$$

where $\alpha = (\alpha_1,\dots,\alpha_p)^T$ are random effects with multi-variate normal distribution

$$
\alpha \sim \mathcal{N}(0,\sigma_a^2 I_p)
$$

and $\epsilon = (\epsilon_1,\dots,\epsilon_N)$ are independent error terms

$$
\epsilon \sim \mathcal{N}(0,\sigma_e^2 I_N).
$$

The random effects are assumed to be independent or the errors. We assumed the random effects are mutually independent, but this can easily be generalized.

Using [properties of the multivariate normal distribution](https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Affine_transformation), the distribution of $\mathbf{y}$ can be seen to be:

$$
p(\mathbf{y}\mid \mathbf{X}) = \mathcal{N}(0, \sigma_a^2 \mathbf{X}\mathbf{X}^T + \sigma_e^2 I_N)
$$

and we see that $\mathbf{X}$, or more precisely the covariance matrix $\mathbf{X}\mathbf{X}^T$ indeed models correlations among the samples in $\mathbf{y}$.

Naturally, fixed and random effects can be combined into a so-called [mixed model](https://en.wikipedia.org/wiki/Mixed_model).

Models of this kind are often used in genetics, where the random effect covariates $\mathbf{X}$ are genetic markers and their covariance matrix $\mathbf{X}\mathbf{X}^T$ expresses the genetic similarity between individuals in the study.

So far we assumed that all random effects had the same variance $\sigma_a^2$. If we assume that each $\alpha_j$ has a different variance $\sigma_j^2$, we would obtain

$$
p(\mathbf{y}\mid \mathbf{X}) = \mathcal{N}(0, \mathbf{K})
$$

where

$$
\mathbf{K} = \sum_{j=1}^p \sigma_j^2 \mathbf{x}_j \mathbf{x}_j^T  + \sigma_e^2 I_N
$$

that is,  $\sigma_j^2$ measures the contribution of covariate $X_j$ to the correlations among samples of $Y$.

To estimate the variance parameters, maximum-likelihood is employed, that is, we find the values of the $\sigma_j^2$ and $\sigma_e^2$ which maximize

$$
\log p(\mathbf{y}) = -\frac12 \log \det (\mathbf{K} ) - \frac12 \mathbf{y}^T \mathbf{K}^{-1} \mathbf{y} - \frac{N}2 \log(2\pi)
$$

This is a difficult, non-convex optimization problem which requires the use of numerical, gradient-based optimizers.

Note that the covariance matrix satisfies

$$
\begin{aligned}
  \mathbb{E}(y_i y_j) = K_{ij}
\end{aligned}
$$
and hence the total variance of $\mathbf{y}$ can be written as

$$
\begin{aligned}
  \mathbb{E}(\mathbf{y}^T\mathbf{y}) = \sum_i \mathbb{E}(y_i^2) = \sum_i K_{ii} = \mathrm{tr}(K)
\end{aligned}
$$

From the definition of $\mathbf{K}$, its trace can be written as

$$
\mathrm{tr}(K) = \sum_{j=1}^p \sigma_j^2 \mathbf{x}_j^T \mathbf{x}_j   + N \sigma_e^2.
$$

Therefore we say that

$$
\frac{\sigma_j^2 \mathbf{x}_j^T \mathbf{x}_j}{\mathrm{tr}(K) }
$$

is the **variance in $\mathbf{y}$ explained by variance component $\mathbf{x}_i$**.

## Kernel-based variance component models

Since the posterior distribution $p(\mathbf{y}\mid \mathbf{X})$ only depends on the matrix $\mathbf{K}$, an immediate generalization is to define variance component models directly in terms of **kernel matrices** instead of starting from a linear, random effect model, namely define

$$
p(\mathbf{y}) = \mathcal{N}(0, \mathbf{K})
$$

For instance, if the elements (samples) $y_i$ of $\mathbf{y}$ are obtained at specific locations (e.g. pixels in an image) $\mathbf{x_i}\in \mathbb{R}^2$ (or $\mathbb{R}^p$ more generally), we can define

$$
K_{ij} = k(\mathbf{x}_i,\mathbf{x}_j)
$$

where the **kernel function** $k$ is a function that measures the similarity or closeness between samples.

Generalizing from before, the kernel matrix  $\mathbf{K}$ can consist of multiple components,

$$
\mathbf{K} = \sum_{j=1}^p \sigma_j^2 \mathbf{K}_j  + \sigma_e^2 I_N
$$

and the **variance explained** by the $j^{\text{th}}$ component is

$$
\frac{\sigma_j^2\mathrm{tr}(\mathbf{K}_j)}{\mathrm{tr}(\mathbf{K})}
$$

## Gaussian processes

In many applications, the "positions" $\mathbf{x}_i$ where samples are obtained are not fixed, but a finite subset of a possibly infinite, continuous range. For instance, when studying a dynamic process, we typically obtain noisy measurements $y_i$ at a finite number of time points $t_i$, and are interested in the entire underlying process $y(t)$ for all times $t$ in some time interval. Similarly, we may have measurements $y_i$ at a finite number of locations $\mathbf{x}_i$, and are interested in the entire function $y(\mathbf{x})$ for all positions $\mathbf{x}$ in some spatial region.

A [Gaussian process](https://en.wikipedia.org/wiki/Gaussian_process) is a [stochastic process](https://en.wikipedia.org/wiki/Stochastic_process), to be understood as a probability distribution over functions $y(\mathbf{x})$ over some continuous domain $D$, such that the set of values of $y(\mathbf{x})$ evaluated at a finite set of points $\mathbf{x}_1,\dots,\mathbf{x}_N\in D$ are jointly normally distributed. A Gaussian process is specified by a kernel function $k(\mathbf{x},\mathbf{x}')$, for $\mathbf{x},\mathbf{x}' \in D$. Writing $\mathbf{y} = (y(\mathbf{x}_1), \dots, y(\mathbf{x}_N))$, the kernel defines the probability distribution

$$
p(\mathbf{y}) = \mathcal{N}(0, \mathbf{K})
$$

where 

$$
K_{ij} = k(\mathbf{x}_i,\mathbf{x}_j)
$$

Hence, for the purposes of parameter estimation, the Gaussian process and variance component viewpoints are identical, which is why the term are often used interchangeably. The main difference lies in the interpretation, where Gaussian process see the data as a finite set of samples from a continuous space. Hence in the Gaussian process viewpoint, we would be particularly interested in **interpolation**, predicting expected values of the function $y(\mathbf{x})$ at locations $\mathbf{x}$ that were not part of the initial (training) data.

Keeping the notation $\mathbf{y}$ for the function values at the training positions $\mathbf{x}_i$, we know that the joint distribution of $\mathbf{y}$ and the value of $y(\mathbf{x})$ at a new position $\mathbf{x}$ is given by

$$
p\left(
\begin{bmatrix}
  \mathbf{y}\\\\
  y(\mathbf{x}) 
\end{bmatrix}
\right) =  \mathcal{N}(0, \mathbf{K}')
$$

where $\mathbf{K}'$ is the $(N+1)\times (N+1)$ matrix

$$
\mathbf{K}' =
\begin{bmatrix}
  \mathbf{K} & \mathbf{k}\\\\
  \mathbf{k}^T & c
\end{bmatrix}
$$

where

$$
\begin{aligned}
  \mathbf{K} &= \left[ k(\mathbf{x}_i,\mathbf{x}_j) \right] \in \mathbb{R}^{N\times N}\\\\
  \mathbf{k} &= k(\mathbf{x}_i,\mathbf{x}) \in \mathbb{R}^{N} \\\\
  c &= k(\mathbf{x},\mathbf{x})\in \mathbb{R}
\end{aligned}
$$

Using [properties of the multivariate normal distribution](https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions), the conditional distribution of the unseen value given the training data is:

$$
\begin{aligned}
  p\left( y(\mathbf{x}) \mid \mathbf{y} \right) &= \mathcal{N}(\mu,\sigma^2)\\\\
  \mu &= \mathbf{k}^T \mathbf{K}^{-1}\mathbf{y}\\\\
  \sigma^2 &= c - \mathbf{k}^T \mathbf{K}^{-1} \mathbf{k}
\end{aligned}
$$

Generalization to interpolation to multiple points is immediate.

## Assignment


{{< alert title="Assignment" >}}


{{< /alert >}}