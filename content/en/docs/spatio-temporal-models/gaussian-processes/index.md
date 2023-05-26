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
p(\mathbf{y}\mid \mathbf{X}) = \mathcal{N}(0, \mathbf{K} + \sigma_e^2 I_N)
$$

where

$$
\mathbf{K} = \sum_{j=1}^p \sigma_j^2 \mathbf{x}_j \mathbf{x}_j^T
$$

that is,  $\sigma_j^2$ measures the contribution of covariate $X_j$ to the correlations among samples of $Y$.

To estimate the variance parameters, maximum-likelihood is employed, that is, we find the values of the $\sigma_j^2$ and $\sigma_e^2$ which maximize

$$
\log p(\mathbf{y}\mid \mathbf{X}) = -\frac12 \log \det \bigl(\mathbf{K} + \sigma_e^2 I_N\bigr) - \frac12 \mathbf{y}^T \bigl(\mathbf{K} + \sigma_e^2 I_N\bigr)^{-1} \mathbf{y} - \frac{N}2 \log(2\pi)
$$

This is a difficult, non-convex optimization problem which requires the use of numerical, gradient-based optimizers.

## Kernel-based variance component models



## Gaussian processes

## Assignment


{{< alert title="Assignment" >}}


{{< /alert >}}