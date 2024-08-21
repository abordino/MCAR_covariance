# Testing Missing Completely at Random (MCAR) using covariance matrices. 
 ## Abstract, based on https://arxiv.org/abs/2401.05256
 We study the problem of testing whether the missing values of a potentially high-dimensional dataset are Missing Completely at Random (MCAR). We reduce the problem of testing MCAR to the problem of testing the compatibility of the sequence of covariance matrices associated to the missingness patterns, motivated by the fact that this procedure is feasible also in high-dimensional settings. We check the compatibility of covariance matrices by checking the consistency of variances and compatibility of correlation matrices. Our first contribution is to define a sensible measure of incompatibility for a sequence of correlation matrices, which can be characterised as the optimal value of a Semi-definite Programming (SDP) problem. More precisely, we define our measure of incompatibility for a sequence $\Sigma_{\mathbb{S}}$ of correlation matrices, where $\mathbb{S}$ is the set of all missingness patterns, as 
```math
     R(\Sigma_{\mathbb{S}})  = \inf \{ \epsilon \in [0,1] : \Sigma_\mathbb{S} \in (1-\epsilon) \mathcal{P}_\mathbb{S}^0 + \epsilon \mathcal{P}_\mathbb{S} \},
```
where $`\mathcal{P}_\mathbb{S}^0`$ is the space of compatible sequences of correlation matrices, and $`\mathcal{P}_\mathbb{S}`$ is the space of all sequences of correlation matrices. We estimate this measure by the natural plug-in estimator and study its concentration properties to introduce suitable critical values. This results in an oracle test under the extra hypothesis of non-singularity of the sequence of correlation matrices, which relies on a novel concentration inequality for the spectral norm of the difference between the Pearson sample correlation matrix and its population version. This leads to a testing procedure whose separation rate is of the order of 
```math
\max_{S \in \mathbb{S}}\sqrt{\frac{|S|\log(|S\|\mathbb{S}|/\alpha)}{ n_S}}.
```
 We prove that this separation is minimax optimal in some examples of interest, and further validate our methodology with some numerical simulations.

## Scope of this repository
This repository contains the code for all the simulations in the paper.
