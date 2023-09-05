# Testing Missing Completely at Random (MCAR) using correlation matrices. 
 ## Abstract
 We study the problem of testing whether the missing values of a potentially high-dimensional dataset are Missing Completely at Random (MCAR). Following a previous work by Berrett and Samworth (2022), we reduce the problem of testing MCAR to the problem of testing the compatibility of the sequence of correlation matrices associated to the missingness patterns, motivated by the fact that this procedure is feasible also in high-dimensional settings. Our first contribution is to define a sensible test statistic for a sequence of correlation matrices, which can be characterised as the optimal value of a Semidefinite Programming (SDP) problem. More precisely, we define our measure of incompatibility for a sequence $`\Sigma_{\mathbb{S}}`$ of correlation matrices, where $`\mathbb{S}`$ is the set of all missingness patterns, as 
```math
     R(\Sigma_{\mathbb{S}})  = \inf \{ \epsilon \in [0,1] : \Sigma_\mathbb{S} \in (1-\epsilon) \tilde{\mathcal{P}}_\mathbb{S}^0 + \epsilon \tilde{\mathcal{P}}_\mathbb{S} \},
```
where $`\tilde{\mathcal{P}}_\mathbb{S}^0`$ is the space of compatible sequences of correlation matrices, and $`\tilde{\mathcal{P}}_\mathbb{S}`$ is the space of sequences of correlation matrices, not necessarily compatible. $`\Sigma_{\mathbb{S}}`$ is compatible if and only if $`R(\Sigma_{\mathbb{S}}) = 0`$, hence, if $`R(\Sigma_{\mathbb{S}}) > 0`$ we can reject the hypothesis of compatibility, and so the hypothesis of MCAR. This results in an oracle test under the extra hypothesis of non-singularity of the sequence of correlation matrices, which relies on a novel concentration inequality for the spectral norm of the difference between the Pearson sample correlation matrix and its population version. More precisely, under the hypothesis of subgaussianity, we prove that 
```math
||\hat{P} - P||_2 \lesssim \sqrt{\frac{d \log d}{n}},
```
with high probability, where $\hat{P}$ is the Pearson sample correlation matrix and $P$ is the associated population version. This leads to a testing procedure whose separation rate is of the order of 
```math
C_{\alpha}\lesssim \frac{1}{c} \max_{S \in \mathbb{S}}\sqrt{\frac{|S|\log(|S||\mathbb{S}|/\alpha)}{n_S}},
```
where $`c > 0`$ is a quantitative measure of non-singularity of $`\Sigma_\mathbb{S}`$,
which we believe to be minimax optimal, at least in a couple of examples. Furthermore, we include a testing procedure that can be applied in practice, which is based on the technique of sample splitting. 

## Scope of this repository
Create an R package to implement this test. 
