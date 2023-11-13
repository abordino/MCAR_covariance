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


1. **Find** $`\hat{Q}_\mathbb{S}`$ from the dual decomposition $`\hat{\Sigma}_\mathbb{S} = (1-R(\hat{\Sigma}_\mathbb{S}))\hat{Q}_\mathbb{S} + R(\hat{\Sigma}_\mathbb{S})\hat{\Sigma}'_\mathbb{S}`$. At the same time, **define** $`R^{(0)} = R(\hat{\Sigma}_\mathbb{S})`$.
2. Rotate the original data $`X_\mathbb{S}`$, i.e. for all $`S \in \mathbb{S}`$, for all $`i \in [n_S]`$ **do** $`\tilde{X}_{S,i} = \hat{Q}_S^{1/2}\hat{\Sigma}_S^{-1/2}X_{S,i}`$.
3. For all $`b \in [B]`$
4.  For all $`S \in \mathbb{S}`$, **bootstrap** from $`\tilde{X}_S`$. 
5.  $`\quad`$ **Compute** $`\Sigma_\mathbb{S}^{(b)} = (\Sigma_S^{(b)})_{S \in \mathbb{S}}`$, where $`\Sigma_S^{(b)} = Corr(\tilde{X}_S^{(b)})`$.
6.  $`\quad`$ **Compute** $`R^{(b)} = R(\Sigma_\mathbb{S}^{(b)})`$.
7. **Reject** $`H_0`$ if and only if
```math
    1 + \sum_{i=1}^B 1 \{R^{(b)} \geq R^{(0)}\} \geq \alpha(1+B).
 ```
