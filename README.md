# Testing Missing Completely at Random (MCAR) using covariance matrices. 
 ## Abstract
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
Develop an R package to implement a bootstrap version of this test. We will use as test statistic $`T = R(\Sigma_\mathbb{S}) + V(\sigma_{\mathbb{S}}^2) + M(\mu_{\mathbb{S}})`$, where 
```math
 M(\mu_{\mathbb{S}}) = \max_{j \in [d]} \max_{S_1, S_2 \in \mathbb{S}_j} |\mu_{S_1,j} - \mu_{S_2,j}|,
```
 and observe that $`M(\mu_{\mathbb{S}}) = 0`$ if and only if $`\mu_{\mathbb{S}}`$ is consistent.


1. Given data $`X_\mathbb{S}`$, compute $`\hat{\sigma}^2_\mathbb{S} = VarX_\mathbb{S}`$, and rescale it such that $`av(\hat{\sigma}^2_{\mathbb{S}_j}) = 1`$ for all $`j \in [d]`$; i.e. replace $`X_j`$ with $`X_j/\sqrt{av(\hat{\sigma}^2_{\mathbb{S}_j})}`$ for all $`j \in [d]`$.
2. Compute $`\hat{\mu}_\mathbb{S} = \mathbb{E}X_\mathbb{S}`$, $`\hat{M}_\mathbb{S} = CovX_\mathbb{S} = diag( \sigma_{\mathbb{S}}^2)^{1/2} \cdot \Sigma_{\mathbb{S}} \cdot diag(\sigma_{\mathbb{S}}^2)^{1/2}`$.
3. Compute $`T^{(0)} = R(\hat{\Sigma}_\mathbb{S}) + V(\hat{\sigma}_{\mathbb{S}}^2) + M(\hat{\mu}_{\mathbb{S}})`$, and compute at the same time the dual decomposition $`\hat{\Sigma}_\mathbb{S} = (1-R(\hat{\Sigma}_\mathbb{S}))\hat{Q}_\mathbb{S} + R(\hat{\Sigma}_\mathbb{S})\hat{\Sigma}'_\mathbb{S}`$.
4. Rotate the original data $`X_\mathbb{S}`$, i.e. for all $`S \in \mathbb{S}`$, for all $`i \in [n_S]`$ do $`\tilde{X}_{S,i} = \hat{Q}_S^{1/2}\hat{M}_S^{-1/2}(X_{S,i}-\hat{\mu}_S + \hat{\mu}_{|S})`$, where $`\hat{\mu}_j = |\mathbb{S}_j|^{-1}\sum_{S\in \mathbb{S}_j} \mu_{S,j}`$.
5. for $`b \in [B]`$:
    5.1. For all $`S \in \mathbb{S}`$, let $`\tilde{X}_{S,i}^{(b)}`$ be a nonparametric bootstrap sample from $`\tilde{X}_{S,i}`$, for $`i \in [n_S]`$. 
    5.2. Compute $`\hat{\mu}_{\mathbb{S},b} = \mathbb{E}X_\mathbb{S}^{(b)}`$, $`\hat{M}_{\mathbb{S},b} = CovX_\mathbb{S}^{(b)} = diag( \hat{\sigma}_{\mathbb{S},b}^2)^{1/2} \cdot \hat{\Sigma}_{\mathbb{S},b} \cdot diag( \hat{\sigma}_{\mathbb{S},b}^2)^{1/2}`$.
    5.3. Compute $`T^{(b)} = R(\hat{\Sigma}_{\mathbb{S},b}) + V(\hat{\sigma}_{\mathbb{S},b}^{2}) + M(\hat{\mu}_{\mathbb{S},b})`$.
6. Reject $`H_0`$ if and only if
 ```math
1 + \sum_{i=1}^B 1\{T^{(b)} \geq T^{(0)}\} \geq \alpha(1+B).
```
