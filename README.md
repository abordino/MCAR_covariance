# MCAR_covariance: Tests for Missing Completely at Random based on covariance matrices

This repository contains the implementation of the tests presented in [this paper](https://arxiv.org/abs/2401.05256) along with code to reproduce all simulations and experiments. 

## Repository Structure

In the main directory, you will find implementations of the algorithms presented in the paper, along with the code for Little's test as introduced in [this other paper](https://www.tandfonline.com/doi/abs/10.1080/01621459.1988.10478722).

- **`MCAR_test/`**  
  Contains auxiliary functions needed for the implementation of our procedures.

- **`simulCycle/`**  
  Compares our tests with Little's methodologies in the case where the missingness pattern is a d-cycle.
  
- **`simulMissMethods/`**  
  Compares our tests with Little's methodologies in the case of a more general missingness mechanism, which is simulated using the R-package [missMethods](https://cran.r-project.org/web/packages/missMethods/index.html).

