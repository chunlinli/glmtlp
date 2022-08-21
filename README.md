# **glmtlp**: An R Package For Truncated Lasso Penalty

<p align="center">
  <img src="GLMTLP.png" alt="glmtlp" width="500"/>
</p>

Efficient procedures for constrained likelihood estimation and inference with truncated lasso penalty (Shen et al., 2010; Zhang 2010) for linear, generalized linear, and Gaussian graphical models. 

Currently this package is in active development and will be released very soon.
## New features [Unreleased]

This version supports regression from summary statistics and out-of-core model fitting using an ultrahigh-dimensional, multi-gigabyte datasets that cannot be loaded into memory. It's highly optimized and much more memory-efficient as compared to existing penalized regression packages like [**glmnet**](https://github.com/cran/glmnet) and [**ncvreg**](https://github.com/pbreheny/ncvreg/). 

- Add regression with summary data input

- Add inference function

- Add OpenMP support

- Add external memory computation support

- Add sparse coefficient matrix output

- Add implementation of Poisson regression

## Highlights

constrained likelihood approach, inference

an improved algorithm 

any GLM

summary data

big data, memory management

visualization



## References

Li, C., Shen, X., & Pan, W. (2021). Inference for a large directed graphical model with interventions. *arXiv preprint* arXiv:2110.03805. <https://arxiv.org/abs/2110.03805>.

Shen, X., Pan, W., & Zhu, Y. (2012). Likelihood-based selection and sharp parameter estimation. *Journal of the American Statistical Association*, 107(497), 223-232. <https://doi.org/10.1080/01621459.2011.645783>.

Shen, X., Pan, W., Zhu, Y., & Zhou, H. (2013). On constrained and regularized high-dimensional regression. *Annals of the Institute of Statistical Mathematics*, 65(5), 807-832. <https://doi.org/10.1007/s10463-012-0396-3>.

Tibshirani, R., Bien, J., Friedman, J., Hastie, T., Simon, N., Taylor, J., & Tibshirani, R. J. (2012). Strong rules for discarding predictors in lasso‐type problems. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 74(2), 245-266. <https://doi.org/10.1111/j.1467-9868.2011.01004.x>.

Yang, Y. & Zou, H. A coordinate majorization descent algorithm for l1 penalized learning. *Journal of Statistical Computation and Simulation* 84.1 (2014): 84-95. <https://doi.org/10.1080/00949655.2012.695374>.

Zhu, Y., Shen, X., & Pan, W. (2020). On high-dimensional constrained maximum likelihood inference. *Journal of the American Statistical Association*, 115(529), 217-230. <https://doi.org/10.1080/01621459.2018.1540986>.

Zhu, Y. (2017). An augmented ADMM algorithm with application to the generalized lasso problem. *Journal of Computational and Graphical Statistics*, 26(1), 195-204. <https://doi.org/10.1080/10618600.2015.1114491>.

Part of the code is adapted from [**glmnet**](https://github.com/cran/glmnet), [**ncvreg**](https://github.com/pbreheny/ncvreg/), and [**biglasso**](https://github.com/YaohuiZeng/biglasso).

**Warm thanks to the authors of above open-sourced softwares.**