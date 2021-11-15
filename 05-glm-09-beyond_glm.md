---
title: "GLM part 9: Beyond GLM"
author: "Nick Green, Kennesaw State University"
output:
  html_document: 
    toc: yes
    toc_float: true
    number_sections: yes
    keep_md: yes
font_size: 16pt
bibliography: stat_refs.bib
csl: ecology.csl
---

# Overview

In other sections we have explored the power and utility of GLMs for biological data. This page describes two methods that can be thought of as extensions of GLM: generalized additive models (GAM) and generalized estimating equations (GEE). GAM will get a more thorough treatment in their own section in the future.

# Generalized additive models (GAM)

**Generalized additive models (GAM)** go one step further than GLMs by relaxing an additional assumption: the assumption of a linear relationship between the response variable on the link scale and the predictors. GAMs fit a curve to data that can vary in complexity. These curves are commonly called **smoothers**. However, the parameters of the curve are not tested individually in the way that the parameters in something like a quadratic model or non-linear model would be. Instead, the curve itself is treated as a regression coefficient. Compare the equation for the linear predictor of a GLM with *k* predictors

$$g\left(E\left(Y\right)\right)=\beta_0+\beta_1X_1+\beta_1X_1+\ldots+\beta_kX_k$$

to the equation of the predictor of a GAM with k predictors:

$$g\left(E\left(Y\right)\right)=\beta_0+f_1\left(X_1\right)+f_2\left(X_2\right)+\ldots+f_k\left(X_k\right)$$

Instead of each predictor being multiplied by a coefficient, a smoothing function of each predictor is estimated. 

When should you use GAM? The most common application is situations where the response is not linearly related to the predictors, or if the response curve takes on a complicated shape. As you might suspect, the smoothers fitted in a GAM can be vulnerable to overfitting just as can be polynomial models or multiple regression models. Researchers usually control for this overfitting of the smoothers by limiting the complexity of the curves; this parameter is usually called the number of “knots”. 

GAMs will get a more thorough treatment later in the context of nonlinear models. Some good references to get you started are @wood2017generalized and @zuur2007analysing. The best-known R package for fitting GAMs is mgcv ([CRAN page](https://cran.r-project.org/web/packages/mgcv/index.html), accessed 2021-09-17).

# Generalized estimating equations (GEE)

All GLMs assume that observations are independent of each other. However, this cannot always be assumed in biological data. Observations may by correlated with each other across time and space, or within treatment or blocking groups. Such violations of the independence assumption can often be handled well by mixed models, which are described in their own pages elsewhere on this site. However, the technique of **generalized estimating equations (GEE)** is another way to analyze data with inherent correlations between observations. Some references to start with if you are interested are @liang1986longitudinal and @hojsgaard2006. Two important packages used with GEE are `gee` [@geepackage] and `geepack` [@hojsgaard2006].

---


[**Go back to main page**](https://greenquanteco.github.io/index.html)

# References

<div id="refs"></div>

# Legal notice

This site is for educational purposes only. This work and its content is released under the [Creative Commons Attribution-ShareAlike 4.0](https://creativecommons.org/licenses/by-sa/4.0/) license. Inclusion of third-party data falls under guidelines of fair use as defined in [section 107 of the US Copyright Act of 1976](https://www.law.cornell.edu/uscode/text/17/107). 
