
<!-- README.md is generated from README.Rmd. Please edit that file -->

# halfsibdesign

<!-- badges: start -->
<!-- badges: end -->

A package containing several algorithms for finding REML estimates for
covariances in the balanced *q*-dimensional half-sib design
*y<sub>ijk</sub> = *μ* + *α*<sub>i</sub> + *β*<sub>ij</sub> +
*ϵ*<sub>ijk</sub>*, *1 ≤ *i* ≤ *I*, 1 ≤ *j* ≤ *J*, 1 ≤ *k* ≤ *K**, where
**α*<sub>i</sub>*  ∼ 𝒩(0,*A*), **β*<sub>ij</sub>*  ∼ 𝒩(0,*B*),
**ϵ*<sub>ijk</sub>*  ∼ 𝒩(0,*E*).

## Installation

You can install the development version of halfsibdesign from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("damian-t-p/halfsibdesign")
```

## Example

Simulate a half-sib experiment with specified parameters:

``` r
set.seed(1)

library(halfsibdesign)

q <- 4 # number of traits

I <- 100 # number of sires
J <- 3 # number of dams
K <- 5 # number of individuals per line

mu <- 1:q

sigma_a <- 5
sigma_b <- 3
sigma_e <- 1

A <- sigma_a^2 * diag(c(0, 0, 1, 1))
B <- sigma_b^2 * diag(q)
E <- sigma_e^2 * diag(q)

df <- rhalfsib(mu, A, I, B, J, E, K)
```

First, perform a MANOVA fit

``` r
manova_2way(df)
#> $S1
#>              [,1]          [,2]          [,3]        [,4]
#> [1,]  1.054245901 -0.0516353150 -0.0042920491 -0.02807315
#> [2,] -0.051635315  0.9371151905  0.0008297587 -0.01024345
#> [3,] -0.004292049  0.0008297587  1.0084415477  0.02166148
#> [4,] -0.028073153 -0.0102434479  0.0216614828  1.08337053
#> 
#> $S2
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  8.2208523  0.2058116  0.9053192 -0.9549915
#> [2,]  0.2058116  9.9534095 -0.2880542  0.6135789
#> [3,]  0.9053192 -0.2880542 10.8004023  0.9599887
#> [4,] -0.9549915  0.6135789  0.9599887 11.3201377
#> 
#> $S3
#>              [,1]       [,2]         [,3]       [,4]
#> [1,]  0.544404629 -0.1397165  0.006040896  0.9287818
#> [2,] -0.139716531 -0.2036637  1.401206022 -0.1107346
#> [3,]  0.006040896  1.4012060 22.991012239 -1.0674308
#> [4,]  0.928781804 -0.1107346 -1.067430766 18.6056761
#> 
#> $reml_crit
#> [1] -8112.177
```

Notice that `S3`, the within-sires estimate is not non-negative
definite, se we compute the correponding REML estimate.

``` r
stepreml_2way(df)
#> $S1
#>              [,1]          [,2]          [,3]        [,4]
#> [1,]  1.054245901 -0.0516353150 -0.0042920491 -0.02807315
#> [2,] -0.051635315  0.9371151905  0.0008297587 -0.01024345
#> [3,] -0.004292049  0.0008297587  1.0084415477  0.02166148
#> [4,] -0.028073153 -0.0102434479  0.0216614828  1.08337053
#> 
#> $S2
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  8.2130940  0.1575627  0.9090954 -0.9564348
#> [2,]  0.1575627  9.6533524 -0.2645701  0.6046033
#> [3,]  0.9090954 -0.2645701 10.7985643  0.9606912
#> [4,] -0.9564348  0.6046033  0.9606912 11.3198692
#> 
#> $S3
#>              [,1]        [,2]         [,3]       [,4]
#> [1,]  0.552215223 -0.09114281  0.002239256  0.9302348
#> [2,] -0.091142809  0.09841398  1.377563800 -0.1016985
#> [3,]  0.002239256  1.37756380 22.992862606 -1.0681380
#> [4,]  0.930234791 -0.10169853 -1.068137977 18.6059464
#> 
#> $reml_crit
#> [1] -8112.313
```
