Malnutritional analysis
================

This model uses code developed by Nicky Welton to estimate the
assocaition between malnutrition and pneumonia
mortality.

``` r
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE)
knitr::opts_knit$set(root.dir = here::here())
```

Packages

# Run model with example data

First check the model runs succesfully in Jags using example data.

## Example data

Need to convert the example data into matrices and vectors for JAGS.

| r1 |   n1 |  r2 |   n2 | r3 |  n3 | g1 | g2 | g3 | na | pi2 |
| -: | ---: | --: | ---: | -: | --: | -: | -: | -: | -: | --: |
|  9 |  140 |  23 |  140 | 10 | 138 |  1 |  2 |  3 |  3 | 1.0 |
| 11 |   78 |  12 |   85 | 29 | 170 |  1 |  2 |  3 |  3 | 1.0 |
| 75 |  731 | 363 |  714 | NA |   1 |  1 |  2 | NA |  2 | 1.0 |
|  2 |  106 |   9 |  205 | NA |   1 |  1 |  2 | NA |  2 | 1.0 |
| 58 |  549 | 237 | 1561 | NA |   1 |  1 |  2 | NA |  2 | 1.0 |
|  0 |   33 |   9 |   48 | NA |   1 |  1 |  2 | NA |  2 | 1.0 |
|  3 |  100 |  31 |   98 | NA |   1 |  1 |  2 | NA |  2 | 1.0 |
|  1 |   31 |  26 |   95 | NA |   1 |  1 |  2 | NA |  2 | 1.0 |
|  6 |   39 |  17 |   77 | NA |   1 |  1 |  2 | NA |  2 | 1.0 |
| 79 |  702 |  77 |  694 | NA |   1 |  2 |  3 | NA |  2 | 1.0 |
| 18 |  671 |  21 |  535 | NA |   1 |  2 |  3 | NA |  2 | 1.0 |
| 64 |  642 | 107 |  761 | NA |   1 | 12 |  3 | NA |  2 | 0.6 |
|  5 |   62 |   8 |   90 | NA |   1 | 12 |  3 | NA |  2 | 0.5 |
| 20 |  234 |  34 |  237 | NA |   1 | 12 |  3 | NA |  2 | 0.4 |
|  0 |   20 |   9 |   20 | NA |   1 |  1 | 23 | NA |  2 | 0.7 |
|  8 |  116 |  19 |  149 | NA |   1 |  1 | 23 | NA |  2 | 0.5 |
| 95 | 1107 | 143 | 1031 | NA |   1 |  1 | 23 | NA |  2 | 0.6 |
| 15 |  187 |  36 |  504 | NA |   1 |  1 | 23 | NA |  2 | 0.5 |
| 78 |  584 |  73 |  675 | NA |   1 |  1 | 23 | NA |  2 | 0.6 |
| 69 | 1177 |  54 |  888 | NA |   1 |  1 | 23 | NA |  2 | 0.3 |
| 20 |   49 |  16 |   43 | NA |   1 |  1 | 23 | NA |  2 | 0.8 |
|  7 |   66 |  32 |  127 | NA |   1 | 12 |  3 | NA |  2 | 0.7 |
| 12 |   76 |  20 |   74 | NA |   1 |  1 |  3 | NA |  2 | 1.0 |
|  9 |   55 |   3 |   26 | NA |   1 | 12 |  3 | NA |  2 | 0.7 |

## Run model

The following is the BUGS/JAGS code for the fixed effects
    model.

    ## model{                                                                        for(i in 1:ns){                                                               # LOOP THROUGH STUDIES     delta[i,1]<-0   mu[i] ~ dnorm(0,.0001)                                              # vague priors for all trial baselines   for (k in 1:na[i]) {                                                       # LOOP THROUGH GROUPS     r[i,k] ~ dbin(p[i,k],n[i,k])                                           # binomial likelihood     logit(p[i,k]) <- mu[i] + delta[i,k]                                   # model for linear predictor      rhat[i,k] <- p[i,k] * n[i,k]                                            # expected value of the numerators      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k]))             #Deviance contribution          + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))   }   resdev[i] <- sum(dev[i,1:na[i]])                        # summed residual deviance contribution for this trial   for (k in 2:na[i]) {                                           # LOOP THROUGH ARMS      delta[i,k] <-  di[i,g[i,k]] - di[i,g[i,1]]             # NMA model   }   for (k in 1:ng){     di[i,k]<-d[k]     }   di[i,12]<-pi2[i]*d[2]   di[i,23]<-pi2[i]*d[2]+(1-pi2[i])*d[3] }  totresdev <- sum(resdev[])                                           #Total Residual Deviance d[1]<- 0                                                                      # group effect is zero for reference group for (k in 2:ng)  { d[k] ~ dnorm(0,.0001)}                           # vague priors for group effects  # pairwise ORs and LORs for all possible pair-wise comparisons for (c in 1:(ng-1)) {  for (k in (c+1):ng) {        or[c,k] <- exp(d[k] - d[c])        lor[c,k] <- (d[k]-d[c])       } }  }                                                                                 # *** PROGRAM ENDS

## Initialise model and return results

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 50
    ##    Unobserved stochastic nodes: 26
    ##    Total graph size: 1371
    ## 
    ## Initializing model

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 5
    ##    Unobserved stochastic nodes: 3
    ##    Total graph size: 40
    ## 
    ## Initializing model

    ## 
    ## Iterations = 2001:4000
    ## Thinning interval = 1 
    ## Number of chains = 1 
    ## Sample size per chain = 2000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##        Mean      SD Naive SE Time-series SE
    ## d[1] 0.0000 0.00000 0.000000       0.000000
    ## d[2] 0.9260 0.06882 0.001539       0.004653
    ## d[3] 0.6201 0.08599 0.001923       0.005042
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##        2.5%    25%    50%    75%  97.5%
    ## d[1] 0.0000 0.0000 0.0000 0.0000 0.0000
    ## d[2] 0.7910 0.8806 0.9271 0.9712 1.0626
    ## d[3] 0.4442 0.5657 0.6208 0.6811 0.7794
