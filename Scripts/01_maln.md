Malnutrition analysis
================

This model uses code developed by Nicky Welton to estimate the
association between malnutrition and pneumonia mortality.

``` r
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, cache = TRUE)
knitr::opts_knit$set(root.dir = here::here())
```

Packages

``` r
library(tidyverse)
library(rjags)
library(mcmcplots)
library(ggplot2)
library(R2OpenBUGS)
```

# Nomenclature

For all analyses the malnutrition categories are interpreted as folows:-

| coding | label              |
|--------|--------------------|
| 1      | None               |
| 2      | Moderate           |
| 3      | Severe             |
| 12     | None or moderate   |
| 23     | Moderate or severe |

# Fixed effects model

## Run model with example data

First check the model runs successfully in Jags using example data.

Need to convert the example data into matrices and vectors for JAGS.

``` r
metadata <- list(ng=3,ns=24)
ex_df <- read_delim("Supporting/example_data.txt", delim = "\t")
matrices <- list(r = NA, n = NA, g = NA)
matrices[] <- map(names(matrices), ~ ex_df %>% 
                  select(starts_with(.x), -na) %>% 
                  as.matrix())
na <- ex_df %>% pull(na)
pi2 <- ex_df %>% pull(pi2)


list_data <- c(matrices, 
               list(na = na, pi2 = pi2),
               metadata)
knitr::kable(ex_df)
```

|  r1 |   n1 |  r2 |   n2 |  r3 |  n3 |  g1 |  g2 |  g3 |  na | pi2 |
|----:|-----:|----:|-----:|----:|----:|----:|----:|----:|----:|----:|
|   9 |  140 |  23 |  140 |  10 | 138 |   1 |   2 |   3 |   3 | 1.0 |
|  11 |   78 |  12 |   85 |  29 | 170 |   1 |   2 |   3 |   3 | 1.0 |
|  75 |  731 | 363 |  714 |  NA |   1 |   1 |   2 |  NA |   2 | 1.0 |
|   2 |  106 |   9 |  205 |  NA |   1 |   1 |   2 |  NA |   2 | 1.0 |
|  58 |  549 | 237 | 1561 |  NA |   1 |   1 |   2 |  NA |   2 | 1.0 |
|   0 |   33 |   9 |   48 |  NA |   1 |   1 |   2 |  NA |   2 | 1.0 |
|   3 |  100 |  31 |   98 |  NA |   1 |   1 |   2 |  NA |   2 | 1.0 |
|   1 |   31 |  26 |   95 |  NA |   1 |   1 |   2 |  NA |   2 | 1.0 |
|   6 |   39 |  17 |   77 |  NA |   1 |   1 |   2 |  NA |   2 | 1.0 |
|  79 |  702 |  77 |  694 |  NA |   1 |   2 |   3 |  NA |   2 | 1.0 |
|  18 |  671 |  21 |  535 |  NA |   1 |   2 |   3 |  NA |   2 | 1.0 |
|  64 |  642 | 107 |  761 |  NA |   1 |  12 |   3 |  NA |   2 | 0.6 |
|   5 |   62 |   8 |   90 |  NA |   1 |  12 |   3 |  NA |   2 | 0.5 |
|  20 |  234 |  34 |  237 |  NA |   1 |  12 |   3 |  NA |   2 | 0.4 |
|   0 |   20 |   9 |   20 |  NA |   1 |   1 |  23 |  NA |   2 | 0.7 |
|   8 |  116 |  19 |  149 |  NA |   1 |   1 |  23 |  NA |   2 | 0.5 |
|  95 | 1107 | 143 | 1031 |  NA |   1 |   1 |  23 |  NA |   2 | 0.6 |
|  15 |  187 |  36 |  504 |  NA |   1 |   1 |  23 |  NA |   2 | 0.5 |
|  78 |  584 |  73 |  675 |  NA |   1 |   1 |  23 |  NA |   2 | 0.6 |
|  69 | 1177 |  54 |  888 |  NA |   1 |   1 |  23 |  NA |   2 | 0.3 |
|  20 |   49 |  16 |   43 |  NA |   1 |   1 |  23 |  NA |   2 | 0.8 |
|   7 |   66 |  32 |  127 |  NA |   1 |  12 |   3 |  NA |   2 | 0.7 |
|  12 |   76 |  20 |   74 |  NA |   1 |   1 |   3 |  NA |   2 | 1.0 |
|   9 |   55 |   3 |   26 |  NA |   1 |  12 |   3 |  NA |   2 | 0.7 |

## Run model

The following is the BUGS/JAGS code for the fixed effects model.

``` r
a <- read_lines("Supporting/FE_model.txt")
print(a)
```

    ##  [1] "model{                                                                       "                                           
    ##  [2] "for(i in 1:ns){                                                               # LOOP THROUGH STUDIES"                    
    ##  [3] "    delta[i,1]<-0"                                                                                                       
    ##  [4] "  mu[i] ~ dnorm(0,.0001)                                              # vague priors for all trial baselines"            
    ##  [5] "  for (k in 1:na[i]) {                                                       # LOOP THROUGH GROUPS"                      
    ##  [6] "    r[i,k] ~ dbin(p[i,k],n[i,k])                                           # binomial likelihood"                        
    ##  [7] "    logit(p[i,k]) <- mu[i] + delta[i,k]                                   # model for linear predictor"                  
    ##  [8] "     rhat[i,k] <- p[i,k] * n[i,k]                                            # expected value of the numerators"         
    ##  [9] "     dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k]))             #Deviance contribution"                          
    ## [10] "         + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))"                                              
    ## [11] "  }"                                                                                                                     
    ## [12] "  resdev[i] <- sum(dev[i,1:na[i]])                        # summed residual deviance contribution for this trial"        
    ## [13] "  for (k in 2:na[i]) {                                           # LOOP THROUGH ARMS"                                    
    ## [14] "     delta[i,k] <-  di[i,g[i,k]] - di[i,g[i,1]]             # NMA model"                                                 
    ## [15] "  }"                                                                                                                     
    ## [16] "  for (k in 1:ng){"                                                                                                      
    ## [17] "    di[i,k]<-d[k]"                                                                                                       
    ## [18] "    }"                                                                                                                   
    ## [19] "  di[i,12]<-pi2[i]*d[2]"                                                                                                 
    ## [20] "  di[i,23]<-pi2[i]*d[2]+(1-pi2[i])*d[3]"                                                                                 
    ## [21] "}"                                                                                                                       
    ## [22] ""                                                                                                                        
    ## [23] "totresdev <- sum(resdev[])                                           #Total Residual Deviance"                           
    ## [24] "d[1]<- 0                                                                      # group effect is zero for reference group"
    ## [25] "for (k in 2:ng)  { d[k] ~ dnorm(0,.0001)}                           # vague priors for group effects"                    
    ## [26] ""                                                                                                                        
    ## [27] "# pairwise ORs and LORs for all possible pair-wise comparisons"                                                          
    ## [28] "for (c in 1:(ng-1)) {  for (k in (c+1):ng) {"                                                                            
    ## [29] "       or[c,k] <- exp(d[k] - d[c])"                                                                                      
    ## [30] "       lor[c,k] <- (d[k]-d[c])"                                                                                          
    ## [31] "      }"                                                                                                                 
    ## [32] "}"                                                                                                                       
    ## [33] ""                                                                                                                        
    ## [34] "}                                                                                 # *** PROGRAM ENDS"

## Check model runs with example data

``` r
mod_example <- jags.model(file = "Supporting/FE_model.txt",
                            data = list_data, n.chains = 1, n.adapt = 1000)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 50
    ##    Unobserved stochastic nodes: 26
    ##    Total graph size: 1074
    ## 
    ## Initializing model

``` r
update(mod_example, 100)

data(LINE)
LINE$recompile()
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 5
    ##    Unobserved stochastic nodes: 3
    ##    Total graph size: 36
    ## 
    ## Initializing model

``` r
LINE.out <- coda.samples(mod_example, c("d"),n.iter = 100)
summary(LINE.out)
```

    ## 
    ## Iterations = 1101:1200
    ## Thinning interval = 1 
    ## Number of chains = 1 
    ## Sample size per chain = 100 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##        Mean      SD Naive SE Time-series SE
    ## d[1] 0.0000 0.00000 0.000000        0.00000
    ## d[2] 0.9420 0.06647 0.006647        0.01660
    ## d[3] 0.5711 0.08785 0.008785        0.01902
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##        2.5%    25%    50%    75%  97.5%
    ## d[1] 0.0000 0.0000 0.0000 0.0000 0.0000
    ## d[2] 0.8414 0.8892 0.9344 0.9906 1.0770
    ## d[3] 0.3986 0.5215 0.5762 0.6209 0.7481

## Process real data into format for analysis

Process real data, in the first instance assume that the proportion in
each category is known for all data.

First rename variables and recode. There should be a maximum of three
groups for each study/malnutrition category type combination.

## Categorise pattern of collapsing for each study/manutrition category combination

Calculate which studies have missing event data, and which have missing
“n” data. In the present analysis, any with complete n (ie count) data
have complete events data.

``` r
mort <- read_csv("Data/Mortality_Numbers.csv")
names(mort) <- str_to_lower(names(mort)) 
names(mort) <- str_replace_all(names(mort), "\\s|\\-", "_")

if(params$restrict == "pre"){
  mort <- mort %>% 
    filter(pre_or_post_2000 == "pre")
}

if(params$restrict == "post"){
  mort <- mort %>% 
    filter(pre_or_post_2000 == "post")
}

mort <- mort %>% 
  select(-pre_or_post_2000)

mort <- mort %>% 
  mutate_at(vars(survived, died, total), as.integer)

mort <- mort %>% 
  group_by(study, measure) %>% 
  mutate(all_missing = all(is.na(total))) %>% 
  ungroup() %>% 
  filter(!all_missing) %>% 
  select(-all_missing)

maln_cat3 <- c("None", "Moderate", "Severe")

mort <- mort %>% 
  group_by(study, malnutrition_category) %>% 
  mutate(all_n =  all(maln_cat3 %in% malnutrition_severity) & !any(is.na(total)),
            all_r =  all(maln_cat3 %in% malnutrition_severity) & !any(is.na(died)),
            col12 = !all(maln_cat3 %in% malnutrition_severity) & any(malnutrition_severity == "All"),
            col23 = !all(maln_cat3 %in% malnutrition_severity) & any(malnutrition_severity == "Non-Severe")) %>% 
  ungroup() 

all_present_smry <- mort %>% 
  distinct(study, malnutrition_category, .keep_all = TRUE) %>% 
  group_by(all_n, all_r, col12, col23) %>% 
  summarise(studies = sum(!duplicated(study)),
            studies_times_categories = n())
col12 <- all_present_smry$studies[all_present_smry$col12]
col23 <- all_present_smry$studies[all_present_smry$col23]
knitr::kable(all_present_smry)
```

| all\_n | all\_r | col12 | col23 | studies | studies\_times\_categories |
|:-------|:-------|:------|:------|--------:|---------------------------:|
| FALSE  | FALSE  | FALSE | TRUE  |       7 |                          7 |
| FALSE  | FALSE  | TRUE  | FALSE |       3 |                          4 |
| TRUE   | TRUE   | FALSE | FALSE |      16 |                         24 |

Of those studies with collapsed data, 3 collapse into none/moderate and
7 collapse into moderate severe.

## Calculate the proportion in group 2 for the collapsed studies

Since all of those with missing event data have missing totals for each
category, we can only examine the proportion in group 2 where there is
complete data. Calculate the proportion in the second category (moderate
malnutrition) within each collapsed category. Having done so collapse
the Ns. Where the proportion is unknown, assume it is the same as the
mean proportion.

First need to classify which are collapsed.

``` r
pi2 <- mort %>% 
  filter(all_n) %>% 
  select(study, malnutrition_category, malnutrition_severity, total) %>% 
  mutate(total = as.integer(total)) %>% 
  spread(malnutrition_severity, total) %>% 
  mutate(pi2_coll12 = Moderate/(Moderate + None),
         pi2_coll23 = Moderate/(Moderate + Severe)) %>% 
  select(study, malnutrition_category, pi2_coll12, pi2_coll23)


pi2_for_mdl <- mort %>% 
  filter(all_n) %>% 
  select(study, malnutrition_category, malnutrition_severity, total) %>% 
  mutate(total = as.integer(total)) %>% 
  spread(malnutrition_severity, total) %>% 
  group_by(study, malnutrition_category) %>% 
  summarise_all(sum) %>% 
  ungroup()

cat_names <- unique(pi2_for_mdl$malnutrition_category) %>%  sort()
pi2_for_mdl <- map(cat_names, ~ pi2_for_mdl %>% 
                     filter(malnutrition_category == .x) %>% 
                     select(None, Moderate, Severe) %>% 
                     as.matrix())
names(pi2_for_mdl) <- cat_names

pi2_smry <- pi2 %>% 
  group_by(malnutrition_category) %>% 
  summarise_at(vars(pi2_coll12, pi2_coll23), mean) %>% 
  ungroup()

pi2_lng <- pi2 %>% 
  gather("collapsed_categories", "value", pi2_coll12, pi2_coll23) 
plot_dist <- ggplot(pi2_lng, aes(x = value, fill = collapsed_categories)) + geom_histogram() +
  facet_wrap(collapsed_categories ~ malnutrition_category)
plot_dist
```

![](01_maln_files/figure-gfm/calculatenproprs-1.png)<!-- -->

The plot above shows the proportion of participants category 2 of 1 + 2
and in category 2 of 2+3 for all three types of measure. The proportion
in 2 (moderate) is sometimes surprisingly low or high.

``` r
mort_slct <- mort %>% 
  select(study, malnutrition_category, n = total, r = died, maln = malnutrition_severity,
         col12, col23) %>% 
  inner_join(pi2_smry) %>%
  mutate(pi2 = case_when(
    col12 ~ pi2_coll12,
    col23 ~ pi2_coll23,
    TRUE ~ 1
  ))
  

grp_lbls <-  mort_slct %>% 
  group_by(study, malnutrition_category) %>% 
  mutate(g = seq_along(study)) %>% 
  ungroup() %>% 
  mutate(g_lbl = case_when(
    maln == "None" ~ 1L,
    maln == "Moderate" ~ 2L,
    maln == "Severe" ~ 3L,
    maln == "Non-Severe" ~ 12L,
    maln == "All" ~ 23L)
  )
```

First run for the “w/a” definition of malnutrition.

This code rearranges the dataframe so that we have a matrix of N’s,
events and group labels as per the earlier structure.

``` r
grp_lbls2 <- grp_lbls %>% 
  select(study, malnutrition_category, r, n, g, g_lbl, pi2) %>% 
  filter(malnutrition_category == "w/a") %>% 
  distinct(study, g, .keep_all = TRUE)

grp_lbls2n <- grp_lbls2 %>% 
  select(-r, -g_lbl) %>% 
  spread(g, n)

grp_lbls2n <- grp_lbls2 %>% 
  select(-r, -g_lbl) %>% 
  spread(g, n)

grp_lbls2r <- grp_lbls2 %>% 
  select(-n, -g_lbl) %>% 
  spread(g, r)

grp_lbls2g <- grp_lbls2 %>% 
  select(-n, -r) %>% 
  spread(g, g_lbl)

grp_lbls2_col12 <- grp_lbls2 %>% 
  group_by(study) %>% 
  summarise(res = any(12 %in% g_lbl)) %>% 
  pull(res)
```

Check that the restructuring has kept the order of the studies. The
following tests should be TRUE.

``` r
identical(grp_lbls2n %>% select(1:2),
          grp_lbls2g %>%  select(1:2))
```

    ## [1] TRUE

``` r
identical(grp_lbls2n %>% select(1:2),
          grp_lbls2r %>%  select(1:2))
```

    ## [1] TRUE

## Try running models on real data

Simplest model, assume that the proportion with moderate malnutrition of
those with none-moderate, and the proportion with moderate of those with
moderate-severe in those studies with collapsed data is the same as the
average for those studies where there data are not collapsed.

``` r
na <- 3 - apply(grp_lbls2n %>% select(`1`, `2`, `3`) %>% as.matrix(), 1, function(x) x %>% 
        as.integer() %>% 
        is.na() %>% 
        sum())
list_data2 <- list(r = grp_lbls2r %>% select(`1`, `2`, `3`) %>% as.matrix(),
                   n = grp_lbls2n %>% select(`1`, `2`, `3`) %>% as.matrix(),
                   g = grp_lbls2g %>% select(`1`, `2`, `3`) %>% as.matrix(),
                   na = na,
                   pi2 = grp_lbls2$pi2,
                   ng = 3,
                   ns = nrow(grp_lbls2n))
```

``` r
mod1 <- jags.model(file = "Supporting/FE_model.txt",
                            data = list_data2, n.chains = 2, n.adapt = 1000)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 61
    ##    Unobserved stochastic nodes: 25
    ##    Total graph size: 1256
    ## 
    ## Initializing model

``` r
update(mod1, 1000)
data(LINE)
LINE$recompile()
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 5
    ##    Unobserved stochastic nodes: 3
    ##    Total graph size: 36
    ## 
    ## Initializing model

``` r
LINE.out <- coda.samples(mod1, c("d"),n.iter = 1000)
summary(LINE.out)
```

    ## 
    ## Iterations = 2001:3000
    ## Thinning interval = 1 
    ## Number of chains = 2 
    ## Sample size per chain = 1000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##        Mean      SD  Naive SE Time-series SE
    ## d[1] 0.0000 0.00000 0.0000000       0.000000
    ## d[2] 0.7844 0.03007 0.0006723       0.001225
    ## d[3] 1.4249 0.02825 0.0006318       0.001287
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##        2.5%    25%    50%   75%  97.5%
    ## d[1] 0.0000 0.0000 0.0000 0.000 0.0000
    ## d[2] 0.7238 0.7645 0.7841 0.806 0.8419
    ## d[3] 1.3692 1.4055 1.4255 1.443 1.4801

``` r
a <- summary(LINE.out)
a <- a$statistics[,"Mean"]
a <- round(a,2)
```

The above tables show that the log odds ratio is 0.78 for moderate
versus none and 1.42 for severe versus none.

## Next step

Next step will be to estimate the proportion of people in category two
rather than assume it is fixed. Can estimate this from the data or can
get input from subject-matter experts about the likely proportion. The
following estimates it from the data.

In order to prepare for this final decision we developed two models. One
model assumes (approach one) that the proportion of none-moderate which
are moderate lies, with equal probability, somewhere between 0 and 1 and
that the proportion of moderate-severe which are moderate lies between 0
and 1 (ie a uniform prior). A second model (approach two) assumes that
the proportion in each category is exchangeable between studies, using a
random effects model to estimate the proportion for those studies where
it is not recorded.

These models (which for speed only run for a small number of samples)
both give similar result.

### Approach One

Add an indicator variables for whether we want proportion in category
two of 1 and 2, or of 2 and 3. Where there is no missing data, the
proportion of two in 2/3 will be calculated, but this is not used in the
code.

In the first instance explore the use of a flat distribution for the
proportion in category two. This is run in both JAGS and OPENBUGS,
giving the same results.

``` r
list_data3 <- list_data2
list_data3[["pi2"]] <- NULL
list_data3[["coll12"]] <- grp_lbls2_col12
mod2 <- jags.model(file = "Supporting/FE_model_estimate_prop.txt",
                            data = list_data3, n.chains = 1, n.adapt = 1000)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 61
    ##    Unobserved stochastic nodes: 27
    ##    Total graph size: 1214
    ## 
    ## Initializing model

``` r
update(mod2, 1000)

data(LINE)
LINE$recompile()
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 5
    ##    Unobserved stochastic nodes: 3
    ##    Total graph size: 36
    ## 
    ## Initializing model

``` r
LINE.out <- coda.samples(mod2, c("d", "theta1", "theta2"),n.iter = 1000)
summary(LINE.out)
```

    ## 
    ## Iterations = 2001:3000
    ## Thinning interval = 1 
    ## Number of chains = 1 
    ## Sample size per chain = 1000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##           Mean      SD  Naive SE Time-series SE
    ## d[1]   0.00000 0.00000 0.0000000       0.000000
    ## d[2]   0.79737 0.02911 0.0009206       0.001695
    ## d[3]   1.41606 0.02964 0.0009374       0.001705
    ## theta1 0.06938 0.06183 0.0019551       0.004580
    ## theta2 0.69633 0.23303 0.0073691       0.013050
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##            2.5%     25%     50%     75%  97.5%
    ## d[1]   0.000000 0.00000 0.00000 0.00000 0.0000
    ## d[2]   0.737691 0.77703 0.79905 0.81737 0.8496
    ## d[3]   1.354987 1.39646 1.41718 1.43757 1.4705
    ## theta1 0.001181 0.02119 0.05075 0.09939 0.2297
    ## theta2 0.155408 0.56000 0.74264 0.88235 0.9875

``` r
inits <- function() {
  list(d = c(NA,-2,1),
       mu = rep(0, nrow(list_data3$r)))
}
list_data3_openbugs <- list_data3
list_data3_openbugs$coll12 <- as.integer(list_data3_openbugs$coll12)
mode2_bugs <- bugs(data = list_data3_openbugs, 
                   parameters.to.save = c("d", "theta1", "theta2"), n.iter = 1000,
                   inits = inits,
                   model.file = "FE_model_estimate_prop_openbugs.txt")
print(mode2_bugs, digits.summary = 3)
```

    ## Inference for Bugs model at "FE_model_estimate_prop_openbugs.txt", 
    ## Current: 3 chains, each with 1000 iterations (first 500 discarded)
    ## Cumulative: n.sims = 1500 iterations saved
    ##             mean    sd    2.5%     25%     50%     75%   97.5%  Rhat n.eff
    ## d[2]       0.798 0.029   0.740   0.779   0.799   0.817   0.854 1.008   250
    ## d[3]       1.416 0.029   1.359   1.397   1.416   1.437   1.476 1.007   450
    ## theta1     0.064 0.058   0.002   0.020   0.047   0.090   0.219 1.001  1500
    ## theta2     0.707 0.215   0.201   0.572   0.752   0.880   0.989 1.017   560
    ## deviance 414.108 7.323 401.147 409.175 413.600 418.500 430.052 1.001  1500
    ## 
    ## For each parameter, n.eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
    ## 
    ## DIC info (using the rule, pD = Dbar-Dhat)
    ## pD = 25.560 and DIC = 439.700
    ## DIC is an estimate of expected predictive error (lower deviance is better).

Show get similar result with Dirichlet prior, just a check on Dirichlet
coding. **Not sure that I (DM) have coded this correctly.**

``` r
list_data3 <- list_data2
list_data3[["pi2"]] <- NULL
list_data3[["coll12"]] <- grp_lbls2_col12
mod2 <- jags.model(file = "Supporting/FE_model_estimate_prop_dirich.txt",
                            data = list_data3, n.chains = 1, n.adapt = 1000)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 61
    ##    Unobserved stochastic nodes: 26
    ##    Total graph size: 1216
    ## 
    ## Initializing model

``` r
update(mod2, 1000)

data(LINE)
LINE$recompile()
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 5
    ##    Unobserved stochastic nodes: 3
    ##    Total graph size: 36
    ## 
    ## Initializing model

``` r
LINE.out <- coda.samples(mod2, c("d", "theta"),n.iter = 1000)
summary(LINE.out)
```

    ## 
    ## Iterations = 2001:3000
    ## Thinning interval = 1 
    ## Number of chains = 1 
    ## Sample size per chain = 1000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##             Mean      SD  Naive SE Time-series SE
    ## d[1]     0.00000 0.00000 0.0000000       0.000000
    ## d[2]     0.79674 0.02879 0.0009103       0.001524
    ## d[3]     1.41650 0.02887 0.0009131       0.001617
    ## theta[1] 0.23673 0.20172 0.0063789       0.034756
    ## theta[2] 0.04801 0.04508 0.0014257       0.005336
    ## theta[3] 0.71526 0.20977 0.0066334       0.036558
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##               2.5%     25%     50%    75%  97.5%
    ## d[1]     0.0000000 0.00000 0.00000 0.0000 0.0000
    ## d[2]     0.7391502 0.77875 0.79532 0.8155 0.8520
    ## d[3]     1.3559338 1.39817 1.41745 1.4360 1.4728
    ## theta[1] 0.0002239 0.07864 0.18998 0.3591 0.7302
    ## theta[2] 0.0024407 0.01368 0.03342 0.0686 0.1555
    ## theta[3] 0.1785696 0.59527 0.75872 0.8795 0.9890

### Approach Two

This is the modelling where we estimate the proportion in the moderate
category using the data on proportions from other studies.

``` r
pi2_choose <- pi2_for_mdl$`w/a`
pi2_choose12 <- rowSums(pi2_choose[, 1:2])
pi2_choose23 <- rowSums(pi2_choose[, 2:3])

list_data4 <- list_data3
list_data4[["n_complete2a"]]   <- pi2_choose[,2]
list_data4[["n_complete2b"]]   <- pi2_choose[,2]
list_data4[["n_complete12"]] <- pi2_choose12
list_data4[["n_complete23"]] <- pi2_choose23
list_data4[["ns_complete"]] <- length(pi2_choose12)


mod3 <- jags.model(file = "Supporting/FE_model_estimate_prop2.txt",
                            data = list_data4, n.chains = 2, n.adapt = 1000)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 91
    ##    Unobserved stochastic nodes: 105
    ##    Total graph size: 1547
    ## 
    ## Initializing model

``` r
update(mod3, 1000)
# mcmcplot(as.mcmc(mod3))
data(LINE)
LINE$recompile()
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 5
    ##    Unobserved stochastic nodes: 3
    ##    Total graph size: 36
    ## 
    ## Initializing model

``` r
LINE.out <- coda.samples(mod3, c("d", "theta1_prop", "theta2_prop"),n.iter = 1000)
summary(LINE.out)
```

    ## 
    ## Iterations = 2001:3000
    ## Thinning interval = 1 
    ## Number of chains = 2 
    ## Sample size per chain = 1000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##               Mean      SD  Naive SE Time-series SE
    ## d[1]        0.0000 0.00000 0.0000000       0.000000
    ## d[2]        0.7918 0.03032 0.0006779       0.001312
    ## d[3]        1.4176 0.02891 0.0006465       0.001219
    ## theta1_prop 0.2711 0.03922 0.0008770       0.001328
    ## theta2_prop 0.5957 0.03697 0.0008266       0.001323
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##               2.5%    25%    50%    75%  97.5%
    ## d[1]        0.0000 0.0000 0.0000 0.0000 0.0000
    ## d[2]        0.7326 0.7719 0.7917 0.8116 0.8538
    ## d[3]        1.3622 1.3995 1.4179 1.4358 1.4749
    ## theta1_prop 0.1939 0.2454 0.2706 0.2974 0.3478
    ## theta2_prop 0.5226 0.5711 0.5955 0.6198 0.6678

# Random effects model - RE for effects but not for proportions

Return to approach 1, but re-run the model assuming random effects model
for the effect of malnutrition on mortality.

``` r
mode3_bugs <- bugs(data = list_data3_openbugs, 
                   parameters.to.save = c("d", "theta1", "theta2"), n.iter = 10000,
                   inits = inits,
                   model.file = "RE_model_estimate_prop_openbugs.txt")
print(mode3_bugs, digits.summary = 3)
```

    ## Inference for Bugs model at "RE_model_estimate_prop_openbugs.txt", 
    ## Current: 3 chains, each with 10000 iterations (first 5000 discarded)
    ## Cumulative: n.sims = 15000 iterations saved
    ##             mean     sd    2.5%     25%     50%     75%   97.5%  Rhat n.eff
    ## d[2]       0.711  0.114   0.487   0.635   0.709   0.783   0.944 1.001 12000
    ## d[3]       1.508  0.113   1.306   1.429   1.501   1.577   1.753 1.001 15000
    ## theta1     0.136  0.129   0.003   0.041   0.097   0.191   0.480 1.001  4000
    ## theta2     0.642  0.242   0.114   0.475   0.680   0.843   0.984 1.002  1600
    ## deviance 364.435 11.458 344.100 356.200 363.900 371.700 388.700 1.003   950
    ## 
    ## For each parameter, n.eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
    ## 
    ## DIC info (using the rule, pD = Dbar-Dhat)
    ## pD = 43.390 and DIC = 407.800
    ## DIC is an estimate of expected predictive error (lower deviance is better).

# RE for effects and for proportions

Combine approaches 1 and 3, where assuming random effects for effect of
malnutrition, and random effects for proportion in each malnutrition
category.

Uncomment `mcmcplot(as.mcmc(mode4_bugs))` to see model diagnostics.Look
reasonable to me without too much autocorrelation. Chains appear to
converge.

``` r
list_data4_openbugs <- c(list_data3_openbugs, list_data4[c("n_complete2a", "n_complete2b", 
                                                         "n_complete12", "n_complete23", 
                                                         "ns_complete")])
mode4_bugs <- bugs(data = list_data4_openbugs, 
                   parameters.to.save = c("d", "theta1_prop", "theta2_prop"), n.iter = 10000,
                   inits = inits,
                   model.file = "RE_model_RE_estimate_prop_openbugs.txt")
# mcmcplot(as.mcmc(mode4_bugs))
print(mode4_bugs, digits.summary = 3)
```

    ## Inference for Bugs model at "RE_model_RE_estimate_prop_openbugs.txt", 
    ## Current: 3 chains, each with 10000 iterations (first 5000 discarded)
    ## Cumulative: n.sims = 15000 iterations saved
    ##                mean     sd    2.5%     25%     50%     75%   97.5%  Rhat n.eff
    ## d[2]          0.698  0.120   0.459   0.621   0.698   0.776   0.938 1.001  6000
    ## d[3]          1.533  0.116   1.320   1.454   1.528   1.605   1.779 1.001  7700
    ## theta1_prop   0.287  0.039   0.213   0.261   0.286   0.311   0.368 1.001  3900
    ## theta2_prop   0.595  0.040   0.513   0.569   0.595   0.620   0.671 1.001 15000
    ## deviance    599.929 14.036 574.400 590.300 599.200 608.800 629.600 1.001  5200
    ## 
    ## For each parameter, n.eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
    ## 
    ## DIC info (using the rule, pD = Dbar-Dhat)
    ## pD = 73.830 and DIC = 673.800
    ## DIC is an estimate of expected predictive error (lower deviance is better).
