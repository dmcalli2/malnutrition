##01_maln
library(tidyverse)
library(rjags)

maln <- read_csv("Data/Mortality_Numbers.csv")

## Run model with example data
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

mod_example <- jags.model(file = "Supporting/FE_model.txt",
                            data = list_data, n.chains = 2, n.adapt = 1000)
update(mod_example, 1000)

data(LINE)
LINE$recompile()
LINE.out <- coda.samples(mod_example, c("d"),n.iter = 2000)
summary(LINE.out)
