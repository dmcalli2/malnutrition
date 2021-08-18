# runopenbugs
library(R2OpenBUGS)
library(tidyverse)

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

inits <- function() {
  list(d = c(NA,-2,1),
       mu = c(0,0,0,0,0,	0,0,0,0,0,	0,0,0,0,0,	0,0,0,0,0,	0,0,0,0))
}


mod1 <- bugs(data = list_data, inits = inits, n.iter = 1000, model.file = "FE_model.txt",
             parameters.to.save = "d")
summary(mod1)

mod2 <- bugs(data = list_data, inits = inits, n.iter = 1000, model.file = "RE_model.txt",
             parameters.to.save = "d")
summary(mod1)
