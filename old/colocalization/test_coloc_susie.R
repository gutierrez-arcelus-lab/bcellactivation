library(tidyverse)
library(coloc)
library(susieR)

data(coloc_test_data)
attach(coloc_test_data)

str(D3)

S3 <- runsusie(D3)
str(S3)
summary(S3)

S4 <- runsusie(D4)


susie.res <- coloc.susie(S3, S4)

print(susie.res$summary)


#Susie
fit3 <- susie_rss(D3$beta / sqrt(D3$varbeta), D3$LD, n = D3$N)

str(fit3)

fit3$lbf_variable[1:10, c(82, 103, 105), drop = FALSE]

fit3$sets
fit3$lbf




#coloc.susie
s1 <- S3
s2 <- S4
cs1 <- s1$sets
cs2 <- s2$sets

idx1 <- cs1$cs_index
idx2 <- cs2$cs_index
bf1 <- s1$lbf_variable[idx1, , drop = FALSE]
bf2 <- s2$lbf_variable[idx2, , drop = FALSE]

ret <- coloc.bf_bf(bf1, bf2)
## renumber index to match
ret$summary[,idx1:=cs1$cs_index[idx1]]
ret$summary[,idx2:=cs2$cs_index[idx2]]
ret
