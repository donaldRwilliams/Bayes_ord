library(Rcpp)
library(BDgraph)

## note that int will take over in division
#z = na.omit(BGGM::bfi[1:1000,1:10])

z = as.matrix(BGGM::ptsd[,1:20] + 1)

sourceCpp("test.cpp")


system.time({
temp =  test(as.matrix(z), samples = 2500,  
             temp = max(z[,1]) - 1,
             int_mat = cor(z))
})





system.time({
  fit <- ordinal_gibbs(z,delta = 10, MH_step  = .01, samples = 2500)
})


system.time({
fit_bd <- bdgraph(z, method = "gcgm")
})

round(fit$mat[1,-1],3)
round(apply(temp$cor_save, c(1,2), mean)[1,-1], 3)


round(fit$mat_sd[1,-1],3)
round(apply(temp$cor_save, c(1,2), sd)[1,-1], 3)









