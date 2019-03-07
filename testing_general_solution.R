rm(list = ls())

devtools::install_github("donaldRwilliams/BGGM")




z <- BGGM::ptsd[,1:4] +1


tz=matrix(0,d,6)x <- rep(1, nrow(z))

d=4

k=1

b<-matrix(0,(d*k))

s= diag(d)

tz=matrix(0,d,6)

tz[,1]=-Inf; tz[,2]=0;  tz[,6]=Inf

zstar=matrix(0,nrow(z),d)



tz[,1]=-Inf; tz[,2]=0;  tz[,6]=Inf

zstar=matrix(0,nrow(z),d)



ctz=tz

acc1=acc2=acctot=0;



thresh <- seq(min(z[,1]) +1, max(z[,1]) - 1)


for(p in 1:d){
  for(t in 1:length(thresh)){
  tz[p,thresh[t] + 1]=qnorm(sum(z[,p]<= thresh[t])/nrow(z), mean=-qnorm(sum(z[,p]==1)/nrow(z),mean=0,sd=1),sd=1)
  }
}


Psi <- diag(d)


tz
r <- NA
r_t <- NA

rm(list = ls())

cor_save <- matrix(0, nrow = 5000, ncol = .5 * d * (d-1))
thresh_save <- matrix(0, nrow = 5000, ncol = 3)



ordinal_gibbs <- function(z, delta, MH_step, samples){
  
  x <- rep(1, nrow(z))
 
  thresh <- seq(min(z[,1]) + 1, max(z[,1]) - 1)
  
  d <- ncol(z)
  
  Psi <- diag(d)
  
  k <- 1
  
  b <-matrix(0,(d*k))
  
    s <- cor(z) 
  # diag(d)
  
  tz <-  matrix(0,d, length(unique(z[,1]))  + 1 )
  
  zstar=matrix(0,nrow(z),d)
  
  cor_save <- matrix(0, nrow = samples, ncol = .5 * d * (d-1))
  
  for(p in 1:d){
    for(t in 1:length(thresh)){
      tz[p,thresh[t] + 1]=qnorm(sum(z[,p]<= thresh[t]) / nrow(z), mean=-qnorm(sum(z[,p]==1)/nrow(z),mean=0,sd=1),sd=1)
    }
  }
  
  tz[,1] <- - Inf 
  tz[,2] <- 0
  tz[,length(unique(z[,1]))  + 1]  <- Inf
  ctz <- tz
mm2 = matrix(0, nrow(z), d)
for(i in 2:samples){
  
  bb= matrix(b,k,d)
  
  m=x%*%bb

  for(j in 1:d){
  
  mm = t(m[,j] + t(s[j,-j])%*%solve(s[-j,-j])%*% t((zstar[,-j]-m[,-j])))
  
  mm2[,j] = t(m[,j] + t(s[j,-j])%*%solve(s[-j,-j])%*% t((zstar[,-j]-m[,-j])))
  
  ss=s[j,j] - t(s[j,-j])%*%solve(s[-j,-j])%*%s[j,-j]
  
  zstar[,j]= qnorm(runif(nrow(z), min= pnorm(tz[j,z[,j]], mm , sqrt(ss)),
                                 
                          max=pnorm(tz[j,(z[,j] + 1)], mm, sqrt(ss))),
                          mean= mm,sd=sqrt(ss))
}

R <- NA
MH_step = 0.01
for(p in 1:d){
  for(t in 1:length(thresh)){
  
 repeat {
   ctz[p, thresh[t] + 1]= qnorm(runif(1,  
                                
                                min = pnorm(0, mean=tz[p, thresh[t] + 1], sd = MH_step), max = 1), 
                                
                                mean=tz[p, thresh[t] + 1], sd = MH_step)
  
  R_top <- (pnorm(ctz[p,z[,p]+1] - m[,p], mean = 0, sd = 1) - pnorm(ctz[p,z[,p]] - m[,p], mean = 0 , sd = 1))
  
  R_bottom <- (pnorm(tz[p, z[,p] + 1]-m[,p],0,1) - pnorm(tz[p,z[,p]]-m[,p],0,1))
  
  R <- sum(log((R_top / R_bottom))) + log((1-pnorm(- tz[p, thresh[t] + 1] / MH_step,0,1))/(1-pnorm(-ctz[p, thresh[t] + 1]/ MH_step,0,1)))             
  
  if(!is.nan(R)) break
  }
  
    if(R > log(runif(1, 0, 1))){

    tz[p,thresh[t] + 1] = ctz[p,thresh[t] + 1]
  }
    
    }
}


# round(solve(solve(cor(z) / length(x))) , 4) == round(  solve(solve(cor(z))%x%(t(x)%*%x)), 4)

vb=solve(solve(s)%x%(t(x)%*%x))

mn=vb%*%(as.vector(t(x)%*%zstar%*%t(solve(s))))

b=mn+t(rnorm((d*k),0,1)%*%chol(vb))

e <-  zstar - t(b %*%  t(x))

v=t(e)%*%e

# S <- t(zstar) %*% zstar

inv <- rWishart(1, delta + nrow(z) - 1, solve(v + Psi, tol  = 1e-20))[,,1]

Psi <- rWishart(1, 1000 + ncol(z)  + delta - 2, solve(inv + diag(1000, d),  tol  = 1e-20))[,,1]

pcs  <- cov2cor(inv)
pcs <- pcs * -1
diag(pcs) <- 1

s <- corpcor::pcor2cor(pcs)

cor_save[i,] <- pcs[upper.tri(pcs)]
}

mat <- mat_sd <- matrix(0, p, p)
mat[upper.tri(mat)] <- round(colMeans(cor_save), 3)
mat_sd[upper.tri(mat_sd)] <- apply(cor_save, 2, sd)

list(mat = mat, mat_sd = mat_sd)

}


zstar

z <- BGGM::ptsd  + 1

b <- BDgraph::bdgraph(z, method = "gcgm")

round(cov2cor(b$K_hat),3)

system.time({
fit <- ordinal_gibbs(z,delta = 10, MH_step  = .05, samples = 500)
})


colMeans(fit)

apply(fit,2, sd)


mat <- matrix(0, 10, 10)
mat[upper.tri(mat)] <- round(colMeans(fit), 3)

mat



corpcor::cor2pcor(cor(z))

hist(fit[,3])

for(i in 1:7){
print(shapiro.test( fit[,i] )$p.value)
}

#corpcor::cor2pcor(cor(z))

colMeans(fit)

colMeans(cor_save)

# for(p in 1:d){  
#  if(r_t[p] > log(runif(1, 0, 1))){
#       for(t in 1:length(thresh)){
#       tz[p,thresh[t] + 1] = ctz[p,thresh[t] + 1]
#     }
#   }
# }



thresh_save

par(mfrow=c(1,3))
nrow(unique(thresh_save))


hist(thresh_save[seq(100, 5000, by = 1),3])

coda::acfplot(as.mcmc(thresh_save[seq(1, 5000, by = 2),2]))

cor(cor_save)

sqrt(1 / (221 - 8 - 3))


apply(cor_save, 2, sd)

hist(cor_save[,8])


colMeans(cor_save)



test <- polycor::polychor(z[,1], z[,5], std.err = T)
sqrt(test$var)

sqrt((1 - 0.03^2) / (221 - 8  -3))

pcs <- psych::polychoric(z)

pcs$rho
pcs$rho[upper.tri(pcs$rho)]

pcs <- corpcor::cor2pcor(pcs$rho)
pcs
pcs[upper.tri(pcs)]

pc_pearson <- cor(z)

plot(pcs$rho[upper.tri(pcs$rho)], colMeans(cor_save[seq(1000, 5000, 2),]), ylab = "Bayes", xlab = "MLE")

plot(pc_pearson[upper.tri(pc_pearson)], colMeans(cor_save), ylab = "Bayes", xlab = "MLE")

