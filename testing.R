#create variables and starting values

x=as.matrix(read.table("mvnprob.dat1")[1:100,1:7])
z=as.matrix(read.table("mvnprob.dat1")[1:221,8])


z <- BGGM::ptsd[,2:5]  + 1

x <- rep(1, 221)


d=4
k=1

b<-matrix(0,(d*k))

s=cs=diag(d)

tz=matrix(0,d,6)

tz[,1]=-Inf; tz[,2]=0;  tz[,6]=Inf

zstar=matrix(0,nrow(z),d)


tz[1,3]=qnorm(sum(z[,1]<=2)/nrow(z), mean=-qnorm(sum(z[,1]==1)/nrow(z),mean=0,sd=1),sd=1)
tz[1,4]=qnorm(sum(z[,1]<=3)/nrow(z), mean=-qnorm(sum(z[,1]==1)/nrow(z),mean=0,sd=1),
              sd=1)

tz[1,5]=qnorm(sum(z[,1]<=4)/nrow(z), mean=-qnorm(sum(z[,1]==1)/nrow(z),mean=0,sd=1),sd=1)



tz[2,3]=qnorm(sum(z[,2]<=2)/nrow(z), mean=-qnorm(sum(z[,2]==1)/nrow(z),mean=0,sd=1), sd=1)

tz[2,4]=qnorm(sum(z[,2]<=3)/nrow(z), mean=-qnorm(sum(z[,2]==1)/nrow(z),mean=0,sd=1), sd=1)
tz[2,5]=qnorm(sum(z[,2]<=4)/nrow(z), mean=-qnorm(sum(z[,2]==1)/nrow(z),mean=0,sd=1), sd=1)


tz[3,3]=qnorm(sum(z[,3]<=2)/nrow(z), mean=-qnorm(sum(z[,3]==1)/nrow(z),mean=0,sd=1), sd=1)
tz[3,4]=qnorm(sum(z[,3]<=3)/nrow(z), mean=-qnorm(sum(z[,3]==1)/nrow(z),mean=0,sd=1), sd=1)
tz[3,5]=qnorm(sum(z[,3]<=4)/nrow(z), mean=-qnorm(sum(z[,3]==1)/nrow(z),mean=0,sd=1), sd=1)



tz[4,3]=qnorm(sum(z[,4]<=2)/nrow(z), mean=-qnorm(sum(z[,4]==1)/nrow(z),mean=0,sd=1), sd=1)
tz[4,4]=qnorm(sum(z[,4]<=3)/nrow(z), mean=-qnorm(sum(z[,4]==1)/nrow(z),mean=0,sd=1), sd=1)
tz[4,5]=qnorm(sum(z[,4]<=4)/nrow(z), mean=-qnorm(sum(z[,4]==1)/nrow(z),mean=0,sd=1), sd=1)





thresh_mat <- matrix(0, 6000, 3)

ctz=tz
acc1=acc2=acctot=0;


s <- psych::polychoric(z)
s <- s$rho

cor_save <- matrix(0, nrow = 6000, ncol = 6)
Psi <- diag(4)

for(i in 2:6000){
  #draw latent data: one-iteration gibbs sampler for tmvn simulation
  bb=matrix(b,k,d)
  m=x%*%bb
  
 
   for(j in 1:d){
     
      mm= m[,j] + t(s[j,-j])%*%solve(s[-j,-j])%*%t((zstar[,-j]-m[,-j]))
      
      
      ss=s[j,j] - t(s[j,-j])%*%solve(s[-j,-j])%*%s[j,-j]
  
    
    
       zstar[,j]= qnorm(runif(nrow(z), min=pnorm(tz[j,z[,j]],mm , sqrt(ss)),
                        max=pnorm(tz[j,(z[,j] + 1)],mm,sqrt(ss))),
                        mean=mm,sd=sqrt(ss))
  }
  
 
 for(l in 3:5){
  ctz[1,l]= qnorm(runif(1,  min=pnorm(0, mean=tz[1,l],sd=.01),max=1), mean=tz[1,l],sd=.01)
 
  
 }
   
  r=as.matrix((pnorm(ctz[1,z[,1]+1]-m[,1],0,1)
                -pnorm(ctz[1,z[,1]]-m[,1],0,1))
             /
                (pnorm(tz[1,z[,1]+1]-m[,1],0,1)
                 -pnorm(tz[1,z[,1]]-m[,1],0,1)))
  r1 <- NA
  
  for(p in 3:5){
  
  r1[p] = t(log(r))%*%matrix(1,nrow(z)) + log((1-pnorm(-tz[1,p]/.01,0,1))/(1-pnorm(-ctz[1,p]/.01,0,1)))
  
  }
  
  
r1 <- sum(r1[3:5])
  
 for(t in 3:5){
    ctz[2,t]= qnorm(runif(1,  min=pnorm(0, mean=tz[2,t],sd=.01),max=1), mean=tz[2,t],sd=.01)
  }
  
  
  r=as.matrix((pnorm(ctz[2,z[,2]+1]-m[,2],0,1)
               -pnorm(ctz[2,z[,2]]-m[,2],0,1))
              /
                (pnorm(tz[2,z[,2]+1]-m[,2],0,1)
                 -pnorm(tz[2,z[,2]]-m[,2],0,1)))
  
  
  r2 <- NA
  for(g in 3:5){
    
    r2[g] = t(log(r))%*%matrix(1,nrow(z)) + log((1-pnorm(-tz[2,g]/.01,0,1))/(1-pnorm(-ctz[2,g]/.01,0,1)))
    
  }
  
  r2 <- sum(r2[3:5])
  
  
  for(t in 3:5){
    ctz[3,t]= qnorm(runif(1,  min=pnorm(0, mean=tz[3,t],sd=.01),max=1), mean=tz[3,t],sd=.01)
  }
  
  
  r=as.matrix((pnorm(ctz[3,z[,3]+1]-m[,3],0,1)
               -pnorm(ctz[3,z[,3]]-m[,3],0,1))
              /
                (pnorm(tz[3,z[,3]+1]-m[,3],0,1)
                 -pnorm(tz[3,z[,3]]-m[,3],0,1)))
 
  for(g in 3:5){
    
    r3[g] = t(log(r))%*%matrix(1,nrow(z)) + log((1-pnorm(-tz[3,g]/.01,0,1))/(1-pnorm(-ctz[3,g]/.01,0,1)))
    
  }
  r3 <- sum(r3[3:5])
  
  
  
  
  
  
  for(t in 3:5){
    ctz[4,t]= qnorm(runif(1,  min=pnorm(0, mean=tz[4,t],sd=.01),max=1), mean=tz[4,t],sd=.01)
  }
  
  
  r=as.matrix((pnorm(ctz[4,z[,4]+1]-m[,4],0,1)
               -pnorm(ctz[4,z[,4]]-m[,4],0,1))
              /
                (pnorm(tz[4,z[,4]+1]-m[,4],0,1)
                 -pnorm(tz[4,z[,4]]-m[,4],0,1)))
  
  
  r4 <- NA
  for(g in 3:5){
    
    r4[g] = t(log(r))%*%matrix(1,nrow(z)) + log((1-pnorm(-tz[4,g]/.01,0,1))/(1-pnorm(-ctz[4,g]/.01,0,1)))
    
  }
  r4 <- sum(r4[3:5])
  
  # sum of sums ?
  if(sum(c(r1, r2, r3, r4)) > log(runif(1,0,1)) ){
    tz[1,3]=ctz[1,3]; tz[1,4]=ctz[1,4]; tz[1,5]=ctz[1,5]; 
    tz[2,3]=ctz[2,3]; tz[2,4]=ctz[2,4]; tz[2,5]=ctz[2,5];  
    tz[3,3]=ctz[3,3]; tz[3,4]=ctz[3,4]; tz[3,5]=ctz[3,5]; 
    tz[4,3]=ctz[4,3]; tz[4,4]=ctz[4,4]; tz[4,5]=ctz[4,5]; acc1=acc1+1}
  
  
  vb=solve(solve(s)%x%(t(x)%*%x))
  
  set.seed(1)
  MASS::mvrnorm(1, mn, vb)
  
  mn=vb%*%(as.vector(t(x)%*%zstar%*%t(solve(s))))
  
  set.seed(1)
  b=mn+t(rnorm((d*k),0,1)%*%chol(vb))
  
  #use metropolis-hastings sampling to draw sigma
  
  #e=matrix((as.vector(zstar)-(diag(d)%x%x%*%b)),nrow(z),d)
  e <-  zstar - t(b %*%  t(x))  
  v=t(e)%*%e
  
  #like=-.5*((d+nrow(z)+1)*log(det(s)) + sum(diag(v%*%solve(s))))
  
  
  s <- cov2cor(rWishart(1, 221, v + Psi)[,,1])
  
  Psi <- rWishart(1, 3 + delta + p - 2, solve(s + diag(4),  tol  = 1e-20))[,,1]
  #cs[upper.tri(cs)] <- s[upper.tri(s)] + rnorm(6,mean=0,sd=.025)
  #
  #
  #cs <- BGGM:::symmteric_mat(cs)

  
  # sigma <- rWishart(1, delta + n - 1, s + Psi)[,,1]
     # 
       # 
       # Psi <- rWishart(1, 1 + delta + p - 2, solve(sigma + diag(4),  tol  = 1e-20))[,,1]
       # 
       # cs <- cov2cor(sigma)
  
 # if( det(cs) > 0 ){
 # 
 #    cslike=-.5*((d+nrow(z)+1)*log(det(cs)) + sum(diag(v%*%solve(cs))))
 # 
 #    if((cslike-like)>log(runif(1,0,1))){
 # 
 #      s = cs; acctot=acctot+1
 # 
 #    }
 #  }
  cor_save[i,] <- s[upper.tri(s)]
  thresh_mat[i,] <- tz[3,3:5]
  
}


Psi <- diag(4)
cor(zstar)

delta = 1
n = 221

samps <- unique(cor_save)


samps

colMeans(cor_save)

colMeans( thresh_mat)

cor_save
coda::traceplot(coda::as.mcmc( cor_save[2:6000 ,1]))


sd(samps[,3])
cor_save[1000:1500,]


colMeans(cor_save[1000:6000,])

hist(samps[,2])

Psi <- diag(3)

 colMeans(samps)


apply(cor_save[1000:6000,], 2, sd)



compare <- psych::polychoric(z)
compare$rho[upper.tri(compare$rho)]
compare$tau

test <- polycor::polychor(z[,1], z[,3], std.err = T)
sqrt(test$var)

colMeans(thresholds)

sd(samps[,3])


sd(samps[,3])

sd(samps[,2])

thresholds <- unique(thresh_mat)
thresholds


hist(samps[,1])
test <-  unique(cor_save)

test[[2]]
length(test)
mat_res <- matrix(0, nrow = length(test), 3)
for(i in 2:length(test)){
mat_res[i,] <- test[[i]][upper.tri(test[[i]])  ]
  
  
  
  
}

colMeans(mat_res)

colMeans(thresholds)

test <- polycor::polychor(z[,2], z[,3], std.err = T)

test$rho
sqrt(test$var)
