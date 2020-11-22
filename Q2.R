#Install the maximum likelihood package
install.packages(maxLik)

#Normal likelihood
dens_trunk_norm <- function(x,r,mu){
  
  r*dnorm(x, mu, sd = 1.5, log = TRUE)+(1-r)*pnorm(x, mu, sd = 1.5, log.p = TRUE)
  
}

loglik_trun_norm<-function(dataa,mu){
  
  x=dataa[,1]
  r=dataa[,2]
  
  nn=length(x)
  
  l_i<-numeric(nn)
  
  for (i in 1:nn){
    l_i[i]=dens_trunk_norm(x[i],r[i],mu)
  }
  
  l=sum(l_i)
  
  return(l)
  
}
  
require(maxLik)
mle <- maxLik(logLik = loglik_trun_norm, dataa = dataex2, start = 0)
summary(mle)

