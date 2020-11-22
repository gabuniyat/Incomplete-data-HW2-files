

p_i_squi=function(p,y,mu,sigm,lambda){
  
  numertr=p*dlnorm(y, meanlog = mu, sdlog = sigm, log = FALSE)
  denomintr=p*dlnorm(y, meanlog = mu, sdlog = sigm, log = FALSE)+(1-p)*dexp(y, rate = lambda, log = FALSE)
  p_i_sqi=numertr/denomintr
  
  return(p_i_sqi)
}


est_em=function(theta_0,y){


#Unpack the variables
p_0=theta_0[1]
mu_0=theta_0[2]
sigm_0=sqrt(theta_0[3])
lambda_0=theta_0[4]
n=length(y)

tolr=1000

LL<-c()
iter=0

while (tolr>=10^(-17)){
  
p_isq=p_i_squi(p_0,y,mu_0,sigm_0,lambda_0)

#Update
p=sum(p_isq)/n
mu=sum(p_isq*log(y))/sum(p_isq)
sigm_sq=sum(p_isq*(log(y)-mu)^2)/sum(p_isq)
lambda=sum(1-p_isq)/sum((1-p_isq)*y)

tolr=abs(p_0-p)+abs(mu_0-mu)+abs(sigm_0-sqrt(sigm_sq))+abs(lambda_0-lambda)

p_0=p
mu_0=mu
sigm_0=sqrt(sigm_sq)
lambda_0=lambda

print(tolr)

L=sum(log(p*dlnorm(y, meanlog = mu, sdlog = sqrt(sigm_sq), log = FALSE)+(1-p)*dexp(y, rate = lambda, log = FALSE)))

LL=append(LL,L)
iter=iter+1
}

return(list(LL=LL,tol_crit=tolr,p=p,mu=mu,sigm=sqrt(sigm_sq),lambda=lambda,iter=iter))

}

theta_0<-c(0.1,1,0.5^2,2)
y=dataex5

outpt=est_em(theta_0,y)

LL=outpt$LL
iter=outpt$iter
p=outpt$p
mu=outpt$mu
sigm=outpt$sigm
lambda=outpt$lambda
iter=outpt$iter

plot(LL[10:iter],main="Evolution of the value of the log-likelihood",
     xlab="Iteration",
     ylab="Value")

plot(LL[900:iter],main="Evolution of the value of the log-likelihood",
     xlab="Iteration",
     ylab="Value")

plot(LL[1000:iter],main="Evolution of the value of the log-likelihood",
     xlab="Iteration",
     ylab="Value")

hist(y, main="Histogram and density",
     xlab = "Values",
     ylab = "Density",
     freq = FALSE)
curve(p*dnorm(x, mu, sigm, log = FALSE) + (1-p)*dexp(x, rate = lambda, log = FALSE),  to = 120, from = 0, add = TRUE, lwd=2, col = "blue2")


hist(y, main="Histogram and density",
     xlab = "Values",
     ylab = "Density",
     freq = F, ylim = c(0,0.35))
curve(p*dnorm(x, mu, sigm, log=FALSE) + (1-p)*dexp(x, rate = lambda, log = FALSE),
      add = TRUE, lwd = 2, col = "blue2")