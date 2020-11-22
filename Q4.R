#Install the maximum likelihood package
install.packages(maxLik)

require(maxLik)


p_i_beta=function(x,beta){
  
  #Unpack the parameters
  beta_0=beta[1]
  beta_1=beta[1]
  
  pi_beta=exp(beta_0+x*beta_1)/(1+exp(beta_0+x*beta_1))
  
return(pi_beta) 
  
}


#M-Step: Update the parameters

ll=function(beta,dataa){
  
  x<-dataa[,1]
  y<-dataa[,2]
  
  beta_0=beta[1]
  beta_1=beta[2]
  
  n=length(y)
  l_i<-numeric(n)
  
  for (i in 1:n){
    el0=log(1+exp(beta_0+beta_1*x[i]))
    l_i[i]=-el0+y[i]*(beta_0+beta_1*x[i])
    
  }
  
  ll=sum(l_i)
  
  return(ll)
  
}

e_m=function(dat_a,beta0){


x<-dat_a[,1]
y<-dat_a[,2]
dum=is.na(y)
pp=p_i_beta(x,beta0)
y_miss0=pp[is.na(y)==TRUE]
tol_r=1000
iter_i=0
LL<-c()

while (tol_r>=10^(-17)){

#E-Step:Replace the missing values by the expectations
y[dum==TRUE]=y_miss0

dataa=cbind(x,y)

#M-Step: Update the parameters
mle <- maxLik(logLik = ll, dataa = dataa, start = c(beta0[1], beta0[2]))
beta1=mle$estimate
pp=p_i_beta(x,beta0)
y_miss1=pp[dum==TRUE]

tol_r=sum(abs(beta0-beta1))
y_miss0=y_miss1
beta0=beta1

LL=append(LL,mle$maximum)

iter_i=iter_i+1
print(tol_r)

}



return (list(mle=mle,iter_i=iter_i,beta=beta1,LL=LL))

}

out=e_m(dataex4,c(0,0))
