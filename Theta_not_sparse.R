
n=500
p=1000 #500/1000

## Not sparse
S=matrix(0.8,p,p)+diag(0.2,p,p)
mu=rep(0,p)

X=mvrnorm(n, mu, S)

c=0.3 #0.1/0.3
W=X+mvrnorm(n,rep(0,p),diag(c,p,p))

H=var(W)%*%(var(W)+diag(1,p,p))^(-2)

Theta.hat1=array(dim = c(p,p,5))
Time=NULL
for (i in 1:5){
  #data generation
  cat("\r",round(i/100*100,2), '%     ')
  set.seed(i)
  #Data generation
  X=mvrnorm(n, mu, S)
  
  W=X+mvrnorm(n,rep(0,p),diag(c,p,p))
  
  #boost.graph
  t1=proc.time()
  result1<-GUEST::boost.graph(data = W,thre = 0.2,ite1=1,ite2 = 0,ite3 = 0,rep = 1,sigma_e = c,cor=T)
  Theta.hat1[,,i]=result1$w
  t2=proc.time()
  t=t2-t1
  
  Time=c(Time, t[3])
  if (i==5)
    print(mean(Time))
}

#take average
Theta1=matrix(0,p,p)
for (i in 1:1){
  Theta1=Theta1+Theta.hat1[,,i]
}
Theta1=Theta1/1

round(Theta1[1:12,1:12],3)
S[1:12,1:12]
Theta1[which(Theta1<0.49)]=0

spe_sen_bias_kl(Theta =Theta1  ,Strue = S,p=1000)

#Theta_guest_p500_c01_notsparse=Theta1;Theta_guest_p500_c03_notsparse=Theta1
#Theta_guest_p1000_c01_notsparse=Theta1;Theta_guest_p1000_c03_notsparse=Theta1

shat=diag(0,p,p)
Time=NULL
for (i in 1:5){
  set.seed(i)
  cat("\r",round(i/100*100,2), '%     ')
  X=mvrnorm(n, mu, S)
  
  W=X+mvrnorm(n,rep(0,p),diag(c,p,p))
  
  t1=proc.time()
  
  #shat=shat+glasso(cov(W),rho=0.05)$wi
  
  #clime1=clime(cov(W),lambda = 0.05,sigma = TRUE)
  #shat=shat+matrix(unlist(clime1$Omegalist),nrow = p)
  
  #shat=shat+glasso(cov(W)-diag(c,p,p),rho=0.05)$wi
  
  #shat=shat+huge(W, nlambda = 11)$beta[[11]]
  
  shat=shat+space.joint(W, lam1 = 0.1)$ParCor 

  #shat=shat+QUIC(var(W), rho = 0.1)$X
  
  #shat=shat+var(W)%*%(var(W)+diag(1,p,p))^(-2)
  
  t2=proc.time()
  
  t=t2-t1
  Time=c(Time, t[3])
  if (i==5)
    print(mean(Time))
  
}
shat=shat/5
shat=matrix(as.numeric(shat),nrow=p) #huge

round(S[1:12,1:12],3)
round(shat[1:12,1:12],3)
shat[which(abs(shat)<0.001)]=0

spe_sen_bias_kl(Theta =shat  ,Strue = S,p=1000)

##p500
#shat_glasso_p500_c01_notsparse=shat;shat_glasso_p500_c03_notsparse=shat
#shat_clime_p500_c01_notsparse=shat;shat_clime_p500_c03_notsparse=shat
#shat_wainwright_p500_c01_notsparse=shat;shat_wainwright_p500_c03_notsparse=shat
#shat_huge_p500_c01_notsparse=shat;shat_huge_p500_c03_notsparse=shat
#shat_space_p500_c01_notsparse=shat;shat_space_p500_c03_notsparse=shat
#shat_QUIC_p500_c01_notsparse=shat;shat_QUIC_p500_c03_notsparse=shat
#shat_H_p500_c01_notsparse=shat;shat_H_p500_c03_notsparse=shat

##p1000
#shat_glasso_p1000_c01_notsparse=shat;shat_glasso_p1000_c03_notsparse=shat
#shat_wainwright_p1000_c01_notsparse=shat;shat_wainwright_p1000_c03_notsparse=shat
#shat_huge_p1000_c01_notsparse=shat;shat_huge_p1000_c03_notsparse=shat
#shat_space_p1000_c01_notsparse=shat;shat_space_p1000_c03_notsparse=shat
#shat_QUIC_p1000_c01_notsparse=shat;shat_QUIC_p1000_c03_notsparse=shat
#shat_H_p1000_c01_notsparse=shat;shat_H_p1000_c03_notsparse=shat

setwd("C:/Users/user/Desktop/R")
save.image(file = "Not_Sparse_1009.RData")
