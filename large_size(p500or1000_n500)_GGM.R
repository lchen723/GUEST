
library(glasso)
library(clime)
library(MASS)
library(huge)
library(GUEST)

n=500
q1 = 12
q2 = 988
p = q1+q2 #500/1000

#Data generation
G = XMRF.Sim(n , q1, model = "GGM", graph.type = "scale-free")
X=t(G$X)
Strue = diag(1,p,p)
Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)
X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
X=cbind(X,X1)

c=0.3
W=X+mvrnorm(n,rep(0,p),diag(c,p,p))

#Independent Data generation
Strue = diag(1,p,p)
X=mvrnorm(n,rep(0,p),diag(1,p,p))

c=0.3
W=X+mvrnorm(n,rep(0,p),diag(c,p,p))

#simulation 100 times
Theta.hat1=array(dim = c(p,p,1))
Strue_for_scale_free=array(dim = c(p,p,1))
Time=NULL
for (i in 1:10){
  #data generation
  cat("\r",round(i/100*100,2), '%     ')
  set.seed(i)
  #Data generation
  G = XMRF.Sim(n , q1, model = "GGM", graph.type = "scale-free")
  X=t(G$X)
  X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
  X=cbind(X,X1)
  Strue = diag(1,p,p)
  Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)
  Strue_for_scale_free[,,i]=Strue
  W=X+mvrnorm(n,rep(0,p),diag(c,p,p))

  #boost.graph
  t1=proc.time()
  result1<-boost.graph(data = W,thre = 0.2,ite1=1,ite2 = 0,ite3 = 0,rep = 1,sigma_e = c,cor=T)
  Theta.hat1[,,i]=result1$w
  t2=proc.time()
  t=t2-t1

  Time=c(Time, t[3])
  if (i==1)
    print(mean(Time))
}

#take average
Theta1=matrix(0,p,p)
for (i in 1:10){
  Theta1=Theta1+Theta.hat1[,,i]
}
Theta1=Theta1/10

round(Theta1[1:12,1:12],3)
Strue[1:12,1:12]
Theta1[which(Theta1<0.293)]=0

spe_sen_bias_kl(Theta =Theta1  ,Strue = Strue,p=1000)

theta.hat=Theta_p500_ind_03[1:12,1:12]
net = theta.hat
net = network::network(net, directed = FALSE)
network::network.vertex.names(net)=paste0("X",network::network.vertex.names(net))
graph = GGally::ggnet2(net,size=3,node.color = "lightgray",label=T,label.size = 3,mode = "circle")
graph

#for scale-free
round(Theta.hat1[1:12,1:12,1],3)
Strue_for_scale_free[1:12,1:12,1]
theta1=Theta.hat1[,,1]
theta1[which(theta1<0.306)]=0
spe_sen_bias_kl(Theta = theta1  ,Strue = Strue_for_scale_free[,,1],p=1000)

##p=500
#Theta_guest_p500_c01_free=theta1;Theta_guest_p500_c03_free=theta1

##p=1000
#Theta_guest_p1000_c01_free=theta1;Theta_guest_p1000_c03_free=theta1

shat=diag(0,p,p)
Strue_for_scale_free=array(dim = c(p,p,1))
Time=NULL
for (i in 1:10){
  set.seed(i)
  cat("\r",round(i/100*100,2), '%     ')
  G = XMRF.Sim(n , q1, model = "GGM", graph.type = "scale-free")
  Strue = diag(1,p,p)
  Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)
  Strue_for_scale_free[,,i]=Strue
  X=t(G$X)
  X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
  X=cbind(X,X1)
  
  ##independent
  #X=mvrnorm(n,rep(0,p),diag(1,p,p))
  
  #W=X+mvrnorm(n,rep(0,p),diag(c,p,p))

  t1=proc.time()
 
  #shat=shat+glasso(cov(W),rho=0.05)$wi
  
  #clime1=clime(cov(W),lambda = 0.05,sigma = TRUE)
  #shat=shat+matrix(unlist(clime1$Omegalist),nrow = p)
  
  #shat=shat+glasso(cov(W)-diag(c,p,p),rho=0.05)$wi
  
  #shat=shat+huge(W, nlambda = 11)$beta[[11]]
  
  shat=shat+space.joint(W, lam1 = 0.1)$ParCor 
 
  #shat=shat+QUIC(var(W), rho = 0.1)$X
  
  t2=proc.time()
  #if (i==1)
  #  print(t2-t1)
  
  t=t2-t1
  Time=c(Time, t[3])
  if (i==1)
    print(mean(Time))
  
}
shat=shat/10
shat=matrix(as.numeric(shat),nrow=p) #huge

round(Strue[1:12,1:12],3)
round(shat[1:12,1:12],3)
shat[which(abs(shat)<0.1)]=0

spe_sen_bias_kl(Theta =shat  ,Strue = Strue,p=1000)

theta.hat=shat_space_p500_c03_ind[1:12,1:12]
net = theta.hat
net = network::network(net, directed = FALSE)
network::network.vertex.names(net)=paste0("X",network::network.vertex.names(net))
graph = GGally::ggnet2(net,size=3,node.color = "lightgray",label=T,label.size = 3,mode = "circle")
graph

#Hub
#shat_wainwright_p500_c01_hub=shat;shat_wainwright_p500_c03_hub=shat
#shat_huge_p500_c01_hub=shat;shat_huge_p500_c03_hub=shat
#shat_space_p500_c01_hub=shat;shat_space_p500_c03_hub=shat
#shat_QUIC_p500_c01_hub=shat;shat_QUIC_p500_c03_hub=shat

#shat_huge_p1000_c01_hub=shat;shat_huge_p1000_c03_hub=shat
#shat_space_p1000_c01_hub=shat;shat_space_p1000_c03_hub=shat
#shat_QUIC_p1000_c01_hub=shat;shat_QUIC_p1000_c03_hub=shat

#Lattice
#shat_wainwright_p500_c01_lat=shat;shat_wainwright_p500_c03_lat=shat
#shat_huge_p500_c01_lat=shat;shat_huge_p500_c03_lat=shat
#shat_space_p500_c01_lat=shat;shat_space_p500_c03_lat=shat
#shat_QUIC_p500_c01_lat=shat;shat_QUIC_p500_c03_lat=shat

#shat_huge_p1000_c01_lat=shat;shat_huge_p1000_c03_lat=shat
#shat_space_p1000_c01_lat=shat;shat_space_p1000_c03_lat=shat
#shat_QUIC_p1000_c01_lat=shat;shat_QUIC_p1000_c03_lat=shat

#Independent
#shat_wainwright_p500_c01_ind=shat;shat_wainwright_p500_c03_ind=shat
#shat_huge_p500_c01_ind=shat;shat_huge_p500_c03_ind=shat
#shat_space_p500_c01_ind=shat;shat_space_p500_c03_ind=shat
#shat_QUIC_p500_c01_ind=shat;shat_QUIC_p500_c03_ind=shat

#shat_huge_p1000_c01_ind=shat;shat_huge_p1000_c03_ind=shat
#shat_space_p1000_c01_ind=shat;shat_space_p1000_c03_ind=shat
#shat_QUIC_p1000_c01_ind=shat;shat_QUIC_p1000_c03_ind=shat

##scale-free
#shat_glasso_p500_c01_free=shat;shat_glasso_p500_c03_free=shat
#shat_clime_p500_c01_free=shat;shat_clime_p500_c03_free=shat
#shat_wainwright_p500_c01_free=shat;shat_wainwright_p500_c03_free=shat
#shat_huge_p500_c01_free=shat;shat_huge_p500_c03_free=shat
#shat_space_p500_c01_free=shat;shat_space_p500_c03_free=shat
#shat_QUIC_p500_c01_free=shat;shat_QUIC_p500_c03_free=shat

#shat_glasso_p1000_c01_free=shat;shat_glasso_p1000_c03_free=shat
#shat_wainwright_p1000_c01_free=shat;shat_wainwright_p1000_c03_free=shat
#shat_huge_p1000_c01_free=shat;shat_huge_p1000_c03_free=shat
#shat_space_p1000_c01_free=shat;shat_space_p1000_c03_free=shat
#shat_QUIC_p1000_c01_free=shat;shat_QUIC_p1000_c03_free=shat

setwd("C:/Users/user/Desktop/R")
save.image(file = "revision_1008.RData")

#space
setwd("C:/Users/user/Desktop/R/revision/space/src")
system("R CMD SHLIB JSRM.c")
dyn.load('C:/Users/user/Desktop/R/revision/space/src/JSRM.dll')

#QUIC
setwd("C:/Users/user/Desktop/R/revision/QUIC/src")
system("R CMD SHLIB QUIC.cpp")
dyn.load('C:/Users/user/Desktop/R/revision/QUIC/src/QUIC.dll')

