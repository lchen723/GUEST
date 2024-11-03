
library(glasso)
library(clime)
library(MASS)
library(huge)
library(GUEST)

#space
dyn.load('C:/Users/user/Desktop/R/revision/space/src/JSRM.dll')

#QUIC
dyn.load('C:/Users/user/Desktop/R/revision/QUIC/src/QUIC.dll')

n=500
q1 = 12
q2 = 988
p = q1+q2

##Data generation
#dependent case
G = XMRF.Sim(n , q1, model = "LPGM", graph.type = "lattice")
X=G$X
Strue = diag(1,p,p)
Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)

#(lambda,pi)=(0.5,0.5)/(0.8,0.5)
Z=matrix(rpois(X,lambda = 0.5),nrow = q1)
W=matrix(rbinom(X,1,0.5),nrow=q1)
X_star=X+Z-W
X_star=t(X_star)

X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
X_star=cbind(X_star,X1)


##independent case
Strue = diag(1,p,p)
X=diag(1,q1,n)
Y=matrix(rpois(X,lambda = 0.5),nrow = q1)
#(lambda,pi)=(0.5,0.5)/(0.8,0.5)
Z=matrix(rpois(Y,lambda = 0.5),nrow = q1)
W=matrix(rbinom(Y,1,0.5),nrow=q1)
X_star=Y+Z-W
X_star=t(X_star)

X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
X_star=cbind(X_star,X1)


Theta.hat1=array(dim = c(p,p,1))
for (i in 1:10){
  #data generation
  cat("\r",round(i/100*100,2), '%     ')
  set.seed(i)
  G = XMRF.Sim(n , q1, model = "LPGM", graph.type = "scale-free")
  Strue = diag(1,p,p)
  Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)
  X=G$X
  
  #(lambda,pi)=(0.5,0.5)/(0.8,0.5)
  Z=matrix(rpois(X,lambda = 0.8),nrow = q1)
  W=matrix(rbinom(X,1,0.5),nrow=q1)
  X_star=X+Z-W
  X_star=t(X_star)
  
  ##independent case
  #Strue = diag(1,p,p)
  #X=diag(1,q1,n)
  #Y=matrix(rpois(X,lambda = 0.5),nrow = q1)
  
  ##(lambda,pi)=(0.5,0.5)/(0.8,0.5)
  #Z=matrix(rpois(Y,lambda = 0.8),nrow = q1)
  #W=matrix(rbinom(Y,1,0.5),nrow=q1)
  #X_star=Y+Z-W
  #X_star=t(X_star)
  
  X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
  X_star=cbind(X_star,X1)
  
  #boost.graph
  t1=proc.time()
  result1<-boost.graph(data = X_star,thre = 0.3, inc = 0.0001, ite1=0,ite2 = 0,ite3 = 5,rep = 1,lambda = 0.8,pi = 0.5,cor=T)
  Theta.hat1[,,i]=result1$w
  t2=proc.time()
  t=t2-t1
  
  if (i==1)
    print(t)
}

#take average
Theta1=matrix(0,p,p)
for (i in 1:10){
  Theta1=Theta1+Theta.hat1[,,i]
}
Theta1=Theta1/10

round(Theta1[1:12,1:12],3)
Strue[1:12,1:12]
Theta1[which(Theta1<0.08)]=0

spe_sen_bias_kl(Theta =Theta1  ,Strue = Strue,p=1000)

#Theta_guest_p500_lambda05_free=Theta1;Theta_guest_p500_lambda08_free=Theta1
#Theta_guest_p1000_lambda05_free=Theta1;Theta_guest_p1000_lambda08_free=Theta1

theta.hat=shat[1:12,1:12]
net = theta.hat
net = network::network(net, directed = FALSE)
network::network.vertex.names(net)=paste0("X",network::network.vertex.names(net))
graph = GGally::ggnet2(net,size=3,node.color = "lightgray",label=T,label.size = 3,mode = "circle")
graph



shat=diag(0,p,p)
#Strue_for_scale_free=array(dim = c(p,p,5))
Time=NULL
for (i in 1:10){
  set.seed(i)
  cat("\r",round(i/100*100,2), '%     ')
  set.seed(i)
  G = XMRF.Sim(n , q1, model = "LPGM", graph.type = "scale-free")
  Strue = diag(1,p,p)
  Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)
  X=G$X
  
  #(lambda,pi)=(0.5,0.5)/(0.8,0.5)
  Z=matrix(rpois(X,lambda = 0.8),nrow = q1)
  W=matrix(rbinom(X,1,0.5),nrow=q1)
  X_star=X+Z-W
  X_star=t(X_star)
  
  #independent case
  #Strue = diag(1,p,p)
  #X=diag(1,q1,n)
  #Y=matrix(rpois(X,lambda = 0.5),nrow = q1)
  
  #(lambda,pi)=(0.5,0.5)/(0.8,0.5)
  #Z=matrix(rpois(Y,lambda = 0.8),nrow = q1)
  #W=matrix(rbinom(Y,1,0.5),nrow=q1)
  #X_star=Y+Z-W
  #X_star=t(X_star)
  
  X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
  X_star=cbind(X_star,X1)
  
  t1=proc.time()
  
  shat=shat+glasso(cov(X_star),rho=0.05)$wi
  
  #clime1=clime(cov(X_star),lambda = 0.05,sigma = TRUE)
  #shat=shat+matrix(unlist(clime1$Omegalist),nrow = p)
  
  #shat=shat+glasso(cov(X_star)-diag(c,p,p),rho=0.05)$wi
  
  #shat=shat+huge(X_star, nlambda = 11)$beta[[11]]
  
  #shat=shat+space.joint(X_star, lam1 = 0.1)$ParCor 
  
  shat=shat+QUIC(var(X_star), rho = 0.1)$X
  
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
shat[which(abs(shat)<0.02)]=0

spe_sen_bias_kl(Theta =shat  ,Strue = Strue,p=1000)

##Hub
#shat_wainwright_p500_lambda05_hub=shat;shat_wainwright_p500_lambda08_hub=shat
#shat_huge_p500_lambda05_hub=shat;shat_huge_p500_lambda08_hub=shat
#shat_space_p500_lambda05_hub=shat;shat_space_p500_lambda08_hub=shat
#shat_QUIC_p500_lambda05_hub=shat;shat_QUIC_p500_lambda08_hub=shat

#shat_huge_p1000_lambda05_hub=shat;shat_huge_p1000_lambda08_hub=shat
#shat_space_p1000_lambda05_hub=shat;shat_space_p1000_lambda08_hub=shat
#shat_QUIC_p1000_lambda05_hub=shat;shat_QUIC_p1000_lambda08_hub=shat

##Lattice
#shat_wainwright_p500_lambda05_lat=shat;shat_wainwright_p500_lambda08_lat=shat
#shat_huge_p500_lambda05_lat=shat;shat_huge_p500_lambda08_lat=shat
#shat_space_p500_lambda05_lat=shat;shat_space_p500_lambda08_lat=shat
#shat_QUIC_p500_lambda05_lat=shat;shat_QUIC_p500_lambda08_lat=shat

#shat_huge_p1000_lambda05_lat=shat;shat_huge_p1000_lambda08_lat=shat
#shat_space_p1000_lambda05_lat=shat;shat_space_p1000_lambda08_lat=shat
#shat_QUIC_p1000_lambda05_lat=shat;shat_QUIC_p1000_lambda08_lat=shat

##Independent
#shat_wainwright_p500_lambda05_ind=shat;shat_wainwright_p500_lambda08_ind=shat
#shat_huge_p500_lambda05_ind=shat;shat_huge_p500_lambda08_ind=shat
#shat_space_p500_lambda05_ind=shat;shat_space_p500_lambda08_ind=shat
#shat_QUIC_p500_lambda05_ind=shat;shat_QUIC_p500_lambda08_ind=shat

#shat_huge_p1000_lambda05_ind=shat;shat_huge_p1000_lambda08_ind=shat
#shat_space_p1000_lambda05_ind=shat;shat_space_p1000_lambda08_ind=shat
#shat_QUIC_p1000_lambda05_ind=shat;shat_QUIC_p1000_lambda08_ind=shat

##Scale-free
#shat_glasso_p500_lambda05_free=shat;shat_glasso_p500_lambda08_free=shat
#shat_clime_p500_lambda05_free=shat;shat_clime_p500_lambda08_free=shat
#shat_huge_p500_lambda05_free=shat;shat_huge_p500_lambda08_free=shat
#shat_space_p500_lambda05_free=shat;shat_space_p500_lambda08_free=shat
#shat_QUIC_p500_lambda05_free=shat;shat_QUIC_p500_lambda08_free=shat

#shat_glasso_p1000_lambda05_free=shat;shat_glasso_p1000_lambda08_free=shat
#shat_huge_p1000_lambda05_free=shat;shat_huge_p1000_lambda08_free=shat
#shat_space_p1000_lambda05_free=shat;shat_space_p1000_lambda08_free=shat
#shat_QUIC_p1000_lambda05_free=shat;shat_QUIC_p1000_lambda08_free=shat

setwd("C:/Users/user/Desktop/R")
save.image(file = "revision_1012.RData")

theta.hat=shat[1:12,1:12]
net = theta.hat
net = network::network(net, directed = FALSE)
network::network.vertex.names(net)=paste0("X",network::network.vertex.names(net))
graph = GGally::ggnet2(net,size=3,node.color = "lightgray",label=T,label.size = 3,mode = "circle")
graph


##hub
#Theta_p500_05=Theta1 (thre=0.3, ite3=1)
#Theta_p500_08=Theta1
#Theta_p1000_05=Theta1 (thre=0.3, inc=0.0001, ite3=5)
#Theta_p1000_08=Theta1

##lattice
#Theta_p500_05_lat=Theta1 (thre=0.3, inc=0.0001, ite3=5)
#Theta_p500_08_lat=Theta1
#Theta_p1000_05_lat=Theta1 (thre=0.3, inc=0.0001, ite3=5)
#Theta_p1000_08_lat=Theta1

##independent
#Theta_p500_05_int=Theta1 (thre=0.3, inc=0.0001, ite3=5)
#Theta_p500_08_int=Theta1
#Theta_p1000_05_int=Theta1 (thre=0.3, inc=0.0001, ite3=5)
#Theta_p1000_08_int=Theta1

