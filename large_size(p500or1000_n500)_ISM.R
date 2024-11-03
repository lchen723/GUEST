
n=500
q1 = 12
q2 = 988
p = q1+q2

##Data generation
#dependent case
G = XMRF.Sim(n , q1, model = "ISM", graph.type = "lattice")
X=t(G$X)
Strue = diag(1,p,p)
Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)
for (j in 1:q2){
 X1=rbinom(n,1,0.9) 
 X=cbind(X,X1)
}

#independent case
Strue = diag(1,p,p)
X=NULL
for (j in 1:p){
  X1=rbinom(n,1,0.9) 
  X=cbind(X,X1)
}

W=data.frame()
for (j in 1:p){
  S=rbinom(n,1,0.85) #s=0.85/0.9
  W[1:n,j]=S*X[,j]+(1-S)*(1-X[,j])
}
W=matrix(unlist(W),ncol = p)


#simulation
Theta.hat1=array(dim = c(p,p,1))
for (i in 1:1:5){
  #data generation
  cat("\r",round(i/100*100,2), '%     ')
  set.seed(1)
  #Data generation
  G = XMRF.Sim(n , q1, model = "ISM", graph.type = "scale-free")
  Strue = diag(1,p,p)
  Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)
  X=t(G$X)
  for (j in 1:q2){
    X1=rbinom(n,1,0.9) 
    X=cbind(X,X1)
  }
  
  ##independent case
  #X=NULL
  #for (j in 1:p){
  #  X1=rbinom(n,1,0.9) 
  #  X=cbind(X,X1)
  #}
  
  W=data.frame()
  for (j in 1:p){
    S=rbinom(n,1,0.9) #s=0.85/0.9
    W[1:n,j]=S*X[,j]+(1-S)*(1-X[,j])
  }
  W=matrix(unlist(W),ncol = p)
  
  #boost.graph
  t1=proc.time()
  result1<-boost.graph(data = W,thre = 0.1,ite1=0,ite2 = 15,ite3 = 0,rep = 1,q=0.9,cor=T)
  Theta.hat1[,,i]=result1$w
  t2=proc.time()
  t=t2-t1
  
  if (i==1)
    print(t)
}

#take average
Theta1=matrix(0,p,p)
for (i in 1:1){
  Theta1=Theta1+Theta.hat1[,,i]
}
Theta1=Theta1/1

round(Theta1[1:12,1:12],3)
Strue[1:12,1:12]
Theta1[which(Theta1<0.1)]=0

spe_sen_bias_kl(Theta =Theta1  ,Strue = Strue,p=1000)

#Theta_guest_p500_s085_free=Theta1;Theta_guest_p500_s09_free=Theta1
#Theta_guest_p1000_s085_free=Theta1;Theta_guest_p1000_s09_free=Theta1

theta.hat=shat[1:12,1:12]
net = theta.hat
net = network::network(net, directed = FALSE)
network::network.vertex.names(net)=paste0("X",network::network.vertex.names(net))
graph = GGally::ggnet2(net,size=3,node.color = "lightgray",label=T,label.size = 3,mode="circle")
graph


shat=diag(0,p,p)
Time=NULL
for (i in 1:1){
  set.seed(1)
  cat("\r",round(i/100*100,2), '%     ')
  #Data generation
  G = XMRF.Sim(n , q1, model = "ISM", graph.type = "scale-free")
  Strue = diag(1,p,p)
  Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)
  X=t(G$X)
  for (j in 1:q2){
    X1=rbinom(n,1,0.9) 
    X=cbind(X,X1)
  }
  
  ##independent case
  #X=NULL
  #for (j in 1:p){
  #  X1=rbinom(n,1,0.9) 
  #  X=cbind(X,X1)
  #}
  
  W=data.frame()
  for (j in 1:p){
    S=rbinom(n,1,0.9) #s=0.85/0.9
    W[1:n,j]=S*X[,j]+(1-S)*(1-X[,j])
  }
  W=matrix(unlist(W),ncol = p)
  
  ##independent
  #X=mvrnorm(n,rep(0,p),diag(1,p,p))
  
  #W=X+mvrnorm(n,rep(0,p),diag(c,p,p))
  
  t1=proc.time()
  
  #shat=shat+glasso(cov(W),rho=0.05)$wi
  
  #clime1=clime(cov(W),lambda = 0.05,sigma = TRUE)
  #shat=shat+matrix(unlist(clime1$Omegalist),nrow = p)
  
  shat=shat+glasso(cov(W)-diag(c,p,p),rho=0.05)$wi
  
  #shat=shat+huge(W, nlambda = 11)$beta[[11]]
  
  #shat=shat+space.joint(W, lam1 = 0.1)$ParCor 
  
  #shat=shat+QUIC(var(W), rho = 0.1)$X
  
  t2=proc.time()
  #if (i==1)
  #  print(t2-t1)
  
  t=t2-t1
  Time=c(Time, t[3])
  if (i==1)
    print(mean(Time))
  
}
shat=shat/5
shat=matrix(as.numeric(shat),nrow=p) #huge

round(Strue[1:12,1:12],3)
round(shat[1:12,1:12],3)
shat[which(abs(shat)<0.1)]=0

spe_sen_bias_kl(Theta =shat  ,Strue = Strue,p=1000)

##Hub
#shat_wainwright_p500_s085_hub=shat;shat_wainwright_s09_hub=shat
#shat_huge_p500_s085_hub=shat;shat_huge_s09_hub=shat
#shat_space_p500_s085_hub=shat;shat_space_s09_hub=shat
#shat_QUIC_p500_s085_hub=shat;shat_QUIC_s09_hub=shat

#shat_wainwright_p1000_s085_hub=shat;shat_wainwright_p1000_s09_hub=shat
#shat_huge_p1000_s085_hub=shat;shat_huge_p1000_s09_hub=shat
#shat_space_p1000_s085_hub=shat;shat_space_p1000_s09_hub=shat
#shat_QUIC_p1000_s085_hub=shat;shat_QUIC_p1000_s09_hub=shat

##Lattice
#shat_wainwright_p500_s085_lat=shat;shat_wainwright_p500_s09_lat=shat
#shat_huge_p500_s085_lat=shat;shat_huge_p500_s09_lat=shat
#shat_space_p500_s085_lat=shat;shat_space_p500_s09_lat=shat
#shat_QUIC_p500_s085_lat=shat;shat_QUIC_p500_s09_lat=shat

#shat_wainwright_p1000_s085_lat=shat;shat_wainwright_p1000_s09_lat=shat
#shat_huge_p1000_s085_lat=shat;shat_huge_p1000_s09_lat=shat
#shat_space_p1000_s085_lat=shat;shat_space_p1000_s09_lat=shat
#shat_QUIC_p1000_s085_lat=shat;shat_QUIC_p1000_s09_lat=shat

##Independent
#shat_wainwright_p500_s085_ind=shat;shat_wainwright_p500_s09_ind=shat
#shat_huge_p500_s085_ind=shat;shat_huge_p500_s09_ind=shat
#shat_space_p500_s085_ind=shat;shat_space_p500_s09_ind=shat
#shat_QUIC_p500_s085_ind=shat;shat_QUIC_p500_s09_ind=shat

#shat_wainwright_p1000_s085_ind=shat;shat_wainwright_p1000_s09_ind=shat
#shat_huge_p1000_s085_ind=shat;shat_huge_p1000_s09_ind=shat
#shat_space_p1000_s085_ind=shat;shat_space_p1000_s09_ind=shat
#shat_QUIC_p1000_s085_ind=shat;shat_QUIC_p1000_s09_ind=shat

##scale-free
#shat_glasso_p500_s085_free=shat;shat_glasso_p500_s09_free=shat
#shat_clime_p500_s085_free=shat;shat_clime_p500_s09_free=shat
#shat_wainwright_p500_s085_free=shat
#shat_huge_p500_s085_free=shat;shat_huge_p500_s09_free=shat
#shat_space_p500_s085_free=shat;shat_space_p500_s09_free=shat
#shat_QUIC_p500_s085_free=shat;shat_QUIC_p500_s09_free=shat

#shat_glasso_p1000_s085_free=shat;shat_glasso_p1000_s09_free=shat
#shat_huge_p1000_s085_free=shat;shat_huge_p1000_s09_free=shat
#shat_space_p1000_s085_free=shat;shat_space_p1000_s09_free=shat
#shat_QUIC_p1000_s085_free=shat;shat_QUIC_p1000_s09_free=shat


setwd("C:/Users/user/Desktop/R")
save.image(file = "revision_1010.RData")

theta.hat=shat[1:12,1:12]
net = theta.hat
net = network::network(net, directed = FALSE)
network::network.vertex.names(net)=paste0("X",network::network.vertex.names(net))
graph = GGally::ggnet2(net,size=3,node.color = "lightgray",label=T,label.size = 3,mode = "circle")
graph

save.image(file = "ISM.RData")

#Package paper
#Theta_p500_085=Theta1 (thre=0.1, ite2=30)
#Theta_p500_09=Theta1 
#Theta_p500_085_lat=Theta1 (thre=0.1, ite2=30)

#Theta_p500_09_lat=Theta1 

#Theta_p1000_085=Theta1 (thre=0.1, ite2=15)
#Theta_p1000_09=Theta1
#Theta_p1000_085_lat=Theta (thre=0.1, ite1=20)
#Theta_p1000_09_lat=Theta (thre=0.1, ite1=15)
#Theta_p500_085_ind=Theta1 (thre=0.1, ite2=30)
#Theta_p500_09_ind=Theta1
#Theta_p1000_085_ind=Theta1 (thre=0.1, ite2=20)
#Theta_p1000_09_ind=Theta1 (thre=0.1, ite2=15)

#Thesis
#Theta_n500_p2000=Theta1 (thre=0.2, ite2=10)
#Theta_n500_p2000_lat=Theta1 (thre=0.1, ite2=10)

