
library(MASS)

n=500
q1=12
q2=988
p=q1+q2 #500/1000

##Data generation
set.seed(1)
G = XMRF.Sim(n , q1, model = "LPGM", graph.type = "lattice")
X=G$X

#(lambda,pi)=(0.5,0.5)/(0.8,0.5)
Z=matrix(rpois(X,lambda = 0.8),nrow = q1)
W=matrix(rbinom(X,1,0.5),nrow=q1)
X_star=X+Z-W
X_star=t(X_star)

X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
X_star=cbind(X_star,X1)

#dependent case
#generate real Y
#rho = exp(colSums(X)) / (1+exp(colSums(X)))
#Y = rbinom(n,1,rho) #GGM
Y = rbinom(n,1,0.5) #ISM/count/not sparse

Z1 = X[which(Y==1),]
Z0 = X[which(Y==0),]

#discriminant
d1_GUEST = NULL
d0_GUEST = NULL
d1_glasso = NULL
d0_glasso = NULL
d1_clime = NULL
d0_clime = NULL

for(i in 1:n) {
  
z = X[i,]

d1_GUEST = c(d1_GUEST, mean(Y) - 0.5 * t(colMeans(Z1)) %*% Theta_p500_08 %*% colMeans(Z1) + t(z) %*% Theta_p500 %*% colMeans(Z1))
d0_GUEST = c(d0_GUEST, mean(Y) - 0.5 * t(colMeans(Z0)) %*% Theta_p500_08 %*% colMeans(Z0) + t(z) %*% Theta_p500 %*% colMeans(Z0))
d1_glasso = c(d1_glasso, mean(Y) - 0.5 * t(colMeans(Z1)) %*% shat%*% colMeans(Z1) + t(z) %*% shat %*% colMeans(Z1))
d0_glasso = c(d0_glasso, mean(Y) - 0.5 * t(colMeans(Z0)) %*% shat%*% colMeans(Z0) + t(z) %*% shat %*% colMeans(Z0))
d1_clime = c(d1_clime, mean(Y) - 0.5 * t(colMeans(Z1)) %*% shat_clime_p500_08 %*% colMeans(Z1) + t(z) %*% shat %*% colMeans(Z1))
d0_clime = c(d0_clime, mean(Y) - 0.5 * t(colMeans(Z0)) %*% shat_clime_p500_08 %*% colMeans(Z0) + t(z) %*% shat %*% colMeans(Z0))

}

r=cbind(Y,1*(d1_GUEST > d0_GUEST),1*(d1_glasso > d0_glasso),1*(d1_clime > d0_clime))

result=GUEST::LDA.boost(X_star, Y, Theta_guest_p1000_lambda08_free, lambda = 0.8) 
#calculation of PRE, REC, F
TP=0;FP=0;FN=0
for (i in 1:length(Y)){
  TP=TP+sum(result$class[i]==0 & Y[i]==0)
  FN=FN+sum(result$class[i]==1 & Y[i]==0)
  FP=FP+sum(result$class[i]==0 & Y[i]==1)
}

PRE=TP/(TP+FP)
REC=TP/(TP+FN)
F=2*((PRE*REC)/(PRE+REC))

PRE;REC;F

##GGM
#Hub
#result_wainwright_p500_c01_hub=c(PRE,REC,F);result_wainwright_p500_c03_hub=c(PRE,REC,F)
#result_huge_p500_c01_hub=c(PRE,REC,F);result_huge_p500_c03_hub=c(PRE,REC,F)
#result_space_p500_c01_hub=c(PRE,REC,F);result_space_p500_c03_hub=c(PRE,REC,F)
#result_QUIC_p500_c01_hub=c(PRE,REC,F);result_QUIC_p500_c03_hub=c(PRE,REC,F)

#result_huge_p1000_c01_hub=c(PRE,REC,F);result_huge_p1000_c03_hub=c(PRE,REC,F)
#result_space_p1000_c01_hub=c(PRE,REC,F);result_space_p1000_c03_hub=c(PRE,REC,F)
#result_QUIC_p1000_c01_hub=c(PRE,REC,F);result_QUIC_p1000_c03_hub=c(PRE,REC,F)

##Lattice
#result_wainwright_p500_c01_lat=c(PRE,REC,F);result_wainwright_p500_c03_lat=c(PRE,REC,F)
#result_huge_p500_c01_lat=c(PRE,REC,F);result_huge_p500_c03_lat=c(PRE,REC,F)
#result_space_p500_c01_lat=c(PRE,REC,F);result_space_p500_c03_lat=c(PRE,REC,F)
#result_QUIC_p500_c01_lat=c(PRE,REC,F);result_QUIC_p500_c03_lat=c(PRE,REC,F)

#result_huge_p1000_c01_lat=c(PRE,REC,F);result_huge_p1000_c03_lat=c(PRE,REC,F)
#result_space_p1000_c01_lat=c(PRE,REC,F);result_space_p1000_c03_lat=c(PRE,REC,F)
#result_QUIC_p1000_c01_lat=c(PRE,REC,F);result_QUIC_p1000_c03_lat=c(PRE,REC,F)

##Independent
#result_wainwright_p500_c01_ind=c(PRE,REC,F);result_wainwright_p500_c03_ind=c(PRE,REC,F)
#result_huge_p500_c01_ind=c(PRE,REC,F);result_huge_p500_c03_ind=c(PRE,REC,F)
#result_space_p500_c01_ind=c(PRE,REC,F);result_space_p500_c03_ind=c(PRE,REC,F)
#result_QUIC_p500_c01_ind=c(PRE,REC,F);result_QUIC_p500_c03_ind=c(PRE,REC,F)

#result_huge_p1000_c01_ind=c(PRE,REC,F);result_huge_p1000_c03_ind=c(PRE,REC,F)
#result_space_p1000_c01_ind=c(PRE,REC,F);result_space_p1000_c03_ind=c(PRE,REC,F)
#result_QUIC_p1000_c01_ind=c(PRE,REC,F);result_QUIC_p1000_c03_ind=c(PRE,REC,F)

##Scale-free
#result_guest_p500_c01_free=c(PRE,REC,F);result_guest_p500_c03_free=c(PRE,REC,F)
#result_glasso_p500_c01_free=c(PRE,REC,F);result_glasso_p500_c03_free=c(PRE,REC,F)
#result_clime_p500_c01_free=c(PRE,REC,F);result_clime_p500_c03_free=c(PRE,REC,F)
#result_wainwright_p500_c01_free=c(PRE,REC,F);result_wainwright_p500_c03_free=c(PRE,REC,F)
#result_huge_p500_c01_free=c(PRE,REC,F);result_huge_p500_c03_free=c(PRE,REC,F)
#result_space_p500_c01_free=c(PRE,REC,F);result_space_p500_c03_free=c(PRE,REC,F)
#result_QUIC_p500_c01_free=c(PRE,REC,F);result_QUIC_p500_c03_free=c(PRE,REC,F)

#result_guest_p1000_c01_free=c(PRE,REC,F);result_guest_p1000_c03_free=c(PRE,REC,F)
#result_glasso_p1000_c01_free=c(PRE,REC,F);result_guest_p1000_c03_free=c(PRE,REC,F)
#result_wainwright_p1000_c01_free=c(PRE,REC,F);result_wainwright_p1000_c03_free=c(PRE,REC,F)
#result_huge_p1000_c01_free=c(PRE,REC,F);result_huge_p1000_c03_free=c(PRE,REC,F)
#result_space_p1000_c01_free=c(PRE,REC,F);result_space_p1000_c03_free=c(PRE,REC,F)
#result_QUIC_p1000_c01_free=c(PRE,REC,F);result_QUIC_p1000_c03_free=c(PRE,REC,F)

##ISM
#Hub
#result_wainwright_p500_s085_hub=c(PRE,REC,F);result_wainwright_p500_s09_hub=c(PRE,REC,F)
#result_huge_p500_s085_hub=c(PRE,REC,F);result_huge_p500_s09_hub=c(PRE,REC,F)
#result_space_p500_s085_hub=c(PRE,REC,F);result_space_p500_s09_hub=c(PRE,REC,F)
#result_QUIC_p500_s085_hub=c(PRE,REC,F);result_QUIC_p500_s09_hub=c(PRE,REC,F)

#result_wainwright_p1000_s085_hub=c(PRE,REC,F);result_wainwright_p1000_s09_hub=c(PRE,REC,F)
#result_huge_p1000_s085_hub=c(PRE,REC,F);result_huge_p1000_s09_hub=c(PRE,REC,F)
#result_space_p1000_s085_hub=c(PRE,REC,F);result_space_p1000_s09_hub=c(PRE,REC,F)
#result_QUIC_p1000_s085_hub=c(PRE,REC,F);result_QUIC_p1000_s09_hub=c(PRE,REC,F)

#Lattice
#result_wainwright_p500_s085_lat=c(PRE,REC,F);result_wainwright_p500_s09_lat=c(PRE,REC,F)
#result_huge_p500_s085_lat=c(PRE,REC,F);result_huge_p500_s09_lat=c(PRE,REC,F)
#result_space_p500_s085_lat=c(PRE,REC,F);result_space_p500_s09_lat=c(PRE,REC,F)
#result_QUIC_p500_s085_lat=c(PRE,REC,F);result_QUIC_p500_s09_lat=c(PRE,REC,F)

#result_wainwright_p1000_s085_lat=c(PRE,REC,F);result_wainwright_p1000_s09_lat=c(PRE,REC,F)
#result_huge_p1000_s085_lat=c(PRE,REC,F);result_huge_p1000_s09_lat=c(PRE,REC,F)
#result_space_p1000_s085_lat=c(PRE,REC,F);result_space_p1000_s09_lat=c(PRE,REC,F)
#result_QUIC_p1000_s085_lat=c(PRE,REC,F);result_QUIC_p1000_s09_lat=c(PRE,REC,F)

#Independent
#result_wainwright_p500_s085_ind=c(PRE,REC,F);result_wainwright_p500_s09_ind=c(PRE,REC,F)
#result_huge_p500_s085_ind=c(PRE,REC,F);result_huge_p500_s09_ind=c(PRE,REC,F)
#result_space_p500_s085_ind=c(PRE,REC,F);result_space_p500_s09_ind=c(PRE,REC,F)
#result_QUIC_p500_s085_ind=c(PRE,REC,F);result_QUIC_p500_s09_ind=c(PRE,REC,F)

#result_wainwright_p1000_s085_ind=c(PRE,REC,F);result_wainwright_p1000_s09_ind=c(PRE,REC,F)
#result_huge_p1000_s085_ind=c(PRE,REC,F);result_huge_p1000_s09_ind=c(PRE,REC,F)
#result_space_p1000_s085_ind=c(PRE,REC,F);result_space_p1000_s09_ind=c(PRE,REC,F)
#result_QUIC_p1000_s085_ind=c(PRE,REC,F);result_QUIC_p1000_s09_ind=c(PRE,REC,F)

##Scale-free
#result_guest_p500_s085_free=c(PRE,REC,F);result_guest_p500_s09_free=c(PRE,REC,F)
#result_glasso_p500_s085_free=c(PRE,REC,F);result_glasso_p500_s09_free=c(PRE,REC,F)
#result_clime_p500_s085_free=c(PRE,REC,F);result_clime_p500_s09_free=c(PRE,REC,F)
#result_wainwright_p500_s085_free=c(PRE,REC,F)
#result_huge_p500_s085_free=c(PRE,REC,F);result_huge_p500_s09_free=c(PRE,REC,F)
#result_space_p500_s085_free=c(PRE,REC,F);result_space_p500_s09_free=c(PRE,REC,F)
#result_QUIC_p500_s085_free=c(PRE,REC,F);result_QUIC_p500_s09_free=c(PRE,REC,F)

#result_guest_p1000_s085_free=c(PRE,REC,F);result_guest_p1000_s09_free=c(PRE,REC,F)
#result_glasso_p1000_s085_free=c(PRE,REC,F);result_glasso_p1000_s09_free=c(PRE,REC,F)
#result_huge_p1000_s085_free=c(PRE,REC,F);result_huge_p1000_s09_free=c(PRE,REC,F)
#result_space_p1000_s085_free=c(PRE,REC,F);result_space_p1000_s09_free=c(PRE,REC,F)
#result_QUIC_p1000_s085_free=c(PRE,REC,F);result_space_p1000_s09_free=c(PRE,REC,F)

##Counts
#Hub
#result_wainwright_p500_lambda05_hub=c(PRE,REC,F);result_wainwright_p500_lambda08_hub=c(PRE,REC,F)
#result_huge_p500_lambda05_hub=c(PRE,REC,F);result_huge_p500_lambda08_hub=c(PRE,REC,F)
#result_space_p500_lambda05_hub=c(PRE,REC,F);result_space_p500_lambda08_hub=c(PRE,REC,F)
#result_QUIC_p500_lambda05_hub=c(PRE,REC,F);result_QUIC_p500_lambda08_hub=c(PRE,REC,F)

#result_huge_p1000_lambda05_hub=c(PRE,REC,F);result_huge_p1000_lambda08_hub=c(PRE,REC,F)
#result_space_p1000_lambda05_hub=c(PRE,REC,F);result_space_p1000_lambda08_hub=c(PRE,REC,F)
#result_QUIC_p1000_lambda05_hub=c(PRE,REC,F);result_QUIC_p1000_lambda08_hub=c(PRE,REC,F)

#Lattice
#result_wainwright_p500_lambda05_lat=c(PRE,REC,F);result_wainwright_p500_lambda08_lat=c(PRE,REC,F)
#result_huge_p500_lambda05_lat=c(PRE,REC,F);result_huge_p500_lambda08_lat=c(PRE,REC,F)
#result_space_p500_lambda05_lat=c(PRE,REC,F);result_space_p500_lambda08_lat=c(PRE,REC,F)
#result_QUIC_p500_lambda05_lat=c(PRE,REC,F);result_QUIC_p500_lambda08_lat=c(PRE,REC,F)

#result_huge_p1000_lambda05_lat=c(PRE,REC,F);result_huge_p1000_lambda08_lat=c(PRE,REC,F)
#result_space_p1000_lambda05_lat=c(PRE,REC,F);result_space_p1000_lambda08_lat=c(PRE,REC,F)
#result_QUIC_p1000_lambda05_lat=c(PRE,REC,F);result_QUIC_p1000_lambda08_lat=c(PRE,REC,F)

##Independent
#result_wainwright_p500_lambda05_ind=c(PRE,REC,F);result_wainwright_p500_lambda08_ind=c(PRE,REC,F)
#result_huge_p500_lambda05_ind=c(PRE,REC,F);result_huge_p500_lambda08_ind=c(PRE,REC,F)
#result_space_p500_lambda05_ind=c(PRE,REC,F);result_space_p500_lambda08_ind=c(PRE,REC,F)
#result_QUIC_p500_lambda05_ind=c(PRE,REC,F);result_QUIC_p500_lambda08_ind=c(PRE,REC,F)

#result_huge_p1000_lambda05_ind=c(PRE,REC,F);result_huge_p1000_lambda08_ind=c(PRE,REC,F)
#result_space_p1000_lambda05_ind=c(PRE,REC,F);result_space_p1000_lambda08_ind=c(PRE,REC,F)
#result_QUIC_p1000_lambda05_ind=c(PRE,REC,F);result_QUIC_p1000_lambda08_ind=c(PRE,REC,F)

##Scale-free
#result_guest_p500_lambda05_free=c(PRE,REC,F);result_guest_p500_lambda08_free=c(PRE,REC,F)
#result_glasso_p500_lambda05_free=c(PRE,REC,F);result_glasso_p500_lambda08_free=c(PRE,REC,F)
#result_clime_p500_lambda05_free=c(PRE,REC,F);result_clime_p500_lambda08_free=c(PRE,REC,F)
#result_huge_p500_lambda05_free=c(PRE,REC,F);result_huge_p500_lambda08_free=c(PRE,REC,F)
#result_space_p500_lambda05_free=c(PRE,REC,F);result_space_p500_lambda08_free=c(PRE,REC,F)
#result_QUIC_p500_lambda05_free=c(PRE,REC,F);result_QUIC_p500_lambda08_free=c(PRE,REC,F)

#result_guest_p1000_lambda05_free=c(PRE,REC,F);result_guest_p1000_lambda08_free=c(PRE,REC,F)
#result_glasso_p1000_lambda05_free=c(PRE,REC,F);result_glasso_p1000_lambda05_free=c(PRE,REC,F)
#result_huge_p1000_lambda05_free=c(PRE,REC,F);result_huge_p1000_lambda08_free=c(PRE,REC,F)
#result_space_p1000_lambda05_free=c(PRE,REC,F);result_space_p1000_lambda08_free=c(PRE,REC,F)
#result_QUIC_p1000_lambda05_free=c(PRE,REC,F);result_QUIC_p1000_lambda08_free=c(PRE,REC,F)

##Not Sparse
#result_glasso_p500_c01_notsparse=c(PRE,REC,F);result_glasso_p500_c03_notsparse=c(PRE,REC,F)
#result_clime_p500_c01_notsparse=c(PRE,REC,F);result_clime_p500_c03_notsparse=c(PRE,REC,F)
#result_guest_p500_c01_notsparse=c(PRE,REC,F);result_guest_p500_c03_notsparse=c(PRE,REC,F)
#result_wainwright_p500_c01_notsparse=c(PRE,REC,F);result_wainwright_p500_c03_notsparse=c(PRE,REC,F)
#result_huge_p500_c01_notsparse=c(PRE,REC,F);result_huge_p500_c03_notsparse=c(PRE,REC,F)
#result_space_p500_c01_notsparse=c(PRE,REC,F);result_space_p500_c03_notsparse=c(PRE,REC,F)
#result_QUIC_p500_c01_notsparse=c(PRE,REC,F);result_QUIC_p500_c03_notsparse=c(PRE,REC,F)
#result_H_p500_c01_notsparse=c(PRE,REC,F);result_H_p500_c03_notsparse=c(PRE,REC,F)

#result_glasso_p1000_c01_notsparse=c(PRE,REC,F);result_glasso_p1000_c03_notsparse=c(PRE,REC,F)
#result_guest_p1000_c01_notsparse=c(PRE,REC,F);result_guest_p1000_c03_notsparse=c(PRE,REC,F)
#result_wainwright_p1000_c01_notsparse=c(PRE,REC,F);result_wainwright_p1000_c03_notsparse=c(PRE,REC,F)
#result_huge_p1000_c01_notsparse=c(PRE,REC,F);result_huge_p1000_c03_notsparse=c(PRE,REC,F)
#result_space_p1000_c01_notsparse=c(PRE,REC,F);result_space_p1000_c03_notsparse=c(PRE,REC,F)
#result_QUIC_p1000_c01_notsparse=c(PRE,REC,F);result_QUIC_p1000_c03_notsparse=c(PRE,REC,F)
#result_H_p1000_c01_notsparse=c(PRE,REC,F);result_H_p1000_c03_notsparse=c(PRE,REC,F)

setwd("C:/Users/user/Desktop/R")
save.image(file = "Classification_free_1012(ISM+Counts).RData")


n=500
q1 = 12
q2 = 988
p = q1+q2 #500/1000

G = XMRF.Sim(n , q1, model = "LPGM", graph.type = "lattice")
X=G$X

#(lambda,pi)=(0.5,0.5)/(0.8,0.5)
Z=matrix(rpois(X,lambda = 0.8),nrow = q1)
W=matrix(rbinom(X,1,0.5),nrow=q1)
X_star=X+Z-W
X_star=t(X_star)

X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
X_star=cbind(X_star,X1)

##independent case
X=diag(1,q1,n)
Y=matrix(rpois(X,lambda = 0.5),nrow = q1)
#(lambda,pi)=(0.5,0.5)/(0.8,0.5)
Z=matrix(rpois(Y,lambda = 0.8),nrow = q1)
W=matrix(rbinom(Y,1,0.5),nrow=q1)
X_star=Y+Z-W
X_star=t(X_star)

X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
X_star=cbind(X_star,X1)


#Data generation
G = XMRF.Sim(n , q1, model = "GGM", graph.type = "hub")
X=t(G$X)
X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
X=cbind(X,X1)

c=0.3
W=X+mvrnorm(n,rep(0,p),diag(c,p,p))

#Independent Data generation
X=mvrnorm(n,rep(0,p),diag(1,p,p))

c=0.3
W=X+mvrnorm(n,rep(0,p),diag(c,p,p))

##Data generation
#dependent case
G = XMRF.Sim(n , q1, model = "ISM", graph.type = "lattice")
X=t(G$X)
for (j in 1:q2){
  X1=rbinom(n,1,0.9) 
  X=cbind(X,X1)
}

#independent case
X=NULL
for (j in 1:p){
  X1=rbinom(n,1,0.9) 
  X=cbind(X,X1)
}

W=data.frame()
for (j in 1:p){
  S=rbinom(n,1,0.9) #s=0.85/0.9
  W[1:n,j]=S*X[,j]+(1-S)*(1-X[,j])
}
W=matrix(unlist(W),ncol = p)


## Not sparse
S=matrix(0.8,p,p)+diag(0.2,p,p)
mu=rep(0,p)

X=mvrnorm(n, mu, S)

c=0.3 #0.1/0.3
W=X+mvrnorm(n,rep(0,p),diag(c,p,p))

