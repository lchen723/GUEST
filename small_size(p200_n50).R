
n=50
q1 = 12
q2 = 188
p = q1+q2 #200

##GGM
Theta.hat1=array(dim = c(p,p,10))
Time=NULL
for (i in 1:10){
  #data generation
  cat("\r",round(i/100*100,2), '%     ')
  set.seed(i)
  #Data generation
  G = XMRF.Sim(n , q1, model = "GGM", graph.type = "hub")
  X=t(G$X)
  X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
  X=cbind(X,X1)
  Strue = diag(1,p,p)
  Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)
  c=0.3
  W=X+mvrnorm(n,rep(0,p),diag(c,p,p))
  
  #boost.graph
  t1=proc.time()
  result1<-boost.graph(data = W,thre = 0.1,ite1=7,ite2 = 0,ite3 = 0,rep = 1,sigma_e = c,cor=T)
  Theta.hat1[,,i]=result1$w
  t2=proc.time()
  t=t2-t1
  
  Time=c(Time, t[3])
  if (i==10)
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
Theta1[which(Theta1<0.13)]=0

#Theta_guest_p200_n50_c01=Theta1
#Theta_guest_p200_n50_c03=Theta1

spe_sen_bias_kl(Theta =Theta_guest_p200_n50_c01  ,Strue = Strue,p=200)
spe_sen_bias_kl(Theta =Theta_guest_p200_n50_c03  ,Strue = Strue,p=200)


##ISM
Theta.hat1=array(dim = c(p,p,10))
Time=NULL
for (i in 1:10){
  #data generation
  cat("\r",round(i/100*100,2), '%     ')
  set.seed(i)
  #Data generation
  G = XMRF.Sim(n , q1, model = "ISM", graph.type = "hub")
  X=t(G$X)
  Strue = diag(1,p,p)
  Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)
  for (j in 1:q2){
    X1=rbinom(n,1,0.9) 
    X=cbind(X,X1)
  }
  
  W=data.frame()
  for (j in 1:p){
    S=rbinom(n,1,0.9) #s=0.85/0.9
    W[1:n,j]=S*X[,j]+(1-S)*(1-X[,j])
  }
  W=matrix(unlist(W),ncol = p)
  
  #boost.graph
  t1=proc.time()
  result1<-boost.graph(data = W,thre = 0.1,ite1=0,ite2 = 50,ite3 = 0,rep = 1,q = 0.9,cor=T)
  Theta.hat1[,,i]=result1$w
  t2=proc.time()
  t=t2-t1
  
  Time=c(Time, t[3])
  if (i==10)
    print(mean(Time))
}

#take average
Theta1=matrix(0,p,p)
for (i in 1:1){
  Theta1=Theta1+Theta.hat1[,,i]
}
Theta1=Theta1/1

round(Theta1[1:12,1:12],3)
Strue[1:12,1:12]
Theta1[which(Theta1<0.001)]=0


#Theta_guest_p200_n50_s085=Theta1
#Theta_guest_p200_n50_s09=Theta1

spe_sen_bias_kl(Theta =Theta_guest_p200_n50_s085  ,Strue = Strue,p=200)
spe_sen_bias_kl(Theta =Theta_guest_p200_n50_s09  ,Strue = Strue,p=200)


##Counts
Theta.hat1=array(dim = c(p,p,10))
Time=NULL
for (i in 1:10){
  #data generation
  cat("\r",round(i/100*100,2), '%     ')
  set.seed(i)
  ##Data generation
  #dependent case
  G = XMRF.Sim(n , q1, model = "LPGM", graph.type = "hub")
  X=G$X
  Strue = diag(1,p,p)
  Strue[1:q1,1:q1]= G$B +diag(max((eigen(Strue))$values+0.1),q1)
  
  #(lambda,pi)=(0.5,0.5)/(0.8,0.5)
  Z=matrix(rpois(X,lambda = 0.8),nrow = q1)
  W=matrix(rbinom(X,1,0.5),nrow=q1)
  X_star=X+Z-W
  X_star=t(X_star)
  
  X1=mvrnorm(n,rep(0,q2),diag(1,q2,q2))
  X_star=cbind(X_star,X1)
  
  #boost.graph
  t1=proc.time()
  result1<-boost.graph(data = X_star,thre = 0.3,ite1=0,ite2 = 0,ite3 = 2,rep = 1,lambda = 0.8,cor=T)
  Theta.hat1[,,i]=result1$w
  t2=proc.time()
  t=t2-t1
  
  Time=c(Time, t[3])
  if (i==10)
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
Theta1[which(Theta1<0.155)]=0


#Theta_guest_p200_n50_lambda05=Theta1
#Theta_guest_p200_n50_lambda08=Theta1


spe_sen_bias_kl(Theta =Theta_guest_p200_n50_lambda05  ,Strue = Strue,p=200)
spe_sen_bias_kl(Theta =Theta_guest_p200_n50_lambda08  ,Strue = Strue,p=200)


#save.image(file = "small_dimension1015")
