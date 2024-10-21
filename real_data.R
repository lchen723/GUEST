# load plsgenomics library
library(plsgenomics)

# load data set
data(SRBCT)

# how many samples and how many genes ?
dim(SRBCT$X)

# how many samples of class 1, 2, 3 and 4, respectively ?
sum(SRBCT$Y==1) #Ewing's sarcoma (EWS): 29 examples (34.9%)
sum(SRBCT$Y==2) #Burkitt's lymphoma (BL): 11 examples (13.3%)
sum(SRBCT$Y==3) #neuroblastoma (NB): 18 examples (21.7%)
sum(SRBCT$Y==4) #rhabdomyosarcoma (RMS): 25 examples (30.1%)

X=SRBCT$X
X_new = scale(X)

Y=SRBCT$Y

#common entries
c1=0.15;c2=0.35;c3=0.55

H=var(X_new)%*%(var(X_new)+diag(1,2308,2308))^(-2)

#plot
#glasso
Strue1=glasso(cov(X_new),rho = 0.05)$wi
round(Strue1[1:10,1:10],3)

temp1=Strue1[1:50,1:50]
net = temp1
net = network(net, directed = FALSE)
network.vertex.names(net)=paste0("X",network.vertex.names(net))
graph1 = ggnet2(net,size=3,node.color = "lightgray",label = T,label.size = 3,mode = "circle")
graph1

#boost.graph c1
result4_1=boost.graph(data = X, ite1 = 1, ite2 = 0, thre = 0.09, sigma_e = c1, rep = 1)
Strue4_1=result4_1$w
round(Strue4_1[1:10,1:10],3)

temp4_1=Strue4_1[1:50,1:50]
net = temp4_1
net = network(net, directed = FALSE)
network.vertex.names(net)=paste0("X",network.vertex.names(net))
graph4_1 = ggnet2(net,size=3,node.color = "lightgray",label = T,label.size = 3,mode = "circle")
graph4_1

#boost.graph c2
result4_2=boost.graph(data = X_new, ite1 = 1, ite2 = 0, thre = 0.09, sigma_e = c2, rep = 1)
Strue4_2=result4_2$w
round(Strue4_2[1:10,1:10],3)
temp4_2=Strue4_2[1:50,1:50]
temp4_2[which(temp4_2<0.2)]=0

net = temp4_2
net = network(net, directed = FALSE)
network.vertex.names(net)=paste0("X",network.vertex.names(net))
graph4_2 = ggnet2(net,size=3,node.color = "lightgray",label = T,label.size = 3,mode = "circle")
graph4_2

#boost.graph c3
result4_3=boost.graph(data = X_new, ite1 = 1, ite2 = 0, thre = 0.09, sigma_e = c3, rep = 1)
Strue4_3=result4_3$w
temp4_3=Strue4_3[1:50,1:50]

net = temp4_3
net = network(net, directed = FALSE)
network.vertex.names(net)=paste0("X",network.vertex.names(net))
graph4_3 = ggnet2(net,size=3,node.color = "lightgray",label = T,label.size = 3,mode = "circle")
graph4_3


#boost.graph c4
c4=0.75
result4_4=GUEST::boost.graph(data = X_new, ite1 = 1, ite2 = 0, thre = 0.09, sigma_e = c4, rep = 1)
Strue4_4=result4_4$w
temp4_4=Strue4_4[1:50,1:50]

net = temp4_4
net = network(net, directed = FALSE)
network.vertex.names(net)=paste0("X",network.vertex.names(net))
graph4_4 = ggnet2(net,size=3,node.color = "lightgray",label = T,label.size = 3,mode = "circle")
graph4_4

#boost.graph c5
c5=0.95
result4_5=GUEST::boost.graph(data = X_new, ite1 = 1, ite2 = 0, thre = 0.09, sigma_e = c5, rep = 1)
Strue4_5=result4_5$w
temp4_5=Strue4_5[1:50,1:50]

net = temp4_5
net = network(net, directed = FALSE)
network.vertex.names(net)=paste0("X",network.vertex.names(net))
graph4_5 = ggnet2(net,size=3,node.color = "lightgray",label = T,label.size = 3,mode = "circle")
graph4_5


#huge
Strue5=huge(X_new, nlambda = 11)$beta[[11]]
Strue5=matrix(as.numeric(Strue5),nrow=2308)
round(Strue5[1:10,1:10],3)

temp5=Strue5[1:50,1:50]
temp5[which(temp5<0.07)]=0
net = temp5
net = network(net, directed = FALSE)
network.vertex.names(net)=paste0("X",network.vertex.names(net))
graph5 = ggnet2(net,size=3,node.color = "lightgray",label = T,label.size = 3,mode = "circle")
graph5

#space
Strue6=space.joint(X_new, lam1 = 0.1)$ParCor 
round(Strue6[1:10,1:10],3)

temp6=Strue6[1:50,1:50]
temp6[which(temp6<0.15)]=0
net = temp6
net = network(net, directed = FALSE)
network.vertex.names(net)=paste0("X",network.vertex.names(net))
graph6 = ggnet2(net,size=3,node.color = "lightgray",label = T,label.size = 3,mode = "circle")
graph6

#QUIC
Strue7=QUIC(var(X_new), rho = 0.1)$X
round(Strue7[1:10,1:10],3)

temp7=Strue7[1:50,1:50]
net = temp7
net = network(net, directed = FALSE)
network.vertex.names(net)=paste0("X",network.vertex.names(net))
graph7 = ggnet2(net,size=3,node.color = "lightgray",label = T,label.size = 3,mode = "circle")
graph7

#H
temp8=H[1:50,1:50]
temp8[which(temp8<10^10)]=0
temp8=temp8+diag(1,50,50)
net = temp8
net = network(net, directed = FALSE)
network.vertex.names(net)=paste0("X",network.vertex.names(net))
graph8 = ggnet2(net,size=3,node.color = "lightgray",label = T,label.size = 3,mode = "circle")
graph8

#predict
temp1=Strue1
result_g=LDA.boost(X_new, Y, temp1)

TP1=0;FN1=0;FP1=0
TP2=0;FN2=0;FP2=0
TP3=0;FN3=0;FP3=0
TP4=0;FN4=0;FP4=0
for (i in 1:length(Y)){
  TP1=TP1+sum(result_g$class[i]==1 & Y[i]==1)
  FN1=FN1+sum(result_g$class[i]!=1 & Y[i]==1)
  FP1=FP1+sum(result_g$class[i]==1 & Y[i]!=1)

  TP2=TP2+sum(result_g$class[i]==2 & Y[i]==2)
  FN2=FN2+sum(result_g$class[i]!=2 & Y[i]==2)
  FP2=FP2+sum(result_g$class[i]==2 & Y[i]!=2)

  TP3=TP3+sum(result_g$class[i]==3 & Y[i]==3)
  FN3=FN3+sum(result_g$class[i]!=3 & Y[i]==3)
  FP3=FP3+sum(result_g$class[i]==3 & Y[i]!=3)

  TP4=TP4+sum(result_g$class[i]==4 & Y[i]==4)
  FN4=FN4+sum(result_g$class[i]!=4 & Y[i]==4)
  FP4=FP4+sum(result_g$class[i]==4 & Y[i]!=4)
}
P1_g=TP1/(TP1+FP1)
R1_g=TP1/(TP1+FN1)

P2_g=TP2/(TP2+FP2)
R2_g=TP2/(TP2+FN2)

P3_g=TP3/(TP3+FP3)
R3_g=TP3/(TP3+FN3)

P4_g=TP4/(TP4+FP4)
R4_g=TP4/(TP4+FN4)

A_g=(TP1+TP2+TP3+TP4)/length(Y)

PRE_g=(P1_g+P2_g+P3_g+P4_g)/4
REC_g=(R1_g+R2_g+R3_g+R4_g)/4
F_g=2*((PRE_g*REC_g)/(PRE_g+REC_g))

A_g;P1_g;P2_g;P3_g;P4_g;R1_g;R2_g;R3_g;R4_g;PRE_g;REC_g;F_g

temp4_1=Strue4_1
temp4_1[which(Strue4_1<1.7)]=0
temp4_1=temp4_1+diag(1,2308,2308)
result_b1=LDA.boost(X_new, Y, temp4_1)

TP1=0;FN1=0;FP1=0
TP2=0;FN2=0;FP2=0
TP3=0;FN3=0;FP3=0
TP4=0;FN4=0;FP4=0
for (i in 1:length(Y)){
  TP1=TP1+sum(result_b1$class[i]==1 & Y[i]==1)
  FN1=FN1+sum(result_b1$class[i]!=1 & Y[i]==1)
  FP1=FP1+sum(result_b1$class[i]==1 & Y[i]!=1)

  TP2=TP2+sum(result_b1$class[i]==2 & Y[i]==2)
  FN2=FN2+sum(result_b1$class[i]!=2 & Y[i]==2)
  FP2=FP2+sum(result_b1$class[i]==2 & Y[i]!=2)

  TP3=TP3+sum(result_b1$class[i]==3 & Y[i]==3)
  FN3=FN3+sum(result_b1$class[i]!=3 & Y[i]==3)
  FP3=FP3+sum(result_b1$class[i]==3 & Y[i]!=3)

  TP4=TP4+sum(result_b1$class[i]==4 & Y[i]==4)
  FN4=FN4+sum(result_b1$class[i]!=4 & Y[i]==4)
  FP4=FP4+sum(result_b1$class[i]==4 & Y[i]!=4)
}
P1_b1=TP1/(TP1+FP1)
R1_b1=TP1/(TP1+FN1)

P2_b1=TP2/(TP2+FP2)
R2_b1=TP2/(TP2+FN2)

P3_b1=TP3/(TP3+FP3)
R3_b1=TP3/(TP3+FN3)

P4_b1=TP4/(TP4+FP4)
R4_b1=TP4/(TP4+FN4)

A_b1=(TP1+TP2+TP3+TP4)/length(Y)

PRE_b1=(P1_b1+P2_b1+P3_b1+P4_b1)/4
REC_b1=(R1_b1+R2_b1+R3_b1+R4_b1)/4
F_b1=2*((PRE_b1*REC_b1)/(PRE_b1+REC_b1))

A_b1;P1_b1;P2_b1;P3_b1;P4_b1;R1_b1;R2_b1;R3_b1;R4_b1;PRE_b1;REC_b1;F_b1

temp4_2=Strue4_2
temp4_2[which(Strue4_2<1.7)]=0
temp4_2=temp4_2+diag(1,2308,2308)
result_b2=LDA.boost(X_new, Y, temp4_2)

TP1=0;FN1=0;FP1=0
TP2=0;FN2=0;FP2=0
TP3=0;FN3=0;FP3=0
TP4=0;FN4=0;FP4=0
for (i in 1:length(Y)){
  TP1=TP1+sum(result_b2$class[i]==1 & Y[i]==1)
  FN1=FN1+sum(result_b2$class[i]!=1 & Y[i]==1)
  FP1=FP1+sum(result_b2$class[i]==1 & Y[i]!=1)

  TP2=TP2+sum(result_b2$class[i]==2 & Y[i]==2)
  FN2=FN2+sum(result_b2$class[i]!=2 & Y[i]==2)
  FP2=FP2+sum(result_b2$class[i]==2 & Y[i]!=2)

  TP3=TP3+sum(result_b2$class[i]==3 & Y[i]==3)
  FN3=FN3+sum(result_b2$class[i]!=3 & Y[i]==3)
  FP3=FP3+sum(result_b2$class[i]==3 & Y[i]!=3)

  TP4=TP4+sum(result_b2$class[i]==4 & Y[i]==4)
  FN4=FN4+sum(result_b2$class[i]!=4 & Y[i]==4)
  FP4=FP4+sum(result_b2$class[i]==4 & Y[i]!=4)
}
P1_b2=TP1/(TP1+FP1)
R1_b2=TP1/(TP1+FN1)

P2_b2=TP2/(TP2+FP2)
R2_b2=TP2/(TP2+FN2)

P3_b2=TP3/(TP3+FP3)
R3_b2=TP3/(TP3+FN3)

P4_b2=TP4/(TP4+FP4)
R4_b2=TP4/(TP4+FN4)

A_b2=(TP1+TP2+TP3+TP4)/length(Y)

PRE_b2=(P1_b2+P2_b2+P3_b2+P4_b2)/4
REC_b2=(R1_b2+R2_b2+R3_b2+R4_b2)/4
F_b2=2*((PRE_b2*REC_b2)/(PRE_b2+REC_b2))

A_b2;P1_b2;P2_b2;P3_b2;P4_b2;R1_b2;R2_b2;R3_b2;R4_b2;PRE_b2;REC_b2;F_b2

temp4_3=Strue4_3
temp4_3[which(Strue4_3<1.7)]=0
temp4_3=temp4_3+diag(1,2308,2308)
result_b3=LDA.boost(X_new, Y, temp4_3)

TP1=0;FN1=0;FP1=0
TP2=0;FN2=0;FP2=0
TP3=0;FN3=0;FP3=0
TP4=0;FN4=0;FP4=0
for (i in 1:length(Y)){
  TP1=TP1+sum(result_b3$class[i]==1 & Y[i]==1)
  FN1=FN1+sum(result_b3$class[i]!=1 & Y[i]==1)
  FP1=FP1+sum(result_b3$class[i]==1 & Y[i]!=1)

  TP2=TP2+sum(result_b3$class[i]==2 & Y[i]==2)
  FN2=FN2+sum(result_b3$class[i]!=2 & Y[i]==2)
  FP2=FP2+sum(result_b3$class[i]==2 & Y[i]!=2)

  TP3=TP3+sum(result_b3$class[i]==3 & Y[i]==3)
  FN3=FN3+sum(result_b3$class[i]!=3 & Y[i]==3)
  FP3=FP3+sum(result_b3$class[i]==3 & Y[i]!=3)

  TP4=TP4+sum(result_b3$class[i]==4 & Y[i]==4)
  FN4=FN4+sum(result_b3$class[i]!=4 & Y[i]==4)
  FP4=FP4+sum(result_b3$class[i]==4 & Y[i]!=4)
}
P1_b3=TP1/(TP1+FP1)
R1_b3=TP1/(TP1+FN1)

P2_b3=TP2/(TP2+FP2)
R2_b3=TP2/(TP2+FN2)

P3_b3=TP3/(TP3+FP3)
R3_b3=TP3/(TP3+FN3)

P4_b3=TP4/(TP4+FP4)
R4_b3=TP4/(TP4+FN4)

A_b3=(TP1+TP2+TP3+TP4)/length(Y)

PRE_b3=(P1_b3+P2_b3+P3_b3+P4_b3)/4
REC_b3=(R1_b3+R2_b3+R3_b3+R4_b3)/4
F_b3=2*((PRE_b3*REC_b3)/(PRE_b3+REC_b3))

A_b3;P1_b3;P2_b3;P3_b3;P4_b3;R1_b3;R2_b3;R3_b3;R4_b3;PRE_b3;REC_b3;F_b3



#predict with c4
temp4_4=Strue4_4
temp4_4[which(Strue4_4<1.85)]=0
temp4_4=temp4_4+diag(1,2308,2308)
result_b4=LDA.boost(X_new, Y, temp4_4)

TP1=0;FN1=0;FP1=0
TP2=0;FN2=0;FP2=0
TP3=0;FN3=0;FP3=0
TP4=0;FN4=0;FP4=0
for (i in 1:length(Y)){
  TP1=TP1+sum(result_b4$class[i]==1 & Y[i]==1)
  FN1=FN1+sum(result_b4$class[i]!=1 & Y[i]==1)
  FP1=FP1+sum(result_b4$class[i]==1 & Y[i]!=1)
  
  TP2=TP2+sum(result_b4$class[i]==2 & Y[i]==2)
  FN2=FN2+sum(result_b4$class[i]!=2 & Y[i]==2)
  FP2=FP2+sum(result_b4$class[i]==2 & Y[i]!=2)
  
  TP3=TP3+sum(result_b4$class[i]==3 & Y[i]==3)
  FN3=FN3+sum(result_b4$class[i]!=3 & Y[i]==3)
  FP3=FP3+sum(result_b4$class[i]==3 & Y[i]!=3)
  
  TP4=TP4+sum(result_b4$class[i]==4 & Y[i]==4)
  FN4=FN4+sum(result_b4$class[i]!=4 & Y[i]==4)
  FP4=FP4+sum(result_b4$class[i]==4 & Y[i]!=4)
}
P1_b4=TP1/(TP1+FP1)
R1_b4=TP1/(TP1+FN1)

P2_b4=TP2/(TP2+FP2)
R2_b4=TP2/(TP2+FN2)

P3_b4=TP3/(TP3+FP3)
R3_b4=TP3/(TP3+FN3)

P4_b4=TP4/(TP4+FP4)
R4_b4=TP4/(TP4+FN4)

A_b4=(TP1+TP2+TP3+TP4)/length(Y)

PRE_b4=(P1_b4+P2_b4+P3_b4+P4_b4)/4
REC_b4=(R1_b4+R2_b4+R3_b4+R4_b4)/4
F_b4=2*((PRE_b4*REC_b4)/(PRE_b4+REC_b4))

A_b4;P1_b4;P2_b4;P3_b4;P4_b4;R1_b4;R2_b4;R3_b4;R4_b4;PRE_b4;REC_b4;F_b4

#predict with c5
temp4_5=Strue4_5
temp4_5[which(Strue4_5<1.85)]=0
temp4_5=temp4_5+diag(1,2308,2308)
result_b5=LDA.boost(X_new, Y, temp4_5)

TP1=0;FN1=0;FP1=0
TP2=0;FN2=0;FP2=0
TP3=0;FN3=0;FP3=0
TP4=0;FN4=0;FP4=0
for (i in 1:length(Y)){
  TP1=TP1+sum(result_b5$class[i]==1 & Y[i]==1)
  FN1=FN1+sum(result_b5$class[i]!=1 & Y[i]==1)
  FP1=FP1+sum(result_b5$class[i]==1 & Y[i]!=1)
  
  TP2=TP2+sum(result_b5$class[i]==2 & Y[i]==2)
  FN2=FN2+sum(result_b5$class[i]!=2 & Y[i]==2)
  FP2=FP2+sum(result_b5$class[i]==2 & Y[i]!=2)
  
  TP3=TP3+sum(result_b5$class[i]==3 & Y[i]==3)
  FN3=FN3+sum(result_b5$class[i]!=3 & Y[i]==3)
  FP3=FP3+sum(result_b5$class[i]==3 & Y[i]!=3)
  
  TP4=TP4+sum(result_b5$class[i]==4 & Y[i]==4)
  FN4=FN4+sum(result_b5$class[i]!=4 & Y[i]==4)
  FP4=FP4+sum(result_b5$class[i]==4 & Y[i]!=4)
}
P1_b5=TP1/(TP1+FP1)
R1_b5=TP1/(TP1+FN1)

P2_b5=TP2/(TP2+FP2)
R2_b5=TP2/(TP2+FN2)

P3_b5=TP3/(TP3+FP3)
R3_b5=TP3/(TP3+FN3)

P4_b5=TP4/(TP4+FP4)
R4_b5=TP4/(TP4+FN4)

A_b5=(TP1+TP2+TP3+TP4)/length(Y)

PRE_b5=(P1_b5+P2_b5+P3_b5+P4_b5)/4
REC_b5=(R1_b5+R2_b5+R3_b5+R4_b5)/4
F_b5=2*((PRE_b5*REC_b5)/(PRE_b5+REC_b5))

A_b5;P1_b5;P2_b5;P3_b5;P4_b5;R1_b5;R2_b5;R3_b5;R4_b5;PRE_b5;REC_b5;F_b5



#huge
temp5=Strue5
temp5[which(temp5<0.1)]=0
result_huge=LDA.boost(X_new, Y, temp5)

TP1=0;FN1=0;FP1=0
TP2=0;FN2=0;FP2=0
TP3=0;FN3=0;FP3=0
TP4=0;FN4=0;FP4=0
for (i in 1:length(Y)){
  TP1=TP1+sum(result_huge$class[i]==1 & Y[i]==1)
  FN1=FN1+sum(result_huge$class[i]!=1 & Y[i]==1)
  FP1=FP1+sum(result_huge$class[i]==1 & Y[i]!=1)
  
  TP2=TP2+sum(result_huge$class[i]==2 & Y[i]==2)
  FN2=FN2+sum(result_huge$class[i]!=2 & Y[i]==2)
  FP2=FP2+sum(result_huge$class[i]==2 & Y[i]!=2)
  
  TP3=TP3+sum(result_huge$class[i]==3 & Y[i]==3)
  FN3=FN3+sum(result_huge$class[i]!=3 & Y[i]==3)
  FP3=FP3+sum(result_huge$class[i]==3 & Y[i]!=3)
  
  TP4=TP4+sum(result_huge$class[i]==4 & Y[i]==4)
  FN4=FN4+sum(result_huge$class[i]!=4 & Y[i]==4)
  FP4=FP4+sum(result_huge$class[i]==4 & Y[i]!=4)
}
P1_huge=TP1/(TP1+FP1)
R1_huge=TP1/(TP1+FN1)

P2_huge=TP2/(TP2+FP2)
R2_huge=TP2/(TP2+FN2)

P3_huge=TP3/(TP3+FP3)
R3_huge=TP3/(TP3+FN3)

P4_huge=TP4/(TP4+FP4)
R4_huge=TP4/(TP4+FN4)

A_huge=(TP1+TP2+TP3+TP4)/length(Y)

PRE_huge=(P1_huge+P2_huge+P3_huge+P4_huge)/4
REC_huge=(R1_huge+R2_huge+R3_huge+R4_huge)/4
F_huge=2*((PRE_huge*REC_huge)/(PRE_huge+REC_huge))

A_huge;P1_huge;P2_huge;P3_huge;P4_huge;R1_huge;R2_huge;R3_huge;R4_huge;PRE_huge;REC_huge;F_huge


temp6=Strue6
temp6[which(temp6<0.15)]=0
result_space=LDA.boost(X_new, Y, temp6)

TP1=0;FN1=0;FP1=0
TP2=0;FN2=0;FP2=0
TP3=0;FN3=0;FP3=0
TP4=0;FN4=0;FP4=0
for (i in 1:length(Y)){
  TP1=TP1+sum(result_space$class[i]==1 & Y[i]==1)
  FN1=FN1+sum(result_space$class[i]!=1 & Y[i]==1)
  FP1=FP1+sum(result_space$class[i]==1 & Y[i]!=1)
  
  TP2=TP2+sum(result_space$class[i]==2 & Y[i]==2)
  FN2=FN2+sum(result_space$class[i]!=2 & Y[i]==2)
  FP2=FP2+sum(result_space$class[i]==2 & Y[i]!=2)
  
  TP3=TP3+sum(result_space$class[i]==3 & Y[i]==3)
  FN3=FN3+sum(result_space$class[i]!=3 & Y[i]==3)
  FP3=FP3+sum(result_space$class[i]==3 & Y[i]!=3)
  
  TP4=TP4+sum(result_space$class[i]==4 & Y[i]==4)
  FN4=FN4+sum(result_space$class[i]!=4 & Y[i]==4)
  FP4=FP4+sum(result_space$class[i]==4 & Y[i]!=4)
}
P1_space=TP1/(TP1+FP1)
R1_space=TP1/(TP1+FN1)

P2_space=TP2/(TP2+FP2)
R2_space=TP2/(TP2+FN2)

P3_space=TP3/(TP3+FP3)
R3_space=TP3/(TP3+FN3)

P4_space=TP4/(TP4+FP4)
R4_space=TP4/(TP4+FN4)

A_space=(TP1+TP2+TP3+TP4)/length(Y)

PRE_space=(P1_space+P2_space+P3_space+P4_space)/4
REC_space=(R1_space+R2_space+R3_space+R4_space)/4
F_space=2*((PRE_space*REC_space)/(PRE_space+REC_space))

A_space;P1_space;P2_space;P3_space;P4_space;R1_space;R2_space;R3_space;R4_space;PRE_space;REC_space;F_space

#QUIC
temp7=Strue7
result_QUIC=LDA.boost(X_new, Y, temp7)

TP1=0;FN1=0;FP1=0
TP2=0;FN2=0;FP2=0
TP3=0;FN3=0;FP3=0
TP4=0;FN4=0;FP4=0
for (i in 1:length(Y)){
  TP1=TP1+sum(result_QUIC$class[i]==1 & Y[i]==1)
  FN1=FN1+sum(result_QUIC$class[i]!=1 & Y[i]==1)
  FP1=FP1+sum(result_QUIC$class[i]==1 & Y[i]!=1)
  
  TP2=TP2+sum(result_QUIC$class[i]==2 & Y[i]==2)
  FN2=FN2+sum(result_QUIC$class[i]!=2 & Y[i]==2)
  FP2=FP2+sum(result_QUIC$class[i]==2 & Y[i]!=2)

  TP3=TP3+sum(result_QUIC$class[i]==3 & Y[i]==3)
  FN3=FN3+sum(result_QUIC$class[i]!=3 & Y[i]==3)
  FP3=FP3+sum(result_QUIC$class[i]==3 & Y[i]!=3)

  TP4=TP4+sum(result_QUIC$class[i]==4 & Y[i]==4)
  FN4=FN4+sum(result_QUIC$class[i]!=4 & Y[i]==4)
  FP4=FP4+sum(result_QUIC$class[i]==4 & Y[i]!=4)
}
P1_QUIC=TP1/(TP1+FP1)
R1_QUIC=TP1/(TP1+FN1)

P2_QUIC=TP2/(TP2+FP2)
R2_QUIC=TP2/(TP2+FN2)

P3_QUIC=TP3/(TP3+FP3)
R3_QUIC=TP3/(TP3+FN3)

P4_QUIC=TP4/(TP4+FP4)
R4_QUIC=TP4/(TP4+FN4)

A_QUIC=(TP1+TP2+TP3+TP4)/length(Y)

PRE_QUIC=(P1_QUIC+P2_QUIC+P3_QUIC+P4_QUIC)/4
REC_QUIC=(R1_QUIC+R2_QUIC+R3_QUIC+R4_QUIC)/4
F_QUIC=2*((PRE_QUIC*REC_QUIC)/(PRE_QUIC+REC_QUIC))

A_QUIC;P1_QUIC;P2_QUIC;P3_QUIC;P4_QUIC;R1_QUIC;R2_QUIC;R3_QUIC;R4_QUIC;PRE_QUIC;REC_QUIC;F_QUIC


#H
temp8=H
result_H=LDA.boost(X_new, Y, temp8)

TP1=0;FN1=0;FP1=0
TP2=0;FN2=0;FP2=0
TP3=0;FN3=0;FP3=0
TP4=0;FN4=0;FP4=0
for (i in 1:length(Y)){
  TP1=TP1+sum(result_H$class[i]==1 & Y[i]==1)
  FN1=FN1+sum(result_H$class[i]!=1 & Y[i]==1)
  FP1=FP1+sum(result_H$class[i]==1 & Y[i]!=1)
  
  TP2=TP2+sum(result_H$class[i]==2 & Y[i]==2)
  FN2=FN2+sum(result_H$class[i]!=2 & Y[i]==2)
  FP2=FP2+sum(result_H$class[i]==2 & Y[i]!=2)
  
  TP3=TP3+sum(result_H$class[i]==3 & Y[i]==3)
  FN3=FN3+sum(result_H$class[i]!=3 & Y[i]==3)
  FP3=FP3+sum(result_H$class[i]==3 & Y[i]!=3)
  
  TP4=TP4+sum(result_H$class[i]==4 & Y[i]==4)
  FN4=FN4+sum(result_H$class[i]!=4 & Y[i]==4)
  FP4=FP4+sum(result_H$class[i]==4 & Y[i]!=4)
}
P1_H=TP1/(TP1+FP1)
R1_H=TP1/(TP1+FN1)

P2_H=TP2/(TP2+FP2)
R2_H=TP2/(TP2+FN2)

P3_H=TP3/(TP3+FP3)
R3_H=TP3/(TP3+FN3)

P4_H=TP4/(TP4+FP4)
R4_H=TP4/(TP4+FN4)

A_H=(TP1+TP2+TP3+TP4)/length(Y)

PRE_H=(P1_H+P2_H+P3_H+P4_H)/4
REC_H=(R1_H+R2_H+R3_H+R4_H)/4
F_H=2*((PRE_H*REC_H)/(PRE_H+REC_H))

A_H;P1_H;P2_H;P3_H;P4_H;R1_H;R2_H;R3_H;R4_H;PRE_H;REC_H;F_H


library(cowplot)
plot_grid(graph1, graph4_1,graph4_2,graph4_3,
          graph5, graph6, graph7, graph8, ncol = 2,labels = "auto")

setwd("C:/Users/user/Desktop/R")
save.image(file = "Realdata_1007.RData")
