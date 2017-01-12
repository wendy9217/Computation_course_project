#install.packages("‘linprog")
require(linprog)
#install.packages("lpSolve")
library(lpSolve)
#install.packages("Rlab")
library("Rlab")
#install.packages("igraph")
library(igraph)



# L1 relaxation optimization function under our setting using integer optimization
optimal.p=function(R,ntilde,number.nodes.active=3){
  f.obj=rep(1,ntilde)%*%R
  f.con=rbind(rep(1,ntilde),diag(1,ntilde),diag(1,ntilde))
  f.dir=rep(c("<=","<=",">="),c(1,ntilde,ntilde))
  f.rhs=rep(c(number.nodes.active,1,0),c(1,ntilde,ntilde))
  return(lp(direction="max",f.obj,f.con,f.dir,f.rhs,int.vec=1:ntilde)$solution)
}
#optimal.p(R,number.nodes,number.nodes.active=4)

#function using L1 relaxation optimization using random generated structure matrix
random.optimal.p=function(number.nodes,alpha=0.1,lambda=0.1,t=0.02,number.nodes.active=4){
  T=matrix(runif(number.nodes^2),number.nodes)
  diag(T)=0
  T[T<=.3]=0
  f=matrix(0,nrow=number.nodes,ncol=number.nodes)
  P=solve((1+lambda)*diag(number.nodes)-T)
  for(i in 1:number.nodes){
    for (j in 1:number.nodes){
      f[i,j]=alpha/P[i,i]*P[j,i]
    }
  }
  R=diag(0,number.nodes)
  R[f>t]=1
  Ct.matrix=matrix(1,number.nodes,number.nodes)
  Ct.matrix[T==0]=0
  return(list(optimal.p(R,number.nodes,number.nodes.active),Ct.matrix))
}

r1=random.optimal.p(number.nodes=5000,alpha=0.1,lambda=0.1,t=0.02,number.nodes.active=300)[[1]]
names(r1)=paste("Individual",1:5)
# visulation using barplot
barplot(r1,col="red",main="Best Individuals Selected with 5 individuals",yaxt="n")
# visulation over network
ad.net=graph.adjacency(r1[[2]])
plot(ad.net,layout=layout.circle,vertex.color=c("yellow","lightblue")[r1[[1]]+1],
     frame=F,vertex.size=18,margin=-.2)

#data & constant
n=100
m=4
C=matrix(runif(n*m), ncol=m)
Ctilde=matrix(0,ncol=m,nrow=n)
thred=0.5
Ia=3
ntilde=30

#Ctilde
Ctilde[C>thred]=1


#r(u)
r=c(rep(0,n))
for(i in 1:n){
  no_intersect=sum(Ctilde[i,1:Ia])
  no_union=sum(Ctilde[1,(Ia+1):m])+Ia
  r[i]=no_intersect/no_union
}

#H(U)
#test
p_u=rbern(n, 0.5)
p=c(1,1,0,0,1)
H0<-function(p,Ia,Ctilde){
  H=0
  for(i in 1:Ia){
    p_hat=(p%*%Ctilde[,i])/sum(p)
    H=H+p_hat*log(p_hat+10^(-5))
  }
  return((-H/log(Ia))[1,1])
}
#greedyOpt
lambda=0.5
#number.nodes=n
#number.active=ntilde
#r.func.value=r
#function for greedy algorithm without recoding each steps
greedyOpt<-function(number.active,number.nodes,r.func.value,Ia,Ctilde,lambda=0.5){
  p_u=rep(FALSE,number.nodes)
  utility.func=0
  for(i in 1:number.active){
    utility.best=-999
    p_u_loc_opt=p_u
    for(j in (1:number.nodes)[!p_u]){
      p_u.temp=p_u
      p_u.temp[j]=TRUE
      utility.func=lambda*sum(r.func.value[p_u.temp])+(1-lambda)*H0(p_u.temp,Ia,Ctilde)
      if(utility.func>utility.best){
        utility.best=utility.func
        p_u_loc_opt=p_u.temp
      }
    }
    p_u=p_u_loc_opt      
  }
  return(p_u)
}
r2=greedyOpt(number.active=3,number.nodes,r.func.value,Ia,Ctilde,lambda=0.5)
r2=as.numeric(r2)
names(r1)=paste("Individual",1:5)
barplot(r1,col="red",main="Best Individuals Selected with 5 individuals",yaxt="n")

#Best_selected_15_L1_relaxation


#function for greedy algorithm with recoding each steps
greedyOpt1<-function(number.active,number.nodes,r.func.value,Ia,Ctilde,lambda=0.5){
  p_u=rep(FALSE,number.nodes)
  utility.func=0
  pp.matrix=as.numeric(p_u)
  for(i in 1:number.active){
    utility.best=-999
    p_u_loc_opt=p_u
    for(j in (1:number.nodes)[!p_u]){
      p_u.temp=p_u
      p_u.temp[j]=TRUE
      utility.func=lambda*sum(r.func.value[p_u.temp])+(1-lambda)*H0(p_u.temp,Ia,Ctilde)
      if(utility.func>utility.best){
        utility.best=utility.func
        p_u_loc_opt=p_u.temp
      }
    }
    pp.matrix=rbind(pp.matrix,as.numeric(p_u_loc_opt))
    p_u=p_u_loc_opt      
  }
  return(list(p_u,pp.matrix))
}
r2=greedyOpt1(number.active=3,number.nodes,r.func.value,Ia,Ctilde,lambda=0.5)[[2]]
colnames(r2)=paste("Individual",1:5)
#generate Gif to show how the algorithm works
saveGIF({
  for(i in 1:4){
    r.temp=r2[i,]
    barplot(r.temp,col="red",main="Greedy Algorithm",yaxt="n",ylim=c(0,1))
  }
},interval=1, ani.width = 1600, ani.height = 800)

#more on visualization

y1=c(0.42,0.67,0.88,0.92,0.96)
y2=c(0.27,0.34,0.75,0.86,0.89)
x=c(5,10,15,20,25)
plot(y1,type="o",col="blue",ylim=c(0,1), xlab="n", ylab="average cover ratio")
lines(y2,type="o", pch=22, lty=2, col="red")  
legend("bottomright",legend=c("without screening","with screening"),col=c("blue","red"),pch=c(2,12),lty=1)

legend(1,0.8, cex=0.8, 
       col=c("blue","red"), pch=21:22, lty=1:2))

y1=c(5.4,5.9,6.6)
y2=c(6.2,6.7,7.0)
y3=c(6.9,7.4,7.9)
y4=c(7.7,8.9,9.7)
dataxx=cbind(y1,y2,y3,y4)
colnames(dataxx)=paste("n=", c(1000,1500,2000,2500), sep = "")
dataxx=as.data.frame(t(dataxx))
cols=terrain.colors(4)
labels <- paste("n=", c(1000,1500,2000,2500), sep = "")
arg=c(1000,1200,1500,2000)
rownames(dataxx)=paste("n=", c(1000,1500,2000,2500), sep = "")
barplot(t(dataxx),ylim=c(0,16),offset=0,axis.lty=1,names.arg=arg,
        beside=TRUE,col=terrain.colors(3))
barplot(dataxx,col=terrain.colors(3),beside=TRUE,ylim=c(0,16),
        ylab="Time/10e-3s")


legend(1,17,"topleft", legend = paste("n/",c(50,80,100),sep=""), 
       fill = terrain.colors(3), box.col = "transparent",cex=0.7)



