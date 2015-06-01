library(doParallel)
library(ggplot2)
init=function(N,sd=1){#generates N particles with x,y,z random coordinates
  x=rnorm(N,sd=sd)
  y=rnorm(N,sd=sd)
  data.frame(x=x,y=y)
}
U=function(r){#scalar function. Returns potential energy of r-type object
  N=dim(r)[1]
  ind=1
  for(i in 2:N){
    ind[i]=N*(i-1)+i
  }
  frame=data.frame(id="0",rx=r$x,ry=r$y)
  x=merge(frame,frame,by="id")
  x=x[-ind,]
  y=(x[,4]-x[,2])^2+(x[,5]-x[,3])^2
  sum(y^(-3))/2+sum(r^2)
}
U3=function(r){
  sum(mapply(function(xi,yi) {
    sum(unlist(mapply(function(xk,yk) {
      if(xi!=xk & yi!=yk){
        ((xi-xk)^2+(yi-yk)^2)^(-3)
      }
    },r$x,r$y)))
  },r$x,r$y))/2+sum(r$x^2)+sum(r$y^2)
}#U rewritten as mapply(). It is not faster!
U4=function(r){
  N=dim(r)[1]
  if(N>2){
    sum(sapply(1:(N-1),function(i) {
      sum(sapply((i+1):N,function(k){
        ((r$x[i]-r$x[k])^2+(r$y[i]-r$y[k])^2)^(-3)
      }))
    }))+sum(r$x^2)+sum(r$y^2)
  }
  
}#U rewritten as sapply(). It is not faster!
U2=function(r){#scalar function. Returns potential energy of r-type object
  N=dim(r)[1]
  sum=0
  for(i in 1:N){
    sum1=0
    if(i!=N){
      for(j in (i+1):N){
        sum1=sum1+sum((r[i,]-r[j,])^2)^(-3)
      }
    }
    sum=sum+sum1+sum(r[i,]^2)
  }
  sum
}#initial one
gU=function(r){#returns gradient of r object in same, N*3, form
  N=dim(r)[1]
  frame=data.frame(x=0,y=0)
  for(i in 1:N){
    for(j in 1:2){
      r2=r[-i,]
      y=(r[i,j]-r2[,j])/(((r2[,1]-r[i,1])^2+(r2[,2]-r[i,2])^2)^4)
      frame[i,j]=2*r[i,j]-6*sum(y)
    }
  }
  frame
}
gU2=function(r){#returns gradient of r object in same, N*3, form
  N=dim(r)[1]
  frame=data.frame(x=0,y=0)
  for(i in 1:N){
    for(j in 1:2){
      sum=0
      for(k in 1:N){
        if(k!=i){
          sum=sum+(r[i,j]-r[k,j])/(sum((r[i,]-r[k,])^2)^4)
        }
      }
      frame[i,j]=2*r[i,j]-6*sum
    }
  }
  frame
}#initial one
reinit<-function(N=17){
  readRDS(file=paste("data/",N,".RDS",sep=""))
}
#temperature function
temp<-function(r, sd=1){
  N=dim(r)[1]
  dx=rnorm(N,sd=sd)
  dy=rnorm(N,sd=sd)
  data.frame(x=r$x-dx,y=r$y-dy,vx=0,vy=0)
}
E<-function(r){
  N=dim(r)[1]
  U(r[,1:2])+sum(r$vx^2)/2+sum(r$vy^2)/2
}
descent=function(r,lambda=0.3,K=1000,print=FALSE){
  l=lambda
  for(k in 1:K){
    gu=gU(r)
    if(U(r-l*gu)<U(r)){
      r=r-l*gu
      l=lambda
    }else{
      l=l/5
    }
    if(print)print(U(r))
  }
  r
}
macro=function(r,K=10000,dt=0.01){    #returns temperature as an average across velocities + 
  #full energy for troubleshooting
  N=dim(r)[1]
  Es=NULL                   #volle Energie
  Vs=NULL                   #Geschwindigkeiten
  Rs=NULL                   #Radius
  phi=NULL
  for(i in 1:K){
    r[,1:2]=r[,1:2]+dt*r[,3:4]      #r step
    r[,3:4]=r[,3:4]-dt*gU(r)        #v step
    Es=c(Es,E(r))
    Vs=c(Vs,sqrt(r$vx^2+r$vy^2))
    Rs=rbind(Rs,sqrt(r$x^2+r$y^2))
    phi=rbind(phi,phiAll(r))
  }
  list(T=mean(Vs^2),sdE=sd(Es),meanE=mean(Es),sdR=mean(apply(Rs,2,sd)),sdPhi=mean(apply(phi,2,sd)))
}
temperature<-function(r,K=1000,dt=0.01){
  #returns temperature as an average across velocities + 
  #full energy for troubleshooting
  N=dim(r)[1]
  Vs=NULL                   #Geschwindigkeiten
  for(i in 1:K){
    r[,1:2]=r[,1:2]+dt*r[,3:4]      #r step
    r[,3:4]=r[,3:4]-dt*gU(r)        #v step
    Vs=c(Vs,sqrt(r$vx^2+r$vy^2))
  }
  list(T=mean(Vs^2))
}
animate.calc<-function(r,K=1000,dt=0.01){
  Es=NULL
  df=array(dim=c(dim(r),K))
  for(i in 1:K){
    r[,1:2]=r[,1:2]+dt*r[,3:4] #r schritt
    r[,3:4]=r[,3:4]-dt*gU(r)   #v schritt
    df[,,i]=as.matrix(r)
  }
  df
}
animate.plot<-function(df,skip=10,xlim,ylim){
  rr=data.frame(df[,,1])
  plot(rr[,1],rr[,2],xlim=xlim,ylim=ylim)
  K=dim(df)[3]
  for(i in 1:(K%/%skip)){
    rr=data.frame(df[,,i*skip])
    points(rr[,1],rr[,2])
  }
}
phi<-function(r,i,j){#berechnet der winkel zwischen zwei vektoren
  acos((r$x[i]*r$x[j]+r$y[i]*r$y[j])/(sqrt(r$x[i]^2+r$y[i]^2))/sqrt(r$x[j]^2+r$y[j]^2))
}
phiAll<-function(r){
  df=cbind(0,r[,1:2],1:dim(r)[1])
  df=merge(df,df,by="0")
  df=df[df[,4]!=df[,7],]
  acos((df$x.x*df$x.y+df$y.x*df$y.y)/(sqrt(df$x.x^2+df$y.x^2))/sqrt(df$x.y^2+df$y.y^2))
}
ro<-function(r){#berechnet der radius von allen partikeln
  sqrt(r$x^2+r$y^2)
}