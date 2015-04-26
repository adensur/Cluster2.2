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
  vx=rnorm(N,sd=sd)
  vx=(vx-sum(vx)/length(vx))
  vy=rnorm(N,sd=sd)
  vy=(vy-sum(vy)/length(vy))
  data.frame(x=r$x,y=r$y,vx=vx,vy=vy)
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
T=function(r,K=10000,dt=0.01){    #returns temperature as an average across velocities + 
                                    #full energy for troubleshooting
  N=dim(r)[1]
  Es=NULL
  Vs=NULL
  for(i in 1:K){
    r[,1:2]=r[,1:2]+dt*r[,3:4]      #r step
    r[,3:4]=r[,3:4]-dt*gU(r)        #v step
    Es=c(Es,E(r))
    Vs=c(Vs,sqrt(r$vx^2+r$vy^2))
  }
  list(T=mean(Vs^2),sdE=sd(Es),meanE=mean(Es))
}

#checking deltaE ~ sd dependance:
sd=seq(0,100,by=0.05)
Es=sapply(sd,function(x) E(temp(r,x)))
plot(sd,Es,type="l")

#gradient descent cycle
for(n in 1:40){
l=0.5         #lambda
K=100000     #number of grad descent steps
N=n      #number of particles
r=init(N)
for(k in 1:K){
  gu=gU(r)
  if(U(r-l*gu)<U(r)){
    r=r-l*gu
    l=0.5
  }else{
    l=l/5
  }
  #print(c(k,U(r),l))
}
filenum=paste("data/",N,".RDS",sep="")
saveRDS(r,filenum)
print(n)
}

#Molecular dynamic
#-----
K=3000    #number of iterations
dt=0.0001       #delta t (should be small for more precision)
dt1=0.001
dt2=0.01
r=reinit(27)
r=temp(r,sd=0.5)
r1=r
r2=r
plot(r$x,r$y)
Es=NULL
Es1=NULL
Es2=NULL
for(i in 1:K){
  r[,1:2]=r[,1:2]+dt*r[,3:4]      #r step
  r[,3:4]=r[,3:4]-dt*gU(r)        #v step
  r1[,1:2]=r1[,1:2]+dt1*r1[,3:4]      #r step
  r1[,3:4]=r1[,3:4]-dt1*gU(r1)        #v step
  r2[,1:2]=r2[,1:2]+dt2*r2[,3:4]      #r step
  r2[,3:4]=r2[,3:4]-dt2*gU(r2)        #v step
  #print(E(r))
  #points(r$x,r$y)
  Es=c(Es,E(r))
  Es1=c(Es1,E(r1))
  Es2=c(Es2,E(r2))
}

  
  
  
  
  
  
  
  
  
  



