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
  data.frame(x=r$x+rnorm(N,sd=sd),y=r$y+rnorm(N,sd=sd))
}

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



