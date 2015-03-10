init=function(N,sd=1){#generates N particles with x,y,z random coordinates
  x=rnorm(N,sd=sd)
  y=rnorm(N,sd=sd)
  data.frame(x=x,y=y)
}
U=function(r){#scalar function. Returns potential energy of r-type object
  N=dim(r)[1]
  sum=0
  for(i in 1:N){
    sum1=0
    if(i!=N){
      for(j in (i+1):N){
        sum1=sum1+sum((r[i,]-r[j,])^2)^(-6)
      }
    }
    sum=sum1+sum(r[i,]^2)
  }
  sum
}
gU=function(r){#returns gradient of r object in same, N*3, form
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
}
#gradient descent cycle
l=1         #lambda
K=1000     #number of grad descent steps
N=17       #number of particles
r=init(N)
for(k in 1:K){
  gu=gU(r)
  if(U(r-l*gu)<U(r)){
    r=r-l*gu
  }
  else{
    l=l/2
  }
  print(c(k,U(r),l))
}