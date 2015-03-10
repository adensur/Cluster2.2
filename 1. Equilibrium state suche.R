init=function(N,sd=1){#generates N particles with x,y,z random coordinates
  x=rnorm(N,sd=sd)
  y=rnorm(N,sd=sd)
  z=rnorm(N,sd=sd)
  data.frame(x=x,y=y,z=z)
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
  frame=data.frame(x=0,y=0,z=0)
  for(i in 1:N){
    for(j in 1:3){
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