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