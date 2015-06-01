source("lib.R")

#checking deltaE ~ sd dependance:
sd=seq(0,100,by=0.05)
Es=sapply(sd,function(x) E(temp(r,x)))
plot(sd,Es,type="l")

#gradient descent cycle
for(n in 58:100){
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
K=300    #number of iterations
dt=0.001       #delta t (should be small for more precision)
r=reinit(15)
r=temp(r,sd=0.5)
plot(r$x,r$y)
Es=NULL
vs=NULL
for(i in 1:K){
  r[,1:2]=r[,1:2]+dt*r[,3:4]      #r step
  r[,3:4]=r[,3:4]-dt*gU(r)        #v step
  #print(E(r))
  #points(r$x,r$y)
  vs=c(vs,sqrt(r$vx^2+r$vy^2))
  Es=c(Es,E(r))
}

  
  
  
  
  
  
  
  
  
  



