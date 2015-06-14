library(doParallel)
cl<-makeCluster(2)
registerDoParallel(cl)


sds=c(seq(1.2,0.3,by=-0.001),seq(0.250,0.01,by=-0.05))
x1=NULL
x2=NULL
x3=NULL
for(s in sds){
  r=reinit(15)
  r=temp(r,sd=s)
  x=T(r,K=10000,dt=0.01)
  x1=c(x1,x$T)
  x2=c(x2,x$sdE)
  x3=c(x3,x$meanE)
  print(s)
}

frame=data.frame(temperature=x1,sdE=x2,meanE=x3)
frame
plot(data=frame,meanE~temperature)


# foreach(s=sds) %dopar% for(s in sds){
#   r=reinit(15)
#   r=temp(r,sd=s)
#   x=T(r,K=10,dt=0.01)
#   #x1=c(x1,x$T)
#   #x2=c(x2,x$sdE)
#   #x3=c(x3,x$meanE)
#   x
# }
# stopCluster(cl)








