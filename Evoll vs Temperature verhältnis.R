library(doParallel)
library(ggplot2)
##diese Scheiße ist zu langsam

# system.time({
# sds=c(seq(1.5,0.5,by=-0.05),seq(0.45,0.05,by=-0.05))
# x1=NULL
# x2=NULL
# x3=NULL
# for(s in sds){
#   r=reinit(15)
#   r=temp(r,sd=s)
#   x=T(r,K=100,dt=0.01)
#   x1=c(x1,x$T)
#   x2=c(x2,x$sdE)
#   x3=c(x3,x$meanE)
# }
# frame=data.frame(temperature=x1,sdE=x2,meanE=x3)
# })
# max(frame$sdE)
# plot(data=frame,meanE~temperature) 


cl=makeCluster(6)
registerDoParallel(cl)

sds=c(seq(1,0.5,by=-0.005),seq(0.45,0.05,by=-0.05))
system.time({
  frame<-foreach(s=sds,.combine=rbind) %dopar% {
    r=reinit(15)
    r=temp(r,sd=s)
    x=macro(r,K=50000,dt=0.003)
    c(x$T,x$sdE,x$meanE,x$sdR)
  }
})
frame=data.frame(frame)
frame=cbind(frame,sds)
colnames(frame)=c("temperature","sdE","meanE","sdR","sd")
max(frame$sdE)
plot(data=frame,meanE~temperature)
#qplot(data=frame,x=temperature,y=meanE,geom="point")
qplot(data=frame,x=temperature,y=sdR,geom="point")




stopCluster(cl)








