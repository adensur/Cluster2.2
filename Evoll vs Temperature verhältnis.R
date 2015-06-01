
source("lib.R")

cl=makeCluster(1)
registerDoParallel(cl)

sds=c(seq(1,0.5,by=-0.005),seq(0.45,0.05,by=-0.05))
system.time({
  frame<-foreach(s=sds,.combine=rbind) %dopar% {
    r=reinit(13)
    r=temp(r,sd=s)
    r=r[order(ro(r)),]
    x=macro(r,K=100,dt=0.003)
    c(x$T,x$sdE,x$meanE,x$sdR,x$sdPhi78)
  }
})
frame=data.frame(frame)
frame=cbind(frame,sds)
colnames(frame)=c("temperature","sdE","meanE","sdR","sdPhi78","sd")
max(frame$sdE)
#plot(data=frame,meanE~temperature)
qplot(data=frame,x=temperature,y=sdPhi78)
#qplot(data=frame,x=temperature,y=meanE,geom="point",main="Зависимость энергии от температуры", ylab="Полная энергия")+geom_smooth()
#qplot(data=frame,x=temperature,y=sdR,geom=c("point","smooth"),
      main="Зависимость среднеквадратичного отклонения по радиусу для 2х-мерного кластера из 15 частиц",
      ylab="Среднеквадратичное отклонение", )
#saveRDS(list(K=50000,dt=0.003,sds=c(seq(1,0.5,by=-0.005),seq(0.45,0.05,by=-0.05)),frame),"frame1.RDS")



stopCluster(cl)








