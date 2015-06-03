source("lib.R")

temperatures2=NULL
radial.sds2=NULL

sds=c(seq(0.001,0.15,by=0.006))



sds=c(seq(0.001,0.15,by=0.001))
#we need adequate sds to cover temperatures from 0 to about 0.5
#that corresponds to sds of:
N=14
K=50000
system.time({
  df=array(dim=c(N,4,K,length(sds)))
  for(s in 1:length(sds)){
    r=reinit(N)
    r=temperature.create(r,sd=sds[s])
    frame=molecular.dynamic(r,K=K,dt=0.01)
    df[,,,s]=frame
    gc()
  }
})
temperatures=apply(df,MARGIN = 4,temperature3)
radial.sds=apply(df,MARGIN=4,sdR3)
shell1.sds=apply(df,MARGIN=4,function(x) shell.sdR(x,inds=c(1,7,9,14)))
shell2.sds=apply(df,MARGIN=4,function(x) shell.sdR(x,inds = c(2,3,4,5,6,8,10,11,12,13)))
phi11=apply(df,MARGIN=4,function(x) sdPhi(x,ind1=1,inds2=c(7,9,14)))
phi12=apply(df,MARGIN=4,function(x) sdPhi(x,ind1=1,inds2=c(2,3,4,5,6,8,10,11,12,13)))
phi21=apply(df,MARGIN=4,function(x) sdPhi(x,ind1=2,inds2=c(7,9,14)))
phi22=apply(df,MARGIN=4,function(x) sdPhi(x,ind1=2,inds2=c(3,4,5,6,8,10,11,12,13)))

#temperatures2=c(temperatures2,temperatures)
#radial.sds2=c(radial.sds2,radial.sds)

qplot(x=temperatures,y=phi12,xlim=c(0,2))
qplot(x=temperatures,y=radial.sds,xlim=c(0,2),ylim=c(0,0.7),main="Температурная зависимость среднеквадратичного отклонения 
      расстояния от частиц до центра; все частицы; кластер с N=14",xlab="Температура", ylab="Среднеквадратичное отклонение расстояния от частиц до центра")
#animate.plot(df=df[,,,1])

frame=data.frame(temperature=temperatures,radial.sd=radial.sds,inner.shell.radial.sds=shell1.sds,outer.shell.radial.sds=shell2.sds,
                 phi11=phi11,phi12=phi12)
write.csv(frame,"frame14.csv")
qplot(data=frame,x=temperatures,y=radial.sd)
saveRDS(df,"df13K50000.RDS")














