source("lib.R")

temperatures2=NULL
radial.sds2=NULL

sds=c(seq(0.001,0.15,by=0.006))



sds=c(seq(0.001,0.15,by=0.001))
#we need adequate sds to cover temperatures from 0 to about 0.5
#that corresponds to sds of:
N=13
K=50000
system.time({
  df=array(dim=c(N,4,K,length(sds)))
  for(s in 1:length(sds)){
    r=reinit(N)
    r=temperature.create(r,sd=sds[s])
    frame=molecular.dynamic(r,K=K,dt=0.01)
    df[,,,s]=frame
  }
})
temperatures=apply(df,MARGIN = 4,temperature3)
radial.sds=apply(df,MARGIN=4,sdR3)
shell1.sds=apply(df,MARGIN=4,function(x) shell.sdR(x,inds=c(2,5,11,12)))
shell2.sds=apply(df,MARGIN=4,function(x) shell.sdR(x,inds = c(1,3,4,6,7,8,9,10,13)))
phi11=apply(df,MARGIN=4,function(x) sdPhi(x,ind1=2,inds2=c(5,11,12)))
phi12=apply(df,MARGIN=4,function(x) sdPhi(x,ind1=2,inds2=c(1,3,4,6,7,8,9,10,13)))
phi21=apply(df,MARGIN=4,function(x) sdPhi(x,ind1=1,inds2=c(5,11,12)))
phi22=apply(df,MARGIN=4,function(x) sdPhi(x,ind1=1,inds2=c(3,4,6,7,8,9,10,13)))

#temperatures2=c(temperatures2,temperatures)
#radial.sds2=c(radial.sds2,radial.sds)

qplot(x=temperatures,y=phi22)
qplot(x=temperatures,y=radial.sds)
#animate.plot(df=df[,,,1])

