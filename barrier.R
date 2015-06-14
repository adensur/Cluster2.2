source('lib.R')

r=reinit(13)
r

system.time({
  r2=search(r,fixed=c(1,2),K=300000)
})
print(U(r2))

system.time({
  Us=NULL
  dphis=seq(0.01,10,by=0.01)
  r0=reinit(13)
  K=100000
  for(d in dphis){
    r=rotate(r0,dphi=d,inds=c(2,5,11,12))
    r2=search(r,fixed=c(1,2),K=K)
    Us=c(Us,U(r2))
  }
})
qplot(x=dphis,y=Us,geom="line")


plot(r0)
for(d in dphis){
  r=rotate(r0,dphi=d,inds=c(2,5,11,12))
  points(r)
}


#sapply
r0=reinit(13)
dphis=seq(0,pi,by=0.05)
K=500000
system.time({
        U2=sapply(dphis,function(d) {
                U(search(rotate(r0,dphi=d,inds=c(2,5,11,12)),fixed=c(1,2),K=K))
        })
})

qplot(x=dphis,y=U2,geom="line")

#manual
Us=NULL
dphis=seq(0.01,2*pi,by=0.1)
r0=reinit(14)



K=50000
d=dphis[4]
        r=rotate(r0,dphi=d,inds=c(2,5,11,12))
        r2=r
        r2=search(r2,fixed=c(1,2),K=K)
        U(r2)
        r2=search(r2,fixed=c(1,2),K=K)
        U(r2)
        r2=search(r2,fixed=c(1,2),K=K)
        U(r2)
        r2=search(r2,fixed=c(1,2),K=K)
        U(r2)
        r2=search(r2,fixed=c(1,2),K=K)
        U(r2)
        r2=search(r2,fixed=c(1,2),K=K)
        U(r2)
        Us=c(Us,U(r2))
dd=c(dd,d)
qplot(x=dd[-7],y=Us[-7],geom="point")


frame=data.frame(dphis=dd,Us=Us)
write.csv(frame,"barrier.csv")

d=dphis[52]
r0=reinit(13)
K=10000

r=rotate(r0,dphi=d,inds=c(2,5,11,12))
r2=r
r2=search(r2,fixed=c(1,2),K=K)
U(r2)
U3=U()