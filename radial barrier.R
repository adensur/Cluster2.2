source('lib.R')

r=reinit(14)
r
plot(r)
r=radial_displace(r,dr=0.1,ind=1)



#
Us=NULL
dd=NULL

drs=seq(0.15,1.3,by=0.1)
r0=reinit(14)
K=10000
d=0
r=radial_displace(r0,dr = d,ind = 1)
r2=r


for(as in 1:60){
        r2=search2(r2,fixed=c(1),K=K)
        print(U(r2))
}

Us=c(Us,U(r2))
dd=c(dd,d)
qplot(x=dd,y=Us,geom="point")

roes=sapply(dd,function(x) {
        ro(radial_displace(r0,dr=x,ind=1))[1]
})

frame=data.frame(rad=roes,U=Us)
write.csv(frame,"radial_barrier14.csv",row.names=FALSE)
