source('lib.R')

r=reinit(13)
r
plot(r)
r=radial_displace(r,dr=0.1,ind=2)



#
Us=NULL

drs=seq(0.15,1.3,by=0.1)
r0=reinit(13)
K=10000
d=drs[5]
r=radial_displace(r0,dr = d,ind = 2)
r2=r


for(as in 1:60){
        r2=search2(r2,fixed=c(2),K=K)
        U(r2)
}

Us=c(Us,U(r2))
dd=c(dd,d)
qplot(x=dd,y=Us,geom="point")

roes=sapply(dd,function(x) {
        ro(radial_displace(r0,dr=x,ind=2))[2]
})

frame=data.frame(rad=roes,U=Us)
write.csv(frame,"radial_barrier13.csv",row.names=FALSE)