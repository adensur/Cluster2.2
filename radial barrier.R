source('lib.R')

r=reinit(13)
r
plot(r)
r=radial_displace(r,dr=0.1,ind=2)



#
Us=NULL

drs=seq(0.1,2,by=0.1)
r0=reinit(13)
K=10000
d=drs[1]
r=radial_displace(r0,dr = d,ind = 2)
r2=r
r2=search2(r2,fixed=c(2),K=K)
U(r2)
r2=search2(r2,fixed=c(2),K=K)
U(r2)
r2=search2(r2,fixed=c(2),K=K)
U(r2)
r2=search2(r2,fixed=c(2),K=K)
U(r2)
r2=search2(r2,fixed=c(2),K=K)
U(r2)
r2=search2(r2,fixed=c(2),K=K)
U(r2)
Us=c(Us,U(r2))
dd=c(dd,d)
qplot(x=dd[-7],y=Us[-7],geom="point")