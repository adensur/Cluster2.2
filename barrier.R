source('lib.R')

Us=NULL
dphis=seq(0.011,2*pi,by=0.03)
r0=reinit(14)
dd=NULL

K=10000
d=dphis[4]
        r=rotate(r0,dphi=d,inds=c(1,7,9,14))
        r2=r

for(as in 1:60){
        r2=search(r2,fixed=c(1,2),K=K)
        print(U(r2))
}
        
        
Us=c(Us,U(r2))
dd=c(dd,d)
qplot(x=dd,y=Us,geom="point")


frame=data.frame(dphis=dd,Us=Us)
#write.csv(frame,"barrier14.csv")

