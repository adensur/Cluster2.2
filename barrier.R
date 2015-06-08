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
