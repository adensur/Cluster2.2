s=0.01

r=reinit(13)
r=temp(r,sd=s)
temperature(r,K=10000)
#now we need T=0.01

df=animate.calc(r,K=1000000)
animate.plot(df,skip=50,xlim=c(-2,2),ylim=c(-2,2))
