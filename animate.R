s=0.15

r=reinit(13)
r=temp(r,sd=s)
temperature(r,K=10000)


df=animate.calc(r,K=100000)
animate.plot(df,skip=50,xlim=c(-2,2),ylim=c(-2,2))
