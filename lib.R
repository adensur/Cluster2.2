library(doParallel)
library(ggplot2)
rotate<-function(r,dphi=0.001,inds){
  #rotates one shell in respect to another;
  #shell is specified by inds
  r2=r
  r2[inds,1]=r[inds,1]*cos(dphi)-r[inds,2]*sin(dphi)
  r2[inds,2]=r[inds,1]*sin(dphi)+r[inds,2]*cos(dphi)
  r2
}
search<-function(r,fixed,K=10000){
  #finds a potential minimum for a cluster, where 'fixed' specifies indices of the particles, that had to be fixed 
  #(their angle is held constant, and only radius varies)
  lambda=1
  N=dim(r)[1]
  U1=U(r)
  for(i in 1:K){
    r2=r
    r2[-fixed,1]=r2[-fixed,1]+lambda*rnorm(N-2)
    r2[-fixed,2]=r2[-fixed,2]+lambda*rnorm(N-2)
    r2[fixed,]=r2[fixed,]*(1+lambda*rnorm(2))
    U2=U(r2)
    if(U2<U1){
      r=r2
      U1=U2
      lambda=1
    }else{
      lambda=lambda*0.9
    }
  }
  r
}

init=function(N,sd=1){#generates N particles with x,y,z random coordinates
  x=rnorm(N,sd=sd)
  y=rnorm(N,sd=sd)
 m=matrix(c(x,y),nrow=N,ncol=2)
  colnames(m)=c("x","y")
 m
}
U=function(r){#scalar function. Returns potential energy of r-type object
  N=dim(r)[1]
  ind=1
  for(i in 2:N){
    ind[i]=N*(i-1)+i
  }
  frame=data.frame(id="0",rx=r[,1],ry=r[,2])
  x=merge(frame,frame,by="id")
  x=x[-ind,]
  y=(x[,4]-x[,2])^2+(x[,5]-x[,3])^2
  sum(y^(-3))/2+sum(r^2)
}
gU=function(r){#returns gradient of r object in same, N*3, form
  N=dim(r)[1]
  frame=matrix(nrow=N,ncol=2)
  for(i in 1:N){
    for(j in 1:2){
      r2=r[-i,]
      y=(r[i,j]-r2[,j])/(((r2[,1]-r[i,1])^2+(r2[,2]-r[i,2])^2)^4)
      frame[i,j]=2*r[i,j]-6*sum(y)
    }
  }
  colnames(frame)=c("x","y")
  frame
}
reinit<-function(N=17){
  as.matrix(readRDS(file=paste("data/",N,".RDS",sep="")))
}
temperature.create<-function(r, sd=1){
  #gives some energy to system by assigning random displacements to the radius-vectors of the particles;
  #finite temperature depends on sd, but has random component as well
  N=dim(r)[1]
  dx=rnorm(N,sd=sd)
  dy=rnorm(N,sd=sd)
  m=matrix(c(r[,1]-dx,r[,2]-dy,rep(0,times=2*N)),nrow=N,ncol=4)
  colnames(m)=c("x","y","vx","vy")
  m
}
E<-function(r){
  #calculates full energy for a N*4 object. If an N*2 object is given, returns an error
  U(r[,1:2])+sum(r[,3]^2)/2+sum(r[,4]^2)/2
}
descent=function(r,lambda=0.3,K=1000,print=FALSE){
  #calculates gradient descent for a given N*2 vector r;
  #returns a N*2 vector of radiuses of system, corresponding (hopefully) to an equilibrium state
  l=lambda
  for(k in 1:K){
    gu=gU(r)
    if(U(r-l*gu)<U(r)){
      r=r-l*gu
      l=lambda
    }else{
      l=l/5
    }
    if(print)print(U(r))
  }
  r
}
molecular.dynamic<-function(r,K=10000,dt=0.01){
  #calculates a 3-dimensional array corresponding to consequent states of system.
  #complete information is preserved
  df=array(dim = c(dim(r),K))
  for(i in 1:K){
    r[,1:2]=r[,1:2]+dt*r[,3:4]      #r step
    r[,3:4]=r[,3:4]-dt*gU(r)        #v step
    df[,,i]=r
  }
  df
}
temperature3<-function(df){
  #3d-array function, calculating temperature across an array df of consequent states of the system,
  #as provided by molecular.dynamic function.
  #if input df has dimensions N*4*K, where N - number of particles, K - number of MD steps,
  #the function returns 1 value of temperature. Bigger K give more precise temperature
  velocity<-function(r){
    sqrt(r[,3]^2+r[,4]^2)
  }
  velocities=unclass(apply(df,MARGIN=3,velocity))
  result=mean(velocities^2)
  names(result)="Temperature"
  result
}
sdE3<-function(df){
  #3d-array function of energy standard deviation;
  #receives a df object of size N*4*K, where N is the number of particles and K is number of MD steps,
  #as provided by molecular.dynamic() function.
  energies=apply(df,MARGIN=3,E)
  result=sd(energies)
  names(result)="energy.sd"
  result
}
animate.plot<-function(df,skip=10,xlim=c(-2,2),ylim=c(-2,2)){
  rr=data.frame(df[,,1])
  plot(rr[,1],rr[,2],xlim=xlim,ylim=ylim)
  K=dim(df)[3]
  for(i in 1:(K%/%skip)){
    rr=data.frame(df[,,i*skip])
    points(rr[,1],rr[,2])
  }
}


sdR3<-function(df){
  ro<-function(r){#berechnet der radius von allen partikeln
    sqrt(r[,1]^2+r[,2]^2)
  }
  radii=apply(df,MARGIN=3,ro)
  sds=apply(radii,MARGIN=1,sd)
  result=mean(sds)
  names(result)="radial.sd"
  result
}
shell.sdR<-function(df,inds){
  #calculates sd of radii of particles within certain shell (specified by indices of the particles - inds)
  radii=apply(df,MARGIN=3,ro)
  sds=apply(radii[inds,],MARGIN=1,sd)
  result=mean(sds)
  names(result)="shell.radial.sd"
  result
}
ro<-function(r){#berechnet der radius von allen partikeln
  sqrt(r[,1]^2+r[,2]^2)
}
sdPhi<-function(df,ind1,inds2){
  phis=apply(df,MARGIN=3,function(x) phi(x,ind1,inds2))
  sds=apply(phis,MARGIN=1,sd)
  result=mean(sds)
  names(result)="shell.phi.sd"
  result
}

phi<-function(r,i,j){#berechnet der winkel zwischen zwei vektoren
  acos((r[i,1]*r[j,1]+r[i,2]*r[j,2])/(sqrt(r[i,1]^2+r[i,2]^2))/sqrt(r[j,1]^2+r[j,2]^2))
}
# phiAll<-function(r){
#   df=cbind(0,r[,1:2],1:dim(r)[1])
#   df=merge(df,df,by="0")
#   df=df[df[,4]!=df[,7],]
#   acos((df$x.x*df$x.y+df$y.x*df$y.y)/(sqrt(df$x.x^2+df$y.x^2))/sqrt(df$x.y^2+df$y.y^2))
# }
