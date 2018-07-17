#' 2D Diffusion-limited Aggregation
#'
#' @param N number of particles to be added to the pre-existing cluster (could be Inf)
#' @param initial a data-frame or matrix with two columns, containing coordinates of a pre-existing cluster.
#' @param abandon.factor a multiplier of the current containing radius of the cluster to get the distance beyond which particles are abandoned
#' @param diameter.tol tolerance fraction of the particle diameter in sticking (penetration)
#' @param write.every a period of regular exporting the current cluster to a file, in addition to writing it last time upon reaching N particles
#' @param write.to path of the file to export the current cluster to. If NULL, there will be no exporting of data.
#' @details Particle diameter is 1. The default value of \code{initial} means that the pre-existing cluster is a single seed particle at the origin.
#'
#' @return a data-frame of cluster coordinates, same as finally written to \code{write.to}
#' @export
#'
dla2d=function(N=200,initial=data.frame(x=0,y=0), abandon.factor=10, diameter.tol=0.01,write.every=3600, write.to=NULL){

tau=2*pi
X=initial[,1]
Y=initial[,2]
R=sqrt(max(X^2+Y^2))

tolerated=1-diameter.tol #tolerated particle diameter
Ab=0

n=0
tic=proc.time()[3]
while (n<N){
    RStart=R+1
    Abandon=RStart*abandon.factor
    phi=runif(1,0,tau)
    x=RStart*cos(phi)
    y=RStart*sin(phi)

    minDist=sqrt(min((X-x)^2+(Y-y)^2))
    Ab=0

    while(minDist>1){
      phi=runif(1,0,tau)
      Rcur=minDist-tolerated
      x=x+Rcur*cos(phi)
      y=y+Rcur*sin(phi)
      minDist=sqrt(min((X-x)^2+(Y-y)^2))


      if(minDist>Abandon){
        Ab=1
        break
      }
    }

    if(Ab==0){
      n=n+1
      X=c(X,x)
      Y=c(Y,y)
      R=max(R,sqrt(x^2+y^2))
    }


    if(proc.time()[3]-tic > write.every){
      if(!is.null(write.to)){
        saveRDS(data.frame(x=X,y=Y),file=write.to)
      }
      tic=proc.time()[3]
    }

  }

df=data.frame(x=X,y=Y)
if(!is.null(write.to)){
  saveRDS(df,file=write.to)
}
return(df)
}

#' 3D Diffusion-limited Aggregation
#'
#' @param N number of particles to be added to the pre-existing cluster (could be Inf)
#' @param initial a data-frame or matrix with two columns, containing coordinates of a pre-existing cluster.
#' @param abandon.factor a multiplier of the current containing radius of the cluster to get the distance beyond which particles are abandoned
#' @param diameter.tol tolerance fraction of the particle diameter in sticking (penetration)
#' @param write.every a period of regular exporting the current cluster to a file, in addition to writing it last time upon reaching N particles
#' @param write.to path of the file to export the current cluster to. If NULL, there will be no exporting of data.
#' @details Particle diameter is 1. The default value of \code{initial} means that the pre-existing cluster is a single seed particle at the origin.
#'
#' @return a data-frame of cluster coordinates, same as finally written to \code{write.to}
#' @export
#'
dla3d=function(N=200,initial=data.frame(x=0,y=0,z=0), abandon.factor=10, diameter.tol=0.01,write.every=3600, write.to=NULL){

  tau=2*pi
  X=initial[,1]
  Y=initial[,2]
  Z=initial[,3]
  R=sqrt(max(X^2+Y^2+Z^2))

  tolerated=1-diameter.tol #tolerated particle diameter
  Ab=0

  n=0
  tic=proc.time()[3]
  while (n<N){
    RStart=R+1
    Abandon=RStart*abandon.factor
    phi=runif(1,0,tau)
    cosTheta=runif(1,-1,1)
    sinTheta=sqrt(1-cosTheta^2)
    x=RStart*cos(phi)*sinTheta
    y=RStart*sin(phi)*sinTheta
    z=RStart*cosTheta

    minDist=sqrt(min((X-x)^2+(Y-y)^2+(Z-z)^2))
    Ab=0

    while(minDist>1){

      Rcur=minDist-tolerated
      phi=runif(1,0,tau)
      cosTheta=runif(1,-1,1)
      sinTheta=sqrt(1-cosTheta^2)
      x=x+Rcur*cos(phi)*sinTheta
      y=y+Rcur*sin(phi)*sinTheta
      z=z+Rcur*cosTheta

      minDist=sqrt(min((X-x)^2+(Y-y)^2+(Z-z)^2))


      if(minDist>Abandon){
        Ab=1
        break
      }
    }

    if(Ab==0){
      n=n+1
      X=c(X,x)
      Y=c(Y,y)
      Z=c(Z,z)
      R=max(R,sqrt(x^2+y^2+z^2))
    }


    if(proc.time()[3]-tic > write.every){
      if(!is.null(write.to)){
        saveRDS(data.frame(x=X,y=Y,z=Z),file=write.to)
      }
      tic=proc.time()[3]
    }

  }

  df=data.frame(x=X,y=Y,z=Z)
  if(!is.null(write.to)){
    saveRDS(df,file=write.to)
  }
  return(df)
}



  # print(proc.time()-ptm)
  # save(X,Y,R,file=paste("dla-",length(X),".RData", sep=""))
  # color=c('black','purple','blue','cyan','green','orange','red')
  #
  # Stage=ceiling(length(X)/NumberOfStagesInPlot)
  # ccc=rep(color[1],Stage)
  # for (k in 2:NumberOfStagesInPlot){
  #   ccc=c(ccc,rep(color[1+(k-1)%%length(color)],Stage))
  #
  # }
  # ccc=ccc[1:length(X)]
  #
  # png(file=paste("DLA-",length(X),"-particles.png", sep=""),width=2000,height=2000)
  #
  # plot(X,Y,col=ccc, pch=".", type='p',asp=1, main=paste("Diffusion-Limited Aggregation of",length(X), "Particles"), xlab='', ylab='',
  #      xaxt='n',yaxt='n', cex=1, bty="n")
  #
  #
  #
