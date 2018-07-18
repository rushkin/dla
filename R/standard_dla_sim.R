#' Generate DLA cluster
#'
#' Diffusion-limited aggregation (DLA) cluster can be generated in 2- or 3-dimensional Euclidean space.
#'
#' @param N number of particles to be added to the initial cluster (can be Inf)
#' @param initial a data-frame or matrix of coordinates of the initial cluster. The default is a single seed particle at the origin of 2-dimensional plane.
#' The number of columns must be 2 or 3 and serves as the dimensionality of space in which the cluster will be grown.
#' Column names do not matter, only their order.
#' @param source a vector of coordinates of the particle source. A single value Inf is understood, but apart from that, must be the same length as the number of columns in \code{initial}.
#' @param abandon.factor a multiplying coefficient controlling when straying particles should be abandoned. The max cluster distance from the origin or the finite source distance from the origin, whichever is greater, are multiplied by this factor to get the abandoning distance.
#' If a particle strays farther than the abandoning distance from the cluster, it is abandoned.
#' @param diameter.tol a fraction of the particle diameter representing the sticking penetration tolerance.
#' @param report.every if not NULL or Inf, growing the cluster will be split into stages of adding this many particles.
#' @param verbose logical, indicating whether the function should print out progress statements at the end of each stage controlled by \code{report.every}.
#' @param write.to if not NULL, it is the file path to which the function will save the cluster at the end of each stage controlled by \code{report.every}.
#'
#' @return A dataframe containing coordinates of cluster particles. The number of rows will be the number of points in the initial cluster plus \code{N}.
#' @export
#'
#' @examples
dla=function(N=200, initial=data.frame(x=0,y=0), source=Inf, abandon.factor=10, diameter.tol=0.01, report.every=NULL, verbose=TRUE, write.to=NULL){
  tic=proc.time()[3]

  d=ncol(initial)

  if(!(d %in% c(2,3))){
    message('Error. Initial cluster data-frame must have 2 or 3 columns: it is the dimensionality of space.\n')
    return(NULL)
  }

  if(is.null(report.every)) report.every=Inf

  export=(!is.null(write.to))

  df=initial
  N=N+nrow(df)

  if(any(is.infinite(source))){

    if(d==2){
      while(nrow(df)<N){
        df=dla2dinf(N=min(report.every,N-nrow(df)),initial=df,abandon.factor = abandon.factor,diameter.tol=diameter.tol)
        if(export) saveRDS(df,file=write.to)
        if(verbose) cat(paste0('Particles: ',nrow(df),'. Elapsed time: ',round(proc.time()[3]-tic),' sec.\n'))
      }
    }
    if(d==3){
      while(nrow(df)<N){
        df=dla3dinf(N=min(report.every,N-nrow(df)),initial=df,abandon.factor = abandon.factor,diameter.tol=diameter.tol)
        if(export) saveRDS(df,file=write.to)
        if(verbose) cat(paste0('Particles: ',nrow(df),'. Elapsed time: ',round(proc.time()[3]-tic),' sec.\n'))
      }
    }


  }else{

    if(d==2){
      while(nrow(df)<N){
        df=dla2dsource(N=min(report.every,N-nrow(df)),initial=df, source=source, abandon.factor = abandon.factor,diameter.tol=diameter.tol)
        if(export) saveRDS(df,file=write.to)
        if(verbose) cat(paste0('Particles: ',nrow(df),'. Elapsed time: ',round(proc.time()[3]-tic),' sec.\n'))
      }
    }
    if(d==3){
      while(nrow(df)<N){
        df=dla3dsource(N=min(report.every,N-nrow(df)),initial=df, source=source, abandon.factor = abandon.factor,diameter.tol=diameter.tol)
        if(export) saveRDS(df,file=write.to)
        if(verbose) cat(paste0('Particles: ',nrow(df),'. Elapsed time: ',round(proc.time()[3]-tic),' sec.\n'))
      }
    }
  }

  return(df)
}

dla2dinf=function(N=200,initial=data.frame(x=0,y=0), abandon.factor=10, diameter.tol=0.01){

  tau=2*pi
  X=initial[,1]
  Y=initial[,2]
  R=sqrt(max(X^2+Y^2))

  tolerated=1-diameter.tol #tolerated particle diameter
  Ab=0

  n=0

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

  }

  df=data.frame(x=X,y=Y)

  return(df)
}

dla2dsource=function(N=200,initial=data.frame(x=0,y=0), source=c(100,0), abandon.factor=10, diameter.tol=0.01){

  tau=2*pi
  X=initial[,1]
  Y=initial[,2]
  R=sqrt(max(X^2+Y^2))

  tolerated=1-diameter.tol #tolerated particle diameter
  Ab=0

  n=0

  Rsource=sqrt(sum(source^2))

  while (n<N){

      Abandon=max(Rsource,R)*abandon.factor
      x=source[1]
      y=source[2]

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

  }

  df=data.frame(x=X,y=Y)

  return(df)
}

dla3d=function(N=200,initial=data.frame(x=0,y=0,z=0), source=Inf, abandon.factor=10, diameter.tol=0.01,write.every=3600, write.to=NULL){

  tau=2*pi
  X=initial[,1]
  Y=initial[,2]
  Z=initial[,3]
  R=sqrt(max(X^2+Y^2+Z^2))

  tolerated=1-diameter.tol #tolerated particle diameter
  Ab=0

  n=0

  from_inf=any(is.infinite(source))
  Rsource=sqrt(sum(source^2))

  tic=proc.time()[3]
  while (n<N){
    if(from_inf){
      RStart=R+1
      Abandon=RStart*abandon.factor
      phi=runif(1,0,tau)
      cosTheta=runif(1,-1,1)
      sinTheta=sqrt(1-cosTheta^2)
      x=RStart*cos(phi)*sinTheta
      y=RStart*sin(phi)*sinTheta
      z=RStart*cosTheta
    }else{
      Abandon=max(Rsource,R)*abandon.factor
      x=source[1]
      y=source[2]
      z=source[3]
    }

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

    if(!is.null(write.to)){
      if(proc.time()[3]-tic > write.every){
        cat('Particles:',length(X),'\n')
        saveRDS(data.frame(x=X,y=Y,z=Z),file=write.to)
        tic=proc.time()[3]
      }
    }

  }

  df=data.frame(x=X,y=Y,z=Z)
  if(!is.null(write.to)){
    saveRDS(df,file=write.to)
  }
  return(df)
}

dla3dinf=function(N=200,initial=data.frame(x=0,y=0,z=0), abandon.factor=10, diameter.tol=0.01){

  tau=2*pi
  X=initial[,1]
  Y=initial[,2]
  Z=initial[,3]
  R=sqrt(max(X^2+Y^2+Z^2))

  tolerated=1-diameter.tol #tolerated particle diameter
  Ab=0

  n=0

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

  }

  df=data.frame(x=X,y=Y,z=Z)

  return(df)
}

dla3dsource=function(N=200,initial=data.frame(x=0,y=0,z=0), source=c(100,0,0), abandon.factor=10, diameter.tol=0.01){

  tau=2*pi
  X=initial[,1]
  Y=initial[,2]
  Z=initial[,3]
  R=sqrt(max(X^2+Y^2+Z^2))

  tolerated=1-diameter.tol #tolerated particle diameter
  Ab=0

  n=0

  Rsource=sqrt(sum(source^2))

  while (n<N){

      Abandon=max(Rsource,R)*abandon.factor
      x=source[1]
      y=source[2]
      z=source[3]

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

  }

  df=data.frame(x=X,y=Y,z=Z)
  return(df)
}

#' Visualize DLA cluster
#'
#'
#' @param data a dataframe or matrix with 2 or 3 columns of particle coordinates (for a 2D or 3D cluster), like the one produced by \code{dla2d} or \code{dla3d}.
#' @param color0 color for the earliest particles
#' @param color1 color for the latest particles
#'
#' @return
#' @export
#'
vis.dla=function(data,color0='green',color1='red'){

  if(ncol(data)==2){

    title=paste("Diffusion-Limited Aggregation of",nrow(data), "Particles")

    plot(data[,1],data[,2],col=colorRampPalette(c(color0,color1))(nrow(data)), pch=".", type='p',asp=1, main=title, xlab='', ylab='',
       xaxt='n',yaxt='n', cex=1, bty="n")
  }

  if(ncol(data)==3){
    title=paste("Diffusion-Limited Aggregation of",nrow(data), "Particles")
    scatterplot3d::scatterplot3d(data[,1],data[,2],data[,3],color=colorRampPalette(c(color0,color1))(nrow(data)),pch=".", asp=1, main=title, xlab="",ylab="", zlab="", grid=FALSE,axis=TRUE, tick.marks=FALSE, cex.symbols=3)
  }

}
