#' Generate DLA cluster
#'
#' Diffusion-limited aggregation (DLA) cluster can be generated in 2- or 3-dimensional Euclidean space.
#'
#' @param N number of particles to be added to the initial cluster (can be Inf)
#' @param initial a data-frame or matrix of coordinates of the initial cluster. For a single seed particle, a coordinate vector will be also accepted. The default is a single seed particle at the origin of 2-dimensional plane.
#' The number of columns must be 2 or 3 and serves as the dimensionality of space in which the cluster will be grown.
#' Column names do not matter, only their order.
#' @param source a vector of coordinates of the particle source. A single value Inf is understood, but apart from that, must be the same length as the number of columns in \code{initial}.
#' @param abandon.factor a multiplying coefficient controlling when straying particles should be abandoned. The max cluster distance from the origin or the finite source distance from the origin, whichever is greater, are multiplied by this factor to get the abandoning distance.
#' If a particle strays farther than the abandoning distance from the cluster, it is abandoned.
#' @param diameter.tol a fraction of the particle diameter representing the sticking penetration tolerance.
#' @param checkin.every if not NULL or Inf, growing the cluster will be split into stages of adding this many particles.
#' At the end of each stage the function will check whether the source has been blocked by the cluster (if yes - will message stop).
#' It may also print a progress message and export cluster to a file, as controlled by the arguments \code{verbose} and \code{write.to}.
#' @param verbose logical, indicating whether the function should print out progress statements at the end of each stage controlled by \code{checkin.every}.
#' @param write.to if not NULL, it is the file path to which the function will export the cluster at the end of each stage controlled by \code{checkin.every}.
#' @param phi,theta two vectors of length 2, denoting the boundaries of the sector of particle propagation. \code{phi} is the azimuthal angle (longitude in 3D) and \code{theta} is the polar angle of spherical coordinates in 3D (the pole being 0).
#' @return A dataframe containing coordinates of cluster particles. The number of rows will be the number of points in the initial cluster plus \code{N}.
#' @export
#'
dla=function(N=5000, initial=data.frame(x=0,y=0), source=Inf, abandon.factor=10, diameter.tol=0.01, checkin.every=5000, verbose=TRUE, write.to=NULL
             ,flow=FALSE
             ,exit_density=NULL
             ,phi=c(0,2*pi),theta=c(0,pi)){
  tic=proc.time()[3]

  if(is.null(dim(initial))){
    initial=t(matrix(initial))
  }
  phi=sort(phi)
  costheta=sort(cos(theta))

  d=ncol(initial)
  if(is.null(checkin.every)) checkin.every=Inf

  if(!(d %in% c(2,3))){
    message('Error. Initial cluster data-frame must have 2 or 3 columns: it is the dimensionality of space.\n')
    return(NULL)
  }

  if(is.infinite(N)&is.infinite(checkin.every)){
    message('Error. Number of particles N is infinite, but checkin.every is not set to a finite number. This would cause the function to grow the cluster infinitely without exporting. \n')
    return(NULL)
  }

  if(is.infinite(N)&is.null(write.to)){
    message('Error. Number of particles N is infinite, but write.to is NULL. This would cause the function to grow the cluster infinitely without exporting. \n')
    return(NULL)
  }

  if(diameter.tol<=0){
    message('Error. The diameter tolerance needs to be positive. \n')
    return(NULL)
  }


  export=(!is.null(write.to))

  df=initial
  N=N+nrow(df)

  if(any(is.infinite(source))){



    if((d==2)&(!flow)){
      while(nrow(df)<N){
        df=dla2dinf(N=min(checkin.every,N-nrow(df)),initial=df,abandon.factor = abandon.factor,diameter.tol=diameter.tol,phi=phi)
        if(export) saveRDS(df,file=write.to)
        if(verbose) cat(paste0('Particles: ',nrow(df),'. Elapsed time: ',round(proc.time()[3]-tic),' sec.\n'))
      }
    }

    if((d==2)&(flow)){
      while(nrow(df)<N){
        df=dla2dinf_flow(N=min(checkin.every,N-nrow(df)),initial=df,abandon.factor = abandon.factor,diameter.tol=diameter.tol,exit_density=exit_density)
        if(export) saveRDS(df,file=write.to)
        if(verbose) cat(paste0('Particles: ',nrow(df),'. Elapsed time: ',round(proc.time()[3]-tic),' sec.\n'))
      }
    }

    if(d==3){
      while(nrow(df)<N){
        df=dla3dinf(N=min(checkin.every,N-nrow(df)),initial=df,abandon.factor = abandon.factor,diameter.tol=diameter.tol,phi=phi,costheta=costheta)
        if(export) saveRDS(df,file=write.to)
        if(verbose) cat(paste0('Particles: ',nrow(df),'. Elapsed time: ',round(proc.time()[3]-tic),' sec.\n'))
      }
    }


  }else{
    sourceblocked=FALSE

    if(d==2){
      while(nrow(df)<N){
        df=dla2dsource(N=min(checkin.every,N-nrow(df)),initial=df, source=source, abandon.factor = abandon.factor,diameter.tol=diameter.tol,phi=phi)
        if(all(df[nrow(df),]==source)){
          cat('Cluster blocked the source.\n')
          df=df[!duplicated(df),]
          sourceblocked=TRUE
        }
        if(export) saveRDS(df,file=write.to)
        if(verbose) cat(paste0('Particles: ',nrow(df),'. Elapsed time: ',round(proc.time()[3]-tic),' sec.\n'))
        if(sourceblocked) break
      }
    }
    if(d==3){
      while(nrow(df)<N){
        df=dla3dsource(N=min(checkin.every,N-nrow(df)),initial=df, source=source, abandon.factor = abandon.factor,diameter.tol=diameter.tol,phi=phi,costheta=costheta)
        if(all(df[nrow(df),]==source)){
          cat('Cluster blocked the source.\n')
          df=df[!duplicated(df),]
          sourceblocked=TRUE
        }
        if(export) saveRDS(df,file=write.to)
        if(verbose) cat(paste0('Particles: ',nrow(df),'. Elapsed time: ',round(proc.time()[3]-tic),' sec.\n'))
        if(sourceblocked) break
      }
    }
  }

  return(df)
}

dla2dinf=function(N=200,initial=data.frame(x=0,y=0), abandon.factor=10, diameter.tol=0.01, phi=c(0,2*pi)){

  phi0=phi[1]
  phi1=phi[2]
  X=initial[,1]
  Y=initial[,2]
  R=sqrt(max(X^2+Y^2))

  tolerated=1-diameter.tol #tolerated particle diameter
  Ab=0

  n=0

  while (n<N){

      RStart=R+1
      Abandon=RStart*abandon.factor
      phi=runif(1,phi0,phi1)+pi
      # phi=tau*(0.5+rbeta(1,2,2))
      x=RStart*cos(phi)
      y=RStart*sin(phi)

    minDist=sqrt(min((X-x)^2+(Y-y)^2))
    Ab=0

    while(minDist>1){

      phi=runif(1,phi0,phi1)
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

dla2dinf_flow=function(N=200,initial=data.frame(x=0,y=0), abandon.factor=10, diameter.tol=0.01, exit_density=NULL){

  if(is.null(exit_density)){
    #Normalize density function in such a way that max(rho)<=1
    exit_density=function(x,k=0.1,r=1){
      return(exp(k*r*(cos(x)-1)))
    }
  }


  # phi0=phi[1]
  # phi1=phi[2]
  phi1=2*pi
  phi0=0
  X=initial[,1]
  Y=initial[,2]
  R=sqrt(max(X^2+Y^2))

  tolerated=1-diameter.tol #tolerated particle diameter
  Ab=0

  n=0

  while (n<N){

    RStart=R+1
    Abandon=RStart*abandon.factor
    phi=runif(phi0,phi1)
    x=RStart*cos(phi)
    y=RStart*sin(phi)

    minDist=sqrt(min((X-x)^2+(Y-y)^2))
    Ab=0

    while(minDist>1){

      Rcur=minDist-tolerated
      phi=rdens(1,exit_density, params=list(r=Rcur))
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


dla2dsource=function(N=200,initial=data.frame(x=0,y=0), source=c(100,0), abandon.factor=10, diameter.tol=0.01, phi=c(0,2*pi)){

  phi0=phi[1]
  phi1=phi[2]
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

      phi=runif(1,phi0,phi1)
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

dla3dinf=function(N=200,initial=data.frame(x=0,y=0,z=0), abandon.factor=10, diameter.tol=0.01, phi=c(0,2*pi),costheta=c(-1,1)){

  phi0=phi[1]
  phi1=phi[2]
  costheta0=costheta[1]
  costheta1=costheta[2]
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
      phi=runif(1,phi0,phi1)+pi
      cosTheta=runif(1,-costheta1,-costheta0)
      sinTheta=sqrt(1-cosTheta^2)
      x=RStart*cos(phi)*sinTheta
      y=RStart*sin(phi)*sinTheta
      z=RStart*cosTheta


    minDist=sqrt(min((X-x)^2+(Y-y)^2+(Z-z)^2))
    Ab=0

    while(minDist>1){

      Rcur=minDist-tolerated
      phi=runif(1,phi0,phi1)
      cosTheta=runif(1,costheta0,costheta1)
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

dla3dsource=function(N=200,initial=data.frame(x=0,y=0,z=0), source=c(100,0,0), abandon.factor=10, diameter.tol=0.01, phi=c(0,2*pi),costheta=c(-1,1)){

  phi0=phi[1]
  phi1=phi[2]
  costheta0=costheta[1]
  costheta1=costheta[2]
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
      phi=runif(1,phi0,phi1)
      cosTheta=runif(1,costheta0,costheta1)
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
#' Particles will be colored in the chronological order using a color scale. The function can  produce a plot and/or return the data frame of the cluster coordinates with the particle colors added.
#'
#' @param data a dataframe or matrix with 2 or 3 columns of particle coordinates (for a 2D or 3D cluster), like the one produced by \code{dla2d} or \code{dla3d}. An additional column named "color" will be tolerated and ignored.
#' @param color0 color for the earliest particles, start color of the color scale
#' @param color1 color for the latest particles, end color of the color scale
#' @param return.data logical, whether the function should return the data frame used for plotting
#' @param plot logical, whether the function should produce a plot
#' @param dupl.rm logical, whether to remove duplicate particles (particles with identical coordinates), leaving the earliest ones. Normally, there should not be duplicate particles, but one exception is when the source of particles was at a finite point and the cluster reached it, so all particles emitted after that remain at the source.
#' @return if \code{return.data} is TRUE, will return the data frame with the extra column "color" containing particle colors
#' @export
#'
vis.dla=function(data,color0='green',color1='red',return.data=FALSE,plot=TRUE, dupl.rm=TRUE){

  data=as.data.frame(data)

  if('color' %in% colnames(data)){
    data$color=NULL
  }

  if(!(ncol(data) %in% c(2,3))){
    message('Error. The number of columns in the data must be 2 or 3.\n')
    return(NULL)
  }

  d=ncol(data)

  if(d==2) {names(data)=c('x','y')}
  if(d==3) {names(data)=c('x','y','z')}

  if(dupl.rm) data=unique(data)

  data$color=colorRampPalette(c(color0,color1))(nrow(data))

  if(plot){
  title=paste("Diffusion-Limited Aggregation of",nrow(data), "Particles")
  if(d==2){
    plot(data$x,data$y,col=data$color, pch=".", type='p',asp=1, main=title, xlab='', ylab='',
       xaxt='n',yaxt='n', cex=1, bty="n")
  }

  if(d==3){
    scatterplot3d::scatterplot3d(data$x,data$y,data$z,color=data$color,pch=".", asp=1, main=title, xlab="",ylab="", zlab="", grid=FALSE,axis=TRUE, tick.marks=FALSE, cex.symbols=3)
  }
  }
  if(return.data) return(data)

}
