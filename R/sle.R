
#Complex square-root, branch applicable to UHP->UHP mapping
sqrtp=function(z){
  modz=sqrt(z[1]^2+z[2]^2)
  re=sqrt((modz+z[1])*0.5)
  im=sqrt((modz-z[1])*0.5)
  if(z[2]<0) re=-re
  return(c(re,im))
}

#Infinitesimal inverse Loewner map (fourdt is 4*dt)
f=function(z,xi,fourdt){
  z=sqrtp(c((z[1]-xi)^2-z[2]^2-fourdt,2*z[2]*(z[1]-xi)))
  z[1]=z[1]+xi
  return(z)
}

#' sle
#'
#' Generate SLE trace driven by Brownian motion, possibly with the addition of Levy-flights
#'
#' @param kappa strength of the 1d Brownian motion in the driving function.
#' @param tmax max time to which the trace is developed.
#' @param a exponent of the Levy flights component in the driving function.
#' @param kappaL strength of the Levy flights component in the driving function.
#' @param nsteps number of steps.
#' @param p_timescaling exponent determining the distribution of time steps. The succession of times will be made uniform in the variable \code{t^p_timescaling}.
#' @param verbose boolean, to print progress statements or not.
#'
#' @note SLE (Stochastic Loewner Evolution) is a generative stochastic process for growing a stocastic curve (trace) out of a boundary of a 2D domain.
#' It uses a continuous family of conformal maps, parametrized by a "time" parameter: w(z,t). The SLE equation is dw(z,t)/dt = 2/(w(z,t)-xi(t)), w(z,0)=z. Here xi(t) is the real-valued driving function of the process, assumed to be a stochastic process.
#' The mapping w(z,t) is from the upper half plane. In the standard SLE, xi(t) is the 1D Brownian motion with diffusion constant kappa (and an intricate connection to conformal field theory exists).
#' As a generalization, we also allow xi(t) to be a sum of 1D Brownian motion and Levy fligts.
#'
#' For more details see this publication and references therein:
#'
#' Rushkin, I., Oikonomou, P., Kadanoff, L.P. and Gruzberg, I.A., 2006. Stochastic Loewner evolution driven by Lévy processes. Journal of Statistical Mechanics: Theory and Experiment, 2006(01), p.P01001.
#'
#' @return List with components: \code{t} - vector of time values, \code{xi} - vector of values of the driving function, \code{t_cross} - crossover time between Brownian and Levy components (NULL if \code{kappaL = 0}), \code{call_params} - list of call parameters, \code{runtime} - elapsed time in seconds.
#' @export
sle=function(kappa=4,tmax=1, a=1, kappaL=0, nsteps=2000, p_timescaling=0.5, verbose=TRUE){
  tic=proc.time()[3]

  t=(seq(0,tmax^p_timescaling,length.out = nsteps+1))^(1/p_timescaling)
  dt=diff(t)

  xi=cumsum(rnorm(nsteps)*sqrt(dt))*sqrt(kappa)
  t_cross=NULL
  if(kappaL>0){
    xi=xi+cumsum(sample(c(1,-1),nsteps,replace=TRUE)*runif(nsteps)^(-1/a)*(dt^(1/a)))*kappaL^(1/a)
    t_cross=((kappa^a)/(kappaL^2))^(1/(2-a))
  }

  fourdt=4*dt

  #Trace coordinates for times after 0
  gamma=as.data.frame(do.call(rbind,
                              lapply(1:length(xi),function(j){
                                temp=c(xi[j],0)
                                for(i in j:1){
                                  temp=f(temp,xi[i],fourdt=fourdt[i])
                                }

                                if(verbose){
                                  if(j %% 1000 == 0){
                                    runt=round(proc.time()[3]-tic)
                                    cat(paste0(j,' steps; ',runt %/% 60,' min. ',runt %% 60,' sec.\n'))
                                  }
                                }

                                return(temp)
                              })
  )
  )

  #Add at the begnning the initial location of the trace
  xi=c(0,xi)
  gamma=rbind(c(0,0),gamma)
  names(gamma)=c('x','y')
  toc=proc.time()[3]

  return(list(t=t, xi=xi, gamma=gamma, t_cross=t_cross, call_params=list(
    kappa=kappa,tmax=tmax, a=a, kappaL=kappaL, nsteps=nsteps, p_timescaling=0.5
  )
  ,runtime=toc-tic)
  )

}
