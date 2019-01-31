

harm=function(curve, k=1){

  #Make sure the curve is closed:
  z=c(z,z[1])

  z=0.5*(curve[-1]+curve[-length(curve)])

  return(sum((z^(-k))*Conj(z)*diff(curve))*complex(real=0, imaginary=-1)/(2*pi))
}






t=2*pi*seq(0,1,1e-3)

z=complex(real=cos(t), imaginary=sin(t))


harm(z,k=0)
