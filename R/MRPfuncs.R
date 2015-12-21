##### Formal Truncated Generalised Gamma Distribution functions for R #####

dtggd = function(x, logh=14, a=-1, b=1, xmin=10, log=FALSE){
  d = b*10^(a*(x-logh))*exp(-10^((x-logh)*b))/(10^logh * gamma_inc((a+1)/b,10^(b*(xmin-logh))))
  if(log){d=log(d)}
  return(d)
}

ptggd = function(q, logh=14, a=-1, b=1, xmin=10, lower.tail=TRUE, log.p=FALSE){
  p = gamma_inc((a+1)/b,10^(b*(q-logh)))/gamma_inc((a+1)/b,10^(b*(xmin-logh)))
  if(lower.tail){p=1-p}
  p[p>1]=1
  p[p<0]=0
  if(log.p){p=log(p)}
  return(p)
}

qtggd = function(p, logh=14, a=-1, b=1, xmin=10, lower.tail=TRUE, log.p=FALSE, res.approx=1e-4){
  mmax = logh + 2.5/b
  logm = seq(xmin, mmax, res.approx)
  cdf = ptggd(q=logm, logh=logh, a=a, b=b, xmin=xmin, lower.tail=lower.tail)
  icdf = approxfun(cdf, logm)
  p[p>1]=1
  p[p<0]=0
  if(log.p){p=exp(p)}
  return(icdf(p))
}

rtggd = function(n, logh=14, a=-1, b=1, xmin=10, res.approx=1e-4){
  return(qtggd(runif(n), logh=logh, a=a, b=b, xmin=xmin, res.approx=res.approx))
}

## densities in logx,
## this is the pdf of log10(m/Hs)
# dtggd_log = function(y,ymin, a,b){
#   z = (a+1)/b
#   x = 10^(y*b)
#   xmin = 10^(ymin*b)
#   return = b * log(10) * 10^(y*(a+1)) * exp(-x) / gamma_inc(z,xmin)
# }
#
# qtggd_log = function(y,ymin,a,b){
#   z = (a+1)/b
#   return = gamma_inc(z,10^(y*b))/gamma_inc(z,10^(ymin*b))
# }
#
# rtggd_log = function(n,ymin,a,b){
#   ymax = 2.5/b
#
#   y = seq(ymin, ymax, 0.0001)
#   cdf = qtggd_log(y,ymin,a,b)
#   icdf = approxfun(cdf, y)
#
#   return = icdf(runif(n))
# }

#########

.getcos=function(ref){
  cosref = NULL
  data('cosref',envir = environment())
  allownames=tolower(as.character(cosref[,'Ref']))
  if(tolower(ref) %in% allownames==FALSE){stop(paste('Provided ref name is not allowed, must be one of',paste(as.character(cosref[,'Ref']),sep='',collapse=', '),' (case insensitive). See ?cosref for details.'))}
  out=as.numeric(cosref[allownames==tolower(ref),])
  names(out)=colnames(cosref)
  return(out)
}

MRP_dNdM=function(H=10, norm=2.258917e-19, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, parm){
	#Utility function for MRP HMF work
	#Differential dN/dM HMF with arbitrary normalisation
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
  x=10^(H-Hs)
  return(norm*beta*(x^(alpha)*exp(-x^beta)))
}

MRP_dNdlnM=function(H=10, norm=2.258917e-19, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, parm){
	#Utility function for MRP HMF work
	#Differential dN/dlnM HMF with arbitrary normalisation
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
  x=10^(H-Hs)
  return(norm*beta*(x^(alpha+1)*exp(-x^beta)))
}

MRP_dNdlog10M=function(H=10, norm=2.258917e-19, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, parm){
	#Utility function for MRP HMF work
	#Differential dN/dlog10M HMF with arbitrary normalisation
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
  x=10^(H-Hs)
  return(norm*beta*(x^(alpha+1)*exp(-x^beta))*log(10))
}

MRP_PDF=function(H=10, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, Hmin=8, parm){
	#Utility function for MRP HMF work
	#HMF forced to form a true PDF down to Hmin
  #This correctly integrates to 1, i.e. this returns the density for a true PDF.
    if(!missing(parm)){
      if(length(parm)==3){
        Hs=parm[1]
        alpha=parm[2]
        beta=parm[3]
      }
    }
  	x=10^(H-Hs)
  	return(beta*(x^(alpha)*exp(-x^beta)/gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+1)/beta)/(10^Hs)))
}

MRPint_N=function(Hmin=8, norm=2.258917e-19, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, parm){
	#Utility function for MRP HMF work
	#Number of halos above Hmin for arbitrary normalisation
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	return(norm*(10^Hs)*gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+1)/beta))
}

MRPint_M=function(Hmin=8, norm=2.258917e-19, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, parm){
	#Utility function for MRP HMF work
	#Mass of halos above Hmin for arbitrary normalisation
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	return(norm*(10^Hs)^2*gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+2)/beta))
}

MRPnorm=function(Hs=14.47256, alpha= -1.863804 , beta=0.7196221, OmegaM = 0.308, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm){
	#Utility function for MRP HMF work
	#Function to normalise a given differential dn/dm HMF to rhomean0
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){Sigma8=as.numeric(params['OmegaR'])}
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	rhomean0=cosgrowRhoMean(z=0, H0=100, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, Dist='Co', Mass='Msun')
	temp=MRPint_M(Hmin=-Inf,norm=1,Hs=Hs,alpha=alpha,beta=beta)
	return(rhomean0/temp)
}

MRPrho_dNdM=function(H=10, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, OmegaM = 0.308, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm, normtype='cos'){
	#Utility function for MRP HMF work
	#Differential dN/dM
  #HMF forced to integrate to the rhomean0
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	if(normtype=='cos'){
	  norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  }
  if(normtype=='A'){
    norm=A
  }
  x=10^(H-Hs)
  return(norm*beta*(x^(alpha)*exp(-x^beta)))
}

MRPrho_dNdlnM=function(H=10, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, OmegaM = 0.308, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm, normtype='cos'){
	#Utility function for MRP HMF work
	#Differential dN/dlnM
  #HMF forced to integrate to the rhomean0
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	if(normtype=='cos'){
	  norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  }
  if(normtype=='A'){
    norm=A
  }
  x=10^(H-Hs)
  return(norm*beta*(x^(alpha+1)*exp(-x^beta)))
}

MRPrho_dNdlog10M=function(H=10, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, OmegaM = 0.308, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm, normtype='cos'){
	#Utility function for MRP HMF work
	#Differential dN/dlog10M
  #HMF forced to integrate to the rhomean0
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	if(normtype=='cos'){
	  norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  }
  if(normtype=='A'){
    norm=A
  }
  x=10^(H-Hs)
  return(norm*beta*(x^(alpha+1)*exp(-x^beta))*log(10))
}

MRPrhoint_N=function(Hmin=8, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, OmegaM = 0.308, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm){
	#Utility function for MRP HMF work
	#Number of halos above Hmin for arbitrary nromalisation
	#HMF forced to integrate to rhomean0
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
	return(norm*(10^Hs)*gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+1)/beta))
}

MRPrhoint_M=function(Hmin=8, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, OmegaM = 0.308, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm, normtype='cos'){
	#Utility function for MRP HMF work
	#Mass of halos above Hmin for arbitrary normalisation
	#HMF forced to integrate to rhomean0
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	if(normtype=='cos'){
	  norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  }
  if(normtype=='A'){
    norm=A
  }
	return(norm*(10^Hs)^2*gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+2)/beta))
}

MRPmode_dMdlnM=function(Hs=14.47256, alpha=-1.863804 , beta=0.7196221, parm){
	#Utility function for MRP HMF work
	#HMF mass mode for dM/dlnM
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
  return(log10((10^Hs)*((alpha+2)/beta)^(1/beta)))
}

MRPmode_dMlog10M=function(Hs=14.47256, alpha=-1.863804 , beta=0.7196221, parm){
	#Utility function for MRP HMF work
	#HMF mass mode for dM/dlog10M
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
  return(log10((10^Hs)*((alpha+2)/beta)^(1/beta)))
}

MRPmean=function(Hmin=8, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, parm){
	#Utility function for MRP HMF work
	#Mean mass of halos Msun/h
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	return(log10((10^Hs)*gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+2)/beta)/gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+1)/beta)))
}

MRPvariance=function(Hmin=8, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, parm){
	#Utility function for MRP HMF work
	#Variance in mass of halos Msun/h
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	return(log10(((10^Hs)^2)*(gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+3)/beta)/gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+1)/beta)-(gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+2)/beta)/gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+1)/beta))^2)))
}

MRPstdev=function(Hmin=8, Hs=14.47256, alpha=-1.863804 , beta=0.7196221, parm){
	#Utility function for MRP HMF work
	#Variance in mass of halos Msun/h
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	return(0.5*log10(((10^Hs)^2)*(gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+3)/beta)/gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+1)/beta)-(gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+2)/beta)/gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+1)/beta))^2)))
}

###Evo stuff:

MRP.B13.Hs=function(z=0, OmegaM=0.308, Sigma8=0.815, mu=12.2158, sigma=1.64125, a=0.058562, b=1.4394, c=0.39111, d=0.11159, e=0.056010, f=0.42444, g=0.90369, h=0.0029417, ref){
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  f = a + b*Sigma8 + c*OmegaM + d*Sigma8*z + e*z^2 + f*Sigma8*OmegaM*z - g*z - h*z^3
  return(mu + sigma*f)
}

MRP.B13.alpha=function(z=0, OmegaM=0.308, Sigma8=0.815, mu=-1.91000, sigma=0.0268194, a=2.6172, b=2.06023, c=1.4791, d=2.2142, e=0.53400, f=2.70981, g=0.19690, ref){
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  f = a*OmegaM + b*Sigma8 + c*d^OmegaM*e^z - f - g*z
  return(mu + sigma*f)
}

MRP.B13.beta=function(z=0, OmegaM=0.308, Sigma8=0.815, mu=0.49961, sigma=0.12913, a=7.5217, b=0.18866, c=0.36891, d=0.071716, e=0.0029092, f=3.4453, g=0.71052, ref){
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  f = a*Sigma8*OmegaM - b - c*z - d*e^z - f*OmegaM*z*g^z
  return(mu + sigma*f)
}

MRP.B13.A=function(z=0, OmegaM=0.308, Sigma8=0.815, mu=-33.268, sigma=7.3593, a=0.0029187, b=0.15541, c=1.4657, d=0.055025, e=0.24068, f=0.33620, ref){
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  f = z + a*z^3 - b - c*Sigma8 - d*z^2 - e*Sigma8*z - f*OmegaM*z
  return(exp(mu + sigma*f))
}

MRP.B13=function(z=0, OmegaM=0.308, Sigma8=0.815, Hs=list(), alpha=list(), beta=list(), A=list(), ref){
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  Hs=do.call('MRP.B13.Hs', c(z=z, OmegaM=OmegaM, Sigma8=Sigma8, Hs))
  alpha=do.call('MRP.B13.alpha', c(z=z, OmegaM=OmegaM, Sigma8=Sigma8, alpha))
  beta=do.call('MRP.B13.beta', c(z=z, OmegaM=OmegaM, Sigma8=Sigma8, beta))
  A=do.call('MRP.B13.A', c(z=z, OmegaM=OmegaM, Sigma8=Sigma8, A))
  return(c(Hs=Hs, alpha=alpha, beta=beta, A=A))
}

MRP.B13rho_dNdM=function(z=0, H=seq(10,15,by=0.01), OmegaM=0.308, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Sigma8=0.815, ref, masses=TRUE, normtype='cos'){
	#Utility function for MRP HMF work
	#Differential dN/dM
  #HMF forced to integrate to rhomean0 if normtype='cos'
  Hs=MRP.B13.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  alpha=MRP.B13.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  beta=MRP.B13.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  if(normtype=='cos'){
	  norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  }
  if(normtype=='A'){
    norm=MRP.B13.A(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  }
  x=10^(H-Hs)
  y=norm*beta*(x^(alpha)*exp(-x^beta))
  if(masses){out=cbind(H, y)}else(out=y)
  return(out)
}

MRP.B13rho_dNdlnM=function(z=0, H=seq(10,15,by=0.01), OmegaM=0.308, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Sigma8=0.815, ref, masses=TRUE, normtype='cos'){
	#Utility function for MRP HMF work
	#Differential dN/dlnM
  #HMF forced to integrate to rhomean0 if normtype='cos'
  Hs=MRP.B13.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  alpha=MRP.B13.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  beta=MRP.B13.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  if(normtype=='cos'){
	  norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  }
  if(normtype=='A'){
    norm=MRP.B13.A(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  }
  x=10^(H-Hs)
  y=norm*beta*(x^(alpha+1)*exp(-x^beta))
  if(masses){out=cbind(H, y)}else(out=y)
  return(out)
}

MRP.B13rho_dNdlog10M=function(z=0, H=seq(10,15,by=0.01), OmegaM=0.308, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Sigma8=0.815, ref, masses=TRUE, normtype='cos'){
	#Utility function for MRP HMF work
	#Differential dN/dlog10M
	#HMF forced to integrate to rhomean0 if normtype='cos'
  Hs=MRP.B13.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  alpha=MRP.B13.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  beta=MRP.B13.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  if(normtype=='cos'){
	  norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  }
  if(normtype=='A'){
    norm=MRP.B13.A(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  }
  x=10^(H-Hs)
  y=norm*beta*(x^(alpha+1)*exp(-x^beta))*log(10)
  if(masses){out=cbind(H, y)}else(out=y)
  return(out)
}

MRP.B13rhoint_N=function(Hmin=8, z=0, OmegaM=0.308, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Sigma8=0.815, ref, normtype='cos'){
	#Utility function for MRP HMF work
	#Number of halos above Hmin
	#HMF forced to integrate to rhomean0 if normtype='cos'
  Hs=MRP.B13.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  alpha=MRP.B13.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  beta=MRP.B13.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  if(normtype=='cos'){
	  norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  }
  if(normtype=='A'){
    norm=MRP.B13.A(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  }
	return(norm*(10^Hs)*gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+1)/beta))
}

MRP.B13rhoint_M=function(Hmin=8, z=0, OmegaM=0.308, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Sigma8=0.815, ref, normtype='cos'){
	#Utility function for MRP HMF work
	#Mass of halos above Hmin
	#HMF forced to integrate to rhomean0 if normtype='cos'
  Hs=MRP.B13.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  alpha=MRP.B13.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  beta=MRP.B13.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  if(normtype=='cos'){
	  norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  }
  if(normtype=='A'){
    norm=MRP.B13.A(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  }
	return(norm*(10^Hs)^2*gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+2)/beta))
}

MRP.B13survey=function(area=179.936, zmin=0, zmax=0.5, lomass=8, himass=15.5, massbin=0.5, inunit='deg2', res=100, OmegaM=0.308, OmegaL=1-OmegaM-OmegaR, OmegaR=0, w0=-1, wprime=0, Sigma8=0.815, ref){
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){Sigma8=as.numeric(params['OmegaR'])}
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  tempHMFbins={}
  zbin=(zmax-zmin)/res
  zvals=seq(zmin,zmax-zbin,length=res)
  for(z in zvals){
    tempvol=cosvol(area=area, zmax=z+zbin, zmin=z, inunit=inunit, H0=100, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)
    tempHMFbins=rbind(tempHMFbins, -diff(MRP.B13rhoint_N(Hmin=seq(lomass,himass,by=massbin), z=tempvol['volmedz']), OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime)*tempvol[1]*1e9)
  }
  out=colSums(tempHMFbins)
  out=cbind(seq(lomass,himass-massbin,by=massbin)+massbin/2,out, cumsum=rev(cumsum(rev(out))))
  return(out)
}

