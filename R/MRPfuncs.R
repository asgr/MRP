.getcos=function(ref){
  cosref = NULL
  data('cosref',envir = environment())
  allownames=tolower(as.character(cosref[,'Ref']))
  if(tolower(ref) %in% allownames==FALSE){stop(paste('Provided ref name is not allowed, must be one of',paste(as.character(cosref[,'Ref']),sep='',collapse=', '),' (case insensitive). See ?cosref for details.'))}
  out=as.numeric(cosref[allownames==tolower(ref),])
  names(out)=colnames(cosref)
  return(out)
}

MRP_dNdM=function(H=10, norm=2.151349-19, Hs=14.4469690, alpha=-1.8643512, beta=0.7268376, parm){
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

MRP_dNdlnM=function(H=10, norm=2.151349-19, Hs=14.4469690, alpha=-1.8643512, beta=0.7268376, parm){
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

MRP_dNdlog10M=function(H=10, norm=2.151349-19, Hs=14.4469690, alpha=-1.8643512, beta=0.7268376, parm){
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

MRP_PDF=function(H=10, Hs=14.4469690, alpha=-1.8643512, beta=0.7268376, Hmin=8, parm){
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

MRPint_N=function(Hmin=8, norm=2.151349-19, Hs=14.4469690, alpha=-1.8643512, beta=0.7268376, parm){
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

MRPint_M=function(Hmin=8, norm=2.151349-19, Hs=14.4469690, alpha=-1.8643512, beta=0.7268376, parm){
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

MRPnorm=function(Hs=14.4469690, alpha= -1.8643512, beta=0.7268376, OmegaM = 0.315, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm){
	#Utility function for MRP HMF work
	#Function to normalise a given differential dn/dm HMF to rhomean0
  H0=100
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    OmegaL=as.numeric(params['OmegaL'])
    if(!is.na(params['OmegaR'])){OmegaR=as.numeric(params['OmegaR'])}
  }
  if(!missing(parm)){
    if(length(parm)==3){
      Hs=parm[1]
      alpha=parm[2]
      beta=parm[3]
    }
  }
	rhomean0=cosgrowRhoMean(z=0, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, Dist='Co', Mass='Msun')
	temp=MRPint_M(Hmin=-Inf,norm=1,Hs=Hs,alpha=alpha,beta=beta)
	return(rhomean0/temp)
}

MRPrho_dNdM=function(H=10, Hs=14.4469690, alpha=-1.8643512, beta=0.7268376, OmegaM = 0.315, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm){
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
	norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  x=10^(H-Hs)
  return(norm*beta*(x^(alpha)*exp(-x^beta)))
}

MRPrho_dNdlnM=function(H=10, Hs=14.4469690, alpha=-1.8643512, beta=0.7268376, OmegaM = 0.315, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm){
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
	norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  x=10^(H-Hs)
  return(norm*beta*(x^(alpha+1)*exp(-x^beta)))
}

MRPrho_dNdlog10M=function(H=10, Hs=14.4469690, alpha=-1.8643512, beta=0.7268376, OmegaM = 0.315, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm){
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
	norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref, parm)
  x=10^(H-Hs)
  return(norm*beta*(x^(alpha+1)*exp(-x^beta))*log(10))
}

MRPrhoint_N=function(Hmin=8, Hs=14.4469690, alpha= -1.8643512, beta=0.7268376, OmegaM = 0.315, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm){
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

MRPrhoint_M=function(Hmin=8, Hs=14.4469690, alpha= -1.8643512, beta=0.7268376, OmegaM = 0.315, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref, parm){
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
	norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
	return(norm*(10^Hs)^2*gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+2)/beta))
}

MRPmode_dMdlnM=function(Hs=14.4469690, alpha=-1.8643512, beta=0.7268376, parm){
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

MRPmode_dMlog10M=function(Hs=14.4469690, alpha=-1.8643512, beta=0.7268376, parm){
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

MRPmean=function(Hmin=8, Hs=14.4469690, alpha= -1.8643512, beta=0.7268376, parm){
	#Utility function for MRP HMF work
	#Mean mass of halos Msun/h
	return(log10((10^Hs)*gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+2)/beta)/gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+1)/beta)))
}

MRPvariance=function(Hmin=8, Hs=14.4469690, alpha= -1.8643512, beta=0.7268376, parm){
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

MRPstdev=function(Hmin=8, Hs=14.4469690, alpha= -1.8643512, beta=0.7268376, parm){
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

#HMFbins=-diff(MRPintNnorm(10^seq(8,15.5,by=0.5),H0=H0))


###Evo stuff:

# MRP.W13.z0.Hs=function(OmegaM=0.315, Sigma8=0.829, mu=14.476, sigma=0.060345, a=0.63008, b=38.3311, c=18.390){
#   f = a + b*log(OmegaM+Sigma8) - c*Sigma8
#   return(mu+f*sigma)
# }
#
# MRP.W13.z0.alpha=function(OmegaM=0.315, Sigma8=0.829, mu=-1.8084, sigma=0.0059875, a=44.793, b=18.526, c=29.173){
#   f = a*OmegaM + b*Sigma8 - c
#   return(mu+f*sigma)
# }
#
# MRP.W13.z0.beta=function(OmegaM=0.315, Sigma8=0.829, mu=0.76583, sigma=0.014030, a=19.0412, b=14.922, c=13.7617){
#   f = a + b*log(Sigma8) + c*log(OmegaM)
#   return(mu+f*sigma)
# }
#
# MRP.W13.z0=function(OmegaM=0.315, Sigma8=0.829, Hs.mu=14.476, Hs.sigma=0.060345, Hs.a=0.63008, Hs.b=38.3311, Hs.c=18.390, alpha.mu=-1.8084, alpha.sigma=0.0059875, alpha.a=44.793, alpha.b=18.526, alpha.c=29.173, beta.mu=0.76583, beta.sigma=0.014030, beta.a=19.0412, beta.b=14.922, beta.c=13.7617){
#   Hs=MRP.W13.z0.Hs(OmegaM=OmegaM, Sigma8=Sigma8, mu=Hs.mu, sigma=Hs.sigma, a=Hs.a, b=Hs.b, c=Hs.c)
#   alpha=MRP.W13.z0.alpha(OmegaM=OmegaM, Sigma8=Sigma8, mu=alpha.mu, sigma=alpha.sigma, a=alpha.a, b=alpha.b, c=alpha.c)
#   beta=MRP.W13.z0.beta(OmegaM=OmegaM, Sigma8=Sigma8, mu=beta.mu, sigma=beta.sigma, a=beta.a, b=beta.b, c=beta.c)
#   return(c(Hs=Hs, alpha=alpha, beta=beta))
# }
#
# MRP.W13.0z6.Hs=function(OmegaM=0.315, Sigma8=0.829, z=1, mu=12.816, sigma=1.1903, a=1.6487, b=0.34575, c=1.62733, d=0.84700){
#   f = a*Sigma8 + b*OmegaM*Sigma8*z^c - d*z
#   return(mu+f*sigma)
# }
#
# MRP.W13.0z6.alpha=function(OmegaM=0.315, Sigma8=0.829, z=1, mu=-1.8981, sigma=0.02605, a=0.37917, b=75.253, c=0.38414, d=24.102, q=0.35229){
#   f = Sigma8 - a*z/Sigma8 + b*OmegaM*c^z - d*q^z
#   return(mu+f*sigma)
# }
#
# MRP.W13.0z6.beta=function(OmegaM=0.315, Sigma8=0.829, z=1, mu=0.60627, sigma=0.10037, a = 5.8319, b = 0.90429){
#   f = a*OmegaM*Sigma8 - z*b^z
#   return(mu+f*sigma)
# }
#
# MRP.W13.0z6=function(z=1, OmegaM=0.315, Sigma8=0.829, Hs.mu=12.816, Hs.sigma=1.1903, Hs.a = 1.6487, Hs.b = 0.34575, Hs.c = 1.62733, Hs.d = 0.84700, alpha.mu=-1.8981, alpha.sigma=0.02605, alpha.a = 0.37917, alpha.b = 75.253, alpha.c = 0.38414, alpha.d = 24.102, alpha.q = 0.35229, beta.mu=0.60627, beta.sigma=0.10037, beta.a = 5.8319, beta.b = 0.90429){
#   Hs=MRP.W13.0z6.Hs(OmegaM=OmegaM, Sigma8=Sigma8, z=z, mu=Hs.mu, sigma=Hs.sigma, a=Hs.a, b=Hs.b, c=Hs.c, d=Hs.d)
#   alpha=MRP.W13.0z6.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, mu=alpha.mu, sigma=alpha.sigma, a=alpha.a, b=alpha.b, c=alpha.c, d=alpha.d, q=alpha.q)
#   beta=MRP.W13.0z6.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, mu=beta.mu, sigma=beta.sigma, a=beta.a, b=beta.b)
#   return(c(Hs=Hs, alpha=alpha, beta=beta))
# }
#
# MRP.W13.6z10.Hs=function(z=1, OmegaM=0.315, Sigma8=0.829, mu=9.3361, sigma=0.56585, a=8.7068, b=8.5622, c=3.3042, d=0.82902){
#   f = a*OmegaM + b*Sigma8 - c - d*z
#   return(mu+f*sigma)
# }
#
# MRP.W13.6z10.alpha=function(z=1, OmegaM=0.315, Sigma8=0.829, mu=-2.1262, sigma=0.0078966, a=4.6376, b=0.69937, c=7.1539, d=13.565){
#   f = a + b*z - c*Sigma8 - d*OmegaM
#   return(mu+f*sigma)
# }
#
# MRP.W13.6z10.beta=function(z=1, OmegaM=0.315, Sigma8= 0.829, mu=0.38326, sigma=0.025932, a=20.762, b=0.78510){
#   f = Sigma8 + a*OmegaM*Sigma8 - b*z
#   return(mu+f*sigma)
# }
#
# MRP.W13.6z10=function(z=1, OmegaM=0.315, Sigma8=0.829, Hs.mu=9.3361, Hs.sigma=0.56585, Hs.a=8.7068, Hs.b=8.5622, Hs.c=3.3042, Hs.d=0.82902, alpha.mu=-2.1262, alpha.sigma=0.0078966, alpha.a=4.6376, alpha.b=0.69937, alpha.c=7.1539, alpha.d=13.565, beta.mu=0.38326, beta.sigma=0.025932, beta.a=20.762, beta.b=0.78510){
#   Hs=MRP.W13.6z10.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, mu=Hs.mu, sigma=Hs.sigma, a=Hs.a, b=Hs.b, c=Hs.c, d=Hs.d)
#   alpha=MRP.W13.6z10.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, mu=alpha.mu, sigma=alpha.sigma, a=alpha.a, b=alpha.b, c=alpha.c, d=alpha.d)
#   beta=MRP.W13.6z10.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, mu=beta.mu, sigma=beta.sigma, a=beta.a, b=beta.b)
#   return(c(Hs=Hs, alpha=alpha, beta=beta))
# }
#
# MRP.W13=function(z=1, OmegaM=0.315, Sigma8=0.829){
#   if(z==0){output=MRP.W13.z0(OmegaM=OmegaM, Sigma8=Sigma8)}
#   if(z>0 & z<=6){output=MRP.W13.0z6(z=z, OmegaM=OmegaM, Sigma8=Sigma8)}
#   if(z>6 & z<=10){output=MRP.W13.6z10(z=z, OmegaM=OmegaM, Sigma8=Sigma8)}
#   if(z>10){output=MRP.W13.6z10(z=z, OmegaM=OmegaM, Sigma8=Sigma8)}
#   return(output)
# }

MRP.B13.Hs=function(z=0, OmegaM=0.315, Sigma8=0.829, mu=12.2158, sigma=1.64125, a=1.69554, b=5.40377, c=0.884497, d=5.76511, ref){
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  f = OmegaM + a*Sigma8 + b*(c^z) - d
  return(mu + sigma*f)
}

MRP.B13.alpha=function(z=0, OmegaM=0.315, Sigma8=0.829, mu=-1.91000, sigma=0.0268194, a=3.04605, b=0.175997, c=-1.88607, d=1.52347, ref){
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  f = a*OmegaM + Sigma8^2 + b*OmegaM*(z^(c*z)) - d*log(1+z)
  return(mu + sigma*f)
}

MRP.B13.beta=function(z=0, OmegaM=0.315, Sigma8=0.829, mu=0.500559, sigma=0.128926, a=6.27010, b=2.01529, c=0.531006, d=1.89763, e=0.567776, ref){
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  f = a*Sigma8*OmegaM + b*c^z - d - e*OmegaM*z
  return(mu + sigma*f)
}

MRP.B13=function(z=0, OmegaM=0.315, Sigma8=0.829, Hs.mu=12.2158, Hs.sigma=1.64125, Hs.a=1.69554, Hs.b=5.40377, Hs.c=0.884497, Hs.d=5.76511, alpha.mu=-1.91000, alpha.sigma=0.0268194, alpha.a=3.04605, alpha.b=0.175997, alpha.c=-1.88607, alpha.d=1.52347, beta.mu=0.500559, beta.sigma=0.128926, beta.a=6.27010, beta.b=2.01529, beta.c=0.531006, beta.d=1.89763, beta.e=0.567776, ref){
  if(!missing(ref)){
    params=.getcos(ref)
    OmegaM=as.numeric(params['OmegaM'])
    if(!is.na(params['Sigma8'])){Sigma8=as.numeric(params['Sigma8'])}
  }
  Hs=MRP.B13.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, mu=Hs.mu, sigma=Hs.sigma, a=Hs.a, b=Hs.b, c=Hs.c, d=Hs.d)
  alpha=MRP.B13.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, mu=alpha.mu, sigma=alpha.sigma, a=alpha.a, b=alpha.b, c=alpha.c, d=alpha.d)
  beta=MRP.B13.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, mu=beta.mu, sigma=beta.sigma, a=beta.a, b=beta.b, c=beta.c, d=beta.d, e=beta.e)
  return(c(Hs=Hs, alpha=alpha, beta=beta))
}

MRP.B13rho_dNdM=function(z=0, H=seq(10,15,by=0.01), OmegaM = 0.315, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, Sigma8=0.829, ref, masses=TRUE){
	#Utility function for MRP HMF work
	#Differential dN/dM HMF forced to integrate to the rhomean0
  Hs=MRP.B13.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  alpha=MRP.B13.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  beta=MRP.B13.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
	norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  x=10^(H-Hs)
  y=norm*beta*(x^(alpha)*exp(-x^beta))
  if(masses){out=cbind(H, y)}else(out=y)
  return(out)
}

MRP.B13rho_dNdlnM=function(z=0, H=seq(10,15,by=0.01), OmegaM = 0.315, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, Sigma8=0.829, ref, masses=TRUE){
	#Utility function for MRP HMF work
	#Differential dN/dM HMF forced to integrate to the rhomean0
  Hs=MRP.B13.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  alpha=MRP.B13.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  beta=MRP.B13.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
	norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  x=10^(H-Hs)
  y=norm*beta*(x^(alpha+1)*exp(-x^beta))
  if(masses){out=cbind(H, y)}else(out=y)
  return(out)
}

MRP.B13rho_dNdlog10M=function(z=0, H=seq(10,15,by=0.01), OmegaM = 0.315, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, Sigma8=0.829, ref, masses=TRUE){
	#Utility function for MRP HMF work
	#Differential dN/dM
  #HMF forced to integrate to the rhomean0
  Hs=MRP.B13.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  alpha=MRP.B13.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  beta=MRP.B13.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
	norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
  x=10^(H-Hs)
  y=norm*beta*(x^(alpha+1)*exp(-x^beta))*log(10)
  if(masses){out=cbind(H, y)}else(out=y)
  return(out)
}

MRP.B13rhoint_N=function(Hmin=8, z=0, OmegaM = 0.315, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, Sigma8=0.829, ref){
	#Utility function for MRP HMF work
	#Number of halos above Hmin for arbitrary nromalisation
	#HMF forced to integrate to rhomean0
  Hs=MRP.B13.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  alpha=MRP.B13.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  beta=MRP.B13.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
	norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
	return(norm*(10^Hs)*gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+1)/beta))
}

MRP.B13rhoint_M=function(Hmin=8, z=0, OmegaM = 0.315, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, Sigma8=0.829, ref){
	#Utility function for MRP HMF work
	#Mass of halos above Hmin for arbitrary normalisation
	#HMF forced to integrate to rhomean0
  Hs=MRP.B13.Hs(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  alpha=MRP.B13.alpha(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
  beta=MRP.B13.beta(z=z, OmegaM=OmegaM, Sigma8=Sigma8, ref=ref)
	norm=MRPnorm(Hs=Hs, alpha=alpha, beta=beta, OmegaM=OmegaM, OmegaL=OmegaL, OmegaR=OmegaR, w0=w0, wprime=wprime, ref=ref)
	return(norm*(10^Hs)^2*gamma_inc((10^(Hmin-Hs))^beta,a=(alpha+2)/beta))
}

