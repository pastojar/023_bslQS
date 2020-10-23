#---------------------------------------------------------------------

# Function implementing a simple standard posterior probability density (up to
# a multiplicative normalizing constant). It is based on a likelihood 
# definition provided by the user and a multivariate normal or
# lognormal distribution at the prior.

sysanal.logposterior.unlim.swmm <- function(par, model, dataCa,
                                            prior.dist="lognormal", prior.mean=1, prior.sd=1,
                                            prior.cor=NA, prior.def=NA,
                                            loglikeli=sysanal.loglikeli,
                                            ...)
{
  #    print(par)##
  logprior  <- sysanal.calcpdf_mv(z=par,dist=prior.dist,mean=prior.mean,
                                  sd=prior.sd,cor=prior.cor,distdef=prior.def)
  #    print(paste("log prior: ", format(logprior, digits=2,scientific=T)))
  
  if ( is.na(logprior) ) {return(NA)}
  
  Loglikeli <- 0
  for (i in 1 : length(dataCa)) {
    L        <- dataCa[[i]][[1]]
    out.data <- dataCa[[i]][[2]]
    inp.file <- dataCa[[i]][[3]]
    Loglikeli <-  Loglikeli + loglikeli(par=par, model=model, L=L, y.obs=out.data[,2], 
                                        #Var.Bs=Var.Bs, sd.Eps=sd.Eps, par.tr = par.tr,  par.fix = par.fix, 
                                        inp.file=inp.file, out.data=out.data, ...)
  }
  
  
  #   print(paste("log likeli: ", format(loglikeli, digits=2,scientific=T)))
  return(logprior + Loglikeli)
}

#---------------------------------------------------------------------

# Loglikeli: unweighted mean of log values ( == geometric mean of the orig. values)
sysanal.logposterior.unlim.swmm.fair2 <- function(par, model, dataCa,
                                                  prior.dist="lognormal", prior.mean=1, prior.sd=1,
                                                  prior.cor=NA, prior.def=NA,
                                                  loglikeli=sysanal.loglikeli,
                                                  ...)
{
  #    print(par)##
  logprior  <- sysanal.calcpdf_mv(z=par,dist=prior.dist,mean=prior.mean,
                                  sd=prior.sd,cor=prior.cor,distdef=prior.def)
  #    print(paste("log prior: ", format(logprior, digits=2,scientific=T)))
  
  if ( is.na(logprior) ) {return(NA)}
  
  Loglikeli <- c(); Lengths <- c()
  for (i in 1 : length(dataCa)) {
    L        <- dataCa[[i]][[1]]
    out.data <- dataCa[[i]][[2]]
    inp.file <- dataCa[[i]][[3]]
    Loglikeli[i] <-  loglikeli(par=par, model=model, L=L, y.obs=out.data[,2], 
                               #Var.Bs=Var.Bs, sd.Eps=sd.Eps, par.tr = par.tr,  par.fix = par.fix, 
                               inp.file=inp.file, out.data=out.data, ...)
    
    #Lengths[i]   <-  diff(range(sysanal.decode(L)[,2]))  # diff(c(20,10)) = -10 !
    #Lengths[i] = 1
  }
  Loglikeli <- mean( Loglikeli )
  
  #   print(paste("log likeli: ", format(loglikeli, digits=2,scientific=T)))
  return(logprior + Loglikeli)
}

#---------------------------------------------------------------------

# Input-dependent bias likelihood

sysanal.loglikeli.bias.inp.JA <- function(par, model, L, y.obs, Var.Bs, inp=rep(0,length(y.obs)), sd.Eps,
                                          par.tr, par.fix=NULL, ...)
{
  if (any (par[!is.na(par)]<0))
  {return(-Inf)}###
  
  
  # decode layout definition: 
  
  L.decoded <- L
  L.encoded <- L
  if ( is.vector(L) )
  {
    L.decoded <- sysanal.decode(L)
  }
  else
  {
    L.encoded <- rownames(L)
  }
  
  # calculate results of deterministic model:
  
  #y.calc        <- model(par, L.encoded, inp.file=inp.file, out.data=out.data ) #changed!! L.decoded
  y.calc        <- model(par, L.encoded, ... ) #changed!! L.decoded
  
  #JA# removes NAs
  if (length(which(is.na(y.obs))) > 0) {
    y.calc    <- y.calc   [- which(is.na(y.obs))]
    L.decoded <- L.decoded[- which(is.na(y.obs)), ]
    L.encoded <- L.encoded[- which(is.na(y.obs))]
    y.obs     <- y.obs    [- which(is.na(y.obs))]
  }
  
  # evaluate likelihood function:
  par.comb      <- c(par,par.fix)
  
  #  new standard value for the input provided (5.12.14)
  #   
  if (sum(inp)>0) {inp = inp[(par.comb["Del.Max"]+1-round(par.comb["Delta"])):(length(inp)-round(par.comb["Delta"]))] }
  
  # transform results and observations:
  
  if (is.na(par.tr["alpha"])) {
    y.calc.trans  <- sysanal.boxcox(y.calc,par.tr["l1"],par.tr["l2"])
    y.obs.trans   <- sysanal.boxcox(y.obs,par.tr["l1"],par.tr["l2"])
    boxcox.deriv  <- sysanal.boxcox.deriv(y.obs,par.tr["l1"],par.tr["l2"])
  }  else {
    y.calc.trans  <- sysanal.logsinh(y.calc,par.tr["alpha"],par.tr["beta"])
    y.obs.trans   <- sysanal.logsinh(y.obs,par.tr["alpha"],par.tr["beta"])
    boxcox.deriv  <- sysanal.logsinh.deriv(y.obs,par.tr["alpha"],par.tr["beta"])
  }
  
  #   Sigma.Bf      <- Var.Bf(par.comb,L.decoded, inp) #fast B cov
  Sigma.Bs      <- Var.Bs(psi=par.comb, L=L.decoded, inp=inp) #slow B cov
  
  Sigma.Eps     <- diag(sd.Eps(par.comb,L.decoded)^2)
  
  Sum.Sigma     <-  Sigma.Bs +Sigma.Eps #2 new covariance matrices
  
  Sum.Sigma.inv <- solve(Sum.Sigma)
  
  
  log.det.Sum.Sigma <- determinant(Sum.Sigma,logarithm=TRUE)
  
  if ( log.det.Sum.Sigma$sign < 0 ) { 
    warning("determinant Sigma.Eps+Sigma.B < 0") 
    loglikeli=-Inf
  }
  else {
    loglikeli <- ( - 0.5 * length(L.encoded) * log(2*pi) -
                     0.5 * as.numeric(log.det.Sum.Sigma$modulus) -
                     0.5 * t(y.obs.trans-y.calc.trans) %*% (Sum.Sigma.inv ) %*%  (y.obs.trans-y.calc.trans) +
                     sum(log(abs(boxcox.deriv)), na.rm = FALSE)
    )
  }
  
  #   print(loglikeli)
  return(loglikeli)
}

#---------------------------------------------------------------------

# function to plot prior and posterior marginals:

sysanal.plot.margs.JA <- function(postsamp,pridist=list(),ncol=NA,mar=NA,
                                  xlim=list(),ymax=list(),
                                  titles=list(),xlab=list(),ylab=list(),
                                  lty=NA)
{
  # transform samples to a list of all samples to plot the marginals of:
  
  postsamp.list <- list()
  if ( is.data.frame(postsamp) )
  {
    postsamp.list[[1]] <- postsamp
  }
  else
  {
    if ( is.matrix(postsamp) )
    {
      postsamp.list[[1]] <- postsamp
    }
    else
    {
      if ( is.list(postsamp) )  # be careful: a data frame is a list!
      {
        postsamp.list <- postsamp
      }
      else
      {
        stop("sysanal.plot.margs: postsamp is of illegal type")
      }
    }
  }
  nsamp <- length(postsamp.list)
  
  # transform prior definitions to a list of priors to plot:
  
  pridist.list <- list()
  if ( length(pridist) > 0 )
  {
    if ( length(names(pridist)) != length(pridist) )
    {
      pridist.list <- pridist
    }
    else
    {
      pridist.list[[1]] <- pridist
    }
  }
  npri <- length(pridist.list)
  
  # get all variable names of the sample: 
  
  var <- colnames(postsamp.list[[1]])
  if ( nsamp > 1 )
  {
    for ( j in 2:nsamp ) 
    {
      var <- c(var,colnames(postsamp.list[[j]]))
    }
    var <- unique(var)
  }
  nvar <- length(var)
  
  # define layout of plot panels:
  
  if ( !is.na(ncol) ) nc <- ncol
  else                nc <- floor(sqrt(nvar))
  nr <- ceiling(nvar/nc)
  marg <- mar
  if ( is.na(marg[1]) ) marg <- c(2.0,2.0,2.5,0.5) # c(bottom, left, top, right)
  
  # define line types:
  
  lty.loc <- lty
  if ( is.na(lty.loc[1]) ) lty.loc <- 1
  if ( length(lty.loc) < nsamp+npri ) lty.loc <- 1:(nsamp+npri)
  
  # plot marginals of samples and priors:
  
  par.def <- par(no.readonly=T)
  par(mfrow=c(nr,nc),xaxs="i",yaxs="i",mar=marg) # c(bottom, left, top, right)
  for ( i in 1:nvar )
  {
    name <- var[i]
    density.max <- NA
    d <- as.list(rep(NA,nsamp))
    for ( j in 1:nsamp )
    {
      ind <- match(var[i],colnames(postsamp.list[[j]]))
      if ( !is.na(ind) )
      {
        bw <- .5*sd(postsamp.list[[j]][,ind])
        d[[j]]      <- density(postsamp.list[[j]][,ind],bw=bw) ### changed!###
        #             d[[j]]      <- density(postsamp.list[[j]][,ind],adjust=1.0)
        density.max <- max(c(density.max,d[[j]]$y),na.rm=TRUE)
      }
    }
    
    xlim.i <- c(min(d[[j]]$x), max(d[[j]]$x))
    if ( xlim.i[1] == xlim.i[2] ) 
    {
      xlim.i[1] <- 0.5*xlim.i[1]
      xlim.i[2] <- 1.5*xlim.i[2]
    }
    if ( length(xlim) > 0 )
    {
      ind <- match(name,names(xlim))
      if ( !is.na(ind) )
      {
        xlim.i <- xlim[[ind]]
      }
    }
    
    if ( is.na(density.max) ) ylim.i <- c(0,1)
    else                      ylim.i <- c(0,1.1*density.max)
    if ( length(ymax) > 0 )
    {
      ind <- match(name,names(ymax))
      if ( !is.na(ind) )
      {
        ylim.i <- c(0,ymax[[ind]])
      }
    }
    
    title.i <- name
    if ( length(titles) > 0 )
    {
      ind <- match(name,names(titles))
      {
      if ( !is.na(ind) )
      {
        title.i <- titles[[ind]]
      }
      }
    }
    
    xlab.i <- name
    if ( length(xlab) > 0 )
    {
      ind <- match(name,names(xlab))
      if ( !is.na(ind) )
      {
        xlab.i <- xlab[[ind]]
      }
    }
    
    ylab.i <- "f"
    if ( length(ylab) > 0 )
    {
      ind <- match(name,names(ylab))
      if ( !is.na(ind) )
      {
        ylab.i <- ylab[[ind]]
      }
    }
    
    plot(numeric(0),numeric(0),main=title.i,xlab=xlab.i,ylab=ylab.i,
         xlim=xlim.i,ylim=ylim.i,type="n")
    
    for ( j in 1:nsamp )
    {
      if ( !is.na(d[[j]])[[1]] )
      {
        #             lines(d[[j]],lty=lty.loc[j])
        polygon(x=c(d[[j]]$x, rev(d[[j]]$x)),
                y=c(d[[j]]$y, rev(rep(0, length(d[[j]]$y)))),### changed!###
                col = "lightgrey", border=NA)
      }
    }
    
    if ( npri > 0 )
    {
      x <- seq(xlim.i[1],xlim.i[2],length=101)  # WHY 101? :D
      y <- numeric(0)
      for ( j in 1:npri )
      {
        if ( length(pridist.list[[j]]) > 0 )
        {
          ind <- sysanal.multmatch(name,names(pridist.list[[j]]))
          if ( !is.na(ind[1]) )
          {
            for ( k in 1:length(ind) )
            {
              y <- sysanal.calcpdf(x=x, distpar=pridist.list[[j]][[ind[k]]])
              lines(x,y,lty=lty.loc[nsamp+j])
            }
            
            if (pridist.list[[1]][[ind[1]]][1] == "NormalTrunc")
            {
              q1.pr <- truncnorm::qtruncnorm(0.05, a = as.numeric(pridist.list[[j]][[ind[k]]][4]), b = as.numeric(pridist.list[[j]][[ind[k]]][5]), 
                                             mean = as.numeric(pridist.list[[j]][[ind[k]]][2]), sd = as.numeric(pridist.list[[j]][[ind[k]]][3]))
              
              q2.pr <- truncnorm::qtruncnorm(0.95, a = as.numeric(pridist.list[[j]][[ind[k]]][4]), b = as.numeric(pridist.list[[j]][[ind[k]]][5]), 
                                             mean = as.numeric(pridist.list[[j]][[ind[k]]][2]), sd = as.numeric(pridist.list[[j]][[ind[k]]][3]))
              
              
              q1.po <- quantile(postsamp.list[[1]][,i], probs = 0.05)
              q2.po <- quantile(postsamp.list[[1]][,i], probs = 0.95)
            }
            #legend("topright", legend = paste("pr VS po: ", round(1-(q2.po-q1.po)/(q2.pr-q1.pr), digits=3)*100, "%", sep=""))
            legend("topright", legend = paste("pr / po = ", round((q2.po-q1.po)/(q2.pr-q1.pr), digits=3), sep=""))
          }
        }
      }
    }
  }
  par(par.def)
}

#---------------------------------------------------------------------

