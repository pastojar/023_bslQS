#---------------------------------------------------------------------

# Function implementing a simple standard posterior probability density (up to
# a multiplicative normalizing constant). It is based on a likelihood 
# definition provided by the user and a multivariate normal or
# lognormal distribution at the prior.

# Loglikeli: unweighted mean of log values ( == geometric mean of the orig. values)
sysanal.logposterior.unlim.swmm.fair2.ErrMod <- function(par, model, dataCa,
                                                         prior.mean=1, prior.sd=1, prior.cor=NA, 
                                                         prior.dist="lognormal", prior.def=NA, loglikeli=sysanal.loglikeli,
                                                         ...)
{
  #    print(par)##
  logprior  <- sysanal.calcpdf_mv(z = par, dist = prior.dist, mean = prior.mean,
                                  sd = prior.sd, cor = prior.cor, distdef = prior.def)
  #    print(paste("log prior: ", format(logprior, digits=2,scientific=T)))
  
  if ( is.na(logprior) ) {return(NA)}
  
  Loglikeli <- c(); Lengths <- c()
  for (i in 1 : length(dataCa)) {
    L        <- dataCa[[i]][[1]]
    out.data <- dataCa[[i]][[2]]
    inp.file <- dataCa[[i]][[3]]
    
    Loglikeli[i] <-  loglikeli( y.obs = out.data[,2], eventID = names(dataCa)[i],
                                par = par, model = model, L = L, inp.file = inp.file, out.data = out.data,
                                #Var.Bs=Var.Bs, sd.Eps=sd.Eps, par.tr = par.tr,  par.fix = par.fix 
                                ...
                                )
    
    #Lengths[i]   <-  diff(range(sysanal.decode(L)[,2]))  # diff(c(20,10)) = -10 !
    #Lengths[i] = 1
  }
  Loglikeli <- mean( Loglikeli )
  
  #   print(paste("log likeli: ", format(loglikeli, digits=2,scientific=T)))
  return(logprior + Loglikeli)
}

#---------------------------------------------------------------------

# Input-dependent bias likelihood

sysanal.loglikeli.bias.inp.JA.ErrMod <- function( par, model, eventID, L, y.obs, Var.Bs, sd.Eps, par.tr, 
                                                  inp = rep(0,length(y.obs)) , 
                                                  par.fix=NULL, ... ) {
  if ( any (par[!is.na(par)]<0) ) {
    return(-Inf)
  } 
  
  par.comb      <- c(par,par.fix)
  
  # decode layout definition:
  L.decoded <- L
  L.encoded <- L
  if ( is.vector(L) ) {
    L.decoded <- sysanal.decode(L)
  }  else {
    L.encoded <- rownames(L)
  }
  
  # calculate results of the deterministic model:
  y.calc <- model(par = par, L = L.encoded, eventID = eventID) 

  
  #  new standard value for the input provided (5.12.14)
  #   
  if (sum(inp)>0) {
    inp <- inp[ (par.comb["Del.Max"]+1-round(par.comb["Delta"])) : (length(inp)-round(par.comb["Delta"])) ] 
  }
  
  # transform results and observations:
  if (is.na(par.tr["alpha"])) {
    y.calc.trans  <- sysanal.boxcox(y.calc, par.tr["l1"], par.tr["l2"])
    y.obs.trans   <- sysanal.boxcox(y.obs, par.tr["l1"], par.tr["l2"])
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
  } else {
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


sysanal.predict.inp.bias.L1.JA <- function(par,model,L1,y.obs,
                                           Var.Bs, inp,sd.Eps,
                                           L2,
                                           y.calc=NA, par.tr,  ...)
{
  # decode likelihood definitions:
  
  L1.decoded <- L1
  L1.encoded <- L1
  if ( is.vector(L1) )
  {
    L1.decoded <- sysanal.decode(L1)
  }
  #   else
  #   {
  #     L1.encoded <- rownames(L1)
  #   }
  n1 <- length(L1.encoded)
  
  L.decoded <- L1.decoded
  L.encoded <- L1.encoded
  n2 <- 0
  n <- n1
  
  # calculate results of deterministic model:
  y.calc    <- model(par = par, L = L1, ...)
  #   if ( is.na(y.calc[1]) )
  #   {   
  #     if ( !is.na(L2[1]) ) y.calc    <- model(par, else y.calc    <- model(par,c(L1),...)
  #   }
  #   else
  #   {
  #     if ( length(y.calc) != n )
  #     {
  #       cat("*** y.calc is not of correct length:",length(y.calc),
  #           "instead of",n,"\n")
  #       return(NA)
  #     }
  #   }
  
  y.calc.L1 <- y.calc[1:n1] 
  y.calc.L2 <- NA
  
  # transform results and observations:
  
  if(is.na(par.tr["alpha"])) 
  {
    y.calc.L1.trans  <- sysanal.boxcox(y.calc.L1,par.tr["l1"],par.tr["l2"])
    if ( !is.na(L2[1]) )  y.calc.L2        <- sysanal.boxcox(y.calc[(n1+1):length(y.calc)],par.tr["l1"],par.tr["l2"])
    y.obs.trans      <- sysanal.boxcox(y.obs,par.tr["l1"],par.tr["l2"])
  }
  else{
    y.calc.L1.trans  <- sysanal.logsinh(y.calc.L1,par.tr["alpha"],par.tr["beta"])
    if ( !is.na(L2[1]) )  y.calc.L2        <- sysanal.logsinh(y.calc[(n1+1):length(y.calc)],par.tr["alpha"],par.tr["beta"])
    y.obs.trans      <- sysanal.logsinh(y.obs,par.tr["alpha"],par.tr["beta"])
  }
  
  # calculate predictions for layout 1:   
  
  if (sum(inp)>0) {
    inp <- inp[ (par.comb["Del.Max"]+1-round(par.comb["Delta"])) : (length(inp)-round(par.comb["Delta"])) ]  
  }
  
  #   Sigma.Bf.L1       <- Var.Bf(par,L1.decoded,inp) #fast B cov
  Sigma.Bs.L1       <- Var.Bs(par,L1.decoded,inp) #slow B cov
  Sigma.Eps.L1      <- diag(sd.Eps(par,L1.decoded)^2)
  
  Sum.Sigma.L1       <- Sigma.Bs.L1 + Sigma.Eps.L1
  Sum.Sigma.L1.inv   <- solve(Sum.Sigma.L1) 
  
  #   Sum.Sigma.ind.L1   <- Sigma.Eps.L1+Sigma.Bf.L1 # SE*
  #   Sum.Sigma.ind.L1.inv <- solve(Sum.Sigma.ind.L1)
  
  Bs.mean.L1 <- as.numeric(Sigma.Bs.L1 %*% Sum.Sigma.L1.inv %*% (y.obs.trans-y.calc.L1.trans))
  Bs.var.L1  <- diag(Sigma.Bs.L1 %*% Sum.Sigma.L1.inv %*% (Sigma.Eps.L1)) #+Sigma.Bf.L1
  
  #   x11();plot(  Bs.var.L1)
  #   dim(Sigma.Bs.L1)
  #   dim(Sum.Sigma.L1.inv)
  #   dim(Sigma.Eps.L1)
  #   x11();image(Sigma.Bs.L1)
  #   x11();image(Sum.Sigma.L1.inv)
  #   Bf.mean.L1 <- as.numeric(Sigma.Bf.L1 %*% Sum.Sigma.L1.inv %*% (y.obs.trans-y.calc.L1.trans))
  #   Bf.var.L1  <- diag(Sigma.Bf.L1 %*% Sum.Sigma.L1.inv %*% (Sigma.Eps.L1)) # +Sigma.Bs.L1
  #   #   Bf.var.L1  <- diag(Sigma.Bf.L1 %*% Sum.Sigma.ind.L1.inv %*% (Sigma.Eps.L1)) # variance
  # #   Bf.mean.L1 <- as.numeric(Sigma.Bf.L1 %*% Sum.Sigma.ind.L1.inv %*% (y.obs.trans-y.calc.L1.trans-Bs.mean.L1)) #Eq. 28
  #   # Err corr
  
  names(Bs.mean.L1) <- L1.encoded 
  names(Bs.var.L1)  <- L1.encoded
  
  
  
  Y.mean.L1 <- y.calc.L1.trans + (Sigma.Bs.L1) %*% Sum.Sigma.L1.inv %*% ( y.obs.trans - y.calc.L1.trans ) 
  Y.var.L1  <- diag((Sigma.Bs.L1) %*% Sum.Sigma.L1.inv %*% (Sigma.Eps.L1)) + diag(Sigma.Eps.L1) 
  
  names(Y.mean.L1) <- L1.encoded
  names(Y.var.L1)  <- L1.encoded
  
  
  
  return(list(y.calc.L1        = y.calc.L1.trans,
              Bs.mean.L1        = Bs.mean.L1,
              Bs.var.L1         = Bs.var.L1,
              Y.mean.L1         = Y.mean.L1,
              Y.var.L1          = Y.var.L1,
              y.calc.L2         = y.calc.L2,
              fsigma_b2.nL1     = Sigma.Bs.L1[nrow(Sigma.Bs.L1),ncol(Sigma.Bs.L1)]
              #               y.calc.L2        = y.calc.L2.trans,
              #               Sum.Sigma.L1.inv = Sum.Sigma.L1.inv,
              #               Sigma.Bs.L1       = Sigma.Bs.L1,
              #               Sigma.Bs.L2L1     = Sigma.Bs.L2L1
  ))
}




# ---------------------------------------------------------------------------
# Input-dependent B and E conditional on parameters 
# ---------------------------------------------------------------------------

sysanal.predict.bias.OU.JA <- function( parsamp.L1, model, L1, y.obs, #CMP
                                        predict.bias.cond,
                                        ppt,
                                        probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),
                                        L2, y.calc=NA, 
                                        par.tr, par.fix=NULL, 
                                        inp , 
                                        ... )  {
  
  
  # decode likelihood definitions:
  
  if ( is.na(L1[1]) ) { 
    L1 <-NA
  }
  
  L1.decoded <- L1
  L1.encoded <- L1
  if ( !is.na(L1) ) {
    L1.decoded <- sysanal.decode(L1)
  } else {
    L1.encoded <- rownames(L1)
  }
  n1 <- length(L1.encoded)
  
  L.decoded <- L1.decoded
  L.encoded <- L1.encoded
  n2 <- 0
  n <- n1
  
  L2.available <- TRUE
  if ( is.vector(L2) ) { 
    if ( is.na(L2[1]) ) { 
      L2.available <- FALSE 
    }
  }   
  if ( L2.available ) { 
    L2.decoded <- L2
    L2.encoded <- L2
    if ( is.vector(L2) ) {
      L2.decoded <- sysanal.decode(L2)
    } else {
      L2.encoded <- rownames(L2)
    }
    n2 <- length(L2.encoded)
    
    if ( !is.na(L1[1]) ) { 
      L.decoded <- rbind(L1.decoded,L2.decoded) 
    } else  { 
      L.decoded <- L2.decoded
    }
    
    if ( !is.na(L1[1]) ) { 
      L.encoded <- c(L1.encoded,L2.encoded)
    } else  { 
      L.encoded <- L2.encoded
    }
    
    if ( !is.na(L1[1]) ) { 
      n <- n1 + n2
    } else  {
      n <- n2
    }
  }
  
  # initialize result arrays:
  
  y.L1.samp <- matrix(nrow=nrow(parsamp.L1),ncol=n1)
  colnames(y.L1.samp) <- L1.encoded
  Bs.L1.samp <- y.L1.samp
  Y.L1.samp <- y.L1.samp
  neg.var.Bs.L1 <- rep(0,n1)
  names(neg.var.Bs.L1) <- L1.encoded 
  neg.var.Y.L1 <- neg.var.Bs.L1
  #   neg.var.Bf.L1 <- rep(0,n1)
  #   names(neg.var.Bf.L1) <- L1.encoded
  
  y.L2.samp    <- NA
  Bs.L2.samp    <- NA
  Y.L2.samp    <- NA
  neg.var.Bs.L2 <- NA 
  neg.var.Y.L2 <- NA
  Bf.L2.samp    <- NA
  neg.var.Bf.L2 <- NA
  #   
  
  if ( !is.na(L2[1]) ) {
    y.L2.samp <- matrix(nrow=nrow(parsamp.L1),ncol=n2)
    colnames(y.L2.samp) <- L2.encoded
    Bs.L2.samp <- y.L2.samp
    Y.L2.samp <- y.L2.samp
    bias.slow.var <- y.L2.samp #!!!!
    neg.var.Bs.L2 <- rep(0,n2)
    names(neg.var.Bs.L2) <- L2.encoded
    neg.var.Y.L2 <- neg.var.Bs.L2
  }
  
  
  # calculate samples:
  
  par.old <- rep(NA,ncol(parsamp.L1))
  num.eval <- 0
  
  for ( j in 1:nrow(parsamp.L1) ) {
    par <- c(parsamp.L1[j,],par.fix)
    
    if ( !is.na(L1[1]) ) {  
      if ( j==1 | sum((par-par.old)^2, na.rm = T) != 0 )  { # do not reevaluate if parameter values stayed the same
        
        par.old <- par
        num.eval <- num.eval + 1      
        res <- predict.bias.cond(par     = par,
                                 model   = model,
                                 L1      = L1, # before: L1.decoded
                                 y.obs   = y.obs,
                                 L2      = L2, 
                                 y.calc  = y.calc,
                                 par.tr  = par.tr,
                                 inp     = inp ,
                                 Var.Bs  = sysanal.Var.Bs ,
                                 sd.Eps  = sysanal.sd.Eps.L ,
                                 ...)
      }
      y.L1.samp[j,] <- res$y.calc.L1; # remember: it's transformed
      
      if ( n2 > 0 ) y.L2.samp[j,] <- res$y.calc.L2; # remember: it's transformed  
      
      
      
    } else { # i.e. if there is just L2
      
      y.calc.L2    <- model(par, L2, ...)
      names(y.calc.L2) <- rownames(L2)
      
      
      # transform results and observations FOR MULTIOUTPUT WITH DIFFERENT TRANSF OF EACH:
      out.vars  <- unique(as.character(L.decoded$var))
      
      
      y.calc.L2.trans <- rep(NA, length(y.calc.L2))
      
      if(!is.na(par.tr[paste("l1")]))  {
        y.calc.L2.trans  <- sysanal.boxcox(y.calc.L2,par.tr[paste("l1")],par.tr["l2"])
      } else {                                                                             #JA#
        #if (!is.na(par.tr[paste("alpha",var, sep="_")])) 
        if (!is.na(par.tr[paste("alpha")])) {
          #y.calc.L2.trans [ind.var.L2]  <- sysanal.logsinh(y.calc.L2[ind.var.L2],par.tr[paste("alpha")],par.tr[paste("beta",var, sep="_")])
          y.calc.L2.trans <- sysanal.logsinh(y.calc.L2, par.tr[paste("alpha")], par.tr[paste("beta")])          
        } else  {                                                                          #JA#
          print("problem with transf param in pred")
        }
      }
      
      y.L2.samp[j,]  <- y.calc.L2.trans 
      
    }
    
    
    
    
    # Draw realization/OU paths of the errors
    if ( !is.na(L1[1]) ) {
      
      for ( i in 1:n1 ) {
        var <- res$Y.var.L1[i]
        if ( var < 0 ) { 
          cat("* Warning: negative variance of Y:",var,"at",L1.encoded[i],"\n") 
          var <- 0; neg.var.Y.L1[i] <- neg.var.Y.L1[i]+1
        }
        Y.L1.samp[j,i] <- rnorm( 1, res$Y.mean.L1[i] , sqrt(var) ) # Draw realizations for Y
      }
      
      if ( ! is.na(res$Bs.mean.L1[1]) ) {
        
        for ( i in 1:n1 ) {   
          var <- res$Bs.var.L1[i]
          if ( var < 0 ) { 
            cat("* Warning: negative variance of Bs:",var,"at",L1.encoded[i],"\n") 
            var <- 0; neg.var.Bs.L1[i] <- neg.var.Bs.L1[i]+1 
          }
          Bs.L1.samp[j,i] <- rnorm( 1, res$Bs.mean.L1[i], sqrt(var) )  # Draw realizations for Bs 
        }
        
        for ( i in 1:n1 ) {   
          var <- res$Bf.var.L1[i]
          #         if ( var < 0 ) 
          #         { 
          #           cat("* Warning: negative variance of Bf:",var,"at",L1.encoded[i],"\n") 
          #           var <- 0; neg.var.Bf.L1[i] <- neg.var.Bf.L1[i]+1 
          #         }
          #         Bf.L1.samp[j,i] <- rnorm(1,res$Bf.mean.L1[i],sqrt(var))  # Draw realizations for Bf 
        }
      }
    }
    
    
  }
  
  
  
  neg.var.Bs.L1 <- neg.var.Bs.L1/n1
  #   neg.var.Bf.L1 <- neg.var.Bf.L1/n1
  neg.var.Y.L1 <- neg.var.Y.L1/n1
  
  
  
  # AR(1) predictions in L2 
  
  if ( n2 > 0 )  {
    Dt        <- L.decoded$val[2]-L.decoded$val[1]
    vars <- unique(L.decoded$var)
    var=vars[1]
    for ( var in vars ) {
      out_count = which(vars==var) # which variable are we obs: i, ii, iii...
      
      for ( j in 1:nrow(parsamp.L1) ) {   # for all (selected) MCMC realizations (of several parameter sets) # Add: for every variable
                                          # draw realizations of the observed system output and the model
        par       <- c(parsamp.L1[j,],par.fix)
        
        beta      <- 1/par["corrlen"]^1; if ( is.na(beta) ) beta <- 0 # !!difference with original!!     
        exp.2dec  <- exp(-2*beta*(Dt))
        exp.dec   <- exp(-beta*(Dt))
        
        name.ks.B <- paste("ks",var,sep="_")
        ks.B      <- par[name.ks.B] ; if (is.na(ks.B)) ks.B <- 0
        
        name.sd.B <- paste("sd.B",var,sep="_")
        sigma_b2  <- par[name.sd.B]^2 ; if ( is.na(sigma_b2) ) sigma_b2  <- 0
        
        name.sd.Eps <- paste("sd.Eps",var,sep="_")
        sd.Eps      <- par[name.sd.Eps] ;  if ( is.na(sd.Eps) ) stop(paste("*** parameter", name.sd.Eps, "not found in sd.Eps.L"))
        
        #         if (!is.na(y.obs))    bs        = Bs.L1.samp[j,n1/length(vars)*out_count] 
        #         else                  bs        = rnorm(1,0,par[name.sd.B]) # discuss with Carlo case where we have no L1
        
        bs        <- rnorm(1,0,par[name.sd.B])
        
        
        if (sum(ppt)>0) {
          ppt_t <- ppt[ (par["Del.Max"]+1-round(par["Delta"])) : (length(ppt)-round(par["Delta"])) ]
        } else {
          ppt_t <- ppt
        }
        
        
        for ( i in 1:(n2/length(vars)) ) # for all points in time L2
        { 
          #       slow_proc_var   = sigma_b2 + (fsigma_b2 -sigma_b2-(ppt[n1+i]*ks.B)^2)*exp.2dec + 1*(ppt[n1+i]*ks.B)^2 # independent variance
          
          slow_jump_var   <- (sigma_b2+ (ppt_t[n1/length(vars)+i]*ks.B)^2)*(1-exp.2dec) # unconditional on previous variance
          
          Bs.L2.samp[j,(i+((n2/length(vars))*(out_count-1)))] <- rnorm( 1, bs*exp.dec, sqrt(slow_jump_var) ) 
          
          E.L2.samp <- rnorm(1,0,sd.Eps)
          
          Y.L2.samp[j,(i+((n2/length(vars))*(out_count-1)))]  <- y.L2.samp[j,(i+((n2/length(vars))*(out_count-1)))] +    # realizations/points in the path
                                                                 Bs.L2.samp[j,(i+((n2/length(vars))*(out_count-1)))] + 
                                                                 E.L2.samp  
          
          #       fsigma_b2       = slow_proc_var
          
          bias.slow.var[j,(i+((n2/length(vars))*(out_count-1)))] = slow_jump_var
          
          bs             = Bs.L2.samp[j,(i+((n2/length(vars))*(out_count-1)))] # newB=act_val (dep on prev val=bs, jump_var=slow_jump_var)
          
          #       print(slow_jump_var)
          #       print(bs_1)
          
        }
      } 
    }
    
  }   
  
  
  
  # derive quantiles:
  
  y.L1.quant <- matrix(nrow=length(probs),ncol=n1)
  colnames(y.L1.quant) <- L1.encoded
  rownames(y.L1.quant) <- probs
  Btot.L1.quant <- y.L1.quant
  yplusBtot.L1.quant <- y.L1.quant
  #   Bs.L1.quant <- y.L1.quant
  #   Bf.L1.quant <- y.L1.quant
  Y.L1.quant  <- y.L1.quant
  
  if ( n1 >0) { 
    for ( i in 1:n1 ) {  # for all time points in inference period and (we have distributions given by different paths)
      
      #     B.L1.quant[,i]      <- quantile(B.L1.samp[,i],probs=probs,na.rm=TRUE) # modif!!
      y.L1.quant[,i]         <- quantile(y.L1.samp[,i],probs=probs,na.rm=TRUE)
      Btot.L1.quant[,i]      <- quantile(Bs.L1.samp[,i],probs=probs,na.rm=TRUE) 
      #     Bs.L1.quant[,i]        <- quantile(Bs.L1.samp[,i],probs=probs,na.rm=TRUE)
      #     Bf.L1.quant[,i]        <- quantile(Bf.L1.samp[,i],probs=probs,na.rm=TRUE)
      yplusBtot.L1.quant[,i] <- quantile(y.L1.samp[,i]+Bs.L1.samp[,i],probs=probs,na.rm=TRUE) # modif!!
      Y.L1.quant[,i]         <- quantile(Y.L1.samp[,i],probs=probs,na.rm=TRUE)
    }
  }
  
  
  y.L2.quant  <- NA
  Btot.L2.quant <- NA
  yplusBtot.L2.quant <- NA
  Y.L2.quant  <- NA
  bias.slow.var <- NA
  
  if ( n2 > 0 ) {
    y.L2.quant <- matrix(nrow=length(probs),ncol=n2)
    colnames(y.L2.quant) <- L2
    rownames(y.L2.quant) <- probs
    #     Bs.L2.quant <- y.L2.quant
    #     Bf.L2.quant <- y.L2.quant
    Btot.L2.quant <- y.L2.quant
    yplusBtot.L2.quant <- y.L2.quant
    Y.L2.quant <- y.L2.quant
    
    for ( i in 1:n2 ) {
      y.L2.quant[,i]         <- quantile(y.L2.samp[,i],probs=probs,na.rm=TRUE)
      Btot.L2.quant[,i]      <- quantile(Bs.L2.samp[,i],probs=probs,na.rm=TRUE)
      #       Bs.L2.quant[,i]        <- quantile(Bs.L2.samp[,i],probs=probs,na.rm=TRUE)
      #       Bf.L2.quant[,i]        <- quantile(Bf.L2.samp[,i],probs=probs,na.rm=TRUE)
      yplusBtot.L2.quant[,i] <- quantile(y.L2.samp[,i]+Bs.L2.samp[,i],probs=probs,na.rm=TRUE)
      Y.L2.quant[,i]         <- quantile(Y.L2.samp[,i],probs=probs,na.rm=TRUE)
    }
  }
  
  
  
  return(list(y.L1.samp       = y.L1.samp,
              Y.L1.samp       = Y.L1.samp,
              Bs.L1.samp      = Bs.L1.samp,
              #               Bf.L1.samp      = Bf.L1.samp,
              Bs.L2.samp      = Bs.L2.samp,
              #               Bf.L2.samp      = Bf.L2.samp,
              Y.L2.samp       = Y.L2.samp,
              y.L2.samp       = y.L2.samp,
              yplusB.L1.quant = yplusBtot.L1.quant,
              yplusB.L2.quant = yplusBtot.L2.quant,
              Y.L1.quant      = Y.L1.quant,
              y.L1.quant      = y.L1.quant,
              Y.L2.quant      = Y.L2.quant,
              y.L2.quant      = y.L2.quant,
              #               Bs.L1.quant     = Bs.L1.quant,
              #               Bf.L1.quant     = Bf.L1.quant,
              #               Bs.L2.quant     = Bs.L2.quant,
              #               Bf.L2.quant     = Bf.L2.quant,
              Btot.L1.quant   = Btot.L1.quant,
              Btot.L2.quant   = Btot.L2.quant, 
              bias.slow.var.L2= bias.slow.var,
              neg.var.Bs.L1   = neg.var.Bs.L1,
              #               neg.var.Bf.L1   = neg.var.Bf.L1,
              neg.var.Y.L1    = neg.var.Y.L1,
              neg.var.Bs.L2   = neg.var.Bs.L2,
              #               neg.var.Bf.L2   = neg.var.Bf.L2,
              neg.var.Y.L2    = neg.var.Y.L2))
}




