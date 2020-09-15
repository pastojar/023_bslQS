#########################################
# J. Pastorek, APR 2016
# based on a script by O. Wani
#
# MCMC Bayesian calibration, prediction and uncertainty analysis
#########################################


###################################################################################################################################

Ca <- function(prodata, par, input, transf, seed, runs) {
  
 
  ##############################################################################
  ## Defines initial  error model parameters
   dataCa <- prodata$Ca

  ref.out   <- c(y.ref_Q = 200)
  par.init  <- c(par,              # initial (SWMM) model parameters defined elsewhere
                 sd.Eps_Q = 0.5,   
                 sd.B_Q   = 50,    # based on anaylsis of previous data
                 corrlen  = 0.5
                ) 
  pri_sd    <- 0.5*par.init 
  par.fix   <- c(Del.Max = 0, 
                 ks_Q    = 0, 
                 Delta   = 0
                 )
  
  par.tr <- transf$par.tr
  
  lim.sx.ks_Q = NA
  
  
  ##############################################################################
  ## Tranforms error model hyperparameters
  
  var="Q"
  
  (if (is.na(par.tr["alpha"])) 
  {  # BC or no transformation
    
    if (!is.na(pri_sd["corrlen"])) {
      par.init[paste("sd.B",var, sep="_")] <- par.init[paste("sd.B",var, sep="_")]*
                                              sysanal.boxcox.deriv(ref.out[paste("y.ref",var, sep="_")],
                                                                   par.tr["l1"], par.tr["l2"])
    }  
    
    par.init[paste("sd.Eps",var, sep="_")] <- par.init[paste("sd.Eps",var, sep="_")]*
                                              sysanal.boxcox.deriv(ref.out[paste("y.ref",var, sep="_")],
                                                                   par.tr["l1"], par.tr["l2"])
    
    if (!is.na(par.init[paste("ks",var, sep="_")])) { # we have input dependence via an heteroskedastik EM
      lim.sx.ks_Q = as.numeric(sysanal.boxcox(0,par.tr["l1"],par.tr["l2"]) / max(imber.vel))
      pri_sd["ks_Q"] = (-lim.sx.ks_Q+(sysanal.boxcox(par.init["ks_Q"],par.tr["l1"],par.tr["l2"]) / max(imber.vel)))
      par.init["ks_Q"] = 0.1*pri_sd["ks_Q"]  # start from a small value
    } 
  }  
  
  else  
  { # log-sinh
    
    if (!is.na(pri_sd["corrlen"])) {
      par.init[paste("sd.B",var, sep="_")] <- par.init[paste("sd.B",var, sep="_")]*
                                              sysanal.logsinh.deriv(ref.out[paste("y.ref",var, sep="_")],
                                                                    par.tr["alpha"],par.tr["beta"])
    }
    
    par.init[paste("sd.Eps",var, sep="_")] <- par.init[paste("sd.Eps",var, sep="_")]*
                                              sysanal.logsinh.deriv(ref.out[paste("y.ref",var, sep="_")], 
                                                                    par.tr["alpha"],par.tr["beta"])
    
    if (!is.na(par.init[paste("ks",var, sep="_")])) {
      lim.sx.ks_Q = as.numeric(sysanal.logsinh(0,par.tr["alpha"],par.tr["beta"]) / max(imber.vel))
      pri_sd["ks_Q"] = (-lim.sx.ks_Q+(sysanal.logsinh(par.init["ks_Q"],par.tr["alpha"],par.tr["beta"]) / max(imber.vel)))
      par.init["ks_Q"] = 0.1*pri_sd["ks_Q"]
    } 
  } 
  )   
  
  ##############################################################################
  ## Defines prior parameter distribution
  prior.pbdis<- list(#mult.imp     =c("NormalTrunc", par.init["mult.imp"], par.init["mult.imp"], 0.8, 1.2 ),
                     #mult.wid     =c("NormalTrunc", par.init["mult.wid"], par.init["mult.wid"], 0.3, 1.7 ),
                     mult.slo     =c("NormalTrunc", par.init["mult.slo"], par.init["mult.slo"], 0.3, 1.7 ),
                     mult.Nim     =c("NormalTrunc", par.init["mult.Nim"], par.init["mult.Nim"], 0.3, 1.7 ),
                     mult.Sim     =c("NormalTrunc", par.init["mult.Sim"], par.init["mult.Sim"], 0.3, 1.7 ),
                     #mult.Spe     =c("NormalTrunc", par.init["mult.Spe"], par.init["mult.Spe"], 0.3, 1.7 ),
                     #mult.Pze     =c("NormalTrunc", par.init["mult.Pze"], par.init["mult.Pze"], 0.3, 1.7 ),
                     #mult.rou     =c("NormalTrunc", par.init["mult.rou"], par.init["mult.rou"], 0.3, 1.7 ),
                     corrlen      =c("NormalTrunc", par.init["corrlen"], pri_sd["corrlen"], 0.01, 3 ), # [hr]; same units as layout
                     sd.Eps_Q     =c("NormalTrunc", par.init["sd.Eps_Q"], pri_sd["sd.Eps_Q"], 0.01, 1.5 ),#param of the observ err model
                     sd.B_Q       =c("NormalTrunc", par.init["sd.B_Q"], pri_sd["sd.B_Q"],0, 1e+08 ),
                     ks_Q         =c("NormalTrunc", lim.sx.ks_Q, pri_sd["ks_Q"], lim.sx.ks_Q, 1e+08 ),
                     Delta        =c("Exponential", par.init["Delta"])
  )
  # print(prior.pbdis)
  
  
  ##############################################################################
  ## Defines objective function
  
  logposterior.unlim.swmm <- function(par)
  { 
    names(par) <- names(par.init)
    out <- sysanal.logposterior.unlim.swmm.fair2( par,
                                                  model         = model.swmm,
                                                  dataCa        = dataCa, 
                                                  prior.dist    = "indep",
                                                  prior.def     = prior.pbdis,
                                                  loglikeli     = sysanal.loglikeli.bias.inp.JA,
                                                  par.fix       = par.fix,
                                                  par.tr        = par.tr,
                                                  Var.Bs        = sysanal.Var.Bs,
                                                  sd.Eps        = sysanal.sd.Eps.L
                                                 )
    
    if (rnorm(1, mean = 0, sd = 1) > 1.9) { # to monitor the progress during iterations
      print(paste("log post: ", format(out, digits=2)))
      print(par)
    }
    return(out)
  }
  
  
  # Checks ability to compute posterior
  
  logposterior.unlim.swmm(par = par.init)
  
  
  # Prepares calibration
  
  nlogposterior <- function(par) {
    out = -1*logposterior.unlim.swmm(par)
    return(out)
  }
  
  
  ############################################################################## 
  ## RUNS CALIBRATION in 2 steps
  
  set.seed(seed)
  
  # i. Optimization
  
  low_ran = rep(0.1,length(par.init)); names(low_ran)=names(par.init)
  up_ran  = par.init*10

  Opt.precal <- GenSA::GenSA( par = par.init,
                              fn    = nlogposterior, 
                              lower =  low_ran,
                              upper =  up_ran,
                              control= list( max.call = runs[1], verbose=T)
                            )

  par.optim.1  <- Opt.precal$par;  names(par.optim.1) <- names(par.init)

  
  # ii. Improves jump distribution and keeps optimizing
  
  RAM      <- adaptMCMC::MCMC( p     = logposterior.unlim.swmm,
                               init  = par.optim.1*rnorm(n=length(par.init),mean=1,sd=0.01),
                               scale = diag((pri_sd/2)^2,length(par.init)),
                               n     = runs[2] + runs[3], 
                               adapt = runs[2],    
                               acc.rate = 0.3, # maybe bit smaller (0.27)           use gamma ?
                               n.start = 100 # maybe larger
                             )
  
  
  # Time taken so far
  end_time=proc.time()
  time_taken=end_time-start_time
  time_taken
  
  return(list(pr.dis=prior.pbdis, Opt.precal=Opt.precal, par.optim.1=par.optim.1, RAM=RAM,
              par.tr=par.tr, par.fix=par.fix))
}




###################################################################################################################################




Pre <- function(prodata, transf, runs, RAM, par.tr, par.fix) {
  
  dataCa  <- prodata$Ca
  dataPre <- prodata$Pre
  
  ##############################################################################
  ## Predictions...   (with uncertainty propagation)

  if (runs[4] < runs[3] ) {
    MCMC.propa <- RAM$samples[ sample( seq((nrow(RAM$samples) - runs[3] + 1), nrow(RAM$samples), by=1), runs[4] ),  ]
  } 
  if (runs[4] == runs[3]) {
    MCMC.propa <- RAM$samples[seq((nrow(RAM$samples) - runs[3] + 1), nrow(RAM$samples), 1), ]
  }
  
  # Calibration Phase
  
  res.swmm.LCa <- list()
  for (i in 1 : length(dataCa)) {
    res.swmm.LCa[[i]] <- CaPre.predict.Ca(evdata = dataCa[[i]], MCMC.propa = MCMC.propa, par.tr = par.tr, par.fix = par.fix)
  }
  
  
  # Validation Phase (predicts for Pre events)
  
  res.swmm.LPre <- list()
  for (i in 1 : length(dataPre)) {
    res.swmm.LPre[[i]] <- CaPre.predict.Pre(evdata = dataPre[[i]], MCMC.propa = MCMC.propa, par.tr = par.tr, par.fix = par.fix)
  }
  
  
  
  
  #########################################################
  ## Back Transform

  #  Ca events
  
  bTr.Ca <- list()
  for (i in 1 : length(dataCa)) {
    bTr.Ca[[i]] <- CaPre.bTr.Ca(transf = transf, res.swmm.LCa = res.swmm.LCa[[i]], L.Ca = dataCa[[i]][[1]])
  }
  
  
  #  Pre events
  
  bTr.Pre <- list()
  for (i in 1 : length(dataPre)) {
    bTr.Pre[[i]] <- CaPre.bTr.Pre(transf = transf, res.swmm.LPre = res.swmm.LPre[[i]], L.Pre = dataPre[[i]][[1]])
  }

  
 
  
  return(list(MCMC.propa=MCMC.propa, bTr.Ca=bTr.Ca, bTr.Pre=bTr.Pre))

}

###################################################################################################################################

