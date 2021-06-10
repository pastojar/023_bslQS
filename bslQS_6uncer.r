
#######################################
## loading the local data
load( paste0( getwd(), "/outputs/bsl.QS_1cal.Rdata" ) )
load( paste0( getwd(), "/outputs/bsl.QS_2rr.Rdata" ) )
devtools::load_all(".")



#######################################
## Defines observed data and their layout for the calibration and prediction
Urquell <- system.file("swmm", "inpfile.inp", package = package) # path to the swmm catchment model

dataCa <- setupSWMMX( eventIDs = eventIDsCa, flow.data.proc = match_with_IDs(rainfall_datfr = flow.data.proc, IDs = eventIDsCa), 
                      Urquell = Urquell, package = package )
for ( i_dataCa in 1:length(dataCa) ) {
  dataCa[[i_dataCa]]$Q_Data[,2] <- imputeTS::na_interpolation( x = dataCa[[i_dataCa]]$Q_Data[,2], option = "linear", maxgap = 4  )
  dataCa[[i_dataCa]]$Q_Data[,3] <- imputeTS::na_interpolation( x = dataCa[[i_dataCa]]$Q_Data[,3], option = "linear", maxgap = 4  )
}

dataPre <- setupSWMMX( eventIDs = eventIDsPre, flow.data.proc = match_with_IDs(rainfall_datfr = flow.data.proc, IDs = eventIDsPre), 
                       Urquell = Urquell, package = package )
for ( i_dataPre in 1:length(dataPre) ) {
  dataPre[[i_dataPre]]$Q_Data[,2] <- imputeTS::na_interpolation( x = dataPre[[i_dataPre]]$Q_Data[,2], option = "linear", maxgap = 5  )
  dataPre[[i_dataPre]]$Q_Data[,3] <- imputeTS::na_interpolation( x = dataPre[[i_dataPre]]$Q_Data[,3], option = "linear", maxgap = 5  )
}


#######################################
## creates the desired rainfall data
CML_bsl <- sup.rain.data( scens = "read CML14_bslQuantSm", periods = periods )

scens <- as.character(c())
for ( i_link in rownames(par_opt_all) ) {
  
  pars_str  <- paste( par_names[1], par_opt_all[i_link,][1], par_names[2], par_opt_all[i_link,][2], sep = "-")
  proc_meth <- paste0(WAA_meth, "-", pars_str, "--AttSpec--AtoR")
  
  scens <- c( scens, paste0( "CML_bsl__single-", i_link, "--", proc_meth ) )
}
CML_Rain <- sup.rain.data( scens = scens )

subset_choice <- list()
subset_choice["arb-3-8-12-15"] <- paste( c(1, 6, 9, 12), collapse = "_" )    # CMLs # 3, 8, 12, 15

myRain <- sup.rain.data( scens = paste0("CML_Rain__subcols_mean-which_cols-", subset_choice["arb-3-8-12-15"]) ,
                         periods = periods )


#######################################
## faking the r-r model
myRuno <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                               data_new  = myRain,
                               package   = package )
colnames(myRuno)[ !colnames(myRuno) %in% c("time", "timestamp", "id", "Qobs", "sd_Qobs") ] <- "Qmod"

pseudoHydro_build <- function( myRuno ) {
  
  pseudoHydro <- function(  par, L, eventID   )  { 

    out <- match_with_IDs(rainfall_datfr = myRuno, IDs = eventID)
    hlp <- out$Qmod
    names(hlp) <- paste0( "Q_",  format(out$timestamp, nsmall = 6 ) )
    out <- hlp
    
    return(out) 
  }
  
  return( pseudoHydro )
}

pseudoHydro <- pseudoHydro_build( myRuno = myRuno )


#-------
save.image( file = paste0(getwd(), "/outputs/", package, "_6uncer.Rdata") )
#-------
load( paste0( getwd(), "/outputs/bsl.QS_1cal.Rdata" ) )
load( paste0( getwd(), "/outputs/bsl.QS_2rr.Rdata" ) )
load( paste0( getwd(), "/outputs/bsl.QS_6uncer.Rdata" ) )
devtools::load_all(".")
#-------


############################################################################## 
## PREPARES THE CALIBRATION

#######################################
## Defines the parameters and their initial values
par.init.untr  <- c( sd.Eps_Q = 2,      # units of Q data, i.e. [l/s]
                     sd.B_Q   = 0.001,  # units of Q data, i.e. [l/s]
                     corrlen  = 0.5 )   # [h]; same units as layout

# par.fix   <- c( Del.Max = 0, 
#                 ks_Q    = 0, 
#                 Delta   = 0 )
par.fix <- NA
lim.sx.ks_Q <- NA

#######################################
## Defines prior parameter distributions
pri.sd.untr    <- c( sd.Eps_Q = 2,              # units of Q data, i.e. [l/s]
                     sd.B_Q   = 25,             # units of Q data, i.e. [l/s]
                     corrlen  = 0.25  )         # [h]; same units as layout

prior.pbdis <- list(
  sd.Eps_Q     = c( "NormalTrunc", par.init.untr["sd.Eps_Q"], pri.sd.untr["sd.Eps_Q"], 0.01, 100 ),   # units of Q data, i.e. [l/s]   
  sd.B_Q       = c( "NormalTrunc", par.init.untr["sd.B_Q"],   pri.sd.untr["sd.B_Q"],   0,    1000 ),  # units of Q data, i.e. [l/s]
  corrlen      = c( "NormalTrunc", par.init.untr["corrlen"],  pri.sd.untr["corrlen"],  0.01, 3 )      # [h]; same units as layout
  # ks_Q         = c( "NormalTrunc", lim.sx.ks_Q,          pri.sd["ks_Q"], lim.sx.ks_Q, 1e+08 ),
  # Delta        = c( "Exponential", par.init.untr["Delta"] )
)
# print(prior.pbdis)



#######################################
## Sets up the output (i.e. Q) transformation method
transf  <- list( transf = "BC", # Box-Cox transformation
                 par.tr = c(l1 = 0.45, l2 = 1) )  # two different parameters of the Box-Cox transformation
par.tr  <- transf$par.tr
ref.out <- c(y.ref_Q = 300)     # units of Q data, i.e. [l/s]


#######################################
## Transforms the prior distributions of the error model parameters 
var <- "Q"

if ( is.na(par.tr["alpha"]) ) {  # BC or no transformation
  
  if ( !is.null(prior.pbdis[["corrlen"]]) ) {
    prior.pbdis[[paste0("sd.B_", var)]][-1] <- as.numeric(prior.pbdis[[paste0("sd.B_", var)]][-1]) *
                                               sysanal.boxcox.deriv( ref.out[paste0("y.ref_", var)],
                                                                     par.tr["l1"], par.tr["l2"] )
  }  
  
  prior.pbdis[[paste0("sd.Eps_", var)]][-1] <- as.numeric(prior.pbdis[[paste0("sd.Eps_", var)]][-1]) *
                                             sysanal.boxcox.deriv( ref.out[paste0("y.ref_", var)],
                                                                   par.tr["l1"], par.tr["l2"] )

} else {   # log-sinh
  
  if ( !is.null(prior.pbdis[["corrlen"]]) ) {
    prior.pbdis[[paste0("sd.B_", var)]][-1] <- as.numeric(prior.pbdis[[paste0("sd.B_", var)]][-1]) *
                                               sysanal.boxcox.deriv( ref.out[paste0("y.ref_", var)],
                                                                     par.tr["alpha"], par.tr["beta"] )
  }  
  
  prior.pbdis[[paste0("sd.Eps_", var)]][-1] <- as.numeric(prior.pbdis[[paste0("sd.Eps_", var)]][-1]) *
                                               sysanal.boxcox.deriv( ref.out[paste0("y.ref_", var)],
                                                                     par.tr["alpha"], par.tr["beta"] )
} 

par.init <- as.numeric( unlist(prior.pbdis)[c(F,T,F,F,F)] )
pri.sd   <- as.numeric( unlist(prior.pbdis)[c(F,F,T,F,F)] )
low_ran  <- as.numeric( unlist(prior.pbdis)[c(F,F,F,T,F)] )
up_ran   <- as.numeric( unlist(prior.pbdis)[c(F,F,F,F,T)] )

names(par.init) <- names(pri.sd) <- names(low_ran) <- names(up_ran) <- names( prior.pbdis )



#######################################
## Defines objective function
logposterior.unlim.swmm <- function(par) {
  
  names(par) <- names(par.init)
  out <- sysanal.logposterior.unlim.swmm.fair2.ErrMod( par        = par,
                                                       model      = pseudoHydro,
                                                       dataCa     = dataCa, 
                                                       prior.dist = "indep",
                                                       prior.def  = prior.pbdis,
                                                       loglikeli  = sysanal.loglikeli.bias.inp.JA.ErrMod,
                                                       Var.Bs     = sysanal.Var.Bs,
                                                       sd.Eps     = sysanal.sd.Eps.L,
                                                       par.fix    = par.fix,
                                                       par.tr     = par.tr
  )
  
  if (rnorm(1, mean = 0, sd = 1) > 1.9) { # to monitor the progress during iterations
    print(paste("log post: ", format(out, digits=2)))
    print(par)
  }
  return(out)
}

#######################################
# Checks the ability to compute posterior
logposterior.unlim.swmm( par = par.init ) 
logposterior.unlim.swmm( par = low_ran )
logposterior.unlim.swmm( par = up_ran )


#######################################
# Prepares the objective function to be minimized
neg.logposterior <- function(par) {
  out <- -1 * logposterior.unlim.swmm(par)
  return(out)
}


############################################################################## 
## RUNS CALIBRATION in 2 steps

runs <- 10*c(5, 8, 2, 2)
seed <- 42
set.seed(seed)

# i. Optimization
Opt.precal <- GenSA::GenSA( par     = par.init,
                            fn      = neg.logposterior, 
                            lower   = low_ran,
                            upper   = up_ran,
                            control = list( max.call = runs[1], verbose=T)
)

par.optim.1  <- Opt.precal$par;  names(par.optim.1) <- names(par.init)


# ii. Improves jump distribution and keeps optimizing
RAM      <- adaptMCMC::MCMC( p        = logposterior.unlim.swmm,
                             init     = par.optim.1 * rnorm(n = length(par.init), mean = 1, sd=0.01),
                             scale    = diag( (pri.sd/2)^2, length(par.init) ),
                             n        = runs[2] + runs[3], 
                             adapt    = runs[2],    
                             acc.rate = 0.3, # maybe bit smaller (0.27)           use gamma ?
                             n.start  = 100  # maybe larger
)


#######################################           
## plots calibration chains, marginal priors and posteriors
out_dir <- file.path( getwd(), "outputs", "uncer" )
if ( dir.exists( out_dir) == F ) {
  dir.create(out_dir)  
}

check <- FALSE
check <- plot.chains.margs( RAM = RAM, pr.dis = prior.pbdis, pack.dir = out_dir, 
                            which_samples = 1:runs[2], plot_name = "2_Adapt" )
if (check==FALSE) { dev.off() } # closes graphic device if plotting fails

check <- FALSE
check <- plot.chains.margs( RAM = RAM, pr.dis = prior.pbdis, pack.dir = out_dir, 
                            which_samples = (runs[2]+1) : (runs[2]+runs[3]), plot_name = "3_Non-Adapt" )
if (check==FALSE) { dev.off() } # closes graphic device if plotting fails


#-------
save.image( file = paste0(getwd(), "/outputs/", package, "_6uncer.Rdata") )
#-------
load( paste0( getwd(), "/outputs/bsl.QS_1cal.Rdata" ) )
load( paste0( getwd(), "/outputs/bsl.QS_2rr.Rdata" ) )
load( paste0( getwd(), "/outputs/bsl.QS_6uncer.Rdata" ) )
devtools::load_all(".")
#-------



##############################################################################
## Predictions with uncertainty propagation

#######################################
# Reuses parameter samples from the last phase of the calibration
if ( runs[4] <= runs[3]  ) {
  set.seed(seed)
  par_samples_Pre_ord <- sample( x = (nrow(RAM$samples)-runs[4]+1) : nrow(RAM$samples), size = runs[4], replace = F )
} else { stop("Wrong numbers of iterations")  }
par_samples_Pre <- RAM$samples[ par_samples_Pre_ord , ]


#######################################           
## plots prediction chains, marginal priors and posteriors
check <- FALSE
check <- plot.chains.margs( RAM = RAM, pr.dis = prior.pbdis, pack.dir = out_dir, 
                            which_samples = par_samples_Pre_ord , plot_name = "4_Predict" )
if (check==FALSE) { dev.off() } # closes graphic device if plotting fails


#######################################
# Prediction and back-transformation for the calibration phase (Ca events) and validation phase (Pre events)
res.LCa <- list(); res.LCa.bTr <- list()
for (i in 1 : length(dataCa)) {
  res.LCa[[i]] <- CaPre.predict( L1 = dataCa[[i]][[1]], L2 = NA, y.obs = dataCa[[i]][[2]][,2], 
                                 inp.file = dataCa[[i]][[3]], eventID = names(dataCa[i]), model = pseudoHydro,
                                 par_samples_Pre = par_samples_Pre, par.tr = par.tr, par.fix = par.fix )
  
  res.LCa.bTr[[i]] <- CaPre.bTr(transf = transf, res.L = res.LCa[[i]] )
}

res.LCa.as.Pre <- list(); res.LCa.as.Pre.bTr <- list()
for (i in 1 : length(dataCa)) {
  res.LCa.as.Pre[[i]] <- CaPre.predict( L1 = NA, L2 = dataCa[[i]][[1]], y.obs = NA, 
                                        inp.file = dataCa[[i]][[3]], eventID = names(dataCa[i]), model = pseudoHydro,
                                        par_samples_Pre = par_samples_Pre, par.tr = par.tr, par.fix = par.fix )
  
  res.LCa.as.Pre.bTr [[i]] <- CaPre.bTr(transf = transf, res.L = res.LCa.as.Pre[[i]] )
}

res.LPre <- list(); res.LPre.bTr <- list()
for (i in 1 : length(dataPre)) {
  res.LPre[[i]] <- CaPre.predict( L1 = NA, L2 = dataPre[[i]][[1]], y.obs = NA, 
                                  inp.file = dataPre[[i]][[3]], eventID = names(dataPre[i]), model = pseudoHydro,
                                  par_samples_Pre = par_samples_Pre, par.tr = par.tr, par.fix = par.fix )
  
  res.LPre.bTr [[i]] <- CaPre.bTr(transf = transf, res.L = res.LPre[[i]] )
}



#-------
save.image( file = paste0(getwd(), "/outputs/", package, "_6uncer.Rdata") )
#-------
load( paste0( getwd(), "/outputs/bsl.QS_1cal.Rdata" ) )
load( paste0( getwd(), "/outputs/bsl.QS_2rr.Rdata" ) )
load( paste0( getwd(), "/outputs/bsl.QS_6uncer.Rdata" ) )
devtools::load_all(".")
#-------



#########################################################
## plots basic hydrographs for Ca events and Pre events
pdf( paste(out_dir, "/5_hydrographs_Ca.pdf", sep="") , height = 6,  width = 7 , fonts = "Times")
for ( i in 1 : length(dataPre) ) { 
  CaPre.plot.new( data_obs = dataCa[[i]], data_mod = res.LCa.bTr[[i]], eventSet = "Ca" )  
}
dev.off()

pdf( paste(out_dir, "/5_hydrographs_Ca_as_Pre.pdf", sep="") , height = 6,  width = 7 , fonts = "Times")
for ( i in 1 : length(dataPre) ) { 
  CaPre.plot.new( data_obs = dataCa[[i]], data_mod = res.LCa.as.Pre.bTr[[i]], eventSet = "Pre" )  
}
dev.off()

pdf( paste(out_dir, "/5_hydrographs_Pre.pdf", sep="") , height = 6,  width = 7 , fonts = "Times")
for ( i in 1 : length(dataPre) ) { 
  CaPre.plot.new( data_obs = dataPre[[i]], data_mod = res.LPre.bTr[[i]], eventSet = "Pre" )  
}
dev.off()



#########################################################
## statistics of the predicted model outputs
stats_CaasPre <- stats_inf( data_mod = res.LCa.as.Pre.bTr, data_obs = dataCa )
stats_Pre <- stats_inf( data_mod = res.LPre.bTr, data_obs = dataPre )



#-------
save.image( file = paste0(getwd(), "/outputs/", package, "_6uncer.Rdata") )
#-------
load( paste0( getwd(), "/outputs/bsl.QS_1cal.Rdata" ) )
load( paste0( getwd(), "/outputs/bsl.QS_2rr.Rdata" ) )
load( paste0( getwd(), "/outputs/bsl.QS_6uncer.Rdata" ) )
devtools::load_all(".")
#-------



#######################################
## plots hydrographs with performance statistics
check <- FALSE
check <- plot_hydro_stats( data_obs = dataCa, data_mod = res.LCa.as.Pre.bTr, stats = stats_CaasPre, eventSet = "Pre",
                           out_dir = out_dir )
if (check==FALSE) { dev.off() } # closes graphic device if plotting fails
check <- FALSE
check <- plot_hydro_stats( data_obs = dataPre, data_mod = res.LPre.bTr, stats = stats_Pre, eventSet = "Pre",
                           out_dir = out_dir )
if (check==FALSE) { dev.off() } # closes graphic device if plotting fails




