######################################
## J. Pastorek, OCT 2020
######################################


#######################################
## setting up the environment
print( getwd() )  # make sure the wd is set correctly (it should be on default if inside an R-project)
pack <- devtools::load_all()  # loading the R package  
package <- environmentName( pack$env )   # name of this package
out.dir <- file.path( getwd(), "outputs" )  # output directory


#######################################
## reads FG and RG statistical overview and CML meta data
## ! Decide here about the statistics files (the periods used)  (it defines the rainfall-runoff events)
uni.data <- read.stats(FG.ov.path  = system.file("extdata", "flow_stats_Q2min_(12_tpl+locRGs_smooth+remRGs)+MP1_H2min.csv", package = package),
                       RG.ov.path  = system.file("extdata", "rainfall_stats_locRGs_smooth_(12_tpl+locRGs_smooth+remRGs)+MP1_H2min_2mm_short.csv", package = package)
)
uni.data[["CML_meta"]] <- read.csv( system.file( "extdata", "meta_25xCML_complet.csv", package = package ), sep = ";", stringsAsFactors = F )


#######################################
## selects events (time periods) to work with by their number;   
## ! See e.g.    View(uni.data$RG.overview)  or   View(uni.data$FG.overview)
all.events  <-  ( 1 : length(uni.data$RG.overview[,1]) )

bad.events  <-  c( c(2:9, 84:88 ), (62), (63:69))                # damaged FG data  || too long  ||  winter  
NA_for_CML  <-  c( 1, 29, 90, 96, 97,     10:18, 37, 47, 56 )   # 10:18, 37, 47, 56 could be simulated, but only with too few CMLs      
NA_for_locRGs <- c( c(19:22, 82:83) , c(75, 79, 93) , c(80:83) )           #  NA loc RGs  ||  weird locRGs || RG1 stuffed by insect 

good.events <-  setdiff(all.events, union(bad.events, union(NA_for_locRGs, NA_for_CML)) )

periods     <-  uni.data$RG.overview[ good.events,  ] [, c("st", "en") ]

set.seed(1); which_events <- rbinom(n = length(good.events) , size = 1, prob = 0.5)
eventIDsCa    <- as.character( periods$st[ !which_events ] )                     # desired events for CALIBRATION  
eventIDsPre   <- setdiff( as.character( periods$st ) , eventIDsCa )   # desired events for PREDICTION


#######################################
## selects the discharge data to work with
## and reads the data for desired time periods
flow.data      <- read_select_data(rain_data_name = "Q_uncert" , periods = periods )
flow.data.proc <- flow.data[,c("time","id", "MP1", "sd_Q")]




#######################################
## overview of available attenuation and rainfall data 
## and of methods for processing the data

      # rain_data_name
"locRGs"                        # local RGs (7 RGs at 5 locations) - all available data
"locRGs_smooth"                 # local RGs (7 RGs at 5 locations) - all available data - smooth
"remRGs_all"                    # remote RGs (all 23 gauges) - data from a 3-year period
"remRGs_mean3"                  # remote RGs (mean of the closest 3 - numbers 10, 13, and 22) - data from a 3-year period

"CML4_kR"                       # CML (25 links) data from a 3-year period with k-R model
"CML5_cor15"                    # CML (19 links) data from a 3-year period corrected by remote RGs (15 min)
"CML5_cor60"                    # CML (19 links) data from a 3-year period corrected by remote RGs (60 min)
"CML6_kR_noWAA"                 # CML (25 links) data from a 3-year period with k-R model, WAA = 0
"CML7_varWAA"                   # CML (25 links) data from a 3-year period with k-R model, various WAA (from 0 to 3 dB)
"CML08_varAlp"                  # CML (25 links) data from a 3-year period with k-R model, various Alpha parameter, WAA=1.4
"CML09_drywet"                  # CML (25 links) data from a 3-year period with k-R model, new baseline B (??) , various WAA (from 0 to 3 dB), Alpha by ICU
"CML10_linWAA"                  # CML (25 links) data from a 3-year period with k-R model, WAA = 1.4  BUT multiplied by a factor growing lineraly with (tot. attenuation - baseline)
"CML11_only_A"                  # TPL3, (minus) baseline accord. to Fenicia, mean of the 2 channels
"CML12_tpl"                     # TPL3, single channels
"CML14_bslQuantSm"              # CML after baseline separation using "baseQuantSmooth", 16 link, mean of the two channels

      # rain_data_proc_meth
"noProc"                   # no processing
"meanAll"                  # mean of all time series
"aggregby-min-5"           # aggregates time series to a coarser time step
"aggregbykeepLin-min-15"   # aggregates time series to a coarser time step but disaggregates it back afterwards using linear interpolation
"single-***"               # selects a single time series
"meanof-***"               # takes the mean of the partial matches of the name specified

"ThPol3"                   # for local RGs - 3 Thiessen polygons for SWMM
"mean3loc"                 # for local RGs - mean of the local RGs at the 3 locations

"keep16paper"              # keeps only data from the 16 CMLs analyzed in the paper (Pastorek et al., 2019)
"keep19paper"              # keeps only data from the 19 CMLs with paperID (Pastorek et al., 2019)
"freq-fr-25"               # keeps only data from the CMLs with circa the same frequency
"mean19"                   # mean of 19 CML time series
"mean16"                   # mean of 16 CML time series
"mean04"                   # mean of  4 CML time series
"meanChan"                 # means of the two channels of the respective CMLs
"basInterp"                # separating baseline by interpolating between the last and the next dry timestep
"basFeni-m-0.00568"        # separating baseline with a low-pass filter parameter m (Fenicia et al., 2012)
"baseQuantSmooth"          # separating baseline calculated as moving quantile window through smooting of hourly data
"WAAconst-WAA-1.57"        # separating WAA as a constant offset, similarly to Overeem et al. (2011)
"WAAKhaRo-C-7-d-0.125"     # separating WAA depending on the total measured A; Kharadly and Ross (2001)
"WAAGaRu-C-7-d-0.55"       # separating WAA depending on the absolute (not specific) rainfall-induced A; see Garcia-Rubia et al. (2011)
"WAAGaRuSp-C-5.5-d-1.5"    # separating WAA depending on the specific rainfall-induced A; similar to Garcia-Rubia et al. (2011)
"WAA3GaRuSp-C-7-d-1.5-z-1" # separating WAA depending on the specific rainfall-induced A; added a 3rd parameter "z"
"WAAKhaRoVal-C-7-d-1.5-z-1"# separating WAA similary like KhaRo, but depending on rainfall intensity R; added a 3rd parameter "z" 
"WAAVal-k-0.68-alp-0.34"   # separating WAA depending on specific rainfall-induced A; relation proposed by Valtr et al. (2019); numerical eq. solver
"WAAValpq-p-1.5-q-0.6"     # a modification of the relation proposed by Valtr et al. (2019)
"WAAlog-a-2-b*5"           # separating WAA depending on specific rainfall-induced A; logarithmic curve, inspired by Valtr et al. (2019)
"WAASchl-Wmax-2.3-tau-15"  # separating WAA using a time-dependent model of  Schleiss et al. (2013)
"AttSpec"                  # total attenuation to specific attenuation (divided by the path length)
"AtoR"                     # specific attenuation to rainfall



#####################################################################################################################
## calibrates using the "Ca" events

#######################################
## reads rainfall data for desired time periods,
## applies the selected processing method and deals with NAs,
scens <- as.character(c())
scens <- c( scens, "read remRGs_mean3__aggregby-min-60" )
refRain_Ca <- sup.rain.data( scens = scens, periods = periods[ which( as.character(periods$st) %in% eventIDsCa ) , ] )

scens <- as.character(c())
scens <- c( scens, "read CML14_bslQuantSm" )
CML_bsl_Ca <- sup.rain.data( scens = scens, periods = periods[ which( as.character(periods$st) %in% eventIDsCa ) , ] )


#######################################
## calibrates the chosen WAA model
time_start <- proc.time()

WAA_meth  <- "WAAVal" 
par_names <- c("k", "alp")
par_init  <- c(0.68, 0.34)   # from Valtr et al., 2019
# par_init <-  c(0.697, 0.502)   # from "022_Adj_04_corr60min_cal_perAll" optimization
lower     <- c( 0.001, 0.01 )
upper     <- c( 3, 3 ) 
max.call  <- 1

par_opt_all <- data.frame(  matrix(vector(), 0, length(par_names))   );  colnames(par_opt_all) <- par_names;   
for ( i_link in colnames( CML_bsl_Ca )[ ! colnames( CML_bsl_Ca ) %in% c("id", "time") ] ) {
  i_link <- strsplit( i_link, "_-_" )[[1]][1]
  
  #######################################
  ## selects data for the given CML
  scens <- as.character(c())
  scens <- c( scens, paste0("CML_bsl_Ca__single-", i_link) )
  CML_bsl_Ca_link <- sup.rain.data( scens = scens )
  
  #######################################
  ## defines the function to be minimized
  WAA_inf <- function(pars) {
    
    scens <- as.character(c())
    
    pars_str <- paste(par_names[1], pars[1], par_names[2], pars[2], sep = "-")
    proc_meth <- paste0(WAA_meth, "-", pars_str, "--AttSpec--AtoR--aggregby-min-60")
    
    scens <- c( scens, paste0( "CML_bsl_Ca_link__", proc_meth ) )
    newRain <- sup.rain.data( scens = scens )
    
    
    #######################################
    ## calculates performance statistics
    out_vec <- c()
    for ( i_col in colnames( newRain[   !colnames(newRain) %in% c("time", "id") ] ) ) {
      
      noNAs <-  !is.na(newRain[i_col]) 
      mod <- newRain[i_col] [ noNAs ] 
      obs <- refRain_Ca[ , ! colnames( refRain_Ca ) %in% c("id", "time") ] [ noNAs ]
      
      out_vec <- c( out_vec,  sqrt( mean( (mod-obs)^2 ) ) )  # RMSE
    }
    
    out <- mean( out_vec )
    
    print(out)
    
    return(out)
  }
  
  
  #######################################
  ## infers the paramaters of the WAA model
  par_init <- par_init
  set.seed(42)
  Opt.precal1 <- GenSA::GenSA( par = par_init,
                               fn    = WAA_inf, 
                               lower =  lower,
                               upper =  upper,
                               control= list( max.call = max.call, verbose = T, simple.function = F )
  )
  
  par_opt_all[i_link,] <-t(Opt.precal1$par)
}

time_end   <- proc.time(); time_taken <- time_end - time_start; time_taken



#####################################################################################################################
## Evaluates using the "Pre" events

#######################################
## defines data to be evaluated
## applies the selected processing method - parameters optimized above
scens <- as.character(c())
scens <- c( scens, "read CML14_bslQuantSm" )
CML_bsl_Pre <- sup.rain.data( scens = scens, periods = periods[ which( as.character(periods$st) %in% eventIDsPre ) , ] )

scens <- as.character(c())
for ( i_link in rownames(par_opt_all) ) {
  
  pars_str  <- paste( par_names[1], par_opt_all[i_link,][1], par_names[2], par_opt_all[i_link,][2], sep = "-")
  proc_meth <- paste0(WAA_meth, "-", pars_str, "--AttSpec--AtoR")
  
  scens <- c( scens, paste0( "CML_bsl_Pre__single-", i_link, "--", proc_meth ) )
}
newRain <- sup.rain.data( scens = scens )


#########################################################################
## Rainfall-Rainfall evaluation

#######################################
## defines data to be used as the reference
scens <- as.character(c())
scens <- c( scens, "read locRGs_smooth__mean3loc--aggregby-min-60" )
refRain_Pre <- sup.rain.data( scens = scens, periods = periods[ which( as.character(periods$st) %in% eventIDsPre ) , ] )


#######################################
## defines data to be evaluated
scens <- as.character(c())
scens <- c( scens, "newRain__aggregby-min-60" )
sup.group.res_rain <- sup.rain.data( scens = scens )

sup.group.res_rain["Qobs"] <- refRain_Pre[ !colnames(refRain_Pre) %in% c( "id", "time" ) ] 
sup.group.res_rain["sd_Qobs"] <- NA[]


#######################################
## calculates timestamps
sup.group.res_rain$id <- as.character(sup.group.res_rain$id)
sup.group.res_rain["timestamp"] <- NA[]
for ( i_id in unique(sup.group.res_rain$id) ) {   # for every selected event
  i_id <- as.POSIXct(i_id, tz="UTC")
  
  event_data <- match_with_IDs(rainfall_datfr = sup.group.res_rain, IDs = i_id) # selects data for the given event using its ID
  
  time_re <- as.POSIXct(event_data$time, format="%d/%m/%Y %H:%M",tz="")
  shifts <- as.numeric( (time_re - time_re[1]) / 3600)
  origo  <- kimisc::hms.to.seconds( format( time_re[1] , format="%H:%M:%S" ) ) /3600 # [h]
  t.grid <- format(origo + shifts, nsmall = 6)
  
  sup.group.res_rain[ sup.group.res_rain$time %in% event_data$time , "timestamp" ] <- as.numeric(t.grid)
}

#######################################
## defines event subsets for statistics
events.subsets <- list(all       = periods$st[ as.character(periods$st) %in% uni.data$RG.overview$id[ which( uni.data$RG.overview$meanRain_Rmax10 > 0 )  ] ] ,
                       strong    = periods$st[ as.character(periods$st) %in% uni.data$RG.overview$id[ which( uni.data$RG.overview$meanRain_Rmax10 > 12 ) ] ] ,
                       #strongOLD = periods$st[ as.character(periods$st) %in% uni.data$RG.overview$id[ c(13, 16, 19, 21, 28) ] ] ,
                       strongest = periods$st[ as.character(periods$st) %in% uni.data$RG.overview$id[ which( uni.data$RG.overview$meanRain_Rmax10 > 20 ) ] ] ,
                       medium    = periods$st[ as.character(periods$st) %in% intersect( uni.data$RG.overview$id[ which( uni.data$RG.overview$meanRain_Rmax10 < 12 ) ],
                                                                                        uni.data$RG.overview$id[ which( uni.data$RG.overview$meanRain_Rmax10 > 5  ) ] ) ] ,
                       light     = periods$st[ as.character(periods$st) %in% uni.data$RG.overview$id[ which( uni.data$RG.overview$meanRain_Rmax10 < 5 ) ] ]
)

#######################################
## calculates performance statistics
statistics_rain_inf <- list()
for ( i_subset in names(events.subsets) ) {
  hlp <- match_with_IDs( rainfall_datfr = sup.group.res_rain, IDs = events.subsets[[i_subset]] ); 
  hlp$id <- as.character(hlp$id)
  statistics_rain_inf[[i_subset]] <- simple.stats.sup.group( sup.group.res = hlp )  
}



#########################################################################
## Rainfall-Runoff simulations and evaluation

#######################################
## runs rainfall-runoff simulations - for all rainfall data specifies above
if (exists("sup.group.res")) rm(sup.group.res); ThPol3_check <- F
for ( i_scen in colnames(newRain)[ !colnames(newRain) %in% c("time", "id") ] ) {
  print( i_scen ); time_start <- proc.time()
  
  #######################################
  ## prepares data for SWMM (or other R-R model)
  Urquell <- system.file("swmm", "inpfile.inp", package = package) # path to the catchment model
  
  prodata <- list(); prodata$Pre <- list();
  prodata$Pre  <- setupSWMMX( eventIDs = eventIDsPre, flow.data.proc = match_with_IDs(rainfall_datfr = flow.data.proc, IDs = eventIDsPre), 
                              Urquell = Urquell, package = package )
  
  #######################################
  ## prepares rainfall files for SWMM (or other R-R model)
  if ( grepl("ThPol3", i_scen) ) {
    if ( ThPol3_check == T ) { next() }
    i_scen <- substr(i_scen, 7, nchar(i_scen))
    ThPol3_check <- T
    prodata_rain_name <- paste0( c("RG1", "RG2", "RG3"), "_-_", i_scen  )
  } else { 
    prodata_rain_name <- i_scen 
  }
  Rain_File_Tab_Pre  <- setupRainFiles( rain.data.proc = newRain[ c("time", "id", prodata_rain_name) ], 
                                        package = package )
  Rain_File_Tab <- Rain_File_Tab_Pre
  
  #######################################
  ## defines SWMM model parameters to be used   
  ## ! When modifying, do not forget to change also parameters in modelSWMM and CaPre and .awk file !
  par      <- c(#mult.imp = 1, #mult.wid = 1, 
    mult.slo = 1, 
    mult.Nim = 1, 
    mult.Sim = 1 
    #mult.Spe = 1, #mult.Pze = 1, #mult.rou = 1
  )
  
  #######################################
  ## runs the model 
  hlpRRsim <- group.run.noInf( par = par, prodata = prodata, Rain_File_Tab = Rain_File_Tab,
                               RRmodel = model.swmm ) # model.swmm   model.1res
  
  #######################################
  ## reshapes the modelling results
  colnames(hlpRRsim)[ colnames(hlpRRsim) %in%  "Qmod" ] <- i_scen
  if ( ! exists(x = "sup.group.res") ) {
    sup.group.res <- hlpRRsim
  } else {
    sup.group.res <- cbind( sup.group.res, hlpRRsim[ ! colnames(hlpRRsim) %in% colnames(sup.group.res) ] )
  }
  time_end   <- proc.time(); time_taken <- time_end - time_start; print(time_taken)
}
sup.group.res$sd_Qobs <- match_with_IDs(rainfall_datfr = flow.data.proc, IDs = eventIDsPre)$sd_Q * 1000    # [m^3/s] --> [l/s]


#######################################
## calculates performance statistics
statistics_inf <- list()
for ( i_subset in names(events.subsets) ) {
  hlp <- match_with_IDs( rainfall_datfr = sup.group.res, IDs = events.subsets[[i_subset]] ); 
  hlp$id <- as.character(hlp$id)
  statistics_inf[[i_subset]] <- simple.stats.sup.group( sup.group.res = hlp )  
}



#####################################################################################################################
## exports the data
for (j in names(events.subsets)) {
  write.table( statistics_rain_inf [[ j ]]$overview_noEv , paste0(out.dir, "/stat_rain_noEv_", j ,".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
  write.table( statistics_inf [[ j ]]$overview_noEv , paste0(out.dir, "/stat_runoff_noEv_", j ,".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
}
for (j in names(statistics_rain_inf$all$overview_ev)) {
  write.table( statistics_rain_inf$all$overview_ev [[j]] , paste0(out.dir, "/stat_rain_perEv_",  j, ".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
  write.table( statistics_inf$all$overview_ev [[j]] , paste0(out.dir, "/stat_runoff_perEv_",  j, ".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
}

save.image( file = paste0(out.dir, "/", package, ".Rdata") )




#####################################################################################################################
## plots hydrographs
sup.group.plot.noInf( mod.scens.to.plot = colnames(newRain)[ !colnames(newRain) %in% c("time", "id") ][c(1,2,4)], 
                      name = paste0( "baseQuantSmooth_vs_baseInter"),
                      sup.group.res = sup.group.res, newRain = newRain, out.dir = out.dir )
dev.off()




