######################################
## J. Pastorek, SEP 2020
######################################


#######################################
## name of this package
Package <- "Adj.04.corr60min"  # stays only a global variable...
pack.dir <- substr(system.file("extdata", "gawk.exe", package = Package), 1,       # path to the package 
                   nchar(system.file("extdata", "gawk.exe", package = Package))-9)

#######################################
## reads FG and RG statistical overview and CML meta data
## ! Decide here about the statistics files (the periods used)  (it defines the rainfall-runoff events)
uni.data <- read.stats(FG.ov.path  = system.file("extdata", "flow_stats_Q2min_(12_tpl+locRGs_smooth+remRGs)+MP1_H2min.csv", package = Package),
                       RG.ov.path  = system.file("extdata", "rainfall_stats_locRGs_smooth_(12_tpl+locRGs_smooth+remRGs)+MP1_H2min_2mm_short.csv", package = Package)
)
uni.data[["CML_meta"]] <- read.csv( system.file( "extdata", "meta_25xCML_complet.csv", package = Package ), sep = ";", stringsAsFactors = F )


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
## overview of initial attenuation or rainfall data to work with
## and of methods for processing the data
rain_data_name <- "locRGs"                        # local RGs (7 RGs at 5 locations) - all available data
rain_data_name <- "locRGs_smooth"                 # local RGs (7 RGs at 5 locations) - all available data - smooth
rain_data_name <- "remRGs_all"                    # remote RGs (all 23 gauges) - data from a 3-year period
rain_data_name <- "remRGs_mean3"                  # remote RGs (mean of the closest 3 - numbers 10, 13, and 22) - data from a 3-year period

rain_data_name <- "CML4_kR"                       # CML (25 links) data from a 3-year period with k-R model
rain_data_name <- "CML5_cor15"                    # CML (19 links) data from a 3-year period corrected by remote RGs (15 min)
rain_data_name <- "CML5_cor60"                    # CML (19 links) data from a 3-year period corrected by remote RGs (60 min)
rain_data_name <- "CML6_kR_noWAA"                 # CML (25 links) data from a 3-year period with k-R model, WAA = 0
rain_data_name <- "CML7_varWAA"                   # CML (25 links) data from a 3-year period with k-R model, various WAA (from 0 to 3 dB)
rain_data_name <- "CML08_varAlp"                  # CML (25 links) data from a 3-year period with k-R model, various Alpha parameter, WAA=1.4
rain_data_name <- "CML09_drywet"                  # CML (25 links) data from a 3-year period with k-R model, new baseline B (??) , various WAA (from 0 to 3 dB), Alpha by ICU
rain_data_name <- "CML10_linWAA"                  # CML (25 links) data from a 3-year period with k-R model, WAA = 1.4  BUT multiplied by a factor growing lineraly with (tot. attenuation - baseline)
rain_data_name <- "CML11_only_A"                  # TPL3, (minus) baseline accord. to Fenicia, mean of the 2 channels
rain_data_name <- "CML12_tpl"                     # TPL3, single channels


rain_data_proc_meth <- "noProc"                   # no processing
rain_data_proc_meth <- "meanAll"                  # mean of all time series
rain_data_proc_meth <- "aggregby-min-5"           # aggregates time series to a coarser time step
rain_data_proc_meth <- "aggregbykeepLin-min-15"   # aggregates time series to a coarser time step but disaggregates it back afterwards using linear interpolation
rain_data_proc_meth <- "single-***"               # selects a single time series
rain_data_proc_meth <- "meanof-***"               # takes the mean of the partial matches of the name specified

rain_data_proc_meth <- "ThPol3"                   # for local RGs - 3 Thiessen polygons for SWMM
rain_data_proc_meth <- "mean3loc"                 # for local RGs - mean of the local RGs at the 3 locations

rain_data_proc_meth <- "keep16paper"              # keeps only data from the 16 CMLs analyzed in the paper (Pastorek et al., 2019)
rain_data_proc_meth <- "keep19paper"              # keeps only data from the 19 CMLs with paperID (Pastorek et al., 2019)
rain_data_proc_meth <- "freq-fr-25"               # keeps only data from the CMLs with circa the same frequency
rain_data_proc_meth <- "mean19"                   # mean of 19 CML time series
rain_data_proc_meth <- "mean16"                   # mean of 16 CML time series
rain_data_proc_meth <- "mean04"                   # mean of  4 CML time series
rain_data_proc_meth <- "meanChan"                 # means of the two channels of the respective CMLs
rain_data_proc_meth <- "basInterp"                # separating baseline by interpolating between the last and the next dry timestep
rain_data_proc_meth <- "basFeni-m-0.00568"        # separating baseline with a low-pass filter parameter m (Fenicia et al., 2012)
rain_data_proc_meth <- "WAAconst-WAA-1.57"        # separating WAA as a constant offset, similarly to Overeem et al. (2011)
rain_data_proc_meth <- "WAAKhaRo-C-7-d-0.125"     # separating WAA depending on the total measured A; Kharadly and Ross (2001)
rain_data_proc_meth <- "WAAGaRu-C-7-d-0.55"       # separating WAA depending on the absolute (not specific) rainfall-induced A; see Garcia-Rubia et al. (2011)
rain_data_proc_meth <- "WAAGaRuSp-C-5.5-d-1.5"    # separating WAA depending on the specific rainfall-induced A; similar to Garcia-Rubia et al. (2011)
rain_data_proc_meth <- "WAA3GaRuSp-C-7-d-1.5-z-1" # separating WAA depending on the specific rainfall-induced A; added a 3rd parameter "z"
rain_data_proc_meth <- "WAAKhaRoVal-C-7-d-1.5-z-1"# separating WAA similary like KhaRo, but depending on rainfall intensity R; added a 3rd parameter "z" 
rain_data_proc_meth <- "WAAVal-k-0.68-alp-0.34"   # separating WAA depending on specific rainfall-induced A; relation proposed by Valtr et al. (2019); numerical eq. solver
rain_data_proc_meth <- "WAAValpq-p-1.5-q-0.6"     # a modification of the relation proposed by Valtr et al. (2019)
rain_data_proc_meth <- "WAAlog-a-2-b*5"           # separating WAA depending on specific rainfall-induced A; logarithmic curve, inspired by Valtr et al. (2019)
rain_data_proc_meth <- "WAASchl-Wmax-2.3-tau-15"  # separating WAA using a time-dependent model of  Schleiss et al. (2013)
rain_data_proc_meth <- "AttSpec"                  # total attenuation to specific attenuation (divided by the path length)
rain_data_proc_meth <- "AtoR"                     # specific attenuation to rainfall



#####################################################################################################################
## adjusts using the "Pre" events

#######################################
## reads rainfall data for desired time periods,
## applies the selected processing method and deals with NAs,
scens <- as.character(c())
scens <- c( scens, "read remRGs_all__single-10--aggregby-min-60" )
remRGs <- sup.rain.data( scens = scens, periods = periods[ which( as.character(periods$st) %in% eventIDsPre ) , ] )

scens <- as.character(c())
scens <- c( scens, "read CML12_tpl__keep16paper--meanChan--basInterp" )
CML_base <- sup.rain.data( scens = scens,  periods = periods[ which( as.character(periods$st) %in% eventIDsPre ) , ] )


#######################################
## adjusts (calibrates per chunks of length N) the chosen WAA model
time_start <- proc.time()

WAA_meth  <- "WAAVal" 
par_names <- c("k", "alp")
par_init  <- c(0.68, 0.34)
lower     <- c( 0.001, 0.01 )
upper     <- c( 3, 3 ) 

N <- 5

par_opt_all <- array( data = 0, dim = c( nrow(remRGs), 16, length(par_names) ) )  # parameter values equal to 0 => zero QPE
dimnames(par_opt_all) <- list( as.character(remRGs$time), paste0("#", as.character(c(3:9, 11:19))), par_names)

dt <- difftime( remRGs$time[2], remRGs$time[1] )
for ( i_ev in unique(remRGs$id) ) {  # for each event
  
  times_ev <- remRGs$time[remRGs$id %in% i_ev]
  first <- TRUE
  for ( i_ord in 1:length(times_ev) ) {  # for each timestep ~~ chunk of length N
    i_i <- rev(times_ev) [i_ord]
    i_i_char <- as.character( rev(times_ev) ) [i_ord]
    
    i_iminusAlmostN <- i_i - (N-1)*dt
    
    if ( !first ) {
      if ( i_iminusAlmostN < min(times_ev) ) {break}  
    }
    first <- FALSE
    
    i_itoAlmostN <- seq( i_i, i_iminusAlmostN, -dt )
    remRGs_chunk <- remRGs[ remRGs$time %in% i_itoAlmostN , ]
    
    times_CML_max <- i_i + dt/2   # corresponding to the aggregby function
    times_CML_min <- i_iminusAlmostN - dt/2
    times_CML_chunk <- CML_base$time[ which( ( CML_base$time < times_CML_max ) * ( CML_base$time > times_CML_min ) == 1 ) ]
    
    CML_base_chunk  <- CML_base[ CML_base$time %in% times_CML_chunk,  ]
    
    # if zero reference, than zero parameters => zero QPE
    if ( sum(remRGs_chunk[ , ! colnames( remRGs_chunk ) %in% c("id", "time") ], na.rm = T) == 0 ) { next }  
    
    for ( i_link in colnames( CML_base )[ ! colnames( CML_base ) %in% c("id", "time") ] ) {  # for each CML
      i_link <- strsplit( i_link, "_-_" )[[1]][1]
      
      CML_base_chunk_link <- CML_base_chunk[ c( "id", "time", colnames(CML_base_chunk)[ grep( paste0(i_link, "_-_"), colnames(CML_base_chunk) ) ] ) ]
      
      # if no CML data, than NA parameter values => NA QPE
      if ( length( which( !is.na( CML_base_chunk_link[ colnames( CML_base_chunk_link )[ ! colnames( CML_base_chunk_link ) %in% c("id", "time") ] ] ) ) ) == 0 ) { 
        par_opt_all[ i_i_char, i_link,  ] <- t(as.numeric(c(NA,NA))) 
      } else {
        
        # if no attenuation observed, than zero parameters => zero QPE
        if ( sum( CML_base_chunk_link[ colnames( CML_base_chunk_link )[ ! colnames( CML_base_chunk_link ) %in% c("id", "time") ] ] , na.rm = T ) == 0 ) { next }
        
        #######################################
        ## defines the function to be minimized
        WAA_inf <- function(pars) {
          
          mod.scens <- as.character(c())
          
          pars_str <- paste(par_names[1], pars[1], par_names[2], pars[2], sep = "-")
          proc_meth <- paste0(WAA_meth, "-", pars_str, "--AttSpec--AtoR--aggregby-min-60")
          
          mod.scens <- c( mod.scens, paste0( "CML_base_chunk_link__", proc_meth ) )
          sup.rain.data <- sup.rain.data( scens = mod.scens )
          
          
          #######################################
          ## calculates performance statistics
          out_vec <- c()
          for ( i_col in colnames( sup.rain.data[   !colnames(sup.rain.data) %in% c("time", "id") ] ) ) {
            
            noNAs <-  !is.na(sup.rain.data[i_col]) 
            mod <- sup.rain.data[i_col] [ noNAs ] 
            obs <- remRGs_chunk[ , ! colnames( remRGs_chunk ) %in% c("id", "time") ] [ noNAs ]
            
            out_vec <- c( out_vec,  sqrt( mean( (mod-obs)^2 ) ) )
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
                                     control= list( max.call = 500, verbose = T, simple.function = F )
        )
        
        par_opt_all[ i_i_char, i_link,  ] <- t(Opt.precal1$par) 
      }
    }
  }
}

time_end   <- proc.time(); time_taken <- time_end - time_start; time_taken



#####################################################################################################################
## Rainfall-Rainfall evaluation

#######################################
## defines data to be evaluated
scens <- as.character(c())
scens <- c( scens, "CML_base__noProc" )
CML_R <- sup.rain.data( scens = scens )

CML_R[ ! colnames(CML_R) %in% c("id", "time") ] <- NA
colnames(CML_R) [ ! colnames(CML_R) %in% c("id", "time") ] <-  sub( "noProc", "WAAValAdj", colnames(CML_R)[ ! colnames(CML_R) %in% c("id", "time") ] )

dt <- difftime( remRGs$time[2], remRGs$time[1] )
for ( i_ev in unique(remRGs$id) ) {  # for each event
  
  times_ev <- remRGs$time[remRGs$id %in% i_ev]
  first <- TRUE
  for ( i_ord in 1:length(times_ev) ) {  # for each reference timestep ~~ chunk of length N
    i_i <- rev(times_ev) [i_ord]
    i_i_char <- as.character( rev(times_ev) ) [i_ord]
    
    i_iminusN <- i_i - (N)*dt
    i_iminusAlmostN <- i_i - (N-1)*dt
    
    if ( !first ) {
      if ( i_iminusAlmostN < min(times_ev) ) {break}  
    }
    first <- FALSE
    
    i_itoN <- seq( i_i, i_iminusN, -dt )
    times_CML_chunk <- CML_base$time[ which( ( CML_base$time <= max(i_itoN) ) * ( CML_base$time > min(i_itoN) ) == 1 ) ]
    
    CML_base_chunk  <- CML_base[ CML_base$time %in% times_CML_chunk,  ]
    CML_R_chunk     <- CML_R[ CML_R$time %in% times_CML_chunk,  ]
    
    for ( i_link in colnames( CML_base )[ ! colnames( CML_base ) %in% c("id", "time") ] ) {  # for each CML
      i_link <- strsplit( i_link, "_-_" )[[1]][1]
  
      CML_base_chunk_link <- CML_base_chunk[ c( "id", "time", colnames(CML_base_chunk)[ grep( paste0(i_link, "_-_"), colnames(CML_base_chunk) ) ] ) ]
      
      CML_R_chunk_link <- CML_R_chunk[ c( "id", "time", colnames(CML_R_chunk)[ grep( paste0(i_link, "_-_"), colnames(CML_R_chunk) ) ] ) ]
      
      # if WAA model parameters not available
      if ( sum( is.na( par_opt_all[i_i_char, i_link,] ) ) != 0 ) {
        CML_R_chunk_link[ ! colnames(CML_R_chunk_link) %in% c("id", "time") ] <- NA
      } else {
        # if WAA model parameters equal zero
        if ( sum(par_opt_all[i_i_char, i_link,]) == 0 ) {
          CML_R_chunk_link[ ! colnames(CML_R_chunk_link) %in% c("id", "time") ] <- 0
        } else {
          
          scens <- as.character(c())
          
          pars_str <- paste(par_names[1], par_opt_all[i_i_char, i_link, 1], par_names[2], par_opt_all[i_i_char, i_link, 2], sep = "-")
          proc_meth <- paste0(WAA_meth, "-", pars_str, "--AttSpec--AtoR")
          
          scens <- c( scens, paste0( "CML_base_chunk_link__", proc_meth ) )
          CML_R_chunk_link <- sup.rain.data( scens = scens )
        }
        
        CML_R_chunk[ , colnames(CML_R_chunk)[ grep( paste0(i_link, "_-_"), colnames(CML_R_chunk) ) ]  ] <- CML_R_chunk_link[ ! colnames(CML_R_chunk_link) %in% c("id", "time") ]          
      }
    }
    
    CML_R[ CML_R$time %in% CML_R_chunk$time,  ] <- CML_R_chunk
  }
  
  # if CML data timesteps after the last (aggregated) reference data timestep, then zero QPEs (end of the event <=> no rainfall)
  leftovers <- CML_R[ CML_R$id %in% i_ev  , "time" ] > max( times_ev )
  CML_R[ CML_R$id %in% i_ev  , ! colnames(CML_R) %in% c("id", "time") ] [ leftovers, ] <- 0
}

scens <- as.character(c())
scens <- c( scens, "CML_R__aggregby-min-60" )
sup.group.res_rain <- sup.rain.data( scens = scens )

#######################################
## defines data to be used as reference
sup.group.res_rain["Qobs"] <- remRGs[ !colnames(remRGs) %in% c( "id", "time" ) ] 
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


#####################################################################################################################
## Rainfall-Runoff simulations and evaluation

#######################################
## reads rainfall data for desired time periods - "Pre" events
## applies the selected processing method and deals with NAs - parameters optimized above
scens <- as.character(c())
scens <- c( scens, "CML_R__noProc" )
sup.rain.data <- sup.rain.data( scens = scens )

#######################################
## runs rainfall-runoff simulations - for all rainfall data specifies above
sup.group.res <- list(); ThPol3_check <- F
for ( i_scen in colnames(sup.rain.data)[ !colnames(sup.rain.data) %in% c("time", "id") ] ) {
  
  #######################################
  ## prepares data for SWMM (or other R-R model)
  Urquell <- system.file("extdata", "inpfile.inp", package = Package) # path to the catchment model
  
  prodata <- list(); prodata$Pre <- list();
  prodata$Pre  <- setupSWMMX( eventIDs = eventIDsPre, flow.data.proc = match_with_IDs(rainfall_datfr = flow.data.proc, IDs = eventIDsPre), 
                              Urquell = Urquell, pack.dir = pack.dir )
  
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
  Rain_File_Tab_Pre  <- setupRainFiles( rain.data.proc = sup.rain.data[ c("time", "id", prodata_rain_name) ], 
                                        pack.dir = pack.dir )
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
  sup.group.res[[ i_scen ]]  <- group.run.noInf(par = par, prodata = prodata, Rain_File_Tab = Rain_File_Tab,
                                                RRmodel = model.swmm) # model.swmm   model.1res
}

#######################################
## reshapes the modelling results
sup.group.res <- reshape.res( sup.group.res = sup.group.res)
sup.group.res$sd_Qobs <- match_with_IDs(rainfall_datfr = flow.data.proc, IDs = eventIDsPre)$sd_Q * 1000    # [m^3/s] --> [l/s]


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
statistics_inf <- list()
for ( i_subset in names(events.subsets) ) {
  hlp <- match_with_IDs( rainfall_datfr = sup.group.res, IDs = events.subsets[[i_subset]] ); 
  hlp$id <- as.character(hlp$id)
  statistics_inf[[i_subset]] <- simple.stats.sup.group( sup.group.res = hlp )  
}



#######################################
## exports the data
for (j in names(events.subsets)) {
  write.table( statistics_rain_inf [[ j ]]$overview_noEv , paste0(pack.dir, "/stat_rain_noEv_", j ,".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
  write.table( statistics_inf [[ j ]]$overview_noEv , paste0(pack.dir, "/stat_runoff_noEv_", j ,".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
}
for (j in names(statistics_rain_inf$all$overview_ev)) {
  write.table( statistics_rain_inf$all$overview_ev [[j]] , paste0(pack.dir, "/stat_rain_perEv_",  j, ".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
  write.table( statistics_inf$all$overview_ev [[j]] , paste0(pack.dir, "/stat_runoff_perEv_",  j, ".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
}

save.image( file = paste0(pack.dir, "/", Package, ".Rdata") )




#######################################
## plots hydrographs
sup.group.plot.noInf( mod.scens.to.plot = colnames(sup.rain.data)[ !colnames(sup.rain.data) %in% c("time", "id") ][c(1,2,4)], 
                      name = paste0( "adj60min"),
                      sup.group.res = sup.group.res, sup.rain.data = sup.rain.data, pack.dir = pack.dir )
dev.off()




