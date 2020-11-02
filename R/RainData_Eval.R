
##################################################################################################################################################

Eval_rain_rain <- function( data_ref, data_new , events.subsets  ){
  ## Rainfall-Rainfall evaluation
  ##
  ## data_ref - data to be used as the reference
  ## data_new - data to be evaluated
  ## events.subsets - event subsets for statistics
  
  data_new["Qobs"] <- data_ref[ !colnames(data_ref) %in% c( "id", "time" ) ] 
  data_new["sd_Qobs"] <- NA[]
  
  
  #######################################
  ## calculates timestamps
  data_new$id <- as.character(data_new$id)
  data_new["timestamp"] <- NA[]
  for ( i_id in unique(data_new$id) ) {   # for every selected event
    i_id <- as.POSIXct(i_id, tz="UTC")
    
    event_data <- match_with_IDs(rainfall_datfr = data_new, IDs = i_id) # selects data for the given event using its ID
    
    time_re <- as.POSIXct(event_data$time, format="%d/%m/%Y %H:%M",tz="")
    shifts <- as.numeric( (time_re - time_re[1]) / 3600)
    origo  <- kimisc::hms.to.seconds( format( time_re[1] , format="%H:%M:%S" ) ) /3600 # [h]
    t.grid <- format(origo + shifts, nsmall = 6)
    
    data_new[ data_new$time %in% event_data$time , "timestamp" ] <- as.numeric(t.grid)
  }
  
  #######################################
  ## calculates performance statistics
  statistics_rain_inf <- list()
  for ( i_subset in names(events.subsets) ) {
    hlp <- match_with_IDs( rainfall_datfr = data_new, IDs = events.subsets[[i_subset]] ); 
    hlp$id <- as.character(hlp$id)
    statistics_rain_inf[[i_subset]] <- simple.stats.sup.group( sup.group.res = hlp )  
  }
  
  out <- list( RainData = data_new, statistics = statistics_rain_inf ) 
  return( out )
}



##################################################################################################################################################

Eval_rain_runoff <- function( data_flow, data_new , package , RRmodel = model.swmm  ){  
  ## Rainfall-Rainfall evaluation
  ##
  ## data_flow - discharge data to be used as the reference
  ## data_new - data to be evaluated
  ## package
  ## RRmodel -  model.swmm  or  model.1res
  
  #######################################
  ## runs rainfall-runoff simulations - for all rainfall data specified above
  ThPol3_check <- F
  for ( i_scen in colnames(data_new)[ !colnames(data_new) %in% c("time", "id") ] ) {
    print( i_scen )
    
    #######################################
    ## prepares data for SWMM (or other R-R model)
    Urquell <- system.file("swmm", "inpfile.inp", package = package) # path to the swmm catchment model
    
    eventIDs <- as.character( unique( data_new$id ) )
    
    prodata <- list(); prodata$Pre <- list();
    prodata$Pre  <- setupSWMMX( eventIDs = eventIDs, flow.data.proc = match_with_IDs(rainfall_datfr = data_flow, IDs = eventIDs), 
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
    Rain_File_Tab_Pre  <- setupRainFiles( rain.data.proc = data_new[ c("time", "id", prodata_rain_name) ], 
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
                                 RRmodel = RRmodel ) # model.swmm   model.1res
    
    #######################################
    ## reshapes the modelling results
    colnames(hlpRRsim)[ colnames(hlpRRsim) %in%  "Qmod" ] <- i_scen
    if ( ! exists(x = "sup.group.res") ) {
      sup.group.res <- hlpRRsim
    } else {
      sup.group.res <- cbind( sup.group.res, hlpRRsim[ ! colnames(hlpRRsim) %in% colnames(sup.group.res) ] )
    }
  }
  sup.group.res$sd_Qobs <- match_with_IDs(rainfall_datfr = data_flow, IDs = eventIDs)$sd_Q * 1000    # [m^3/s] --> [l/s]
  
  
  #######################################
  ## calculates performance statistics
  statistics_inf <- list()
  for ( i_subset in names(events.subsets) ) {
    hlp <- match_with_IDs( rainfall_datfr = sup.group.res, IDs = events.subsets[[i_subset]] ); 
    hlp$id <- as.character(hlp$id)
    statistics_inf[[i_subset]] <- simple.stats.sup.group( sup.group.res = hlp )  
  }
  
  
  out <- list( RainData = data_new, FlowData = sup.group.res, statistics = statistics_inf ) 
  return( out )
}



##################################################################################################################################################

Merge_Eval_rain_rain <- function( ...  ) { 
  
  vars <- list(...)
  
  for ( i_n in 1:(length(vars)-1) ) {
    
    if ( i_n == 1 ) { 
      dat1 <- vars[[1]] 
    } else {
      dat1 <- out
    }
    dat2 <- vars[[1+i_n]] 
    
    
    out <- dat1
    
    out$RainData <- merge(dat1$RainData, dat2$RainData)
    
    for ( i_subs in names(out$statistics) ) {
      out$statistics[[i_subs]]$overview_noEv <- rbind( dat1$statistics[[i_subs]]$overview_noEv, dat2$statistics[[i_subs]]$overview_noEv )
      
      for ( i_metr in names(  out$statistics[[i_subs]]$overview_ev ) ) {
        out$statistics[[i_subs]]$overview_ev[[i_metr]] <- cbind( dat1$statistics[[i_subs]]$overview_ev[[i_metr]] , dat2$statistics[[i_subs]]$overview_ev[[i_metr]] )
        
      }
    }
  }
  
  return(out)
}



##################################################################################################################################################

Merge_Eval_rain_runoff <- function( ...  ) { 

  vars <- list(...)
  
  for ( i_n in 1:(length(vars)-1) ) {
    
    if ( i_n == 1 ) { 
      dat1 <- vars[[1]] 
    } else {
      dat1 <- out
    }
    dat2 <- vars[[1+i_n]] 
    
    
    out <- dat1
    
    out$RainData <- merge(dat1$RainData, dat2$RainData)
    out$FlowData <- merge(dat1$FlowData, dat2$FlowData)
    
    for ( i_subs in names(out$statistics) ) {
      out$statistics[[i_subs]]$overview_noEv <- rbind( dat1$statistics[[i_subs]]$overview_noEv, dat2$statistics[[i_subs]]$overview_noEv )
      
      for ( i_metr in names(  out$statistics[[i_subs]]$overview_ev ) ) {
        out$statistics[[i_subs]]$overview_ev[[i_metr]] <- cbind( dat1$statistics[[i_subs]]$overview_ev[[i_metr]] , dat2$statistics[[i_subs]]$overview_ev[[i_metr]] )
        
      }
    }
  }
  
  return(out)
}


