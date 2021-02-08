
##################################################################################################################################################

Eval_rain_runo <- function( data_rain_ref, data_rain_new , data_Q = NULL, events.subsets  ){
  ## Rainfall-Rainfall  and  Rainfall-Runoff  evaluation
  ##
  ## data_rain_ref - rainfall data to be used as the reference
  ## data_rain_new - rainfall data to be evaluated
  ## data_Q        - runoff data, to be evaluated and as well reference
  ## events.subsets - event subsets for statistics
  
  data_rain_new["Qobs"] <- data_rain_ref[ !colnames(data_rain_ref) %in% c( "id", "time" ) ] 
  data_rain_new["sd_Qobs"] <- NA[]
  
  
  #######################################
  ## calculates timestamps
  data_rain_new$id <- as.character(data_rain_new$id)
  data_rain_new["timestamp"] <- NA[]
  for ( i_id in unique(data_rain_new$id) ) {   # for every selected event
    i_id <- as.POSIXct(i_id, tz="UTC")
    
    event_data <- match_with_IDs(rainfall_datfr = data_rain_new, IDs = i_id) # selects data for the given event using its ID
    
    time_re <- as.POSIXct(event_data$time, format="%d/%m/%Y %H:%M",tz="")
    shifts <- as.numeric( (time_re - time_re[1]) / 3600)
    origo  <- kimisc::hms.to.seconds( format( time_re[1] , format="%H:%M:%S" ) ) /3600 # [h]
    t.grid <- format(origo + shifts, nsmall = 6)
    
    data_rain_new[ data_rain_new$time %in% event_data$time , "timestamp" ] <- as.numeric(t.grid)
  }
  
  #######################################
  ## calculates performance statistics
  statistics_rain <- list()
  statistics_runo <- list()
  for ( i_subset in names(events.subsets) ) {
    hlp_rain <- match_with_IDs( rainfall_datfr = data_rain_new, IDs = events.subsets[[i_subset]] ); 
    hlp_rain$id <- as.character(hlp_rain$id)
    statistics_rain[[i_subset]] <- simple.stats.sup.group( sup.group.res = hlp_rain )  
    
    if ( is.null(data_Q) ) {  statistics_runo[[i_subset]] <- NULL } else {
      hlp_runo <- match_with_IDs( rainfall_datfr = data_Q, IDs = events.subsets[[i_subset]] ); 
      hlp_runo$id <- as.character(hlp_runo$id)
      statistics_runo[[i_subset]] <- simple.stats.sup.group( sup.group.res = hlp_runo )  
    }
  }
  
  out <- list( RainData = data_rain_new, FlowData = data_Q, statistics_rain = statistics_rain, statistics_runo = statistics_runo ) 
  return( out )
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
    
    for ( i_subs in names(out$statistics_rain) ) {
      out$statistics_rain[[i_subs]]$overview_noEv <- rbind( dat1$statistics_rain[[i_subs]]$overview_noEv, dat2$statistics_rain[[i_subs]]$overview_noEv )
      
      for ( i_metr in names(  out$statistics_rain[[i_subs]]$overview_ev ) ) {
        out$statistics_rain[[i_subs]]$overview_ev[[i_metr]] <- cbind( dat1$statistics_rain[[i_subs]]$overview_ev[[i_metr]] , dat2$statistics_rain[[i_subs]]$overview_ev[[i_metr]] )
        
      }
    }
    
    for ( i_subs in names(out$statistics_runo) ) {
      out$statistics_runo[[i_subs]]$overview_noEv <- rbind( dat1$statistics_runo[[i_subs]]$overview_noEv, dat2$statistics_runo[[i_subs]]$overview_noEv )
      
      for ( i_metr in names(  out$statistics_runo[[i_subs]]$overview_ev ) ) {
        out$statistics_runo[[i_subs]]$overview_ev[[i_metr]] <- cbind( dat1$statistics_runo[[i_subs]]$overview_ev[[i_metr]] , dat2$statistics_runo[[i_subs]]$overview_ev[[i_metr]] )
        
      }
    }
    
  }
  
  
  return(out)
}


