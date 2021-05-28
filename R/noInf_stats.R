
overview.stats <- function(statistics, events.subsets) {
  ## this function calculates overview statistics from a data frame with statistics (the results of the simple.stats function)
  ##
  ## inputs:
  ##    statistics - a data frame with statistics (the results of the simple.stats function)
  ##    events.subsets - list with vectors of event IDs
  ##
  ## outputs:
  ##    out - a list of data frames with statistical overview (one data frame per event subset)
  
  out <- list()
  for ( i_subset in names(events.subsets) ) {
    hlp <- as.data.frame( matrix(NA, ncol = length(statistics)*2, nrow = length(statistics[[1]][1,]) ) )
    rownames(hlp) <- names( statistics[[1]][1,] )
    hlp2 <- paste0( "E(", substr( names(statistics), 7 , nchar(names(statistics)) ), ")"  )
    hlp2 <- c(hlp2, paste0( "sd(", substr( names(statistics), 7 , nchar(names(statistics)) ), ")"  ) )
    colnames(hlp) <- hlp2
    
    for ( i_stat in colnames(hlp) ) {
      if ( substr(i_stat, 1, 2) == "sd" ) {
        stat_name <- substr( i_stat, 4, nchar(i_stat)-1 )  
        table_sel <- statistics[[ paste0("table_", stat_name) ]] [ which( row.names(statistics[[ paste0("table_", stat_name) ]]) %in% as.character(events.subsets[[i_subset]]) ), ]
        table_sel <- as.matrix(table_sel)
        
        i_stats_vector <- round( apply( table_sel, 2, sd, na.rm = T), 3)
      }
      
      if ( substr(i_stat, 1, 1) == "E" ) {
        stat_name <- substr( i_stat, 3, nchar(i_stat)-1 )
        table_sel <- data.frame(statistics[[ paste0("table_", stat_name) ]] [ which( row.names(statistics[[ paste0("table_", stat_name) ]]) %in% as.character(events.subsets[[i_subset]]) ), ])
        names(table_sel) <- colnames(statistics[[ paste0("table_", stat_name) ]])
        
        i_stats_vector <- round( apply( table_sel, 2, mean, na.rm = T), 3)
      }
      
      hlp[, i_stat] <- i_stats_vector  
    }
  
    out[[i_subset]] <- hlp
  }
  
  return(out)
}  

#---------------------------------------------------------------------

simple.stats.sup.group <- function(sup.group.res) {
  ## this function calculates basic statistics for a data frame with simulation results
  ##
  ## inputs:
  ##    sup.group.res  - a data frame with simulation results for many scenarios
  ##
  ## outputs:
  ##    a list of data frames (one data frame per statistic)
  
  
  event_IDs <- unique(sup.group.res$id)
  
  scen_names <- names(sup.group.res)[ -which(names(sup.group.res) %in% c("time", "id", "timestamp", "Qobs", "sd_Qobs") ) ]
  
  statistics <- list(); first1 <- T
  for ( i_scen in scen_names ) {
    print(i_scen)
    
    group.res <- sup.group.res[, names(sup.group.res) %in% c("time", "id", "timestamp", "Qobs", "sd_Qobs") ]
    group.res$Qmod <- sup.group.res[, i_scen]
    
    statistics_single_scen <- simple.stats(group.res = group.res)
    
    statistics[[ i_scen ]] <- list( event.table = statistics_single_scen$event.table, 
                                    group.overview.ev = statistics_single_scen$group.overview.ev, 
                                    group.overview.noEv = statistics_single_scen$group.overview.noEv )
    
    if ( is.na(statistics [[ i_scen ]]) ) {next}
  }
  
  for ( i_scens in 1:length(scen_names) )  {
    if ( i_scens == 1 ) {
      overview_noEv <- as.data.frame( t( statistics[[ scen_names[i_scens] ]]$group.overview.noEv ) )
    } else {
      overview_noEv <- rbind( overview_noEv, statistics[[ scen_names[i_scens] ]]$group.overview.noEv )  
    }
    rownames(overview_noEv)[i_scens ] <- scen_names[i_scens]
  }
  
  statistics_new <- list()
  for ( i_stats in names( statistics[[1]][[1]] ) ) {
    
    for ( i_scens in 1:length(scen_names) )  {
      
      if ( i_scens == 1 ) {
        hlp <-  statistics[[scen_names[i_scens]]]$event.table[i_stats]
      } else {
        hlp <- cbind( hlp, statistics[[scen_names[i_scens]]]$event.table[i_stats] )  
      }
      
      names(hlp)[ which(names(hlp) == i_stats ) ] <- scen_names[i_scens]
    }
    
    statistics_new[[ paste0("table_", i_stats) ]] <-  hlp
  }
  
  return( list( overview_ev = statistics_new , overview_noEv = overview_noEv ))
}


#---------------------------------------------------------------------


simple.stats <- function(group.res) {
  ## this function calculates basic statistics for a data frame with simulation results
  ##
  ## inputs:
  ##    group.res  - a data frame with simulation results 
  ##
  ## outputs:
  ##    a list of data frames (one data frame per statistic)
  
  
  event_IDs <- unique(group.res$id)
  
  first2 <- T
  for ( EventID in event_IDs ) {
    Q_all <- match_with_IDs( group.res, EventID )
    
    Qdata <- data.frame( Qobs = Q_all$Qobs, Qmod = Q_all$Qmod, timestamp = Q_all$timestamp, Q_obs_sd = Q_all$sd_Qobs )
    
    stats <- simple.stats.core( Qdata = Qdata )
    
    if (first2 == T) {
      event.table <-  data.frame( matrix(NA, ncol = length(stats), nrow = length(event_IDs) ), row.names = event_IDs )  
      names(event.table) <- names(stats)
      
      first2 <- F
    }
    
    event.table[EventID,] <- stats
  }
  
  group.overview.ev <- round( c( apply(event.table, 2, mean, na.rm = T) , apply(event.table, 2, sd,   na.rm = T) ), digits = 3)
  names(group.overview.ev) <- c( paste("E(", names(stats), ")", sep = "") , paste("sd(", names(stats), ")", sep = "") )
  
  Qdata <- data.frame( Qobs = group.res$Qobs, Qmod = group.res$Qmod, timestamp = group.res$timestamp, Q_obs_sd = group.res$sd_Qobs )
  group.overview.noEv <- simple.stats.core( Qdata = Qdata )
  group.overview.noEv <- group.overview.noEv[ !names(group.overview.noEv) %in% c("dQmax", "shift_Qmax") ]
 
  return(list( event.table = event.table, group.overview.ev = group.overview.ev,  group.overview.noEv =  group.overview.noEv ))
}


#------------------------------------------------------------------------------------------------------------------------


simple.stats.core <- function(Qdata) {
  ## this function calculates basic statistics for a data frame with simulation results
  ##
  ## inputs:
  ##    Qdata - a data frame with simulation results (Qobs, Qmod, timestamp)
  ##
  ## outputs:
  ##    stats_out - a vector of statistics defined below
  
  
  stats_out_names <- c("dV", "dQmax", "shift_Qmax", "NSE", "NNSE", "PCC", "CR", "RMSE", "SCC")
  
  if ( sum(is.na(Qdata$Qmod)) == length(Qdata$Qmod) ) {
    stats_out <- rep(NA, length(stats_out_names))
    names(stats_out) <- stats_out_names
    return(stats_out)
  }
  
  if ( length(which(is.na(Qdata$Qobs))) > 0 ) {
    Qdata_noNA <- Qdata[ -which(is.na(Qdata$Qobs)) , ]
  } else {
    Qdata_noNA <- Qdata
  }
  if ( length(which(is.na(Qdata_noNA$Qmod))) > 0 ) {
    Qdata_noNA <- Qdata_noNA[ -which(is.na(Qdata_noNA$Qmod)) , ]
  } else {
    Qdata_noNA <- Qdata_noNA
  }
  
  
  if ( nrow(Qdata_noNA) < 2 ) { 
    stats_out <- rep(NA, length(stats_out_names))
    names(stats_out) <- stats_out_names
    return(stats_out)
  } else {
    
    Qobs = Qdata_noNA$Qobs; Qmod = Qdata_noNA$Qmod; timestamp = Qdata_noNA$timestamp; Qobs_sd = Qdata_noNA$Q_obs_sd
    
    if (length(Qobs) != length(Qmod)) {stop("Error: different lengths of observed and modelled data.")}
    if (length(Qobs) != length(timestamp)) {stop("Error: different lengths of data and timestamps.")}
    
    
    # indivdual statistics
    
    stats_mod_names <- c("V", "Vpeak", stats_out_names) 
    
    V.obs       <- Vtot (Qdata = Qobs , timestep = timestamp  )
    Vpeak.obs   <- Vpeak(Qdata = Qobs , timestep = timestamp  )
    
    V      <- Vtot (Qdata = Qmod , timestep = timestamp  )
    Vpeak  <- Vpeak(Qdata = Qmod , timestep = timestamp  )
    
    
    dV     <- round( (V - V.obs ) / V.obs , 3 )  
    dQmax  <- round( (Vpeak - Vpeak.obs ) / Vpeak.obs , 3 )
    shift  <- time.shift.1h.window(mod = Qmod, obs = Qobs, timestamp = timestamp)
    NS     <- enesko(Qmod, Qobs)
    NNS    <- 1 / (2 - NS)
    PCC    <- cor(Qobs, Qmod, method = "pearson")
    CR     <- cont.ratio(Qobs = Qobs, Qmod = Qmod, Qobs_sd = Qobs_sd)
    RMSE   <- RMSE(Qmod = Qmod, Qobs = Qobs)
    SCC    <- SCC(Qmod = Qmod, Qobs = Qobs)
    
    stats_mod <- c(V, Vpeak, dV, dQmax, shift, NS, NNS, PCC, CR, RMSE, SCC)
    names(stats_mod) <- stats_mod_names
    
    stats_out <- stats_mod[stats_out_names]  
  }
  
  return(stats_out)
}


#---------------------------------------------------------------------

#  total discharged volume
Vtot <- function(Qdata, timestep) {
  Vtot <- 0
  
  uniq <- unique(diff(timestep))
  timestep_mode <- uniq[which.max(tabulate(match(diff(timestep), uniq)))] # determines the mode
  
  for (i in 1 : (length(timestep)-1) ) {
    if ( (timestep[i+1] - timestep[i]) == timestep_mode ) {  # to avoid calculating between various events
      Vtot <-  Vtot + ( mean(c(Qdata[i+1], Qdata[i])) * (timestep[i+1] - timestep[i])*3600 )   #    l/s * s = l  
    } else {
      next()
    }
  }
  
  ret <- round(Vtot, 3)
  return(ret)
}

#---------------------------------------------------------------------

# 4 min (2x 2-min time step) maximal discharge
Vpeak <- function(Qdata, timestep) {
  
  Qmax <- max(Qdata[-c(1, length(Qdata))])
  max.timestep <- which(Qdata==Qmax)
  if (length(max.timestep) > 1) { max.timestep <- max.timestep[2] }
  V <- 0
  for (i in (max.timestep-1) : (max.timestep)) {
    V <- V + ( mean(c(Qdata[i+1], Qdata[i])) * (timestep[i+1] - timestep[i])*3600 )   #    l/s * s = l
  }
  
  ret <- round(V, 3)
  return(ret)
}

#---------------------------------------------------------------------

# time shift of the maximum
time.shift <- function(series1, series2, timestep) {
  
  max1 <- max(series1[-c(1, length(series1))])
  max1.timestep <- which(series1==max1)
  if (length(max1.timestep) > 1) { max1.timestep <- max1.timestep[2] }
  
  max2 <- max(series2[-c(1, length(series2))])
  max2.timestep <- which(series2==max2)
  if (length(max2.timestep) > 1) { max2.timestep <- max2.timestep[2] }
  
  shift <- timestep[max1.timestep] - timestep[max2.timestep]
  ret <- round(shift, 3)
  return(ret)
}

#---------------------------------------------------------------------

# time shift of the maximum - 1 hour window
time.shift.1h.window <- function(mod, obs, timestamp) {
  
  frame <- data.frame(mod = mod, obs = obs, timestamp = timestamp)
  
  max_obs <- max(frame$obs[ -c(1, length(frame$obs)) ])
  max_obs_timestamp <- frame$timestamp[ which(frame$obs == max_obs) ]
  if (length(max_obs_timestamp) > 1) { max_obs_timestamp <- max_obs_timestamp[2] }
  
  search_window <- intersect( which( frame$timestamp > (max_obs_timestamp - 1.0001) ),    # searches only in a window +- 1 hour
                              which( frame$timestamp < (max_obs_timestamp + 1.0001) ) )   # relative to the observed maximmum
  
  max_mod <- max( frame$mod[ search_window ] )
  max_mod_timestamp <- frame$timestamp[ intersect( which( frame$mod == max_mod ) , search_window ) ]
  if (length(max_mod_timestamp) > 1) { max_mod_timestamp <- max_mod_timestamp[2] }
  
  
  shift <- max_mod_timestamp - max_obs_timestamp
  ret <- round(shift, 3)
  return(ret)
}

#---------------------------------------------------------------------

# NS efficiency without new simulations
enesko <- function (mod, obs) {
  if (length(mod) != length(obs)) { stop("error - different length of arguments") }
  
  cit <- c(); men <- c()
  for (i in 1:length(mod)) {
    cit[i] <- (obs[i] - mod[i]   )^2
    men[i] <- (obs[i] - mean(obs))^2
  }
  ret <- 1 - (sum(cit) / sum(men))
  ret <- round(ret, 3)
  
  return(ret)
}

#---------------------------------------------------------------------

# containing ratio

cont.ratio <- function(Qobs, Qmod, Qobs_sd) {

  Q_plus_2sd  <- Qobs + (2*Qobs_sd)
  Q_minus_2sd <- Qobs - (2*Qobs_sd); Q_minus_2sd[ Q_minus_2sd < 0 ] <- 0
  
  over  <- which( Qmod > Q_plus_2sd )
  under <- which( Qmod < Q_minus_2sd )
  
  out <-  1 - ( length(over) + length(under) ) / length(Qobs)
  
  return(out)
}

#---------------------------------------------------------------------

# RMSE

RMSE <- function(Qobs, Qmod) {
 
  out <-  sqrt( mean( (Qmod-Qobs)^2 ) )
  
  return(out)
}

#---------------------------------------------------------------------

# Pearson's rank correlation coefficient

SCC <- function(Qobs, Qmod) {
  
  rho <- cor.test( x = Qobs,
                   y = Qmod,
                   method = 'spearman')
  
  out <- round( rho$estimate, 5 )
  
  return(out)
}




