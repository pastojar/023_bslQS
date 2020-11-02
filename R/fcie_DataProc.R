########################################################################################################################

sup.rain.data <- function ( scens, periods = NULL ) {
  ## prepares rainfall data for r-r simulations
  ##
  ## reads rainfall data for desired time periods,
  ## applies the selected processing method and deals with NAs,
  #
  # inputs:
  #   scens - names of the rainfall scenarios
  #
  # outputs:
  #   out - a data frame of rainfall data processed by the selected method
  
  sup.rain.data <- data.frame()
  for ( i in 1:length(scens) ) {
    print(paste0("-------------", scens[i], "------------------------"))
    
    if ( substr(scens[i], 1, 5) == "read " ) {
      scen_name <- substr(scens[i], 6, nchar(scens[i]))
      
      rain_data_name <- strsplit( scen_name, "__" ) [[1]][1]
      # reading the data
      rain.data <- read_select_data( rain_data_name = rain_data_name , periods = periods )
     
    } else {
      scen_name <- scens[i]
      
      rain_data_name <- strsplit( scen_name, "__" ) [[1]][1]
      # reading the data
      rain.data <- get(rain_data_name)
    }
    
    # processing the data
    if ( !is.na(strsplit( scen_name, "__" ) [[1]][2]) ) {
      rain.data <- process_rainfall_data( rain.data = rain.data, rain_data_proc_meth = strsplit( scen_name, "__" ) [[1]][2] )
    }

    rain.data.proc <- manage_NAs(time.series = rain.data, time.limit = 5, rem = F) # time limit in minutes
    
    # naming the columns
    # if ( ncol(rain.data.proc) == 5  ) {
    #   if ( "RG1" %in% names(rain.data.proc) ) {
    #     names(rain.data.proc)[names(rain.data.proc) %in% c("RG1", "RG2", "RG3")] <- paste(scen_name, c("RG1", "RG2", "RG3"), sep="_" )  
    #   }
    # } else {
      names(rain.data.proc)[ !names(rain.data.proc) %in% c("id", "time") ] <-
        paste0( names(rain.data.proc)[ !names(rain.data.proc) %in% c("id", "time") ] ,
                "_-_", scen_name )
    # } 
    
    if ( i == 1 ) {
      sup.rain.data <- rain.data.proc[ c( names(rain.data.proc)[ names(rain.data.proc) %in% c("time", "id") ], names(rain.data.proc)[ !names(rain.data.proc) %in% c("time", "id") ]  ) ]
    } else {
      sup.rain.data <- cbind (sup.rain.data, rain.data.proc[ !names(rain.data.proc) %in% c("time", "id") ] )    
    }
  }
  
  out <- sup.rain.data
  return(out)
}

########################################################################################################################

process_rainfall_data <- function (rain.data, rain_data_proc_meth) {
  # this function reads rainfall data and processes them using the specified method
  #
  # inputs:
  #   rain_data - name of the rain data to be processed
  #   rain_data_proc_meth - the processing method to be used
  #
  # outputs:
  #   out - a data frame of rainfall data processed by the selected method

  
  # if ( substr(rain_data_name, 1, 8) == "multiple"  ) {
  #   
  #   rain_data_combs <- eval(parse(text = substr(rain_data_name, 10, nchar(rain_data_name))))
  #   for ( i_rain_data in 1:length(rain_data_combs[,1]) ) {
  #    
  #     rain_data_name_IN      <- rain_data_combs$rain_data_names[ i_rain_data]
  #     rain_data_proc_meth_IN <- paste0("single_", rain_data_combs$names_CML[i_rain_data])
  #     
  #     hlp <- process_rainfall_data( rain_data_name =  rain_data_name_IN, 
  #                                   rain_data_proc_meth = rain_data_proc_meth_IN )
  #     if ( i_rain_data == 1 ) {
  #       maatrix <- hlp
  #       names(maatrix) [ which( names(maatrix) == "rain" ) ] <- substr(rain_data_proc_meth_IN, 8, nchar(rain_data_proc_meth_IN) )
  #     } else {
  #       maatrix <- cbind( maatrix, hlp$rain )
  #     }
  #     names(maatrix) [ which( names(maatrix) == "hlp$rain" ) ] <- substr(rain_data_proc_meth_IN, 8, nchar(rain_data_proc_meth_IN) )
  #     
  #   }
  #   
  #   out <- eval(parse(text = rain_data_proc_meth)) (rain.data = maatrix)
  #   return(out)
  # } 
 
  
  # if multiple processing methods specifed, takes one at a time
  if ( length( strsplit( rain_data_proc_meth, "--" )[[1]] ) > 1 ) {
    
    methods <- strsplit( rain_data_proc_meth, "--" )[[1]]
    proc_meth_hlp <- strsplit( methods[1] , "-" )[[1]]
    proc_meth <- proc_meth_hlp[1]   
    proc_meth_par <- as.numeric( proc_meth_hlp[-1][ c(FALSE, TRUE) ] )
    names(proc_meth_par) <- proc_meth_hlp[-1][ c(TRUE, FALSE) ]
    
    assign( "new_out",   eval(parse(text = proc_meth)) (rain.data = rain.data, proc_meth_par = proc_meth_par),  envir=globalenv() )

    rain_data_proc_meth_IN <-  sub( pattern = paste0(methods[1], "--"), replacement = "", x = rain_data_proc_meth )
    
    out <- process_rainfall_data( rain.data = new_out, 
                                  rain_data_proc_meth = rain_data_proc_meth_IN )
    rm( list = c("new_out"), envir = globalenv() )
    
  } else { # if only one method 
    
   
      proc_meth_hlp <- strsplit( rain_data_proc_meth , "-" )[[1]]
      proc_meth <- proc_meth_hlp[1]   
      proc_meth_par <- as.numeric( proc_meth_hlp[-1][ c(FALSE, TRUE) ] )
      names(proc_meth_par) <- proc_meth_hlp[-1][ c(TRUE, FALSE) ]
      
      out <- eval(parse(text = proc_meth)) (rain.data = rain.data, proc_meth_par = proc_meth_par)  
     
   
  }
  
  return(out)
}



########################################################################################################################

# processing method functions --------------------------------------------------------------

# no processing --------------------------
noProc <- function(rain.data, proc_meth_par) {
  
  out      <- rain.data[ colnames(rain.data) %in% c( "time", "id" ) ]
  out      <- cbind( out, rain.data[ names(rain.data)[ !(names(rain.data) %in% names(out)) ] ] )
  
  return(out)
}

# keeps only data from the 19 CMLs with paperID (Pastorek et al., 2019) --------------------------
keep19paper <- function(rain.data, proc_meth_par) {
  
  out      <- rain.data[ colnames(rain.data) %in% c( "time", "id" ) ]
  
  to.keep  <- c()
  to.keep <-  colnames(rain.data) [ grep( "#",  colnames(rain.data) ) ]
 
  out      <- cbind( out, rain.data[to.keep] )
  
  return(out)
}

# keeps only data from the 16 CMLs analyzed in the paper (Pastorek et al., 2019) --------------------------
keep16paper <- function(rain.data, proc_meth_par) {
  
  out      <- rain.data[ colnames(rain.data) %in% c( "time", "id" ) ]
  
  to.keep  <- c()
  
  if (  length( grep( "-A",  colnames(rain.data)[ !colnames(rain.data) %in% c("time", "id") ] ) ) > 0  ) {
    for ( i_colnames in colnames(rain.data) ) {
      if ( strsplit(i_colnames, "-")[[1]][1] %in% c("#3", "#4", "#5", "#6", "#7", "#8", "#9", "#11", "#12", "#13", "#14", "#15",  "#16", "#17", "#18", "#19") ) {
        to.keep <- c(to.keep, i_colnames)
      }
    }
  } else {
    for ( i_colnames in colnames(rain.data) ) {
      if ( strsplit(i_colnames, "_")[[1]][1] %in% c("#3", "#4", "#5", "#6", "#7", "#8", "#9", "#11", "#12", "#13", "#14", "#15",  "#16", "#17", "#18", "#19") ) {
        to.keep <- c(to.keep, i_colnames)
      }
    }
  }
  
  out      <- cbind( out, rain.data[to.keep] )
  
  return(out)
}

# keeps only data from the CMLs with circa the same frequency --------------------------
freq <- function(rain.data, proc_meth_par) {
  
  out      <- rain.data[ colnames(rain.data) %in% c( "time", "id" ) ]
  
  i_fr <- as.numeric( proc_meth_par )
  
  to.keep  <- c()
 
  colnames_short <- unlist( strsplit( colnames(rain.data)[ !colnames(rain.data) %in% c("id", "time") ], "_-_" ) )[ c(T,F) ]
  freq_vec <- apply( uni.data$CML_meta[ uni.data$CML_meta$paperNo %in% colnames_short , c("freqA", "freqB") ] ,
                     MARGIN = 1, FUN = mean, na.rm = T )
  names(freq_vec) <-  uni.data$CML_meta$paperNo[ uni.data$CML_meta$paperNo %in% colnames_short ]
  
  i_fr_ids <- names( which( abs(freq_vec - i_fr) < 1 ) )
  
  to.keep  <- colnames(rain.data)[ !colnames(rain.data) %in% c("id", "time") ] [ which( colnames_short %in% i_fr_ids  ) ]
  
  out      <- cbind( out, rain.data[to.keep] )
  
  return(out)
}

# means of the two channels of the respective CMLs --------------------------
meanChan <- function(rain.data, proc_meth_par) {
  
  out      <- rain.data[ colnames(rain.data) %in% c( "time", "id" ) ]
  
  CMLs <- substr( colnames(rain.data)[!colnames(rain.data) %in% c("time", "id")] , 1 , nchar( colnames(rain.data)[!colnames(rain.data) %in% c("time", "id")] ) -2 )
  for ( i_CML in unique(CMLs) ) {
    hlp <- apply( rain.data[ colnames(rain.data)[!colnames(rain.data) %in% c("time", "id")] [which(CMLs==i_CML)] ] , 1, mean, na.rm = T )
    out <- cbind(out, hlp)
    colnames(out)[ colnames(out) %in% "hlp" ] <- i_CML
  }
  
  return(out)
}


# aggregates time series to a coarser time step --------------------------
aggregby <- function( rain.data, proc_meth_par){
  
  by_step <- proc_meth_par['min']
  
  hlp <- rain.data[ !(colnames(rain.data) %in% c("time", "id")) ]

  t.num <- as.numeric( rain.data$time )
  t.num <- round( t.num / (by_step*60) ) * by_step*60
  tim.agr <-  as.POSIXct(t.num, origin="1970-01-01 00:00:00 UTC", tz="UTC") #time indexes for aggregtion
  
  y <- aggregate( x = hlp , by = list(time = tim.agr), FUN = mean, na.rm = T )
  y <- y[ y$time %in% rain.data$time,  ]  # removes timesteps which are not included in the input data frame
  
  if ( ! is.null(rain.data$id)  ) {
    id <- rain.data$id[ rain.data$time %in% y$time  ]
    out <- cbind( id, y  )  
  } else {
    out <- y
  }
  
  return(out)
}

# aggregates time series to a coarser time step but disaggregates it back afterwards using linear interpolation --------------------------
aggregbykeepLin <- function( rain.data, proc_meth_par){
  
  by_step <- proc_meth_par['min']

  hlp <- rain.data[ !(colnames(rain.data) %in% c("time", "id")) ]  
  
  t.num <- as.numeric( rain.data$time )
  t.num <- round( t.num / (by_step*60) ) * by_step*60
  tim.agr <-  as.POSIXct(t.num, origin="1970-01-01 00:00:00 UTC", tz="UTC") #time indexes for aggregtion
  
  y_short <- aggregate( x = hlp , by = list(time = tim.agr), FUN = mean, na.rm = T )
  y_short <- y_short[ y_short$time %in% rain.data$time,  ]  # removes timesteps which are not included in the input data frame
  
  
  # adds zero observations at the starts and ends of the events
  y_short_new <- y_short[0,]
  
  time_diff <- difftime( time1 =  y_short$time[ -1 ] ,
                         time2 =  y_short$time[ -nrow(y_short) ], 
                         units = c("mins") )
  
  if ( length( which( time_diff != by_step ) ) > 0 ) {
    event_ends <-  which( time_diff != by_step )
  } else {
    event_ends <- nrow(y_short)
  }
  
  i_e_e_old <- 0
  for ( i_e_e in event_ends ) {
    y_short_new_ev <- rbind( y_short[(i_e_e_old +1), ] , 
                             y_short[ (i_e_e_old +1) : i_e_e , ] , 
                             y_short[i_e_e, ] )
    
    y_short_new_ev$time[ 1 ] <- y_short_new_ev$time[ 1 ] - unname(by_step)*60 
    y_short_new_ev$time[ nrow(y_short_new_ev) ] <- y_short_new_ev$time[ nrow(y_short_new_ev) ] + unname(by_step)*60 
    y_short_new_ev[ c(1, nrow(y_short_new_ev)), !(colnames( y_short_new_ev) %in% c("time", "id")) ] <- 0
    
    y_short_new <- rbind( y_short_new, y_short_new_ev )
    
    i_e_e_old <- i_e_e
  }

  y <- rain.data
  y[ !(colnames(y) %in% c("time", "id")) ] <- approx( x = y_short_new$time, y = y_short_new[, !(colnames(y_short_new) %in% c("time", "id")) ] , xout = y$time, method = "linear"  ) [["y"]]
  y[ which( is.na(y[ !(colnames(y) %in% c("time", "id")) ]) ) ,  !(colnames(y) %in% c("time", "id")) ] <- 0
  
  out <- y
  
  return(out)
}

# aggregates time series to a coarser time step but disaggregates it back afterwards using a constant value --------------------------
aggregbykeepCon <- function( rain.data, proc_meth_par){
  
  by_step <- proc_meth_par['min']
  
  hlp <- rain.data[ !(colnames(rain.data) %in% c("time", "id")) ]
  
  t.num <- as.numeric( rain.data$time )
  t.num <- round( t.num / (by_step*60) ) * by_step*60
  tim.agr <-  as.POSIXct(t.num, origin="1970-01-01 00:00:00 UTC", tz="UTC") #time indexes for aggregtion
  
  y_short <- aggregate( x = hlp , by = list(time = tim.agr), FUN = mean, na.rm = T )
  y_short <- y_short[ y_short$time %in% rain.data$time,  ]  # removes timesteps which are not included in the input data frame
  
  
  # adds zero observations at the starts and ends of the events
  y_short_new <- y_short[0,]
  
  if ( length( which( diff(y_short$time, ) != by_step/60 ) ) > 0 ) {
    event_ends <-  which( diff(y_short$time, ) != by_step/60 )
  } else {
    event_ends <- nrow(y_short)
  }
  
  i_e_e_old <- 0
  for ( i_e_e in event_ends ) {
    y_short_new_ev <- rbind( y_short[(i_e_e_old +1), ] , 
                             y_short[ (i_e_e_old +1) : i_e_e , ] , 
                             y_short[i_e_e, ] )
    
    y_short_new_ev$time[ 1 ] <- y_short_new_ev$time[ 1 ] - unname(by_step)*60 
    y_short_new_ev$time[ nrow(y_short_new_ev) ] <- y_short_new_ev$time[ nrow(y_short_new_ev) ] + unname(by_step)*60 
    y_short_new_ev[ c(1, nrow(y_short_new_ev)), !(colnames( y_short_new_ev) %in% c("time", "id")) ] <- 0
    
    y_short_new <- rbind( y_short_new, y_short_new_ev )
    
    i_e_e_old <- i_e_e
  }
  
  y <- rain.data
  y[ !(colnames(y) %in% c("time", "id")) ] <- approx( x = y_short_new$time, y = y_short_new[, !(colnames(y_short_new) %in% c("time", "id")) ] , xout = y$time, method = "constant"  ) [["y"]]
  y[ which( is.na(y[ !(colnames(y) %in% c("time", "id")) ]) ) ,  !(colnames(y) %in% c("time", "id")) ] <- 0
  
  out <- y
  
  return(out)
}

# separating baseline with a low-pass filter parameter m=0.00145 (for details see eq. 3 and 17 in Fenicia et al., 2012) --------------------------
# originally coded by Martin Fencl

# Inputs: tpl - vector with total path loss [dB]
#         m   - filter parameter [-]

basFeni <- function(rain.data, proc_meth_par) {

  Atten.data <- rain.data
  
  m <- proc_meth_par['m']
  #m <- 0.00145
  #m <- 0.00568
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    B_allEv <- c()
    for( i_event in as.character(unique(Atten.data$id)) ) {
      
      tpl <- Atten.data[ which(Atten.data$id == as.POSIXct(i_event, tz = "UTC")) , i_colnames]
      
      na.ids <- which(is.na(tpl))
      if ( length(na.ids) == length(tpl) ) {
        B_allEv <- c(B_allEv, tpl)
        next()
      }
      
      Blong <- rep(NA, length(tpl))
      
      if ( length(na.ids) == 0 ) {
        Bshort <- tpl
      } else {
        Bshort <- tpl[-na.ids]
      }
      
      # estimate baseline
      for ( i in 1:length(Bshort) ) {
        if ( i == 1 ) {
          b0 <- Bshort[1]
        } else {
          Bshort[i] <- min( (1-m)*b0 + m*Bshort[i] , Bshort[i] )
          b0 <- Bshort[i]
        }
      }
      
      # implement filtered NA values
      if ( length(na.ids) == 0 ) {
        Blong <- Bshort
      } else {
        Blong[-na.ids] <- Bshort    
      }
      
      B_allEv <- c(B_allEv, Blong)
    }
    
    out[ , i_colnames] <-  Atten.data[ , i_colnames] - B_allEv
  }
  
  return(out)
}


# separating baseline calculated as moving quantile window through smooting of hourly data  --------------------------
baseQuantSmooth <- function (rain.data, q = .5, win = 7 * 24, ...) {
  # adapted from baseline_Qsmoothing from https://github.com/fenclmar/Rcmlrain/blob/master/Rcmlrain/fun_CML.r
  #
  # baseline calculated as moving quantile window through smooting of hourly data using first
  # predefined function and then quantile. The window is moving is moving by hourly steps.
  # The window length should be set up considering typical duratio of rain
  # events to ensure sufficienlty high ratio of dry weather records is within
  # window range (in each step). 
  #
  # rain.data - tpl data
  # q - quantile (0-1) of tpl hourly subset (within moving window) used for
  #            basline.
  # win - smoothing window size in hours (default is one week)
  
  tpl <- rain.data
  by_step <- 60; names(by_step) = "min"  

  tpl_h_unzoo <- aggregby( rain.data = tpl, proc_meth_par = by_step )
  tpl_h <- zoo::zoo( x = tpl_h_unzoo[ ! colnames(tpl_h_unzoo) %in% c("time", "id") ], order.by = tpl_h_unzoo$time )
  
  b_h  <- zoo::rollapply(tpl_h, win, quantile, probs = q, na.rm = T,
                         align = 'center', by = 1, partial = T)
  b_h2 <- zoo::rollapply(b_h, win, mean, na.rm = T, align = 'center',
                         by = 1, partial = T)
  
  tpl_zoo <- zoo::zoo( x = tpl[ ! colnames(tpl) %in% c("time", "id") ], order.by = tpl$time )
  
  bsl <- tpl_zoo
  bsl[ ] <- NA
  bsl[ zoo::index(b_h2), ] <- b_h2
  bsl <- zoo::na.approx(bsl, method = 'linear', rule = 2, f = .5)
  
  hlp <- tpl_zoo - bsl

  out <- zoo::fortify.zoo(model = hlp)
  out <- out[ ! colnames(out) %in% "Index" ]
  out[ out < 0 ]  <- 0
  out <- cbind( zoo::index(hlp), out )
  colnames(out)[1] <- "time"
   
  return (out)
}


# separating baseline by interpolating between the last and the next dry timestep --------------------------
basInterp <- function(rain.data, proc_meth_par) {
  
  # 1. identifies wet periods

  # wet_periods <- wet_periods_CML( rain.data = rain.data )
  
  wet_periods <- periods  #  periods is only a global variable !
  
  
  # 2. calculates the baseline for all wet periods identified above 
  
  # out - only zeros and NAs initially
  out <- rain.data
  out[ ! colnames(out) %in% c("time", "id") ] <- 0[]
  out[ ! colnames(out) %in% c("time", "id") ] <- rain.data[ ! colnames(rain.data) %in% c("time", "id") ] * out[ ! colnames(out) %in% c("time", "id") ]
  
  for ( i_event in 1:nrow(wet_periods) ) {
    
    data_event <- match_with_periods( raindata =  rain.data , periods = wet_periods[i_event,]  )
    
    for ( i_colnames in colnames(data_event)[ !(colnames(data_event) %in% c("time", "id")) ] ) {
     
      tpl <- data_event[[ i_colnames ]]
      
      # filters NA values
      na.ids <- which(is.na(tpl))
      if ( length(na.ids) == length(tpl) ) {
        data_event[[ i_colnames ]] <- tpl
        next()
      }
      
      Blong <- rep(NA, length(tpl))
      
      if ( length(na.ids) == 0 ) {
        Bshort <- tpl
      } else {
        Bshort <- tpl[-na.ids]
      }
      
      # estimates baseline
      dry_last <- Bshort[1]
      dry_next <- Bshort[length(Bshort)]
      Blong    <- seq( dry_last, dry_next, along.with = Blong )
      
      hlp <- tpl - Blong
      hlp[ hlp < 0 ] <- 0
      
      data_event[[ i_colnames ]] <- hlp
    }
    
    out[ which( out$time %in% data_event$time ) , !colnames( out ) %in% c("time", "id") ] <- data_event[ , !colnames( data_event ) %in% c("time", "id") ]

  }
  
  return(out)
}


# separating WAA using CCDFs of RG and CML observations --------------------------
WAA_CCDFs <- function(rain.data, proc_meth_par) {
  
  data_corr_pre <- rain.data
  knowledge <- get( names(proc_meth_par) )
  
  data_corr_post <- data_corr_pre[ c("id", "time") ]
  for ( i_col in colnames(data_corr_pre)[ !colnames(data_corr_pre) %in% c("id", "time") ] ) {
    
    CCDF_Ameas  <- knowledge[[ grep( i_col, names(knowledge) ) ]]$CCDF_Ameas
    Ar_teor     <- knowledge[[ grep( i_col, names(knowledge) ) ]]$Ar_teor
    Ameas_min   <- knowledge[[ grep( i_col, names(knowledge) ) ]]$Ameas_min
    
    A_meas <- data_corr_pre[[i_col]]
    ind_noNA <- !is.na(A_meas)
    ind_supercrit <- which( A_meas[ ind_noNA ] > Ameas_min )
    A_meas_supercrit <-  A_meas[ ind_noNA ] [ ind_supercrit ]

    Ameas_probs <- CCDF_Ameas( A_meas_supercrit )
    
    A_pred  <- rep( as.numeric(NA), length(A_meas) )
    A_pred[ ind_noNA ] <- 0
    A_pred[ ind_noNA ] [ ind_supercrit ] <- unname(quantile( x = Ar_teor ,  probs = Ameas_probs, na.rm = T ))
    
    data_corr_post [ i_col ] <- A_pred
  }
  
  return(data_corr_post)
}



# separating WAA using a time-dependent model of  Schleiss et al. (2013) --------------------------
WAASchl  <- function(rain.data, proc_meth_par) {

  Atten.data <- rain.data
  
  i_Wmax <- proc_meth_par["Wmax"] # [dB]
  i_tau  <- proc_meth_par["tau"]  # [min]
  # i_Wmax = 2.3  # [dB]
  # i_tau = 15  # [min]
  WAA <- Atten.data; WAA[, -which(colnames(WAA) %in% c("time", "id") ) ] <- NA
  out <- WAA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    A_tot_colnames <- Atten.data[ , c("time", i_colnames) ]
    w_colnames <- c()
    
    for ( i_events in unique(WAA$id) ) {

      A_tot_events <- A_tot_colnames[ which(Atten.data$id == i_events) , ]
      w_event <- c()
      for ( i in 1:nrow(A_tot_events) ) {

        A_tot <- A_tot_events[i, i_colnames]

        if ( i == 1 ) { # first timestep of the event

          w <- A_tot
          if ( is.na(w) ) { w <- 0 }

        } else {        # all the other timesteps

          if ( is.na(A_tot) ) {
            w <- as.numeric(NA)
          } else {

            w_new <- w0  +  (i_Wmax-w0) *3* ( as.numeric(A_tot_events$time[i]) - as.numeric(A_tot_events$time[i-1]) )/60 / i_tau

            w <- min( A_tot, i_Wmax , w_new, na.rm = TRUE  )
          }
        }

        w_event <- c(w_event, w)
        w0 <- w
        w  <- as.numeric(NULL)
      }
      w_event <- w_event * (w_event >= 0)  # sets negative values to zero

      w_colnames <- c(w_colnames, w_event)
    }
    
    WAA[ , i_colnames ] <- w_colnames
  }
  
  out[ !(colnames(Atten.data) %in% c("time", "id")) ] <- Atten.data[ !(colnames(Atten.data) %in% c("time", "id")) ] - WAA[ !(colnames(WAA) %in% c("time", "id")) ]

  return(out)
}

# separating WAA depending on specific rainfall-induced A; logarithmic curve, inspired by Valtr et al. (2019) --------------------------
WAAlog <- function(rain.data, proc_meth_par) {
  
  Atten.data <- rain.data
  
  i_a = proc_meth_par["a"]; i_b = proc_meth_par["b"]
  # i_a = 2; i_b = 5
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    if ( substr(i_colnames, 1, 1) == "#"  ) {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    } else {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    }
    
    A_r <- out[i_colnames]
    
    hlp <- Atten.data[ , i_colnames ] [which( !is.na(Atten.data[ , i_colnames]) ) ]   # without NAs, because of the Lambert function
    A_r[ which( !is.na(Atten.data[ ,i_colnames]) ), i_colnames ] <- ( ( i_a * i_b * gsl::lambert_W0( ( exp(L/(i_a*i_b)+hlp/i_a)*L ) / (i_a*i_b) ) ) - L ) / ( i_b ) 
    
    A_r <- A_r * (A_r >= 0)  # sets negative values to zero
    out[ , i_colnames] <-  A_r
  }
  
  return(out)
}

# separating WAA depending on specific rainfall-induced A; relation proposed by Valtr et al. (2019); numerical eq. solver --------------------------
WAAVal <- function(rain.data, proc_meth_par) {
  
  Atten.data <- rain.data
  
  k = proc_meth_par["k"]; alp = proc_meth_par["alp"]
  # k = 0.68; alp = 0.34
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    if ( substr(i_colnames, 1, 1) == "#"  ) {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
      alpha <- uni.data$CML_meta$alpha[ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ]
      beta <- uni.data$CML_meta$beta[ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ]
    } else {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
      alpha <- uni.data$CML_meta$alpha[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ]
      beta <- uni.data$CML_meta$beta[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ]
    }

    A_r <- out[i_colnames]
    
    hlp_in  <- Atten.data[ , i_colnames ] [ !is.na(Atten.data[ , i_colnames]) ]   # without NAs
    hlp_out <- A_r[ !is.na(Atten.data[ , i_colnames]), ]
    
    zeros <- which( hlp_in == 0 )  # deals with zeros without solving numerical eq.
    if ( length(zeros) == length(hlp_in) ) { 
      hlp_out <- rep(0, length(hlp_out))
    } else {
      if ( length(zeros) > 0 ) {
        hlp_out[ zeros ] <- 0
        rest_in <- hlp_in [ -zeros  ] 
      } else {
        rest_in  <- hlp_in  
      }
        
      rest_out <- c()
      rest_in_cut <- seq( 1, length(rest_in), by = min(10, length(rest_in)) )  # cuts the attenuation data to smaller pieces 
      for ( i_cut in 1:length(rest_in_cut)  ) {
        if ( i_cut == length(rest_in_cut) ) {
          At <- rest_in[ rest_in_cut[i_cut] : length(rest_in) ]  
        } else {
          At <- rest_in[ rest_in_cut[i_cut] : (rest_in_cut[i_cut+1]-1) ]    
        }
        
        WAAfn <- function(Ar) ( At - Ar - 2*k*( alpha*(Ar/L)^beta )^alp )
        lol <- rootSolve::multiroot( f = WAAfn, start = rep(1, length(At)), positive = T )$root
        
        rest_out <- c(rest_out, lol)
      }
      if ( length(zeros) > 0 ) {
        hlp_out[ -zeros ] <- rest_out
      } else {
        hlp_out <- rest_out
      }
    }

    A_r[ !is.na(Atten.data[ , i_colnames]), ] <- hlp_out
    
    A_r <- A_r * (A_r >= 0)  # sets negative values to zero
    out[ , i_colnames] <-  A_r
  }
  
  return(out)
}


# a modification of the relation proposed by Valtr et al. (2019); numerical eq. solver --------------------------
WAAValpq <- function(rain.data, proc_meth_par) {
  
  Atten.data <- rain.data
  
  p = proc_meth_par["p"]; q = proc_meth_par["q"]
  #p = 1.5; q = 0.6
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    if ( substr(i_colnames, 1, 1) == "#"  ) {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    } else {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    }
    
    A_r <- out[i_colnames]
    
    hlp_in  <- Atten.data[ , i_colnames ] [ !is.na(Atten.data[ , i_colnames]) ]   # without NAs
    hlp_out <- A_r[ !is.na(Atten.data[ , i_colnames]), ]
    
    zeros <- which( hlp_in == 0 )  # deals with zeros without solving numerical eq.
    hlp_out[ zeros ] <- 0
    if ( length(zeros) > 0 ) { 
      rest_in  <- hlp_in [ -zeros  ] 
    } else {
      rest_in  <- hlp_in
    }
    
    rest_out <- c()
    rest_in_cut <- seq(1, length(rest_in), by=10)  # cuts the attenuation data to smaller pieces 
    for ( i_cut in 1:length(rest_in_cut)  ) {
      if ( i_cut == length(rest_in_cut) ) {
        At <- rest_in[ rest_in_cut[i_cut] : length(rest_in) ]  
      } else {
        At <- rest_in[ rest_in_cut[i_cut] : (rest_in_cut[i_cut+1]-1) ]    
      }
      
      WAAfn <- function(Ar) ( At - Ar - 2*p*(Ar/L)^q )
      lol <- rootSolve::multiroot( f = WAAfn, start = rep(1, length(At)), positive = T )$root
      
      rest_out <- c(rest_out, lol)
    }
    if ( length(zeros) > 0 ) {
      hlp_out[ -zeros ] <- rest_out
    } else {
      hlp_out <- rest_out
    }
    
    A_r[ !is.na(Atten.data[ , i_colnames]), ] <- hlp_out
    
    A_r <- A_r * (A_r >= 0)  # sets negative values to zero
    out[ , i_colnames] <-  A_r
  }
  
  return(out)
}


# separating WAA depending on specific rainfall-induced A; similar to Garcia-Rubia et al. (2011) --------------------------
WAAGaRuSp <- function(rain.data, proc_meth_par) {
  
  Atten.data <- rain.data
  
  i_C = proc_meth_par["C"]; i_d = proc_meth_par["d"]
  # i_C = 5.5; i_d = 1.5
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    if ( substr(i_colnames, 1, 1) == "#"  ) {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    } else {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    }
    
    A_r <- out[i_colnames]
    
    hlp <- Atten.data[ , i_colnames ] [which( !is.na(Atten.data[ , i_colnames]) ) ]   # without NAs, because of the Lambert function
    A_r[ which( !is.na(Atten.data[ ,i_colnames]) ), i_colnames ] <- ( L * gsl::lambert_W0( i_C*i_d*exp(i_d*(i_C - hlp)/L) / L) + hlp*i_d - i_C*i_d ) / ( i_d )
    
    A_r <- A_r * (A_r >= 0)  # sets negative values to zero
    out[ , i_colnames] <-  A_r
  }
  
  return(out)
}

# separating WAA depending on the absolute (not specific) rainfall-induced A; see Garcia-Rubia et al. (2011) --------------------------
WAAGaRu <- function(rain.data, proc_meth_par) {
  
  Atten.data <- rain.data
  
  i_C = proc_meth_par["C"]; i_d = proc_meth_par["d"]
  # i_C = 7; i_d = 0.55
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    if ( substr(i_colnames, 1, 1) == "#"  ) {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    } else {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    }
    
    A_r <- out[i_colnames]
    
    hlp <- Atten.data[ , i_colnames ] [which( !is.na(Atten.data[ , i_colnames]) ) ]   # without NAs, because of the Lambert function
    A_r[ which( !is.na(Atten.data[ ,i_colnames]) ), i_colnames ] <- ( gsl::lambert_W0( i_C*i_d*exp(i_d*(i_C - hlp))) + hlp*i_d - i_C*i_d ) / ( i_d )
    
    
    A_r <- A_r * (A_r >= 0)  # sets negative values to zero
    out[ , i_colnames] <-  A_r
  }
  
  return(out)
}

# separating WAA depending on the total measured A; Kharadly and Ross (2001) --------------------------
WAAKhaRo <- function(rain.data, proc_meth_par) {
  
  Atten.data <- rain.data
  
  i_C = proc_meth_par["C"]; i_d = proc_meth_par["d"]
  # i_C = 7; i_d = 0.125
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {

    A_r <- Atten.data[ ,i_colnames] - i_C * (1 - exp(-i_d * Atten.data[ ,i_colnames]))
    
    A_r <- A_r * (A_r >= 0)  # sets negative values to zero
    out[ , i_colnames] <-  A_r
  }
  
  return(out)
}

# separating WAA depending on specific rainfall-induced A; added a 3rd parameter "z" --------------------------
WAA3GaRuSp <- function(rain.data, proc_meth_par) {
  
  Atten.data <- rain.data
  
  i_C =  proc_meth_par["C"]
  i_d =  proc_meth_par["d"]
  i_z =  proc_meth_par["z"]
  # i_C = 5.5; i_d = 1.5
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    if ( substr(i_colnames, 1, 1) == "#"  ) {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    } else {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    }
    
    A_r <- out[i_colnames]
    
    hlp_in  <- Atten.data[ , i_colnames ] [ !is.na(Atten.data[ , i_colnames]) ]   # without NAs
    hlp_out <- A_r[ !is.na(Atten.data[ , i_colnames]), ]
    
    zeros <- which( hlp_in == 0 )  # deals with zeros without solving numerical eq.
    hlp_out[ zeros ] <- 0
    
    rest_in  <- hlp_in [ -zeros  ]
    rest_out <- c()
    rest_in_cut <- seq(1, length(rest_in), by=10)  # cuts the attenuation data to smaller pieces 
    for ( i_cut in 1:length(rest_in_cut)  ) {
      if ( i_cut == length(rest_in_cut) ) {
        At <- rest_in[ rest_in_cut[i_cut] : length(rest_in) ]  
      } else {
        At <- rest_in[ rest_in_cut[i_cut] : (rest_in_cut[i_cut+1]-1) ]    
      }
      
      WAAfn <- function(Arsp) ( At - Arsp*L - i_C*(1 - exp(-i_d * Arsp^i_z )) )
      lol <- rootSolve::multiroot( f = WAAfn, start = rep(1, length(At)), positive = T )$root
      lol <- lol*L  # Arsp --> Ar
      
      rest_out <- c(rest_out, lol)
    }
    hlp_out[ -zeros ] <- rest_out
    
    A_r[ !is.na(Atten.data[ , i_colnames]), ] <- hlp_out
    
    A_r <- A_r * (A_r >= 0)  # sets negative values to zero
    out[ , i_colnames] <-  A_r
   
  }
  
  return(out)
}

# separating WAA depending on the absolute (not specific) rainfall-induced A; see Garcia-Rubia et al. (2011) --------------------------
WAA3GaRu <- function(rain.data, proc_meth_par) {
  
  Atten.data <- rain.data
  
  i_C =  proc_meth_par["C"]
  i_d =  proc_meth_par["d"]
  i_z =  proc_meth_par["z"]
  # i_C = 5.5; i_d = 1.5
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    if ( substr(i_colnames, 1, 1) == "#"  ) {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    } else {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
    }
    
    A_r <- out[i_colnames]
    
    hlp_in  <- Atten.data[ , i_colnames ] [ !is.na(Atten.data[ , i_colnames]) ]   # without NAs
    hlp_out <- A_r[ !is.na(Atten.data[ , i_colnames]), ]
    
    zeros <- which( hlp_in == 0 )  # deals with zeros without solving numerical eq.
    hlp_out[ zeros ] <- 0
    
    rest_in  <- hlp_in [ -zeros  ]
    rest_out <- c()
    rest_in_cut <- seq(1, length(rest_in), by=10)  # cuts the attenuation data to smaller pieces 
    for ( i_cut in 1:length(rest_in_cut)  ) {
      if ( i_cut == length(rest_in_cut) ) {
        At <- rest_in[ rest_in_cut[i_cut] : length(rest_in) ]  
      } else {
        At <- rest_in[ rest_in_cut[i_cut] : (rest_in_cut[i_cut+1]-1) ]    
      }
      
      WAAfn <- function(Ar) ( At - Ar - i_C*(1 - exp(-i_d * Ar^i_z )) )
      lol <- rootSolve::multiroot( f = WAAfn, start = rep(1, length(At)), positive = T )$root
      
      rest_out <- c(rest_out, lol)
    }
    hlp_out[ -zeros ] <- rest_out
    
    A_r[ !is.na(Atten.data[ , i_colnames]), ] <- hlp_out
    
    A_r <- A_r * (A_r >= 0)  # sets negative values to zero
    out[ , i_colnames] <-  A_r
    
  }
  
  return(out)
}

# separating WAA depending on the total measured A; Kharadly and Ross (2001) --------------------------
WAA3KhaRo <- function(rain.data, proc_meth_par) {
  
  Atten.data <- rain.data
  
  i_C =  proc_meth_par["C"]
  i_d =  proc_meth_par["d"]
  i_z =  proc_meth_par["z"]
  # i_C = 7; i_d = 0.125
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    A_r <- Atten.data[ ,i_colnames] - i_C * (1 - exp(-i_d * Atten.data[ ,i_colnames]^i_z))
    
    A_r <- A_r * (A_r >= 0)  # sets negative values to zero
    out[ , i_colnames] <-  A_r
  }
  
  return(out)
}

# separating WAA similary like KhaRo, but depending on rainfall intensity R; added a 3rd parameter "z" --------------------------
WAAKhaRoVal <- function(rain.data, proc_meth_par) {
  
  Atten.data <- rain.data
  
  i_C =  proc_meth_par["C"]
  i_d =  proc_meth_par["d"]
  i_z =  proc_meth_par["z"]
  # i_C = 7; i_d = 0.1; i_z = 1
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    if ( substr(i_colnames, 1, 1) == "#"  ) {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
      alpha <- uni.data$CML_meta$alpha[ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ]
      beta  <- uni.data$CML_meta$beta [ which(uni.data$CML_meta$paperNo == strsplit(i_colnames, "_")[[1]][1]) ]
    } else {
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 )  # [km]
      alpha <- uni.data$CML_meta$alpha[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ]
      beta  <- uni.data$CML_meta$beta [ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ]
    }
    a <- (1/alpha)^(1/beta)
    b <- (1/beta)
    
    A_r <- out[i_colnames]
    
    hlp_in  <- Atten.data[ , i_colnames ] [ !is.na(Atten.data[ , i_colnames]) ]   # without NAs
    hlp_out <- A_r[ !is.na(Atten.data[ , i_colnames]), ]
    
    zeros <- which( hlp_in == 0 )  # deals with zeros without solving numerical eq.
    hlp_out[ zeros ] <- 0
    
    rest_in  <- hlp_in [ -zeros  ]
    rest_out <- c()
    rest_in_cut <- seq(1, length(rest_in), by=10)  # cuts the attenuation data to smaller pieces 
    for ( i_cut in 1:length(rest_in_cut)  ) {
      if ( i_cut == length(rest_in_cut) ) {
        At <- rest_in[ rest_in_cut[i_cut] : length(rest_in) ]  
      } else {
        At <- rest_in[ rest_in_cut[i_cut] : (rest_in_cut[i_cut+1]-1) ]    
      }
      
      WAAfn <- function(R) ( At - L*a*(R^b) - i_C*(1 - exp(-i_d * R ^i_z )) )
      lol <- rootSolve::multiroot( f = WAAfn, start = rep(1, length(At)), positive = T )$root
      
      lol <- a*lol^b  # R --> Arsp
      lol <- lol*L    # Arsp --> Ar
      
      rest_out <- c(rest_out, lol)
    }
    hlp_out[ -zeros ] <- rest_out
    
    A_r[ !is.na(Atten.data[ , i_colnames]), ] <- hlp_out
    
    A_r <- A_r * (A_r >= 0)  # sets negative values to zero
    out[ , i_colnames] <-  A_r
    
  }
  
  return(out)
}

# separating WAA as a constant offset, similarly to Overeem et al. (2011) --------------------------
WAAconst <- function(rain.data, proc_meth_par) {
  
  Atten.data <- rain.data
  
  WAA = proc_meth_par["WAA"];
  # WAA <- 1.57
  out <- Atten.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(Atten.data)[ !(colnames(Atten.data) %in% c("time", "id")) ] ) {
    
    A_r <- Atten.data[ ,i_colnames] - WAA 
    
    A_r <- A_r * (A_r >= 0)  # sets negative values to zero
    out[ , i_colnames] <-  A_r
  }
  
  return(out)
}

# attenuation to specific attenuation (divided by the path length) --------------------------
AttSpec <- function(rain.data, proc_meth_par) {
  
  out <- rain.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(rain.data)[ !(colnames(rain.data) %in% c("time", "id")) ] ) {
    
    if ( grepl(pattern = "#", x = i_colnames) ) {
      i_paperNo <- paste0( "#", strsplit( strsplit( i_colnames, "#" )[[1]][2] , "_-_" )[[1]][1] )
      L <- uni.data$CML_meta$length[ which(uni.data$CML_meta$paperNo ==  i_paperNo) ] / 1000 # [km]
    } else {
      L <- uni.data$CML_meta$length[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ] / 1000 # [km]
    }
    
    Ar_spec <- out[i_colnames]
    Ar_spec[, i_colnames ] <- rain.data[, i_colnames]  /  ( L ) 
    Ar_spec <- Ar_spec * (Ar_spec >= 0)  # sets negative values to zero
    
    out[ , i_colnames] <-  Ar_spec
  }
  
  return(out)
}

# specific attenuation to rainfall --------------------------
AtoR <- function(rain.data, proc_meth_par) {

  out <- rain.data; out[, -which(colnames(out) %in% c("time", "id") ) ] <- NA
  for ( i_colnames in colnames(rain.data)[ !(colnames(rain.data) %in% c("time", "id")) ] ) {
    
    if ( grepl(pattern = "#", x = i_colnames) ) {
      i_paperNo <- paste0( "#", strsplit( strsplit( i_colnames, "#" )[[1]][2] , "_-_" )[[1]][1] )
      alpha <- uni.data$CML_meta$alpha[ which(uni.data$CML_meta$paperNo == i_paperNo) ]
      beta <-  uni.data$CML_meta$beta [ which(uni.data$CML_meta$paperNo == i_paperNo) ]
    } else {
      alpha <- uni.data$CML_meta$alpha[ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ]
      beta <-  uni.data$CML_meta$beta [ which(uni.data$CML_meta$id == strsplit(i_colnames, "_")[[1]][1]) ]
    }
    
    out[ , i_colnames] <-  alpha * ( rain.data[ , i_colnames] ^ beta )  
    
    out[ , i_colnames] <- out[ , i_colnames] * (out[ , i_colnames] >= 0)  # sets negative values to zero
    out[ , i_colnames] <- round( out[ , i_colnames] , digits = 3 )
  }
  
  return(out)
}

# for local RGs - 3 Thiessen polygons for SWMM --------------------------
ThPol3 <- function(rain.data, proc_meth_par) {
  
  out      <- data.frame(time = rain.data$time, id = rain.data$id)
  out$RG1  <- rain.data$RG1_PS
  out$RG2  <- apply(rain.data[ , c("RG2_SC", "RG5_SC") ], 1, mean, na.rm = T  ) # mean of 2 RGs at Sport Centrum
  out$RG3  <- apply(rain.data[ , c("RG3_GP", "RG6_GP") ], 1, mean, na.rm = T  ) # mean of 2 RGs at Green Park
  
  out[ c("RG1", "RG2", "RG3") ] <- round( out[ c("RG1", "RG2", "RG3") ], digits = 3 )
  
  return(out)
}

# for local RGs - the same as ThPl3, only a different name --------------------------
keep3 <- function(rain.data, proc_meth_par) {
  
  out      <- data.frame(time = rain.data$time, id = rain.data$id)
  out$RG1  <- rain.data$RG1_PS
  out$RG2  <- apply(rain.data[ , c("RG2_SC", "RG5_SC") ], 1, mean, na.rm = T  ) # mean of 2 RGs at Sport Centrum
  out$RG3  <- apply(rain.data[ , c("RG3_GP", "RG6_GP") ], 1, mean, na.rm = T  ) # mean of 2 RGs at Green Park
  
  out[ c("RG1", "RG2", "RG3") ] <- round( out[ c("RG1", "RG2", "RG3") ], digits = 3 )
  
  return(out)
}


########################################################################################################################

# processing method functions with a single vector of rainfall as output

# selects a single time series;  --------------------------
single <- function(rain.data, proc_meth_par) {
  
  time.series <- names(proc_meth_par)
  
  out      <- rain.data[ colnames(rain.data) %in% c( "time", "id" ) ]
  
  if ( length( which( names(rain.data)==time.series ) ) == 1 ) {   # if a column matches the name exactly
    out[ , time.series ] <- rain.data[ , time.series ]
  } else {
    out[ , time.series ] <- rain.data[ , grep(time.series, names(rain.data)) ]    #  TEMPORARY, REMOVE SOON !  
  }
  
  out[ , time.series ] <- round( out[ , time.series ], digits = 3 )
  
  return(out)
}


# takes the mean of the partial matches of the name specified --------------------------
meanof <- function(rain.data, proc_meth_par) {
  
  time.series <- names(proc_meth_par)
  
  out      <- data.frame(time = rain.data$time, id = rain.data$id)
  
  cols_to_take <- grep(time.series, names(rain.data))
  if ( length( cols_to_take ) >= 1 ) {
    
    hlp <- apply( rain.data[ , cols_to_take ], 1, mean, na.rm = T )
    out[ , time.series ] <- hlp
    
  } else {
    break()
  }
  
  out[ , time.series ] <- round( out[ , time.series ], digits = 3 )
  
  return(out)
}


# mean of the RGs at the 3 locations --------------------------
mean3loc <- function(rain.data, proc_meth_par) {
  
  hlp      <- data.frame(time = rain.data$time, id = rain.data$id)
  
  hlp$RG1  <- rain.data$RG1_PS
  hlp$RG2  <- apply(rain.data[ , c("RG2_SC", "RG5_SC") ], 1, mean, na.rm = T  ) # mean of 2 RGs at Sport Centrum
  hlp$RG3  <- apply(rain.data[ , c("RG3_GP", "RG6_GP") ], 1, mean, na.rm = T  ) # mean of 2 RGs at Green Park
  
  out      <- data.frame(time = rain.data$time, id = rain.data$id)
  out$rain <- apply( hlp[ , c("RG1", "RG2", "RG3") ], 1, mean, na.rm = T  )
  out$rain <- round( out$rain, digits = 3 )
  
  return(out)
}

# mean of all time series--------------------------
meanAll <- function(rain.data, proc_meth_par) {
  
  out <- data.frame(time = rain.data$time, id = rain.data$id)
  
  if ( sum( !colnames(rain.data) %in% c("time", "id") ) == 1 ) {
    out$meanAll <- rain.data[ , !colnames(rain.data) %in% c("time", "id")]
  } else {
    out$meanAll <- apply(rain.data[ , !colnames(rain.data) %in% c("time", "id")], 1, mean, na.rm = T  )
  }
  
  out$meanAll <- round( out$meanAll, digits = 3 )
  
  return(out)
}

# mean of 19 CML time series--------------------------
mean19 <- function(rain.data, proc_meth_par) {
  
  out      <- data.frame(time = rain.data$time, id = rain.data$id)
  to.keep  <- colnames(rain.data) %in% c("#1", "#2", "#3", "#4", "#5", "#6", "#7", "#8", "#9", "#10", "#11", "#12", "#13", "#14", "#15",  "#16", "#17", "#18", "#19")
  out$rain <- apply(rain.data[ , to.keep ], 1, mean, na.rm = T  )
  out$rain <- round( out$rain, digits = 3 )
  
  return(out)
}

# mean of 16 CML time series--------------------------
mean16 <- function(rain.data, proc_meth_par) {
  
  out      <- data.frame(time = rain.data$time, id = rain.data$id)
  to.keep  <- colnames(rain.data) %in% c("#3", "#4", "#5", "#6", "#7", "#8", "#9", "#11", "#12", "#13", "#14", "#15",  "#16", "#17", "#18", "#19")
  out$rain <- apply(rain.data[ , to.keep ], 1, mean, na.rm = T  )
  out$rain <- round( out$rain, digits = 3 )
  
  return(out)
}

# mean of 4 CML time series--------------------------
mean04 <- function(rain.data, proc_meth_par) {
  
  out      <- data.frame(time = rain.data$time, id = rain.data$id)
  to.keep  <- colnames(rain.data) %in% c("#1", "#2", "#5", "#7")
  out$rain <- apply(rain.data[ , to.keep ], 1, mean, na.rm = T  ) 
  out$rain <- round( out$rain, digits = 3 )
  
  return(out)
}



########################################################################################################################

manage_NAs <- function(time.series, time.limit, rem) {
  # this is a roof function for the "manage_NAs_core" function;
  # it just cuts the input data time series into single events
  #
  # inputs:
  #   time.series - a data frame of  time series (e.g. rainfall data - time , data... , id )
  #   time.limit  - the maximium allowed length of the gap in minutes
  #   rem         - boolean, Should the events with too much NAs be removed?
  #
  # outputs:
  #   out - a data frame of time series
  
  if ( typeof(time.series[1,1]) != "double" ) { stop("The first row is not of type double.") }
  out <- time.series[ -(1 : length(time.series[,1])), ]
  
  uncorrected <- c()
  
  if ( ! "id" %in% colnames(time.series) ) {
    fc_core <- manage_NAs_core(time.series = time.series, time.limit = time.limit, rem = rem)
    out <- fc_core[["out"]]
  } else {
    
    for ( event_ID in unique(as.character(time.series$id)) ) {
      event_ID <- as.POSIXct(event_ID, tz="UTC")
      time.series.event <- match_with_IDs(time.series, event_ID)
      
      fc_core <- manage_NAs_core(time.series = time.series.event, time.limit = time.limit, rem = rem)
      out_core <- fc_core[["out"]]
      if ( fc_core[["uncorrected"]] == T ) {
        uncorrected <- c(uncorrected, as.character(event_ID))
      }
      
      out <- rbind(out, out_core)
    }
  }
  
  # if ( length(uncorrected) > 0 ) {            # if there are some time steps with NA rainfall data
  #   
  #   print( "manage_NAs - These events contain too long gaps (NAs or NaNs) and therefore have not been corrected.") 
  #   print(paste0("manage_NAs - ", uncorrected))
  #   if (rem == T) {
  #     print("manage_NAs - These events were removed from data: ")
  #   }
  #   else {
  #     print("manage_NAs - Data for these events still contain NAs or NaNs: ")
  #   } 
  # }
  
  
  return(out)
}


########################################################################################################################

manage_NAs_core <- function(time.series, time.limit, rem) {
  # this function - removes NAs from a time series by interpolation if the gap is not too long;
  #               - if the gap is too long, it can return an empty data frame or ignore this gap
  #
  # inputs:
  #   time.series - a data frame of  time series (e.g. rainfall data - time , data... , id )
  #   time.limit  - the maximium allowed length of the gap in minutes
  #   rem         - boolean; If there is a too long gap, should only an empty data frame be returned ?
  #
  # outputs:
  #   out           - a data frame of time series
  #   uncorrected   - boolean, indicating wheather the data frame was (OK or succesfully corrected) or (there have been a too long gap)
  
  if ( typeof(time.series[1,1]) != "double" ) { stop("The first row is not of type double.") }
  out <- time.series
  
  for ( i in colnames(time.series)[ ! colnames(time.series) %in% c("id", "time") ] ) {  # for all columns except the first one (that should be time)
    NApos <- which(is.na(time.series[,i]))
    if ( length( NApos ) > 0 ) {
      
      OK_last_pos <- NA; OK_next_pos <- NA
      for ( j in 1 : length( NApos ) ) {
        if ( is.na(OK_last_pos) ) {
          if ( NApos[j] > 1 && !is.na(time.series[NApos[j]-1, i]) ) {
            OK_last_pos <- NApos[j]-1
          }
          else { next }  
        }
        
        if ( is.na( time.series[NApos[j]+1, i] ) ) {next}
        else { 
          OK_next_pos <- NApos[j]+1
          #time_diff <- time.series[ (OK_next_pos-1) , 1 ] - time.series[ (OK_last_pos-1) , 1 ]    # 19.12.2018
          time_diff <- time.series$time[ (OK_next_pos) ] - time.series$time[ (OK_last_pos) ]
          
          if ( (as.numeric(time_diff , units = "secs" ) - 1*60) <= time.limit*60 ) {
            hlp <- seq(from = time.series[ OK_last_pos , i ] , to = time.series[ OK_next_pos , i ] , 
                       length.out = (OK_next_pos - OK_last_pos + 1) )
            hlp <- hlp[ -c(1,length(hlp)) ]
            out[ (OK_last_pos+1) : (OK_next_pos-1) , i ] <- hlp
          }
          
          OK_last_pos <- NA; OK_next_pos <- NA  
        }
      }
      
    }
    
    uncorrected <- F
    if ( length( which(is.na(out[,i])) ) > 0 ) {
      uncorrected <- T
      if ( rem == T) {        
        out <- out[ -(1 : length(out[,1])), ]  
      }
    }
    
  }
  
  ret <- list(out = out, uncorrected = uncorrected)
  
  return(ret)
}
  
  
  