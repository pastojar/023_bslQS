
########################################################################################################################


read_periods <- function(fil.ev) {
  ## this function reads periods of interest (starts and ends)
  ##
  ## inputs:
  ## fil.ev  - (one or more) path to the file with starts and ends
  ##
  ## outputs:
  ## event - a list of data frames with two columns (st and end as POSIXct)
  
  event <- list()
  for (i in 1:length(fil.ev)) { 
    event[[i]] <- read.table(fil.ev[i], sep=";", header=T, stringsAsFactors=F)
    event[[i]][ ,1] <- as.POSIXct(event[[i]][ ,1], tz="UTC")
    event[[i]][ ,2] <- as.POSIXct(event[[i]][ ,2], tz="UTC")
  }
  
  return(event)
}


########################################################################################################################


read_connect_RG_files_XY <- function(path.rainfall){
  ## function to read and connect RG files (regular time step data) available in a given folder 
  ##
  ## inputs:
  ## path.rainfall - (one or more) path to a folder RG files
  ##
  ## outputs:
  ## dat_list - a list of data frames with two columns (time in UTC, R in mm/h)
  
  dat_list <- list()
  for (j in 1:length(path.rainfall)) {
    path <- path.rainfall[j]
    f.nam <- dir(path)
    id <- which(substr(f.nam, nchar(f.nam)-2, nchar(f.nam))=="dat")
    
    #basic checks on input parameters
    if(length(id)==0){stop("no RG dat files available in the folder")}
    
    dat <- read.table(paste(path[1], "//", f.nam[1], sep=""), header=T) # 1st iteration
    if(length(id)>1){
      for(i in 2:length(id)){                                           # other iterations
        dat.i <- read.table(paste(path, "//", f.nam[i], sep=""), header=T)
        dat <- rbind(dat, dat.i)
      }
    }
    
    dat[,1] <- as.POSIXct(dat[,1], tz="UTC")
    
    dat_list[[j]] <- dat
    names(dat_list)[j] <- substr( path.rainfall[j], nchar(path.rainfall[j])-5,  nchar(path.rainfall[j]) )
  }
  
  return(dat_list)
}


########################################################################################################################


list_to_dataframe <- function( RG, timestep ) {
  ## function to create a unified data frame  
  ##
  ## inputs:
  ## RG  - rainfall data in the format of the output of read_connect_RG_files_XY function (a list of data frames)
  ## timestep - desired time step for the output [seconds] 
  ##
  ## outputs:
  ## RG.total - a data frame based on the input list RG
  ##                                         - with a unified time span (over a total observed period)
  ##                                         
  
  
  # defines total observed period
  starts <- c()   ; ends <- c()
  for (i in 1:length(RG)) {
    starts[i] <- as.character(RG[[i]] [1,1])
    ends[i]   <- as.character(RG[[i]] [length(RG[[i]][,1]),1])
  }
  starts <- min( as.POSIXct( starts, origin="1970-01-01 00:00:00 UTC", tz="UTC") )
  ends   <- max( as.POSIXct( ends,   origin="1970-01-01 00:00:00 UTC", tz="UTC") )
  
  T.ser <- seq( starts, ends , timestep ) 
#  T.ser <- seq( RG[[ longest_timeseries ]] [1, "time"] , RG[[ longest_timeseries ]] [length(RG[[ longest_timeseries ]][,"time"]), "time"] , 60 ) 
  
  
  # creates a data frame with data from all RGs with the same time stemps (over the total observed period)
  RG.total <- as.data.frame( matrix(NA, ncol = length(RG)+1, nrow = length(T.ser)) )  
  colnames(RG.total) <- c("time", names(RG) )
  RG.total$time <- T.ser
  for (i in 1:length(RG)) {             
    print(paste("time series from ", names(RG)[i], sep=""))
    
    id <- which( is.element(T.ser, RG[[i]][,1]) == T ) 
    RG.total[id, i+1] <- RG[[i]][,2]
  }
  
  return( RG.total ) 
}


########################################################################################################################


add_mean_rainfall <- function( rainfall ) {
  ## function to calculate mean areal rainfall 
  ##
  ## inputs:
  ## rainfall  - rainfall data in the format of the output of list_to_dataframe function (a data frame)
  ##
  ## outputs:
  ## out - a data frame like the one in input, but with added areal (mean) rainfall

  out <- rainfall
  
  #  makes space for the menan time series in the 2nd column
  out_names <- names(out)[2:(length(out[1,])+0)]
  out[ , length(out[1,]) + 1 ] <- NA
  names(out)[   3:length(out[1,]) ] <- out_names
  out       [ , 3:length(out[1,]) ] <- out[ , 2:(length(out[1,])-1) ]
  out [,2] <- NA
  names(out)[2] <- "meanRain"
  
  # calulates the mean
  rain_time_series <- which( names(out) != ( "time") & names(out) != "id" & names(out) != "meanRain" )
  out[, 2] <-  apply( out[ ,  rain_time_series] , 1 , mean , na.rm=T )
  
  return( out )
}

########################################################################################################################
########################################################################################################################

match_with_periods <- function( raindata, periods ) {
  ## the function selects rainfall data based on periods of interest
  ##
  ## inputs:
  ## raindata - rainfall data in the format of the output of areal_rainfall function (a data frame)
  ## periods  - a data frame with two columns (st and end as POSIXct)
  ##
  ## outputs:
  ## raindata_events - a data frame (like the raindata input) but only for selected periods,
  ##                   including a column with event IDs ("st" of every event)
  
  if ( is.character( periods[1,1] ) ) {
    periods$st <- as.POSIXct(periods$st, tz = "UTC" )
    periods$en <- as.POSIXct(periods$en, tz = "UTC" )
  }
  
  out <- raindata[0,]
 
  for ( i_period in 1:length(periods$st) ) {
    
    i_id <- periods$st[i_period]
    
    hlp  <- intersect(which( raindata$time <= periods$en[i_period] ) , which( raindata$time >= periods$st[i_period] ) )
    hlp2 <- hlp[1]-1
    hlp2 <- c( hlp2, hlp )
    
    i_raindata <- raindata[ hlp2 ,  ]
    i_raindata$id <- i_id
    
    out <- rbind( out, i_raindata )
  }
  
  
  return(out)
}


zz_match_with_periods <- function( raindata, periods ) {
  ## the function selects rainfall data based on periods of interest
  ##
  ## inputs:
  ## raindata - rainfall data in the format of the output of areal_rainfall function (a data frame)
  ## periods  - a data frame with two columns (st and end as POSIXct)
  ##
  ## outputs:
  ## raindata_events - a data frame (like the raindata input) but only for selected periods,
  ##                   including a column with event IDs ("st" of every event)
  
  
  
  ## removes periods which start befor the data start
    if ( length( which( periods$st < raindata[1,1] ) )  > 0 ) {
      periods <- periods[ -which( periods$st < raindata[1,1] ), ]
      print("Warning: Some periods of interest start before the data start. Data could not be loaded for these periods.")
    }
  ## creates a vector with time steps for the events and a vector with IDs of the events
    st_en_ev <- data.frame( "st_ev" = as.POSIXct(periods[,1], tz="UTC"), "en_ev" = as.POSIXct(periods[,2], tz="UTC"))
    ev_lengths <- c()  
    step <- raindata[2,1] - raindata[1,1]
    
    # 1st iteration
    help <- seq(st_en_ev[1,1], st_en_ev[1,2], step)  # a vector of time steps for the event
    ev_lengths[1] <- length(help)                    # remembers the length of the event
    TimE <- help
    TimID <- rep( st_en_ev[1,1] , ev_lengths[1])     # creates a vector with rain event's IDs
    
    # other iterations
    if ( length(st_en_ev[,1]) > 1 ) {
      for (i in 2:length(st_en_ev[,1])) {
        help <- seq(st_en_ev[i,1],st_en_ev[i,2], step) # a vector of time steps for the event
        ev_lengths[i] <- length(help)                  # remembers the length of each event     
        TimE <- c(TimE, help)
        TimID <- c(TimID, rep( st_en_ev[i,1] , ev_lengths[i]))
      }
      
      TimE  <- as.POSIXlt( TimE , tz = "UTC"); TimE  <- as.POSIXct( TimE , tz = "UTC")
      TimID <- as.POSIXlt( TimID, tz = "UTC"); TimID <- as.POSIXct( TimID, tz = "UTC")
    }
    
  ## treats read data with irregular timing          WHEN DOES THIS HAPPEN ??
    matchNApos <- which( is.na(match(TimE, raindata[,1])) == T )
    if ( length(matchNApos) > 0 ) {
      for ( ii in 1:length(matchNApos) ) {
        remain <- as.numeric(TimE[matchNApos[ii]]) %% as.numeric(step*60)
        if ( matchNApos[ii]==1  ) {
          TimE[matchNApos[ii]] <- as.POSIXct(TimE[matchNApos[ii]] -remain, origin="1970-01-01 00:00:00 UTC", tz="UTC")
        }else{
          if ( TimE[matchNApos[ii]-1] < TimE[matchNApos[ii]] - remain ) {   # either decreases time to the closest step
            TimE[matchNApos[ii]] <- as.POSIXct(TimE[matchNApos[ii]] -remain, origin="1970-01-01 00:00:00 UTC", tz="UTC")
          }else {
            TimE <- TimE [ -matchNApos[ii] ]    # or deletes the data if the closest lower step already exists
            # SHOULDNT HERE BE TimID AS WELL ??
          }
        }
      }
    }
  
  ## selects only desired (within the observed periods) data     
    match_pos2 <- match(TimE, raindata[,1])
    raindata_events <- raindata[match_pos2,]
    raindata_events$id <- TimID 

  return(raindata_events)
}


########################################################################################################################


match_with_IDs <- function ( rainfall_datfr, IDs  ) {
  ## the function selects rainfall data based on event IDs
  ##
  ## inputs:
  ## rainfall_datfr - a data frame with a column "id"
  ## IDs  - a vector of IDs to be matched
  ##
  ## outputs:
  ## out - a data frame (like the rainfall_datfr) but only rows where the id was matching
  
  
  wrk <- rainfall_datfr
  
  if ( typeof(wrk$id)=="character" && typeof(IDs)!="character" ) {
    wrk$id <- as.POSIXct(wrk$id, tz = "UTC")
  }
  
  if ( typeof(wrk$id)!="character" && typeof(IDs)=="character" ) {
    IDs <- as.POSIXct(IDs, tz = "UTC")
  }
  
  out <- wrk [ which( wrk$id %in% IDs ) ,   ]
  
  return(out)
}


########################################################################################################################


basic_statistics <- function( rainfall_datfr,  periods = NA  ) {
  ## the function creates a table with basic statistics for a rainfall data set
  ##
  ## inputs:
  ##    rainfall_datfr - rainfall data in the format of the output of areal_rainfall function (a data frame)
  ##                     Use data matched with periods only if the periods used for matching are the same as the
  ##    periods  - a data frame with two columns (st and end as POSIXct)
  ##
  ## outputs:
  ##    ev.info - a data frame with basic statistics of the events
  ## 
  
  
  # checks what kinds of data are being evaluated
  noIDs <- TRUE
  if ( "id" %in% names(rainfall_datfr) ) { noIDs <- FALSE }
  
  noPeriods <- FALSE
  if ( length(periods) == 1 ) {  
    if ( is.na(periods) ) { noPeriods <- TRUE }
    else { stop("not valid periods") }
  }
  
  if ( noIDs & noPeriods) { stop("no periods specified and no IDs in the rain data") }
  
  if (noIDs) { # if only periods are available
    rainfall_datfr <- match_with_periods( raindata = rainfall_datfr, periods = periods   ) # adds IDs
    noIDs <- FALSE
  }
  
  IDs    <- unique( rainfall_datfr$id )
  n_rows <- length( IDs )
  
  if (noPeriods) { # if only IDs are available
    st  <- rep(NA, n_rows)
    en  <- rep(NA, n_rows)
    duration <- rep(NA, n_rows)
  }
  if ( !noIDs & !noPeriods) {  # if both are available
    st  <- periods[,1]
    en  <- periods[,2]
    duration <- difftime(en, st, units="mins")
    if ( FALSE %in% (st == IDs)  ) { stop("periods mismatch with IDs") } 
  }
 
  # determines how many rainfall data time series are included in the data frame
  names_time_series <- names(rainfall_datfr) [ which( names(rainfall_datfr) != ( "time") & names(rainfall_datfr) != "id" ) ]
  
  # creats a mask for the output
  stat.names <- c("height", "Rmax", "Rmax10")
  stat.names_all <- c() 
  for ( i in 1 : length(names_time_series) ) {
    stat.names_i   <- paste(names_time_series[i], stat.names, sep="_")
    stat.names_all <- c(stat.names_all, stat.names_i) 
  }
  stat_mask <- as.data.frame( matrix(NA, ncol = length(stat.names_all) , nrow = n_rows) )
  names(stat_mask) <- stat.names_all
  ev.info <- cbind( "id" = IDs , "st" = st , "en" = en , "duration" = duration , stat_mask ) 
  
  # calculates statistics for one or more specified rainfall data series
  for (j in 1 : n_rows) {  # for all events
    ID <- IDs [j]
    datfr <- match_with_IDs( rainfall_datfr = rainfall_datfr , IDs = ID ) 
    
    for ( i in names_time_series ) {    # for all rainfall data time series
      columns <- which( substr(stat.names_all, 1, nchar(i)) == i ) + 4
        
      simple_datfr <- data.frame( time = datfr[,1] , rain = datfr[,i] )
      
      ev.info[ j, columns ] <-  RG.stat.interval_core( RG.i = simple_datfr , nas = T ) [2:4] # basic statistics without duration
    }
  }

  return ( ev.info )
}

# # if the input was a list of data frames, each with a different event periods 
# 
# res.all <- data.frame("duration"=NA, "height"=NA, "Rmax"=NA, "Rmax10"=NA)
# res.i   <- data.frame("st"=NA, "en"=NA, 
#                       "duration"=NA, "height"=NA, "Rmax"=NA, "Rmax10"=NA)
# res.i$st <- as.character(res.i$st); res.i$en <- as.character(res.i$en) 
# ev.info <- cbind(event[[8]][1:length( event[[8]][,1]), ], 
#                  res.all, res.i, res.i, res.i, res.i, res.i, res.i, res.i) # event[[ RGs+1 ]]
# 
# for(j in 1:length(ev.info[,1])){
#   #calculate statistics for areal rainfall
#   ev.info[j, 3:6] <- RG.stat.interval(RG[[8]], event[[8]]$st[j], event[[8]]$en[j], T)  # RGs+1
#   #calculate statistics for point RG rainfalls
#   for(i in 1:7){
#     id <- which(event[[i]][,1] >= event[[8]][j,1] & event[[i]][,2] <= event[[8]][j,2]) # RGs+1
#     if(length(id)==0){next}else{
#       ev.info[j, (i*6+1)] <- as.character(event[[i]]$st[id[1]])  # start of the first one
#       ev.info[j, (i*6+2)] <- as.character(event[[i]]$en[id[length(id)]]) # end of the last one
#       ev.info[j, (i*6+3):(i*6+6)] <-  RG.stat.interval(RG[[i]], as.POSIXct(ev.info[j, (i*6+1)], tz="UTC"), 
#                                                        as.POSIXct(ev.info[j, (i*6+2)], tz="UTC"), T)
#     }
#     if(length(id)>1){print(c(i, "___", id))}
#   }
# }


########################################################################################################################


RG.stat.interval_core <- function(RG.i, nas) {
  ## function to make selected statistics for the given data frame
  ##
  ## Inputs:  RG.i  -  data.frame of two columns (time, R[mm/h])
  ##          nas - logical, na.rm = nas
  ##
  ## Outputs: res - vector with basic rainfall statistics
  
  res <- c("duration"=NA, "height"=NA, "Rmax"=NA, "Rmax10"=NA)
  
  if ( length(RG.i[,2])==0 | ! ( FALSE %in% is.na(RG.i[,2]) )  ) {
    return(res)
  }

  res[1] <- NA                               #duration of rainfall
  res[2] <- sum(RG.i[,2], na.rm=nas)/60      #total height of rainfall [mm]
  res[3] <- max(RG.i[,2], na.rm=nas)         #max rain rate [mm/h]
    
  if (length(RG.i[,1]) < 10) {               #max 10min rain rate [mm/h]
    res[4] <- sum(RG.i[,2], na.rm=nas)/10
  } 
  else {
    RG.i.mtrx <- matrix(c(RG.i[,2], rep(0,10)), nrow=10, ncol= length(RG.i[,2])+9, byrow=T)
    res[4] <- max(apply(RG.i.mtrx, 2, mean, na.rm=nas), na.rm=nas)
  }
  
  
  return(round(res, 1))
}  


########################################################################################################################


read_regular_FGdata <- function(path, metric){
  ## function to read and connect data from a regular FG file (.csv)
  ##
  ## inputs:
  ## path  - path to the FG file
  ##
  ## outputs:
  ## dat - data.frame with 2 columns (time, Q)
  ##                  
  
  
  if (substr(path, nchar(path)-2, nchar(path)) != "csv"){
    stop("wrong path to the file")  #basic checks on input parameters
  }
  
  dat_help<- read.csv(path, header=T, sep=";")
  dat_help[,1] <- as.POSIXct( dat_help[,1] , tz = "UTC")
  
  
  dat <- dat_help[c("time", metric)]
  
  return(dat)
}


########################################################################################################################


basic_statistics_flows <- function( flow_datfr,  periods = NA  ) {
  ## the function creates a table with basic statistics for a discharge data set
  ##
  ## inputs:
  ##    flow_datfr - discharge data (a data frame)
  ##                     Use data matched with periods only if the periods used for matching are the same as the
  ##    periods  - a data frame with two columns (st and end as POSIXct)
  ##
  ## outputs:
  ##    ev.info - a data frame with basic statistics of the events
  ## 
  
  
  # checks what kinds of data are being evaluated
  noIDs <- TRUE
  if ( "id" %in% names(flow_datfr) ) { noIDs <- FALSE }
  
  noPeriods <- FALSE
  if ( length(periods) == 1 ) {  
    if ( is.na(periods) ) { noPeriods <- TRUE }
    else { stop("not valid periods") }
  }
  
  if ( noIDs & noPeriods) { stop("no periods specified and no IDs in the rain data") }
  
  if (noIDs) { # if only periods are available
    flow_datfr <- match_with_periods( raindata = flow_datfr, periods = periods   ) # adds IDs
    noIDs <- FALSE
  }
  
  IDs    <- unique( flow_datfr$id )
  n_rows <- length( IDs )
  
  if (noPeriods) { # if only IDs are available
    st  <- rep(NA, n_rows)
    en  <- rep(NA, n_rows)
    duration <- rep(NA, n_rows)
  }
  if ( !noIDs & !noPeriods) {  # if both are available
    st  <- periods[,1]
    en  <- periods[,2]
    duration <- difftime(en, st, units="mins")
    if ( FALSE %in% (st == IDs)  ) { stop("periods mismatch with IDs") } 
  }
  
  # determines how many flow data time series are included in the data frame
  names_time_series <- names(flow_datfr) [ which( names(flow_datfr) != ( "time") & names(flow_datfr) != "id" ) ]
  
  # creats a mask for the output
  stat.names <- c("Vtot" , "Qmax", "TQmax")
  stat.names_all <- c() 
  for ( i in 1 : length(names_time_series) ) {
    stat.names_i   <- paste(names_time_series[i], stat.names, sep="_")
    stat.names_all <- c(stat.names_all, stat.names_i) 
  }
  stat_mask <- as.data.frame( matrix(NA, ncol = length(stat.names_all) , nrow = n_rows) )
  names(stat_mask) <- stat.names_all
  ev.info <- cbind( "id" = IDs , "st" = st , "en" = en , "duration" = duration , stat_mask )
  
  # calculates statistics for one or more specified flow data series
  for (j in 1 : n_rows) {  # for all events
    ID <- IDs [j]
    datfr <- match_with_IDs( rainfall_datfr = flow_datfr , IDs = ID ) 
    
    for ( i in names_time_series ) {    # for all time series
      columns <- which( substr(stat.names_all, 1, nchar(i)) == i ) + 4
      
      simple_datfr <- data.frame( time = datfr[,1] , flow = datfr[,i] )
      
      ev.info[ j, columns ] <-  FG.stat.interval_core( FG.i = simple_datfr , nas = T ) [2:4] # basic statistics without duration
      
      ev.info[ j, columns[3] ] <- as.character( format.POSIXct( simple_datfr$time[ as.numeric(ev.info[ j, columns[3] ]) ] , 
                                                                format = "%Y-%m-%d %H:%M:%S" ) )
    }
  }
  
  return ( ev.info )
}


########################################################################################################################


FG.stat.interval_core <- function(FG.i, nas) {
  ## function to make selected statistics for the given data frame
  ##
  ## Inputs:  FG.i  -  data.frame of two columns (time, Q)
  ##          nas - logical, na.rm = nas
  ##
  ## Outputs: res - vector with basic discharge statistics
  
  res <- c("duration"=NA, "Vtot" = NA, "Qmax" = NA, "TQmax" = NA)
  
  if ( length(FG.i[,2])==0 | ! ( FALSE %in% is.na(FG.i[,2]) )  ) {
    return(res)
  }
  
  stps <- as.numeric( FG.i$time[2] - FG.i$time[1] )
  
  res[1] <- NA                                           #duration of rainfall
  
  res[2] <- sum(FG.i[,2], na.rm = nas) * 60 * stps       # total runoff in [m^3]
  res[3] <- max(FG.i[,2], na.rm = nas)                   # max discharge in [m^3/s]
  res[4] <- which(FG.i[,2]==res[3]) [1]                  # time step with the max discharge
  
  return(res)
}

