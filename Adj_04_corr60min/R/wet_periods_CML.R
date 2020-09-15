
########################################################################################################################

wet_periods_CML <- function ( rain.data ) {
  # identifies wet periods from the mean of all time series
  

  R <- rain.data
  
  # calculates the mean          # instad of the mean, it is possible to use a single time series (e.g. a single RG)
  mean.rain <- R[ , c("time", colnames(R)[ !colnames(R) %in% c("time", "id") ] [1] ) ]
  names(mean.rain)[2] <- "R_mean"
  mean.rain$R_mean <- NA[]
  mean.rain$R_mean <- apply(R[, -1], 1, mean, na.rm = TRUE)  
  dat  <- mean.rain
  
  
  # standard deviation of TPL calculated for moving
  # window of "width" minutes aligned to the center
  sigma <- dat; sigma[!names(dat) %in% "time"] <- NA[]
  width <- 60
  for( i_timeser in names(dat)[ !names(dat) %in% "time"] ) {
    
    hlp <- zoo::zoo( x = dat[[i_timeser]] , order.by = dat[["time"]] )
    sd_tpl <- zoo::rollapply(hlp, width = width, sd, by = 1, align = 'center' )
    
    first_steps <- head(hlp, n = width/2 - 1); first_steps[] <- NA
    last_steps  <- tail(hlp, n = width/2 ); last_steps[]  <- NA
    
    sigma[i_timeser] <- zoo::coredata( rbind(first_steps, sd_tpl, last_steps) )
  }
  
  
  # -----------------------------------------------------------------------------------------------------------------------------
  # identify strange periods (sigma = 0 for longer than 6 hours (1/4 day))
  
  iszero <- function(x) {
    
    zero <- F
    if ( length(which(is.na(x))) == length(x) ) {
      zero <- NA
    } else {
      zero <- sum(x, na.rm = T) == 0
    }
    
    return (zero)
  }
  
  identify_no_fluctuations <- function (sigma, width = 60 * 6, by = 60 * 6) {
    # Function to identify longer periods with no fluctuations
    # sigma - floating window standard deviation appied to TPL
    # width - windows size  (elements)
    
    nosignal <- zoo::rollapply( data = sigma, width = width, FUN = iszero, by = by, align = 'right' )
    
    nosignal2 <- sigma; nosignal2[] <- NA
    
    nosignal2[ zoo::index(nosignal) ] <- nosignal
    nosignal2 <- zoo::na.approx(nosignal2, method = 'constant', f = 1, rule = 2)
    
    return(nosignal2)
    
  }
  
  replace_no_signal <- function (sigma, width = 60 * 6, by = 60 * 6) {
    # Function to identify longer periods with no fluctuations and set them
    # to NAs
    # sigma - floating window standard deviation appied to TPL
    # width - windows size  (elements)
    # by - step by which is window moved [elements]
    
    nosignal <- identify_no_fluctuations( sigma = sigma  )
    sigma[nosignal == 1] <- NA
    
    return(sigma)
    
  }
  # -----------------------------------------------------------------------------------------------------------------------------
  
  
  # classify periods with the highest sigma (corresponding to the specified quantile) as wet
  first <- TRUE
  quant <- 0.9
  for ( i_timeser in names(sigma)[ !names(sigma) %in% "time"] ) {
    
    hlp <- zoo::zoo( x = sigma[[i_timeser]] , order.by = sigma[["time"]] )
    hlp <- replace_no_signal( sigma = hlp ) # sets longer periods with no fluctuations to NA
    
    wet <- hlp[ hlp > quantile(hlp, quant, na.rm = T) ]
    
    if (first == TRUE) {
      first <- FALSE
      wet_mtx <- data.frame( zoo::index(wet) , zoo::coredata(wet) )
      colnames(wet_mtx) <- c("time", i_timeser)
    } else {
      wet_mtx <- cbind(wet_mtx, zoo::coredata(wet))
      colnames(wet_mtx)[ which( colnames(wet_mtx)== "zoo::coredata(wet)" ) ] <- i_timeser
    }
  }
  
  # -----------------------------------------------------------------------------------------------------------------------------
  rainy.periods <- function(tim, gap.lim){
    ## function to identify dry and rainy periods using raw data
    ## (time stamps of tips from rain gauges)
    ##
    ## Inputs:   tim     - vector of time stamps of RG tips
    ##           gap.lim - length of gap between two tips [minutes]
    ##                     which idicates end of rainfall
    ## Outputs:  data frame with begginings and ends of rainfall periods
    
    gap <- difftime(tim[-1], tim[-length(tim)], units="mins")
    end.ind <- c(which(gap > gap.lim), length(tim))
    beg.ind <- c(1, which(gap > gap.lim) + 1)
    
    return(data.frame("st" =tim[beg.ind], "en"=tim[end.ind]))
  }
  # -----------------------------------------------------------------------------------------------------------------------------
  
  # define starts and ends of wet periods
  events <- rainy.periods(tim = wet_mtx$time, gap.lim = 60) # gap.lim - length of max. gap between two tips [min]
  events$en <- events$en + 60*60  # adds one hour to the end of the periods; to make sure it is dry at the event end
  
  # rounds to minutes and add one minute to the beggining and end of rainfall
  t.num1 <- floor(as.numeric(events[,1])/60)*60 
  t.num2 <- ceiling(as.numeric(events[,2])/60)*60 
  tim1 <-  as.POSIXct(t.num1, origin="1970-01-01 00:00:00 UTC", tz="UTC") - 60
  tim2 <-  as.POSIXct(t.num2, origin="1970-01-01 00:00:00 UTC", tz="UTC") + 60
  
  events [,1] <- tim1
  events [,2] <- tim2

  
  return( events )
}


