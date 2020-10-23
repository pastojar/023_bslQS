
########################################################################################################################
# runs the SWMM R-R model and remebers the results

group.run.noInf <- function(par, prodata, Rain_File_Tab, RRmodel) {
  
  hlp0 <- list()
  for (j in 1:length(prodata)) {
    for (i in names(prodata[[j]])) {
      hlp0[[ i ]] <- prodata[[j]] [[ i ]]
    }  
  }
  prodata <- hlp0
  
  
  invalid <- c(); results <- list(); Qmod <- list(); Qobs <- list()
  for (i in names(prodata) ) {
    if ( !i %in% row.names(Rain_File_Tab)[which(Rain_File_Tab$Rain_File_Availabs==F)]  ) {
      
      ## runs the model
        hlp <- RRmodel( par = par, L = prodata[[i]][[1]], inp.file = prodata[[i]][[3]] )
        out <- unname( hlp )
      
      ## remembers the time stamps
        timestamp <- sysanal.decode(names(hlp))$val
        time <- prodata[[i]][[2]][,1]
        time <- as.POSIXct(time, format="%d/%m/%Y %H:%M", origin="1970-01-01 00:00:00 UTC", tz="UTC")
      
      ## saves the modelled and observed Q
        results[[i]]  <-  data.frame(timestamp = timestamp, Qmod = out, Qobs = prodata[[i]][[2]][,2], time = time)
    }
    
    if ( i %in% row.names(Rain_File_Tab)[which(Rain_File_Tab$Rain_File_Availabs==F)]  ) {
      
      timestamp <- sysanal.decode( prodata[[i]]$Layout )$val
      time <- prodata[[i]][[2]][,1]
      time <- as.POSIXct(time, format="%d/%m/%Y %H:%M", origin="1970-01-01 00:00:00 UTC", tz="UTC")
      
      results[[i]] <- data.frame(timestamp = timestamp, 
                                 Qmod = NA[], 
                                 Qobs = prodata[[i]][[2]][,2], 
                                 time = time)
    }
  }
  
  # puts the results in a dataframe
  group.Qres <- data.frame()
  event_IDs <- names(results)
  for ( i1 in 1:length(event_IDs) ) {

    hlp <- data.frame( time = results[[i1]]$time, timestamp = results[[i1]]$timestamp, 
                       id = event_IDs[i1], Qobs = results[[i1]]$Qobs , Qmod = results[[i1]]$Qmod, stringsAsFactors = F)
    
    group.Qres <- rbind(group.Qres, hlp)
  }

  
  return( group.Qres )
}


########################################################################################################################

group.plot.noInf <- function(group.res, rain.data, out.dir, name) {
  
  pdf( paste(out.dir, "/hydrographs__", name, ".pdf", sep="") )
  
  for (EventID in names(group.res)) {
    if ( !is.na(group.res[[EventID]]) ) {
      
      Rain <- match_with_IDs(rain.data, as.POSIXct(EventID, tz="UTC"))
      
      timestamp_Rain_1st <- kimisc::hms.to.seconds(format(Rain$time[1], format="%H:%M:%S")) /3600 #[h]
      dt <- min(abs( ( kimisc::hms.to.seconds(format(Rain$time[-1],          format="%H:%M:%S")) /3600 ) -
                     ( kimisc::hms.to.seconds(format(Rain$time[-nrow(Rain)], format="%H:%M:%S")) /3600 ) ) )
      timestamp_Rain <- seq( from = timestamp_Rain_1st, by =dt , length.out = length(Rain$time) ) 
      timestamp_Q    <- group.res[[EventID]]$timestamp
      timestamp      <- sort(union(timestamp_Rain, timestamp_Q))
      
      timestamp_Rain <- timestamp_Rain - timestamp[1]
      timestamp_Q    <- timestamp_Q - timestamp[1]
      timestamp      <- timestamp - timestamp[1]
      
      Qmod <- group.res[[EventID]]$Qmod
      Qobs <- group.res[[EventID]]$Qobs
      
      ymax <-  1.37 * max( max(Qmod), max(Qobs) )
      par(mar = c(3.2, 3.2, 1, 3.2))
      plot(x= NA, y = NA, xlim = range(timestamp), ylim = c(0, ymax), ylab = NA, xlab = NA, type = "n")
        mtext(side = 2, line = 2, "Discharge [l/s]")
        mtext(side = 1, line = 2, "Time [h]")  
      
        # observed
        points(x = timestamp_Q, y = Qobs)
        # modelled
        lines(x = timestamp_Q, y = Qmod)
        
      # rain
      Rain_toPrint <- apply( Rain[ names(Rain)[ !(names(Rain) %in% c("time", "id")) ] ], 1, mean )
      par(new = T)
      plot(x = NA, y = NA, xlim = range(timestamp), ylim = c(4*max(Rain_toPrint), min(Rain_toPrint)),   # only RG1 now !!!!!!!!!!!!!!!!!
             axes = F, xlab = NA, ylab = NA, type = "n")
        lines(x = timestamp_Rain, y = Rain_toPrint, col = "blue", pch = 20)
        axis(side = 4, labels = round(seq(from = 0, to = 1.2*max(Rain_toPrint), length.out = 4 ), digits = 0) , 
                           at = round(seq(from = 0, to = 1.2*max(Rain_toPrint), length.out = 4 ), digits = 0) ) 
        mtext(side = 4, line = 2, "Precipitation [mm/h]", adj = 1)
        
        # legend("topright", legend = c("modelled", "observed"), pch = c(NA, 1), lty = c(1, NA), horiz = T)
        # legend("topleft", legend = paste( "Event: ", EventID , sep = "" ) )
     
    }
    
    else {
      plot(x =50, y = 50, type = "n")
      legend("topleft", legend = paste("Event: ", EventID, "  - rainfall data NA", sep=""))
    }
  }
  
  dev.off() 
}


########################################################################################################################

sup.group.plot.noInf <- function(mod.scens.to.plot, name, sup.group.res, newRain, out.dir) {
  
  # manages the number of scenarios to be plotted
  if ( length(mod.scens.to.plot) > 4) { stop("Too much data sets.") }
  
  if ( length( grep("ThPol3", mod.scens.to.plot) ) > 0 ) { 
    mod.scens.to.plot_rain  <- paste(  c("RG1", "RG2", "RG3"), mod.scens.to.plot , sep="_-_" ) 
  } else {
    mod.scens.to.plot_rain <- mod.scens.to.plot
  }
  newRain <- newRain[ c("time", "id", mod.scens.to.plot_rain) ]
  
  sup.group.res <- sup.group.res[ c("timestamp", "id", "Qobs", "sd_Qobs", mod.scens.to.plot) ]
  
  
  
  pdf( paste0(out.dir, "/hydro+hyeto2__", name, ".pdf"),
       pointsize = 8, paper = "a4", height = 11.69 , width =  8.27 )
  
  layout(mat = matrix(c(1,1,2,3,4,5,6,7,8,9), 5, 2, byrow = T), heights = c(0.8,2,2,2,2) )
  
  
  # defines colors for scenarios
  scen_colors <- data.frame( scen_name = mod.scens.to.plot_rain, rain = NA, Q = NA )
  scen_colors$rain[ grep( "_RG",  scen_colors$scen_name  ) ] <- c( rgb(red = 0.99, green = 0, blue = 0, alpha = 0.6), 
                                                                   rgb(red = 1, green = 0.6, blue = 0, alpha = 0.6), 
                                                                   rgb(red = 1, green = 0, blue = 0.6, alpha = 0.6) )
  scen_colors$Q[ grep( "_RG",  scen_colors$scen_name  ) ] <- rgb(red = 0.99, green = 0, blue = 0, alpha = 0.5) 
  
  scen_colors$rain[ !grepl( "_RG",  scen_colors$scen_name  ) ] <- c( rgb(red = 0, green = 1, blue = 0, alpha = 0.6), 
                                                                     rgb(red = 0, green = 0, blue = 1, alpha = 0.6), 
                                                                     rgb(red = 0, green = 0.6, blue = 0.7, alpha = 0.6),
                                                                     rgb(red = 0.4, green = 0.2, blue = 0.5, alpha = 0.6))
  scen_colors$Q[ !grepl( "_RG",  scen_colors$scen_name  ) ] <- scen_colors$rain[  !grepl( "_RG",  scen_colors$scen_name  ) ]
  
  
  i_plot <- 0
  EventIDs <- unique(as.character(newRain$id))
  for ( EventID in as.character( EventIDs[ order( uni.data$RG.overview$meanRain_Rmax10[ uni.data$RG.overview$id %in% EventIDs ] ) ] ) ) {
    i_plot <- i_plot + 1

    # selects data for the event
    Rain.all <- match_with_IDs( newRain, as.POSIXct(EventID, tz="UTC") )
    Q.all    <- match_with_IDs( sup.group.res, EventID )
    Q.all    <- Q.all[ - which( is.na(Q.all$Qobs) ),  ]
    Qobs     <- Q.all$Qobs
    sd_Qobs  <- Q.all$sd_Qobs
    
    # defines y axis range
    Rain.all_pure <- Rain.all[] [ !(names(Rain.all) %in% c("time", "id")) ]
    ylim_R <- c(4*max(Rain.all_pure, na.rm = T), min(Rain.all_pure, na.rm = T))
    Q.all_pure <- Q.all[ , !(names(Q.all) %in% c("time", "timestamp", "id", "sd_Qobs")) ]
    Q.all_pure$Qobs_plus2sd <- Q.all_pure$Qobs + 2*Q.all$sd_Qobs
    ylim_Q <- c(0, 1.37 * max( Q.all_pure, na.rm = T ))  
    
    # defines x axis range
    timestamp_Rain_1st <- kimisc::hms.to.seconds(format(Rain.all$time[1], format="%H:%M:%S")) 
    dt <- min(abs( kimisc::hms.to.seconds(format(Rain.all$time[-1], format="%H:%M:%S"))   -
                   kimisc::hms.to.seconds(format(Rain.all$time[-nrow(Rain.all)], format="%H:%M:%S")) ) )
    timestamp_Rain <- seq( from = timestamp_Rain_1st, by = dt , length.out = length(Rain.all$time) ) 
    
    timestamp_Q    <- round( Q.all$timestamp *60*60 )
    timestamp      <- sort(union(timestamp_Rain, timestamp_Q))
    
    timestamp_Rain <- timestamp_Rain - timestamp[1]
    timestamp_Q    <- timestamp_Q - timestamp[1]
    timestamp      <- timestamp - timestamp[1]
    
    
    # plotting hydrographs
    
    # legend
    if ( i_plot %% 8 == 1  ) {
      par(mar = c(0,0,0,0))
      plot.new()
      legend("center", legend = c( mod.scens.to.plot_rain , "Q observed", "95% confidence interval of Q observed" ), 
             ncol = 2, lty = c( rep(1, length(mod.scens.to.plot_rain) + 1 ), 0 ),
             col = c(scen_colors$rain,
                     rgb(red = 0.2, green = 0.2, blue = 0.2, alpha = 0.99),
                     rgb(red = 0.6, green = 0.6, blue = 0.6, alpha = 0.5)  ),
             lwd = c( rep(0.7, length(mod.scens.to.plot_rain) + 1 ), NA),
             pch = c( rep(NA, length(mod.scens.to.plot_rain) + 1  ), 15 ), pt.cex = 2  )
    }
    
    # Q observed
    par(mar = c(3.2, 3.2, 2.5, 3.2))
    # plot(x= NA, y = NA, xlim = range(timestamp), ylim = ylim_Q, 
    #      ylab = NA, xlab = NA, type = "n", cex.axis = 0.9)
    plot(x= NA, y = NA, xlim = range(timestamp), ylim = ylim_Q, 
         ylab = NA, xlab = NA, type = "n", xaxt = "n"  )
    timestamp_hrs <- timestamp / (60*60)
    axis(side = 1, labels = seq( 0, max(timestamp_hrs), by = 2 ), 
         at = seq( 0, max(timestamp), by = 2*60*60 ) )
    lines(x = timestamp_Q, y = Qobs, col = rgb(red = 0.2, green = 0.2, blue = 0.2, alpha = 0.99), lwd = 0.7)
    
    # Q observed  +/- 2*sd
    Q_plus_2sd  <- Qobs + 2*sd_Qobs
    Q_minus_2sd <- Qobs - 2*sd_Qobs 
    polygon( c(timestamp_Q,rev(timestamp_Q)),
             c( Q_minus_2sd , rev(Q_plus_2sd) ), 
             col = rgb(red = 0.6, green = 0.6, blue = 0.6, alpha = 0.5), border=NA ) 
    
    # iterates data sets / scenarios
    i_rain_scen <- 0
    for ( i_datasets_rain in mod.scens.to.plot_rain ) {
      i_rain_scen <- i_rain_scen + 1
      
      # selects the specified data set
      Rain <- Rain.all[ c("time", "id", i_datasets_rain) ]
      Rain_to_print <- Rain[ , !(names(Rain) %in% c("time", "id")) ]
      
      if ( grepl( "ThPol3", i_datasets_rain ) ) {
        Qmod <- Q.all[ , grep( "ThPol3", names(Q.all) ) ]
      } else {
        Qmod <- Q.all[ , i_datasets_rain ]
      }
      
      # rain
      par(new = T)
      plot(x = NA, y = NA, xlim = range(timestamp_Rain), ylim = ylim_R,
           axes = F, xlab = NA, ylab = NA, type = "n")
      
      if ( i_rain_scen == 1 ) {
        title( paste0( EventID, ",  Rmax10 = ",  uni.data$RG.overview$meanRain_Rmax10[ uni.data$RG.overview$id %in% EventIDs ] [match(EventID, as.character(EventIDs))], " mm/h,  ") ) 
        
        axis(side = 4, labels = round(seq(from = 0, to = 1.2*max(Rain.all_pure, na.rm = T), length.out = 4 ), digits = 0) , 
             at = round(seq(from = 0, to = 1.2*max(Rain.all_pure, na.rm = T), length.out = 4 ), digits = 0) ) 
        mtext(side = 4, line = 2, "Precipitation [mm/h]", adj = 1, cex = 0.8)
      }
      
      lines(x = timestamp_Rain, y = Rain_to_print,  col = scen_colors$rain[ which( scen_colors$scen_name == i_datasets_rain ) ], lwd = 0.5)
     
      # Q modelled
      par(new = T)
      plot(x= NA, y = NA, xlim = range(timestamp), ylim = ylim_Q, 
           axes = F, ylab = NA, xlab = NA, type = "n")
      if ( i_rain_scen == 1 ) {
        mtext(side = 2, line = 2, "Discharge [l/s]", cex = 0.8)
        mtext(side = 1, line = 2, "Time [h]", cex = 0.8)
      }

      lines(x = timestamp_Q, y = Qmod,  col = scen_colors$Q[ which( scen_colors$scen_name == i_datasets_rain ) ], lwd = 0.7)
      
    }
  }
  
  dev.off()
  
  
}




