######################################
# J. Pastorek, JAN 2018
######################################

#---------------------------------------------------------------------

group.run.plot <- function(dir, name, par, ev.data) {
  hlp <- list()
  for (i in 1 : length(ev.data$Ca)) {
    hlp[[i]] <- ev.data$Ca[[i]]
  }
  for (i in ( length(ev.data$Ca) + (1 : length(ev.data$Pre)) ) ) {
    hlp[[i]] <- ev.data$Pre[[i-length(ev.data$Ca)]]
  }
  ev.data <- hlp
  
  invalid <- c()
  pdf( paste(dir, "/hydrographs_", name, ".pdf", sep="") )
  for (i in 1 : length(ev.data)) {
    if (!is.null(ev.data[[i]]) && ev.data[[i]][[1]]!=42) {
      hlp <-  model.swmm(par = par, inp.file = ev.data[[i]][[3]], 
                         out.data = ev.data[[i]][[2]], L = ev.data[[i]][[1]])
      timestep <- sysanal.decode(names(hlp))$val
      out <- unname(hlp)
      plot(x = c(min(timestep), max(timestep)), y = c(0, max(max(out), max(ev.data[[i]][[2]][,2]))), 
           ylab = "Discharge [l/s]", xlab = "Timestep [h]", type = "n")
      # observed
      points(x = timestep, y = ev.data[[i]][[2]][,2])
      # modelled
      lines(x = timestep, y = out)
      
      legend("topright", legend = c("observed", "modelled"), pch = c(1,NA), lty = c(NA, 1))
      legend("topleft", legend = paste( "Event: ", substr(ev.data[[i]][[3]], nchar(ev.data[[i]][[3]])-22, 
                                                          nchar(ev.data[[i]][[3]])-4), sep="" ) )
    }
    if (!is.null(ev.data[[i]]) && ev.data[[i]][[1]]==42) {
      plot(x =50, y = 50, type = "n")
      legend("topright", legend = paste("event #", i, ": setupSWMM.RG: error, missing RG data", sep=""))
      invalid <- c(invalid, i)
    }
  }
  if (length(invalid) > 0) {print(paste("missing RG data: ", invalid))}
  dev.off()
}


#---------------------------------------------------------------------

plot.chains.margs <- function(RAM, pr.dis, which_samples, pack.dir, plot_name) {
  
  RAM_to_plot <- RAM
  RAM_to_plot$samples <- cbind(RAM$samples[ which_samples , ], log.post = RAM_to_plot$log.p[ which_samples ])
  RAM_to_plot <- adaptMCMC::convert.to.coda(RAM_to_plot)

  pdf( paste0(pack.dir, "/", plot_name, "_chains+margs.pdf") )
  plot(RAM_to_plot)
  dev.off()
  
  pdf( paste0(pack.dir, "/", plot_name, "_chains.pdf") )
  coda::cumuplot(RAM_to_plot)
  dev.off()
  
  pdf( paste0(pack.dir, "/", plot_name, "_margs+prior.pdf") )
  sysanal.plot.margs.JA(postsamp=RAM_to_plot, pridist=pr.dis)
  dev.off()
  
  # pdf( paste(pack.dir, "/3_MCMC.nonAdapt_chains+margs.pdf", sep="") )
  # plot(RAM_nonAdapt)
  # dev.off()
  # 
  # pdf( paste(pack.dir, "/3_MCMC.nonAdapt_chains.pdf", sep="") )
  # coda::cumuplot(RAM_nonAdapt)
  # dev.off()
  #   
  # pdf( paste(pack.dir, "/3_MCMC.nonAdapt_margs+prior.pdf", sep="") )
  # sysanal.plot.margs.JA(postsamp=RAM_nonAdapt, pridist=pr.dis)
  # dev.off()
  
  return(TRUE)
}


#---------------------------------------------------------------------

CaPre.predict <- function(L1, L2, y.obs, inp.file, eventID, model, par_samples_Pre, par.tr, par.fix) {
  
  if (!is.na(L1)) { L <- L1 }
  if (!is.na(L2)) { L <- L2 }
  
  input    <- inp.file
  
  ret <- sysanal.predict.bias.OU.JA (
    parsamp.L1        = par_samples_Pre,
    model             = model,
    eventID           = eventID,
    #inp.file          = input,
    #out.data          = out.data,
    inp               = rep(0,length(L)),
    ppt               = rep(0,length(L)),
    L1                = L1,
    #y.obs             = out.data[2:nrow(out.data),2], # WHYYYYY to leave out the first one???
    y.obs             = y.obs,
    L2                = L2,
    predict.bias.cond = sysanal.predict.inp.bias.L1.JA,
    par.tr            = par.tr,
    par.fix           = par.fix )
  
  return(ret)
}


#---------------------------------------------------------------------

CaPre.bTr <- function(transf, res.L) {
  par.tr <- transf$par.tr
  
  if (transf$transf == "LogSinh") {
    inv.fcia <- sysanal.logsinh.inv
    Par.tr <- c(par.tr["alpha"], par.tr["beta"])
  }
  if (transf$transf == "BC") {
    inv.fcia <- sysanal.boxcox.inv
    Par.tr <- c(par.tr["l1"], par.tr["l2"])
  }
  
  ret <- list()
  for ( i_names in names(res.L) ) {
    if ( length(res.L[[i_names]]) == 0 ) { next }
    if ( is.na(res.L[[i_names]]) ) { next }
    ret[[i_names]] <- inv.fcia( res.L[[i_names]], Par.tr[1], Par.tr[2] )
  }
  
  return(ret)
}


#---------------------------------------------------------------------

CaPre.VerInd.Pre <- function(data_obs, data_mod) {
  
  #o=data_obs[2:nrow(data_obs),2] # WHYYYYY to leave out the first one???
  o <- data_obs[1:nrow(data_obs),2]
  t.pre <- seq(1,length(o),by = 1)
 
  up <- data_mod[(row.names(data_mod)=="0.95"),]
  lo <- data_mod[(row.names(data_mod)=="0.05"),]
  ob <- o
 
  
  red_points <- red_points( obs = ob, up = up, lo = lo )

  
  if ( length(which(is.na(ob))) > 0 ) { 
    avail.data <- ob[ - which(is.na(ob)) ] 
  } else {
    avail.data <- ob 
  }
  
  reliab <-  ( 1 - length(which(!is.na(red_points))) / length(avail.data) ) * 100  # prediction reliability [%]
  
  ABW <- mean(up - lo)                        # Average Band Width
  normABW <- mean( (up - lo) / mean(ob) )     # Average Band Width normalized
  
  MIS   <- MIS  (low = lo, upp = up, obs = ob)          # Mean of Interval Scores
  MISS  <- MISS (low = lo, upp = up, obs = ob)          
  NMISS <- nMISS(low = lo, upp = up, obs = ob)
 
  ret <- c( ABW = ABW, normABW = normABW, rlb.prcnt = reliab, MIS = MIS, MISS = MISS, NMISS = NMISS )
  return(ret)
}


#---------------------------------------------------------------------

red_points <- function( obs, up, lo ) {
  avail.data <- which(!is.na(obs))
  red_points <- rep(NA, length(obs)) # data outside of the predicted band
  for (i in avail.data) {
    if ( (up[i] < obs[i])
         ||
         (lo[i] > obs[i])
    ) {
      red_points[i] <- obs[i]  
    }
  }
  
  return(red_points)
}

#---------------------------------------------------------------------

# compute mean of the interval scores as defined by Gneiting and Raftery (2007) and used by Breinholt et al. (2012), 
MIS <-function ( low, upp, obs, conf = 0.1 ) {  #  0.1 ~ 90%
  sh_vec <- upp-low   #sharpness
  under <- (low-obs)*(low>obs)
  over  <- (obs-upp)*(obs>upp)
  
  scores <- sh_vec + 2/conf*(under+over)
  
  MIS <- mean( scores , na.rm = TRUE)  
  return(MIS)
}

# MISS - Mean Interval Skill Score 
MISS <-function ( low, upp, obs, conf = 0.1 ) {  #  0.1 ~ 90%   
  
  # standard MIS
  MIS_stand <- MIS( low = low, upp = upp, obs = obs, conf = conf )
  
  # reference MIS where the 90% (95%, 80%...) range of observations is used instead of the predicted intervals
  MIS_ref <- MIS( low = rep( quantile( obs, conf/2 ),   length(obs) ) , 
                  upp = rep( quantile( obs, 1-conf/2 ), length(obs) ) , 
                  obs = obs, conf = conf )
  
  MISS <- 1 - MIS_stand / MIS_ref   
  return(MISS)
}

# NMISS - normalized Mean Interval Skill Score
nMISS <-function ( low, upp, obs, conf = 0.1 ) {  #  0.1 ~ 90%   
  
  MISS <- MISS( low = low, upp = upp, obs = obs, conf = conf )
  
  nMISS <- 1 / ( 2 - MISS )
  
  return(nMISS)
}

#---------------------------------------------------------------------

# calculates statistics for prediction events
stats_inf <- function( data_mod,  data_obs ) {
  
  stats_it   <- list()
  stats_qntl <- list()
  stats_band <- list()
  
  for ( i_ev in 1:length(data_mod) ) {      #  for each event
    
    Qdata <- data.frame( Qobs      = data_obs[[i_ev]]$Q_Data[,2] ,
                         timestamp = sysanal.decode( L = data_obs[[i_ev]]$Layout )$val )
    
    # statistics for each iteration
    data_mod_ev <- data.frame( data_mod[[i_ev]]$Y.L2.samp )

    for ( i_it in 1:nrow(data_mod_ev) ) {
      Qdata$Qmod <- as.numeric(data_mod_ev[i_it,])

      stat_ev_it <- simple.stats.core( Qdata = Qdata )

      if ( i_it == 1 ) {
        stat_ev <- data.frame( t(stat_ev_it) )
      } else {
        stat_ev <- rbind(stat_ev, stat_ev_it)
      }
    }

    stats_it[[ names(data_obs)[i_ev] ]] <- stat_ev
    
    
    # statistics for prediction band quantiles
    data_mod_ev_qntls <- apply( X = data_mod_ev, FUN = quantile, MARGIN = 2, probs = seq(0.005, 0.995, 0.005) )
    rownames( data_mod_ev_qntls ) <- seq(0.005, 0.995, 0.005)
    
    for ( i_qntl in 1:nrow(data_mod_ev_qntls) ) {
      Qdata$Qmod <- as.numeric(data_mod_ev_qntls[i_qntl,])
      
      stat_ev_qntl <- simple.stats.core( Qdata = Qdata )
      
      if ( i_qntl == 1 ) {
        stat_ev <- data.frame( t(stat_ev_qntl) )
      } else {
        stat_ev <- rbind(stat_ev, stat_ev_qntl)
      }
    }
    rownames( stat_ev ) <- rownames( data_mod_ev_qntls )
    
    stats_qntl[[ names(data_obs)[i_ev] ]] <- stat_ev
    
    
    # statistics for prediction bands of each event - Verification Indicies (ABW, reliab, MIS...)
    stats_band_ev <- CaPre.VerInd.Pre( data_obs = data_obs[[i_ev]]$Q_Data, data_mod = data_mod[[i_ev]]$Y.L2.quant )
    
    if ( i_ev == 1 ) {
      stats_band <- data.frame( t(stats_band_ev) )
    } else {
      stats_band <- rbind(stats_band, stats_band_ev)
    }
    rownames( stats_band )[i_ev] <- names( data_obs )[i_ev]
    
  }
 
  return( list( stats_it = stats_it, stats_qntl = stats_qntl, stats_band = stats_band ) )    
} 


#---------------------------------------------------------------------

# plots basic hydrographs
CaPre.plot.new <- function( data_obs, data_mod, eventSet ) {
  
  out.data_obs <- data_obs[[2]]
  
  if ( eventSet == "Ca"  ) { L <- "L1" }
  if ( eventSet == "Pre" ) { L <- "L2" }
  
  y.quant   <- data_mod[[ paste0("y.", L, ".quant") ]]
  y_B.quant <- data_mod[[ paste0("yplusB.", L, ".quant") ]]
  Y.quant   <- data_mod[[ paste0("Y.", L, ".quant") ]]
  timestep  <- sysanal.decode( colnames(y.quant) )[,2]
  
  
  plot( x = timestep, 
        y = out.data_obs[1:nrow(out.data_obs),2], 
        bty = "o", xaxt = "n",
        ylab = "Discharge [l/s]", xlab = "Timestep [h]", 
        ylim = c( 0 , max( c(Y.quant[(row.names(Y.quant)=="0.95")], out.data_obs[1:nrow(out.data_obs),2] ) ) ) )

  polygon( c( timestep,rev( timestep)),
           c(Y.quant[(row.names(Y.quant)=="0.05")],
             rev(Y.quant[(row.names(Y.quant)=="0.95")])), #set the limits (1st and last quantiles)
           col=gray(0.1), 
           border=NA )
  
  polygon( c( timestep,rev( timestep)),
           c(y_B.quant[(row.names(y_B.quant)=="0.05")],
             rev(y_B.quant[(row.names(y_B.quant)=="0.95")])), #set the limits (1st and last quantiles)
           col=gray(0.5), 
           border=NA )
  
  polygon( c( timestep,rev( timestep)),
           c(y.quant[(row.names(y.quant)=="0.05")],
             rev(y.quant[(row.names(y.quant)=="0.95")])), #set the limits (1st and last quantiles)
           # col=gray(0.8),
           col = "magenta",
           border=NA )
  
  lines(x =  timestep, y = Y.quant[(row.names(Y.quant)=="0.5")], lty = "66")
  
  points(x =  timestep, y = out.data_obs[1:nrow(out.data_obs),2], col="blue")
  
  VerInd <- CaPre.VerInd.Pre( data_obs = data_obs$Q_Data, data_mod = Y.quant )
  points(x =  timestep, y = VerInd$red_points, col="red")

  
}


#---------------------------------------------------------------------

# plots hydrographs with statistics overview
plot_hydro_stats <- function( data_obs, data_mod, stats, eventSet, out_dir ) {
 
  
  # OVERVIEW OF THE PERFORMANCE STATISTICS FOR ALL EVENTS TOGETHER
  
  for ( i_ev in 1:length(stats$stats_it) ) {
  
    if ( i_ev == 1 ) {
      stats_it_mergedEv   <- stats$stats_it[[i_ev]]
      stats_qntl_mergedEv <- data.frame( stats$stats_qntl[[i_ev]], qntl = rownames(stats$stats_qntl[[i_ev]]), stringsAsFactors = F )
    } else {
      stats_it_mergedEv   <- rbind( stats_it_mergedEv,   stats$stats_it[[i_ev]] )
      stats_qntl_mergedEv <- rbind( stats_qntl_mergedEv, data.frame( stats$stats_qntl[[i_ev]], qntl = rownames(stats$stats_qntl[[i_ev]]), stringsAsFactors = F ) )
    }
    
  }

  png( paste0(out_dir, "/7_stats_overview.png") ,
       type="cairo", units = "in", width = 6*4, height = 2, res = 150 )

    par( mar = c(2, 1, 1, 1), mfrow = c(1, 6), cex = 1.2 )
  
    for ( i_stat in c("dV", "dQmax", "NNSE", "SCC") ) {
      
      vioplot::vioplot( stats_it_mergedEv[, i_stat] ,
                        ylab = NA ,
                        plotCentre = "line", col = gray(0.7), range = 0, horizontal = T )
      mtext(side = 3, line = 0, text = i_stat, cex = 1.2)
      
      # lightens the distribution extremes
      if ( i_stat %in% c("dV", "dQmax") ) {
        polygon( y = c(0.5, 0.5, 1.5, 1.5),
                 x = c( quantile(stats_it_mergedEv[, i_stat], c(0, 0.05)) ,
                        rev( quantile(stats_it_mergedEv[, i_stat], c(0, 0.05)) ) ),
                 col = rgb(1,1,1,0.6), border = rgb(1,1,1,0.6) )
        polygon( y = c(0.5, 0.5, 1.5, 1.5),
                 x = c( quantile(stats_it_mergedEv[, i_stat], c(0.95, 1)) ,
                        rev( quantile(stats_it_mergedEv[, i_stat], c(0.95, 1)) ) ),
                 col = rgb(1,1,1,0.6), border = rgb(1,1,1,0.6) )
      } else {
        polygon( y = c(0.5, 0.5, 1.5, 1.5),
                 x = c( quantile(stats_it_mergedEv[, i_stat], c(0, 0.1)) ,
                        rev( quantile(stats_it_mergedEv[, i_stat], c(0, 0.1)) ) ),
                 col = rgb(1,1,1,0.6), border = rgb(1,1,1,0.6) )
      }
      
      # lines for median predictions
      for ( i_ev in 1:length( stats_qntl_mergedEv[ stats_qntl_mergedEv$qntl == "0.5" , i_stat])  ) {
        lines( y = c(0.9, 1.1) , 
               x = rep( stats_qntl_mergedEv[ stats_qntl_mergedEv$qntl == "0.5" , i_stat] [i_ev] , 2) ,
               col = "purple" )    
      }
      
    }
    
    labels <- c( "reliab.", "NMISS" ); names(labels) <- c( "rlb.prcnt", "NMISS" )
    for ( i_stat in names(labels)  ) {
      
      if ( i_stat == "NMISS"   )   { ylim <- range(stats$stats_band[, i_stat]) }
      if ( i_stat == "rlb.prcnt" ) { stats$stats_band[, i_stat] <- stats$stats_band[, i_stat] /100
                                     ylim <- range(stats$stats_band[, i_stat]) }
      
      boxplot( stats$stats_band[, i_stat], range = 0, ylim = ylim,
               col = gray(0.7), horizontal = T, ylab = NA )
      mtext(side = 3, line = 0, text = labels[i_stat], cex = 1.2)
    }
  
  dev.off()
  
  
  # HYDROGRAPHS FOR INDIVIDUAL EVENTS WITH PERFORMANCE STATISTICS
  
  if ( eventSet == "Ca"  ) { L <- "L1" }
  if ( eventSet == "Pre" ) { L <- "L2" }
  
  pdf( paste0(out_dir, "/6_hydro+stats.pdf"),
       pointsize = 8, paper = "a4", height = 11.69 , width =  8.27, )
  
  layout( mat = matrix( c( 1,1,1,     1,1,1,
                           2,2,2,     9,9,9,
                           3,4,5,     10,11,12,
                           6,7,8,     13,14,15,
                           16,16,16,  23,23,23,
                           17,18,19,  24,25,26,
                           20,21,22,  27,28,29,
                           30,30,30,  37,37,37,
                           31,32,33,  38,39,40,
                           34,35,36,  41,42,43
                          ) ,
                        nrow = 10, ncol =  6, byrow = T ), 
          heights =  c(1,4,0.6,0.6,4,0.6,0.6,4,0.6,0.6) )
  
  
  for ( i_ev in 1:length(data_obs) ) {
    
    if ( i_ev %% 6 == 1  ) {
      par(mar = c(0,0,0,0))
      plot.new()
      legend("center", x.intersp = 0.2, text.width = 0.25, box.lwd = 0.5,
             legend = c( "Q observed - within predicted bounds", "Q observed - out of predicted bounds", 
                         "Q predicted - median", "90% bounds (Q predicted or performance metrics)" ),
             ncol = 2, 
             lty = c( 0, 0, 1, 0 ),
             col = c( "steelblue1",
                      "firebrick1",
                      "purple",
                      gray(0.6) ),
             lwd = c( NA, NA, 1, NA),
             pch = c( 1, 1, NA, 15 ), pt.cex = c(1, 1, 1.5, 2) )
    }
    
    
    ### plots hydrographs
    
    # selects data for the event
    data_obs_ev <- data_obs[[i_ev]]$Q_Data
    data_mod_ev <- data_mod[[i_ev]]
    
    Qobs <- data_obs_ev[ , 2 ]
    y    <- data_mod_ev[[ paste0("y.", L, ".quant") ]]
    y_B  <- data_mod_ev[[ paste0("yplusB.", L, ".quant") ]]
    Y    <- data_mod_ev[[ paste0("Y.", L, ".quant") ]]
    timestep <- sysanal.decode( colnames(y) )[,2]    
    timestep <- timestep - timestep[1] 
    
    # plots the hydrograph data
    par(mar = c(3.2, 3.2, 3, 1))
    plot( x = timestep, 
          y = Qobs, 
          type = "n", xaxt = "n", axes = F,
          ylab = NA, xlab = NA, 
          ylim = c( 0 , max( c( Y[(row.names(Y)=="0.95")] , Qobs ) ) ), )
    box (lwd = 0.5); axis(side = 1, lwd = 0.5); axis(side = 2, lwd = 0.5)
    
    
    title( paste0( names(data_obs)[i_ev], ",  Rmax10 = ",  
                   uni.data$RG.overview$meanRain_Rmax10[ uni.data$RG.overview$id %in% names(data_obs) ] [match(names(data_obs)[i_ev], as.character(names(data_obs)))], " mm/h") ) 
    mtext(side = 2, line = 2, "Discharge [l/s]", cex = 0.8)
    mtext(side = 1, line = 2, "Time [h]", cex = 0.8)
    
    
    polygon( c(timestep, rev( timestep)),
             c(Y[(row.names(Y)=="0.05")],
               rev(Y[(row.names(Y)=="0.95")])), #set the limits (1st and last quantiles)
             col=gray(0.6), 
             border=NA )
    
    # polygon( c(timestep,rev( timestep)),
    #          c(y_B[(row.names(y_B)=="0.05")],
    #            rev(y_B[(row.names(y_B)=="0.95")])), #set the limits (1st and last quantiles)
    #          col=gray(0.5), 
    #          border=NA )
    # 
    # polygon( c(timestep,rev( timestep)),
    #          c(y[(row.names(y)=="0.05")],
    #            rev(y[(row.names(y)=="0.95")])), #set the limits (1st and last quantiles)
    #          # col=gray(0.8),
    #          col = "magenta",
    #          border=NA )
    
    # median predictions
    lines( x = timestep, y = Y[(row.names(Y)=="0.5")], col = "purple")
    
    
    red_points <- red_points( obs = Qobs, 
                              up  = Y[(row.names(Y)=="0.95")], 
                              lo  = Y[(row.names(Y)=="0.05")] )
    blu_points <- Qobs
    blu_points[ !is.na(red_points) ] <- NA
    
    points(x = timestep, y = blu_points, col="steelblue1", pch = 1 )
    points(x = timestep, y = red_points, col="firebrick1",     pch = 1 )
    

    ### plots stats
    for ( i_stat in c("dV", "dQmax", "NNSE", "SCC") ) {
      
      if ( i_stat %in% c("dV", "dQmax") ) { ylim <- c(-1, 1) }
      if ( i_stat %in% c("NNSE") )        { ylim <- c(0.3,1) }
      if ( i_stat %in% c("SCC") )         { ylim <- c(0,  1) }
      
      par(mar = c(2, 1.8, 0.2, 0.2))
      vioplot::vioplot( stats$stats_it[[i_ev]][, i_stat] ,
                        xlab = NA , yaxt = "n", axes = F, ylim = ylim, border = NA,
                        plotCentre = "line", col = gray(0.6), range = 0, horizontal = T, lwd = 0.5 )
      
      # lightens the distribution extremes
      if ( i_stat %in% c("dV", "dQmax") ) {
        polygon( y = c(0.5, 0.5, 1.5, 1.5),
                 x = c( quantile(stats$stats_it[[i_ev]][, i_stat], c(0, 0.05)) ,
                        rev( quantile(stats$stats_it[[i_ev]][, i_stat], c(0, 0.05)) ) ),
                 col = rgb(1,1,1,0.6), border = rgb(1,1,1,0.6) )
        polygon( y = c(0.5, 0.5, 1.5, 1.5),
                 x = c( quantile(stats$stats_it[[i_ev]][, i_stat], c(0.95, 1)) ,
                        rev( quantile(stats$stats_it[[i_ev]][, i_stat], c(0.95, 1)) ) ),
                 col = rgb(1,1,1,0.6), border = rgb(1,1,1,0.6) )
      } else {
        polygon( y = c(0.5, 0.5, 1.5, 1.5),
                 x = c( quantile(stats$stats_it[[i_ev]][, i_stat], c(0, 0.1)) ,
                        rev( quantile(stats$stats_it[[i_ev]][, i_stat], c(0, 0.1)) ) ),
                 col = rgb(1,1,1,0.6), border = rgb(1,1,1,0.6) )
      }
      
      # lines for median predictions
      lines( y = c(0.8, 1.2) , 
             x = rep( stats$stats_qntl[[i_ev]][ "0.5" , i_stat], 2 ) ,
             col = "purple" )    
      
      # adds the axis and label
      axis(1, lwd = 0.5); abline( h = 0.46, lwd = 0.5, )
      if ( i_stat %in% c("dV", "dQmax") ) { abline( v = 0, lwd = 0.5, lty = "dashed") }
      mtext(side = 2, line = 0, text = i_stat, cex = 0.8, )
    }
    
    labels <- c( "reliab.", "NMISS" ); names(labels) <- c( "rlb.prcnt", "NMISS" )
    for ( i_stat in names(labels)  ) {
      
      if ( i_stat == "NMISS"   )   { ylim <- range(stats$stats_band[, i_stat]) }
      if ( i_stat == "rlb.prcnt" ) { # stats$stats_band[, i_stat] <- stats$stats_band[, i_stat] /100  # [%] --> [-]
                                     ylim <- range(stats$stats_band[, i_stat]) }
      
      boxplot( stats$stats_band[, i_stat], horizontal = T, range = 0, ylim = ylim,
               border = (gray(0.7)), axes = F, lwd = 0.7 )
      axis(1, lwd = 0.5); abline( h = 0.46 , lwd = 0.5, )
      mtext(side = 2, line = 0, text = labels[i_stat], cex = 0.8, )
      
      lines( y = c(0.65, 1.35) , 
             x = rep( stats$stats_band[i_ev, i_stat] , 2) ,
             col = "black" )
      points( x = stats$stats_band[i_ev, i_stat],
              y = 1, pch = 10,  col = "black", cex = 0.7  )
    }

  }
  dev.off()
  
  
  return(TRUE)   
}
  

# #---------------------------------------------------------------------
# 
# # plots statistics overview
# plot.Pre.res <- function( dataCa, dataPre, 
#                           bTr.Ca, bTr.Pre,
#                           pack.dir, statistics) {
# 
#   for (ii in 1 : length(to.plot.list[[1]]$statistics)) {
# 
# 
#     # boxplots C  ( all iterations )
#     pdf( paste(pack.dir, "/00_Pre_", names(to.plot.list[[1]]$statistics)[ii], "_stats_bxpltC.pdf", sep="") , height = 6,  width = 7)
#       nValues <- length(to.plot.list[[1]]$statistics[[ii]]$all.iterations[1,]) * length(to.plot.list[[1]]$statistics[[ii]]$all.iterations[,1])    
#       BoxPlot <- data.frame(matrix(NA, ncol=3, nrow= nValues * length(to.plot.list)))
#       colnames(BoxPlot) <- c("value", "metric", "datSrc")
#       
#       reliab <- c()
#       for (i in 1:length(to.plot.list)) {
#         for (j in 1:length(to.plot.list[[i]]$statistics[[ii]]$all.iterations[1,])) {
#           pos <- (i-1) * nValues + 
#             (j-1) * (length(to.plot.list[[i]]$statistics[[ii]]$all.iterations[,1])) + 
#             (1 : length(to.plot.list[[i]]$statistics[[ii]]$all.iterations[,1]))
#           BoxPlot$value [pos]  <- to.plot.list[[i]]$statistics[[ii]]$all.iterations[,j]
#           BoxPlot$metric[pos]  <- colnames(to.plot.list[[i]]$statistics[[ii]]$all.iterations)[j]
#           BoxPlot$datSrc[pos]  <- names(to.plot.list)[i]
#         }
#         reliab[i] <- to.plot.list[[i]]$statistics[[ii]]$ret.all["reliab"]
#         names(reliab)[i] <- names(to.plot.list)[[i]]
#       }
#       
#       par(tcl=-0.5, family="serif", omi=c(0,0,0,0), cex.lab=0.8, cex.axis=0.7, cex.main=0.8,
#           mai=c(0.3, 0.5, 0.2, 0.1), mgp=c(1.3, 0.6, 0) )
#       split.screen(c(2, 1))       # splits display into two screens
#       split.screen(c(1, 2), screen = 2) # splits the bottom half into 2
#       
#       # plot up 
#       screen(1) 
#         Plot1 <- c( which( BoxPlot$metric == paste(intToUtf8(0x03B4), "V(Y)", sep="") ),          # delta V
#                     which( BoxPlot$metric == paste(intToUtf8(0x03B4), "Vpeak(Y)", sep="") ),      # delta Vpeak
#                     which( BoxPlot$metric == "NS(Y)")                                             # NSE
#                    )     
#         boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n',
#                 col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7, 9, 10, 11),  ylab = "[-]",
#                 main = "C - all events and iterations",
#                 outline = F  # no outlyers!
#                 #ylim = c(-max(abs(BoxPlot$value[Plot1]), na.rm = T), max(abs(BoxPlot$value[Plot1]), na.rm = T)) 
#                 )
#         legend(x = "top", legend = names( to.plot.list)[order(names(to.plot.list))],
#                fill = c('white', 'gray80', 'gray30'), cex = 0.6, horiz = T )
#         axis(side = 1, at = c(2, 6, 10), line = -0.8, lwd = 0,
#              labels = c(paste(intToUtf8(0x03B4), "V", sep=""), paste(intToUtf8(0x03B4), "Vpeak", sep=""), "NSE") )
#         lines(x=c(0,12), y=c(0,0), lty=2 )
#       close.screen(1)      
#       
#       # plot down left
#       screen(3)
#         Plot1 <- c( which( BoxPlot$metric == "shift(Qmax(Y))" ) )     
#         boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1],  xaxt='n',
#                 col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3), ylab = "[h]", outline = F  # no outlyers!
#                 )
#         axis(side = 1, at = c(2), line = -0.8, lwd = 0,
#              labels = c( "shift(Qmax(Y))" ) )
#         lines(x=c(0,4), y=c(0,0), lty=2)
#       close.screen(3) 
#       
#       # plot down right
#       screen(4)
#         barCenters <- barplot(height = as.matrix(reliab[order(names(reliab))]), 
#                               ylim = c(0,1), beside = T, space=0.2,
#                               border = "black", axes = TRUE, ylab = "[-]",
#                               col = c('white', 'gray80', 'gray30'),
#                               main = "reliability")
#       close.screen(4) 
#       
#     close.screen(all = TRUE)
#     dev.off()
#   
#   }  
#   
#   return(TRUE)
# }  
# 
# 
# #---------------------------------------------------------------------
# 
# CaPre.plot.Pre <- function(plotdata1) {
#   timestepPre <- plotdata1$timestepPre;  
#   y_B.quant <- plotdata1$bct$bct.yplusB.L2.quant.Pre;  
#   y.quant   <- plotdata1$bct$bct.y.L2.quant.Pre; 
#   Y.quant   <- plotdata1$bct$bct.Y.L2.quant.Pre; Y.samp <- plotdata1$bct$bct.Y.L2.samp.Pre
#   bind.ret  <- plotdata1$bind.ret;  
#   VerInd    <- plotdata1$VerInd
#   transf <- plotdata1$transf
#   data.Pre <- plotdata1$data.Pre; out.data.Pre <- data.Pre[[2]]
#   
#   
#   par(mfrow=c(1,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0.1,0))
#   
#   # left plot
#   par(mai=c(0.4, 0.4, 0, 0))
#   plot(x = timestepPre, y = out.data.Pre[1:nrow(out.data.Pre),2], 
#        ylab = "Discharge [l/s]", xlab = "Timestep [h]",
#        ylim = c( 0 , max(Y.quant[(row.names(Y.quant)=="0.95")])*1.1 )
#       )
#   
#   polygon(c(timestepPre,rev(timestepPre)),
#           c(Y.quant[(row.names(Y.quant)=="0.05")],
#             rev(Y.quant[(row.names(Y.quant)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.1), border=NA)
#   
#   polygon(c(timestepPre,rev(timestepPre)),
#           c(y_B.quant[(row.names(y_B.quant)=="0.05")],
#             rev(y_B.quant[(row.names(y_B.quant)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.5), border=NA)
#   
#   polygon(c(timestepPre,rev(timestepPre)),
#           c(y.quant[(row.names(y.quant)=="0.05")],
#             rev(y.quant[(row.names(y.quant)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.8), border=NA)
#   
#   lines(x = timestepPre, y = apply(Y.samp, 2, mean), lty = "66")
#   
#   points(x = timestepPre, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep", col="blue")
#   
#   points(x = timestepPre, y = VerInd[2,], col="red")
#   
#   legend("topright", inset=0.00, legend=plotdata1$data.source,
#          pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=0.8)
#   
# #   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
# #                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
# #                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
# #                                    "; Pre event ", 
# # 
# #   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
# #                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
# #                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
# #          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
# #   
#   # right plot
#   par(mai=c(0.2, 0.0, 0.2, 0.0))
#   plot.new()
#   legend("top", inset=0.00, pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1,
#          legend = c(paste("event ", substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), sep="")
#          )
#   )
#   legend("top", inset=0.1, pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=0.8, ncol=2,
#                                         # "NS(E(Y))"  "delta E(V)"     "delta E(Vpeak)" "MIS"
#          legend = c( " ", names(plotdata1$bind.ret[c(2, 5, 8, 9, 12, 19)]),
#                      plotdata1$data.source, plotdata1$bind.ret[c(2, 5, 8, 9, 12, 19)]
#          )
#   )
#   
#   mtext("Timestep [h]", side=1, outer=T, at=0.25)
#   mtext("Discharge [l/s]", side=2, outer=T, at=0.5)
# }
# 
# #---------------------------------------------------------------------
# 
# CaPre.plot.2.Pre <- function(plotdata1, plotdata2) {
#   
#   (if (plotdata1$data.Pre[[3]] == plotdata2$data.Pre[[3]]) 
#   {
#     data.Pre <- plotdata1$data.Pre; out.data.Pre <- data.Pre[[2]]
#   }
#   else {stop("event mismatch")}
#   )
#   
#   timestepPre1 <- plotdata1$timestepPre-plotdata1$timestepPre[1]; 
#   y_B.quant1   <- plotdata1$bct$bct.yplusB.L2.quant.Pre;  
#   y.quant1     <- plotdata1$bct$bct.y.L2.quant.Pre; 
#   Y.quant1     <- plotdata1$bct$bct.Y.L2.quant.Pre; 
#   Y.samp1      <- plotdata1$bct$bct.Y.L2.samp.Pre
#   bind.ret1    <- plotdata1$bind.ret;  
#   VerInd1      <- plotdata1$VerInd
#   transf1      <- plotdata1$transf
#   
#   timestepPre2 <- plotdata2$timestepPre-plotdata2$timestepPre[1]; 
#   y_B.quant2   <- plotdata2$bct$bct.yplusB.L2.quant.Pre;  
#   y.quant2     <- plotdata2$bct$bct.y.L2.quant.Pre; 
#   Y.quant2     <- plotdata2$bct$bct.Y.L2.quant.Pre; 
#   Y.samp2      <- plotdata2$bct$bct.Y.L2.samp.Pre
#   bind.ret2    <- plotdata2$bind.ret;  
#   VerInd2      <- plotdata2$VerInd
#   transf2      <- plotdata2$transf
#   
#   par(mfrow=c(1,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0.2,0))
#   
#   # plot left
#   par(mai=c(0.4,0.4,0.02,0))
#   
#   plot(x = timestepPre1, y = out.data.Pre[1:nrow(out.data.Pre),2], bty="o",
#        ylab = " ", xlab = " ",
#        xlim = c(0, floor(max(timestepPre1)+0.7)),
#        ylim = c(0 , max(Y.quant1[(row.names(Y.quant1)=="0.95")])*1.1 )
#   )
#   
#   polygon(c(timestepPre1,rev(timestepPre1)),
#           c(Y.quant1[(row.names(Y.quant1)=="0.05")],
#             rev(Y.quant1[(row.names(Y.quant1)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.1), border=NA)
#   
#   polygon(c(timestepPre1,rev(timestepPre1)),
#           c(y_B.quant1[(row.names(y_B.quant1)=="0.05")],
#             rev(y_B.quant1[(row.names(y_B.quant1)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.5), border=NA)
#   
#   polygon(c(timestepPre1,rev(timestepPre1)),
#           c(y.quant1[(row.names(y.quant1)=="0.05")],
#             rev(y.quant1[(row.names(y.quant1)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.8), border=NA)
#   
#   lines(x = timestepPre1, y = apply(Y.samp1, 2, mean), lty = "66")
#   
#   points(x = timestepPre1, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep",  col="blue")
#   
#   points(x = timestepPre1, y = VerInd1[2,], col="red")
#  
#   legend("topright", inset=0.00, legend=plotdata3$data.source,       # legend, name of the rain data source 
#          pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1)
#   
# #   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
# #                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
# #                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
# #                                    "; Pre event ", 
# #                                    substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), ")", sep=""))
# #   
# #   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
# #                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
# #                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
# #          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
#   
#   # plot right
#   par(mai=c(0.4,0.2,0.02,0.2))
#   
#   plot(x = timestepPre2, y = out.data.Pre[1:nrow(out.data.Pre),2], bty="o", yaxt="n",
#        ylab = " ", xlab = " ",
#        xlim = c(0, floor(max(timestepPre2)+0.7)),
#        ylim = c( 0 , max(Y.quant2[(row.names(Y.quant2)=="0.95")])*1.1 )
#   )
#   
#   polygon(c(timestepPre2,rev(timestepPre2)),
#           c(Y.quant2[(row.names(Y.quant2)=="0.05")],
#             rev(Y.quant2[(row.names(Y.quant2)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.1), border=NA)
#   
#   polygon(c(timestepPre2,rev(timestepPre2)),
#           c(y_B.quant2[(row.names(y_B.quant2)=="0.05")],
#             rev(y_B.quant2[(row.names(y_B.quant2)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.5), border=NA)
#   
#   polygon(c(timestepPre2,rev(timestepPre2)),
#           c(y.quant2[(row.names(y.quant2)=="0.05")],
#             rev(y.quant2[(row.names(y.quant2)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.8), border=NA)
#   
#   lines(x = timestepPre2, y = apply(Y.samp2, 2, mean), lty = "66")
#   
#   points(x = timestepPre2, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep",  col="blue")
#   
#   points(x = timestepPre2, y = VerInd2[2,], col="red")
#   
#   legend("topright", inset=0.00, legend=plotdata3$data.source,
#          pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1)
#   
# #   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
# #                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
# #                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
# #                                    "; Pre event ", 
# #                                    substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), ")", sep=""))
# #   
# #   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
# #                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
# #                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
# #          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
#   
#   mtext("Timestep [h]", side=1, outer=T, at=0.5)
#   mtext("Discharge [l/s]", side=2, outer=T, at=0.5)
#   mtext(paste("event ", substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), sep=""),
#           side=3, outer=T, at=0.5)
#  
# }
# 
# #---------------------------------------------------------------------
# 
# CaPre.plot.3.Pre <- function(plotdata1, plotdata2, plotdata3) {
#   
#   (if ( (plotdata1$data.Pre[[3]] == plotdata2$data.Pre[[3]]) && (plotdata2$data.Pre[[3]] == plotdata3$data.Pre[[3]]) ) 
#   {
#     data.Pre <- plotdata1$data.Pre; out.data.Pre <- data.Pre[[2]]
#   }
#   else {stop("event mismatch")}
#   )
#   
#   timestepPre1 <- plotdata1$timestepPre-plotdata1$timestepPre[1]; 
#   y_B.quant1   <- plotdata1$bct$bct.yplusB.L2.quant.Pre;  
#   y.quant1     <- plotdata1$bct$bct.y.L2.quant.Pre; 
#   Y.quant1     <- plotdata1$bct$bct.Y.L2.quant.Pre; 
#   Y.samp1      <- plotdata1$bct$bct.Y.L2.samp.Pre
#   bind.ret1    <- plotdata1$bind.ret;  
#   VerInd1      <- plotdata1$VerInd
#   transf1      <- plotdata1$transf
#   
#   timestepPre2 <- plotdata2$timestepPre-plotdata2$timestepPre[1]; 
#   y_B.quant2   <- plotdata2$bct$bct.yplusB.L2.quant.Pre;  
#   y.quant2     <- plotdata2$bct$bct.y.L2.quant.Pre; 
#   Y.quant2     <- plotdata2$bct$bct.Y.L2.quant.Pre; 
#   Y.samp2      <- plotdata2$bct$bct.Y.L2.samp.Pre
#   bind.ret2    <- plotdata2$bind.ret;  
#   VerInd2      <- plotdata2$VerInd
#   transf2      <- plotdata2$transf
#   
#   timestepPre3 <- plotdata3$timestepPre-plotdata3$timestepPre[1]; 
#   y_B.quant3   <- plotdata3$bct$bct.yplusB.L2.quant.Pre;  
#   y.quant3     <- plotdata3$bct$bct.y.L2.quant.Pre; 
#   Y.quant3     <- plotdata3$bct$bct.Y.L2.quant.Pre; 
#   Y.samp3      <- plotdata3$bct$bct.Y.L2.samp.Pre
#   bind.ret3    <- plotdata3$bind.ret;  
#   VerInd3      <- plotdata3$VerInd
#   transf3      <- plotdata3$transf
#   
#   par(mfrow=c(2,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0.2,0))
#   
#   # plot up left 
#   par(mai=c(0.2, 0.4, 0.2, 0.1))
#   
#   plot(x = timestepPre1, y = out.data.Pre[1:nrow(out.data.Pre),2], bty="o", xaxt="n",
#        ylab = " ", xlab = " ",
#        xlim = c(0, floor(max(timestepPre1)+0.7)),
#        ylim = c(0 , max(c(Y.quant1[(row.names(Y.quant1)=="0.95")], Y.quant2[(row.names(Y.quant2)=="0.95")], 
#                           Y.quant3[(row.names(Y.quant3)=="0.95")], out.data.Pre[1:nrow(out.data.Pre),2]))
#                 *1.05 )
#   )
#   
#   polygon(c(timestepPre1,rev(timestepPre1)),
#           c(Y.quant1[(row.names(Y.quant1)=="0.05")],
#             rev(Y.quant1[(row.names(Y.quant1)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.1), border=NA)
#   
#   polygon(c(timestepPre1,rev(timestepPre1)),
#           c(y_B.quant1[(row.names(y_B.quant1)=="0.05")],
#             rev(y_B.quant1[(row.names(y_B.quant1)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.5), border=NA)
#   
#   polygon(c(timestepPre1,rev(timestepPre1)),
#           c(y.quant1[(row.names(y.quant1)=="0.05")],
#             rev(y.quant1[(row.names(y.quant1)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.8), border=NA)
#   
#   lines(x = timestepPre1, y = apply(Y.samp1, 2, mean), lty = "66")
#   
#   points(x = timestepPre1, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep",  col="blue")
#   
#   points(x = timestepPre1, y = VerInd1[2,], col="red")
#   
#   legend("topright", inset=0.00, legend =plotdata1$data.source,        # legend, name of the rain data source 
#          pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1)
#   
#   #   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
#   #                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
#   #                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
#   #                                    "; Pre event ", 
#   #                                    substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), ")", sep=""))
#   #   
#   #   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
#   #                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
#   #                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
#   #          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
#   
#   # plot up right
#   par(mai=c(0.2, 0.0, 0.2, 0.0))
#   plot.new()
#   legend("top", inset=0.00, pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1.2,
#          legend = c(paste("event ", substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), sep="")
#          )
#   )
#   legend("top", inset=0.15, pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=0.8, ncol=4,
#          legend = c( " ", names(plotdata1$bind.ret[c(8, 9, 13, 14, 18, 19, 20, 21, 22, 24, 25)]),
#                      plotdata1$data.source, plotdata1$bind.ret[c(8, 9, 13, 14, 18, 19, 20, 21, 22, 24, 25)],
#                      plotdata2$data.source, plotdata2$bind.ret[c(8, 9, 13, 14, 18, 19, 20, 21, 22, 24, 25)],
#                      plotdata3$data.source, plotdata3$bind.ret[c(8, 9, 13, 14, 18, 19, 20, 21, 22, 24, 25)]
#          )
#   )             
#   
#   # plot down left
#   par(mai=c(0.4, 0.4, 0, 0.1))
#   
#   plot(x = timestepPre2, y = out.data.Pre[1:nrow(out.data.Pre),2], bty="o",
#        ylab = " ", xlab = " ",
#        xlim = c(0, floor(max(timestepPre2)+0.7)),
#        ylim = c(0 , max(c(Y.quant1[(row.names(Y.quant1)=="0.95")], Y.quant2[(row.names(Y.quant2)=="0.95")], 
#                           Y.quant3[(row.names(Y.quant3)=="0.95")], out.data.Pre[1:nrow(out.data.Pre),2]))
#                 *1.05 )
#   )
#   
#   polygon(c(timestepPre2,rev(timestepPre2)),
#           c(Y.quant2[(row.names(Y.quant2)=="0.05")],
#             rev(Y.quant2[(row.names(Y.quant2)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.1), border=NA)
#   
#   polygon(c(timestepPre2,rev(timestepPre2)),
#           c(y_B.quant2[(row.names(y_B.quant2)=="0.05")],
#             rev(y_B.quant2[(row.names(y_B.quant2)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.5), border=NA)
#   
#   polygon(c(timestepPre2,rev(timestepPre2)),
#           c(y.quant2[(row.names(y.quant2)=="0.05")],
#             rev(y.quant2[(row.names(y.quant2)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.8), border=NA)
#   
#   lines(x = timestepPre2, y = apply(Y.samp2, 2, mean), lty = "66")
#   
#   points(x = timestepPre2, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep",  col="blue")
#   
#   points(x = timestepPre2, y = VerInd2[2,], col="red")
#   
#   legend("topright", inset=0.00, legend=plotdata2$data.source,
#          pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1)
#   
#   #   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
#   #                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
#   #                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
#   #                                    "; Pre event ", 
#   #                                    substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), ")", sep=""))
#   #   
#   #   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
#   #                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
#   #                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
#   #          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
#   
#   # plot down right
#   par(mai=c(0.4, 0.1, 0, 0.4))
#   
#   plot(x = timestepPre3, y = out.data.Pre[1:nrow(out.data.Pre),2], bty="o", yaxt="n",
#        ylab = " ", xlab = " ",
#        xlim = c(0, floor(max(timestepPre3)+0.7)),
#        ylim = c(0 , max(c(Y.quant1[(row.names(Y.quant1)=="0.95")], Y.quant2[(row.names(Y.quant2)=="0.95")], 
#                           Y.quant3[(row.names(Y.quant3)=="0.95")], out.data.Pre[1:nrow(out.data.Pre),2]))
#                 *1.05 )
#   )
#   
#   polygon(c(timestepPre3,rev(timestepPre3)),
#           c(Y.quant3[(row.names(Y.quant3)=="0.05")],
#             rev(Y.quant3[(row.names(Y.quant3)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.1), border=NA)
#   
#   polygon(c(timestepPre3,rev(timestepPre3)),
#           c(y_B.quant3[(row.names(y_B.quant3)=="0.05")],
#             rev(y_B.quant3[(row.names(y_B.quant3)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.5), border=NA)
#   
#   polygon(c(timestepPre3,rev(timestepPre3)),
#           c(y.quant3[(row.names(y.quant3)=="0.05")],
#             rev(y.quant3[(row.names(y.quant3)=="0.95")])), #set the limits (1st and last quantiles)
#           col=gray(0.8), border=NA)
#   
#   lines(x = timestepPre3, y = apply(Y.samp3, 2, mean), lty = "66")
#   
#   points(x = timestepPre3, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep",  col="blue")
#   
#   points(x = timestepPre3, y = VerInd3[2,], col="red")
#   
#   legend("topright", inset=0.00, legend=plotdata3$data.source,
#          pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1)
#   
#   #   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
#   #                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
#   #                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
#   #                                    "; Pre event ", 
#   #                                    substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), ")", sep=""))
#   #   
#   #   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
#   #                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
#   #                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
#   #          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
#   
#   mtext("Timestep [h]", side=1, outer=T, at=0.5)
#   mtext("Discharge [l/s]", side=2, outer=T, at=0.5)
#   
# }
# 
# #---------------------------------------------------------------------