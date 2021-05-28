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
 
  up = data_mod[(row.names(data_mod)=="0.95"),]
  lo = data_mod[(row.names(data_mod)=="0.05"),]
  ob = o
  Out <- matrix(c(lo,up,ob), ncol = length(ob), nrow = 3, byrow = T,
                    dimnames = list(c("lo","up","ob"),paste("Q",round(t.pre, digits = 5),sep="_")) )
  
  red_points <- red_points( obs = ob, up = up, lo = lo )

  
  if ( length(which(is.na(ob))) > 0 ) { 
    avail.data <- ob[ - which(is.na(ob)) ] 
  } else {
    avail.data <- ob 
  }
  
  reliab <-  ( 1 - length(which(!is.na(red_points))) / length(avail.data) ) * 100  # prediction reliability [%]
  
  ABW <- mean(up - lo)                 # Average Band Width
  relABW <- as.numeric(ABW) / sd(ob)   # Average Band Width relative to standard deviation of observations
  
  QuSc.pre <- apply(Out, 2, quscore)     # Interval Scores
  MIS <- mean( QuSc.pre , na.rm = TRUE)  # Mean of Interval Scores
  relMIS <- MIS / sd(ob)                 # Mean of Interval Scores relative to standard deviation of observations
  
  MIS_ABW <- MIS / ABW        
  
  ret <- c( ABW = ABW, relABW = relABW, rlb.prcnt = reliab, MIS = MIS, relMIS = relMIS, MIS_ABW = MIS_ABW )
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

# compute quantile score as defined by Gneiting 2007, use input from a single time step as vector, 
# use apply to evaluate time series ( res=apply(data,2,quscore,conf=0.05) )
# conf = confidence level, corresponding to (1-conf)% interval
quscore = function(x, conf=0.1) {
  low=x[1] #lower bound
  upp=x[2] #upper bound
  obs=x[3] #observations
  sh=upp-low #sharpness
  undersh=(low-obs)*(low>obs)
  oversh=(obs-upp)*(obs>upp)
  score=sh+2/conf*(undersh+oversh)
  return(score)
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


# 
# statist.CaPre.res <- function(data_mod, data_obs, skip) {   # skip - number of event ignored when calculating overall stats (ret.all)
#   
#   Y_quant <- names(data_mod[[1]])[ grepl("Y.L..quant", names(data_mod[[1]]) ) ]
#   
#   # Verification Indicies (ABW, reliab, MIS, red points)
#   VerInd.statistics <- c("ABW", "relABW", "reliab", "MIS", "relMIS", "n.timesteps", "n.red_points")
#   VerInd <- list();  VerInd.stat <-  data.frame(matrix(NA, ncol=length(VerInd.statistics), nrow=length(data_obs)))
#   names(VerInd.stat) <- VerInd.statistics
#   for (i in 1 : length(data_mod)) {
#     hlp <- CaPre.VerInd.Pre(data_obs = data_obs[[i]][[2]], data_mod = data_mod[[i]][[Y_quant]])
#     VerInd.stat[i,] <- as.numeric( c(hlp$ABW, hlp$relABW, hlp$reliab, hlp$MIS, hlp$relMIS, 
#                                    length(hlp$red_points), length(which(is.na(hlp$red_points) == TRUE)) )
#                                   )  
#     VerInd[[i]] <- data.frame(matrix(NA, ncol=length(hlp$QuSc), nrow=2)); names(VerInd[[i]]) <- data_obs[[i]][[1]]
#     VerInd[[i]][1,] <- hlp$QuSc; VerInd[[i]][2,] <- hlp$red_points
#   }
#   
#   
#   # NSE, Vtot and Vpeak
#   my.stats   <- c("id", "NSE(E(Y))", "NSE(Y_95)", "NSE(Y_05)", 
#                   
#                   paste(intToUtf8(0x03B4), "V(E(Y))", sep=""), paste(intToUtf8(0x03B4), "V(Y_95)", sep=""), 
#                                                                paste(intToUtf8(0x03B4), "V(Y_05)", sep=""),
#                   
#                   paste( "E(", intToUtf8(0x03B4), "V)", sep=""),  paste( "sd(", intToUtf8(0x03B4), "V)", sep=""),
#                   
#                   paste(intToUtf8(0x03B4), "Vpeak(E(Y))", sep=""), paste(intToUtf8(0x03B4), "Vpeak(Y_95)", sep=""), 
#                                                                    paste(intToUtf8(0x03B4), "Vpeak(Y_05)", sep=""),
#                   
#                   paste( "E(", intToUtf8(0x03B4), "Vpeak)", sep=""),  paste( "sd(", intToUtf8(0x03B4), "Vpeak)", sep=""),
#                   
#                   "shift(Qmax(E(Y)))", "shift(Qmax(Y_95))", "shift(Qmax(Y_05))",
#                   
#                   paste( "E(shift(Qmax))", sep=""),  paste( "sd(shift(Qmax))", sep=""),
#                   
#                   paste( "E(NSE)", sep=""),  paste( "sd(NSE)", sep="") )
#   
#   ret <- data.frame(matrix(NA, ncol=length(my.stats), nrow=length(data_obs)))
#   names(ret) <- my.stats
#   
#   for ( i in 1 : length(data_obs) ) {
#     id <- substr( data_obs[[i]][[3]], nchar(data_obs[[i]][[3]])-22, nchar(data_obs[[i]][[3]])-4 ) # event id (starting time)
#     
#     bct <- data_mod[[i]]
#     timestep  <- sysanal.decode( colnames( data_mod[[i]][[Y_quant]] ) )[,2]
#     obs <- data_obs[[i]][[2]][,2]
#     Vobs      <- Vtot (Qdata = obs, timestep = timestep)
#     Vpeak.obs <- Vpeak(Qdata = obs, timestep = timestep)
#     
#     # statistics for all predicted data
#     if (i==1) {
#       my.stats.event <- c("NS(Y)", "V(Y)", "Vpeak(Y)", "shift(Qmax(Y))", 
#                           paste(intToUtf8(0x03B4), "V(Y)", sep=""), paste(intToUtf8(0x03B4), "Vpeak(Y)", sep="") )
#       event.table.all <-  data.frame(matrix(NA, ncol=length(my.stats.event), nrow=length(bct[[Y_quant]][,1])*length(data_obs) ))
#       names(event.table.all) <- my.stats.event
#     }
#     event.table <-  data.frame(matrix(NA, ncol=length(my.stats.event), nrow=length(bct[[Y_quant]][,1])))
#     names(event.table) <- my.stats.event
#     
#     for (j in 1 : length(bct[[Y_quant]][,1]) ) {
#       Qmod  <- bct[[Y_quant]][j,]
#       
#       NS     <- enesko(mod = Qmod, obs = obs)
#       V      <- Vtot(Qdata = Qmod, timestep = timestep)
#       Vpeak  <- Vpeak(Qdata = Qmod, timestep = timestep)
#       shift  <- timestep[which(Qmod==max(Qmod))] - timestep[which(obs==max(obs))]; shift <- round( shift, 3)
#       dV     <- round( (V - Vobs) / Vobs, 3 )
#       dVpeak <- round( (Vpeak - Vpeak.obs) / Vpeak.obs, 3 )
# 
#       event.table[j,] <- c(NS, V, Vpeak, shift, dV, dVpeak)
#     }
#     
#     ignore <- F
#     if (length(skip) > 0 ) {          # checks whether to ignore the given event for the overall statistics  
#       for (j in 1 : length(skip) ) {
#         if (i == skip[j]) {
#           ignore <- T
#         }
#       }  
#     }
#     if (ignore == F) {
#       event.table.all[ ((i-1)*length(bct[[Y_quant]][,1]) + 1) : (i*length(bct[[Y_quant]][,1])), ] <- event.table
#     }
#     
#     
#     # statistics for E(Y), Y_05 and Y_95
#     Qmod.med  <- bct[[Y_quant]][(row.names(bct[[Y_quant]])=="0.5"),]
#     Qmod.95   <- bct[[Y_quant]][(row.names(bct[[Y_quant]])=="0.95"),]
#     Qmod.05   <- bct[[Y_quant]][(row.names(bct[[Y_quant]])=="0.05"),]
#     
#     # NS efficiency
#     NS.med <- enesko(mod = Qmod.med, obs = obs)
#     NS.95   <- enesko(mod = Qmod.95, obs = obs)
#     NS.05   <- enesko(mod = Qmod.05, obs = obs)
#     
#     # relative errors delta for total V
#     Vmod.med <- Vtot(Qdata = Qmod.med, timestep = timestep)
#     Vmod.95   <- Vtot(Qdata = Qmod.95  , timestep = timestep)
#     Vmod.05   <- Vtot(Qdata = Qmod.05  , timestep = timestep)
#     
#     deltaV.med <-  round( (Vmod.med - Vobs) / Vobs, 3 )
#     deltaV.95   <-  round( (Vmod.95   - Vobs) / Vobs, 3 )
#     deltaV.05   <-  round( (Vmod.05   - Vobs) / Vobs, 3 )
#     
#     # relative errors delta for peak V
#     Vpeak.mod.med <- Vpeak(Qdata = Qmod.med, timestep = timestep)
#     Vpeak.mod.95   <- Vpeak(Qdata = Qmod.95  , timestep = timestep)
#     Vpeak.mod.05   <- Vpeak(Qdata = Qmod.05  , timestep = timestep)
#     
#     deltaVpeak.med <-  round( (Vpeak.mod.med - Vpeak.obs) / Vpeak.obs, 3 )
#     deltaVpeak.95   <-  round( (Vpeak.mod.95   - Vpeak.obs) / Vpeak.obs, 3 )
#     deltaVpeak.05   <-  round( (Vpeak.mod.05   - Vpeak.obs) / Vpeak.obs, 3 )
#     
#     # Qmax time shifts
#     shift.med <- time.shift(series1 = Qmod.med, series2 = obs, timestep = timestep)
#     shift.95   <- time.shift(series1 = Qmod.95,   series2 = obs, timestep = timestep)
#     shift.05   <- time.shift(series1 = Qmod.05,   series2 = obs, timestep = timestep)
#     
#     # interval scores for total V and peak V
#     IS.V     <- quscore(x=c(Vmod.05, Vmod.95, Vobs), conf=0.1); IS.V <- round(IS.V, 0) # [l]
#     IS.Vpeak <- quscore(x=c(Vpeak.mod.05, Vpeak.mod.95, Vpeak.obs), conf=0.1); IS.Vpeak <- round(IS.Vpeak, 0) # [l]
#     
#     
#     ret[i,] <- c(id, NS.med, NS.95, NS.05, 
#                  deltaV.med, deltaV.95, deltaV.05, 
#                  round(mean(event.table[,5]), 3), round(sd(event.table[,5]), 3), # dV
#                  deltaVpeak.med, deltaVpeak.95, deltaVpeak.05,
#                  round(mean(event.table[,6]), 3), round(sd(event.table[,6]), 3), # dVpeak
#                  shift.med, shift.95, shift.05,
#                  round(mean(event.table[,4]), 3), round(sd(event.table[,4]), 3), # time shift
#                  round(mean(event.table[,1]), 3), round(sd(event.table[,1]), 3)  # NS
#     )
#   }
#   
#   bind.ret <- cbind(ret, VerInd.stat)
#   for (i in 2:length(bind.ret[1,])) {          # converts numeric values back to numeric
#     bind.ret[,i] <- as.numeric(bind.ret[,i])
#   }
#   
#   # E and sd for all events together; it is more informative to use the boxplot overview
#   ret.all <- c(round(mean(event.table.all[,5], na.rm=T), 3), round(sd(event.table.all[,5], na.rm=T), 3),  # dV
#                round(mean(event.table.all[,6], na.rm=T), 3), round(sd(event.table.all[,6], na.rm=T), 3),  # dVpeak
#                round(mean(event.table.all[,4], na.rm=T), 3), round(sd(event.table.all[,4], na.rm=T), 3),  # time shift
#                round(mean(event.table.all[,1], na.rm=T), 3), round(sd(event.table.all[,1], na.rm=T), 3), # NS
#                round( sum(as.numeric(bind.ret$n.red_points[!1:length(bind.ret$n.red_points) %in% skip])) / 
#                       sum(as.numeric(bind.ret$n.timesteps [!1:length(bind.ret$n.timesteps) %in% skip]))  , 3) # reliab
#               )
#   names(ret.all) <- c(my.stats[c(8, 9, 13, 14, 18, 19, 20, 21)], "reliab")
#   
#   # data for boxplots A and B
#   boxplot.data <- (bind.ret[ , -c(1, 27, 28)])
#   if (length(skip) > 0 ) {          # checks which events to ignore for the overall statistics  
#     boxplot.data <- boxplot.data[-skip,]  
#   }
#   
#   # data for boxplot C
#   all.iterations <- event.table.all[, -c(2,3)]
#   
#   return(list(bind.ret = bind.ret, VerInd=VerInd, ret.all=ret.all, boxplot.data=boxplot.data, all.iterations=all.iterations))
# }


#---------------------------------------------------------------------

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

# plots statistics overview
plot_stats <- function( stats, data_obs ) {
 
  for ( i_ev in 1:length(stats$stats_qntl) ) {
    
    for ( i_stat in colnames( stats$stats_qntl[[1]] ) ) {
      stats_ev <- stats$stats_qntl[[i_ev]] [ rownames(stats$stats_qntl[[i_ev]]) %in% as.character(seq(0.05, 0.95, 0.005)) , ]
      assign( paste0( i_stat, "_ev" ) , stats_ev[ , i_stat ] )
      
      if ( i_ev == 1 ) {
        assign( paste0( i_stat, "_all_ev" ) , eval( parse( text = paste0( i_stat, "_ev" ) ) ) ) 
      } else {
        assign( paste0( i_stat, "_all_ev" ) , c( eval( parse( text = paste0( i_stat, "_all_ev" ) ) ), 
                                                 eval( parse( text = paste0( i_stat, "_ev"     ) ) ) )   )
      }
      
    }
    
  }
  
  png( paste0(out_dir, "/stats_qntl.png") ,
       type="cairo", units = "in", width = 3, height = 5*4, res = 150 )
  
    par( mar = c(1, 4, 1, 1), mfrow = c(4, 1), cex = 1.5 )
    
    for ( i_stat in c("dV", "dQmax", "NNSE", "SCC") ) {
      boxplot( eval( parse( text = paste0( i_stat, "_all_ev" ) ) ) ,
               ylab = i_stat ,
               range = 0  )
    }
  dev.off()
  
  
  #-----------------------------------------------
  
  for ( i_ev in 1:length(stats$stats_it) ) {
    
    for ( i_stat in colnames( stats$stats_it[[1]] ) ) {
      stats_ev <- stats$stats_it[[i_ev]] 
      assign( paste0( i_stat, "_ev" ) , stats_ev[ , i_stat ] )
      
      if ( i_ev == 1 ) {
        assign( paste0( i_stat, "_all_ev" ) , eval( parse( text = paste0( i_stat, "_ev" ) ) ) ) 
      } else {
        assign( paste0( i_stat, "_all_ev" ) , c( eval( parse( text = paste0( i_stat, "_all_ev" ) ) ), 
                                                 eval( parse( text = paste0( i_stat, "_ev"     ) ) ) )   )
      }
      
    }
    
  }
  
  png( paste0(out_dir, "/stats_it.png") ,
       type="cairo", units = "in", width = 3, height = 5*4, res = 150 )
  
  par( mar = c(1, 4, 1, 1), mfrow = c(4, 1), cex = 1.5 )
  
  for ( i_stat in c("dV", "dQmax", "NNSE", "SCC") ) {
    boxplot( eval( parse( text = paste0( i_stat, "_all_ev" ) ) ) ,
             ylab = i_stat ,
             range = 0  )
  }
  dev.off()
  
  return(TRUE)   
}
  

#---------------------------------------------------------------------

# plots statistics overview
plot.Pre.res <- function( dataCa, dataPre, 
                          bTr.Ca, bTr.Pre,
                          pack.dir, statistics) {

  for (ii in 1 : length(to.plot.list[[1]]$statistics)) {


    # boxplots C  ( all iterations )
    pdf( paste(pack.dir, "/00_Pre_", names(to.plot.list[[1]]$statistics)[ii], "_stats_bxpltC.pdf", sep="") , height = 6,  width = 7)
      nValues <- length(to.plot.list[[1]]$statistics[[ii]]$all.iterations[1,]) * length(to.plot.list[[1]]$statistics[[ii]]$all.iterations[,1])    
      BoxPlot <- data.frame(matrix(NA, ncol=3, nrow= nValues * length(to.plot.list)))
      colnames(BoxPlot) <- c("value", "metric", "datSrc")
      
      reliab <- c()
      for (i in 1:length(to.plot.list)) {
        for (j in 1:length(to.plot.list[[i]]$statistics[[ii]]$all.iterations[1,])) {
          pos <- (i-1) * nValues + 
            (j-1) * (length(to.plot.list[[i]]$statistics[[ii]]$all.iterations[,1])) + 
            (1 : length(to.plot.list[[i]]$statistics[[ii]]$all.iterations[,1]))
          BoxPlot$value [pos]  <- to.plot.list[[i]]$statistics[[ii]]$all.iterations[,j]
          BoxPlot$metric[pos]  <- colnames(to.plot.list[[i]]$statistics[[ii]]$all.iterations)[j]
          BoxPlot$datSrc[pos]  <- names(to.plot.list)[i]
        }
        reliab[i] <- to.plot.list[[i]]$statistics[[ii]]$ret.all["reliab"]
        names(reliab)[i] <- names(to.plot.list)[[i]]
      }
      
      par(tcl=-0.5, family="serif", omi=c(0,0,0,0), cex.lab=0.8, cex.axis=0.7, cex.main=0.8,
          mai=c(0.3, 0.5, 0.2, 0.1), mgp=c(1.3, 0.6, 0) )
      split.screen(c(2, 1))       # splits display into two screens
      split.screen(c(1, 2), screen = 2) # splits the bottom half into 2
      
      # plot up 
      screen(1) 
        Plot1 <- c( which( BoxPlot$metric == paste(intToUtf8(0x03B4), "V(Y)", sep="") ),          # delta V
                    which( BoxPlot$metric == paste(intToUtf8(0x03B4), "Vpeak(Y)", sep="") ),      # delta Vpeak
                    which( BoxPlot$metric == "NS(Y)")                                             # NSE
                   )     
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n',
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7, 9, 10, 11),  ylab = "[-]",
                main = "C - all events and iterations",
                outline = F  # no outlyers!
                #ylim = c(-max(abs(BoxPlot$value[Plot1]), na.rm = T), max(abs(BoxPlot$value[Plot1]), na.rm = T)) 
                )
        legend(x = "top", legend = names( to.plot.list)[order(names(to.plot.list))],
               fill = c('white', 'gray80', 'gray30'), cex = 0.6, horiz = T )
        axis(side = 1, at = c(2, 6, 10), line = -0.8, lwd = 0,
             labels = c(paste(intToUtf8(0x03B4), "V", sep=""), paste(intToUtf8(0x03B4), "Vpeak", sep=""), "NSE") )
        lines(x=c(0,12), y=c(0,0), lty=2 )
      close.screen(1)      
      
      # plot down left
      screen(3)
        Plot1 <- c( which( BoxPlot$metric == "shift(Qmax(Y))" ) )     
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1],  xaxt='n',
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3), ylab = "[h]", outline = F  # no outlyers!
                )
        axis(side = 1, at = c(2), line = -0.8, lwd = 0,
             labels = c( "shift(Qmax(Y))" ) )
        lines(x=c(0,4), y=c(0,0), lty=2)
      close.screen(3) 
      
      # plot down right
      screen(4)
        barCenters <- barplot(height = as.matrix(reliab[order(names(reliab))]), 
                              ylim = c(0,1), beside = T, space=0.2,
                              border = "black", axes = TRUE, ylab = "[-]",
                              col = c('white', 'gray80', 'gray30'),
                              main = "reliability")
      close.screen(4) 
      
    close.screen(all = TRUE)
    dev.off()
  
  }  
  
  return(TRUE)
}  


#---------------------------------------------------------------------

CaPre.plot.Pre <- function(plotdata1) {
  timestepPre <- plotdata1$timestepPre;  
  y_B.quant <- plotdata1$bct$bct.yplusB.L2.quant.Pre;  
  y.quant   <- plotdata1$bct$bct.y.L2.quant.Pre; 
  Y.quant   <- plotdata1$bct$bct.Y.L2.quant.Pre; Y.samp <- plotdata1$bct$bct.Y.L2.samp.Pre
  bind.ret  <- plotdata1$bind.ret;  
  VerInd    <- plotdata1$VerInd
  transf <- plotdata1$transf
  data.Pre <- plotdata1$data.Pre; out.data.Pre <- data.Pre[[2]]
  
  
  par(mfrow=c(1,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0.1,0))
  
  # left plot
  par(mai=c(0.4, 0.4, 0, 0))
  plot(x = timestepPre, y = out.data.Pre[1:nrow(out.data.Pre),2], 
       ylab = "Discharge [l/s]", xlab = "Timestep [h]",
       ylim = c( 0 , max(Y.quant[(row.names(Y.quant)=="0.95")])*1.1 )
      )
  
  polygon(c(timestepPre,rev(timestepPre)),
          c(Y.quant[(row.names(Y.quant)=="0.05")],
            rev(Y.quant[(row.names(Y.quant)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.1), border=NA)
  
  polygon(c(timestepPre,rev(timestepPre)),
          c(y_B.quant[(row.names(y_B.quant)=="0.05")],
            rev(y_B.quant[(row.names(y_B.quant)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.5), border=NA)
  
  polygon(c(timestepPre,rev(timestepPre)),
          c(y.quant[(row.names(y.quant)=="0.05")],
            rev(y.quant[(row.names(y.quant)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.8), border=NA)
  
  lines(x = timestepPre, y = apply(Y.samp, 2, mean), lty = "66")
  
  points(x = timestepPre, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep", col="blue")
  
  points(x = timestepPre, y = VerInd[2,], col="red")
  
  legend("topright", inset=0.00, legend=plotdata1$data.source,
         pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=0.8)
  
#   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
#                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
#                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
#                                    "; Pre event ", 
# 
#   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
#                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
#                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
#          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
#   
  # right plot
  par(mai=c(0.2, 0.0, 0.2, 0.0))
  plot.new()
  legend("top", inset=0.00, pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1,
         legend = c(paste("event ", substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), sep="")
         )
  )
  legend("top", inset=0.1, pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=0.8, ncol=2,
                                        # "NS(E(Y))"  "delta E(V)"     "delta E(Vpeak)" "MIS"
         legend = c( " ", names(plotdata1$bind.ret[c(2, 5, 8, 9, 12, 19)]),
                     plotdata1$data.source, plotdata1$bind.ret[c(2, 5, 8, 9, 12, 19)]
         )
  )
  
  mtext("Timestep [h]", side=1, outer=T, at=0.25)
  mtext("Discharge [l/s]", side=2, outer=T, at=0.5)
}

#---------------------------------------------------------------------

CaPre.plot.2.Pre <- function(plotdata1, plotdata2) {
  
  (if (plotdata1$data.Pre[[3]] == plotdata2$data.Pre[[3]]) 
  {
    data.Pre <- plotdata1$data.Pre; out.data.Pre <- data.Pre[[2]]
  }
  else {stop("event mismatch")}
  )
  
  timestepPre1 <- plotdata1$timestepPre-plotdata1$timestepPre[1]; 
  y_B.quant1   <- plotdata1$bct$bct.yplusB.L2.quant.Pre;  
  y.quant1     <- plotdata1$bct$bct.y.L2.quant.Pre; 
  Y.quant1     <- plotdata1$bct$bct.Y.L2.quant.Pre; 
  Y.samp1      <- plotdata1$bct$bct.Y.L2.samp.Pre
  bind.ret1    <- plotdata1$bind.ret;  
  VerInd1      <- plotdata1$VerInd
  transf1      <- plotdata1$transf
  
  timestepPre2 <- plotdata2$timestepPre-plotdata2$timestepPre[1]; 
  y_B.quant2   <- plotdata2$bct$bct.yplusB.L2.quant.Pre;  
  y.quant2     <- plotdata2$bct$bct.y.L2.quant.Pre; 
  Y.quant2     <- plotdata2$bct$bct.Y.L2.quant.Pre; 
  Y.samp2      <- plotdata2$bct$bct.Y.L2.samp.Pre
  bind.ret2    <- plotdata2$bind.ret;  
  VerInd2      <- plotdata2$VerInd
  transf2      <- plotdata2$transf
  
  par(mfrow=c(1,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0.2,0))
  
  # plot left
  par(mai=c(0.4,0.4,0.02,0))
  
  plot(x = timestepPre1, y = out.data.Pre[1:nrow(out.data.Pre),2], bty="o",
       ylab = " ", xlab = " ",
       xlim = c(0, floor(max(timestepPre1)+0.7)),
       ylim = c(0 , max(Y.quant1[(row.names(Y.quant1)=="0.95")])*1.1 )
  )
  
  polygon(c(timestepPre1,rev(timestepPre1)),
          c(Y.quant1[(row.names(Y.quant1)=="0.05")],
            rev(Y.quant1[(row.names(Y.quant1)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.1), border=NA)
  
  polygon(c(timestepPre1,rev(timestepPre1)),
          c(y_B.quant1[(row.names(y_B.quant1)=="0.05")],
            rev(y_B.quant1[(row.names(y_B.quant1)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.5), border=NA)
  
  polygon(c(timestepPre1,rev(timestepPre1)),
          c(y.quant1[(row.names(y.quant1)=="0.05")],
            rev(y.quant1[(row.names(y.quant1)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.8), border=NA)
  
  lines(x = timestepPre1, y = apply(Y.samp1, 2, mean), lty = "66")
  
  points(x = timestepPre1, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep",  col="blue")
  
  points(x = timestepPre1, y = VerInd1[2,], col="red")
 
  legend("topright", inset=0.00, legend=plotdata3$data.source,       # legend, name of the rain data source 
         pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1)
  
#   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
#                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
#                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
#                                    "; Pre event ", 
#                                    substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), ")", sep=""))
#   
#   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
#                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
#                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
#          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
  
  # plot right
  par(mai=c(0.4,0.2,0.02,0.2))
  
  plot(x = timestepPre2, y = out.data.Pre[1:nrow(out.data.Pre),2], bty="o", yaxt="n",
       ylab = " ", xlab = " ",
       xlim = c(0, floor(max(timestepPre2)+0.7)),
       ylim = c( 0 , max(Y.quant2[(row.names(Y.quant2)=="0.95")])*1.1 )
  )
  
  polygon(c(timestepPre2,rev(timestepPre2)),
          c(Y.quant2[(row.names(Y.quant2)=="0.05")],
            rev(Y.quant2[(row.names(Y.quant2)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.1), border=NA)
  
  polygon(c(timestepPre2,rev(timestepPre2)),
          c(y_B.quant2[(row.names(y_B.quant2)=="0.05")],
            rev(y_B.quant2[(row.names(y_B.quant2)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.5), border=NA)
  
  polygon(c(timestepPre2,rev(timestepPre2)),
          c(y.quant2[(row.names(y.quant2)=="0.05")],
            rev(y.quant2[(row.names(y.quant2)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.8), border=NA)
  
  lines(x = timestepPre2, y = apply(Y.samp2, 2, mean), lty = "66")
  
  points(x = timestepPre2, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep",  col="blue")
  
  points(x = timestepPre2, y = VerInd2[2,], col="red")
  
  legend("topright", inset=0.00, legend=plotdata3$data.source,
         pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1)
  
#   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
#                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
#                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
#                                    "; Pre event ", 
#                                    substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), ")", sep=""))
#   
#   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
#                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
#                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
#          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
  
  mtext("Timestep [h]", side=1, outer=T, at=0.5)
  mtext("Discharge [l/s]", side=2, outer=T, at=0.5)
  mtext(paste("event ", substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), sep=""),
          side=3, outer=T, at=0.5)
 
}

#---------------------------------------------------------------------

CaPre.plot.3.Pre <- function(plotdata1, plotdata2, plotdata3) {
  
  (if ( (plotdata1$data.Pre[[3]] == plotdata2$data.Pre[[3]]) && (plotdata2$data.Pre[[3]] == plotdata3$data.Pre[[3]]) ) 
  {
    data.Pre <- plotdata1$data.Pre; out.data.Pre <- data.Pre[[2]]
  }
  else {stop("event mismatch")}
  )
  
  timestepPre1 <- plotdata1$timestepPre-plotdata1$timestepPre[1]; 
  y_B.quant1   <- plotdata1$bct$bct.yplusB.L2.quant.Pre;  
  y.quant1     <- plotdata1$bct$bct.y.L2.quant.Pre; 
  Y.quant1     <- plotdata1$bct$bct.Y.L2.quant.Pre; 
  Y.samp1      <- plotdata1$bct$bct.Y.L2.samp.Pre
  bind.ret1    <- plotdata1$bind.ret;  
  VerInd1      <- plotdata1$VerInd
  transf1      <- plotdata1$transf
  
  timestepPre2 <- plotdata2$timestepPre-plotdata2$timestepPre[1]; 
  y_B.quant2   <- plotdata2$bct$bct.yplusB.L2.quant.Pre;  
  y.quant2     <- plotdata2$bct$bct.y.L2.quant.Pre; 
  Y.quant2     <- plotdata2$bct$bct.Y.L2.quant.Pre; 
  Y.samp2      <- plotdata2$bct$bct.Y.L2.samp.Pre
  bind.ret2    <- plotdata2$bind.ret;  
  VerInd2      <- plotdata2$VerInd
  transf2      <- plotdata2$transf
  
  timestepPre3 <- plotdata3$timestepPre-plotdata3$timestepPre[1]; 
  y_B.quant3   <- plotdata3$bct$bct.yplusB.L2.quant.Pre;  
  y.quant3     <- plotdata3$bct$bct.y.L2.quant.Pre; 
  Y.quant3     <- plotdata3$bct$bct.Y.L2.quant.Pre; 
  Y.samp3      <- plotdata3$bct$bct.Y.L2.samp.Pre
  bind.ret3    <- plotdata3$bind.ret;  
  VerInd3      <- plotdata3$VerInd
  transf3      <- plotdata3$transf
  
  par(mfrow=c(2,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0.2,0))
  
  # plot up left 
  par(mai=c(0.2, 0.4, 0.2, 0.1))
  
  plot(x = timestepPre1, y = out.data.Pre[1:nrow(out.data.Pre),2], bty="o", xaxt="n",
       ylab = " ", xlab = " ",
       xlim = c(0, floor(max(timestepPre1)+0.7)),
       ylim = c(0 , max(c(Y.quant1[(row.names(Y.quant1)=="0.95")], Y.quant2[(row.names(Y.quant2)=="0.95")], 
                          Y.quant3[(row.names(Y.quant3)=="0.95")], out.data.Pre[1:nrow(out.data.Pre),2]))
                *1.05 )
  )
  
  polygon(c(timestepPre1,rev(timestepPre1)),
          c(Y.quant1[(row.names(Y.quant1)=="0.05")],
            rev(Y.quant1[(row.names(Y.quant1)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.1), border=NA)
  
  polygon(c(timestepPre1,rev(timestepPre1)),
          c(y_B.quant1[(row.names(y_B.quant1)=="0.05")],
            rev(y_B.quant1[(row.names(y_B.quant1)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.5), border=NA)
  
  polygon(c(timestepPre1,rev(timestepPre1)),
          c(y.quant1[(row.names(y.quant1)=="0.05")],
            rev(y.quant1[(row.names(y.quant1)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.8), border=NA)
  
  lines(x = timestepPre1, y = apply(Y.samp1, 2, mean), lty = "66")
  
  points(x = timestepPre1, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep",  col="blue")
  
  points(x = timestepPre1, y = VerInd1[2,], col="red")
  
  legend("topright", inset=0.00, legend =plotdata1$data.source,        # legend, name of the rain data source 
         pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1)
  
  #   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
  #                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
  #                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
  #                                    "; Pre event ", 
  #                                    substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), ")", sep=""))
  #   
  #   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
  #                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
  #                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
  #          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
  
  # plot up right
  par(mai=c(0.2, 0.0, 0.2, 0.0))
  plot.new()
  legend("top", inset=0.00, pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1.2,
         legend = c(paste("event ", substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), sep="")
         )
  )
  legend("top", inset=0.15, pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=0.8, ncol=4,
         legend = c( " ", names(plotdata1$bind.ret[c(8, 9, 13, 14, 18, 19, 20, 21, 22, 24, 25)]),
                     plotdata1$data.source, plotdata1$bind.ret[c(8, 9, 13, 14, 18, 19, 20, 21, 22, 24, 25)],
                     plotdata2$data.source, plotdata2$bind.ret[c(8, 9, 13, 14, 18, 19, 20, 21, 22, 24, 25)],
                     plotdata3$data.source, plotdata3$bind.ret[c(8, 9, 13, 14, 18, 19, 20, 21, 22, 24, 25)]
         )
  )             
  
  # plot down left
  par(mai=c(0.4, 0.4, 0, 0.1))
  
  plot(x = timestepPre2, y = out.data.Pre[1:nrow(out.data.Pre),2], bty="o",
       ylab = " ", xlab = " ",
       xlim = c(0, floor(max(timestepPre2)+0.7)),
       ylim = c(0 , max(c(Y.quant1[(row.names(Y.quant1)=="0.95")], Y.quant2[(row.names(Y.quant2)=="0.95")], 
                          Y.quant3[(row.names(Y.quant3)=="0.95")], out.data.Pre[1:nrow(out.data.Pre),2]))
                *1.05 )
  )
  
  polygon(c(timestepPre2,rev(timestepPre2)),
          c(Y.quant2[(row.names(Y.quant2)=="0.05")],
            rev(Y.quant2[(row.names(Y.quant2)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.1), border=NA)
  
  polygon(c(timestepPre2,rev(timestepPre2)),
          c(y_B.quant2[(row.names(y_B.quant2)=="0.05")],
            rev(y_B.quant2[(row.names(y_B.quant2)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.5), border=NA)
  
  polygon(c(timestepPre2,rev(timestepPre2)),
          c(y.quant2[(row.names(y.quant2)=="0.05")],
            rev(y.quant2[(row.names(y.quant2)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.8), border=NA)
  
  lines(x = timestepPre2, y = apply(Y.samp2, 2, mean), lty = "66")
  
  points(x = timestepPre2, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep",  col="blue")
  
  points(x = timestepPre2, y = VerInd2[2,], col="red")
  
  legend("topright", inset=0.00, legend=plotdata2$data.source,
         pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1)
  
  #   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
  #                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
  #                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
  #                                    "; Pre event ", 
  #                                    substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), ")", sep=""))
  #   
  #   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
  #                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
  #                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
  #          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
  
  # plot down right
  par(mai=c(0.4, 0.1, 0, 0.4))
  
  plot(x = timestepPre3, y = out.data.Pre[1:nrow(out.data.Pre),2], bty="o", yaxt="n",
       ylab = " ", xlab = " ",
       xlim = c(0, floor(max(timestepPre3)+0.7)),
       ylim = c(0 , max(c(Y.quant1[(row.names(Y.quant1)=="0.95")], Y.quant2[(row.names(Y.quant2)=="0.95")], 
                          Y.quant3[(row.names(Y.quant3)=="0.95")], out.data.Pre[1:nrow(out.data.Pre),2]))
                *1.05 )
  )
  
  polygon(c(timestepPre3,rev(timestepPre3)),
          c(Y.quant3[(row.names(Y.quant3)=="0.05")],
            rev(Y.quant3[(row.names(Y.quant3)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.1), border=NA)
  
  polygon(c(timestepPre3,rev(timestepPre3)),
          c(y_B.quant3[(row.names(y_B.quant3)=="0.05")],
            rev(y_B.quant3[(row.names(y_B.quant3)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.5), border=NA)
  
  polygon(c(timestepPre3,rev(timestepPre3)),
          c(y.quant3[(row.names(y.quant3)=="0.05")],
            rev(y.quant3[(row.names(y.quant3)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.8), border=NA)
  
  lines(x = timestepPre3, y = apply(Y.samp3, 2, mean), lty = "66")
  
  points(x = timestepPre3, y = out.data.Pre[1:nrow(out.data.Pre),2], ylab = "Discharge", xlab = "Timestep",  col="blue")
  
  points(x = timestepPre3, y = VerInd3[2,], col="red")
  
  legend("topright", inset=0.00, legend=plotdata3$data.source,
         pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1)
  
  #   legend("topleft", legend = paste("90% quantiles for y+B+E (", 
  #                                    names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
  #                                    names(transf$par.tr[2]),"=", transf$par.tr[2], 
  #                                    "; Pre event ", 
  #                                    substr(data.Pre[[3]], nchar(data.Pre[[3]])-22, nchar(data.Pre[[3]])-4), ")", sep=""))
  #   
  #   legend("topright", inset=.05, legend =c(paste("reliab = ", bind.ret$reliab, "%", sep=""),
  #                                           paste("ABW =", bind.ret$ABW), paste("ABW/sd =", bind.ret$relABW), 
  #                                           paste("MIS =", bind.ret$MIS), paste("MIS/sd =", bind.ret$relMIS)),   
  #          pch=c(NA, NA, NA),col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n",cex=1.2)
  
  mtext("Timestep [h]", side=1, outer=T, at=0.5)
  mtext("Discharge [l/s]", side=2, outer=T, at=0.5)
  
}

#---------------------------------------------------------------------