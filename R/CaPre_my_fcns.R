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

CaPre.VerInd.Pre <- function(out.data.Pre, bct.Y.L2.quant.Pre) {
  
  #o=out.data.Pre[2:nrow(out.data.Pre),2] # WHYYYYY to leave out the first one???
  o=out.data.Pre[1:nrow(out.data.Pre),2]
  t.pre=seq(1,length(o),by = 1)
 
  up.L2 = bct.Y.L2.quant.Pre[(row.names(bct.Y.L2.quant.Pre)=="0.95"),]
  lo.L2 = bct.Y.L2.quant.Pre[(row.names(bct.Y.L2.quant.Pre)=="0.05"),]
  ob.L2 = o
  Out.L2 <- matrix(c(lo.L2,up.L2,ob.L2), ncol = length(ob.L2), nrow = 3, byrow = T,
                    dimnames = list(c("lo.L2","up.L2","ob.L2"),paste("Q",round(t.pre, digits = 5),sep="_")) )
  
  
  avail.data <- which(!is.na(ob.L2))
  red.points <- rep(NA, length(ob.L2)) # data outside of the predicted band
  for (i in avail.data) {
      if ( (up.L2[i] < ob.L2[i])
           ||
           (lo.L2[i] > ob.L2[i])
          ) {
        red.points[i] <- ob.L2[i]  
      }
  }
  
  if ( length(which(is.na(ob.L2))) > 0 ) { avail.data <- ob.L2[ - which(is.na(ob.L2))] }
  reliab <-  format(( 1 - length(which(!is.na(red.points))) / length(avail.data) ) * 100, digits=3) # prediction reliability [%]
  
  ABW <- format(mean(up.L2 - lo.L2), digits=3) # Average Band Width
  relABW <- format(as.numeric(ABW) / sd(ob.L2), digits=3) # Average Band Width relative to standard deviation of observations
  
  QuSc.pre <- apply(Out.L2, 2, quscore)
  MIS.L2 <- format((mean( QuSc.pre , na.rm = TRUE)), digits=3) # Mean of Interval Scores
  relMIS.L2 <- format(as.numeric(MIS.L2) / sd(ob.L2), digits=3) # Mean of Interval Scores relative to standard deviation of observations
  
  ret <- list(ABW=ABW, relABW=relABW, reliab=reliab, MIS=MIS.L2, relMIS=relMIS.L2, QuSc=QuSc.pre, red.points=red.points)
  return(ret)
}

#---------------------------------------------------------------------


# calculates statistics for prediction events
statist.CaPre.res <- function(Pre.res, prodata, skip) {   # skip - number of event ignored when calculating overall stats (ret.all)
  
  dataPre <- prodata$Pre
  
  # Verification Indicies (ABW, reliab, MIS, red points)
  VerInd.statistics <- c("ABW", "relABW", "reliab", "MIS", "relMIS", "n.timesteps", "n.red.points")
  VerInd <- list();  VerInd.stat <-  data.frame(matrix(NA, ncol=length(VerInd.statistics), nrow=length(dataPre)))
  names(VerInd.stat) <- VerInd.statistics
  for (i in 1 : length(Pre.res$bTr.Pre)) {
    hlp <- CaPre.VerInd.Pre(out.data.Pre = dataPre[[i]][[2]], bct.Y.L2.quant.Pre = Pre.res$bTr.Pre[[i]]$bct.Y.L2.quant.Pre)
    VerInd.stat[i,] <- as.numeric( c(hlp$ABW, hlp$relABW, hlp$reliab, hlp$MIS, hlp$relMIS, 
                                   length(hlp$red.points), length(which(is.na(hlp$red.points) == TRUE)) )
                                  )  
    VerInd[[i]] <- data.frame(matrix(NA, ncol=length(hlp$QuSc), nrow=2)); names(VerInd[[i]]) <- dataPre[[i]][[1]]
    VerInd[[i]][1,] <- hlp$QuSc; VerInd[[i]][2,] <- hlp$red.points
  }
  
  
  # NSE, Vtot and Vpeak
  my.stats   <- c("id", "NSE(E(Y))", "NSE(Y_95)", "NSE(Y_05)", 
                  
                  paste(intToUtf8(0x03B4), "V(E(Y))", sep=""), paste(intToUtf8(0x03B4), "V(Y_95)", sep=""), 
                                                               paste(intToUtf8(0x03B4), "V(Y_05)", sep=""),
                  
                  paste( "E(", intToUtf8(0x03B4), "V)", sep=""),  paste( "sd(", intToUtf8(0x03B4), "V)", sep=""),
                  
                  paste(intToUtf8(0x03B4), "Vpeak(E(Y))", sep=""), paste(intToUtf8(0x03B4), "Vpeak(Y_95)", sep=""), 
                                                                   paste(intToUtf8(0x03B4), "Vpeak(Y_05)", sep=""),
                  
                  paste( "E(", intToUtf8(0x03B4), "Vpeak)", sep=""),  paste( "sd(", intToUtf8(0x03B4), "Vpeak)", sep=""),
                  
                  "shift(Qmax(E(Y)))", "shift(Qmax(Y_95))", "shift(Qmax(Y_05))",
                  
                  paste( "E(shift(Qmax))", sep=""),  paste( "sd(shift(Qmax))", sep=""),
                  
                  paste( "E(NSE)", sep=""),  paste( "sd(NSE)", sep="") )
  
  ret <- data.frame(matrix(NA, ncol=length(my.stats), nrow=length(dataPre)))
  names(ret) <- my.stats
  
  for (i in 1 : length(dataPre)) {
    id <- substr( dataPre[[i]][[3]], nchar(dataPre[[i]][[3]])-22, nchar(dataPre[[i]][[3]])-4 ) # event id (starting time)
    
    bct <- Pre.res$bTr.Pre[[i]]
    timestep <- Pre.res$bTr.Pre[[i]]$timestepPre
    obs <- dataPre[[i]][[2]][,2]
    Vobs      <- Vtot (Qdata = obs, timestep = timestep)
    Vpeak.obs <- Vpeak(Qdata = obs, timestep = timestep)
    
    # statistics for all predicted data
    if (i==1) {
      my.stats.event <- c("NS(Y)", "V(Y)", "Vpeak(Y)", "shift(Qmax(Y))", 
                          paste(intToUtf8(0x03B4), "V(Y)", sep=""), paste(intToUtf8(0x03B4), "Vpeak(Y)", sep="") )
      event.table.all <-  data.frame(matrix(NA, ncol=length(my.stats.event), nrow=length(bct$bct.Y.L2.samp.Pre[,1])*length(dataPre) ))
      names(event.table.all) <- my.stats.event
    }
    event.table <-  data.frame(matrix(NA, ncol=length(my.stats.event), nrow=length(bct$bct.Y.L2.samp.Pre[,1])))
    names(event.table) <- my.stats.event
    
    for (j in 1 : length(bct$bct.Y.L2.samp.Pre[,1]) ) {
      Qmod  <- bct$bct.Y.L2.samp.Pre[j,]
      
      NS     <- enesko(mod = Qmod, obs = obs)
      V      <- Vtot(Qdata = Qmod, timestep = timestep)
      Vpeak  <- Vpeak(Qdata = Qmod, timestep = timestep)
      shift  <- timestep[which(Qmod==max(Qmod))] - timestep[which(obs==max(obs))]; shift <- round( shift, 3)
      dV     <- round( (V - Vobs) / Vobs, 3 )
      dVpeak <- round( (Vpeak - Vpeak.obs) / Vpeak.obs, 3 )

      event.table[j,] <- c(NS, V, Vpeak, shift, dV, dVpeak)
    }
    
    ignore <- F
    if (length(skip) > 0 ) {          # checks whether to ignore the given event for the overall statistics  
      for (j in 1 : length(skip) ) {
        if (i == skip[j]) {
          ignore <- T
        }
      }  
    }
    if (ignore == F) {
      event.table.all[((i-1)*length(bct$bct.Y.L2.samp.Pre[,1]) + 1) : (i*length(bct$bct.Y.L2.samp.Pre[,1])),] <- event.table
    }
    
    
    # statistics for E(Y), Y_05 and Y_95
    Qmod.mean <- apply(bct$bct.Y.L2.samp.Pre, 2, mean)
    Qmod.95   <- bct$bct.Y.L2.quant.Pre[(row.names(bct$bct.Y.L2.quant.Pre)=="0.95"),]
    Qmod.05   <- bct$bct.Y.L2.quant.Pre[(row.names(bct$bct.Y.L2.quant.Pre)=="0.05"),]
    
      # NS efficiency
      NS.mean <- enesko(mod = Qmod.mean, obs = obs)
      NS.95   <- enesko(mod = Qmod.95, obs = obs)
      NS.05   <- enesko(mod = Qmod.05, obs = obs)
    
      # relative errors delta for total V
      Vmod.mean <- Vtot(Qdata = Qmod.mean, timestep = timestep)
      Vmod.95   <- Vtot(Qdata = Qmod.95  , timestep = timestep)
      Vmod.05   <- Vtot(Qdata = Qmod.05  , timestep = timestep)
      
      deltaV.mean <-  round( (Vmod.mean - Vobs) / Vobs, 3 )
      deltaV.95   <-  round( (Vmod.95   - Vobs) / Vobs, 3 )
      deltaV.05   <-  round( (Vmod.05   - Vobs) / Vobs, 3 )
      
      # relative errors delta for peak V
      Vpeak.mod.mean <- Vpeak(Qdata = Qmod.mean, timestep = timestep)
      Vpeak.mod.95   <- Vpeak(Qdata = Qmod.95  , timestep = timestep)
      Vpeak.mod.05   <- Vpeak(Qdata = Qmod.05  , timestep = timestep)
      
      deltaVpeak.mean <-  round( (Vpeak.mod.mean - Vpeak.obs) / Vpeak.obs, 3 )
      deltaVpeak.95   <-  round( (Vpeak.mod.95   - Vpeak.obs) / Vpeak.obs, 3 )
      deltaVpeak.05   <-  round( (Vpeak.mod.05   - Vpeak.obs) / Vpeak.obs, 3 )
    
      # Qmax time shifts
      shift.mean <- time.shift(series1 = Qmod.mean, series2 = obs, timestep = timestep)
      shift.95   <- time.shift(series1 = Qmod.95,   series2 = obs, timestep = timestep)
      shift.05   <- time.shift(series1 = Qmod.05,   series2 = obs, timestep = timestep)
    
      # interval scores for total V and peak V
      IS.V     <- quscore(x=c(Vmod.05, Vmod.95, Vobs), conf=0.1); IS.V <- round(IS.V, 0) # [l]
      IS.Vpeak <- quscore(x=c(Vpeak.mod.05, Vpeak.mod.95, Vpeak.obs), conf=0.1); IS.Vpeak <- round(IS.Vpeak, 0) # [l]
    
        
    ret[i,] <- c(id, NS.mean, NS.95, NS.05, 
                 deltaV.mean, deltaV.95, deltaV.05, 
                 round(mean(event.table[,5]), 3), round(sd(event.table[,5]), 3), # dV
                 deltaVpeak.mean, deltaVpeak.95, deltaVpeak.05,
                 round(mean(event.table[,6]), 3), round(sd(event.table[,6]), 3), # dVpeak
                 shift.mean, shift.95, shift.05,
                 round(mean(event.table[,4]), 3), round(sd(event.table[,4]), 3), # time shift
                 round(mean(event.table[,1]), 3), round(sd(event.table[,1]), 3)  # NS
                )
  }
  
  bind.ret <- cbind(ret, VerInd.stat)
  for (i in 2:length(bind.ret[1,])) {          # converts numeric values back to numeric
    bind.ret[,i] <- as.numeric(bind.ret[,i])
  }
  
  # E and sd for all events together; it is more informative to use the boxplot overview
  ret.all <- c(round(mean(event.table.all[,5], na.rm=T), 3), round(sd(event.table.all[,5], na.rm=T), 3),  # dV
               round(mean(event.table.all[,6], na.rm=T), 3), round(sd(event.table.all[,6], na.rm=T), 3),  # dVpeak
               round(mean(event.table.all[,4], na.rm=T), 3), round(sd(event.table.all[,4], na.rm=T), 3),  # time shift
               round(mean(event.table.all[,1], na.rm=T), 3), round(sd(event.table.all[,1], na.rm=T), 3), # NS
               round(sum(as.numeric(bind.ret$n.red.points[-skip])) / sum(as.numeric(bind.ret$n.timesteps[-skip]))  , 3) # reliab
              )
  names(ret.all) <- c(my.stats[c(8, 9, 13, 14, 18, 19, 20, 21)], "reliab")
  
  # data for boxplots A and B
  boxplot.data <- (bind.ret[ , -c(1, 27, 28)])
  if (length(skip) > 0 ) {          # checks which events to ignore for the overall statistics  
    boxplot.data <- boxplot.data[-skip,]  
  }
  
  # data for boxplot C
  all.iterations <- event.table.all[, -c(2,3)]
  
  return(list(bind.ret = bind.ret, VerInd=VerInd, ret.all=ret.all, boxplot.data=boxplot.data, all.iterations=all.iterations))
}

#---------------------------------------------------------------------

plot.Pre.res <- function(pack.dir, prodata, to.plot.list, Ca.res) {
  
  for (j in 1 : length(to.plot.list)) {
    MCMC.propa <- to.plot.list[[j]]$Pre.res$MCMC.propa
    
    pdf( paste(pack.dir, "/4_predict_chains+margs_", names(to.plot.list)[j], ".pdf", sep="") )
      RAMtoplot3 <- Ca.res$RAM
      RAMtoplot3$samples <- MCMC.propa
      RAMtoplot3 <- adaptMCMC::convert.to.coda(RAMtoplot3)
      plot(RAMtoplot3)
    dev.off()
    
    pdf( paste(pack.dir, "/4_predict_cumuplot_", names(to.plot.list)[j], ".pdf", sep="") )
      coda::cumuplot(RAMtoplot3)
    dev.off()
  }
  
  # prepares data for hydrographs 
  hydro.list <- list()
  for (j in 1 : length(to.plot.list)) {
    hydrographsCa  <- list()
    hydrographsPre <- list()
    
    Pre.res    <- to.plot.list[[j]]$Pre.res
    statistics <- to.plot.list[[j]]$statistics[[1]] 
    transf     <- to.plot.list[[j]]$transf
      
    for (i in 1 : length(prodata$Ca)) {
      hydrographsCa[[i]] <- list( timestepCa = Pre.res$bTr.Ca[[i]]$timestepCa, 
                                  data.Ca = prodata$Ca[[i]], 
                                  bct = Pre.res$bTr.Ca[[i]],
                                  transf = transf,
                                  data.source = names(to.plot.list)[j]
                                ) 
    }
      
    for (i in 1 : length(prodata$Pre)) {
      hydrographsPre[[i]] <- list( timestepPre = Pre.res$bTr.Pre[[i]]$timestepPre, 
                                   data.Pre = prodata$Pre[[i]], 
                                   bct = Pre.res$bTr.Pre[[i]],
                                   bind.ret = statistics$bind.ret[i,],
                                   VerInd = statistics$VerInd[[i]],
                                   transf = transf,
                                   data.source = names(to.plot.list)[j]
                                 ) 
    }
    hydro.list[[j]] <- list(hydrographsCa = hydrographsCa, hydrographsPre = hydrographsPre)
  }

  # plots hydrographs
      # Ca events
  for (i in 1 : length(hydro.list)) {
    pdf( paste(pack.dir, "/hydrographs_Ca_", hydro.list[[i]]$hydrographsCa[[1]]$data.source, ".pdf", sep="") )
      for (ii in 1 : length(hydro.list[[i]]$hydrographsCa)) { 
        CaPre.plot.Ca(plotdata1 = hydro.list[[i]]$hydrographsCa[[ii]])
      }
    dev.off()
  }
      # Pre events
  if (length(hydro.list) == 1) {
    pdf( paste(pack.dir, "/00 hydrographs_Pre.pdf", sep="") , height = 3.5,  width = 7)
    for (i in 1 : length(hydro.list[[1]]$hydrographsPre)) { 
      CaPre.plot.Pre(plotdata1 = hydro.list[[1]]$hydrographsPre[[i]])
    }
  }
  
  if (length(hydro.list) == 2) {
    pdf( paste(pack.dir, "/00 hydrographs_Pre.pdf", sep="") , height = 3.5,  width = 7)
    for (i in 1 : length(hydro.list[[1]]$hydrographsPre)) { 
      CaPre.plot.2.Pre(plotdata1 = hydro.list[[1]]$hydrographsPre[[i]], 
                       plotdata2 = hydro.list[[2]]$hydrographsPre[[i]] )
    }
  }
  
  if (length(hydro.list) == 3) {
    pdf( paste(pack.dir, "/00 hydrographs_Pre.pdf", sep="") , height = 6,  width = 7 , fonts = "Times")
    # statistics per event
    for (i in 1 : length(hydro.list[[1]]$hydrographsPre)) { 
      CaPre.plot.3.Pre(plotdata1 = hydro.list[[1]]$hydrographsPre[[i]], 
                       plotdata2 = hydro.list[[2]]$hydrographsPre[[i]],
                       plotdata3 = hydro.list[[3]]$hydrographsPre[[i]] )   
    }
  }
    
  dev.off()
  
  
  # prints statistics overview
 
  for (ii in 1 : length(to.plot.list[[1]]$statistics)) {
    
    # table   
    pdf( paste(pack.dir, "/00_Pre_", names(to.plot.list[[1]]$statistics)[ii], "_stats_table.pdf", sep="") , height = 3.5,  width = 7)
      par(mai=c(0.2, 0.0, 0.2, 0.0))
      plot.new()
      legend("top", inset=0.00, pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=1.2,
             legend = c(paste("overall statistics", sep=""))
      )
      Legend <- c( " ", names(to.plot.list[[1]]$statistics[[ii]]$ret.all))
      for (i in 1 : length(hydro.list)) {
        Legend <- c(Legend,  names(to.plot.list)[[i]], to.plot.list[[i]]$statistics[[ii]]$ret.all)
      }
      
      legend("top", inset=0.15, pch=c(NA, NA, NA), col=c(1,rgb(69,22,198, maxColorValue = 255)),  bty = "n", cex=0.8, ncol=length(hydro.list)+1,
             legend = Legend
      )
    dev.off()
    
    # barplots
    pdf( paste(pack.dir, "/00_Pre_", names(to.plot.list[[1]]$statistics)[ii], "_stats_barplot.pdf", sep="") , height = 6,  width = 7)
      statMatr <- matrix(NA, nrow = 3, ncol = 9); 
      colnames(statMatr) <- names(to.plot.list[[1]]$statistics[[ii]]$ret.all)
      rownames(statMatr) <- names(to.plot.list) 
      for (i in 1:length(to.plot.list)) {
        statMatr[i,] <- to.plot.list[[i]]$statistics[[ii]]$ret.all
      }
      
      par(tcl=-0.5, family="serif", omi=c(0,0,0,0), cex.lab=0.8, cex.axis=0.7 )
      split.screen(c(2, 1))       # splits display into two screens
      split.screen(c(1, 2), screen = 2) # splits the bottom half into 2
      # plot up
        screen(1) 
        par(mai=c(0.6, 0.5, 0.3, 0.1))
        barCenters <- barplot(height = statMatr[,c(1,3,7)], ylim = c(-1, 1),
                              beside = T, axes = T, legend.text = T,
                              names.arg = c(NA,NA,NA),
                              main = "dimensionless statitics: E + sd",
                              border = "black",
                              args.legend = list(title = "rainfall data", 
                                                 x = "bottomright",
                                                 cex = .7))
        
        mtext(side=1, text=names(statMatr[1,c(1,3,7)]), at=barCenters[2,], line = 0)
        
        arrows(barCenters, statMatr[,c(1,3,7)] - 0.5 * statMatr[,c(2,4,8)], 
               barCenters, statMatr[,c(1,3,7)] + 0.5 * statMatr[,c(2,4,8)],
               lwd = 1.5, angle = 90,
               code = 3, length = 0.05)
        
        lines(x=c(1,15), y=c(-1,-1), lty=2); lines(x=c(1,15), y=c(1,1), lty=2); lines(x=c(1,15), y=c(0,0), lty=2)
      
      # plot down left
        screen(3) 
        par(mai=c(0.3, 1.1, 0.6, 0.1))
        barCenters <- barplot(height = as.matrix(statMatr[,5]), 
                              ylim = c(-max(statMatr[,5] + (0.5 * statMatr[,6]))-0.5, max(statMatr[,5] + (0.5 * statMatr[,6]))+0.5),
                              beside = T, 
                              border = "black", axes = TRUE,
                              ylab = "shift(Qmax) [h]",
                              main = "shift(Qmax): E + sd")
        arrows(barCenters, statMatr[,5] - 0.5 * statMatr[,6], 
               barCenters, statMatr[,5] + 0.5 * statMatr[,6],
               lwd = 1.5, angle = 90,
               code = 3, length = 0.05)
       
      # plot down right
        screen(4)
        par(mai=c(0.5, 1.1, 0.6, 0.1))
        barCenters <- barplot(height = as.matrix(statMatr[,9]), 
                              ylim = c(0,1),
                              beside = T, 
                              border = "black", axes = TRUE,
                              main = "reliability")

    close.screen(all = TRUE)
    dev.off()
    
    
    # boxplots A 
    pdf( paste(pack.dir, "/00_Pre_", names(to.plot.list[[1]]$statistics)[ii], "_stats_bxpltA.pdf", sep="") , height = 6,  width = 7)
      nValues <- length(to.plot.list[[1]]$statistics[[ii]]$boxplot.data[1,]) * 
                 length(to.plot.list[[1]]$statistics[[ii]]$boxplot.data[,1])    
      BoxPlot <- data.frame(matrix(NA, ncol=3, nrow= nValues * length(to.plot.list)))
      colnames(BoxPlot) <- c("value", "metric", "datSrc")
      
      for (i in 1:length(to.plot.list)) {
        for (j in 1:length(to.plot.list[[i]]$statistics[[ii]]$boxplot.data[1,])) {
          pos <- (i-1) * nValues + 
                 (j-1) * (length(to.plot.list[[i]]$statistics[[ii]]$boxplot.data[,1])) + 
                 (1 : length(to.plot.list[[i]]$statistics[[ii]]$boxplot.data[,1]))
          BoxPlot$value [pos]  <- to.plot.list[[i]]$statistics[[ii]]$boxplot.data[,j]
          BoxPlot$metric[pos]  <- colnames(to.plot.list[[i]]$statistics[[ii]]$boxplot.data)[j]
          BoxPlot$datSrc[pos]  <- names(to.plot.list)[i]
        }
      }
      
      par(tcl=-0.5, family="serif", omi=c(0,0,0,0), cex.lab=0.8, cex.axis=0.7, cex.main=0.8,
          mai=c(0.3, 0.5, 0.2, 0.1), mgp=c(1.3, 0.6, 0) )
      split.screen(c(2, 1))       # splits display into two screens
      split.screen(c(1, 3), screen = 1) # splits the top    half into 3
      split.screen(c(1, 3), screen = 2) # splits the bottom half into 3
      
      # plot up left
      screen(3) 
        Plot1 <- c( which( BoxPlot$metric == paste("E(", intToUtf8(0x03B4), "V)", sep="") ),        # delta V
                    which( BoxPlot$metric == paste("E(", intToUtf8(0x03B4), "Vpeak)", sep="") )     # delta Vpeak
                   )    
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n',
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7),  ylab = "[-]",
                ylim = c(-max(abs(BoxPlot$value[Plot1])), max(abs(BoxPlot$value[Plot1]))) )
        legend(x = "top", legend = names( to.plot.list)[order(names(to.plot.list))],
               fill = c('white', 'gray80', 'gray30'), cex = 0.6, horiz = T )
        axis(side = 1, at = c(2, 6), line = -0.8, lwd = 0,
             labels = c(paste("E(", intToUtf8(0x03B4), "V)", sep=""), paste("E(", intToUtf8(0x03B4), "Vpeak)", sep="")) )
        lines(x=c(0,8), y=c(0,0), lty=2 )
      close.screen(3)      
      
      # plot up middle
      screen(4)
        Plot1 <- c( which( BoxPlot$metric == paste("sd(", intToUtf8(0x03B4), "V)", sep="") ),        # sd V   
                    which( BoxPlot$metric == paste("sd(", intToUtf8(0x03B4), "Vpeak)", sep="") )     # sd Vpeak
        ) 
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1],  xaxt='n',
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7), ylab = "[-]",
                ylim = c(0, max(abs(BoxPlot$value[Plot1]))),
                main = "A - variations among events")     # title of the pdf file
        axis(side = 1, at = c(2, 6), line = -0.8, lwd = 0,
             labels=c(paste("sd(", intToUtf8(0x03B4), "V)", sep=""), paste("sd(", intToUtf8(0x03B4), "Vpeak)", sep="")) )
        lines(x=c(0,8), y=c(0,0), lty=2)
      close.screen(4) 
      
      # plot up right
      screen(5)
        Plot1 <- c( which( BoxPlot$metric == "reliab" ) )
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n', ylab = "[%]",
                col = c('white', 'gray80', 'gray30') , ylim = c(0, 100))
        axis(side = 1, at = c(2), line = -0.8, lwd = 0,
             labels=c("reliab") )
      close.screen(5) 
      
      # plot down left
      screen(6) 
        Plot1 <- c( which( BoxPlot$metric ==  "E(NSE)" ),        
                    which( BoxPlot$metric == "sd(NSE)" )     
        )    
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n', ylab = "[-]",,
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7) )
        axis(side = 1, at = c(2, 6), line = -0.8, lwd = 0,
             labels=c("E(NSE)", "sd(NSE)") )
        lines(x=c(0,8), y=c(0,0), lty=2)
      close.screen(6) 
      
      # plot down middle
      screen(7) 
        Plot1 <- c( which( BoxPlot$metric ==  "E(shift(Qmax))" ),        
                    which( BoxPlot$metric == "sd(shift(Qmax))" )     
        )    
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n', ylab = "[h]",
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7) )
        axis(side = 1, at = c(2, 6), line = -0.8, lwd = 0,
             labels=c("E(shift(Qmax))", "sd(shift(Qmax))") )
        lines(x=c(0,8), y=c(0,0), lty=2)
      close.screen(7) 
      
      # plot down right
      screen(8) 
        Plot1 <- c( which( BoxPlot$metric == "relABW" ),        
                    which( BoxPlot$metric == "relMIS" )     
        )    
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n', ylab = "[?]",
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7) )
        axis(side = 1, at = c(2, 6), line = -0.8, lwd = 0,
             labels=c("relABW", "relMIS") )
        lines(x=c(0,8), y=c(0,0), lty=2)
      close.screen(8) 
      
    close.screen(all = TRUE)
    dev.off()
    
    
    # boxplots B (   E(dV)   vs   dV(E(Y))   )
    pdf( paste(pack.dir, "/00_Pre_", names(to.plot.list[[1]]$statistics)[ii], "_stats_bxpltB.pdf", sep="") , height = 6,  width = 7)
      nValues <- length(to.plot.list[[1]]$statistics[[ii]]$boxplot.data[1,]) * 
                 length(to.plot.list[[1]]$statistics[[ii]]$boxplot.data[,1])    
      BoxPlot <- data.frame(matrix(NA, ncol=3, nrow= nValues * length(to.plot.list)))
      colnames(BoxPlot) <- c("value", "metric", "datSrc")
      
      for (i in 1:length(to.plot.list)) {
        for (j in 1:length(to.plot.list[[i]]$statistics[[ii]]$boxplot.data[1,])) {
          pos <- (i-1) * nValues + 
            (j-1) * (length(to.plot.list[[i]]$statistics[[ii]]$boxplot.data[,1])) + 
            (1 : length(to.plot.list[[i]]$statistics[[ii]]$boxplot.data[,1]))
          BoxPlot$value [pos]  <- to.plot.list[[i]]$statistics[[ii]]$boxplot.data[,j]
          BoxPlot$metric[pos]  <- colnames(to.plot.list[[i]]$statistics[[ii]]$boxplot.data)[j]
          BoxPlot$datSrc[pos]  <- names(to.plot.list)[i]
        }
      }
      
      par(tcl=-0.5, family="serif", omi=c(0,0,0,0), cex.lab=0.8, cex.axis=0.7, cex.main=0.8,
          mai=c(0.3, 0.5, 0.2, 0.1), mgp=c(1.3, 0.6, 0) )
      split.screen(c(2, 1))       # splits display into two screens
      split.screen(c(1, 3), screen = 1) # splits the top    half into 3
      split.screen(c(1, 3), screen = 2) # splits the bottom half into 3
      
      # plot up left
      screen(3) 
        Plot1 <- c( which( BoxPlot$metric == paste("E(", intToUtf8(0x03B4), "V)", sep="") ),        # delta V
                    which( BoxPlot$metric == paste("E(", intToUtf8(0x03B4), "Vpeak)", sep="") )     # delta Vpeak
                   )     
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n',
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7),  ylab = "[-]",
                ylim = c(-max(abs(BoxPlot$value[Plot1])), max(abs(BoxPlot$value[Plot1]))) )
        legend(x = "top", legend = names( to.plot.list)[order(names(to.plot.list))],
               fill = c('white', 'gray80', 'gray30'), cex = 0.6, horiz = T )
        axis(side = 1, at = c(2, 6), line = -0.8, lwd = 0,
             labels = c(paste("E(", intToUtf8(0x03B4), "V)", sep=""), paste("E(", intToUtf8(0x03B4), "Vpeak)", sep="")) )
        lines(x=c(0,8), y=c(0,0), lty=2 )
      close.screen(3)      
      
      # plot up middle
      screen(4)
        Plot1 <- c( which( BoxPlot$metric == paste(intToUtf8(0x03B4), "V(E(Y))", sep="") ),        # delta V(E(Y)
                    which( BoxPlot$metric == paste(intToUtf8(0x03B4), "Vpeak(E(Y))", sep="") )     # delta Vpeak(E(Y)
                   ) 
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1],  xaxt='n',
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7), ylab = "[-]",
                ylim = c(-max(abs(BoxPlot$value[Plot1])), max(abs(BoxPlot$value[Plot1]))),
                main = "B - E(X) vs. X(E(Y))")     # title of the pdf file
        axis(side = 1, at = c(2, 6), line = -0.8, lwd = 0,
             labels=c(paste(intToUtf8(0x03B4), "V(E(Y))", sep=""), paste(intToUtf8(0x03B4), "Vpeak(E(Y))", sep="")) )
        lines(x=c(0,8), y=c(0,0), lty=2)
      close.screen(4) 
      
      # plot up right
      screen(5)
        Plot1 <- c( which( BoxPlot$metric == "reliab" ) )
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n', ylab = "[%]",
                col = c('white', 'gray80', 'gray30') , ylim = c(0, 100))
        axis(side = 1, at = c(2), line = -0.8, lwd = 0,
             labels=c("reliab") )
      close.screen(5) 
      
      # plot down left
      screen(6) 
        Plot1 <- c( which( BoxPlot$metric ==  "E(NSE)" ),        
                    which( BoxPlot$metric == "NSE(E(Y))" )     
        )    
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n', ylab = "[-]",,
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7) )
        axis(side = 1, at = c(2, 6), line = -0.8, lwd = 0,
             labels=c("E(NSE)", "NSE(E(Y))") )
        lines(x=c(0,8), y=c(0,0), lty=2)
      close.screen(6) 
      
      # plot down middle
      screen(7) 
        Plot1 <- c( which( BoxPlot$metric ==  "E(shift(Qmax))" ),        
                    which( BoxPlot$metric == "shift(Qmax(E(Y)))" )     
        )    
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n', ylab = "[h]",
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7) )
        axis(side = 1, at = c(2, 6), line = -0.8, lwd = 0,
             labels=c("E(shift(Qmax))", "shift(Qmax(E(Y)))") )
        lines(x=c(0,8), y=c(0,0), lty=2)
      close.screen(7) 
      
      # plot down right
      screen(8) 
        Plot1 <- c( which( BoxPlot$metric == "relABW" ),        
                    which( BoxPlot$metric == "relMIS" )     
        )    
        boxplot(BoxPlot$value[Plot1] ~ BoxPlot$datSrc[Plot1] + BoxPlot$metric[Plot1], xaxt='n', ylab = "[?]",
                col = c('white', 'gray80', 'gray30'), at = c(1, 2, 3, 5, 6, 7) )
        axis(side = 1, at = c(2, 6), line = -0.8, lwd = 0,
             labels=c("relABW", "relMIS") )
        lines(x=c(0,8), y=c(0,0), lty=2)
      close.screen(8) 
    
    close.screen(all = TRUE)
    dev.off()
    
    
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

CaPre.plot.Ca <- function(plotdata1) {
  timestepCa  <- plotdata1$timestepCa; data.Ca <- plotdata1$data.Ca; 
  y_B.quant <- plotdata1$bct$bct.yplusB.L1.quant.Ca; 
  y.quant <- plotdata1$bct$bct.y.L1.quant.Ca
  Y.quant <- plotdata1$bct$bct.Y.L1.quant.Ca; Y.samp <- plotdata1$bct$bct.Y.L1.samp.Ca
  out.data.Ca <- data.Ca[[2]]
  transf <- plotdata1$transf
  
  plot(x = timestepCa, y = out.data.Ca[1:nrow(out.data.Ca),2], ylab = "Discharge [l/s]", xlab = "Timestep [h]", 
       ylim = c( 0 , max(Y.quant[(row.names(Y.quant)=="0.95")])*1.1 ) 
  )
  
  polygon(c(timestepCa,rev(timestepCa)),
          c(Y.quant[(row.names(Y.quant)=="0.05")],
            rev(Y.quant[(row.names(Y.quant)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.1), border=NA)
  
  polygon(c(timestepCa,rev(timestepCa)),
          c(y_B.quant[(row.names(y_B.quant)=="0.05")],
            rev(y_B.quant[(row.names(y_B.quant)=="0.95")])), #set the limits (1st and last quantiles)
          col=gray(0.5), border=NA)
  
  polygon(c(timestepCa,rev(timestepCa)),
          c(y.quant[(row.names(y.quant)=="0.05")],
            rev(y.quant[(row.names(y.quant)=="0.95")])),#set the limits (1st and last quantiles)
          col=gray(.8), border=NA)
  
  legend("topleft", legend = paste("90% quantiles for y+B+E (", 
                                   names(transf$par.tr[1]),"=", transf$par.tr[1], ", ",
                                   names(transf$par.tr[2]),"=", transf$par.tr[2],
                                   "; Ca event ", substr(data.Ca[[3]], nchar(data.Ca[[3]])-22, nchar(data.Ca[[3]])-4), ")", sep=""))
  
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