boxplot.noInf <-  function (pack.dir, statistics,  rainfall_data_name) {

  for (ii in 1: length(statistics)) {
  
    pdf( paste(pack.dir, "/00_bxplt_",  rainfall_data_name, "_", names(statistics)[[ii]], ".pdf", sep="") , height = 6,  width = 7)
    
    nValues <- length(statistics[[ii]]$event.table.mod[1,]) * length(statistics[[ii]]$event.table.mod[,1])     
    BoxPlot <- data.frame(matrix(NA, ncol=2, nrow=nValues ))
    colnames(BoxPlot) <- c("value", "metric")
    
    
    for (j in 1:length(statistics[[ii]]$event.table.mod[1,])) {
      pos <-        (j-1) * (length(statistics[[ii]]$event.table.mod[,1])) + 
                    (1 : length(statistics[[ii]]$event.table.mod[,1]))
      BoxPlot$value [pos]  <- statistics[[ii]]$event.table.mod[,j]
      BoxPlot$metric[pos]  <- colnames(statistics[[ii]]$event.table.mod)[j]
    }
    
    
    par(tcl=-0.5, family="serif", omi=c(0,0,0,0), cex.lab=0.8, cex.axis=0.7, cex.main=0.8,
        mai=c(0.3, 0.5, 0.2, 0.1), mgp=c(1.3, 0.6, 0) )
    split.screen(c(2, 1))       # splits display into two screens
    split.screen(c(1, 2), screen = 2) # splits the bottom half into 2
    
    # plot up 
    screen(1) 
      Plot1 <- c( which( BoxPlot$metric == "dV" ),          # delta V
                  which( BoxPlot$metric == "dVpeak" )      # delta Vpeak
      )     
      boxplot(BoxPlot$value[Plot1] ~ BoxPlot$metric[Plot1], xaxt='n',
              at = c(1, 3),  ylab = "[-]",
              main = paste(rainfall_data_name, "_", names(statistics)[[ii]], sep=""),
              outline = T,  # with outliers!
      )
      # legend(x = "top", legend = names(statistics)[order(names(statistics))],
      #        fill = c('white', 'gray80', 'gray30'), cex = 0.6, horiz = T )
      
      axis(side = 1, at = c(1, 3), labels = c("dV", "dVpeak") )
      
      lines(x=c(0,12), y=c(0,0), lty=2 )
    close.screen(1)      
    
    # plot down left
    screen(3)
      Plot2 <- c( which( BoxPlot$metric == "NS")           
      )     
      boxplot(BoxPlot$value[Plot2] ~  BoxPlot$metric[Plot2], xaxt='n',
              ylab = "[-]",
              outline = T  # with outliers!
      )
      # legend(x = "top", legend = names( statistics)[order(names(statistics))],
      #        fill = c('white', 'gray80', 'gray30'), cex = 0.6, horiz = T )
      
      axis(side = 1, at =1, labels = c( "NSE") )
      
      lines(x=c(0,12), y=c(0,0), lty=2 )
    close.screen(3) 
    
    # plot down right
    screen(4)
      Plot3 <- c( which( BoxPlot$metric == "shift(Qmax)" ) )     
      boxplot(BoxPlot$value[Plot3] ~ BoxPlot$metric[Plot3],  xaxt='n',
              ylab = "[h]", 
              outline = T  # with outliers!
      )
      axis(side = 1, at = 1, labels = c( "shift(Qmax)" ) )
      
      lines(x=c(0,4), y=c(0,0), lty=2)
    close.screen(4) 
      
    close.screen(all = TRUE)
    dev.off()
      
  }

}

