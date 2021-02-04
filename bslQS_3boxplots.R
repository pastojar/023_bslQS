
load( paste0( getwd(), "/outputs/bsl.QS_2eval.Rdata" ) )
devtools::load_all(".")

plot_name <- "Cal perLink - loc RGs mean - 60 min"
RGrefName <- "rain_-_locRGs_smooth__mean3loc--aggregbykeepLin-min-60"


out_dir <- file.path( getwd(), "outputs", "boxplots" )
if ( dir.exists( out_dir) == F ) {
  dir.create(out_dir)  
}


# # # RG reference
# dir_ref <- "D:/OneDrive - České vysoké učení technické/sim_results/023_bslQS/023_bslQS_03_ref"
# Rdata_name <- "bsl.QS.Rdata"
# load( file.path(dir_ref, Rdata_name) )
# 
# # hlp <- statistics_inf [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ "locRGs_smooth__ThPol3" ]
# # names(hlp) <- "locRGs"
# # data_to_plot <- hlp
# # hlp <- statistics_inf [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [  RGrefName ]
# # names(hlp) <- "remRGs"
# # data_to_plot <- cbind( data_to_plot, hlp )
# 
# 
# 
# # merged_RainRain   <- Merge_Eval_rain_rain(   newRain_RainRain , newRainCombNSE_RainRain, newRainCombSCC_RainRain, newRainCombdV_RainRain, newRainArb_RainRain  )
# # merged_RainRunoff <- Merge_Eval_rain_runoff( newRain_RainRunoff , newRainCombNSE_RainRunoff, newRainCombSCC_RainRunoff, newRainCombdV_RainRunoff, newRainArb_RainRunoff  )






freq <- apply( 
  X =  uni.data$CML_meta[ uni.data$CML_meta$paperNo %in% c(  "#3",  "#4",  "#5",  "#6",  "#7",  "#8",  "#9", "#11", "#12", "#13", "#14", "#15", "#16", "#17", "#18", "#19" ) , c("freqA", "freqB")  ]  , 
  MARGIN = 1, FUN = mean, na.rm = T  
)
names( freq ) <- c(  "#3",  "#4",  "#5",  "#6",  "#7",  "#8",  "#9", "#11", "#12", "#13", "#14", "#15", "#16", "#17", "#18", "#19" )
col_CML <- rep( as.character(NA), length(freq) )
col_CML[ which( freq < 28 ) ] <- gray(0.92)
col_CML[ which( freq > 28 & freq < 35 ) ] <- gray(0.6)
col_CML[ which( freq > 35 ) ] <- gray(0.3)


ylim <- as.data.frame( matrix( nrow = 4,   c( c(0, 1) ,  c(0, 100), c(0.6, 1), 100*c(-0.5, 0.5) ) , byrow = T  ), row.names = c("NSE", "RMSE", "SCC", "dV") )
for ( j_subset in c("all", "strong", "medium", "light") ) {
  
  for ( j_metric in c("NSE", "RMSE", "SCC", "dV") ) {

    
    data_to_plot <- list()
    
    # single CMLs
    data_hlp <- merged_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
      grep( pattern = "single" , x = names( merged_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
    colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_-_" ) ) [c(T,F)]
    data_hlp["noEv",] <- merged_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] [ 
                                 grep( pattern = "single" , x = rownames( merged_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
    data_to_plot[["single"]] <- data_hlp
    # hlp <- c()
    # for ( i_CML in c( "03", "04", "05", "06", "07", "08", "09", "11", "12", "13", "14", "15", "16", "17", "18", "19") ) {
    #   hlp <-  merged_RainRunoff$statistics [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ grep( pattern = paste0( "#", as.numeric(i_CML)) , x = names(merged_RainRunoff$statistics [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
    #   names(hlp) <- i_CML
    #   data_to_plot <- cbind(data_to_plot, hlp); 
    # }
   
    # CML combinations - dV
    data_hlp <- merged_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
      grep( pattern = "dV" , x = names( merged_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
    colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_" ) ) [c(F,T)]
    data_hlp["noEv",] <- merged_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] [ 
      grep( pattern = "dV" , x = rownames( merged_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
    data_to_plot[["comb - best dV"]] <- data_hlp
    # hlp <- c()
    # for ( i_n in 1:16 ) {
    #         i_n_char <- paste0( "best", as.character(i_n ), "_" )
    #         hlp <-  merged_RainRunoff$statistics [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ grep( pattern = i_n_char , x = names(merged_RainRunoff$statistics [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
    #         names(hlp) <- substr( i_n_char, 1, nchar(i_n_char)-1)
    #         data_to_plot <- cbind(data_to_plot, hlp); 
    # }
    
    # CML combinations - NSE
    data_hlp <- merged_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
      grep( pattern = "NSE" , x = names( merged_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
    colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_" ) ) [c(F,T)]
    data_hlp["noEv",] <- merged_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] [ 
      grep( pattern = "NSE" , x = rownames( merged_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
    data_to_plot[["comb - best NSE"]] <- data_hlp
    
    # CML combinations - SCC
    data_hlp <- merged_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
      grep( pattern = "SCC" , x = names( merged_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
    colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_" ) ) [c(F,T)]
    data_hlp["noEv",] <- merged_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] [ 
      grep( pattern = "SCC" , x = rownames( merged_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
    data_to_plot[["comb - best SCC"]] <- data_hlp
    


    png( paste0( out_dir, "/", j_subset, "_", j_metric , ".png") ,
         type="cairo", units = "in", width = 7*length(data_to_plot), height = 4.5, res = 150 )
    
      par(mar=c(3.5,2,2,1), mfrow = c(1, length(data_to_plot)) )
      
      for ( i_data in 1:length(data_to_plot) ) {
        data_i <- data_to_plot[[i_data]]
        
        if ( j_metric == "dV" ) {  # [-] --> [%]       
          data_i <-  data_i * 100 
        }
        
        if ( names(data_to_plot)[i_data] == "single" ) {
          col_plot <- col_CML   
        } else { 
          col_plot <- NA
        }
        
        boxplot( data_i[ !rownames(data_i) %in% "noEv" , ]  , outline = T, range = 1.5,  # default - 1.5 , to extremes - 0
                 names = F, horizontal = F, border = gray(0.2),  cex.axis = 1.2,
                 ylim = as.numeric(ylim[j_metric,]), col = col_plot  )
        points( x = 1:ncol(data_i), y = data_i["noEv",], pch = "-", cex = 4, col = "red" )
        
        # mtext(side = 3, line = 0.2, text = paste0( plot_name, " - ", j_subset, " - ", j_metric ), cex = 1.1)
        mtext(side = 3, line = 0.2, text = names(data_to_plot)[i_data], cex = 1.1)
        mtext(side = 1, line = 0.7, text = colnames(data_i), cex =1.1, at = 1:ncol(data_i), las = 2 )
        
        
        abline(h = seq( from = ylim[j_metric,1], to = ylim[j_metric,2], by = (max(ylim[j_metric,]) - min(ylim[j_metric,])) /20 ), 
               col = gray(0.55), lwd = 0.15, lty = 2 ) 
        
        abline(v = c(2, 18)+0.5, col = gray(0.8), lwd = 0.3 )
      }
    
    dev.off()
    
  }
  
}


