
#######################################
## setting up the environment
devtools::load_all(".")

out_dir <- file.path( getwd(), "outputs", "boxplots" )
if ( dir.exists( out_dir) == F ) {
  dir.create(out_dir)  
}


#######################################
## loading the local data
load( paste0( getwd(), "/outputs/bsl.QS_1cal.Rdata" ) )
load( paste0( getwd(), "/outputs/bsl.QS_3stats.Rdata" ) )

refRain_Ca_name <- names(refRain_Ca)[3]
if ( grepl( "aggregby-", refRain_Ca_name ) ) {
  refRain_Ca_name <- sub( "aggregby-", "aggregbykeepLin-", refRain_Ca_name  )  
} 
# plot_name <- "Cal perLink - loc RGs mean - 60 min"


#######################################
## loading external data
## and merging with the local data
dir_ref <- "D:/OneDrive - České vysoké učení technické/sim_results/023_bslQS/023_bslQS_06_ref"
Rdata_name <- "bsl.QS_3stats.Rdata"
load( file.path(dir_ref, Rdata_name) )

mergedref_stats <- Merge_Eval_rain_runoff( merged_stats, ref_stats )

stats_to_plot <- mergedref_stats


#######################################
## setting up boxplot colors
freq <- apply( 
  X =  uni.data$CML_meta[ uni.data$CML_meta$paperNo %in% c(  "#3",  "#4",  "#5",  "#6",  "#7",  "#8",  "#9", "#11", "#12", "#13", "#14", "#15", "#16", "#17", "#18", "#19" ) , c("freqA", "freqB")  ]  , 
  MARGIN = 1, FUN = mean, na.rm = T  
)
names( freq ) <- c(  "#3",  "#4",  "#5",  "#6",  "#7",  "#8",  "#9", "#11", "#12", "#13", "#14", "#15", "#16", "#17", "#18", "#19" )
col_CML <- rep( as.character(NA), length(freq) )
col_CML[ which( freq < 28 ) ] <- gray(0.92)
col_CML[ which( freq > 28 & freq < 35 ) ] <- gray(0.6)
col_CML[ which( freq > 35 ) ] <- gray(0.3)

myCols <- colorRampPalette( colors = fields::tim.colors(n = 10, alpha = 0.8), alpha = T )


#######################################
## plotting

# event classification schemes
event_class_sta <- subset( uni.data$RG.overview, subset = id %in% eventIDsPre , 
                           select = c( "id", "meanRain_Rmax10", "var_locNrem_5",   "var_locNrem_15",  "var_locNrem_60"  )  )
event_class_sta <- merge( event_class_sta, QlocRG_var )
rownames(event_class_sta) <- event_class_sta$id
event_class_sta <- event_class_sta[ !colnames(event_class_sta) %in% "id" ]  
event_class_col <- event_class_sta[,0]


# boxplot lower and upper boundaries
ylim <- as.data.frame( matrix( nrow = 5,   c( c(0, 1), c(0, 1),  c(0, 100), c(0.6, 1), 100*c(-1, 1) ) , byrow = T  ), 
                       row.names = c("NSE", "NNSE", "RMSE", "SCC", "dV") )


# for all event classification schemes, event subsets, and metrics
for ( i_col in 1:ncol(event_class_sta) ) {
  
  # color according to the event type
  event_class_col[ colnames(event_class_sta)[i_col] ] <- round( ( log(event_class_sta[i_col]*100) - log(min(event_class_sta[i_col]*100))  ) / 
                                                                  ( log(max(event_class_sta[i_col]*100 , na.rm = T)) - log(min(event_class_sta[i_col]*100))  )
                                                                * nrow(event_class_sta[i_col]) 
                                                               ) + 1
  
  #event_class_col <- event_class_col + 3  # shifts the colors towards red
  
  if (colnames(event_class_sta)[i_col] == "meanRain_Rmax10") { subsetets_i <- c("all", "strong", "medium", "light") }
  if (colnames(event_class_sta)[i_col] == "var_locNrem_5"  ) { subsetets_i <- c("all", "var_5_1", "var_5_2", "var_5_3") }
  if (colnames(event_class_sta)[i_col] == "var_locNrem_15" ) { subsetets_i <- c("all", "var_15_1", "var_15_2", "var_15_3") }
  if (colnames(event_class_sta)[i_col] == "var_locNrem_60" ) { subsetets_i <- c("all", "var_60_1", "var_60_2", "var_60_3") }
  if (colnames(event_class_sta)[i_col] == "varQ_5"         ) { subsetets_i <- c("all", "varQ_5_1", "varQ_5_2", "varQ_5_3") }
  if (colnames(event_class_sta)[i_col] == "varQ_15"        ) { subsetets_i <- c("all", "varQ_15_1", "varQ_15_2", "varQ_15_3") }
  if (colnames(event_class_sta)[i_col] == "varQ_60"        ) { subsetets_i <- c("all", "varQ_60_1", "varQ_60_2", "varQ_60_3") }
  
  for ( j_subset in subsetets_i ) {
    
    for ( j_metric in c("NSE", "NNSE", "RMSE", "SCC", "dV") ) {
      
      ##############################################################################
      ############# individual CMLs and "smart" combinations
      data_to_plot <- list()
      
      # single CMLs + ref RGs
      data_hlp <- stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
        grep( pattern = "single-#" , x = names( stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
      colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_-_" ) ) [c(T,F)]
      data_hlp["noEv",] <- stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[j_metric] [
        grep( pattern = "single-#" , x = rownames( stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
      
      data_hlp$loc <- c( stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [["locRGs_smooth__ThPol3"]] ,
                         stats_to_plot$statistics_runo [[j_subset]]$overview_noEv["locRGs_smooth__ThPol3", j_metric] )
      
      data_hlp$cal <- c( stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [[refRain_Ca_name]] ,
                         stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[refRain_Ca_name, j_metric] )
      
      data_to_plot[["single"]] <- data_hlp
      
      
      # CML combinations - dV
      data_hlp <- stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
        grep( pattern = "dV" , x = names( stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
      colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_" ) ) [c(F,T)]
      data_hlp["noEv",] <- stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[j_metric] [ 
        grep( pattern = "dV" , x = rownames( stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
      data_to_plot[["comb - best dV"]] <- data_hlp
      
      
      # CML combinations - NSE
      data_hlp <- stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
        grep( pattern = "NSE" , x = names( stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
      colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_" ) ) [c(F,T)]
      data_hlp["noEv",] <- stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[j_metric] [ 
        grep( pattern = "NSE" , x = rownames( stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
      data_to_plot[["comb - best NSE"]] <- data_hlp

      
      # CML combinations - SCC
      data_hlp <- stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
        grep( pattern = "SCC" , x = names( stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
      colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_" ) ) [c(F,T)]
      data_hlp["noEv",] <- stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[j_metric] [ 
        grep( pattern = "SCC" , x = rownames( stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
      data_to_plot[["comb - best SCC"]] <- data_hlp
      
      
      # CML combinations - "arbitrary"
      data_hlp <- stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
        grep( pattern = "arb" , x = names( stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
      colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "arb-" ) ) [c(F,T)]
      data_hlp["noEv",] <- stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[j_metric] [ 
        grep( pattern = "arb" , x = rownames( stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
      data_to_plot[["comb - arbitrary"]] <- data_hlp
      
      
      
      # plotting
      plot_dir <- file.path( out_dir, "ind+smart" )
      if ( dir.exists( plot_dir) == F ) {
        dir.create(plot_dir)  
      }
      png( paste0( plot_dir, "/", colnames(event_class_sta)[i_col], "__", j_subset, "_", j_metric , ".png") ,
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
        
        boxplot( data_i[ !rownames(data_i) %in% "noEv" , ]  , outline = T, range = 0,  # default - 1.5 , to extremes - 0
                 names = F, horizontal = F, border = gray(0.2),  cex.axis = 1.2,
                 ylim = as.numeric(ylim[j_metric,]), col = c(col_plot, NA, NA)  )
        
        for ( i_num in 1:ncol(data_i) ) {
          # xes <- jitter(rep(i_num, length(data_i[ !rownames(data_i) %in% "noEv", i_scen])), amount = 0.2 ) 
          xes <- event_class_col[ rownames(data_i)[ !rownames(data_i )%in% "noEv" ], colnames(event_class_col)[i_col]  ]
          xes <- xes - min(xes) 
          xes <- ( ( xes / max(xes) - 0.5 ) *0.8 ) + i_num
          
          i_scen <- colnames(data_i)[i_num]
          set.seed(1)
          points( x = xes,
                  y = data_i[ !rownames(data_i) %in% "noEv" , i_scen],
                  # col =  fields::tim.colors(n = round(nrow(event_class_col)*1.2), alpha = 0.8) [ event_class_col[ rownames(event_class_col) %in% rownames(data_i) , colnames(event_class_col)[i_col] ]  ] ,
                  col = myCols( n = round(nrow(event_class_col)*1.1) ) [  event_class_col[ rownames(event_class_col) %in% rownames(data_i) , colnames(event_class_col)[i_col] ]  ] ,
                  pch = 19,
                  cex = 1.8,
                  #ylim = c(y_lim[i_metric, "low"], y_lim[i_metric, "upp"])
          )
        }
        
        points( x = 1:ncol(data_i), y = data_i["noEv",], pch = "-", cex = 5, col = "magenta" )
        
        # mtext(side = 3, line = 0.2, text = paste0( plot_name, " - ", j_subset, " - ", j_metric ), cex = 1.1)
        mtext(side = 3, line = 0.2, text = names(data_to_plot)[i_data], cex = 1.1)
        mtext(side = 1, line = 0.7, text = colnames(data_i), cex = 0.9, at = 1:ncol(data_i), las = 2 )
        
        
        abline(h = seq( from = ylim[j_metric,1], to = ylim[j_metric,2], by = (max(ylim[j_metric,]) - min(ylim[j_metric,])) /20 ), 
               col = gray(0.55), lwd = 0.15, lty = 2 ) 
        if ( i_data == 1 ) {
          loc_perf <- data_i[ "noEv", "loc" ] 
        }
        abline(h = loc_perf, 
               col = "magenta", lwd = 0.5, lty = 2 )
        abline(h = 0,
               col = gray(0.25), lwd = 0.35, lty = 2 )
        
        abline(v = c(16, 18)+0.5, col = gray(0.8), lwd = 0.3 )
      }
      
      dev.off()
      
      
      ##############################################################################
      ############# combinations of two CMLs
      data_to_plot <- list()
      
      for ( i_CMlno in 1:16 ) {
        i_CMLpaperNo <- c( "#3", "#4", "#5", "#6", "#7", "#8", "#9", "#11", "#12", "#13", "#14", "#15", "#16", "#17", "#18", "#19" )[i_CMlno]
        
        indices <- grep( paste0( "^", as.character(i_CMlno), "_" ) , colnames(  stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] ) [65:184] ) 
        indices <- c( indices, grep( paste0( "_", as.character(i_CMlno), "_" ) , colnames(  stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] ) [65:184] ) )
        data_hlp <- stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [,65:184] [ indices ]
        
        colnames(data_hlp) <- unlist( strsplit( colnames(data_hlp), "_-_" ) ) [c(T,F)]
        colnames(data_hlp) <- paste( c( "#3", "#4", "#5", "#6", "#7", "#8", "#9", "#11", "#12", "#13", "#14", "#15", "#16", "#17", "#18", "#19" ) [as.numeric(unlist(strsplit(colnames(data_hlp), "_")))][c(T,F)] ,
                                     c( "#3", "#4", "#5", "#6", "#7", "#8", "#9", "#11", "#12", "#13", "#14", "#15", "#16", "#17", "#18", "#19" ) [as.numeric(unlist(strsplit(colnames(data_hlp), "_")))][c(F,T)] ,
                                     sep = "_" )
        
        data_hlp["noEv",] <- stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[j_metric] [65:184,] [ indices ]
        
                
        data_hlp[,i_CMLpaperNo] <- c( stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [[ 
                                        grep( pattern = paste0("single-", i_CMLpaperNo) , x = names( stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]] , 
                                      stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[ 
                                        grep( pattern = paste0("single-", i_CMLpaperNo) , x = rownames( stats_to_plot$statistics_runo [[j_subset]]$overview_noEv ) ), j_metric ] )

        data_hlp$loc <- c( stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [["locRGs_smooth__ThPol3"]] ,
                           stats_to_plot$statistics_runo [[j_subset]]$overview_noEv["locRGs_smooth__ThPol3", j_metric] )
        data_hlp$cal <- c( stats_to_plot$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [[refRain_Ca_name]] ,
                           stats_to_plot$statistics_runo [[j_subset]]$overview_noEv[refRain_Ca_name, j_metric] )
        
        data_to_plot[[ i_CMLpaperNo ]] <- data_hlp
      }
     
      
      # plotting
      plot_dir <- file.path( out_dir, "combOF2" )
      if ( dir.exists( plot_dir) == F ) {
        dir.create(plot_dir)  
      }
      png( paste0( plot_dir, "/", colnames(event_class_sta)[i_col], "__", j_subset, "_", j_metric , ".png") ,
           type="cairo", units = "in", width = 7*4, height = 4.5*4, res = 150 )
      
      par( mar = c(5,2.5,2,0.5), mfrow = c(4,4) )
      
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
        
        boxplot( data_i[ !rownames(data_i) %in% "noEv" , ]  , outline = T, range = 0,  # default - 1.5 , to extremes - 0
                 names = F, horizontal = F, border = gray(0.2),  cex.axis = 1.2,
                 ylim = as.numeric(ylim[j_metric,]), col = c(col_plot, NA, NA)  )
        
        for ( i_num in 1:ncol(data_i) ) {
          # xes <- jitter(rep(i_num, length(data_i[ !rownames(data_i) %in% "noEv", i_scen])), amount = 0.2 ) 
          xes <- event_class_col[ rownames(data_i)[ !rownames(data_i )%in% "noEv" ], colnames(event_class_col)[i_col]  ]
          xes <- xes - min(xes) 
          xes <- ( ( xes / max(xes) - 0.5 ) *0.8 ) + i_num
          
          i_scen <- colnames(data_i)[i_num]
          set.seed(1)
          points( x = xes,
                  y = data_i[ !rownames(data_i) %in% "noEv" , i_scen],
                  # col =  fields::tim.colors(n = round(nrow(event_class_col)*1.2), alpha = 0.8) [ event_class_col[ rownames(event_class_col) %in% rownames(data_i) , colnames(event_class_col)[i_col] ]  ] ,
                  col = myCols( n = round(nrow(event_class_col)*1.1) ) [  event_class_col[ rownames(event_class_col) %in% rownames(data_i) , colnames(event_class_col)[i_col] ]  ] ,
                  pch = 19,
                  cex = 1.8,
                  #ylim = c(y_lim[i_metric, "low"], y_lim[i_metric, "upp"])
          )
        }
        
        points( x = 1:ncol(data_i), y = data_i["noEv",], pch = "-", cex = 5, col = "magenta" )
        
        # mtext(side = 3, line = 0.2, text = paste0( plot_name, " - ", j_subset, " - ", j_metric ), cex = 1.1)
        mtext(side = 3, line = 0.2, text = names(data_to_plot)[i_data], cex = 1.1)
        mtext(side = 1, line = 0.7, text = colnames(data_i), cex =0.9, at = 1:ncol(data_i), las = 2 )
        
        
        abline(h = seq( from = ylim[j_metric,1], to = ylim[j_metric,2], by = (max(ylim[j_metric,]) - min(ylim[j_metric,])) /20 ), 
               col = gray(0.55), lwd = 0.15, lty = 2 ) 
        if ( i_data == 1 ) {
          loc_perf <- data_i[ "noEv", "loc" ] 
        }
        abline(h = loc_perf, 
               col = "magenta", lwd = 0.5, lty = 2 )
        abline(h = 0,
               col = gray(0.25), lwd = 0.35, lty = 2 )
        
        abline(v = c(15, 16, 18)+0.5, col = gray(0.8), lwd = 0.3 )
      }
      
      dev.off()
      
      
    }
    
  }
  
  
}


