
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
dir_ref <- "D:/OneDrive - České vysoké učení technické/sim_results/023_bslQS/023_bslQS_05_ref"
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


#######################################
## plotting

# boxplot lower and upper boundaries
ylim <- as.data.frame( matrix( nrow = 4,   c( c(0, 1) ,  c(0, 100), c(0.6, 1), 100*c(-1, 1) ) , byrow = T  ), row.names = c("NSE", "RMSE", "SCC", "dV") )

# event classification schemes
event_class_sta <- subset(uni.data$RG.overview, subset = id %in% eventIDsPre , 
                     select = c( "meanRain_Rmax10", "var_locNrem_5",   "var_locNrem_15",  "var_locNrem_30",  "var_locNrem_60"  )  )
event_class_sta$locRG_Qvar <- locRG_NSE_tab$var
rownames(event_class_sta) <- eventIDsPre

# for all event classification schemes, event subsets, and metrics
event_class_col <- event_class_sta[,0]
for ( i_col in 1:ncol(event_class_sta) ) {
  
  if ( colnames(event_class_sta)[i_col] == "locRG_Qvar" ) { 
    event_class_col[ colnames(event_class_sta)[i_col] ] <- NA
    event_class_col[ colnames(event_class_sta)[i_col] ] [order(event_class_sta[i_col]), ] <- nrow(event_class_sta[i_col]):1
  } else {
    event_class_col[ colnames(event_class_sta)[i_col] ] <- round( ( log(event_class_sta[i_col]*100) - log(min(event_class_sta[i_col]*100))  ) / 
                                                                    ( log(max(event_class_sta[i_col]*100 , na.rm = T)) - log(min(event_class_sta[i_col]*100))  )
                                                                  * nrow(event_class_sta[i_col]) 
                                                                 ) 
  }
  event_class_col <- event_class_col + 3  # shifts the colors towards red
  
  if (i_col == 1 ) { subsetets_i <- c("all", "strong", "medium", "light") }
  if (i_col == 2 ) { subsetets_i <- c("all", "var_5_1", "var_5_2", "var_5_3") }
  if (i_col == 3 ) { subsetets_i <- c("all", "var_15_1", "var_15_2", "var_15_3") }
  if (i_col == 4 ) { subsetets_i <- c("all", "var_30_1", "var_30_2", "var_30_3") }
  if (i_col == 5 ) { subsetets_i <- c("all", "var_60_1", "var_60_2", "var_60_3") }
  if (i_col == 6 ) { subsetets_i <- c("all", "locRG_Qvar_60_1", "locRG_Qvar_60_2", "locRG_Qvar_60_3") }
  
  for ( j_subset in subsetets_i ) {
    
    for ( j_metric in c("NSE", "RMSE", "SCC", "dV") ) {
      
      
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
      
      
      # plotting
      png( paste0( out_dir, "/", colnames(event_class_sta)[i_col], "__", j_subset, "_", j_metric , ".png") ,
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
                  col =  fields::tim.colors(n = round(nrow(event_class_col)*1.2), alpha = 0.8) [ event_class_col[ rownames(event_class_col) %in% rownames(data_i) , colnames(event_class_col)[i_col] ]  ] ,
                  pch = 19,
                  cex = 1.8,
                  #ylim = c(y_lim[i_metric, "low"], y_lim[i_metric, "upp"])
          )
        }
        
        points( x = 1:ncol(data_i), y = data_i["noEv",], pch = "-", cex = 5, col = "magenta" )
        
        # mtext(side = 3, line = 0.2, text = paste0( plot_name, " - ", j_subset, " - ", j_metric ), cex = 1.1)
        mtext(side = 3, line = 0.2, text = names(data_to_plot)[i_data], cex = 1.1)
        mtext(side = 1, line = 0.7, text = colnames(data_i), cex =1.1, at = 1:ncol(data_i), las = 2 )
        
        
        abline(h = seq( from = ylim[j_metric,1], to = ylim[j_metric,2], by = (max(ylim[j_metric,]) - min(ylim[j_metric,])) /20 ), 
               col = gray(0.55), lwd = 0.15, lty = 2 ) 
        # abline(h = 0, 
        #        col = gray(0.25), lwd = 0.35, lty = 2 )
        if ( i_data == 1 ) {
          loc_perf <- data_i[ "noEv", "loc" ] 
        }
        abline(h = loc_perf, 
               col = gray(0.25), lwd = 0.5, lty = 2 )
        
        abline(v = c(16, 18)+0.5, col = gray(0.8), lwd = 0.3 )
      }
      
      dev.off()
      
    }
    
  }
  
  
}



