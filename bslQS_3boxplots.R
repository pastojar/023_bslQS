
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
load( paste0( getwd(), "/outputs/bsl.QS_2eval.Rdata" ) )

RGcalName <- names(refRain_Ca)[3]
if ( grepl( "aggregby-", RGcalName ) ) {
  RGcalName <- sub( "aggregby-", "aggregbykeepLin-", RGcalName  )  
}
# plot_name <- "Cal perLink - loc RGs mean - 60 min"


#######################################
## loading external data
## and merging with the local data
dir_ref <- "D:/OneDrive - České vysoké učení technické/sim_results/023_bslQS/023_bslQS_04_ref"
Rdata_name <- "bsl.QS_2eval.Rdata"
load( file.path(dir_ref, Rdata_name) )

mergedref_stats <- Merge_Eval_rain_runoff( merged_stats, ref_stats )


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
## for all subsets and metrics
ylim <- as.data.frame( matrix( nrow = 4,   c( c(0, 1) ,  c(0, 100), c(0.6, 1), 100*c(-0.5, 0.5) ) , byrow = T  ), row.names = c("NSE", "RMSE", "SCC", "dV") )
for ( j_subset in c("all", "strong", "medium", "light") ) {
  
  for ( j_metric in c("NSE", "RMSE", "SCC", "dV") ) {

    
    data_to_plot <- list()
    
    # single CMLs + ref RGs
    data_hlp <- mergedref_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
      grep( pattern = "single-#" , x = names( mergedref_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
    colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_-_" ) ) [c(T,F)]
    data_hlp["noEv",] <- mergedref_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] [
      grep( pattern = "single-#" , x = rownames( mergedref_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
    
    data_hlp$loc <- c( mergedref_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [["locRGs_smooth__ThPol3"]] ,
                          mergedref_stats$statistics_runo [[j_subset]]$overview_noEv["locRGs_smooth__ThPol3", j_metric] )
    
    data_hlp$cal <- c( mergedref_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [[RGcalName]] ,
                          mergedref_stats$statistics_runo [[j_subset]]$overview_noEv[RGcalName, j_metric] )
    
    data_to_plot[["single"]] <- data_hlp
   
   
    # CML combinations - dV
    data_hlp <- mergedref_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
      grep( pattern = "dV" , x = names( mergedref_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
    colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_" ) ) [c(F,T)]
    data_hlp["noEv",] <- mergedref_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] [ 
      grep( pattern = "dV" , x = rownames( mergedref_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
    data_to_plot[["comb - best dV"]] <- data_hlp
    
    
    # CML combinations - NSE
    data_hlp <- mergedref_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
      grep( pattern = "NSE" , x = names( mergedref_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
    colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_" ) ) [c(F,T)]
    data_hlp["noEv",] <- mergedref_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] [ 
      grep( pattern = "NSE" , x = rownames( mergedref_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
    data_to_plot[["comb - best NSE"]] <- data_hlp
    
    
    # CML combinations - SCC
    data_hlp <- mergedref_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]] [ 
      grep( pattern = "SCC" , x = names( mergedref_stats$statistics_runo [[j_subset]]$overview_ev[[ paste0("table_", j_metric) ]]) ) ]
    colnames( data_hlp ) <- unlist( strsplit( colnames( data_hlp ), "_" ) ) [c(F,T)]
    data_hlp["noEv",] <- mergedref_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] [ 
      grep( pattern = "SCC" , x = rownames( mergedref_stats$statistics_runo [[j_subset]]$overview_noEv[j_metric] ) ) , ]
    data_to_plot[["comb - best SCC"]] <- data_hlp
    

    # plotting
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
                 ylim = as.numeric(ylim[j_metric,]), col = c(col_plot, NA, NA)  )
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


