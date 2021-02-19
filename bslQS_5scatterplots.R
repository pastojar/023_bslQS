
#######################################
## setting up the environment
devtools::load_all(".")

out_dir <- file.path( getwd(), "outputs", "scatterplots" )
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
## color according to the event type
locRG_NSE_tab$col [order(locRG_NSE_tab$var)] <- nrow(locRG_NSE_tab):1
# locRG_NSE_tab$col <- log( -locRG_NSE_tab$var +1 )
# locRG_NSE_tab$col <- nrow(locRG_NSE_tab) * ( ( locRG_NSE_tab$col - min(locRG_NSE_tab$col ) ) / max( locRG_NSE_tab$col - min(locRG_NSE_tab$col )  ) )  + 1


#######################################
## plotting
# myCols <- colorRampPalette( colors = RColorBrewer::brewer.pal(11, "Spectral" )[11:1] , alpha = 0.9  )
# myCols <- colorRampPalette( colors = c("red", "yellow", "green", "blue")[4:1] , alpha = 0.9  )
# alpha = 0.1myCols <- colorRampPalette( colors = c( rgb(1,0,0,1), rgb(1,0,0,0) ), alpha = T )
myCols <- colorRampPalette( colors = fields::tim.colors(n = 10, alpha = 0.8), alpha = T )

# for ( i_cml in colnames(stats_to_plot$FlowData)[ grepl( "single-#", colnames(stats_to_plot$FlowData) ) ] ) {
for ( i_cml in colnames(stats_to_plot$FlowData)[ 6:ncol(stats_to_plot$FlowData) ] ) {
  
  data_i <- stats_to_plot$FlowData[ , c("id", "Qobs", i_cml) ]

  png( paste0( out_dir, "/", i_cml, "__vs__Qobs",  ".png") ,
       type="cairo", units = "in", width = 7, height = 7, res = 250 )
  
  par(mar=c(3.5,2,2,1), bg = grey(0.8))
  
    plot( x = data_i$Qobs ,
          y = data_i[,i_cml] ,
          log = "xy" ,
          xlim = c(5.5, max(Qobs, na.rm = T) ) ,
          ylim = c(5.5, max(Qobs, na.rm = T) ) ,
          type = "n"  )
    
    for ( i_ev in unique(stats_to_plot$FlowData$id) ) {
      data_ij <- match_with_IDs(data_i, i_ev)
      
      points( x = data_ij$Qobs ,
              y = data_ij[,i_cml] ,
              #pch = 16, cex = 0.4,
              pch = 19, cex = 0.2,
              # col =  [ locRG_NSE_tab$col[ rownames(locRG_NSE_tab) %in% i_ev ]  ] 
              col = myCols( n = round(nrow(locRG_NSE_tab)*1.1) ) [ locRG_NSE_tab$col[ rownames(locRG_NSE_tab) %in% i_ev ]  ]
      )
    }
    
    lines( x = c(1,2000) ,
           y = c(1,2000) ,
           col = gray(0.5) ,
           lwd = 1.5 )
    
  dev.off()
}







