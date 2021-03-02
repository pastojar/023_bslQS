
#######################################
## loading the local data
load( paste0( getwd(), "/outputs/bsl.QS_1cal.Rdata" ) )
load( paste0( getwd(), "/outputs/bsl.QS_2rr.Rdata" ) )
devtools::load_all(".")


#######################################
## loading external data
dir_ref <- "D:/OneDrive - České vysoké učení technické/sim_results/023_bslQS/023_bslQS_06_ref"
Rdata_name <- "bsl.QS_2rr.Rdata"
load( file.path(dir_ref, Rdata_name) )



#######################################
## defines event subsets for statistics
RG_overview_Pre <- uni.data$RG.overview[ uni.data$RG.overview$id %in% eventIDsPre , ]
eventIDs <- as.POSIXct( RG_overview_Pre$id , tz = "UTC")

# Rmax10
events.subsets <- list(all       = eventIDs[ as.character(eventIDs) %in% RG_overview_Pre$id[ which( RG_overview_Pre$meanRain_Rmax10 > 0 )  ] ] ,
                       strong    = eventIDs[ as.character(eventIDs) %in% RG_overview_Pre$id[ which( RG_overview_Pre$meanRain_Rmax10 > 12 ) ] ] ,
                       strongest = eventIDs[ as.character(eventIDs) %in% RG_overview_Pre$id[ which( RG_overview_Pre$meanRain_Rmax10 > 20 ) ] ] ,
                       medium    = eventIDs[ as.character(eventIDs) %in% intersect( RG_overview_Pre$id[ which( RG_overview_Pre$meanRain_Rmax10 < 12 ) ],
                                                                                        RG_overview_Pre$id[ which( RG_overview_Pre$meanRain_Rmax10 > 5  ) ] ) ] ,
                       light     = eventIDs[ as.character(eventIDs) %in% RG_overview_Pre$id[ which( RG_overview_Pre$meanRain_Rmax10 < 5 ) ] ]
)

# rainfall variability
for ( i_agg in c(5, 15, 60) ) {
  events.subsets [[ paste0( "var_", i_agg, "_1" ) ]]  <-  eventIDs[ as.character(eventIDs)   %in%
                                                                        RG_overview_Pre$id[ which( RG_overview_Pre[ , paste0("var_locNrem_", i_agg) ] <=
                                                                                                          quantile( RG_overview_Pre[ , paste0("var_locNrem_", i_agg) ], 0.33 ) )  ] ]

  events.subsets [[ paste0( "var_", i_agg, "_3" ) ]]  <-  eventIDs[ as.character(eventIDs)   %in%
                                                                        RG_overview_Pre$id[ which( RG_overview_Pre[ , paste0("var_locNrem_", i_agg) ] >=
                                                                                                          quantile( RG_overview_Pre[ , paste0("var_locNrem_", i_agg) ], 0.66 ) )  ] ]

  events.subsets [[ paste0( "var_", i_agg, "_2" ) ]]  <-  eventIDs[ ! eventIDs %in% c( events.subsets [[ paste0( "var_", i_agg, "_1" ) ]],
                                                                                           events.subsets [[ paste0( "var_", i_agg, "_3" ) ]] ) ]
}


# # variability among the runoff from the local RGs measured by NNSE
# locRG_stats <- list()
# for ( i_locRG in 1:3 ) {
#   locRG_names <- colnames(refRain_Q)[ grep( "locRGs_smooth__keep3--single-RG", colnames(refRain_Q) ) ] [c(T,F,F)]
#   i_locRG_name <- locRG_names[i_locRG]
# 
#   hlp_runo <- refRain_Q[ , c("time", "timestamp", "id", locRG_names[ !locRG_names %in%  i_locRG_name ] ) ]
#   hlp_runo["sd_Qobs"] <- NA
#   hlp_runo["Qobs"] <-  refRain_Q[,i_locRG_name]
# 
#   locRG_stats[[i_locRG_name]] <- simple.stats.sup.group( sup.group.res = hlp_runo )
# 
#   hlp <- locRG_stats[[i_locRG]][["overview_ev"]][["table_NNSE"]]
#   colnames(hlp) <- paste0(i_locRG_name, "__", locRG_names[ !locRG_names %in%  i_locRG_name ]  )
#   if ( i_locRG == 1 ) {
#     locRG_NNSE_tab <- hlp
#   } else {
#     locRG_NNSE_tab <- cbind( locRG_NNSE_tab, hlp )
#   }
# }
# locRG_NNSE_tab$var <- apply( locRG_NNSE_tab, 1, mean )
# 
# events.subsets [[ paste0( "Qvar_NNSE_1" ) ]] <- periods$st[ as.character(periods$st) %in% rownames( locRG_NNSE_tab[ order(locRG_NNSE_tab$var, decreasing = T), ] ) [ 1:8 ] ]
# events.subsets [[ paste0( "Qvar_NNSE_2" ) ]] <- periods$st[ as.character(periods$st) %in% rownames( locRG_NNSE_tab[ order(locRG_NNSE_tab$var, decreasing = T), ] ) [ 9:16 ] ]
# events.subsets [[ paste0( "Qvar_NNSE_3" ) ]] <- periods$st[ as.character(periods$st) %in% rownames( locRG_NNSE_tab[ order(locRG_NNSE_tab$var, decreasing = T), ] ) [ 17:24 ] ]


# variability among the runoff from the local and 3 remote RGs measured by the normalized variance
locRG_names <- colnames(refRain_Q)[ grep( "locRGs_smooth__keep3--single-RG", colnames(refRain_Q) ) ] [c(T,F,F)]
locRG_names <- c( locRG_names, colnames(refRain_Q)[ grep( "10_-_remRGs_all__single-10", colnames(refRain_Q) ) ] [c(T,F,F)] )
locRG_names <- c( locRG_names, colnames(refRain_Q)[ grep( "13_-_remRGs_all__single-13", colnames(refRain_Q) ) ] [c(T,F,F)] )
locRG_names <- c( locRG_names, colnames(refRain_Q)[ grep( "22_-_remRGs_all__single-22", colnames(refRain_Q) ) ] [c(T,F,F)] )

QlocRG_var <- data.frame( id = unique(refRain_Q$id), stringsAsFactors = F  )
for ( i_min in c(5, 15, 60) ) {
  par <- i_min; names(par) <- "min"
  
  QlocRG_aggreg <- aggregby( refRain_Q[  c("time", "timestamp", "id", locRG_names) ], par )
  for ( i_ev in unique(QlocRG_aggreg$id) ) {
    data_ev <- subset( QlocRG_aggreg, subset = as.character(id) == i_ev, 
                       select = colnames(QlocRG_aggreg)[ !colnames(QlocRG_aggreg) %in% c("id", "time", "timestamp") ] )
    
    mean_var <- round( mean( apply( X = data_ev, MARGIN = 1, FUN = sd, na.rm = T ) / 
                               apply( X = data_ev, MARGIN = 1, FUN = mean, na.rm = T ) , 
                             na.rm = T ) , 
                       digits = 3 )
    QlocRG_var [ QlocRG_var$id %in% i_ev , paste0( "varQ_", i_min ) ] <-  mean_var
  }  
  
  events.subsets [[ paste0( "varQ_", i_min, "_1" ) ]]  <- eventIDs[ as.character(eventIDs)   %in%
                                                                      QlocRG_var$id[ which( QlocRG_var[ , paste0("varQ_", i_min) ] <=
                                                                                              quantile( QlocRG_var[ , paste0("varQ_", i_min) ], 0.33 ) ) ]  ]
  
  events.subsets [[ paste0( "varQ_", i_min, "_3" ) ]]  <- eventIDs[ as.character(eventIDs)    %in%
                                                                      QlocRG_var$id[ which( QlocRG_var[ , paste0("varQ_", i_min) ] >=
                                                                                              quantile( QlocRG_var[ , paste0("varQ_", i_min) ], 0.66 ) ) ]  ] 
  
  events.subsets [[ paste0( "varQ_", i_min, "_2" ) ]]  <- eventIDs[ !as.character(eventIDs) %in%  c( as.character(events.subsets [[ paste0( "varQ_", i_min, "_1" ) ]]),
                                                                                                     as.character(events.subsets [[ paste0( "varQ_", i_min, "_3" ) ]]) ) ]
}




#######################################
## Rainfall-Rainfall  and  Rainfall-Runoff   evaluation
refRain_eval_name <- "locRGs_smooth__mean3loc--aggregby-min-60"
merged_stats <-  Eval_rain_runo( data_rain_ref  = sup.rain.data( scens = paste0("read ", refRain_eval_name), 
                                                                 periods = periods[ eventIDsPre, ] ), 
                                 data_rain_new  = sup.rain.data( scens = "mergedRain__aggregby-min-60" ),
                                 data_Q         = mergedRuno,
                                 events.subsets = events.subsets ) 



#####################################################################################################################
## exports the data
for (j in names(events.subsets)) {
  write.table( merged_stats$statistics_rain [[ j ]]$overview_noEv , paste0(getwd(), "/outputs/", "/stat_rain_noEv_", j ,".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
  write.table( merged_stats$statistics_runo [[ j ]]$overview_noEv , paste0(getwd(), "/outputs", "/stat_runoff_noEv_", j ,".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
}
for (j in names(merged_stats$statistics_rain$all$overview_ev)) {
  write.table( merged_stats$statistics_rain$all$overview_ev [[j]] , paste0(getwd(), "/outputs", "/stat_rain_perEv_",  j, ".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
  write.table( merged_stats$statistics_runo$all$overview_ev [[j]] , paste0(getwd(), "/outputs", "/stat_runoff_perEv_",  j, ".csv"),
               sep = ";", col.names = NA, row.names = TRUE )
}


save( merged_stats, QlocRG_var, file = paste0(getwd(), "/outputs/", package, "_3stats.Rdata") )

