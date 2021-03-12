
load( paste0( getwd(), "/outputs/bsl.QS_1cal.Rdata" ) )
devtools::load_all(".")


#####################################################################################################################
## Evaluates using the "Pre" events 

#######################################
## defines data to be evaluated - individual CMLs
## applies the selected processing method - parameters optimized above
scens <- as.character(c())
scens <- c( scens, "read CML14_bslQuantSm" )
CML_bsl_Pre <- sup.rain.data( scens = scens, periods = periods[ which( as.character(periods$st) %in% eventIDsPre ) , ] )

scens <- as.character(c())
for ( i_link in rownames(par_opt_all) ) {
  
  pars_str  <- paste( par_names[1], par_opt_all[i_link,][1], par_names[2], par_opt_all[i_link,][2], sep = "-")
  proc_meth <- paste0(WAA_meth, "-", pars_str, "--AttSpec--AtoR")
  
  scens <- c( scens, paste0( "CML_bsl_Pre__single-", i_link, "--", proc_meth ) )
}
newRain <- sup.rain.data( scens = scens )


#######################################
## prepares and runs rainfall-runoff simulations
newRain_Q <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                                  data_new  = newRain,
                                  package   = package )
mergedRain <- newRain
mergedRuno <- newRain_Q


#######################################
## defines event subsets for statistics
events.subsets <- list( all = periods$st[ as.character(periods$st) %in% uni.data$RG.overview$id[ which( uni.data$RG.overview$meanRain_Rmax10 > 0 )] ] )

#######################################
## Rainfall-Rainfall  and  Rainfall-Runoff   evaluation
refRain_eval_name <- "locRGs_smooth__mean3loc--aggregby-min-60"
newRain_stats <- Eval_rain_runo( data_rain_ref  = sup.rain.data( scens = paste0("read ", refRain_eval_name), 
                                                                 periods = periods[ eventIDsPre, ] ), 
                                 data_rain_new  = sup.rain.data( scens = "newRain__aggregby-min-60" ),
                                 data_Q         = newRain_Q,
                                 events.subsets = events.subsets ) 




#####################################################################################################################
## Combining the CMLs - random combinations

#######################################
## defines data to be evaluated - CML subsets - random combinations of 2
## applies the selected processing method 
subset_choice <- t( combn(1:16, 2) )  
for ( i_row in 1 : nrow(subset_choice) ) {
  row_i <- subset_choice[i_row,]
  which_cols <- paste( row_i, collapse = "_" ) 
  
  lol <- sup.rain.data( scens = paste0("newRain__subcols_mean-which_cols-", which_cols) )
  colnames(lol)[ !colnames(lol) %in% c("time", "id") ] <- paste0( "twos_", paste0(c(  "#3",  "#4",  "#5",  "#6",  "#7",  "#8",  "#9", "#11", "#12", "#13", "#14", "#15", "#16", "#17", "#18", "#19" )[row_i], collapse = "" )  )
  
  if ( i_row == 1 ) {
    newRainComb_rand2 <- lol
  } else {
    newRainComb_rand2 <- merge(newRainComb_rand2,lol)
  }
}

#######################################
## prepares and runs rainfall-runoff simulations
newRainComb_rand2_Q <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                                            data_new  = newRainComb_rand2,
                                            package   = package )

mergedRain <- merge(mergedRain, newRainComb_rand2)
mergedRuno <- merge(mergedRuno, newRainComb_rand2_Q)



#######################################
## defines data to be evaluated - CML subsets - random combinations of 3
## applies the selected processing method 
subset_choice <- t( combn(1:16, 3) )  
set.seed(42);
subset_choice <- subset_choice[ as.logical( rbinom( n = nrow(subset_choice), size = 1, prob = 0.5 ) ) , ]
for ( i_row in 1 : nrow(subset_choice) ) {
  row_i <- subset_choice[i_row,]
  which_cols <- paste( row_i, collapse = "_" ) 
  
  lol <- sup.rain.data( scens = paste0("newRain__subcols_mean-which_cols-", which_cols) )
  colnames(lol)[ !colnames(lol) %in% c("time", "id") ] <- paste0( "threes_", paste0(c(  "#3",  "#4",  "#5",  "#6",  "#7",  "#8",  "#9", "#11", "#12", "#13", "#14", "#15", "#16", "#17", "#18", "#19" )[row_i], collapse = "" )  )
  
  if ( i_row == 1 ) {
    newRainComb_rand3 <- lol
  } else {
    newRainComb_rand3 <- merge(newRainComb_rand3,lol)
  }
}

#######################################
## prepares and runs rainfall-runoff simulations
newRainComb_rand3_Q <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                                            data_new  = newRainComb_rand3,
                                            package   = package )

mergedRain <- merge(mergedRain, newRainComb_rand3)
mergedRuno <- merge(mergedRuno, newRainComb_rand3_Q)



#####################################################################################################################
## Combining the CMLs - runoff statistics

#######################################
## defines data to be evaluated - CML subsets - NSE
## applies the selected processing method 
for ( i_nCML in 1:nrow(newRain_stats$statistics_runo$all$overview_noEv) ) {
  lol <- newRain[ colnames(newRain) %in% c( "time", "id", 
                                            rownames( newRain_stats$statistics_runo$all$overview_noEv[ order(newRain_stats$statistics_runo$all$overview_noEv$NSE, decreasing = T ) , ] ) [1:i_nCML]
                                           ) ] 
  if ( i_nCML == 1 ) {
    newRainCombNSE <- meanAll(lol)  
    colnames(newRainCombNSE)[3] <- paste0("bestNSE_", i_nCML)
  } else {
    newRainCombNSE[ paste0("bestNSE_", i_nCML) ] <- meanAll(lol)[3]  
  }
}

#######################################
## prepares and runs rainfall-runoff simulations
newRainCombNSE_Q <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                                         data_new  = newRainCombNSE,
                                         package   = package )

mergedRain <- merge(mergedRain, newRainCombNSE)
mergedRuno <- merge(mergedRuno, newRainCombNSE_Q)


#######################################
## defines data to be evaluated - CML subsets - SCC
## applies the selected processing method 
for ( i_nCML in 1:nrow(newRain_stats$statistics_runo$all$overview_noEv) ) {
  lol <- newRain[ colnames(newRain) %in% c( "time", "id", 
                                            rownames( newRain_stats$statistics_runo$all$overview_noEv[ order(newRain_stats$statistics_runo$all$overview_noEv$SCC, decreasing = T ) , ] ) [1:i_nCML]
                                           ) ] 
  if ( i_nCML == 1 ) {
    newRainCombSCC <- meanAll(lol)  
    colnames(newRainCombSCC)[3] <- paste0("bestSCC_", i_nCML)
  } else {
    newRainCombSCC[ paste0("bestSCC_", i_nCML) ] <- meanAll(lol)[3]  
  }
}

#######################################
## prepares and runs rainfall-runoff simulations
newRainCombSCC_Q <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                                         data_new  = newRainCombSCC,
                                         package   = package )

mergedRain <- merge(mergedRain, newRainCombSCC)
mergedRuno <- merge(mergedRuno, newRainCombSCC_Q)



#######################################
## defines data to be evaluated  - CML subsets - dV
## applies the selected processing method 
for ( i_nCML in 1:nrow(newRain_stats$statistics_runo$all$overview_noEv) ) {
  lol <- newRain[ colnames(newRain) %in% c( "time", "id", 
                                            rownames( newRain_stats$statistics_runo$all$overview_noEv[ order(abs(newRain_stats$statistics_runo$all$overview_noEv$dV), decreasing = F ) , ] ) [1:i_nCML]
                                           ) ] 
  if ( i_nCML == 1 ) {
    newRainCombdV <- meanAll(lol)  
    colnames(newRainCombdV)[3] <- paste0("bestdV_", i_nCML)
  } else {
    newRainCombdV[ paste0("bestdV_", i_nCML) ] <- meanAll(lol)[3]  
  }
}

#######################################
## prepares and runs rainfall-runoff simulations
newRainCombdV_Q <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                                        data_new  = newRainCombdV,
                                        package   = package )

mergedRain <- merge(mergedRain, newRainCombdV)
mergedRuno <- merge(mergedRuno, newRainCombdV_Q)



#####################################################################################################################
## Combining the CMLs - other arbitrary subsets

#######################################
## defines data to be evaluated 
## applies the selected processing method 
subset_choice <- list()
subset_choice["arb-3-4-7-8"] <- paste( c(1, 2, 5, 6), collapse = "_" )       # CMLs # 3, 4, 7, 8
subset_choice["arb-3-4-7-8-15"] <- paste( c(1, 2, 5, 12), collapse = "_" )   # CMLs # 3, 4, 7, 8, 15
subset_choice["arb-3-8-12-15"] <- paste( c(1, 6, 9, 12), collapse = "_" )    # CMLs # 3, 8, 12, 15
subset_choice["arb-3-4-5-6-7-8"] <- paste( c(1,2,3,4,5,6), collapse = "_" )  # CMLs # 3, 4, 5, 6, 7, 8
subset_choice["arb-length1"] <- paste( c(1, 2, 3, 4, 5), collapse = "_" )    # CMLs # 3, 4, 5, 6, 7
subset_choice["arb-length2"] <- paste( c(6,7,8,9,10,11), collapse = "_" )    # CMLs # 8, 9, 11, 12, 13, 14 
subset_choice["arb-length3"] <- paste( c(12,13,14,15,16), collapse = "_" )   # CMLs # 15, 16, 17, 18, 19

for ( i_row in 1 : length(subset_choice) ) {
  which_cols <- subset_choice[[i_row]]
  
  lol <- sup.rain.data( scens = paste0("newRain__subcols_mean-which_cols-", which_cols) )
  colnames(lol)[3] <- names(subset_choice)[[i_row]]
  
  if ( i_row == 1 ) {
    newRainArb <- lol
  } else {
    newRainArb <- merge(newRainArb,lol)
  }
}


#######################################
## prepares and runs rainfall-runoff simulations
newRainArb_Q <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                                     data_new  = newRainArb,
                                     package   = package )
mergedRain <- mergedRain[ , !colnames(mergedRain) %in% c( "meanAll_-_newRain_arb1__meanAll", 
                                                          "meanAll_-_newRain_arb2__meanAll", 
                                                          "meanAll_-_newRain_arb3__meanAll" ) ]
mergedRuno <- mergedRuno[ , !colnames(mergedRuno) %in% c( "meanAll_-_newRain_arb1__meanAll", 
                                                          "meanAll_-_newRain_arb2__meanAll", 
                                                          "meanAll_-_newRain_arb3__meanAll" ) ] 

mergedRain <- merge(mergedRain, newRainArb)
mergedRuno <- merge(mergedRuno, newRainArb_Q)






#####################################################################################################################
## plots hydrographs
sup.group.plot.noInf( name = strsplit( names(refRain_Ca)[3], "_-_" )[[1]][2],
                      mod.scens.to.plot = colnames(mergedRain)[ !colnames(mergedRain) %in% c("time", "id") ][c(1,2,20,21)],
                      sup.group.res = mergedRuno, 
                      newRain = mergedRain,
                      out.dir = file.path( getwd(), "outputs") )
dev.off()



save( mergedRain, mergedRuno, file = paste0(getwd(), "/outputs/", package, "_2rr.Rdata") )

