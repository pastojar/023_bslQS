
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



#####################################################################################################################
## Combining the CMLs - rainfall statistics

#######################################
## defines data to be evaluated - CML subsets - NSE
## applies the selected processing method 
for ( i_nCML in 1:nrow(newRain_stats$statistics_rain$all$overview_noEv) ) {
  lol <- newRain[ colnames(newRain) %in% c( "time", "id", 
                                            unlist( strsplit( rownames( newRain_stats$statistics_rain$all$overview_noEv[ order(newRain_stats$statistics_rain$all$overview_noEv$NSE, decreasing = T ) , ] ) [1:i_nCML] ,
                                                              "_-_newRain_" ) )[c(T,F)]  ) ] 
  if ( i_nCML == 1 ) {
    newRainCombNSE_rain <- meanAll(lol)  
    colnames(newRainCombNSE_rain)[3] <- paste0("bestNSE_rain_", i_nCML)
  } else {
    newRainCombNSE_rain[ paste0("bestNSE_rain_", i_nCML) ] <- meanAll(lol)[3]  
  }
}

#######################################
## prepares and runs rainfall-runoff simulations
newRainCombNSE_rain_Q <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                                              data_new  = newRainCombNSE_rain,
                                              package   = package )



#######################################
## defines data to be evaluated - CML subsets - SCC
## applies the selected processing method 
for ( i_nCML in 1:nrow(newRain_stats$statistics_rain$all$overview_noEv) ) {
  lol <- newRain[ colnames(newRain) %in% c( "time", "id", 
                                            unlist( strsplit( rownames( newRain_stats$statistics_rain$all$overview_noEv[ order(newRain_stats$statistics_rain$all$overview_noEv$NSE, decreasing = T ) , ] ) [1:i_nCML] ,
                                                              "_-_newRain_" ) )[c(T,F)]  ) ] 
  if ( i_nCML == 1 ) {
    newRainCombSCC_rain <- meanAll(lol)  
    colnames(newRainCombSCC_rain)[3] <- paste0("bestSCC_rain_", i_nCML)
  } else {
    newRainCombSCC_rain[ paste0("bestSCC_rain_", i_nCML) ] <- meanAll(lol)[3]  
  }
}

#######################################
## prepares and runs rainfall-runoff simulations
newRainCombSCC_rain_Q <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                                              data_new  = newRainCombSCC_rain,
                                              package   = package )



#######################################
## defines data to be evaluated  - CML subsets - dV
## applies the selected processing method 
for ( i_nCML in 1:nrow(newRain_stats$statistics_rain$all$overview_noEv) ) {
  lol <- newRain[ colnames(newRain) %in% c( "time", "id", 
                                            unlist( strsplit( rownames( newRain_stats$statistics_rain$all$overview_noEv[ order(newRain_stats$statistics_rain$all$overview_noEv$NSE, decreasing = T ) , ] ) [1:i_nCML] ,
                                                              "_-_newRain_" ) )[c(T,F)]  ) ] 
  if ( i_nCML == 1 ) {
    newRainCombdV_rain <- meanAll(lol)  
    colnames(newRainCombdV_rain)[3] <- paste0("bestdV_rain_", i_nCML)
  } else {
    newRainCombdV_rain[ paste0("bestdV_rain_", i_nCML) ] <- meanAll(lol)[3]  
  }
}

#######################################
## prepares and runs rainfall-runoff simulations
newRainCombdV_rain_Q <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                                             data_new  = newRainCombdV,
                                             package   = package )




#######################################
## defines data to be evaluated - other arbitrary subsets
## applies the selected processing method 
newRain_arb1 <- newRain[ c(1 ,2, 3, 4, 7, 8) ]  # CMLs # 3, 4, 7, 8
newRain_arb2 <- newRain[ c(1 ,2, 3, 4, 5, 6, 7) ]  # CMLs # 3, 4, 5, 6, 7
newRain_arb3 <- newRain[ c(1 ,2, 3, 4, 7, 8, 14) ]  # CMLs # 3, 4, 7, 8, 15

arbs <- 1:3
scens <- as.character(c())
scens <- c( scens, paste0("newRain_arb", arbs, "__meanAll") )
newRainArb <- sup.rain.data(scens = scens)

#######################################
## prepares and runs rainfall-runoff simulations
newRainArb_Q <- PrepNRunRain_runoff( data_flow = flow.data.proc, 
                                     data_new  = newRainArb,
                                     package   = package )





#####################################################################################################################
## merges the rain and runoff data
mergedRain <- merge( merge( merge( merge( newRain, 
                                                            newRainArb ),
                                                     newRainCombdV),
                                              newRainCombNSE), 
                                       newRainCombSCC)

# mergedRain <- merge(merge(merge(merge( merge( merge( merge( newRain, 
#                                                             newRainArb ),
#                                                      newRainCombdV),
#                                               newRainCombNSE), 
#                                        newRainCombSCC),
#                                 newRainCombdV_rain),
#                           newRainCombNSE_rain), 
#                     newRainCombSCC_rain)

mergedRuno <-merge( merge( merge( merge( newRain_Q, 
                                                               newRainArb_Q ),
                                                        newRainCombdV_Q),
                                                 newRainCombNSE_Q), 
                                          newRainCombSCC_Q)

# mergedRuno <- merge( merge( merge( merge( merge( merge( merge( newRain_Q, 
#                                                                newRainArb_Q ),
#                                                         newRainCombdV_Q),
#                                                  newRainCombNSE_Q), 
#                                           newRainCombSCC_Q),
#                                    newRainCombdV_rain_Q),
#                             newRainCombNSE_rain_Q), 
#                      newRainCombSCC_rain_Q)



#####################################################################################################################
## plots hydrographs
sup.group.plot.noInf( name = strsplit( names(refRain_Ca)[3], "_-_" )[[1]][2],
                      mod.scens.to.plot = colnames(mergedRain)[ !colnames(mergedRain) %in% c("time", "id") ][c(1,2,20,21)],
                      sup.group.res = mergedRuno, 
                      newRain = mergedRain,
                      out.dir = file.path( getwd(), "outputs") )
dev.off()



save( mergedRain, mergedRuno, file = paste0(getwd(), "/outputs/", package, "_2rr.Rdata") )

