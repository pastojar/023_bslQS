######################################################################################
# 1 cascade reservoir model
#
# J. Pastorek, JUL 2019
# based on a script by D. DelGiudice
######################################################################################

model.1res <- function(par, L, inp.file) 
{
  # model:
  # ------
  #
  #   flow through a 1 cascade reservoir model:
  #
  #   dS1                                             S1
  #   ----  =  Qin - Qout ;     Qin= A*P  ;   Qout =  -- 
  #    dt                                             K
  #
  #
  #   S1   = S1_0*exp(-Dt/K) + Qin[t]*K - Qin[t]*K*exp(-Dt/K)     #//  what is this exactely ???
  # 
  #                     
  #   flow  =   Qout
  #                 
  #
  # state variables:
  # ----------------
  #
  #   S1      cascade state 1
  #
  # parameters:
  # -----------
  #
  #   A       catchment area [m2]
  #   K       reservoir residence time [hr]
  #   Rain    Inp [m/hr]
  #
  #   S1_0    initial cascade state  
  


  
  # set original parameter values:
  # --------------
  par.orig <- c(A = (1/3) * 1.3 * 1000 * 1000,  # [m^2]
                k = 0.3,   # [h]
                S1_0 = 4   # initial cascade state
                )
  
  # mulitply the original parameters:
  # --------------
  param <- par[ names(par.orig) ]
  names(param) <- names(par.orig)
  param[ which(is.na(param)) ] <- par.orig[ which(is.na(param)) ]
  
  par.mod <-  param
  
    
  # get rain data:
  # --------------
  rain.path <- paste(substr(inp.file, 1, nchar(inp.file)-35), "raindata/",
                     substr(inp.file, nchar(inp.file)-22, nchar(inp.file)-4 ), 
                     "_rain.dat", sep = "")
  rain.all  <- read.csv(rain.path, sep=";", header=T, stringsAsFactors = F)   
  
  Rain     <- as.data.frame( rain.all["time"] )
  Rain[, "rain"] <- apply( rain.all[,c("RG1", "RG2", "RG3")], MARGIN = 1, mean  )   #// using the mean of available RGs  
  

  Rain[,"rain"] <- Rain[,"rain"] / 1000 # mm/h -> m/h

    
  # manage time layout:
  # --------------
  Rain$L <- as.POSIXct(Rain[,1], format="%m/%d/%Y %H:%M:%S",tz="UCT")
  shifts <- as.numeric( (Rain$L - Rain$L[1]) )
  origo_rain  <- kimisc::hms.to.seconds( format( Rain$L[1] , format="%H:%M:%S" ) )  
  Rain$L <- origo_rain + shifts
  ult_rain <- as.numeric( Rain$L[ nrow(Rain) ] )

  Dt <-  ( as.numeric(as.POSIXct(Rain$time[2], format="%m/%d/%Y %H:%M:%S",tz="UCT")) - 
           as.numeric(as.POSIXct(Rain$time[1], format="%m/%d/%Y %H:%M:%S",tz="UCT")) ) 


  flow_L <- sysanal.decode(L)$val *60*60
  origo_flow <- as.numeric( flow_L[1] )
  origo      <- min( origo_flow, origo_rain )

  ult_flow <- as.numeric( flow_L[ length(flow_L) ] )
  ult.cal  <- max( ult_flow, ult_rain )
  
  t.grid  <- round( seq(origo, ult.cal + 0.1, Dt) )
  Lm      <- paste( "Q" , format(t.grid / (60*60), nsmall = 6 ) , sep="_" )
  
  
  # STEP I: start simulations a bit earlier and end later
  # burn-in phase
  #---------------
  inustus_temp_orig <- 20   # [time steps]
  inustus_temp_ult  <- 5    # [time steps]
  
  origo_elon   <- origo   - (inustus_temp_orig * Dt)
  ult.cal_elon <- ult.cal + (inustus_temp_ult  * Dt)
  
  hlp <- round( seq(origo_elon, ult.cal_elon + 0.1, by = Dt) )
  Rain.elon <- data.frame( rain = 0, L = hlp ) 
  Rain.elon$rain[ which( as.character(Rain.elon$L) %in% as.character(Rain$L) ) ]  <- Rain$rain

  
  # calculate results:
  # ------------------
  Qin <- Rain.elon; Qin$rain <-  par.mod["A"] * Qin$rain  # m^3/h
  
  SWF <- Qin; names(SWF)[ which( names(SWF) == "rain") ] <- "Q"
  SWF$Q <- linres_anal(Qin = Qin$rain, K = par.mod["k"], S_0 = par.mod["S1_0"], Dt = Dt) # m^3/h
  SWF$Q <- SWF$Q / 3.6 # l/s
   
  flow <- SWF$Q[ which( SWF$L %in% t.grid ) ] + 6.3  # BASEFLOW (6.3 l/s); 
  res  <- as.numeric(flow)
  res[ which(res<0) ]  <- 1e-7 # eliminate neg val
  
  shortLm <- substr(Lm, 1, nchar(Lm)-1 )
  shortL  <- substr(L , 1, nchar(L )-1 )
  matching <- match(shortL, shortLm)
  if ( length( which(is.na(matching)) ) != 0 ) {
    stop("timestamp mismatch")
  }
  
  yM <- as.numeric( as.character(res[matching]) ) # removes data from timesteps when there is NA in observed data
  names(yM) = L
  
  return(yM)
}


######################################################################################
#-----------------------------------------------------------------------------
#
# Qout = Qin - dS/dt        # I think this way is worse than the one below
#
Qout_Qin_dSdt <- function(Qin, K, S_0, Dt) {
  
  S   <- rep(NA,length(Qin))
  SWF <- rep(NA,length(Qin))
  for (t in 1:length(Qin))
  {
    S[t]  <- S_0*exp(-1/K*Dt) + Qin[t]*K - Qin[t]*K*exp(-1/K*Dt) # analytical sol linear reservoir
    
    SWF[t]<- Qin[t] - (S[t]-S_0)/Dt  
    
    S_0    <- S[t] 
  }
  
  return(SWF)
}


######################################################################################
#-----------------------------------------------------------------------------
#
# Analytical solution         # almost the same as above, but I think more precise
#
linres_anal <- function(Qin, K, S_0, Dt) {
  
  S    <- rep( NA, length(Qin) )
  Qout <- rep( NA, length(Qin) )
  for ( t in 1:length(Qin) )
  {
    S[t]  <- S_0*exp(-Dt/K) + Qin[t]*K - Qin[t]*K*exp(-Dt/K) 
    
    Qout[t] <- S[t] / K
    
    S_0    <- S[t] 
  }
  
  return(Qout)
}

