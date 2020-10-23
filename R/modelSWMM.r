######################################################################################
# runs SWMM model
#
# J. Pastorek, JAN 2018
######################################################################################

model.swmm <- function(par, L, inp.file)
{
  # o. path to the package files
  
  pack.dir <- strsplit( inp.file, "/swmm")[[1]][1]
                     
  
  # i. sets model parameters
  
  to.run.inp.file <- paste(substr(inp.file, 1, nchar(inp.file)-4), "_par.inp", sep="")
  
  shell( paste( file.path( pack.dir, "extdata", "gawk.exe" ), 
              #" -v imperviousness=",par["mult.imp"],
              #" -v width=", par["mult.wid"],
              " -v slope=", par["mult.slo"],
              " -v nimperv=", par["mult.Nim"],
              " -v simperv=", par["mult.Sim"],
              #" -v sperv=", par["mult.Spe"],
              #" -v pctzero=", par["mult.Pze"],
              #" -v manning=",par["mult.rou"],      
              " -f ", file.path( pack.dir, "awk", "change-3par_v2.awk" ),
              " ", inp.file, " > ",
              to.run.inp.file,
              sep="")
  )
  
  # ii. runs swmm
  
  system(paste( file.path( pack.dir, "extdata", "swmm5.exe" ),
               " ", to.run.inp.file,
               " ", file.path( pack.dir, "swmm", "yM.out" )
  ),                 
  wait=T)
  
  hlp <- shell(paste( file.path( pack.dir, "extdata", "gawk.exe" ),
                     " -f ", file.path( pack.dir, "awk", "get-outlet.awk" ),  
                     " ", file.path( pack.dir,"swmm", "yM.out" ),
                     sep=""
  ),
  intern = TRUE
  )
  
  hlp <- unlist(strsplit(as.character(hlp), ";"))
  
  yM <- data.frame("Qm"   = hlp[c(T,F)],
                   "time" = hlp[c(F,T)])
  
  yM = yM[2:(length(yM[,1])-1), ] # to compensate for the early start (FG - 2*2 min) and late end (FG + 1*2 min) of SWMM
  
  
  # iii. formats swmm output
  
  # creates time layouts for modelled data (Lm)
  if ( as.numeric(yM[2,2])-as.numeric(yM[1,2]) > 0 ) {    # calculates time step of the flow data in hours
    dt.mod = (as.numeric(as.POSIXct(yM[2,2], format="%H:%M:%S",tz="UCT"))
              -as.numeric(as.POSIXct(yM[1,2], format="%H:%M:%S",tz="UCT")) ) / 3600 
  } else {
    dt.mod = (as.numeric(as.POSIXct(yM[3,2], format="%H:%M:%S",tz="UCT"))
              -as.numeric(as.POSIXct(yM[2,2], format="%H:%M:%S",tz="UCT")) ) / 3600  
  }
  origo     = kimisc::hms.to.seconds(format(as.POSIXct(yM[1,2], format="%H:%M:%S",tz="UCT"), format="%H:%M:%S"))/3600#hr
  ult.cal   = (nrow(yM)-1)*dt.mod + origo #hr -- end of event 1
  t.grid  = format( seq(origo, ult.cal, dt.mod), nsmall=6 )
  Lm      = paste("Q",t.grid,sep="_")
  
  matching <- match(L, Lm)
  yM <- as.numeric(as.character(yM$Qm [matching])) # removes data from timesteps when there is NA in observed data
  
  yM <- yM + 6.3   # BASEFLOW (6.3 l/s); implemented here instead of the .inp file !
  
  names(yM) = L
  
  return(yM)
}
