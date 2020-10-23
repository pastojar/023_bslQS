######################################
# J. Pastorek, JUL 2019
######################################


setupSWMMX <- function(eventIDs, flow.data.proc, Urquell, pack.dir)  {
  
  FG.data <- flow.data.proc
  
  prodata <- list(); invalid <- c()
  eventIDs <- as.character(eventIDs)
  for (Ev.id in eventIDs) {   # for every selected event
    
    # selects data for the given event using its ID
    Ev.id <- as.POSIXct(Ev.id, tz="UTC")
    
    
    Ev.FG.data <- match_with_IDs(rainfall_datfr = FG.data, IDs = Ev.id)
    
    
    
    # creates time layouts for observed flow data in hours
    # -------------------------
    Ev.FG.data.t <- as.POSIXct(Ev.FG.data[,1], format="%d/%m/%Y %H:%M",tz="")
    shifts <- as.numeric( (Ev.FG.data.t - Ev.FG.data.t[1]) / 3600)
    origo  <- kimisc::hms.to.seconds( format( Ev.FG.data.t[1] , format="%H:%M:%S" ) ) /3600 # [h]
    t.grid <- format(origo + shifts, nsmall = 6)
    # t = unlist(strsplit(as.character(Ev.FG.data[,1]), ":"))
    # t = matrix(t, ncol = 2, byrow = T)
    # if ( as.numeric(t[2,2])-as.numeric(t[1,2]) > 0 ) {
    #   dt.mod = (as.numeric(t[2,2])-as.numeric(t[1,2]))/60 
    # } else {
    #   dt.mod = (as.numeric(t[3,2])-as.numeric(t[2,2]))/60
    # }
    # 
    # origo     = kimisc::hms.to.seconds(format(as.POSIXct(Ev.FG.data[1,1], format="%d/%m/%Y %H:%M",tz=""), format="%H:%M:%S"))/3600 #[h]
    # ult.cal   = (nrow(Ev.FG.data)-1)*dt.mod + origo # end of event 1 [h]
    # t.grid = format(seq(origo,ult.cal,dt.mod), nsmall=6)
    L <- paste("Q", t.grid, sep="_")

    
    
    # process FG data
    # -------------------------
    Ev.FG.data  = Ev.FG.data[, - which(names(Ev.FG.data)=="id") ]
    
    Ev.FG.data[,2] = Ev.FG.data[,2]*1000   # m^3/s  ---> l/s
    
    Ev.FG.data[,1] = format(Ev.FG.data[,1], "%d/%m/%Y %H:%M")
    
    #Ev.FG.data  =Ev.FG.data[seq(1,nrow(Ev.FG.data),by=2),] # 4 minutes time step
    Ev.FG.data  = Ev.FG.data                                # 2 minutes time step now !
    
    colnames(Ev.FG.data) = c("Date.Time","Flow.L.s.")

    
      
    # modifies inp file for a given event (sets proper time details and adds path to rain time series) 
    # and remembers the path to the inp file
    # -------------------------
    out.name <- format(Ev.id, format = "%Y.%m.%d_%H.%M.%S")
    RG.path  <- paste(pack.dir, "/data/", out.name, "_rainRG", sep="")
    
    Ev.FG.timestep <- as.POSIXct(Ev.FG.data$Date.Time, format="%d/%m/%Y %H:%M", origin="1970-01-01 00:00:00 UTC", tz="UTC")
    
    START_TIME_full <- head(Ev.FG.timestep, n=1)   
    START_TIME_full <- START_TIME_full - 2*2*60 # 2*2 mins before FG data starts
    START_TIME <- format(START_TIME_full, format = "%H:%M:%S")
    START_DATE <- format(START_TIME_full, format = "%m/%d/%Y")
    
    END_TIME_full <- tail(Ev.FG.timestep, n=1)   
    END_TIME_full <- END_TIME_full  + 2*60    # 2 mins after FG data ends  
    END_TIME <- format(END_TIME_full, format = "%H:%M:%S")
    END_DATE <- format(END_TIME_full, format = "%m/%d/%Y")
    
    t.sim     = c(START_DATE, START_TIME, END_DATE, END_TIME, WET_STEP="00:02:00", REPORT_STEP="00:02:00")  # report step 2 mins
    
    inp.path <- paste(pack.dir, "/", out.name, ".inp", sep="")
    shell( paste(system.file("extdata", "gawk.exe", package = Package),
                 " -v START_DATE=", t.sim[1],
                 " -v START_TIME=", t.sim[2],
                 " -v END_DATE=", t.sim[3],
                 " -v END_TIME=", t.sim[4], 
                 " -v WET_STEP=", t.sim[5],
                 " -v REPORT_STEP=", t.sim[6],
                 " -v RAIN_PATH=", RG.path,  # the name of the time series has to be "RAINDAT" (else change time_modif.awk)
                 " -f ", system.file("extdata", "time_modif.awk", package = Package),  
                 " ", Urquell, " > ",
                 inp.path,
                 sep="")
    )
    
    
    # creates output and checks for errors
    # -------------------------------------
    prodata[[as.character(Ev.id)]] <- list(Layout = L, Q_Data = Ev.FG.data, inp.path = inp.path)
  }

  
  return( prodata )
}  


######################################################################################
######################################################################################


setupRainFiles <- function( rain.data.proc, pack.dir )  {
  
  # checks the input format first
  # -------------------------------
  if (length(rain.data.proc[1,]) != 5 && length(rain.data.proc[1,]) != 3) {
    stop("Wrong format of the rainfall data time series.")
  }
  
  if (length(rain.data.proc[1,]) == 5) {
    if ( sum( substr(names(rain.data.proc)[3:5], 1, 3) == c("RG1", "RG2", "RG3") ) != 3 ) {
      stop("Wrong names of the rainfall data time series.")
    } else {
      names(rain.data.proc)[3:5] <- c("RG1", "RG2", "RG3")
    }
  }
  
  if (length(rain.data.proc[1,]) == 3) {
    names(rain.data.proc)[3] <- "RG1"
    rain.data.proc$RG2 <- rain.data.proc$RG1
    rain.data.proc$RG3 <- rain.data.proc$RG1
  }
  
  
  Rain_File_Availabs <- c(); invalid <- c()
  eventIDs <- as.character( unique( rain.data.proc$id ) )
  for (Ev.id in eventIDs) {   # for every given event
    
    # selects data for the given event using its ID
    # -------------------------
    
    if ( is.character(Ev.id) ) {
      Ev.id <- as.POSIXct(Ev.id, tz="UTC")
    }
    Ev.RG.data <- match_with_IDs(rainfall_datfr = rain.data.proc, IDs = Ev.id)
    
    
    # process RG data
    # -------------------------
    
    # checks if there are any data at all or any NAs in the data for the given event in the input rainfall data frame
    if ( length(Ev.RG.data[,1]) < 1 || length(which(is.na(Ev.RG.data))) > 0 ) {
      invalid <- c(invalid, as.character(Ev.id))            
      Rain_File_Availab <- F
      names(Rain_File_Availab) <- as.character(Ev.id)
    } 
    
    else {
      Ev.RG.data.SWMM <- Ev.RG.data
      Ev.RG.data.SWMM[,1] <- format(Ev.RG.data.SWMM[,1], "%m/%d/%Y %H:%M:%S")
      out.name <- format(Ev.id, format = "%Y.%m.%d_%H.%M.%S")
      
      # saves rain time series file for a given event and remembers the path to the file
      # -------------------------
      write.table(Ev.RG.data.SWMM, paste(pack.dir, "/data/", out.name, "_rain.dat", sep=""), col.names=T, row.names=F, sep=";",quote=F)
      rain.path <- c()
      
      for (i in 1:3) {        # the current SWMM model works with 3 RGs
        rain.path[i] <- paste(pack.dir, "/data/", out.name, "_rainRG", i, ".dat", sep="")
        # which rain data to write
        write.table(Ev.RG.data.SWMM[ , c("time", c("RG1", "RG2", "RG3")[i]) ], rain.path[i],  col.names=F, row.names=F, sep=" ",quote=F) 
      }
      
      Rain_File_Availab <- T
      names(Rain_File_Availab) <- as.character(Ev.id)
    }
    
    Rain_File_Availabs <- c(Rain_File_Availabs, Rain_File_Availab)
  }
  
  Rain_File_Availabs_Tab <- data.frame(Rain_File_Availabs = Rain_File_Availabs)
  row.names(Rain_File_Availabs_Tab) <- names(Rain_File_Availabs)
  
  # report message
  # ----------------
  # if ( length(invalid) > 0 ) {            # if there are some time steps with NA rainfall data
  #   print( paste(deparse(substitute(eventIDs)), ": Rainfall data not available for events: ", sep="") ) 
  #   print(invalid)
  # }
  # if ( length(invalid) == 0 ) {
  #   print(paste(deparse(substitute(eventIDs)), ": Rainfall data available for all time steps of all events.", sep=""))
  # }
  
  
  return( Rain_File_Availabs_Tab )
}  


