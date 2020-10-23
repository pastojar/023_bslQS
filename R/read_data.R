######################################
# J. Pastorek, JUL 2017
######################################


read.stats <- function(FG.ov.path, RG.ov.path) {
  
  # reads FG stats
  FG.overview <- read.csv(FG.ov.path, sep=";", header=T)
  
  TQmax_names <- c()
  for (i in 1:length(names(FG.overview))) {
    if ( length(strsplit(names(FG.overview)[i], "_")[[1]]) == 2 ) {
      if ( strsplit(names(FG.overview)[i], "_")[[1]][2] == "TQmax" ) {
        TQmax_names <- c(TQmax_names, names(FG.overview)[i])  
      }
    }
  }
  
  for (i in c("st", "en", TQmax_names)) {
    FG.overview[,i] <- as.POSIXct( FG.overview[,i], origin="1970-01-01 00:00:00 UTC", tz="UTC" )
  }
  FG.overview$id <- as.character(FG.overview$id)
  
  # reads RG stats
  RG.overview <- read.csv(RG.ov.path, sep=";", header=T)
  for (i in c("st", "en")) {
    RG.overview[i] <- as.POSIXct( RG.overview[,i], origin="1970-01-01 00:00:00 UTC", tz="UTC" )
  }
  RG.overview$id <- as.character(RG.overview$id)
 
  # checks if statistscs for rain and flow data were created using the same events
  matchRGFG <- match(RG.overview$id, FG.overview$id)
  if ( length(which(is.na(matchRGFG))) != 0 ) {
    stop("RG FG data timestamp mismatch")
  }
  
  
  return(list(FG.overview = FG.overview, RG.overview = RG.overview))
}

######################################################################################
######################################################################################

read_select_data <- function (rain_data_name, periods) {
  # this function reads standardized rainfall data and matches them with periods of interest
  #
  # inputs:
  #   rain_data_name - the name of the RData file containing the R data frame
  #   periods - a data frame with desired "st" and "en" times
  #
  # outputs:
  #   rain.data - a data frame of rainfall data from desired periods with IDs for every time step
  
  
  
  rain.dat.path <- file.path( getwd(), "inst", "extdata", paste0(rain_data_name, ".RData") )
  
  # loads the rainfall data ( a data frame named "R")
  load( rain.dat.path ) 
  
  # matches the data and adds IDs
  if ( as.character(periods) == "whole" ) {
    rain.data <- R
  } else {
    rain.data <- match_with_periods( raindata = R, periods = periods )  
  }
  
  # renaming the columns
  
  # 1)  for CML data
  if ( substr(rain_data_name, 1, 3) == "CML" ) {
    cols.to.replace <- colnames(rain.data) %in% uni.data$CML_meta$id[!is.na(uni.data$CML_meta$paperNo)]
    if ( length(which(cols.to.replace == TRUE)) != 0 ) {
      # for data combining both CML channels
      colnames(rain.data)[ cols.to.replace ] [ match( uni.data$CML_meta$id[!is.na(uni.data$CML_meta$paperNo)], colnames(rain.data)[ cols.to.replace ] ) ] <- uni.data$CML_meta$paperNo[!is.na(uni.data$CML_meta$paperNo)]   
    } else {  
      # for single channels
      for ( i_col in 1 : length(colnames(rain.data)) ) {
        if ( length( grep( colnames(rain.data)[i_col],  uni.data$CML_meta$id[!is.na(uni.data$CML_meta$paperNo)] ) ) > 0 ) {
          cols.to.replace[i_col] <- TRUE
        }
      }
      if ( length(which(cols.to.replace == TRUE)) != 0 ) {
        for ( i_colnames in colnames(rain.data)[ !colnames(rain.data) %in% c("time", "id") ] ) {
          #cmlID   <- uni.data$CML_meta$id[!is.na(uni.data$CML_meta$paperNo)] [ grep( i_colnames,  uni.data$CML_meta$id[!is.na(uni.data$CML_meta$paperNo)] ) ]
          for ( i_id in uni.data$CML_meta$id ) {
            if ( i_colnames %in% strsplit( i_id, "_" )[[1]] ) {
              cmlID <- i_id    
              hlp <- grep( i_colnames, strsplit( i_id, "_" )[[1]] )
              break()
            }
          }
          if ( cols.to.replace[ which( colnames(rain.data) == i_colnames ) ] == TRUE ) {
            paperNo <- uni.data$CML_meta$paperNo[!is.na(uni.data$CML_meta$paperNo)] [ which( uni.data$CML_meta$id[!is.na(uni.data$CML_meta$paperNo)] == cmlID ) ]
            colnames(rain.data)[ which(colnames(rain.data)==i_colnames) ] <- paste(paperNo, c("A", "B")[hlp], sep = "_")  
          } else {
            colnames(rain.data)[ which(colnames(rain.data)==i_colnames) ] <- paste(cmlID, c("A", "B")[hlp], sep = "_")  
          }
        }
      }
    }
    
    # check for data combining both CML channels
    if ( sum( c("415_416", "78_NA",   "NA_384",  "385_386", "59_60",   "545_547", "79_105",  "57_58",   "76_77",   "74_NA",   "73_NA",   "63_64",   "67_68",   "69_70",   "NA_56",
                "71_72",   "375_376", "65_66",   "61_62" ) %in% colnames(rain.data) ) > 0 ) {
      cols.to.replace <- colnames(rain.data) %in% uni.data$CML_meta$id[!is.na(uni.data$CML_meta$paperNo)]
      colnames(rain.data)[ cols.to.replace ] [ match( uni.data$CML_meta$id[!is.na(uni.data$CML_meta$paperNo)], colnames(rain.data)[ cols.to.replace ] ) ] <- uni.data$CML_meta$paperNo[!is.na(uni.data$CML_meta$paperNo)] 
    }
  } 
  # 2)  for other data 
  else {  

    if ( length( colnames(rain.data)[ !colnames(rain.data) %in% c("time", "id") ] )  == 1 ) {
      names(rain.data)[ !names(rain.data) %in% c("time", "id") ] <- rain_data_name
    }
    
  }

  return(rain.data)
}

