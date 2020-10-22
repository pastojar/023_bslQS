
########################################################################################################################

WAA_clim_learn <- function ( data_reference , data_proxy ) {
  
  list_frame <- list()
  for ( i_col in colnames(data_proxy)[ !colnames(data_proxy) %in% c("id", "time") ] ) {
    
    Ameas <- data_proxy[ c("time", i_col) ]
    Ameas_noNA <- Ameas[ ! is.na(Ameas[[i_col]]) , ]
    
    data_reference_noNA <- data_reference[ which( data_reference$time %in% Ameas_noNA$time ) , ]
    data_reference_no0  <- data_reference_noNA[ which( data_reference$time %in% Ameas_noNA$time ) , ] [
      which( data_reference_noNA[ !colnames(data_reference_noNA) %in% c("id", "time") ] != 0 ) , ]
    rate0 <-  1  -  nrow( data_reference_no0 )  /  nrow( data_reference_noNA ) 
    
    Ameas_min <- quantile( x = Ameas_noNA[[i_col]] , probs = rate0, type = 1 )
    Ameas_supercrit <- Ameas_noNA[ which( Ameas_noNA[ i_col ] > Ameas_min ) , ]
    
    CCDF_Ameas <- ecdf( Ameas_supercrit[[i_col]] )
    
    
    
    if ( grepl(pattern = "#", x = i_col) ) {
      i_paperNo <- paste0( "#", strsplit( strsplit( i_col, "#" )[[1]][2] , "_-_" )[[1]][1] )
      alpha <- uni.data$CML_meta$alpha[ which(uni.data$CML_meta$paperNo == i_paperNo) ]
      beta <-  uni.data$CML_meta$beta [ which(uni.data$CML_meta$paperNo == i_paperNo) ]
      L <- ( uni.data$CML_meta$length[ which(uni.data$CML_meta$paperNo == i_paperNo) ] / 1000 )  # [km]
    } else {
      stop("no CML specified explicitely")
    }
    a <- (1/alpha)^(1/beta)
    b <- 1/beta
    
    Ar_teor <- a* data_reference_no0[ , !colnames(data_reference_no0) %in% c("id", "time") ] ^b *L
    
    
    
    list_col <- list( CCDF_Ameas = CCDF_Ameas , Ar_teor =  Ar_teor, Ameas_min = Ameas_min )
    list_frame[[i_col]] <- list_col
  }   
    
  return(list_frame)
}


