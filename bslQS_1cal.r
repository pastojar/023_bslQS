######################################
## J. Pastorek, Feb 2021
######################################


#######################################
## setting up the environment
print( getwd() )  # make sure the wd is set correctly (it should be on default if inside an R-project)
pack <- devtools::load_all()  # loading the R package  
package <- environmentName( pack$env )   # name of this package


#######################################
## reads FG and RG statistical overview and CML meta data
## ! Decide here about the statistics files (the periods used)  (it defines the rainfall-runoff events)
uni.data <- read.stats(FG.ov.path  = system.file("extdata", "flow_stats_Q2min_(12_tpl+locRGs_smooth+remRGs)+MP1_H2min.csv", package = package),
                       RG.ov.path  = system.file("extdata", "rainfall_stats_locRGs_smooth_(12_tpl+locRGs_smooth+remRGs)+MP1_H2min_2mm_short.csv", package = package)
)
uni.data[["CML_meta"]] <- read.csv( system.file( "extdata", "meta_25xCML_complet.csv", package = package ), sep = ";", stringsAsFactors = F )


#######################################
## selects events (time periods) to work with by their number;   
## ! See e.g.    View(uni.data$RG.overview)  or   View(uni.data$FG.overview)
all.events  <-  ( 1 : length(uni.data$RG.overview[,1]) )

bad.events  <-  c( c(2:9, 79:83, 42), (59), (60:65))                # damaged FG data  || too long  ||  winter  
NA_for_CML  <-  c( 1, 28, 85, 91, 92,     10:18, 35, 46, 53 )   # 10:18, 35, 46, 53 could be simulated, but only with too few CMLs      
NA_for_locRGs <- union( c(19:21, 77:78) , union( c(71, 74, 88) , c(75:78) ) )           #  NA loc RGs  ||  weird locRGs || RG1 stuffed by insect 

good.events <-  setdiff(all.events, union(bad.events, union(NA_for_locRGs, NA_for_CML)) )

periods     <-  uni.data$RG.overview[ good.events,  ] [, c("st", "en") ]

set.seed(288); which_events <- rbinom(n = length(good.events) , size = 1, prob = 0.5)
eventIDsCa    <- as.character( periods$st[ !which_events ] )                     # desired events for CALIBRATION  
eventIDsPre   <- setdiff( as.character( periods$st ) , eventIDsCa )   # desired events for PREDICTION
rownames(periods) [ as.character(periods$st) %in% eventIDsCa ]  <- eventIDsCa
rownames(periods) [ as.character(periods$st) %in% eventIDsPre ] <- eventIDsPre


#######################################
## selects the discharge data to work with
## and reads the data for desired time periods
flow.data      <- read_select_data(rain_data_name = "Q_uncert" , periods = periods )
flow.data.proc <- flow.data[,c("time","id", "MP1", "sd_Q")]






#######################################
## overview of available attenuation data and rainfall data 
## and of methods for processing the data
##
##        SEE    R\\fcie_DataProc.R    !      # sup.rain.data()





#####################################################################################################################
## calibrates using the "Ca" events

#######################################
## reads rainfall data for desired time periods,
## applies the selected processing method and deals with NAs,
scens <- as.character(c())
scens <- c( scens, paste0("read ", "locRGs_smooth__mean3loc--aggregby-min-60" ) )
refRain_Ca <- sup.rain.data( scens = scens, periods = periods[ which( as.character(periods$st) %in% eventIDsCa ) , ] )

scens <- as.character(c())
scens <- c( scens, "read CML14_bslQuantSm" )
CML_bsl_Ca <- sup.rain.data( scens = scens, periods = periods[ which( as.character(periods$st) %in% eventIDsCa ) , ] )


#######################################
## calibrates the chosen WAA model
time_start <- proc.time()

WAA_meth  <- "WAAVal" 
par_names <- c("k", "alp")
par_init  <- c(0.68, 0.34)   # from Valtr et al., 2019
# par_init <-  c(0.697, 0.502)   # from "022_Adj_04_corr60min_cal_perAll" optimization
lower     <- c( 0.001, 0.01 )
upper     <- c( 3, 3 ) 
max.call  <- 5000

par_opt_all <- data.frame(  matrix(vector(), 0, length(par_names))   );  colnames(par_opt_all) <- par_names;   
for ( i_link in colnames( CML_bsl_Ca )[ ! colnames( CML_bsl_Ca ) %in% c("id", "time") ] ) {
  i_link <- strsplit( i_link, "_-_" )[[1]][1]
  
  #######################################
  ## selects data for the given CML
  scens <- as.character(c())
  scens <- c( scens, paste0("CML_bsl_Ca__single-", i_link) )
  CML_bsl_Ca_link <- sup.rain.data( scens = scens )
  
  #######################################
  ## defines the function to be minimized
  WAA_inf <- function(pars) {
    
    scens <- as.character(c())
    
    pars_str <- paste(par_names[1], pars[1], par_names[2], pars[2], sep = "-")
    proc_meth <- paste0(WAA_meth, "-", pars_str, "--AttSpec--AtoR--aggregby-min-60")
    
    scens <- c( scens, paste0( "CML_bsl_Ca_link__", proc_meth ) )
    newRain <- sup.rain.data( scens = scens )
    
    
    #######################################
    ## calculates performance statistics
    out_vec <- c()
    for ( i_col in colnames( newRain[   !colnames(newRain) %in% c("time", "id") ] ) ) {
      
      noNAs <-  !is.na(newRain[i_col]) 
      mod <- newRain[i_col] [ noNAs ] 
      obs <- refRain_Ca[ , ! colnames( refRain_Ca ) %in% c("id", "time") ] [ noNAs ]
      
      out_vec <- c( out_vec,  sqrt( mean( (mod-obs)^2 ) ) )  # RMSE
    }
    
    out <- mean( out_vec )
    
    print(out)
    
    return(out)
  }
  
  
  #######################################
  ## infers the parameters of the WAA model
  par_init <- par_init
  set.seed(42)
  Opt.precal1 <- GenSA::GenSA( par = par_init,
                               fn    = WAA_inf, 
                               lower =  lower,
                               upper =  upper,
                               control= list( max.call = max.call, verbose = T, simple.function = F )
  )
  
  par_opt_all[i_link,] <-t(Opt.precal1$par)
}

time_end   <- proc.time(); time_taken <- time_end - time_start; time_taken

save.image( file = paste0(getwd(), "/outputs/", package, "_1cal.Rdata") )


