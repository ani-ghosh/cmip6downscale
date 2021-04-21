# -------------------------------------------------- #
# CMIP6 data download
# A. Ghosh, R. Hijmans
# -------------------------------------------------- #

source("1b_0_search_CMIP6_functions.R")

vars <- c("pr","tas","tasmax","tasmin")

# list of all models; or remove source-id argument from the search functions
models <- c("ACCESS-CM2","ACCESS-ESM1-5","AWI-CM-1-1-MR","AWI-ESM-1-1-LR",
            "BCC-CSM2-HR","BCC-CSM2-MR","BCC-ESM1","CAMS-CSM1-0","CanESM5",
            "CanESM5-CanOE","CAS-ESM2-0","CESM1-1-CAM5-CMIP5","CESM1-WACCM-SC",
            "CESM2","CESM2-FV2","CESM2-WACCM","CESM2-WACCM-FV2","CIESM","CMCC-CM2-HR4",
            "CMCC-CM2-SR5","CMCC-CM2-VHR4","CMCC-ESM2","CNRM-CM6-1","CNRM-CM6-1-HR",
            "CNRM-ESM2-1","E3SM-1-0","E3SM-1-1","E3SM-1-1-ECA","EC-Earth3","EC-Earth3-AerChem",
            "EC-Earth3-CC","EC-Earth3-LR","EC-Earth3-Veg","EC-Earth3-Veg-LR","EC-Earth3P",
            "EC-Earth3P-HR","EC-Earth3P-VHR","ECMWF-IFS-HR","ECMWF-IFS-LR","FGOALS-f3-H",
            "FGOALS-f3-L","FGOALS-g3","FIO-ESM-2-0","GFDL-AM4","GFDL-CM4","GFDL-CM4C192",
            "GFDL-ESM2M","GFDL-ESM4","GFDL-OM4p5B","GISS-E2-1-G","GISS-E2-1-G-CC",
            "GISS-E2-1-H","GISS-E2-2-G","HadGEM3-GC31-HM","HadGEM3-GC31-LL","HadGEM3-GC31-LM",
            "HadGEM3-GC31-MM","IITM-ESM","INM-CM4-8","INM-CM5-0","INM-CM5-H","IPSL-CM5A2-INCA",
            "IPSL-CM6A-ATM-HR","IPSL-CM6A-LR","IPSL-CM6A-LR-INCA","KACE-1-0-G","KIOST-ESM",
            "MCM-UA-1-0","MIROC-ES2H","MIROC-ES2L","MIROC6","MPI-ESM-1-2-HAM","MPI-ESM1-2-HR",
            "MPI-ESM1-2-LR","MPI-ESM1-2-XR","MRI-AGCM3-2-H","MRI-AGCM3-2-S","MRI-ESM2-0","NESM3",
            "NorCPM1","NorESM1-F","NorESM2-LM","NorESM2-MM","SAM0-UNICON","TaiESM1","UKESM1-0-LL")

varmod <- expand.grid(vars, models)
names(varmod) <- c("vars", "models")

############################################################################################
# historical
dh <- list()

for (i in 1:nrow(varmod)){
  var <- varmod$vars[i]
  model <- varmod$models[i]
  dh[[i]] <- try(getMetaCMIP6(offset = 0,
                              limit = 10000,
                              activity_id="CMIP",
                              experiment_id = "historical",
                              frequency = "day",
                              member_id = "r1i1p1f1",
                              variable_id = var,
                              source_id = model,
                              mip_era = "CMIP6"))
}

# remove any unsuccessful attempts
dhc <- lapply(dh, function(x){if(inherits(x, "data.table")){return(x)}else{NULL}})
dhist <- data.table::rbindlist(dhc, fill = TRUE)


####################################################################################
# future
df <- list()

for (i in 1:nrow(varmod)){
  var <- varmod$vars[i]
  model <- varmod$models[i]
  df[[i]] <- try(getMetaCMIP6(offset = 0,
                              limit = 10000,
                              activity_id="ScenarioMIP",
                              experiment_id = "ssp585",
                              member_id = "r1i1p1f1",
                              frequency = "day",
                              variable_id = var,
                              source_id = model,
                              mip_era = "CMIP6"))
}

# remove any unsuccessful attempts
dfc <- lapply(df, function(x){if(inherits(x, "data.table")){return(x)}else{NULL}})
dfut <- data.table::rbindlist(dfc, fill = TRUE)

# combine both results
dd <- rbind(dhist, dfut)
data.table::fwrite(dd, paste0("data/cmip6_filter_index_", Sys.Date(), ".csv"), row.names = FALSE)

# model selections based on data availability?

############################################################################################
# now download
options(timeout=3600)

# downloader function
getDataCMIP6 <- function(i, idx, downdir, silent=FALSE){
  d <- idx[i,]
  # print something
  if (silent) print(d$file_url); flush.console()
  
  # specify where to save
  # fstr <- strsplit(d$file_url, "/CMIP6/|/cmip6/")[[1]][2]
  flocal <- file.path(downdir, basename(d$file_url))
  d$localfile <- flocal
  
  # where to save
  dir.create(dirname(flocal), FALSE, TRUE)
  
  # should we download?
  if (!file.exists(flocal)){
    # try downloading
    try(download.file(d$file_url, flocal, mode = "wb", quiet=silent))
  }
  return(NULL)
}

##########################################################################################
# change the data directory as needed
downdir <- "~/data/input/climate/CMIP6/daily"

idx <- read.csv("data/cmip6_filter_index_2021-03-24.csv", stringsAsFactors = FALSE)

downloadParallel <- FALSE

if (downloadParallel){
  library(future.apply)
  plan(multiprocess, workers = 12)
  future_lapply(1:nrow(idx), getDataCMIP6, idx, downdir, silent=FALSE, future.seed = TRUE)
} else {
  # download files
  lapply(1:nrow(idx), getDataCMIP6, idx, downdir, silent=FALSE)
}
