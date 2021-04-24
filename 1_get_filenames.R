

get_meta <- function(model, experiment, freq, path=".", vars=c("pr", "tasmax", "tasmin")) {
	fname <- paste0("CMIP6_", model, "_", experiment, "_", freq, "_", substr(Sys.Date(), 1, 7), ".csv")
	fname <- file.path(path, fname)
	cat(fname, "\n"); flush.console()
	if (file.exists(fname)) return(read.csv(fname))
	
	activity_id = ifelse(experiment=="historical", "CMIP", "ScenarioMIP")
	dh <- list()
	for (i in 1:length(vars)){
		qargs = list(offset = 0,
			 	  limit = 10000,
				  activity_id=activity_id, 
				  experiment_id = experiment,
				  frequency = freq,
				  #member_id = "r1i1p1f1",
				  variable_id = vars[i],
				  source_id = model,
				  mip_era = "CMIP6")
							
		dh[[i]]  <- getMetaCMIP6(qargs)
	}
	# remove any unsuccessful attempts
	dhc <- lapply(dh, function(x){if(inherits(x, "data.table")){return(x)}else{NULL}})
	dd <- data.table::rbindlist(dhc, fill = TRUE)
	if (NROW(dd) == 0) return(NULL)
	x <- as.data.frame(dd)
	for (i in 1:ncol(x)) {
		if (inherits(x[[i]], "list")) {
			x[[i]] <- vapply(x[[i]], paste, collapse = ", ", "")			
		}
	}
	try(write.csv(x, fname, row.names = FALSE))
	x
}


# list of all models; or remove source-id argument from the search functions
models <- c("ACCESS-CM2","ACCESS-ESM1-5","AWI-CM-1-1-MR","AWI-ESM-1-1-LR",
            #"BCC-CSM2-HR","BCC-CSM2-MR","BCC-ESM1",
			"CAMS-CSM1-0","CanESM5",
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


setwd("C:/github/cmip6downscale")
source("0_CMIP6download_functions.R")
dir.create("data", FALSE, FALSE)

all <- list()
for (exper in c(126, 245, 370, 585, "historical")) {
	x <- list()
	for (i in 1:length(models)) {
		cat(i, ": ")
		R.utils::withTimeout(
			try(x[[i]] <- get_meta(models[i], experiment=exper, freq="mon", "data")), 
			timeout = 60
		)
	}
	all <- c(all, x)
}


