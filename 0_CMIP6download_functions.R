# -------------------------------------------------- #
# CMIP6 data search
# A. Ghosh, R. Hijmans
# -------------------------------------------------- #

library(xml2)
library(httr)

# helper function to find matches in list
getAccess <- function(nd, pattern = ".nc$"){
	  urlp <- xml2::xml_attr(nd, "urlpath")
	  if(!is.na(urlp)){
		dataset <- xml2::xml_attr(nd, "name")
		if (grepl(pattern, dataset)){
		  xd <- xml_contents(nd)
		  xdd <- dplyr::bind_rows(lapply(xml_attrs(xd), function(x) data.frame(as.list(x), stringsAsFactors=FALSE)))
		  xdd <- xdd[,c("name", "value")]
		  xdd <- xdd[complete.cases(xdd),]
		  d <- data.frame(t(xdd$value))
		  colnames(d) <- xdd$name
		  d <- data.frame(file_name = dataset, file_url = urlp, d)
		  return(d)}
	  }
}

findMetaCMIP6 <- function(dots){
  
  # check validity of args supplied
  
  args <- names(dots)
  # args <- c("project", "activity_id", names(dots))
  
  # work in progress: not the final list 
  argsList <- c("offset","limit","activity_id","model_cohort","product","source_id",
                "institution_id","member_id","source_type","nominal_resolution","experiment_id",
                "sub_experiment_id","variant_label","grid_label","table_id","frequency",
                "realm","variable_id","cf_standard_name","data_node","project", "mip_era")
  
  # stop if non-matching argument supplied; pmatch/match.call is a better alternative to avoid errors
  nomatch <- args[!(args %in% argsList)]
  if(length(nomatch) > 0){stop(paste(paste(nomatch, collapse=","), " is not a recognized query parameter"))}
  
  # query structure to create GET request
  query = list(
    type = "Dataset",
    replica = "false",
    latest = "true",
    facets = "mip_era", # which cmip_era
    format = "application/solr+json"
  )
  
  # construct the query
  query <- c(dots, query)
  
  # evaluate query 
  # query <- eval(query)
  
  # llnl node no longer works, updated to dkrz
  # baseurl <- "https://esgf-node.llnl.gov/esg-search/search/"
  baseurl <- "https://esgf-data.dkrz.de/esg-search/search/"
  
  # get request
  request <- httr::GET(url = baseurl, query = query)
  httr::stop_for_status(request)
  
  response <- httr::content(request, as = "text")
  dj <- jsonlite::fromJSON(response)
  return(dj$response$docs)
}

getURLs <- function(d) {

  uud <- list()
  for (i in 1:nrow(d)) {
    ds <- d[i,]
    # read XML metadata
    ur <- unlist(ds$url)[1]
    if (length(ur) == 0) return(uud)
    u <- try(xml2::read_html(ur), silent = TRUE) 
    # if xml content is missing
    if (inherits(u, "try-error")) {
		if (substr(u[1], 1, 24) == "Error in open.connection") return(uud)
		next
	}
    # get dataset nodes
    nds <- xml2::xml_find_all(u, ".//dataset")
    # find the ones with files
    uu <- lapply(nds, getAccess)
	
    #invalid subscript type 'list'
	#uu[sapply(uu, is.null)] <- NULL
    uu <- data.table::rbindlist(uu, fill = TRUE)
    # construct url for data access
    uu$file_url <- file.path("http:/",ds$data_node, "thredds/fileServer", uu$file_url)
    uud[[i]] <- data.frame(ds[rep(seq_len(nrow(ds)), each = nrow(uu)), ], uu)
  }
  return(uud)
}	


getMetaCMIP6 <- function(qargs) {
  
  # cat("Parsing page ", i, "-----------------\n"); flush.console()
  # 10000 is the maximum limit in one query
	xdd <- try(findMetaCMIP6(qargs)) 
	if (NROW(xdd) < 1) return( NULL )

	ldd <- getURLs(xdd)
	
	sdd <- data.table::rbindlist(ldd, fill = TRUE)
	if (NROW(sdd) < 1) return( NULL )

	names(sdd)[names(sdd) == 'size.1'] <- 'file_size'
  # are there any duplicates
	sdd <- unique(sdd, by = "file_id")
	
## no need for this at this point, and it creates errors 	
	#ftime <- try(sapply(strsplit(sdd$file_name, "_"), function(x)grep(".nc", x, value = TRUE)))
	#if (inherits(ftime, "try-error")) {
#		print(sdd$file_name) #null?
#	} else {
#		ftime <- strsplit(gsub(".nc", "", ftime), "-")
#		stime <- sapply(ftime, "[[", 1)
#		etime <- sapply(ftime, "[[", 2)
#		sdd$file_start_date <- as.Date(stime, "%Y%m%d")
#		sdd$file_end_date <- as.Date(etime, "%Y%m%d")
#	}	 
# "file_start_date","file_end_date",
	tokeep <- c("id","version", "activity_id","cf_standard_name","variable","variable_units","data_node","experiment_id", "frequency","index_node","institution_id","member_id","mip_era","file_name", "file_url",
              "nominal_resolution", "number_of_aggregations",
              "number_of_files","pid","project","source_id","source_type","sub_experiment_id", 
              "url","xlink","realm","replica","latest","geo","geo_units","grid","mod_time", "file_size" ,"checksum","checksum_type",
              "grid_label","data_specs_version","tracking_id","citation_url", "further_info_url", "retracted")
	sdd <- sdd[, ..tokeep]
   return(sdd)
}

# complete list of files
# additional cases
# d <- getMetaCMIP6(nominal_resolution = "100km"))
# may be not honored in the search
# d <- try(findMetaCMIP6(offset = offset, limit = limit, experiment_id = "historical"))


# get important parameters from uu: urlPath and dataset
getDataCMIP6 <- function(d, downdir, silent=FALSE){
  if (!silent) print(d$furl); flush.console()
  # specify where to save
  flocal <- file.path(downdir, d$source_id, basename(d$furl)) #uu[1])
  if (file.exists(flocal)) return()
  dir.create(dirname(flocal), FALSE, TRUE)
  # start downloading
  download.file(d$furl, flocal, mode = "wb", quiet=silent)
  # alternate
  # curl::curl_fetch_disk(fileurl, file.path(localdir, dataset))
}


##########################################################################################
checkDownloadStatus <- function(i, idx, downdir){
  d <- idx[i,]
  
  d$download_status <- "FAIL"
  d$localfilechek <- "FAIL"
  
  # localfile exists?
  # fstr <- strsplit(d$file_url, "/CMIP6/|/cmip6/|/cmip6_data/")[[1]][2]
  flocal <- file.path(downdir, basename(d$file_url))
  d$localfile <- flocal
  
  localfilecheck <- file.size(flocal) >= d$file_size  
  
  if(file.exists(flocal) & localfilecheck){
    d$download_status <- "PASS"
    d$localfilecheck <- "PASS"
  }
  return(d)
}


