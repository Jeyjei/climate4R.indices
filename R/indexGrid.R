#     indexGrid.R Climate Indices in Climate4R
#
#     Copyright (C) 2018 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see <http://www.gnu.org/licenses/>.

#' @title Climate Indices in Climate4R
#' @description Calculation of indices.
#' @param tn A climate4R dataset of daily minimum temperature (degrees C)
#' @param tx A climate4R dataset of daily maximum temperature (degrees C)
#' @param tm A climate4R dataset of daily mean temperature (degrees C)
#' @param pr A climate4R dataset of daily precipitation (mm)
#' @param any A climate4R dataset of any variable.
#' @param baseline Optional climate4R dataset. Only used if \code{index.code = "P"}, for calculating the relevant percentiles.
#' @param index.code Character string, indicating the specific code of the index (see Details).
#' @param time.resolution Output time resolution. Choices are "month", "year" (default) and "climatology".
#' @param ... Optional. A list of arguments internally passed to the functions displayed by \code{\link{indexShow}}.
#' @template templateParallelParams
#' @import transformeR
#' @importFrom parallel stopCluster
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom utils head
#' @details \code{\link{indexShow}} will display on screen the full list of available indices and their codes.
#' The names of the internal functions calculating each index is also displayed, whose help files can aid in
#' the definition of index-specific arguments.
#'
#' @template templateParallel
#' @author M. Iturbide
#' @export
#' @examples
#' require(climate4R.datasets)
#' data("EOBS_Iberia_tas")
#' data("CFS_Iberia_tas")
#' fd <- indexGrid(
#'   tn = EOBS_Iberia_tas,
#'   time.resolution = "year",
#'   index.code = "FD"
#' )
#' per1 <- indexGrid(
#'   tn = EOBS_Iberia_tas,
#'   time.resolution = "year",
#'   index.code = "P",
#'   percent = 90
#' )
#' per2 <- indexGrid(
#'   tn = CFS_Iberia_tas,
#'   time.resolution = "year",
#'   index.code = "P",
#'   baseline = CFS_Iberia_tas,
#'   percent = 90
#' )
#' hdd <- indexGrid(
#'   tn = CFS_Iberia_tas,
#'   tx = CFS_Iberia_tas,
#'   tm = CFS_Iberia_tas,
#'   time.resolution = "year",
#'   index.code = "HDD"
#' )
indexGrid <- function(tn = NULL,
                      tx = NULL,
                      tm = NULL,
                      pr = NULL,
                      any = NULL,
                      baseline = NULL,
                      index.code,
                      time.resolution = "year",
                      ...,
                      parallel = FALSE,
                      max.ncores = 16,
                      ncores = NULL) {
  index.arg.list <- list(...)
  
  # Define a list of possible index codes
  choices <- c(
    "FD", "TNth", "TXth", "GDD", "CDD", "HDD", "P", "dt_st_rnagsn", "nm_flst_rnagsn",
    "dt_fnst_rnagsn", "dt_ed_rnagsn", "dl_agsn", "dc_agsn", "rn_agsn",
    "avrn_agsn", "dc_rnlg_agsn", "tm_agsn", "dc_txh_agsn", "dc_tnh_agsn",
    "gsl", "avg", "nd_thre", "nhw", "dr", "prcptot", "nrd", "lds", "sdii", "prcptot_thre", "ns",
    "TM_cold", "TM_warm", "FD_first", "FD_last", "FD_num", "Tth_first", "GST", "FD_prob", "Tth_prob", 
    "EHE", "GDD", "HI", "BEDD"
  )
  # Check if the index code is valid, otherwise throw an error message
  if (!index.code %in% choices) stop("Non valid index selected: Use indexShow() to select an index.")
  # If the index code is "FD", set the "th" argument to 0 and display a message
  if (index.code == "FD") {
    index.arg.list[["th"]] <- 0
    message("[", Sys.time(), "] th = 0 for index FD. Use index.code = 'TNth' to set a different threshold")
  }
  
  # Match the time resolution argument to one of three possible values
  time.resolution <- match.arg(time.resolution,
    choices = c("month", "year", "climatology")
  )
  
  # Check if the input temperature and precipitation data are daily data, otherwise throw an error message
  if (!is.null(tn)) {
    if (getTimeResolution(tn) != "DD") stop("Daily data is required as input", call. = FALSE)
  }
  if (!is.null(tx)) {
    if (getTimeResolution(tx) != "DD") stop("Daily data is required as input", call. = FALSE)
  }
  if (!is.null(tm)) {
    if (getTimeResolution(tm) != "DD") stop("Daily data is required as input", call. = FALSE)
  }
  if (!is.null(pr)) {
    if (getTimeResolution(pr) != "DD") stop("Daily data is required as input", call. = FALSE)
  }
  
  # Check if the baseline data is daily data and if the index code is "P", otherwise throw an error or warning message
  if (!is.null(baseline)) {
    if (!index.code %in% c("P")) {
      # If the index code is not "P", issue a warning and set baseline to NULL
      warning("Index.code is not 'P', baseline ignored")
      baseline <- NULL
    } else {
      # If the index code is "P", check if the baseline data is daily data, otherwise throw an error
      if (getTimeResolution(baseline) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
  }
  
  # Read metadata from the master file
  aux <- read.master()
  metadata <- aux[grep(paste0("^", index.code, "$"), aux$code, fixed = FALSE), ]
  
  # Check which input variables are present
  a <- c(!is.null(tn), !is.null(tx), !is.null(tm), !is.null(pr), !is.null(any)) %>% as.numeric()
  
  # If the index code is not "P", check if all required input variables are present, otherwise throw an error
  if (!index.code %in% c("P")) {
    b <- metadata[, 4:8] %>% as.numeric()
    if (any(b - a > 0)) {
      stop("The required input variable(s) for ", index.code,
        " index calculation are missing\nType \'?",
        metadata$indexfun, "\' for help",
        call. = FALSE
      )
    }
  } else {
	# If the index code is "P", check if only one input variable is present, otherwise throw an error
    b <- a
    if (sum(b) > 1) stop(index.code, " is applied to single variable.")
  }
  
  # Create a list with only the variables that were provided as input
  grid.list <- list("tn" = tn, "tx" = tx, "tm" = tm, "pr" = pr, "any" = any)[which(as.logical(b))]
  namesgridlist <- names(grid.list)
  
  # Operations for the consistency of the grids
  locs <- lapply(grid.list, isRegular)
  if (!any(sum(unlist(locs)) != 0, sum(unlist(locs)) != length(grid.list))) stop("Regular and Irregular grids can not be combined. See function interpGrid")
  
  # If there is more than one grid, intersect them and interpolate them to a common grid
  if (length(grid.list) > 1) {
    grid.list <- intersectGrid(grid.list, type = "temporal", which.return = 1:length(grid.list))
    names(grid.list) <- namesgridlist
    grid.list <- suppressMessages(lapply(grid.list, function(i) interpGrid(i, getGrid(grid.list[[1]]))))
  }
  
  # Remove any dimensions with a length of 1 from the grids
  grid.list <- lapply(grid.list, function(r) redim(r, drop = TRUE))
  grid.list <- lapply(grid.list, function(r) redim(r, loc = !unique(unlist(locs))))
  if (!unique(unlist(locs))) stop("The implementation of indexGrid for irregular grids is under development.")
  
  # Member loop preparation - Prepare the loop over the ensemble members
  ns.mem <- lapply(grid.list, function(r) getShape(r)[["member"]])
  if (sum(unlist(ns.mem) - rep(ns.mem[[1]], length(ns.mem))) != 0) stop("Number of members is different")
  n.mem <- unique(unlist(ns.mem))
  
  # If there is more than one member, check if parallel processing is possible and set up the function to apply to each member
  if (n.mem > 1) {
    parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
    apply_fun <- selectPar.pplyFun(parallel.pars, .pplyFUN = "lapply")
    if (parallel.pars$hasparallel) on.exit(parallel::stopCluster(parallel.pars$cl))
  } else {
    if (isTRUE(parallel)) message("NOTE: Parallel processing was skipped (unable to parallelize one single member)")
    apply_fun <- lapply
  }
  
  # Member loop
  # This loop is iterating through each member of the ensemble and calculating the corresponding index.
  message("[", Sys.time(), "] Calculating ", index.code, " ...")
  out.m <- apply_fun(1:n.mem, function(m) {
    if (sum(b) == 1 & is.null(baseline) & metadata$indexfun != "agroindexFAO" &
      metadata$indexfun != "agroindexFAO_tier1" & metadata$indexfun != "agroindex_agroclim") {
      # Indices from a single variable
      # This section of code is for indices that are calculated from a single variable.
      # The aggr.arg and fun.call variables determine how to aggregate the data based on the time resolution
      aggr.arg <- switch(time.resolution,
        "month" = "aggr.m",
        "year" = "aggr.y",
        "climatology" = "clim.fun"
      )
      fun.call <- switch(time.resolution,
        "month" = "aggregateGrid",
        "year" = "aggregateGrid",
        "climatology" = "climatology"
      )
      input.arg.list <- list()
      input.arg.list[["grid"]] <- subsetGrid(grid.list[[1]], members = m)
      input.arg.list[[aggr.arg]] <- c(list("FUN" = metadata$indexfun), index.arg.list)
      suppressMessages(do.call(fun.call, input.arg.list))
    } else {
      # Indices from multiple variables or for baseline methods
	  # Indices from multiple variables or for baseline methods
      # This section of code is for indices that are calculated from multiple variables or for baseline methods.
      # The grid.list.aux, months, and years variables are used to determine the data to use in the calculation.
      grid.list.aux <- lapply(grid.list, function(x) subsetGrid(x, members = m))
      months <- switch(time.resolution,
        "month" = as.list(getSeason(grid.list.aux[[1]])),
        "year" = list(getSeason(grid.list.aux[[1]])),
        "climatology" = list(getSeason(grid.list.aux[[1]]))
      )
      years <- switch(time.resolution,
        "month" = as.list(unique(getYearsAsINDEX(grid.list.aux[[1]]))),
        "year" = as.list(unique(getYearsAsINDEX(grid.list.aux[[1]]))),
        "climatology" = list(unique(getYearsAsINDEX(grid.list.aux[[1]])))
      )
	  
	  # Baseline loop
      if (!is.null(baseline)) {
	    # Subsetting the baseline grid by member
        baseline.sub <- suppressWarnings(subsetGrid(baseline, members = m))
		
		# Check if percentile or value was specified for index arguments
        if (is.null(index.arg.list[["percent"]]) & is.null(index.arg.list[["value"]])) stop("Baseline provided but percent or value not specified.")
        
		# If both percentile and value were specified, issue a warning and ignore value
		if (!is.null(index.arg.list[["percent"]]) & !is.null(index.arg.list[["value"]])) {
          warning("Values were given to both percentile and value... value will be ignored (set to NULL)")
          index.arg.list[["value"]] <- NULL
        }
		
		# If percentile was specified, calculate the value using climatology function with the percentile argument
        if (!is.null(index.arg.list[["percent"]])) {
          index.arg.list[["value"]] <- suppressMessages(
            redim(climatology(baseline.sub,
              clim.fun = list(FUN = percentile, percent = index.arg.list[["percent"]])
            ), drop = TRUE)[["Data"]]
          )
          index.arg.list[["percent"]] <- NULL
        } 
		# If value was specified, calculate the percentile using climatology function with the value argument
		else if (!is.null(index.arg.list[["value"]])) {
          index.arg.list[["percent"]] <- suppressMessages(
            redim(climatology(baseline.sub,
              clim.fun = list(FUN = percentile, value = index.arg.list[["value"]])
            ), drop = TRUE)[["Data"]]
          )
          index.arg.list[["value"]] <- NULL
        }
      }
	  
      # EXCEPTION for FAO INDICES and AGROCLIM INDICES (require lat, dates, and NO temporal subsetting)
      if (metadata$indexfun %in% c("agroindexFAO", "agroindexFAO_tier1", "agroindex_agroclim")) {
	    # If the index function is one of the FAO indices, aggregate the input grid to yearly resolution
        if (time.resolution != "year") message(index.code, " is calculated yaear by year by definition. argument time.resolution ignored.")
        out.aux <- suppressMessages(aggregateGrid(grid.list.aux[[1]], aggr.y = list(FUN = "mean", na.rm = TRUE)))
        
		# Get the input arguments for the index function and the dates
		input.arg.list <- lapply(grid.list.aux, function(d) d[["Data"]])
        datess <- as.Date(grid.list.aux[[1]][["Dates"]][["start"]])
        datess <- cbind(as.numeric(format(datess, "%Y")), as.numeric(format(datess, "%m")), as.numeric(format(datess, "%d")))
        lats <- grid.list.aux[[1]][["xyCoords"]][["y"]]
        
		# Add the necessary input arguments to the index argument list
		index.arg.list[["dates"]] <- datess
        index.arg.list[["index.code"]] <- index.code
        
		# Loop over latitudes and longitudes, calling the index function for each point
		latloop <- lapply(1:length(lats), function(l) {
          lonloop <- lapply(1:getShape(grid.list.aux[[1]])["lon"], function(lo) {
            index.arg.list[["lat"]] <- lats[l]
            do.call(metadata$indexfun, c(lapply(input.arg.list, function(z) z[, l, lo]), index.arg.list))
          })
          do.call("abind", list(lonloop, along = 0))
        })
		
		# Merge the results back into a grid
        out.aux[["Data"]] <- unname(aperm(do.call("abind", list(latloop, along = 0)), c(3, 1, 2)))
        attr(out.aux[["Data"]], "dimensions") <- c("time", "lat", "lon")
        out.aux
      } else {
        yg <- lapply(years, function(yi) {
          mg <- lapply(months, function(mi) {
		    # Compute the climatology for the first grid in the list and suppress any messages.
            out.aux <- suppressMessages(climatology(grid.list.aux[[1]]))
			
			# Subset the grid list for the current year and month
            grid.list.sub <- lapply(grid.list.aux, function(x) subsetGrid(x, years = yi, season = mi))
            
			# Extract the data component of each grid and store in a list
			input.arg.list <- lapply(grid.list.sub, function(d) d[["Data"]])
            
			# If the index code is "P", add a name to the input argument list
			if (index.code == "P") names(input.arg.list) <- "var"
			
			# Combine the data list and the index argument list
            input.arg.list <- c(input.arg.list, index.arg.list)
			
			# Compute the index using the specified function and input argument list
            out.aux[["Data"]] <- unname(do.call(metadata$indexfun, input.arg.list))
			
			# Set the dimensions attribute of the index data to "lat" and "lon"
            attr(out.aux[["Data"]], "dimensions") <- c("lat", "lon")
			
			# Return the output for the current month
            out.aux
          })
		  # Combine the output for all months into a single grid along the time dimension
          tryCatch(
            {
              bindGrid(mg, dimension = "time")
            },
            error = function(err) {
			  # If an error occurs, return the output for each month separately
              unlist(mg, recursive = FALSE)
            }
          )
        })
		# Combine the output for all years into a single grid along the time dimension
        tryCatch(
          {
            bindGrid(yg, dimension = "time")
          },
          error = function(err) {
            unlist(yg, recursive = FALSE)
          }
        )
      }
    }
  })
  # Bind the output data along the member dimension
  out <- suppressMessages(suppressWarnings(bindGrid(out.m, dimension = "member")))
  # Add the variable name and level information to the output data's Variable field
  out[["Variable"]] <- list(
    "varName" = index.code,
    "level" = out[["Variable"]][["level"]]
  )
  # Add metadata to the Variable field
  attr(out[["Variable"]], "description") <- metadata$description
  attr(out[["Variable"]], "units") <- metadata$units
  attr(out[["Variable"]], "longname") <- metadata$longname
  message("[", Sys.time(), "] Done")
  
  # Return the processed output data
  return(out)
}


#' @title List all available Indices
#' @description Print a table with a summary of the available indices
#' @return Print a table on the screen with the following columns:
#' \itemize{
#' \item \strong{code}: Code of the index. This is the character string used as input value
#' for the argument \code{index.code} in \code{\link{indexGrid}}
#' \item \strong{longname}: Long description of the index
#' \item \strong{index.fun}: The name of the internal function used to calculate it
#' \item \strong{tn, tx, tm, pr}: A logical value (0/1) indicating the input variables required for index calculation
#' \item \strong{units}: The units of the index (when different from those of the input variable)
#' }
#' @author J. Bedia, M. Iturbide
#' @export

indexShow <- function() {
  read.master()
}



#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom utils read.table

read.master <- function() {
  system.file("master", package = "climate4R.indices") %>% read.table(
    header = TRUE,
    sep = ";",
    stringsAsFactors = FALSE,
    na.strings = ""
  )
}
