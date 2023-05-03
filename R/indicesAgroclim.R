#' @title Wrapper to call agroclim index calculation functions.
#' @param index.code To call to the atomic function of the same name
#' @param ... Other parameters that can be passed to the selected index function (see Details).
#' @details Index calculation is done by different functions that receive the same name 
#' as the corresponding index code. Therefore, type \code{?index.code} to know the arguments
#' that are used in each case, e.g. \code{?HI}.
#' @section Index codes:
#'    \strong{TM_cold} 
#'    
#'    \strong{TM_warm}  
#'    
#'    \strong{FD_first}
#'  
#'    \strong{FD_last}
#'    
#'    \strong{FD_num}
#'    
#'    \strong{Tth_first}
#'    
#'    \strong{GST}
#'    
#'    \strong{FD_prob}
#'    
#'    \strong{Tth_prob}
#'    
#'    \strong{EHE} 
#'    
#'    \strong{GDD} 
#'    
#'    \strong{HI}
#'    
#'    \strong{BEDD}
#'    
#' @author JJ. Velasco
#' @export

agroindex_agroclim <- function(index.code, ...) {
  choices <- c("TM_cold", "TM_warm", "FD_first", "FD_last", "FD_num", 
               "Tth_first", "GST", "FD_prob", "Tth_prob", "EHE",
               "GDD", "HI", "BEDD")
  if (!index.code %in% choices) stop("Non valid index selected: Use indexShow() to select an index.")
  
  # Remove "lat" argument
  input.arg.list <- list(...)
  if (!index.code %in% c("GDD", "GST", "HI", "BEDD"))  input.arg.list$lat <- null
  
  # Apply the atomic function
  do.call(index.code, input.arg.list)
}


################################################################################
## agroclim indices (https://doi.org/10.1016/j.compag.2020.105473)            ##
################################################################################

#########
## HI ##
#########
#' @title Huglin Heliothermal Index (HI)
#' @description Huglin Heliothermal Index (HI). Useful as a zoning tool (Huglin 1978).
#' @return Number of Huglin index (per year)
#' @param tn Vector with daily minimum temperature
#' @param tx Vector with daily maximum temperature
#' @param dates Matrix containing the full range of dates corresponding to "tx" and "tn" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param lat Numeric value indicating the latitude of location.
#' @param year Vector with years of interest (e.g. 1990:1995)
#' @param year.start Vector of dates [in "YYYY-MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season)
#' @param year.end Vector of dates [in "YYYY-MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season)
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @details Depending on the latitude, the function detects the hemisphere and considers growing season from 1st April to 30rd September (northern hemisphere) or from 1st October to 31st March (southern hemisphere).
#' @references Huglin P. (1978) Nouveau mode d'evaluation des possibilities heliothermiques d'un milieu viticole. CR Acad Agr 64: 1117–1126
#' @author JJ. Velasco
#' @examples
#'
#'   index = HI(tm, tn, dates, 26)  # call to the function
#'   index = HI(tm, tn, dates, 26,  # call to the function
#'              year = 1994:2018, 
#'              year.start = AS.dates[[1]]$sdate$days[, 15], 
#'              year.end = AS.dates[[1]]$edate$days[, 15])  
#'
#' @export

HI <- function(tn, tx, dates, lat, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  if(!is.null(lat)) warning("This index doesn't use latitude information.")
  
  if (is.null(year)) {
    year = unique(dates[, 1])  # years of analysis
  }
  
  # initializing output
  index = rep(NA, 1, length(year))  
  
  for (iyear in year) {
    
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year = yearStartEnd(dates, iyear, year.start = year.start[year == iyear], year.end = year.end[year == iyear])  # bounding dates defining the portion of year of interest
    } else {
      ind.year = yearStartEnd(dates, iyear, year.start = NULL, year.end = NULL)  # bounding dates defining the year of interest
    }
    
    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {
        
        # Seleccionamos el intervalo de los datos correspondiente
        tx.year = tx[ind.year$start:ind.year$end]
        tn.year = tn[ind.year$start:ind.year$end]
        
        # Seleccionamos el intervalo de fechas correspondiente
        dates.matrix.year = dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"), 
                                   format = "%Y-%m-%d")
        
        if (sum(is.na(tx.year)) < 0.01*pnan*length(tx.year)) {  # asking for a minimum of pnan (%) of non-missing days 
          
          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] = agroclim::hi(mx = tx.year, 
                                           mn = tn.year, 
                                           dates = dates.date.year, 
                                           lati = lat)
        }
      }
    }
  }
  return(index)
}
# index = HI(tm, tn, dates, 26)  # call to the function
# index = HI(tm, tn, dates, 26,  # call to the function
#            year = 1994:2018, 
#            year.start = AS.dates[[1]]$sdate$days[, 15], 
#            year.end = AS.dates[[1]]$edate$days[, 15])  

