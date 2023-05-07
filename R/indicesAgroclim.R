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
  choices <- c(
    "TM_cold", "TM_warm", "FD_first", "FD_last", "FD_num",
    "Tth_first", "GST", "FD_prob", "Tth_prob", "EHE",
    "GDD", "HI", "BEDD"
  )
  if (!index.code %in% choices) stop("Non valid index selected: Use indexShow() to select an index.")

  # Remove "lat" argument when index.code is not in c("GDD", "GST", "HI", "BEDD")
  input.arg.list <- list(...)
  if (!index.code %in% c("GDD", "GST", "HI", "BEDD")) input.arg.list$lat <- null

  # Apply the atomic function
  do.call(index.code, input.arg.list)
}

#########################
## auxiliary functions ##
#########################

##################
## year_StartEnd ##
##################
#' @title Function to find the position marking the start and the end of a given year (or a user-defined portion of the year)
#' @return Indices marking the start and the end of a given year (or a user-defined portion of the year). The search is done within the "dates" matrix
#' @param dates Matrix containing the full range of dates (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param year Year of interest (e.g. 1995)
#' @param year.start (Optional) User-defined start of the year [in "MM-DD" format].
#' @param year.end (Optional) User-defined end of the year [in "MM-DD" format].
#' @author JJ. Velasco
#' @export

year_StartEnd <- function(dates, year, year.start = NULL, year.end = NULL) {
  if (!is.null(year.start) & !is.null(year.end)) {
    # Add year to year.start and year.end
    year.start.c <- paste(as.character(year), year.start, sep = "-")
    year.end.c <- paste(as.character(year), year.start, sep = "-")

    ind.start <- which(dates[, 1] == as.numeric(substr(year.start.c, 1, 4)) & # start of the year (as defined by the user)
      dates[, 2] == as.numeric(substr(year.start.c, 6, 7)) &
      dates[, 3] == as.numeric(substr(year.start.c, 9, 10)))
    ind.end <- which(dates[, 1] == as.numeric(substr(year.end.c, 1, 4)) & # end of the year (as defined by the user)
      dates[, 2] == as.numeric(substr(year.end.c, 6, 7)) &
      dates[, 3] == as.numeric(substr(year.end.c, 9, 10)))
  } else {
    ind.start <- which(dates[, 1] == year & dates[, 2] == 1 & dates[, 3] == 1) # year's definition by default (1-Jan to 31-Dec)
    ind.end <- which(dates[, 1] == year & dates[, 2] == 12 & dates[, 3] == 31) # year's definition by default (1-Jan to 31-Dec)
  }
  out <- list()
  out$start <- ind.start
  out$end <- ind.end
  return(out)
}


################################################################################
## agroclim indices (https://doi.org/10.1016/j.compag.2020.105473)            ##
################################################################################

#############
## TM_cold ##
#############
#' @title Coldest month of the year.
#' @description Calculates the mean temperature of the coldest month of the year.
#' @return Depending on argument type, the output will be a numeric value with the month (type_output = "month") or with temperature (type_output = "temp") per year.
#' @param tn Vector with daily minimum temperature.
#' @param dates Matrix containing the full range of dates corresponding to "tn" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param type_output Type of output. It can be "temp" for temperature or "month" for the number of the coldest month. The default value is "temp".
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995).
#' @param year.start (Optional) Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
#' @param year.end (Optional) Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored. The default value is 25.
#' @details Depending on argument type_output, the output will be a numeric value with the month (type_output = "month") or with temperature (type_output = "temp") per year.
#' @author JJ. Velasco
#' @examples
#'
#' index <- TM_cold(tn, dates, "temp") # call to the function
#' index <- TM_cold(tn, dates, "month", # call to the function
#'   year = 1994:2018,
#'   year.start = "03-21",
#'   year.end = "10-15"
#' )
#' @export

TM_cold <- function(tn, dates, type_output = "temp", year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # if (!is.null(lat)) warning("This index doesn't use latitude information.")

  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tn.year <- tn[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tn.year)) < 0.01 * pnan * length(tn.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::coldMonth(
            mn = tn.year,
            dates = dates.date.year,
            type = type_output
          )
        }
      }
    }
  }
  return(index)
}
# index <- TM_cold(tn, dates, "temp") # call to the function
# index <- TM_cold(tn, dates, "month", # call to the function
#   year = 1994:2018,
#   year.start = "03-21",
#   year.end = "10-15"
# )



#############
## TM_warm ##
#############
#' @title Warmest month of the year.
#' @description Calculates the mean temperature of the warmest month of the year.
#' @return Depending on argument type, the output will be a numeric value with the month (type_output = "month") or with temperature (type_output = "temp") per year.
#' @param tx Vector with daily maximum temperature.
#' @param dates Matrix containing the full range of dates corresponding to "tx" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param type_output Type of output. It can be "temp" for temperature or "month" for the number of the warmest month. The default value is "temp".
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995).
#' @param year.start (Optional) Vector of dates [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
#' @param year.end (Optional) Vector of dates [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored. The default value is 25.
#' @details Depending on argument type_output, the output will be a numeric value with the month (type_output = "month") or with temperature (type_output = "temp") per year.
#' @author JJ. Velasco
#' @examples
#'
#' index <- TM_warm(tx, dates, "temp") # call to the function
#' index <- TM_warm(tx, dates, "month", # call to the function
#'   year = 1994:2018,
#'   year.start = "03-21",
#'   year.end = "10-15"
#' )
#' @export

TM_warm <- function(tx, dates, type_output = "temp", year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # if (!is.null(lat)) warning("This index doesn't use latitude information.")

  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tx.year <- tx[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tx.year)) < 0.01 * pnan * length(tx.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::warmMonth(
            mx = tx.year,
            dates = dates.date.year,
            type = type_output
          )
        }
      }
    }
  }
  return(index)
}
# index <- TM_warm(tx, dates, "temp") # call to the function
# index <- TM_warm(tx, dates, "month", # call to the function
#   year = 1994:2018,
#   year.start = "03-21",
#   year.end = "10-15"
# )



##############
## FD_first ##
##############
#' @title Mean first frost day.
#' @description Calculates the first frost day of each year.
#' @return Depending on argument type_output, the output will be a numeric vector of julian days (type_output = "doy") or a vector of characters with dates (type_output = "date").
#' @param tn Vector with daily minimum temperature.
#' @param dates Matrix containing the full range of dates corresponding to "tn" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param type_output Type of output. It can be "doy" for day of the year (julian day) or "date" for data format ("dd-mm"). The default value is "doy".
#' @param threshold Temperature threshold considered to trigger frost occurrence. The default value is 0 (ºC).
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995).
#' @param year.start (Optional) Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
#' @param year.end (Optional) Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored. The default value is 25.
#' @details Depending on argument type_output, the output will be a numeric value with the month (type_output = "month") or with temperature (type_output = "temp") per year.
#' @author JJ. Velasco
#' @examples
#'
#' index <- FD_first(tx, dates, "doy", 33) # call to the function
#' index <- FD_first(tx, dates, "date", 33, # call to the function
#'   year = 1994:2018,
#'   year.start = "03-21",
#'   year.end = "10-15"
#' )
#' @export

FD_first <- function(tn, dates, type_output = "doy", threshold = 0, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # if (!is.null(lat)) warning("This index doesn't use latitude information.")

  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tn.year <- tn[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tn.year)) < 0.01 * pnan * length(tn.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::firstFrost(
            mn = tn.year,
            dates = dates.date.year,
            iniday = format(dates.date.year[1], "%m-%d"),
            endday = format(dates.date.year[length(dates.date.year)], "%m-%d"),
            type = type_output,
            thres = threshold
          )
        }
      }
    }
  }
  return(index)
}
# index <- FD_first(tx, dates, "doy", 33) # call to the function
# index <- FD_first(tx, dates, "date", 33, # call to the function
#   year = 1994:2018,
#   year.start = "03-21",
#   year.end = "10-15"
# )


##############
## FD_last ##
##############
#' @title Mean last frost day.
#' @description Calculates the last frost day of each year.
#' @return Depending on argument type_output, the output will be a numeric vector of julian days (type_output = "doy") or a vector of characters with dates (type_output = "date").
#' @param tn Vector of daily minimum temperature.
#' @param dates Matrix containing the full range of dates corresponding to "tn" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param type_output Type of output. It can be "doy" for day of the year (julian day) or "date" for data format ("dd-mm"). The default value is "doy".
#' @param threshold Temperature threshold considered to trigger frost occurrence. The default value is 0 (ºC).
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995).
#' @param year.start (Optional) Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
#' @param year.end (Optional) Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored. The default value is 25.
#' @details Depending on argument type_output, the output will be a numeric value with the month (type_output = "month") or with temperature (type_output = "temp") per year.
#' @author JJ. Velasco
#' @examples
#'
#' index <- FD_last(tn, dates, "doy", 2) # call to the function
#' index <- FD_last(tn, dates, "date", 0, # call to the function
#'   year = 1994:2018,
#'   year.start = "03-21",
#'   year.end = "10-15"
#' )
#' @export

FD_last <- function(tn, dates, type_output = "doy", threshold = 0, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # if (!is.null(lat)) warning("This index doesn't use latitude information.")

  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tn.year <- tn[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tn.year)) < 0.01 * pnan * length(tn.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::lastFrost(
            mn = tn.year,
            dates = dates.date.year,
            iniday = format(dates.date.year[1], "%m-%d"),
            endday = format(dates.date.year[length(dates.date.year)], "%m-%d"),
            type = type_output,
            thres = threshold
          )
        }
      }
    }
  }
  return(index)
}
# index <- FD_last(tx, dates, "doy", 2) # call to the function
# index <- FD_last(tx, dates, "date", 0, # call to the function
#   year = 1994:2018,
#   year.start = "03-21",
#   year.end = "10-15"
# )


##############
## FD_num ##
##############
#' @title Number of frost days.
#' @description Calculates the number of frost days of each year.
#' @return A numeric vector with the annual number of frost days is returned.
#' @param tn Vector with daily minimum temperature.
#' @param dates Matrix containing the full range of dates corresponding to "tn" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param threshold Temperature threshold considered to trigger frost occurrence. The default value is 0 (ºC).
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995).
#' @param year.start (Optional) Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
#' @param year.end (Optional) Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored. The default value is 25.
#' @details Despite the logical threshold of temperature is 0 ºC to determine frost occurrence, the argument "threshold" is open to change in case of different units of temperature.
#' @author JJ. Velasco
#' @examples
#'
#' index <- FD_num(tn, dates, 2) # call to the function
#' index <- FD_num(tn, dates, 0, # call to the function
#'   year = 1994:2018,
#'   year.start = "03-21",
#'   year.end = "10-15"
#' )
#' @export

FD_num <- function(tn, dates, threshold = 0, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # if (!is.null(lat)) warning("This index doesn't use latitude information.")

  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tn.year <- tn[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tn.year)) < 0.01 * pnan * length(tn.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::lastFrost(
            mn = tn.year,
            dates = dates.date.year,
            iniday = format(dates.date.year[1], "%m-%d"),
            endday = format(dates.date.year[length(dates.date.year)], "%m-%d"),
            thres = threshold
          )
        }
      }
    }
  }
  return(index)
}
# index <- FD_num(tn, dates, 2) # call to the function
# index <- FD_num(tn, dates, 0, # call to the function
#   year = 1994:2018,
#   year.start = "03-21",
#   year.end = "10-15"
# )



###############
## Tth_first ##
###############
#' @title First day in the year where P(tmax>threshold) >= threshold_prob
#' @description Calculates the first day in the year where the probability of temperature over a threshold is higher than a predefined threshold.
#' @return Depending on argument type_output, the output will be a numeric vector of julian days (type_output = "doy") or a vector of characters with dates (type_output = "date").
#' @param tx Vector with daily maximum temperature.
#' @param dates Matrix containing the full range of dates corresponding to "tn" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param type_output Type of output. It can be "doy" for day of the year (julian day) or "date" for data format ("dd-mm"). The default value is "doy".
#' @param threshold Temperature threshold considered to trigger occurrence. The default value is 35 (ºC).
#' @param threshold_prob Probability threshold indicating the minimum probability. The default value is 0.1
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995).
#' @param year.start (Optional) Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
#' @param year.end (Optional) Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored. The default value is 25.
#' @details Depending on argument type_output, the output will be a numeric vector of julian days (type_output = "doy") or a vector of characters with dates (type_output = "date").
#' @author JJ. Velasco
#' @examples
#'
#' index <- Tth_first(tx, dates, "doy", 32, 0.25) # call to the function
#' index <- Tth_first(tx, dates, "date", 35, 0.2, # call to the function
#'   year = 1994:2018,
#'   year.start = "03-21",
#'   year.end = "10-15",
#' )
#' @export

Tth_first <- function(tx, dates, type_output = "doy", threshold = 35, threshold_prob = 0.1, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # if (!is.null(lat)) warning("This index doesn't use latitude information.")

  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tx.year <- tx[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tx.year)) < 0.01 * pnan * length(tx.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::tempDayprob(
            mx = tx.year,
            dates = dates.date.year,
            iniday = format(dates.date.year[1], "%m-%d"),
            endday = format(dates.date.year[length(dates.date.year)], "%m-%d"),
            type = type_output,
            thres = threshold,
            prob = threshold_prob
          )
        }
      }
    }
  }
  return(index)
}
# index <- Tth_first(tx, dates, "doy", 32, 0.25) # call to the function
# index <- Tth_first(tx, dates, "date", 35, 0.2, # call to the function
#   year = 1994:2018,
#   year.start = "03-21",
#   year.end = "10-15",
# )


#########
## GST ##
#########
#' @title Growing Season Temperature (GST)
#' @description Growing Season Temperature (GST). Mean daily temperature in growing season.
#' @return A numeric vector with annual values is returned.
#' @param tn Vector with daily minimum temperature
#' @param tx Vector with daily maximum temperature
#' @param dates Matrix containing the full range of dates corresponding to "tx" and "tn" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param lat Numeric value indicating the latitude of location.
#' @param year Vector with years of interest (e.g. 1990:1995)
# #' @param year.start Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
# #' @param year.end Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @details Depending on the latitude, the function detects the hemisphere and considers growing season from 1st April to 31st October (northern hemisphere) or from 1st October to 30rd April (southern hemisphere).
#' @references Jones G, Duff A, Hall A, Myers J (2010) Spatial Analysis of Climate in Winegrape Growing Regions in the Western United States. Am. J. Enol. Vitic. 61:3.
#' @author JJ. Velasco
#' @examples
#'
#' index <- GST(tn, tx, dates, 26) # call to the function
#' index <- GST(tn, tx, dates, 26, # call to the function
#'   year = 1994:2018
#' )
#' @export

GST <- function(tn, tx, dates, lat, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tx.year <- tx[ind.year$start:ind.year$end]
        tn.year <- tn[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tx.year)) < 0.01 * pnan * length(tx.year) & sum(is.na(tn.year)) < 0.01 * pnan * length(tn.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::gst(
            mx = tx.year,
            mn = tn.year,
            dates = dates.date.year,
            lati = lat
          )
        }
      }
    }
  }
  return(index)
}
# index = GST(tn, tx, dates, 26)  # call to the function
# index = GST(tn, tx, dates, 26,  # call to the function
#             year = 1994:2018
# )



###############
## FD_prob ##
###############
#' @title First day in the year where P(tmin<threshold) <= threshold_prob
#' @description Calculates the first day in the year where the probability of temperature below a threshold is below than a predefined threshold.
#' @return Depending on argument type_output, the output will be a numeric vector of julian days (type_output = "doy") or a vector of characters with dates (type_output = "date").
#' @param tn Vector with daily minimum temperature.
#' @param dates Matrix containing the full range of dates corresponding to "tn" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param type_output Type of output. It can be "doy" for day of the year (julian day) or "date" for data format ("dd-mm"). The default value is "doy".
#' @param threshold Temperature threshold considered to trigger frost occurrence. The default value is 0 (ºC).
#' @param threshold_prob Probability threshold indicating the maximum probability. The default value is 0.1
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995).
#' @param year.start (Optional) Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
#' @param year.end (Optional) Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored. The default value is 25.
#' @details Depending on argument type_output, the output will be a numeric vector of julian days (type_output = "doy") or a vector of characters with dates (type_output = "date").
#' @author JJ. Velasco
#' @examples
#'
#' index <- FD_prob(tn, dates, "doy", 0, 0.25) # call to the function
#' index <- FD_prob(tn, dates, "date", -2, 0.2, # call to the function
#'   year = 1994:2018,
#'   year.start = "03-21",
#'   year.end = "10-15",
#' )
#' @export

FD_prob <- function(tn, dates, type_output = "doy", threshold = 0, threshold_prob = 0.1, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # if (!is.null(lat)) warning("This index doesn't use latitude information.")

  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tn.year <- tn[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tn.year)) < 0.01 * pnan * length(tn.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::frostProb(
            mn = tn.year,
            dates = dates.date.year,
            iniday = format(dates.date.year[1], "%m-%d"),
            endday = format(dates.date.year[length(dates.date.year)], "%m-%d"),
            type = type_output,
            thres = threshold,
            prob = threshold_prob
          )
        }
      }
    }
  }
  return(index)
}
# index <- FD_prob(tn, dates, "doy", 0, 0.25) # call to the function
# index <- FD_prob(tn, dates, "date", -2, 0.2, # call to the function
#   year = 1994:2018,
#   year.start = "03-21",
#   year.end = "10-15",
# )


##############
## Tth_prob ##
##############
#' @title Probability of exceed a predefined temperature value.
#' @description Calculates the probability of exceed a predefined temperature value.
#' @return A numeric vector with annual values is returned.
#' @param tx Vector with daily (usually maximum) temperature.
#' @param dates Matrix containing the full range of dates corresponding to "tx" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param threshold Temperature threshold considered to trigger occurrence. The default value is 20 (ºC).
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995).
#' @param year.start (Optional) Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
#' @param year.end (Optional) Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored. The default value is 25.
#' @details Despite the logical threshold of temperature is 20 ºC to determine frost occurrence, the argument "threshold" is open to change in case of different units of temperature.
#' @references Trnka M, Rotter RP, Ruiz-Ramos M, Kersebaum KC, Olesen JE, Zalud Z, Semenov MA (2014) Adverse weather conditions for European wheat production will become more frequent with climate change. Nature Climate Change volume 4, pages 637–643.
#' @author JJ. Velasco
#' @examples
#' index <- Tth_prob(tx, dates, 20) # call to the function
#' index <- Tth_prob(tx, dates, 25, # call to the function
#'   year = 1994:2018,
#'   year.start = "03-21",
#'   year.end = "10-15"
#' )
#' @export

Tth_prob <- function(tx, dates, threshold = 20, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # if (!is.null(lat)) warning("This index doesn't use latitude information.")

  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tx.year <- tx[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tx.year)) < 0.01 * pnan * length(tx.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::tempProb(
            mx = tx.year,
            dates = dates.date.year,
            thres = threshold,
            month = NULL
          )
        }
      }
    }
  }
  return(index)
}
# index <- Tth_prob(tx, dates, 20) # call to the function
# index <- Tth_prob(tx, dates, 25, # call to the function
#   year = 1994:2018,
#   year.start = "03-21",
#   year.end = "10-15"
# )



#########
## EHE ##
#########
#' @title Extreme Heat Exposure (EHE)
#' @description Extreme Heat Exposure (EHE). Useful for climatic risks assessement on wheat and barley.
#' @return If op = "first", the function returns the first day (date format) when the first event is triggered. If op =='doy', julian day is returned. If op = "number", the funciton returns the number of events occurred in the year.
#' @param tx Vector with daily (usually maximum) temperature.
#' @param dates Matrix containing the full range of dates corresponding to "tx" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param op Type of output. If op = "first", the function returns the first day (date format) when the first event is triggered. If op =='doy', julian day is returned. If op = "number", the funciton returns the number of events occurred in the year. The default value is "first".
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995).
#' @param year.start (Optional) Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
#' @param year.end (Optional) Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored. The default value is 25.
#' @details Adapted from Trnka et al. (2014). Event is triggered when the Tmax is above +35°C for at least three days during the period from five days after anthesis (supposed to be May-1st) to maturity (suposed to be July-31st). The minimum daily temperature is usually measured 2 m above ground; thus, the actual crop temperature might be even lower.
#' @author JJ. Velasco
#' @examples
#'
#' index <- EHE(tx, dates, op = "first") # call to the function
#' index <- EHE(tx, dates,
#'   op = "doy", # call to the function
#'   year = 1994:2018,
#'   year.start = "03-21",
#'   year.end = "10-15"
#' )
#' @export

EHE <- function(tx, dates, op = "first", year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  # if (!is.null(lat)) warning("This index doesn't use latitude information.")

  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tx.year <- tx[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tx.year)) < 0.01 * pnan * length(tx.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::ehe(
            mx = tx.year,
            dates = dates.date.year,
            op = op
          )
        }
      }
    }
  }
  return(index)
}
# index <- EHE(tx, dates, op = "first") # call to the function
# index <- EHE(tx, dates, op = "doy", # call to the function
#   year = 1994:2018,
#   year.start = "03-21",
#   year.end = "10-15"
# )


#########
## GDD ##
#########
#' @title Growing Degree Days (GDD)
#' @description Growing Degree Day (GDD) or Winkler index. Useful as a zoning tool to differentiate between grape varieties and climate (Winkler et al. 1974).
#' @return Number of Huglin index (per year)
#' @param tn Vector with daily minimum temperature
#' @param tx Vector with daily maximum temperature
#' @param dates Matrix containing the full range of dates corresponding to "tx" and "tn" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param lat Numeric value indicating the latitude of location.
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995)
# #' @param year.start (Optional) Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
# #' @param year.end (Optional) Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @details Depending on the latitude, the function detects the hemisphere and considers growing season from 1st April to 30rd September (northern hemisphere) or from 1st October to 31st March (southern hemisphere).
#' @references Huglin P. (1978) Nouveau mode d'evaluation des possibilities heliothermiques d'un milieu viticole. CR Acad Agr 64: 1117–1126
#' @author JJ. Velasco
#' @examples
#'
#' index <- GDD(tn, tx, dates, lat = 26) # call to the function
#' index <- GDD(tn, tx, dates,
#'   lat = 26, # call to the function
#'   year = 1994:2018
#' )
#' @export

GDD <- function(tn, tx, dates, lat, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tx.year <- tx[ind.year$start:ind.year$end]
        tn.year <- tn[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tx.year)) < 0.01 * pnan * length(tx.year) & sum(is.na(tn.year)) < 0.01 * pnan * length(tn.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::gdd(
            mx = tx.year,
            mn = tn.year,
            dates = dates.date.year,
            lati = lat,
            iniday = format(dates.date.year[1], "%m-%d"),
            endday = format(dates.date.year[length(dates.date.year)], "%m-%d"),
          )
        }
      }
    }
  }
  return(index)
}
# index <- GDD(tn, tx, dates, 26) # call to the function
# index <- GDD(tn, tx, dates, 26, # call to the function
#   year = 1994:2018
# )



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
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995)
# #' @param year.start (Optional) Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
# #' @param year.end (Optional) Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @details Depending on the latitude, the function detects the hemisphere and considers growing season from 1st April to 30rd September (northern hemisphere) or from 1st October to 31st March (southern hemisphere).
#' @references Huglin P. (1978) Nouveau mode d'evaluation des possibilities heliothermiques d'un milieu viticole. CR Acad Agr 64: 1117–1126
#' @author JJ. Velasco
#' @examples
#'
#' index <- HI(tn, tx, dates, lat = 26) # call to the function
#' index <- HI(tn, tx, dates,
#'   lat = 26, # call to the function
#'   year = 1994:2018
#' )
#' @export

HI <- function(tn, tx, dates, lat, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tx.year <- tx[ind.year$start:ind.year$end]
        tn.year <- tn[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tx.year)) < 0.01 * pnan * length(tx.year) & sum(is.na(tn.year)) < 0.01 * pnan * length(tn.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::hi(
            mx = tx.year,
            mn = tn.year,
            dates = dates.date.year,
            lati = lat
          )
        }
      }
    }
  }
  return(index)
}
# index = HI(tn, tx, dates, 26)  # call to the function
# index = HI(tn, tx, dates, 26,  # call to the function
#            year = 1994:2018
# )


##########
## BEDD ##
##########
#' @title Biologically effective degree-days (BEDD)
#' @description Biologically effective degree-days (BEDD). Gladstones, J. (1992).
#' @return The sum of degree-days (BEDD) is returned as a numeric value.
#' @param tn Vector with daily minimum temperature
#' @param tx Vector with daily maximum temperature
#' @param dates Matrix containing the full range of dates corresponding to "tx" and "tn" (ndates x 3 size); e.g. rbind(c(1995, 3, 1), c(1995, 3, 2), ...)
#' @param lat Numeric value indicating the latitude of location.
#' @param year (Optional) Vector with years of interest (e.g. 1990:1995)
# #' @param year.start (Optional) Day [in "MM-DD" format] defining the beginning of a portion of interest within each year (e.g., the agronomic season).
# #' @param year.end (Optional) Day [in "MM-DD" format] defining the end of a portion of interest within each year (e.g., the agronomic season).
#' @param pnan Any year with a percentage of NA data above "pnan" will be ignored
#' @details Depending on the latitude, the function detects the hemisphere and considers growing season from 1st April to 30rd September (northern hemisphere) or from 1st October to 31st March (southern hemisphere).
#' @references Gladstones, J. (1992) Viticulture and environment (Winetitles: Adelaide).
#' @author JJ. Velasco
#' @examples
#'
#' index <- BEDD(tn, tx, dates, lat = 26) # call to the function
#' index <- BEDD(tn, tx, dates,
#'   lat = 26, # call to the function
#'   year = 1994:2018
#' )
#' @export

BEDD <- function(tn, tx, dates, lat, year = NULL, year.start = NULL, year.end = NULL, pnan = 25) {
  if (is.null(year)) {
    year <- unique(dates[, 1]) # years of analysis
  }

  # initializing output
  index <- rep(NA, 1, length(year))

  for (iyear in year) {
    if (!is.null(year.start) & !is.null(year.end)) {
      ind.year <- year_StartEnd(dates, iyear, year.start = year.start, year.end = year.end) # bounding dates defining the portion of year of interest
    } else {
      ind.year <- year_StartEnd(dates, iyear, year.start = NULL, year.end = NULL) # bounding dates defining the year of interest
    }

    if (length(ind.year$start) != 0 & length(ind.year$end) != 0) {
      if (!is.na(ind.year$start) & !is.na(ind.year$end)) {

        # Select the corresponding data range
        tx.year <- tx[ind.year$start:ind.year$end]
        tn.year <- tn[ind.year$start:ind.year$end]

        # Select the corresponding date range
        dates.matrix.year <- dates[ind.year$start:ind.year$end, ]
        dates.date.year <- as.Date(apply(dates.matrix.year, 1, paste, collapse = "-"),
          format = "%Y-%m-%d"
        )

        # asking for a minimum of pnan (%) of non-missing days
        if (sum(is.na(tx.year)) < 0.01 * pnan * length(tx.year) & sum(is.na(tn.year)) < 0.01 * pnan * length(tn.year)) {

          # Calculate the agroclim index
          library(agroclim)
          index[year == iyear] <- agroclim::bedd(
            mx = tx.year,
            mn = tn.year,
            dates = dates.date.year,
            lati = lat
          )
        }
      }
    }
  }
  return(index)
}
# index = BEDD(tn, tx, dates, 26)  # call to the function
# index = BEDD(tn, tx, dates, 26,  # call to the function
#            year = 1994:2018
# )
