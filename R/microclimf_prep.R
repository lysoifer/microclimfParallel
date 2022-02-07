#' seasonalRast_to_hourlyArray
#'
#' Takes a spatRaster with each layer representing monthly seasonality of a different variable, replicates monthly seasonality to produce hourly, and converts to array.
#'
#' @param seasonalRast raster with nlyrs = 12. Each layer represents one month
#' @param tme vector of POSIX objects for the desired time period
#'
#' @return array of hourly seasonal values in array format where each layer corresponds to one hour.
#' @export
#'
#'
seasonalRast_to_hourlyArray = function (seasonalRast, tme) {
  year1 = as.numeric(format(tme[1], format = "%Y"))
  year2 = as.numeric(format(tme[length(tme)], format = "%Y"))
  years = seq(year1, year2, 1)

  jan = 744
  feb = 672
  feb_leap = 696
  march = 744
  april = 720
  may = 744
  june = 720
  july = 744
  aug = 744
  sep = 720
  oct = 744
  nov = 720
  dec = 744

  months = c(jan, feb, march, april, may, june, july, aug, sep, oct, nov, dec)
  months_leap = c(jan, feb_leap, march, april, may, june, july, aug, sep, oct, nov, dec)

  year_hours = list()
  for (y in 1:length(years)) {
    # initiate list to store hourly pad values in for one year
    mon = list()
    # if a leap year
    if (years[y]%%4 == 0) {
      for (m in 1:length(months)) {
        mon[[m]] = rep(seasonalRast[[m]], each = months_leap[m])
      }
    } else {
      for (m in 1:length(months)) {
        mon[[m]] = rep(seasonalRast[[m]], each = months[m])
      }
    }
    mon = terra::rast(mon)
    #print(dim(mon))
    year_hours = c(year_hours, mon)
  }
  year_hours = terra::rast(year_hours)
  year_hours = as.array(year_hours)
  return(year_hours)
}



#' seasonalRast_to_dailyArray
#'
#'Takes a spatRaster with each layer representing monthly seasonality of a different variable, replicates monthly seasonality to produce hourly, and converts to array.
#'
#' @param seasonalRast raster with nlyrs = 12. Each layer represents one month
#' @param tme vector of POSIX objects for the desired time period
#'
#' @return array of hourly seasonal values in array format where each layer corresponds to one hour.
#' @export
#'
#'
seasonalRast_to_dailyArray = function (seasonalRast, tme) {

  year1 = as.numeric(format(tme[1], format = "%Y"))
  year2 = as.numeric(format(tme[length(tme)], format = "%Y"))
  years = seq(year1, year2, 1)

  jan = 31
  feb = 28
  feb_leap = 29
  march = 31
  april = 30
  may = 31
  june = 30
  july = 31
  aug = 31
  sep = 30
  oct = 31
  nov = 30
  dec = 31

  months = c(jan, feb, march, april, may, june, july, aug, sep, oct, nov, dec)
  months_leap = c(jan, feb_leap, march, april, may, june, july, aug, sep, oct, nov, dec)

  year_days = list()
  for (y in 1:length(years)) {
    # initiate list to store hourly pad values in for one year
    mon = list()
    # if a leap year
    if (years[y]%%4 == 0) {
      for (m in 1:length(months)) {
        mon[[m]] = rep(seasonalRast[[m]], each = months_leap[m])
      }
    } else {
      for (m in 1:length(months)) {
        mon[[m]] = rep(seasonalRast[[m]], each = months[m])
      }
    }
    mon = terra::rast(mon)
    #print(dim(mon))
    year_days = c(year_days, mon)
  }
  year_days = terra::rast(year_days)
  year_days = as.array(year_days)
  return(year_days)
}
