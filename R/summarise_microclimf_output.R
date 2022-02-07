#' microout_to_raster
#'
#' converts output of microclimf models into raster format
#'
#' @param microout output of run_micro or run_micro_dy
#' @param r raster tile with NA buffer of 1 row/column around the edge
#' @param epsg epsg code for output
#' @param climvars vector of climate variables to extract from microout (c('Tz', 'tleaf', 'T0', 'soilm', 'relhum', 'windspeed', 'raddir', 'raddif', 'radlw'))
#'
#' @return named list of spatRasters
#' @export
#'
#'
microout_to_raster = function (microout, r, epsg, climvars) {

  r = terra::rast(r)
  microout_r = list()
  for (i in climvars) {
    m = microout[[i]]
    m = m[-c(1, nrow(m)), -c(1,ncol(m)),]
    rast = terra::rast(m)
    terra::ext(rast) = terra::ext(raster::trim(r))
    terra::crs(rast) = paste0('epsg:', epsg)
    #rast = crop(x = rast, y = trim(r))
    #print('completed crop')
    microout_r[[i]] = rast
  }
  names(microout_r) = climvars
  return(microout_r)
}


#' summarise_microclimf
#'
#' summarizes clim_rast by designated time period. The function first subsets the data by months and hours
#' if provided, and then summarizes by day. All daily summaries are then averaged.
#' For example fun = max, will provide average daily maximum temperature
#'
#' @param clim_rast spatRaster of a climate variable. Layers must represent times used to calculate microclimates
#' @param climdata climate dataframe that was used as input in microclimf
#' @param fun string: summary function to apply to clim_array (see fun options in terra::app)
#' @param months Months over which to apply summary function (1-12). If NA, summary function applied over all months
#' @param hours hours in day over which to apply summary function (0-23). If NA, summary function applied over all hours
#'
#' @return summary raster
#' @export
#'
#'
summarise_microclimf = function(clim_rast, climdata, fun, months = NA, hours = NA) {

  # check if clim_rast is entirely NA values
  if (terra::freq(clim_rast, value = NA)[3]==terra::ncell(clim_rast)) {
    return(clim_rast)
  }

  # name third dimension of array with time stamp
  timenames = format(strptime(climdata$obs_time, format = "%Y-%m-%d %H:%M:%S", tz = 'GMT'), format = "%Y-mo%m-%d hr%H")
  names(clim_rast) = timenames

  # subset desired months if specified
  if (!is.na(months)[1]) {
    m = sprintf("%02d", months)
    m = paste0('mo', m, collapse = '|')
    m.index = grep(m, names(clim_rast)) # identify indices of designated hours to include in summary
    clim_rast = subset(clim_rast, m.index) # subset designated hours
  }
  if (!is.na(hours)[1]) {
    h = sprintf("%02d", hours)
    h = paste0('hr', h, collapse = '|')
    h.index = grep(h, names(clim_rast)) # identify indices of designated hours to include in summary
    clim_rast = subset(clim_rast, h.index) # subset designated hours
  }
  # rename by year and julian day
  timenames2 = format(strptime(names(clim_rast), format = "%Y-mo%m-%d hr%H", tz = 'GMT'), format = "%Y-%j")
  names(clim_rast) = timenames2

  # assign values to unique year_julianday combinations
  i = base::match(names(clim_rast), unique(names(clim_rast)))

  # apply specified function for each day
  day_summ = terra::tapp(clim_rast, index = i, fun = fun)
  names(day_summ) = unique(names(clim_rast))

  # # identify unique days for each year
  # yr_day = unique(names(clim_rast))
  # summ = vector(mode = 'list', length = length(yr_day))
  # names(summ) = yr_day
  #
  # # summarise for each day
  # for (y in yr_day) {
  #   sub = clim_rast[y]
  #   # calculate designated fun for each year
  #   summ[[y]] = app(sub, cores = 4, fun, na.rm = T)
  # }

  # average over all years/days in the list
  #yrs_summ = rast(summ)
  yrs_summ = terra::app(day_summ, cores = 4, fun = 'mean', na.rm = T)
  return(yrs_summ)

}
