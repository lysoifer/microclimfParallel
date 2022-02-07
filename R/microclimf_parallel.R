#' run_microclimf_parallel
#'
#' makes single spatraster with multiple bands for different summaries
#'
#' @param tile_list list of raster tiles to iterate over
#' @param n_cores integer: number of cores to use for parallel processing
#' @param pai spatRaster of PAI for entire study area
#' @param dem spatRaster of DEM for entire study area
#' @param chm spatRaster of CHM for entire study area
#' @param habitats spatRaster of habitat types (habitat type key in microclimf package) for entire study area
#' @param soil spatRaster of soil types for entire study area (soil type key in microclimf)
#' @param pai_seasonal array of PAI values as returned by lai_seasonality()
#' @param paia_seasonal array of PAI values above reqhgt as returned by paia_for_microf()
#' @param folden_seasonal array of foliage density at reqhgt as returned by folden_for_microf()
#' @param weather dataframe of hourly weather as returned by mcera5::extract_clim
#' @param precip vector of daily precipitation as returned by mcera5::extract_precip
#' @param ll dataframe of lat and lon values (see latlongfromSpatRaster())
#' @param reqhgt height at which to model microclimate
#' @param epsg epsg code
#' @param n_years number of years represented in weather data
#'
#'
#' @return summ_rasts: spatRaster: each layer represents one summary variable
#' @export
#'
#' @importFrom foreach %dopar% foreach
#' @importFrom iterators icount
#' @importFrom doParallel registerDoParallel
#'
#'
run_microclimf_parallel = function(tile_list, pai, dem, chm, habitats, soil, n_cores, pai_seasonal, paia_seasonal, folden_seasonal, weather, precip, ll, reqhgt, epsg, n_years) {

  # Construct cluster
  cl = parallel::makeCluster(n_cores)

  # After the function is run, shutdown the cluster.
  on.exit(parallel::stopCluster(cl))

  # Register parallel backend
  doParallel::registerDoParallel(cl)   # Modify with any do*::registerDo*()

  # Compute estimates
  summ_rasts = foreach (i = c(1:length(tile_list)), .packages = c("terra", "raster", "microclimf", "microctools", "microclimfParallel"),
                        .export = c('tile_list', 'pai', 'dem', 'chm', 'habitats', 'soil', 'pai_seasonal', 'paia_seasonal', 'folden_seasonal', 'weather', 'precip', 'll', 'reqhgt', 'epsg', 'n_years')) %dopar% {
    tile = tile_list[[i]]
    raster::values(tile) = 1

    # extend rasters to prevent NA values between tiles due to edge effect calculations
    tile = raster::extend(tile, c(1,1))
    e = terra::ext(tile)

    pai.tile = raster::raster(terra::crop(pai, e))
    dem.tile = raster::raster(terra::crop(dem, e))
    chm.tile = raster::raster(terra::crop(chm, e))
    habitats.tile = raster::raster(terra::crop(habitats, terra::rast(tile)))
    soil.tile = raster::raster(terra::crop(soil, terra::rast(tile)))
    # soil data missing in some areas. Since northern range is mostly clay, make NA clay value (11)
    if(raster::freq(soil.tile, value = NA) == raster::ncell(soil.tile)) {
      raster::values(soil.tile) = 11
    }

    # test for missing data
    pai_test = raster::freq(pai.tile, value = NA) == raster::ncell(pai.tile)
    dem_test = raster::freq(dem.tile, value = NA) == raster::ncell(dem.tile)
    chm_test = raster::freq(chm.tile, value = NA) == raster::ncell(chm.tile)
    habitat_test = raster::freq(habitats.tile, value = NA) == raster::ncell(habitats.tile)

    # check if area is entirely NA (i.e. in the ocean)
    # if area is entirely in ocean or there is missing data, add an NA tile to list of climate variables
    if (pai_test | dem_test | chm_test | habitat_test) {
      placeholder = terra::rast(raster::trim(tile))
      placeholder[placeholder==1] = NA

      # assign placeholder raster is NA
      avg_daily_mean_temp = avg_daytime_max_temp = avg_daytime_min_temp = avg_daytime_mean_temp = avg_night_max_temp=
      avg_night_min_temp= avg_night_mean_temp= avg_daily_mean_hum= avg_daytime_max_hum= avg_daytime_min_hum = avg_daytime_mean_hum =
      avg_night_max_hum = avg_night_min_hum = avg_night_mean_hum = avg_daily_mean_groundtemp = avg_daytime_max_groundtemp =
      avg_daytime_min_groundtemp = avg_daytime_mean_groundtemp = avg_night_max_groundtemp = avg_night_min_groundtemp =
      avg_night_mean_groundtemp = avg_daily_raddir = avg_daily_raddif = avg_daily_radlw = avg_daily_mean_temp_dry = avg_daytime_max_temp_dry =
      avg_daytime_min_temp_dry = avg_daytime_mean_temp_dry = avg_night_max_temp_dry = avg_night_min_temp_dry = avg_night_mean_temp_dry =
      avg_daily_mean_hum_dry = avg_daytime_max_hum_dry = avg_daytime_min_hum_dry = avg_daytime_mean_hum_dry = avg_night_max_hum_dry =
      avg_night_min_hum_dry = avg_night_mean_hum_dry = avg_daily_mean_groundtemp_dry = avg_daytime_max_groundtemp_dry = avg_daytime_min_groundtemp_dry =
      avg_daytime_mean_groundtemp_dry = avg_night_max_groundtemp_dry = avg_night_min_groundtemp_dry = avg_night_mean_groundtemp_dry =
      avg_daily_mean_temp_wet = avg_daytime_max_temp_wet = avg_daytime_min_temp_wet = avg_daytime_mean_temp_wet = avg_night_max_temp_wet =
      avg_night_min_temp_wet = avg_night_mean_temp_wet = avg_daily_mean_hum_wet = avg_daytime_max_hum_wet = avg_daytime_min_hum_wet = avg_daytime_mean_hum_wet =
      avg_night_max_hum_wet = avg_night_min_hum_wet = avg_night_mean_hum_wet = avg_daily_mean_groundtemp_wet = avg_daytime_max_groundtemp_wet =
      avg_daytime_min_groundtemp_wet = avg_daytime_mean_groundtemp_wet = avg_night_max_groundtemp_wet = avg_night_min_groundtemp_wet = avg_night_mean_groundtemp_wet = placeholder


      print(paste(i, 'skipping ocean tile'))
    } else {

      # crop seasonal spatRasters separately so they don't have to be converted to rasterStacks since they will become arrays anyway
      pai_seasonal.tile = terra::crop(pai_seasonal, terra::rast(tile))
      paia_seasonal.tile = terra::crop(paia_seasonal, terra::rast(tile))
      folden_seasonal.tile = terra::crop(folden_seasonal, terra::rast(tile))

      # resample rasters so resolution matches pai resolution (soil calculated at pai resolution, so no need to resample)
      #dem.tile = raster::resample(dem.tile, pai.tile, method = "bilinear")
      #chm.tile = raster::resample(chm.tile, pai.tile, method = "bilinear")

      # prepare inputs to modelin
      climdata = microctools::hourlyncep_convert(weather, lat = ll$lat, long = ll$long)

      #rainfall = precip[[i]]

      # make pai and paia array
      pai.array = as.array(pai_seasonal.tile)

      # replicate array for n_years
      pai.array = array(rep(pai.array, n_years), dim = c(dim(pai.array)[1], dim(pai.array)[2], dim(pai.array)[3]*n_years))

      # estimate vegetation parameters from chm and pai
      vegp = microclimf::vegpfromhab(habitats.tile, hgts = chm.tile, pai = pai.array)

      # soil parameters
      soilc = microclimf::soilcfromtype(soil.tile)

      # model inputs
      micro = microclimf::modelin_dy(climdata, precip, vegp, soilc, dem.tile)

      # repeat seasonal pai_a and folden for number of days
      paia.array = seasonalRast_to_dailyArray(paia_seasonal.tile, climdata$obs_time)
      folden.array = seasonalRast_to_dailyArray(folden_seasonal.tile, climdata$obs_time)

      # run model

      model_dy = microclimf::runmicro_dy(micro, reqhgt, expand = T,
                             pai_a = paia.array, folden = folden.array)
      print(paste(i, 'model finished'))


      mout_r = microout_to_raster(model_dy, tile, epsg, climvars = c('Tz', 'relhum', 'T0', 'raddir', 'raddif', 'radlw'))

      print(paste(i, 'microout_to_raster finished'))

      #avg_daily_mean_temp = merge_tiles(avg_daily_mean_temp)
      avg_daily_mean_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean')
      avg_daytime_max_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(6,18))
      avg_daytime_min_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(6,18))
      avg_daytime_mean_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(6,18))
      avg_night_max_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(0:5,19:23))
      avg_night_min_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(0:5,19:23))
      avg_night_mean_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(0:5,19:23))

      # Avg. daily max/min/mean relative humidity
      avg_daily_mean_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean')
      avg_daytime_max_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(6,18))
      avg_daytime_min_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(6,18))
      avg_daytime_mean_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(6,18))
      avg_night_max_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(0:5,19:23))
      avg_night_min_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(0:5,19:23))
      avg_night_mean_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(0:5,19:23))

      # Avg. daily max/min/mean  GROUND TEMP (TO)
      avg_daily_mean_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean')
      avg_daytime_max_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(6,18))
      avg_daytime_min_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(6,18))
      avg_daytime_mean_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(6,18))
      avg_night_max_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(0:5,19:23))
      avg_night_min_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(0:5,19:23))
      avg_night_mean_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(0:5,19:23))

      # AVG. TOTAL DAILY DIRECT RADIATION
      avg_daily_raddir = summarise_microclimf(mout_r$raddir, climdata = climdata, fun = 'sum')
      # AVG. TOTAL DAILY DIFFUSE RADIATION
      avg_daily_raddif = summarise_microclimf(mout_r$raddif, climdata = climdata, fun = 'sum')
      # AVG. TOTAL DAILY LONGWAVE RADIATION
      avg_daily_radlw = summarise_microclimf(mout_r$radlw, climdata = climdata, fun = 'sum')

      # DRY SEASON
      # Avg. daily max/min/mean temp
      avg_daily_mean_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', months = c(1:5))
      avg_daytime_max_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(6,18), months = c(1:5))
      avg_daytime_min_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(6,18), months = c(1:5))
      avg_daytime_mean_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(1:5))
      avg_night_max_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(1:5))
      avg_night_min_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(1:5))
      avg_night_mean_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(1:5))

      # Avg. daily max/min/mean relative humidity
      avg_daily_mean_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', months = c(1:5))
      avg_daytime_max_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(6,18), months = c(1:5))
      avg_daytime_min_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(6,18), months = c(1:5))
      avg_daytime_mean_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(1:5))
      avg_night_max_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(1:5))
      avg_night_min_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(1:5))
      avg_night_mean_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(1:5))

      # Avg. daily max/min/mean relative GROUND TEMP (TO)
      avg_daily_mean_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', months = c(1:5))
      avg_daytime_max_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(6,18), months = c(1:5))
      avg_daytime_min_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(6,18), months = c(1:5))
      avg_daytime_mean_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(1:5))
      avg_night_max_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(1:5))
      avg_night_min_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(1:5))
      avg_night_mean_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(1:5))

      # WET SEASON
      avg_daily_mean_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', months = c(6:12))
      avg_daytime_max_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(6,18), months = c(6:12))
      avg_daytime_min_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(6,18), months = c(6:12))
      avg_daytime_mean_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(6:12))
      avg_night_max_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(6:12))
      avg_night_min_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(6:12))
      avg_night_mean_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(6:12))

      # Avg. daily max/min/mean relative humidity
      avg_daily_mean_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', months = c(6:12))
      avg_daytime_max_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(6,18), months = c(6:12))
      avg_daytime_min_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(6,18), months = c(6:12))
      avg_daytime_mean_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(6:12))
      avg_night_max_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(6:12))
      avg_night_min_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(6:12))
      avg_night_mean_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(6:12))

      # Avg. daily max/min/mean relative GROUND TEMP (TO)
      avg_daily_mean_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', months = c(6:12))
      avg_daytime_max_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(6,18), months = c(6:12))
      avg_daytime_min_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(6,18), months = c(6:12))
      avg_daytime_mean_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(6:12))
      avg_night_max_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(6:12))
      avg_night_min_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(6:12))
      avg_night_mean_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(6:12))

    }
    # store summary rasters as layered spatRaster
    summ_rasts = c(avg_daily_mean_temp, avg_daytime_max_temp, avg_daytime_min_temp, avg_daytime_mean_temp, avg_night_max_temp,
                   avg_night_min_temp, avg_night_mean_temp, avg_daily_mean_hum, avg_daytime_max_hum, avg_daytime_min_hum, avg_daytime_mean_hum,
                   avg_night_max_hum, avg_night_min_hum, avg_night_mean_hum, avg_daily_mean_groundtemp, avg_daytime_max_groundtemp,
                   avg_daytime_min_groundtemp, avg_daytime_mean_groundtemp, avg_night_max_groundtemp, avg_night_min_groundtemp,
                   avg_night_mean_groundtemp, avg_daily_raddir, avg_daily_raddif, avg_daily_radlw, avg_daily_mean_temp_dry, avg_daytime_max_temp_dry,
                   avg_daytime_min_temp_dry, avg_daytime_mean_temp_dry, avg_night_max_temp_dry, avg_night_min_temp_dry, avg_night_mean_temp_dry,
                   avg_daily_mean_hum_dry, avg_daytime_max_hum_dry, avg_daytime_min_hum_dry, avg_daytime_mean_hum_dry, avg_night_max_hum_dry,
                   avg_night_min_hum_dry, avg_night_mean_hum_dry, avg_daily_mean_groundtemp_dry, avg_daytime_max_groundtemp_dry, avg_daytime_min_groundtemp_dry,
                   avg_daytime_mean_groundtemp_dry, avg_night_max_groundtemp_dry, avg_night_min_groundtemp_dry, avg_night_mean_groundtemp_dry,
                   avg_daily_mean_temp_wet, avg_daytime_max_temp_wet, avg_daytime_min_temp_wet, avg_daytime_mean_temp_wet, avg_night_max_temp_wet,
                   avg_night_min_temp_wet, avg_night_mean_temp_wet, avg_daily_mean_hum_wet, avg_daytime_max_hum_wet, avg_daytime_min_hum_wet, avg_daytime_mean_hum_wet,
                   avg_night_max_hum_wet, avg_night_min_hum_wet, avg_night_mean_hum_wet, avg_daily_mean_groundtemp_wet, avg_daytime_max_groundtemp_wet,
                   avg_daytime_min_groundtemp_wet, avg_daytime_mean_groundtemp_wet, avg_night_max_groundtemp_wet, avg_night_min_groundtemp_wet, avg_night_mean_groundtemp_wet)

    names(summ_rasts) = c('avg_daily_mean_temp', 'avg_daytime_max_temp', 'avg_daytime_min_temp', 'avg_daytime_mean_temp', 'avg_night_max_temp',
                          'avg_night_min_temp', 'avg_night_mean_temp', 'avg_daily_mean_hum', 'avg_daytime_max_hum', 'avg_daytime_min_hum', 'avg_daytime_mean_hum',
                          'avg_night_max_hum', 'avg_night_min_hum', 'avg_night_mean_hum', 'avg_daily_mean_groundtemp', 'avg_daytime_max_groundtemp',
                          'avg_daytime_min_groundtemp', 'avg_daytime_mean_groundtemp', 'avg_night_max_groundtemp', 'avg_night_min_groundtemp',
                          'avg_night_mean_groundtemp', 'avg_daily_raddir', 'avg_daily_raddif', 'avg_daily_radlw', 'avg_daily_mean_temp_dry', 'avg_daytime_max_temp_dry',
                          'avg_daytime_min_temp_dry', 'avg_daytime_mean_temp_dry', 'avg_night_max_temp_dry', 'avg_night_min_temp_dry', 'avg_night_mean_temp_dry',
                          'avg_daily_mean_hum_dry', 'avg_daytime_max_hum_dry', 'avg_daytime_min_hum_dry', 'avg_daytime_mean_hum_dry', 'avg_night_max_hum_dry',
                          'avg_night_min_hum_dry', 'avg_night_mean_hum_dry', 'avg_daily_mean_groundtemp_dry', 'avg_daytime_max_groundtemp_dry', 'avg_daytime_min_groundtemp_dry',
                          'avg_daytime_mean_groundtemp_dry', 'avg_night_max_groundtemp_dry', 'avg_night_min_groundtemp_dry', 'avg_night_mean_groundtemp_dry',
                          'avg_daily_mean_temp_wet', 'avg_daytime_max_temp_wet', 'avg_daytime_min_temp_wet', 'avg_daytime_mean_temp_wet', 'avg_night_max_temp_wet',
                          'avg_night_min_temp_wet', 'avg_night_mean_temp_wet', 'avg_daily_mean_hum_wet', 'avg_daytime_max_hum_wet', 'avg_daytime_min_hum_wet', 'avg_daytime_mean_hum_wet',
                          'avg_night_max_hum_wet', 'avg_night_min_hum_wet', 'avg_night_mean_hum_wet', 'avg_daily_mean_groundtemp_wet', 'avg_daytime_max_groundtemp_wet',
                          'avg_daytime_min_groundtemp_wet', 'avg_daytime_mean_groundtemp_wet', 'avg_night_max_groundtemp_wet', 'avg_night_min_groundtemp_wet', 'avg_night_mean_groundtemp_wet')
    print(paste(i, 'finished'))
    return(summ_rasts)
  }
  return(summ_rasts)
}







#' run_microclimf_parallel_update
#'
#' makes single spatraster with multiple bands for different summaries
#'
#' @param tile list of raster tiles to iterate over
#' @param pai spatRaster of PAI for entire study area
#' @param dem spatRaster of DEM for entire study area
#' @param chm spatRaster of CHM for entire study area
#' @param habitats spatRaster of habitat types (habitat type key in microclimf package) for entire study area
#' @param soil spatRaster of soil types for entire study area (soil type key in microclimf)
#' @param pai_seasonal array of PAI values as returned by lai_seasonality()
#' @param paia_seasonal array of PAI values above reqhgt as returned by paia_for_microf()
#' @param folden_seasonal array of foliage density at reqhgt as returned by folden_for_microf()
#' @param weather dataframe of hourly weather as returned by mcera5::extract_clim
#' @param precip vector of daily precipitation as returned by mcera5::extract_precip
#' @param ll dataframe of lat and lon values (see latlongfromSpatRaster())
#' @param reqhgt height at which to model microclimate
#' @param epsg epsg code
#' @param n_years number of years represented in weather data
#'
#'
#' @return summ_rasts: spatRaster: each layer represents one summary variable
#' @export
#'
#' @importFrom foreach %dopar% foreach
#' @importFrom iterators icount
#' @importFrom doParallel registerDoParallel
#'
#'
run_microclimf_parallel_update = function(tile, pai, dem, chm, habitats, soil, pai_seasonal, paia_seasonal, folden_seasonal, weather, precip, ll, reqhgt, epsg, n_years) {

  # Compute estimates
  raster::values(tile) = 1

  # extend rasters to prevent NA values between tiles due to edge effect calculations
  tile = raster::extend(tile, c(1,1))
  e = terra::ext(tile)

  pai.tile = raster::raster(terra::crop(pai, e))
  dem.tile = raster::raster(terra::crop(dem, e))
  chm.tile = raster::raster(terra::crop(chm, e))
  habitats.tile = raster::raster(terra::crop(habitats, terra::rast(tile)))
  soil.tile = raster::raster(terra::crop(soil, terra::rast(tile)))
  # soil data missing in some areas. Since northern range is mostly clay, make NA clay value (11)
  if(raster::freq(soil.tile, value = NA) == raster::ncell(soil.tile)) {
    raster::values(soil.tile) = 11
  }

  # test for missing data
  pai_test = raster::freq(pai.tile, value = NA) == raster::ncell(pai.tile)
  dem_test = raster::freq(dem.tile, value = NA) == raster::ncell(dem.tile)
  chm_test = raster::freq(chm.tile, value = NA) == raster::ncell(chm.tile)
  habitat_test = raster::freq(habitats.tile, value = NA) == raster::ncell(habitats.tile)

  # check if area is entirely NA (i.e. in the ocean)
  # if area is entirely in ocean or there is missing data, add an NA tile to list of climate variables
  if (pai_test | dem_test | chm_test | habitat_test) {
    placeholder = terra::rast(raster::trim(tile))
    placeholder[placeholder==1] = NA

    # assign placeholder raster is NA
    avg_daily_mean_temp = avg_daytime_max_temp = avg_daytime_min_temp = avg_daytime_mean_temp = avg_night_max_temp=
      avg_night_min_temp= avg_night_mean_temp= avg_daily_mean_hum= avg_daytime_max_hum= avg_daytime_min_hum = avg_daytime_mean_hum =
      avg_night_max_hum = avg_night_min_hum = avg_night_mean_hum = avg_daily_mean_groundtemp = avg_daytime_max_groundtemp =
      avg_daytime_min_groundtemp = avg_daytime_mean_groundtemp = avg_night_max_groundtemp = avg_night_min_groundtemp =
      avg_night_mean_groundtemp = avg_daily_raddir = avg_daily_raddif = avg_daily_radlw = avg_daily_mean_temp_dry = avg_daytime_max_temp_dry =
      avg_daytime_min_temp_dry = avg_daytime_mean_temp_dry = avg_night_max_temp_dry = avg_night_min_temp_dry = avg_night_mean_temp_dry =
      avg_daily_mean_hum_dry = avg_daytime_max_hum_dry = avg_daytime_min_hum_dry = avg_daytime_mean_hum_dry = avg_night_max_hum_dry =
      avg_night_min_hum_dry = avg_night_mean_hum_dry = avg_daily_mean_groundtemp_dry = avg_daytime_max_groundtemp_dry = avg_daytime_min_groundtemp_dry =
      avg_daytime_mean_groundtemp_dry = avg_night_max_groundtemp_dry = avg_night_min_groundtemp_dry = avg_night_mean_groundtemp_dry =
      avg_daily_mean_temp_wet = avg_daytime_max_temp_wet = avg_daytime_min_temp_wet = avg_daytime_mean_temp_wet = avg_night_max_temp_wet =
      avg_night_min_temp_wet = avg_night_mean_temp_wet = avg_daily_mean_hum_wet = avg_daytime_max_hum_wet = avg_daytime_min_hum_wet = avg_daytime_mean_hum_wet =
      avg_night_max_hum_wet = avg_night_min_hum_wet = avg_night_mean_hum_wet = avg_daily_mean_groundtemp_wet = avg_daytime_max_groundtemp_wet =
      avg_daytime_min_groundtemp_wet = avg_daytime_mean_groundtemp_wet = avg_night_max_groundtemp_wet = avg_night_min_groundtemp_wet = avg_night_mean_groundtemp_wet = placeholder

  } else {

    # crop seasonal spatRasters separately so they don't have to be converted to rasterStacks since they will become arrays anyway
    pai_seasonal.tile = terra::crop(pai_seasonal, terra::rast(tile))
    paia_seasonal.tile = terra::crop(paia_seasonal, terra::rast(tile))
    folden_seasonal.tile = terra::crop(folden_seasonal, terra::rast(tile))

    # resample rasters so resolution matches pai resolution (soil calculated at pai resolution, so no need to resample)
    #dem.tile = raster::resample(dem.tile, pai.tile, method = "bilinear")
    #chm.tile = raster::resample(chm.tile, pai.tile, method = "bilinear")

    # prepare inputs to modelin
    climdata = microctools::hourlyncep_convert(weather, lat = ll$lat, long = ll$long)

    #rainfall = precip[[i]]

    # make pai and paia array
    pai.array = as.array(pai_seasonal.tile)

    # replicate array for n_years
    pai.array = array(rep(pai.array, n_years), dim = c(dim(pai.array)[1], dim(pai.array)[2], dim(pai.array)[3]*n_years))

    # estimate vegetation parameters from chm and pai
    vegp = microclimf::vegpfromhab(habitats.tile, hgts = chm.tile, pai = pai.array)

    # soil parameters
    soilc = microclimf::soilcfromtype(soil.tile)

    # model inputs
    micro = microclimf::modelin_dy(climdata, precip, vegp, soilc, dem.tile)

    # repeat seasonal pai_a and folden for number of days
    paia.array = seasonalRast_to_dailyArray(paia_seasonal.tile, climdata$obs_time)
    folden.array = seasonalRast_to_dailyArray(folden_seasonal.tile, climdata$obs_time)

    # run model

    model_dy = microclimf::runmicro_dy(micro, reqhgt, expand = T,
                                       pai_a = paia.array, folden = folden.array)
    print(paste(i, 'model finished'))


    mout_r = microout_to_raster(model_dy, tile, epsg, climvars = c('Tz', 'relhum', 'T0', 'raddir', 'raddif', 'radlw'))

    print(paste(i, 'microout_to_raster finished'))

    #avg_daily_mean_temp = merge_tiles(avg_daily_mean_temp)
    avg_daily_mean_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean')
    avg_daytime_max_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(6,18))
    avg_daytime_min_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(6,18))
    avg_daytime_mean_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(6,18))
    avg_night_max_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(0:5,19:23))
    avg_night_min_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(0:5,19:23))
    avg_night_mean_temp = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(0:5,19:23))

    # Avg. daily max/min/mean relative humidity
    avg_daily_mean_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean')
    avg_daytime_max_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(6,18))
    avg_daytime_min_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(6,18))
    avg_daytime_mean_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(6,18))
    avg_night_max_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(0:5,19:23))
    avg_night_min_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(0:5,19:23))
    avg_night_mean_hum = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(0:5,19:23))

    # Avg. daily max/min/mean  GROUND TEMP (TO)
    avg_daily_mean_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean')
    avg_daytime_max_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(6,18))
    avg_daytime_min_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(6,18))
    avg_daytime_mean_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(6,18))
    avg_night_max_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(0:5,19:23))
    avg_night_min_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(0:5,19:23))
    avg_night_mean_groundtemp = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(0:5,19:23))

    # AVG. TOTAL DAILY DIRECT RADIATION
    avg_daily_raddir = summarise_microclimf(mout_r$raddir, climdata = climdata, fun = 'sum')
    # AVG. TOTAL DAILY DIFFUSE RADIATION
    avg_daily_raddif = summarise_microclimf(mout_r$raddif, climdata = climdata, fun = 'sum')
    # AVG. TOTAL DAILY LONGWAVE RADIATION
    avg_daily_radlw = summarise_microclimf(mout_r$radlw, climdata = climdata, fun = 'sum')

    # DRY SEASON
    # Avg. daily max/min/mean temp
    avg_daily_mean_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', months = c(1:5))
    avg_daytime_max_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(6,18), months = c(1:5))
    avg_daytime_min_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(6,18), months = c(1:5))
    avg_daytime_mean_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(1:5))
    avg_night_max_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(1:5))
    avg_night_min_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(1:5))
    avg_night_mean_temp_dry = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(1:5))

    # Avg. daily max/min/mean relative humidity
    avg_daily_mean_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', months = c(1:5))
    avg_daytime_max_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(6,18), months = c(1:5))
    avg_daytime_min_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(6,18), months = c(1:5))
    avg_daytime_mean_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(1:5))
    avg_night_max_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(1:5))
    avg_night_min_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(1:5))
    avg_night_mean_hum_dry = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(1:5))

    # Avg. daily max/min/mean relative GROUND TEMP (TO)
    avg_daily_mean_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', months = c(1:5))
    avg_daytime_max_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(6,18), months = c(1:5))
    avg_daytime_min_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(6,18), months = c(1:5))
    avg_daytime_mean_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(1:5))
    avg_night_max_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(1:5))
    avg_night_min_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(1:5))
    avg_night_mean_groundtemp_dry = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(1:5))

    # WET SEASON
    avg_daily_mean_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', months = c(6:12))
    avg_daytime_max_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(6,18), months = c(6:12))
    avg_daytime_min_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(6,18), months = c(6:12))
    avg_daytime_mean_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(6:12))
    avg_night_max_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(6:12))
    avg_night_min_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(6:12))
    avg_night_mean_temp_wet = summarise_microclimf(mout_r$Tz, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(6:12))

    # Avg. daily max/min/mean relative humidity
    avg_daily_mean_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', months = c(6:12))
    avg_daytime_max_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(6,18), months = c(6:12))
    avg_daytime_min_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(6,18), months = c(6:12))
    avg_daytime_mean_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(6:12))
    avg_night_max_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(6:12))
    avg_night_min_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(6:12))
    avg_night_mean_hum_wet = summarise_microclimf(mout_r$relhum, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(6:12))

    # Avg. daily max/min/mean relative GROUND TEMP (TO)
    avg_daily_mean_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', months = c(6:12))
    avg_daytime_max_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(6,18), months = c(6:12))
    avg_daytime_min_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(6,18), months = c(6:12))
    avg_daytime_mean_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(6,18), months = c(6:12))
    avg_night_max_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'max', hours = c(0:5,19:23), months = c(6:12))
    avg_night_min_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'min', hours = c(0:5,19:23), months = c(6:12))
    avg_night_mean_groundtemp_wet = summarise_microclimf(mout_r$T0, climdata = climdata, fun = 'mean', hours = c(0:5,19:23), months = c(6:12))

  }
  # store summary rasters as layered spatRaster
  summ_rasts = c(avg_daily_mean_temp, avg_daytime_max_temp, avg_daytime_min_temp, avg_daytime_mean_temp, avg_night_max_temp,
                 avg_night_min_temp, avg_night_mean_temp, avg_daily_mean_hum, avg_daytime_max_hum, avg_daytime_min_hum, avg_daytime_mean_hum,
                 avg_night_max_hum, avg_night_min_hum, avg_night_mean_hum, avg_daily_mean_groundtemp, avg_daytime_max_groundtemp,
                 avg_daytime_min_groundtemp, avg_daytime_mean_groundtemp, avg_night_max_groundtemp, avg_night_min_groundtemp,
                 avg_night_mean_groundtemp, avg_daily_raddir, avg_daily_raddif, avg_daily_radlw, avg_daily_mean_temp_dry, avg_daytime_max_temp_dry,
                 avg_daytime_min_temp_dry, avg_daytime_mean_temp_dry, avg_night_max_temp_dry, avg_night_min_temp_dry, avg_night_mean_temp_dry,
                 avg_daily_mean_hum_dry, avg_daytime_max_hum_dry, avg_daytime_min_hum_dry, avg_daytime_mean_hum_dry, avg_night_max_hum_dry,
                 avg_night_min_hum_dry, avg_night_mean_hum_dry, avg_daily_mean_groundtemp_dry, avg_daytime_max_groundtemp_dry, avg_daytime_min_groundtemp_dry,
                 avg_daytime_mean_groundtemp_dry, avg_night_max_groundtemp_dry, avg_night_min_groundtemp_dry, avg_night_mean_groundtemp_dry,
                 avg_daily_mean_temp_wet, avg_daytime_max_temp_wet, avg_daytime_min_temp_wet, avg_daytime_mean_temp_wet, avg_night_max_temp_wet,
                 avg_night_min_temp_wet, avg_night_mean_temp_wet, avg_daily_mean_hum_wet, avg_daytime_max_hum_wet, avg_daytime_min_hum_wet, avg_daytime_mean_hum_wet,
                 avg_night_max_hum_wet, avg_night_min_hum_wet, avg_night_mean_hum_wet, avg_daily_mean_groundtemp_wet, avg_daytime_max_groundtemp_wet,
                 avg_daytime_min_groundtemp_wet, avg_daytime_mean_groundtemp_wet, avg_night_max_groundtemp_wet, avg_night_min_groundtemp_wet, avg_night_mean_groundtemp_wet)

  names(summ_rasts) = c('avg_daily_mean_temp', 'avg_daytime_max_temp', 'avg_daytime_min_temp', 'avg_daytime_mean_temp', 'avg_night_max_temp',
                        'avg_night_min_temp', 'avg_night_mean_temp', 'avg_daily_mean_hum', 'avg_daytime_max_hum', 'avg_daytime_min_hum', 'avg_daytime_mean_hum',
                        'avg_night_max_hum', 'avg_night_min_hum', 'avg_night_mean_hum', 'avg_daily_mean_groundtemp', 'avg_daytime_max_groundtemp',
                        'avg_daytime_min_groundtemp', 'avg_daytime_mean_groundtemp', 'avg_night_max_groundtemp', 'avg_night_min_groundtemp',
                        'avg_night_mean_groundtemp', 'avg_daily_raddir', 'avg_daily_raddif', 'avg_daily_radlw', 'avg_daily_mean_temp_dry', 'avg_daytime_max_temp_dry',
                        'avg_daytime_min_temp_dry', 'avg_daytime_mean_temp_dry', 'avg_night_max_temp_dry', 'avg_night_min_temp_dry', 'avg_night_mean_temp_dry',
                        'avg_daily_mean_hum_dry', 'avg_daytime_max_hum_dry', 'avg_daytime_min_hum_dry', 'avg_daytime_mean_hum_dry', 'avg_night_max_hum_dry',
                        'avg_night_min_hum_dry', 'avg_night_mean_hum_dry', 'avg_daily_mean_groundtemp_dry', 'avg_daytime_max_groundtemp_dry', 'avg_daytime_min_groundtemp_dry',
                        'avg_daytime_mean_groundtemp_dry', 'avg_night_max_groundtemp_dry', 'avg_night_min_groundtemp_dry', 'avg_night_mean_groundtemp_dry',
                        'avg_daily_mean_temp_wet', 'avg_daytime_max_temp_wet', 'avg_daytime_min_temp_wet', 'avg_daytime_mean_temp_wet', 'avg_night_max_temp_wet',
                        'avg_night_min_temp_wet', 'avg_night_mean_temp_wet', 'avg_daily_mean_hum_wet', 'avg_daytime_max_hum_wet', 'avg_daytime_min_hum_wet', 'avg_daytime_mean_hum_wet',
                        'avg_night_max_hum_wet', 'avg_night_min_hum_wet', 'avg_night_mean_hum_wet', 'avg_daily_mean_groundtemp_wet', 'avg_daytime_max_groundtemp_wet',
                        'avg_daytime_min_groundtemp_wet', 'avg_daytime_mean_groundtemp_wet', 'avg_night_max_groundtemp_wet', 'avg_night_min_groundtemp_wet', 'avg_night_mean_groundtemp_wet')

  return(summ_rasts)
}
