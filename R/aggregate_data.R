#' Aggregating netCDF dataset.
#'
#' Aggregate netCDF dataset to monthly or seasonal.
#' 
#' @param in_time_res character, the time resolution of the input netCDF dataset. Must be one of \code{'daily'}, \code{'dekadal'} and \code{'monthly'}.
#' @param out_time_res character, the time resolution of the aggregated data. Must be one of \code{'monthly'}, \code{'seasonal'} and \code{'annual'}.
#' @param start_date character, the start date of the input netCDF from which the data will be aggregate. The date has the following format:
#' \itemize{
#' \item \code{'daily'}: \code{'yyyy-mm-dd'}, where \code{yyyy} is the year in 4 digits, \code{mm} the month in 2 digits, and \code{dd} the day in 2 digits.
#' \item \code{'dekadal'}: \code{'yyyy-mm-d'}, where \code{d} is representing the dekad, it must take the value 1, 2 or 3.
#' \item \code{'monthly'}: \code{'yyyy-mm'}
#' }
#' @param end_date character, the end date of the input netCDF to which the data will be aggregate. Same format as \code{start_date}. 
#' @param in_nc_path_dir character, full path to the folder containing the input netCDF dataset.
#' @param in_nc_file_fromat character, format of the input netCDF file name, the year, month, dekad and ady must be replaced by \code{\%s}.
#' The date format of the netCDF file names must be in the following order: year goes first, then month, and finally dekad or day.
#' Example: for \emph{ersst.v5.202401.nc} the file name format should be \emph{'ersst.v5.\%s\%s.nc'}.
#' @param out_nc_path_dir character, full path to the folder to store the aggregated data.
#' @param out_nc_file_prefix character, a string to be used as a prefix for the aggregated netCDF file names. Default \code{NULL}, the variable name will be used as file name prefix.
#' @param out_nc_variable a named list to define the output netCDF variable, a list of arguments to be passed on the function \code{\link[ncdf4]{ncvar_def}} of the package \pkg{ncdf4}.
#' Default \code{NULL}, the output netCDF variable will be created from the input netCDF variable.
#' @param aggregate_fun a function to be used to aggregate the data.
#' @param ... other arguments that may be passed on to the function \code{aggregate_fun}.
#' @param season_length integer, the length of the season to be aggregated in case of \code{out_time_res = 'seasonal'}. It must be an integer range from 2 to 12.
#' @param min_frac float, the minimum fraction of non missing data that must be present within each output time resolution.
#' @param var_name character, name of the variable to be read from the input netCDF dataset. If unspecified or left \code{NA},
#' the first variable available in the netCDF data will be taken.
#' @param lon_name character, the name of the longitude dimension from the input netCDF. If unspecified or left \code{NA}, it will be detected. 
#' @param lat_name character, the name of the latitude dimension from the input netCDF. Same as \code{lon_name}.
#' @param extra_dim a named list of pairs of arguments, providing the names and values of any extra dimensions to be extracted.
#' The name of the elements of the list should be in the format \code{<information about the dimension>_name} for the name of the dimension
#' and \code{<information about the dimension>_value} for the value to be extracted from that dimension.
#' Example: if your netCDF data has a pressure levels dimension named \code{"pres"} having the values from 1000 to 10 hPa,
#' if you want to extract the data at the pressure level 850 hPa, 
#' \code{extra_dim} would be \code{extra_dim = list(plev_name = "pres", plev_value = 850)}.
#' 
#' @examples
#' \dontrun{
#' # Computing seasonal (3-months) average SST from 1981 to 2024
#' 
#' in_nc_path_dir = '/home/data/ERSSTv5'
#' in_nc_file_fromat = 'ersst.v5.%s%s.nc'
#' out_nc_path_dir = '/home/data/ERSSTv5_seasonal'
#' 
#' aggregate_netcdf_data('monthly', 'seasonal',
#'                       '1981-01', '2024-12',
#'                       in_nc_path_dir, in_nc_file_fromat,
#'                       out_nc_path_dir,
#'                       aggregate_fun = mean, na.rm = TRUE,
#'                       season_length = 3,
#'                       min_frac = 0.95)
#' 
#' # Computing the annual number of rainy days from 1991 to 2020
#' 
#' in_nc_path_dir = '/home/data/CHIRPSv2_daily' 
#' in_nc_file_fromat = 'chirps_%s%s%s.nc'
#' out_nc_path_dir = '/home/data/CHIRPSv2_rainy_days'
#' out_nc_file_prefix = 'rain_days'
#' out_nc_variable = list(name = 'rainy_days', units = 'days',
#'                        longname = 'Annual number of rainy days',
#'                        missval = -99, prec = 'short')
#' 
#' number_of_rainy_days <- function(x, threshold){
#'     nb_days <- sum(x >= threshold, na.rm = TRUE)
#'     return(nb_days)
#' }
#' 
#' aggregate_netcdf_data('daily', 'annual',
#'                       '1991-01-01', '2020-12-31',
#'                       in_nc_path_dir, in_nc_file_fromat,
#'                       out_nc_path_dir, out_nc_file_prefix,
#'                       aggregate_fun = number_of_rainy_days, threshold = 1.,
#'                       out_nc_variable = out_nc_variable)
#' 
#' # Computing seasonal (3-months) average for u-wind at the pressure level 850 hPa from 1981 to 2024
#' in_nc_path_dir = '/home/data/MONTHLY/UGRD'
#' in_nc_file_fromat = 'u_%s%s.nc'
#' out_nc_path_dir = '/home/data/UGRD_SEASONAL'
#' 
#' aggregate_netcdf_data('monthly', 'seasonal',
#'                       '1981-01', '2024-12',
#'                       in_nc_path_dir, in_nc_file_fromat,
#'                       out_nc_path_dir,
#'                       aggregate_fun = mean, na.rm = TRUE,
#'                       season_length = 3,
#'                       min_frac = 0.95,
#'                       extra_dim = list(plev_name = 'P', plev_value = 850))
#' }
#' 
#' @export

aggregate_netcdf_data <- function(in_time_res, out_time_res,
                                  start_date, end_date,
                                  in_nc_path_dir, in_nc_file_fromat,
                                  out_nc_path_dir, out_nc_file_prefix = NULL,
                                  out_nc_variable = NULL,
                                  aggregate_fun = sum, ...,
                                  season_length = 3, min_frac = 0.95,
                                  var_name = NA, lon_name = NA, lat_name = NA,
                                  extra_dim = NULL)
{
    if(!in_time_res %in% c('daily', 'dekadal', 'monthly'))
        stop('Unknown input time resolution.')
    if(!out_time_res %in% c('monthly', 'seasonal', 'annual'))
        stop('Unknown output time resolution.')
    if(in_time_res == 'monthly' && out_time_res == 'monthly')
        stop('Invalid input and output time resolutions.')

    pattern <- format_pattern(in_nc_file_fromat)
    pattern <- gsub('%s', '.+', pattern)
    ncfiles <- list.files(in_nc_path_dir, pattern = pattern)
    if(length(ncfiles) == 0) stop('No netCDF files found.')

    dates <- extract_filename_dates(ncfiles, in_nc_file_fromat)
    if(in_time_res == 'daily'){
        dates <- substr(dates, 1, 8)
    }else if(in_time_res == 'dekadal'){
        dates <- substr(dates, 1, 7)
    }else{
        dates <- paste0(substr(dates, 1, 6), '15')
        start_date <- paste0(start_date, '-01')
        tmp_date <- as.Date(paste0(end_date, '-01'))
        nd <- nbdays_of_months(tmp_date)
        end_date <- paste0(end_date, '-', nd)
    }
    dates <- as.Date(dates, '%Y%m%d')
    ina <- is.na(dates)
    if(all(ina))
        stop('Unknown netCDF file names date format.')
    dates <- dates[!ina]
    ncfiles <- ncfiles[!ina]

    start_date <- as.Date(start_date, '%Y-%m-%d')
    end_date <- as.Date(end_date, '%Y-%m-%d')
    if(out_time_res == 'seasonal'){
        end_date_s <- end_date
        end_date <- dates_add_months(end_date, season_length - 1)
    }
    it <- dates >= start_date & dates <= end_date
    dates <- dates[it]
    ncfiles <- ncfiles[it]

    if(out_time_res == 'monthly'){
        index <- split(seq_along(dates), format(dates, '%Y%m'))
        odates <- names(index)
        nc_dates <- as.Date(paste0(odates, '15'), '%Y%m%d')
        if(in_time_res == 'daily'){
            nl <- nbdays_of_months(nc_dates)
        }else{
            nl <- rep(3, length(odates))
        }
    }else if(out_time_res == 'annual'){
        index <- split(seq_along(dates), format(dates, '%Y'))
        odates <- names(index)
        nc_dates <- as.Date(paste0(odates, '0101'), '%Y%m%d')
        if(in_time_res == 'daily'){
            nl <- nbdays_of_years(nc_dates)
        }else if(in_time_res == 'dekadal'){
            nl <- rep(36, length(odates))
        }else{
            nl <- rep(12, length(odates))
        }
    }else{
        seasons <- season_define_dates(start_date, end_date_s, season_length)
        index <- lapply(seq_along(seasons$start), function(j){
            ix <- dates >= seasons$start[j] & dates <= seasons$end[j]
            if(!any(ix)) return(NULL)
            which(ix)
        })
        inul <- sapply(index, is.null)
        index <- index[!inul]
        seasons <- seasons[!inul, ]
        seas1 <- format(seasons$start, '%Y-%m')
        seas2 <- format(seasons$end, '%Y-%m')
        odates <- paste0(seas1, '_', seas2)
        nc_dates <- season_middle_dates(odates)
        if(in_time_res == 'daily'){
            nl <- difftime(seasons$end, seasons$start, units = 'days')
            nl <- as.numeric(nl) + 1
        }else if(in_time_res == 'dekadal'){
            nl <- rep(3 * season_length, length(odates))
        }else{
            nl <- rep(season_length, length(odates))
        }
    }

    nba <- sapply(index, length)
    ix <- nba/nl >= min_frac
    if(!any(ix))
        stop('Not enough data to aggregate.')
    index <- index[ix]
    odates <- odates[ix]
    nc_dates <- nc_dates[ix]
    nl <- nl[ix]

    nc <- ncdf4::nc_open(file.path(in_nc_path_dir, ncfiles[1]))
    info <- get_dim_infos_aggregate(nc, var_name, lon_name, lat_name, extra_dim)
    ncvar <- list(name = info$varid,
                  units = nc$var[[info$varid]]$units,
                  longname = nc$var[[info$varid]]$longname,
                  missval = -9999,
                  prec = nc$var[[info$varid]]$prec,
                  compression = 6)
    ncdf4::nc_close(nc)

    dx <- ncdf4::ncdim_def('lon', 'degree_east', info$lon, longname = 'longitude')
    dy <- ncdf4::ncdim_def('lat', 'degree_north', info$lat, longname = 'latitude')

    if(!is.null(out_nc_variable)){
        if(!is.list(out_nc_variable))
            stop('"out_nc_variable" must be a named list.')
        ncvar_args <- methods::formalArgs(ncdf4::ncvar_def)
        ncvar_in <- out_nc_variable[names(out_nc_variable) %in% ncvar_args]
        ncvar <- init_default_list_args(ncvar_in, ncvar)
    }

    fun_args <- list(...)
    in_args <- methods::formalArgs(aggregate_fun)
    if(!is.null(in_args)){
        if('na.rm' %in% in_args){
            fun_args$na.rm <- TRUE
        }
    }else{
        fun_args$na.rm <- TRUE
    }
    fun_args$FUN <- aggregate_fun
    fun_args$MARGIN <- c(1, 2)

    if(is.null(out_nc_file_prefix)){
        nc_prfx <- paste0(info$varid, '_')
    }else{
        nc_prfx <- paste0(out_nc_file_prefix, '_')
    }
    if(out_time_res == 'seasonal'){
        seas_len <- paste0('-', season_length)
    }else{
        seas_len <- ''
    }
    out_dir <- paste0('aggregated_nc_data_', out_time_res, seas_len)
    out_dir <- file.path(out_nc_path_dir, out_dir)
    if(!dir.exists(out_dir))
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    for(j in seq_along(index)){
        t_units <- 'days since 1970-01-01'
        dt <- ncdf4::ncdim_def('time', t_units, as.numeric(nc_dates[j]))
        ncvar$dim <- list(dx, dy, dt)
        nc_grd <- do.call(ncdf4::ncvar_def, ncvar)

        ncdat <- lapply(index[[j]], function(i){
            nc <- ncdf4::nc_open(file.path(in_nc_path_dir, ncfiles[i]))
            z <- ncdf4::ncvar_get(nc, varid = info$varid, collapse_degen = FALSE)
            ncdf4::nc_close(nc)
            extract_nc_data(z, info)
        })

        ncdat <- simplify2array(ncdat)
        nna <- apply(!is.na(ncdat), c(1, 2), sum)
        ina <- nna/nl[j] < min_frac
        nc_args <- c(list(X = ncdat), fun_args)
        ncdat <- do.call(apply, nc_args)
        ncdat[ina] <- ncvar$missval
        ncdat[is.na(ncdat) | is.nan(ncdat)] <- ncvar$missval
        dim(ncdat) <- c(info$nx, info$ny, 1)

        out_ncfile <- paste0(nc_prfx, odates[j], '.nc')
        out_ncfile <- file.path(out_dir, out_ncfile)
        nc <- ncdf4::nc_create(out_ncfile, nc_grd)
        ncdf4::ncvar_put(nc, nc_grd, ncdat)
        ncdf4::nc_close(nc)
    }

    return(0)
}
