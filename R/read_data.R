#' Read monthly netCDF dataset.
#'
#' Read monthly netCDF data for a specific month over a defined region.
#' 
#' @param month integer, the month to be read, must be from 1 to 12.
#' @param nc_path_dir character, full path to the folder containing the netCDF files.
#' @param nc_file_fromat character, format of the netCDF file name, the year and month must be replaced by \code{\%s}.
#' Example: for \emph{ersst.v5.202401.nc} the file name format should be \emph{'ersst.v5.\%s\%s.nc'}.
#' @param years_range an integer vector of length 2, the start and end year of the data to be extracted.
#' Default \code{NULL}, get the whole period available in the dataset.
#' @param region_bbox a numeric vector of length 4, the bounding box of the region to be extracted,
#' with longitude min, latitude min, longitude max and latitude max of the region (West, South, East, North).
#' Default \code{NULL}, no extraction performed, get the entire region from the dataset.
#' @param var_name character, name of the variable to be read from the netCDF data. If unspecified or left \code{NA},
#' the first variable available in the netCDF data will be taken.
#' @param lon_name character, the name of the longitude dimension. If unspecified or left \code{NA}, it will be detected. 
#' @param lat_name character, the name of the latitude dimension. Same as \code{lon_name}.
#' @param ... pairs of arguments, providing the names and values of any extra dimensions to be extracted.
#' The arguments should be in the format \code{<information about the dimension>_name} for the name of the dimension
#' and \code{<information about the dimension>_value} for the value to be extracted from that dimension.
#' Example: if your netCDF data has a pressure levels dimension named \code{"pres"} having the values from 1000 to 10 hPa,
#' if you want to extract the data at the pressure level 850 hPa, 
#' the pair of arguments would be \code{plev_name = "pres"} and \code{plev_value = 850}.
#' 
#' @return This returns an object of class \code{scda_data}. The object has the following elements:
#' \itemize{
#' \item \code{$lon}: a numeric vector of the longitude with length \code{N}.
#' \item \code{$lat}: a numeric vector of the latitude with length \code{M}.
#' \item \code{$date}: a vector of class \code{Date} with length \code{L}.
#' \item \code{$data}: a list of matrices with dimension \code{N x M}, the length of the list is \code{L}.
#' }
#' 
#' @examples
#' \dontrun{
#' # Extracting the SST data for July for the period 1981-2020 
#' # over the region: longitude [170 W - 120 W] and latitude [5 S - 5 N]
#' 
#' nc_path_dir = '/home/data/ERSSTv5'
#' nc_file_fromat = 'ersst.v5.%s%s.nc'
#' years_range = c(1981, 2020)
#' region_bbox = c(-170, -5, -120, 5)
#' 
#' grd_data = get_month_region_netcdf_data(7, nc_path_dir, nc_file_fromat,
#'                                         years_range, region_bbox)
#' 
#' # Extracting the U-wind data for July for the period 1981-2020 at 850 hPa
#' # over the region: longitude [170 W - 120 W] and latitude [5 S - 5 N]
#' 
#' nc_path_dir = '/home/data/MONTHLY/UGRD'
#' nc_file_fromat = 'u_%s%s.nc'
#' years_range = c(1981, 2020)
#' region_bbox = c(-170, -5, -120, 5)
#' 
#' ugrd_data = get_month_region_netcdf_data(7, nc_path_dir, nc_file_fromat,
#'                                          years_range, region_bbox,
#'                                          alt_name = 'P', alt_value = 850)
#' }
#' 
#' @export

get_month_region_netcdf_data <- function(month, nc_path_dir, nc_file_fromat,
                                         years_range = NULL, region_bbox = NULL,
                                         var_name = NA, lon_name = NA, lat_name = NA,
                                         ...)
{
    month <- as.numeric(month)
    if(month < 1 || month > 12) stop('Invalid month.')

    pattern <- format_pattern(nc_file_fromat)
    pattern <- gsub('%s', '.+', pattern)
    ncfiles <- list.files(nc_path_dir, pattern = pattern)
    if(length(ncfiles) == 0) stop('No netCDF files found.')

    dates <- extract_filename_dates(ncfiles, nc_file_fromat)
    dates <- as.Date(paste0(substr(dates, 1, 6), '15'), '%Y%m%d')
    months <- format(dates, '%m')
    years <- format(dates, '%Y')
    nmonth <- as.numeric(months)

    if(!is.null(years_range)){
        nyear <- as.numeric(years)
        imy <- nmonth == month & nyear >= years_range[1] & nyear <= years_range[2]
    }else{
        imy <- nmonth == month
    }
    ncfiles <- sprintf(nc_file_fromat, years[imy], months[imy])
    dates <- dates[imy]

    nc <- ncdf4::nc_open(file.path(nc_path_dir, ncfiles[1]))
    info <- get_dim_infos_extract(nc, var_name, lon_name, lat_name, region_bbox, ...)
    ncdf4::nc_close(nc)

    don <- lapply(seq_along(ncfiles), function(i){
        nc <- ncdf4::nc_open(file.path(nc_path_dir, ncfiles[i]))
        z <- ncdf4::ncvar_get(nc, varid = info$varid, collapse_degen = FALSE)
        ncdf4::nc_close(nc)
        extract_nc_data(z, info)
    })
    names(don) <- format(dates, '%Y-%m')

    out <- list(lon = info$lon, lat = info$lat, date = dates, data = don,
                time_res = 'monthly', month = month, bbox = info$bbox)
    assign_class(out, 'scda_data')
}

#' Read seasonal netCDF dataset.
#'
#' Read seasonal netCDF data for a specific season over a defined region.
#' 
#' @param season character, the season to be read, in the format \code{"start_month-end_month"}, e.g. "12-02" for DJF.
#' @param nc_path_dir character, full path to the folder containing the netCDF files.
#' @param nc_file_fromat character, format of the netCDF file name, the year and month must be replaced by \code{\%s}.
#' Example: for \emph{ersstv5_2024-01_2024-03.nc} the file name format should be \emph{'ersstv5_\%s-\%s_\%s-\%s.nc'}.
#' @param years_range an integer vector of length 2, the start and end year of the data to be extracted.
#' Default \code{NULL}, get the whole period available in the dataset.
#' @param region_bbox a numeric vector of length 4, the bounding box of the region to be extracted,
#' with longitude min, latitude min, longitude max and latitude max of the region (West, South, East, North).
#' Default \code{NULL}, no extraction performed, get the entire region from the dataset.
#' @param var_name character, name of the variable to be read from the netCDF data. If unspecified or left \code{NA},
#' the first variable available in the netCDF data will be taken.
#' @param lon_name character, the name of the longitude dimension. If unspecified or left \code{NA}, it will be detected. 
#' @param lat_name character, the name of the latitude dimension. Same as \code{lon_name}.
#' @param ... pairs of arguments, providing the names and values of any extra dimensions to be extracted.
#' The arguments should be in the format \code{<information about the dimension>_name} for the name of the dimension
#' and \code{<information about the dimension>_value} for the value to be extracted from that dimension.
#' Example: if your netCDF data has a pressure levels dimension named \code{"pres"} having the values from 1000 to 10 hPa,
#' if you want to extract the data at the pressure level 850 hPa, 
#' the pair of arguments would be \code{plev_name = "pres"} and \code{plev_value = 850}.
#' 
#' @return This returns an object of class \code{scda_data}. The object has the following elements:
#' \itemize{
#' \item \code{$lon}: a numeric vector of the longitude with length \code{N}.
#' \item \code{$lat}: a numeric vector of the latitude with length \code{M}.
#' \item \code{$date}: a vector of class \code{Date} with length \code{L}.
#' \item \code{$data}: a list of matrices with dimension \code{N x M}, the length of the list is \code{L}.
#' }
#' 
#' @examples
#' \dontrun{
#' # Extracting the SST data for the season DJF for the period 1981-2020 
#' # over the region: longitude [170 W - 120 W] and latitude [5 S - 5 N]
#' 
#' nc_path_dir = '/home/data/ERSSTv5_seasonal'
#' nc_file_fromat = 'ersstv5_%s-%s_%s-%s.nc'
#' years_range = c(1981, 2020)
#' region_bbox = c(-170, -5, -120, 5)
#' 
#' grd_data = get_season_region_netcdf_data("12-02", nc_path_dir, nc_file_fromat,
#'                                         years_range, region_bbox)
#' }
#' 
#' @export

get_season_region_netcdf_data <- function(season, nc_path_dir, nc_file_fromat,
                                          years_range = NULL, region_bbox = NULL,
                                          var_name = NA, lon_name = NA, lat_name = NA,
                                          ...)
{
    season <- as.numeric(strsplit(season, '-')[[1]])
    if(any(season < 1 | season > 12)) stop('Invalid season.')

    pattern <- format_pattern(nc_file_fromat)
    pattern <- gsub('%s', '.+', pattern)
    ncfiles <- list.files(nc_path_dir, pattern = pattern)
    if(length(ncfiles) == 0) stop('No netCDF files found.')

    dates <- extract_filename_dates(ncfiles, nc_file_fromat)
    syear <- as.numeric(substr(dates, 1, 4))
    smonth <- as.numeric(substr(dates, 5, 6))
    eyear <- as.numeric(substr(dates, 7, 10))
    emonth <- as.numeric(substr(dates, 11, 12))

    seasons <- paste0(smonth, '_', emonth)
    seas <- paste0(season[1], '_', season[2])
    iseas <- seasons == seas

    if(!is.null(years_range)){
        isy <- iseas & syear >= years_range[1] & eyear <= years_range[2]
    }else{
        isy <- iseas
    }
    ncfiles <- ncfiles[isy]
    dates <- dates[isy]
    lname <- season_filename_dates(dates)
    dates <- season_middle_dates(dates)

    nc <- ncdf4::nc_open(file.path(nc_path_dir, ncfiles[1]))
    info <- get_dim_infos_extract(nc, var_name, lon_name, lat_name, region_bbox, ...)
    ncdf4::nc_close(nc)

    don <- lapply(seq_along(ncfiles), function(i){
        nc <- ncdf4::nc_open(file.path(nc_path_dir, ncfiles[i]))
        z <- ncdf4::ncvar_get(nc, varid = info$varid, collapse_degen = FALSE)
        ncdf4::nc_close(nc)
        extract_nc_data(z, info)
    })
    names(don) <- lname

    start <- season[1]
    len <- (season[2] - season[1]) %% 12 + 1
    season <- sprintf('%02d-%02d', season[1], season[2])
    out <- list(lon = info$lon, lat = info$lat, date = dates, data = don,
                time_res = 'seasonal', season = season,
                season_start = start, season_length = len, bbox = info$bbox)
    assign_class(out, 'scda_data')
}

#' Read annual netCDF dataset.
#'
#' Read annual netCDF data over a defined region.
#' 
#' @param nc_path_dir character, full path to the folder containing the netCDF files.
#' @param nc_file_fromat character, format of the netCDF file name, the year must be replaced by \code{\%s}.
#' Example: for \emph{chirps_2024.nc} the file name format should be \emph{'chirps_\%s.nc'}.
#' @param years_range an integer vector of length 2, the start and end year of the data to be extracted.
#' Default \code{NULL}, get the whole period available in the dataset.
#' @param region_bbox a numeric vector of length 4, the bounding box of the region to be extracted,
#' with longitude min, latitude min, longitude max and latitude max of the region (West, South, East, North).
#' Default \code{NULL}, no extraction performed, get the entire region from the dataset.
#' @param var_name character, name of the variable to be read from the netCDF data. If unspecified or left \code{NA},
#' the first variable available in the netCDF data will be taken.
#' @param lon_name character, the name of the longitude dimension. If unspecified or left \code{NA}, it will be detected. 
#' @param lat_name character, the name of the latitude dimension. Same as \code{lon_name}.
#' @param ... pairs of arguments, providing the names and values of any extra dimensions to be extracted.
#' The arguments should be in the format \code{<information about the dimension>_name} for the name of the dimension
#' and \code{<information about the dimension>_value} for the value to be extracted from that dimension.
#' Example: if your netCDF data has a pressure levels dimension named \code{"pres"} having the values from 1000 to 10 hPa,
#' if you want to extract the data at the pressure level 850 hPa, 
#' the pair of arguments would be \code{plev_name = "pres"} and \code{plev_value = 850}.
#' 
#' @return This returns an object of class \code{scda_data}. The object has the following elements:
#' \itemize{
#' \item \code{$lon}: a numeric vector of the longitude with length \code{N}.
#' \item \code{$lat}: a numeric vector of the latitude with length \code{M}.
#' \item \code{$date}: a vector of class \code{Date} with length \code{L}.
#' \item \code{$data}: a list of matrices with dimension \code{N x M}, the length of the list is \code{L}.
#' }
#' 
#' @examples
#' \dontrun{
#' # Extracting CHIRPSv2.0 annual rainfall from global data for the period 1991-2020 
#' # over Madagascar: longitude [42 E - 52 E] and latitude [26 S - 11.5 S]
#' 
#' nc_path_dir = '/home/data/CHIRPSv2.0_annual'
#' nc_file_fromat = 'chirps_%s.nc'
#' years_range = c(1991, 2020)
#' region_bbox = c(42, -26, 52, -11.5)
#' 
#' grd_data = get_annual_region_netcdf_data(nc_path_dir, nc_file_fromat,
#'                                         years_range, region_bbox)
#' }
#' 
#' @export

get_annual_region_netcdf_data <- function(nc_path_dir, nc_file_fromat,
                                          years_range = NULL, region_bbox = NULL,
                                          var_name = NA, lon_name = NA, lat_name = NA,
                                          ...)
{
    pattern <- format_pattern(nc_file_fromat)
    pattern <- gsub('%s', '.+', pattern)
    ncfiles <- list.files(nc_path_dir, pattern = pattern)
    if(length(ncfiles) == 0) stop('No netCDF files found.')

    dates <- extract_filename_dates(ncfiles, nc_file_fromat)
    dates <- as.Date(paste0(substr(dates, 1, 4), '0101'), '%Y%m%d')
    years <- format(dates, '%Y')

    if(!is.null(years_range)){
        nyear <- as.numeric(years)
        iy <- nyear >= years_range[1] & nyear <= years_range[2]
    }else{
        iy <- rep(TRUE, length(years))
    }
    ncfiles <- sprintf(nc_file_fromat, years[iy])
    dates <- dates[iy]

    nc <- ncdf4::nc_open(file.path(nc_path_dir, ncfiles[1]))
    info <- get_dim_infos_extract(nc, var_name, lon_name, lat_name, region_bbox, ...)
    ncdf4::nc_close(nc)

    don <- lapply(seq_along(ncfiles), function(i){
        nc <- ncdf4::nc_open(file.path(nc_path_dir, ncfiles[i]))
        z <- ncdf4::ncvar_get(nc, varid = info$varid, collapse_degen = FALSE)
        ncdf4::nc_close(nc)
        extract_nc_data(z, info)
    })
    names(don) <- format(dates, '%Y')

    out <- list(lon = info$lon, lat = info$lat, date = dates, data = don,
                time_res = 'annual', bbox = info$bbox)
    assign_class(out, 'scda_data')
}