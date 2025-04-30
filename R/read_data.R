#' Read monthly netCDF dataset.
#'
#' Read monthly netCDF data for a specific month over a defined region.
#' 
#' @param month integer, the month to be read, must be from 1 to 12.
#' @param path_dir_data character, full path to the folder containing the netCDF files.
#' @param nc_file_fromat character, format of the netCDF file name, the year and month must be replaced by \code{\%s}.
#' Example: for \emph{ersst.v5.202401.nc} the file name format should be \emph{ersst.v5.\%s\%s.nc}.
#' @param years_range an integer vector of length 2, the start and end year of the data to be extracted.
#' Default \code{NULL}, using the whole period available in the dataset.
#' @param region_bbox a numeric vector of length 4, the bounding box of the region to be extracted,
#' with longitude min, latitude min, longitude max and latitude max of the region (West, South, East, North).
#' Default \code{NULL}, no extraction performed, using the entire region from the dataset.
#' @param var_name character, name of the variable to be read from the netCDF data. If unspecified or left \code{NA},
#' the first variable available in the netCDF data will be taken.
#' @param lon_name character, the name of the longitude dimension. If unspecified or left \code{NA}, it will be detected. 
#' @param lat_name character, the name of the latitude dimension. Same as \code{lon_name}.
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
#' # over the region longitude: [170 W - 120 W] and latitude [5 S - 5 N]
#' 
#' path_dir_data = '/home/data/ERSSTv5'
#' nc_file_fromat = 'ersst.v5.%s%s.nc'
#' years_range = c(1981, 2020)
#' region_bbox = c(-170, -5, -120, 5)
#' 
#' grd_data = get_month_region_netcdf_data(7, path_dir_data, nc_file_fromat,
#'                                         years_range, region_bbox)
#' }
#' 
#' @export

get_month_region_netcdf_data <- function(month, path_dir_data, nc_file_fromat,
                                         years_range = NULL, region_bbox = NULL,
                                         var_name = NA, lon_name = NA,
                                         lat_name = NA)
{
    month <- as.numeric(month)
    if(month < 1 || month > 12) stop('Invalid month.')

    pattern <- format_pattern(nc_file_fromat)
    pattern <- gsub('%s', '.+', pattern)
    ncfiles <- list.files(path_dir_data, pattern = pattern)
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

    nc <- ncdf4::nc_open(file.path(path_dir_data, ncfiles[1]))
    if(is.na(lon_name) || is.na(lat_name)){
        xy_name <- get_latlon_dimnames(nc)
        lon_name <- xy_name$lon
        lat_name <- xy_name$lat
    }

    lon <- nc$dim[[lon_name]]$vals
    lon <- format_longitude(lon)
    lat <- nc$dim[[lat_name]]$vals

    if(is.na(var_name)){
        var_name <- nc$var[[1]]$name
    }
    xy_order <- get_latlon_order(nc, var_name, lon_name, lat_name)
    ncdf4::nc_close(nc)

    if(!is.null(region_bbox)){
        region_bbox[1] <- format_longitude(region_bbox[1])
        region_bbox[3] <- format_longitude(region_bbox[3])
        ix <- lon >= region_bbox[1] & lon <= region_bbox[3]
        iy <- lat >= region_bbox[2] & lat <= region_bbox[4]
        lon <- lon[ix]
        lat <- lat[iy]
    }else{
        ix <- rep(TRUE, length(lon))
        iy <- rep(TRUE, length(lat))
        rx <- range(lon)
        ry <- range(lat)
        region_bbox <- c(rx[1], ry[1], rx[2], ry[2])
    }

    ox <- order(lon)
    oy <- order(lat)
    lon <- lon[ox]
    lat <- lat[oy]

    don <- lapply(seq_along(ncfiles), function(i){
        nc <- ncdf4::nc_open(file.path(path_dir_data, ncfiles[i]))
        z <- ncdf4::ncvar_get(nc, varid = var_name)
        ncdf4::nc_close(nc)
        if(xy_order$ilon < xy_order$ilat){
            z <- z[ix, iy, drop = FALSE]
            z <- z[ox, oy, drop = FALSE]
        }else{
            z <- z[iy, ix, drop = FALSE]
            z <- z[oy, ox, drop = FALSE]
            z <- t(z)
        }
        z
    })
    names(don) <- format(dates, '%Y-%m')

    out <- list(lon = lon, lat = lat, date = dates, data = don,
                time_res = 'monthly', month = month, bbox = region_bbox)

    assign_class(out, 'scda_data')
}

#' Read seasonal netCDF dataset.
#'
#' Read seasonal netCDF data for a specific season over a defined region.
#' 
#' @param season character, the season to be read, in the format \code{"start_month-end_month"}, e.g. "12-02" for DJF.
#' @param path_dir_data character, full path to the folder containing the netCDF files.
#' @param nc_file_fromat character, format of the netCDF file name, the year and month must be replaced by \code{\%s}.
#' Example: for \emph{ersstv5_2024-01_2024-03.nc} the file name format should be \emph{ersstv5_\%s-\%s_\%s-\%s.nc}.
#' @param years_range an integer vector of length 2, the start and end year of the data to be extracted.
#' Default \code{NULL}, using the whole period available in the dataset.
#' @param region_bbox a numeric vector of length 4, the bounding box of the region to be extracted,
#' with longitude min, latitude min, longitude max and latitude max of the region (West, South, East, North).
#' Default \code{NULL}, no extraction performed, using the entire region from the dataset.
#' @param var_name character, name of the variable to be read from the netCDF data. If unspecified or left \code{NA},
#' the first variable available in the netCDF data will be taken.
#' @param lon_name character, the name of the longitude dimension. If unspecified or left \code{NA}, it will be detected. 
#' @param lat_name character, the name of the latitude dimension. Same as \code{lon_name}.
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
#' # over the region longitude: [170 W - 120 W] and latitude [5 S - 5 N]
#' 
#' path_dir_data = '/home/data/ERSSTv5_seasonal'
#' nc_file_fromat = 'ersstv5_%s-%s_%s-%s.nc'
#' years_range = c(1981, 2020)
#' region_bbox = c(-170, -5, -120, 5)
#' 
#' grd_data = get_season_region_netcdf_data("12-02", path_dir_data, nc_file_fromat,
#'                                         years_range, region_bbox)
#' }
#' 
#' @export

get_season_region_netcdf_data <- function(season, path_dir_data, nc_file_fromat,
                                          years_range = NULL, region_bbox = NULL,
                                          var_name = NA, lon_name = NA,
                                          lat_name = NA)
{
    season <- as.numeric(strsplit(season, '-')[[1]])
    if(any(season < 1 | season > 12)) stop('Invalid season.')

    pattern <- format_pattern(nc_file_fromat)
    pattern <- gsub('%s', '.+', pattern)
    ncfiles <- list.files(path_dir_data, pattern = pattern)
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

    nc <- ncdf4::nc_open(file.path(path_dir_data, ncfiles[1]))
    if(is.na(lon_name) || is.na(lat_name)){
        xy_name <- get_latlon_dimnames(nc)
        lon_name <- xy_name$lon
        lat_name <- xy_name$lat
    }

    lon <- nc$dim[[lon_name]]$vals
    lon <- format_longitude(lon)
    lat <- nc$dim[[lat_name]]$vals

    if(is.na(var_name)){
        var_name <- nc$var[[1]]$name
    }
    xy_order <- get_latlon_order(nc, var_name, lon_name, lat_name)
    ncdf4::nc_close(nc)

    if(!is.null(region_bbox)){
        region_bbox[1] <- format_longitude(region_bbox[1])
        region_bbox[3] <- format_longitude(region_bbox[3])
        ix <- lon >= region_bbox[1] & lon <= region_bbox[3]
        iy <- lat >= region_bbox[2] & lat <= region_bbox[4]
        lon <- lon[ix]
        lat <- lat[iy]
    }else{
        ix <- rep(TRUE, length(lon))
        iy <- rep(TRUE, length(lat))
        rx <- range(lon)
        ry <- range(lat)
        region_bbox <- c(rx[1], ry[1], rx[2], ry[2])
    }

    ox <- order(lon)
    oy <- order(lat)
    lon <- lon[ox]
    lat <- lat[oy]

    don <- lapply(seq_along(ncfiles), function(i){
        nc <- ncdf4::nc_open(file.path(path_dir_data, ncfiles[i]))
        z <- ncdf4::ncvar_get(nc, varid = var_name)
        ncdf4::nc_close(nc)
        if(xy_order$ilon < xy_order$ilat){
            z <- z[ix, iy, drop = FALSE]
            z <- z[ox, oy, drop = FALSE]
        }else{
            z <- z[iy, ix, drop = FALSE]
            z <- z[oy, ox, drop = FALSE]
            z <- t(z)
        }
        z
    })
    names(don) <- lname

    start <- season[1]
    len <- (season[2] - season[1]) %% 12 + 1
    season <- sprintf('%02d-%02d', season[1], season[2])
    out <- list(lon = lon, lat = lat, date = dates, data = don,
                time_res = 'seasonal', season = season,
                season_start = start, season_length = len,
                bbox = region_bbox)

    assign_class(out, 'scda_data')
}
