#' Computing climatology.
#'
#' Compute climatology for a specific month or season.
#' 
#' @param data object of class \code{scda_data}.
#' @param base_period an integer vector of length 2, giving the start and end year of the base period to compute the climatology.
#' 
#' @return This returns an object of class \code{scda_climato}. The object has the following elements:
#' \itemize{
#' \item \code{$lon}: a numeric vector of the longitude with length \code{N}.
#' \item \code{$lat}: a numeric vector of the latitude with length \code{M}.
#' \item \code{$base_period}: the base period used to compute the climatology.
#' \item \code{$data}: a list of 2 matrices of the \strong{mean} and \strong{standard deviation}, with dimension \code{N x M}.
#' }
#' And other information from the input data.
#' 
#' @export

compute_climatology <- function(data, base_period = c(1991, 2020)){
    if(!inherits(data, 'scda_data'))
        stop("'data' is not an object of class 'scda_data'.")

    years <- as.numeric(format(data$date, '%Y'))
    it <- years >= base_period[1] & years <= base_period[2]

    tmp <- simplify2array(data$data)
    mtmp <- apply(tmp, c(1, 2), mean, na.rm = TRUE)
    stmp <- apply(tmp, c(1, 2), stats::sd, na.rm = TRUE)
    dims <- sapply(data[c('lon', 'lat')], length)
    dim(mtmp) <- dims
    dim(stmp) <- dims

    out <- unclass(data)
    names(out)[names(out) == 'date'] <- 'base_period'
    out[['base_period']] <- base_period
    out[['data']] <- list('mean' = mtmp, 'sd' = stmp)

    assign_class(out, 'scda_climato')
}

#' Computing anomalies.
#'
#' Compute anomalies for a specific month or season.
#' 
#' @param data object of class \code{scda_data}.
#' @param clim_data object of class \code{scda_climato}, the climatology data computed from the function \code{\link[simpleCDA]{compute_climatology}}.
#' @param years_range an integer vector of length 2, the start and end year of the period the anomalies will be computed.
#' Default \code{NULL}, using the whole period available in the dataset.
#' @param anomaly_type character, it specifies the type of anomalies to be calculated 
#' and must be one of \code{"difference"} and \code{"standardized"}.
#' 
#' @return This returns an object of class \code{scda_anomalies} and \code{scda_data}. The object has the following elements:
#' \itemize{
#' \item \code{$lon}: a numeric vector of the longitude with length \code{N}.
#' \item \code{$lat}: a numeric vector of the latitude with length \code{M}.
#' \item \code{$date}: a vector of class \code{Date} with length \code{L}.
#' \item \code{$data}: a list of matrices with dimension \code{N x M}, the length of the list is \code{L}.
#' }
#' And other information from the input data.
#' 
#' @examples
#' \dontrun{
#' # Standardized anomaly of the SST for the season DJF for the period 1981-2020.
#' 
#' # read the data
#' path_dir_data = '/home/data/ERSSTv5_seasonal'
#' nc_file_fromat = 'ersstv5_%s-%s_%s-%s.nc'
#' grd_data = get_season_region_netcdf_data('12-02', path_dir_data, nc_file_fromat,
#'                                          years_range = c(1981, 2020))
#' 
#' # compute climatology for the period 1991 to 2020
#' clim_data = compute_climatology(grd_data, base_period = c(1991, 2020))
#' 
#' # compute the standardized anomaly
#' anom_data = compute_anomalies(grd_data, clim_data, 
#'                               years_range = c(1981, 2020),
#'                               anomaly_type = 'standardized')
#' }
#' 
#' @export

compute_anomalies <- function(data, clim_data, years_range = NULL,
                              anomaly_type = c('difference', 'standardized'))
{
    if(!inherits(data, 'scda_data'))
        stop("'data' is not an object of class 'scda_data'.")

    if(!inherits(clim_data, 'scda_climato'))
        stop("'clim_data' is not an object of class 'scda_climato'.")

    ddims <- sapply(data[c('lon', 'lat')], length)
    cdims <- sapply(clim_data[c('lon', 'lat')], length)
    if(!isTRUE(all.equal(ddims, cdims)))
        stop("'data' and 'clim_data' are not the same dimensions.")

    anomaly_type <- match.arg(anomaly_type)

    if(!is.null(years_range)){
        years <- as.numeric(format(data$date, '%Y'))
        iy <- years >= years_range[1] & years <= years_range[2]
    }else{
        iy <- rep(TRUE, length(data$date))
    }
    data$date <- data$date[iy]
    data$data <- data$data[iy]
    data$data <- lapply(data$data, function(x){
        z <- x - clim_data$data$mean
        if(anomaly_type == 'standardized'){
            z <- z / clim_data$data$sd
        }
        z
    })

    assign_class(data, 'scda_anomalies')
}
