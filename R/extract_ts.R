#' Extracting time series.
#'
#' Extract time series for a specific month or season over a defined region.
#' 
#' @param data object of class \code{scda_data}.
#' @param output_type character, the type of output time series to be extracted.
#' The available options are \code{"raw"} (default), \code{"anomaly"} and \code{"pcs"}.
#' \itemize{
#' \item \code{"raw"}: return a time series of the raw data spatially averaged over the region defined from \code{data}.
#' \item \code{"anomaly"}: return a time series of the anomalies computed from the data spatially averaged over the region defined from \code{data}.
#' \item \code{"pcs"}: return a time series (PCs) over the region defined from \code{data}.
#' }
#' @param anomaly_type character, if \code{output_type = "anonamly"}, it specifies the type of anomalies to be calculated 
#' and must be one of \code{"difference"} and \code{"standardized"}.
#' @param base_period an integer vector of length 2, if \code{output_type = "anonamly"}, giving the start and end year of the base period 
#' to compute the climatology to calculate the anomalies. 
#' @param eof_number integer, if \code{output_type = "pcs"}, number of \emph{Empirical Orthogonal Functions} (EOF) requested.
#' @param nth_pcs integer, if \code{output_type = "pcs"}, the \code{nth} PCs to be extracted.
#' @param ... other arguments that may be passed on to the function \code{\link[RSpectra]{eigs_sym}} of the package \pkg{RSpectra}, in case of \code{output_type = "pcs"}.
#' 
#' @return This returns a \code{data.frame} with 2 columns, with header names: \code{date} of class \code{Date} and \code{ts} the time series.
#' 
#' @examples
#' \dontrun{
#' # Extracting the first PCs of SST for July for the period 1981-2020 
#' # averaged over the region longitude: [170 W - 120 W] and latitude [5 S - 5 N]
#' 
#' path_dir_data = '/home/data/ERSSTv5'
#' nc_file_fromat = 'ersst.v5.%s%s.nc'
#' years_range = c(1981, 2020)
#' region_bbox = c(-170, -5, -120, 5)
#' 
#' grd_data = get_month_region_netcdf_data(7, path_dir_data, nc_file_fromat,
#'                                         years_range, region_bbox)
#' ts = extract_region_timeSeries(grd_data, output_type = 'pcs', nth_pcs = 1)
#' 
#' }
#' 
#' @export

extract_region_timeSeries <- function(data, output_type = c('raw', 'anomaly', 'pcs'),
                                      anomaly_type = c('difference', 'standardized'),
                                      base_period = c(1991, 2020),
                                      eof_number = 4, nth_pcs = 1, ...)
{
    if(!inherits(data, 'scda_data'))
        stop("'data' is not an object of class 'scda_data'.")

    output_type <- match.arg(output_type)
    anomaly_type <- match.arg(anomaly_type)
    if(output_type %in% c('raw', 'anomaly')){
        val <- lapply(data$data, mean, na.rm = TRUE)
        val <- do.call(c, val)
        val[is.nan(val)] <- NA

        if(output_type == 'anomaly'){
            years <- as.numeric(format(data$date, '%Y'))
            it <- years >= base_period[1] &
                  years <= base_period[2]
            vmn <- mean(val[it], na.rm = TRUE)
            val <- val - vmn
            if(anomaly_type == 'standardized'){
                vsd <- stats::sd(val[it], na.rm = TRUE)
                val <- val/vsd
            }
        }
    }else{
        eigen_data <- eigen_analysis(data, eof_number, ...)
        val <- get_PCs_timeSeries(eigen_data, nth_pcs)
        val <- val$pc
    }

    data.frame(date = data$date, ts = val)
}
