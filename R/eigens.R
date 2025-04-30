#' Performing Principal Component Analysis (PCA).
#'
#' Perform eigenanalysis of covariance matrix.
#' 
#' @param data object of class \code{scda_data}.
#' @param eof_number integer, number of \emph{Empirical Orthogonal Functions} (EOF) requested.
#' @param ... other arguments that may be passed on to the function \code{\link[RSpectra]{eigs_sym}} of the package \pkg{RSpectra}.
#' 
#' @return This returns an object of class \code{scda_eigen}.
#' In addition to the output from the function \code{\link[RSpectra]{eigs_sym}} of the package \pkg{RSpectra}.
#' the object has the additional elements:
#' \itemize{
#' \item \code{$lon}: a numeric vector of the longitude with length \code{N} from \code{data}.
#' \item \code{$lat}: a numeric vector of the latitude with length \code{M} from \code{data}.
#' \item \code{$date}: a vector of class \code{Date} with length \code{L} from \code{data}.
#' \item \code{$anomaly}: a matrix of the standardized anomalies with dimension \code{L x (N*M)}.
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
#' 
#' eigen_data = eigen_analysis(grd_data)
#' }
#'
#' @export

eigen_analysis <- function(data, eof_number = 4, ...){
    if(!inherits(data, 'scda_data'))
        stop("'data' is not an object of class 'scda_data'.")

    # standardize all pixels
    tmp <- simplify2array(data$data)
    mtmp <- mean(tmp, na.rm = TRUE)
    stmp <- stats::sd(tmp, na.rm = TRUE)
    anom_grd <- lapply(data$data, function(x) (x - mtmp)/stmp)
    anom_grd <- lapply(anom_grd, c)
    anom_grd <- do.call(rbind, anom_grd)

    # set NA to zero
    tmp <- anom_grd
    tmp[is.na(tmp) | is.nan(tmp)] <- 0

    # calculate covariance matrix
    nc <- ncol(tmp)
    size <- if(nc < 2000) nc else 2000
    cvm <- propagate::bigcor(tmp, fun = 'cov', size = size, verbose = FALSE)

    # eigen analysis of covariance matrix
    eigen_data <- RSpectra::eigs_sym(cvm[], eof_number, ...)
    dims <- c('lon', 'lat', 'date')
    eigen_data[dims] <- data[dims]
    eigen_data$anomaly <- anom_grd

    assign_class(eigen_data, 'scda_eigen')
}

#' Extracting the principal component (PCs).
#'
#' Extract a specific time series (PCs).
#' 
#' @param eigen_data object of class \code{scda_eigen}.
#' @param k integer, the \code{nth} PCs to be extracted.
#' 
#' @return This returns an object of class \code{scda_eigen_pcs}.
#' The object has the following elements:
#' \itemize{
#' \item \code{$date}: a vector of class \code{Date}.
#' \item \code{$pc}: a numeric vector of the \code{nth} time series (PCs).
#' }
#' 
#' @examples
#' \dontrun{
#' ## grd_data: object of class "scda_data"
#' 
#' ## perform eigenanalysis
#' eigen_data = eigen_analysis(grd_data)
#' 
#' ## get the second PCs
#' pc2 = get_PCs_timeSeries(eigen_data, 2)
#' }
#'
#' @export

get_PCs_timeSeries <- function(eigen_data, k){
    if(!inherits(eigen_data, 'scda_eigen'))
        stop("'eigen_data' is not an object of class 'scda_eigen'.")

    ev <- -1 * eigen_data$vectors
    y <- eigen_data$anomaly %*% ev[, k, drop = FALSE]
    y <- (y - mean(y))/stats::sd(y)
    out <- list(date = eigen_data$date, pc = as.numeric(y[, 1]), k = k)
    assign_class(out, 'scda_eigen_pcs')
}

#' Extracting the spatial patterns (EOFs).
#'
#' Extract a specific spatial patterns (EOFs).
#' 
#' @param eigen_data object of class \code{scda_eigen}.
#' @param k integer, the \code{nth} EOFs to be extracted.
#' 
#' @return This returns an object of class \code{scda_eigen_eofs}.
#' The object has the following elements:
#' \itemize{
#' \item \code{$lon}: a numeric vector of the longitude with length \code{N}.
#' \item \code{$lat}: a numeric vector of the latitude with length \code{M}.
#' \item \code{$eof}: a matrix of the \code{nth} spatial patterns (EOFs) with dimension \code{N x M}.
#' }
#' 
#' @examples
#' \dontrun{
#' ## grd_data: object of class "scda_data"
#' 
#' ## perform eigenanalysis
#' eigen_data = eigen_analysis(grd_data)
#' 
#' ## get the first EOFs
#' eof1 = get_EOFs_spatialPatterns(eigen_data, 1)
#' }
#' 
#' @export

get_EOFs_spatialPatterns <- function(eigen_data, k){
    if(!inherits(eigen_data, 'scda_eigen'))
        stop("'eigen_data' is not an object of class 'scda_eigen'.")

    eof <- -1 * eigen_data$vectors[, k]
    dim(eof) <- c(length(eigen_data$lon), length(eigen_data$lat))
    out <- list(lon = eigen_data$lon, lat = eigen_data$lat, eof = eof, k = k)
    assign_class(out, 'scda_eigen_eofs')
}

#' Display the spatial patterns (EOFs).
#'
#' Plot a specific spatial patterns (EOFs) on map.
#' 
#' @param eigen_data object of class \code{scda_eigen}.
#' @param k integer, the \code{nth} EOFs to be displayed.
#' @param zlim a numeric vector of length 2, giving the minimum and maximum limit of data to display.
#' The values under or above these limits are masked with the corresponding colors assigned to these limits.
#' @param col_palette the color palette function to be used. User defined color palette or a predefined color palettes from the package \pkg{grDevices} or other packages.
#' @param ... other arguments that may be passed on to the function \code{\link[fields]{image.plot}} of the package \pkg{fields}.
#' @param add_map logical, default \code{TRUE} add world map to the plot.
#' @param map_options a named list, if \code{add_map = TRUE}, a list of options to be passed on the function \code{\link[maps]{map}} of the package \pkg{maps}.
#'
#' @examples
#' \dontrun{
#' path_dir_data = '/home/data/ERSSTv5'
#' nc_file_fromat = 'ersst.v5.%s%s.nc'
#' years_range = c(1981, 2020)
#' grd_data = get_month_region_netcdf_data(12, path_dir_data, nc_file_fromat, years_range)
#' eigen_data = eigen_analysis(grd_data)
#' plot_EOFs_spatialPatterns(eigen_data, 2)
#' }
#'  
#' @export

plot_EOFs_spatialPatterns <- function(eigen_data, k, zlim = c(-0.05, 0.05),
                                      col_palette = color_Brewer_Fun, ...,
                                      add_map = TRUE,
                                      map_options = list(col = 'gray95', fill = TRUE,
                                                         lwd = 0.2, border = 'black'))
{
    if(!inherits(eigen_data, 'scda_eigen'))
        stop("'eigen_data' is not an object of class 'scda_eigen'.")

    eof <- get_EOFs_spatialPatterns(eigen_data, k)
    brks <- seq(zlim[1], zlim[2], length.out = 100)
    kolor <- col_palette(length(brks) - 1)
    eof$eof[eof$eof > zlim[2]] <- zlim[2]
    eof$eof[eof$eof < zlim[1]] <- zlim[1]

    data <- eof[c('lon', 'lat', 'eof')]
    names(data) <- c('x', 'y', 'z')
    args <- list(...)
    args <- args[!names(args) %in% c('breaks', 'nlevel', 'col')]
    args$breaks <- brks
    args$col <- kolor
    if(is.null(args$xlab)) args$xlab <- ''
    if(is.null(args$ylab)) args$ylab <- ''
    if(is.null(args$xaxt)){
        args$xaxt <- 'n'
        x_axis <- TRUE
    }
    if(is.null(args$yaxt)){
        args$yaxt <- 'n'
        y_axis <- TRUE
    }
    do.call(fields::image.plot, c(data, args))

    if(x_axis){
        xTck <- graphics::axTicks(1)
        xLab <- format_x_axis_labels(xTck)
        graphics::axis(1, at = xTck, labels = xLab)
    }
    if(y_axis){
        yTck <- graphics::axTicks(2)
        yLab <- format_y_axis_labels(yTck)
        graphics::axis(2, at = yTck, labels = yLab, las = 1)
    }
    graphics::box(bty = 'l')

    if(add_map){
        opts <- list(database = 'world2', col = 'gray95',
                     fill = TRUE, bg = 'white', wrap = c(-180, 180),
                     lwd = 0.2, border = 'black')
        map_options <- init_default_list_args(map_options, opts)
        map_options$add <- TRUE
        do.call(maps::map, map_options)
    }
}
