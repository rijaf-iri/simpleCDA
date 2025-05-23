% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eigens.R
\name{eigen_analysis}
\alias{eigen_analysis}
\title{Performing Principal Component Analysis (PCA).}
\usage{
eigen_analysis(data, eof_number = 4, ...)
}
\arguments{
\item{data}{object of class \code{scda_data}.}

\item{eof_number}{integer, number of \emph{Empirical Orthogonal Functions} (EOF) requested.}

\item{...}{other arguments that may be passed on to the function \code{\link[RSpectra]{eigs_sym}} of the package \pkg{RSpectra}.}
}
\value{
This returns an object of class \code{scda_eigen}.
In addition to the output from the function \code{\link[RSpectra]{eigs_sym}} of the package \pkg{RSpectra}.
the object has the additional elements:
\itemize{
\item \code{$lon}: a numeric vector of the longitude with length \code{N} from \code{data}.
\item \code{$lat}: a numeric vector of the latitude with length \code{M} from \code{data}.
\item \code{$date}: a vector of class \code{Date} with length \code{L} from \code{data}.
\item \code{$anomaly}: a matrix of the standardized anomalies with dimension \code{L x (N*M)}.
}
}
\description{
Perform eigenanalysis of covariance matrix.
}
\examples{
\dontrun{
# Extracting the SST data for July for the period 1981-2020 
# over the region: longitude [170 W - 120 W] and latitude [5 S - 5 N]

nc_path_dir = '/home/data/ERSSTv5'
nc_file_fromat = 'ersst.v5.\%s\%s.nc'
years_range = c(1981, 2020)
region_bbox = c(-170, -5, -120, 5)

grd_data = get_month_region_netcdf_data(7, nc_path_dir, nc_file_fromat,
                                        years_range, region_bbox)

eigen_data = eigen_analysis(grd_data)
}

}
