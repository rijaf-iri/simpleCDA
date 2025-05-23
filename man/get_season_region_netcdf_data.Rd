% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_data.R
\name{get_season_region_netcdf_data}
\alias{get_season_region_netcdf_data}
\title{Read seasonal netCDF dataset.}
\usage{
get_season_region_netcdf_data(
  season,
  nc_path_dir,
  nc_file_fromat,
  years_range = NULL,
  region_bbox = NULL,
  var_name = NA,
  lon_name = NA,
  lat_name = NA,
  ...
)
}
\arguments{
\item{season}{character, the season to be read, in the format \code{"start_month-end_month"}, e.g. "12-02" for DJF.}

\item{nc_path_dir}{character, full path to the folder containing the netCDF files.}

\item{nc_file_fromat}{character, format of the netCDF file name, the year and month must be replaced by \code{\%s}.
Example: for \emph{ersstv5_2024-01_2024-03.nc} the file name format should be \emph{'ersstv5_\%s-\%s_\%s-\%s.nc'}.}

\item{years_range}{an integer vector of length 2, the start and end year of the data to be extracted.
Default \code{NULL}, get the whole period available in the dataset.}

\item{region_bbox}{a numeric vector of length 4, the bounding box of the region to be extracted,
with longitude min, latitude min, longitude max and latitude max of the region (West, South, East, North).
Default \code{NULL}, no extraction performed, get the entire region from the dataset.}

\item{var_name}{character, name of the variable to be read from the netCDF data. If unspecified or left \code{NA},
the first variable available in the netCDF data will be taken.}

\item{lon_name}{character, the name of the longitude dimension. If unspecified or left \code{NA}, it will be detected.}

\item{lat_name}{character, the name of the latitude dimension. Same as \code{lon_name}.}

\item{...}{pairs of arguments, providing the names and values of any extra dimensions to be extracted.
The arguments should be in the format \code{<information about the dimension>_name} for the name of the dimension
and \code{<information about the dimension>_value} for the value to be extracted from that dimension.
Example: if your netCDF data has a pressure levels dimension named \code{"pres"} having the values from 1000 to 10 hPa,
if you want to extract the data at the pressure level 850 hPa, 
the pair of arguments would be \code{plev_name = "pres"} and \code{plev_value = 850}.}
}
\value{
This returns an object of class \code{scda_data}. The object has the following elements:
\itemize{
\item \code{$lon}: a numeric vector of the longitude with length \code{N}.
\item \code{$lat}: a numeric vector of the latitude with length \code{M}.
\item \code{$date}: a vector of class \code{Date} with length \code{L}.
\item \code{$data}: a list of matrices with dimension \code{N x M}, the length of the list is \code{L}.
}
}
\description{
Read seasonal netCDF data for a specific season over a defined region.
}
\examples{
\dontrun{
# Extracting the SST data for the season DJF for the period 1981-2020 
# over the region: longitude [170 W - 120 W] and latitude [5 S - 5 N]

nc_path_dir = '/home/data/ERSSTv5_seasonal'
nc_file_fromat = 'ersstv5_\%s-\%s_\%s-\%s.nc'
years_range = c(1981, 2020)
region_bbox = c(-170, -5, -120, 5)

grd_data = get_season_region_netcdf_data("12-02", nc_path_dir, nc_file_fromat,
                                        years_range, region_bbox)
}

}
