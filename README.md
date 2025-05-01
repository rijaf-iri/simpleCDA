# simpleCDA - Simple Climate Data Analysis

## General
`simpleCDA` is a R package for a simple climate data analysis.


## Installation

```r
library(devtools)
install_github("rijaf-iri/simpleCDA")
```

## Quick start

```r
library(simpleCDA)

# Extracting the first PCs of SST for July for the period 1981-2020 
# over the region: longitude [170 W - 120 W] and latitude [5 S - 5 N]

path_dir_data = '/home/data/ERSSTv5'
nc_file_fromat = 'ersst.v5.%s%s.nc'
years_range = c(1981, 2020)
region_bbox = c(-170, -5, -120, 5)

# reading the data
grd_data = get_month_region_netcdf_data(7, path_dir_data, nc_file_fromat, years_range, region_bbox)

# computing the first time series PCs
ts = extract_region_timeSeries(grd_data, output_type = 'pcs', nth_pcs = 1)
```