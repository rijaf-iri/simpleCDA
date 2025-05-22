get_latlon_dimnames <- function(nc){
    dim_names <- sapply(nc$dim, '[[', 'name')
    dims <- tolower(dim_names)
    lo <- grep('lon', dims)
    if(length(lo) == 0)
        lo <- which(dims == 'x')
    if(length(lo) == 0)
        stop('Unknown longitude dimension name.')
    lon <- dim_names[lo[1]]

    la <- grep('lat', dims)
    if(length(la) == 0)
        la <- which(dims == 'y')
    if(length(la) == 0)
        stop('Unknown latitude dimension name.')
    lat <- dim_names[la[1]]

    list(lon = lon, lat = lat)
}

get_dim_extras <- function(args_n, args_v){
    args <- lapply(seq_along(args_n), function(i){
        x <- as.character(args_n[[i]])
        x <- strsplit(x, '_')[[1]]
        ix <- x[length(x)]
        iv <- paste0(x[-length(x)], collapse = '_')
        list(dim = iv, type = ix, value = args_v[[i]])
    })

    xdims <- sapply(args, '[[', 'dim')
    args_f <- lapply(unique(xdims), function(x){
        v <- args[xdims == x]
        if(length(v) == 1){
            type <- if(v[[1]]$type == 'name') 'value' else 'name'
            type1 <- paste0('"', v[[1]]$dim, '_', type, '"')
            type2 <- paste0('"', v[[1]]$dim, '_', v[[1]]$type, '".')
            stop(paste('No', type1, 'found for', type2))
        }
        z <- lapply(v, '[[', 'value')
        names(z) <- sapply(v, '[[', 'type')
        z
    })

    return(args_f)
}

get_dim_ellipsis <- function(...){
    args_v <- list(...)
    if(length(args_v) == 0) return(NULL)
    args_n <- match.call(expand.dots = FALSE)$...
    xname <- names(args_n)
    lname <- xname != ""
    args_n[lname] <- xname[lname]
    get_dim_extras(args_n, args_v)
}

get_dim_list <- function(list_dim){
    if(is.null(list_dim)) return(NULL)
    args_v <- list_dim
    args_n <- as.list(names(list_dim))
    get_dim_extras(args_n, args_v)
}

get_dim_order <- function(nc, var_name, lon_name, lat_name, dim_dots){
    dim_var <- sapply(nc$var[[var_name]]$dim, '[[', 'name')
    ilon <- which(dim_var == lon_name)
    ilat <- which(dim_var == lat_name)

    dim_order <- list(ilon = ilon, ilat = ilat)
    dim_extra <- dim_var[!dim_var %in% c(lon_name, lat_name)]

    if(!is.null(dim_dots)){
        for(i in seq_along(dim_dots)){
            dname <- dim_dots[[i]]$name
            if(dname %in% dim_var){
                dim_order[dname] <- which(dim_var == dname)
            }else{
                stop(paste0('There is no dimension named "', dname, '"'))
            }
            dim_extra <- dim_extra[dim_extra != dname]
        }
    }

    if(length(dim_extra) > 0){
        for(v in dim_extra){
           dim_order[v] <- which(dim_var == v)
        }
    }

    return(dim_order)
}

get_dim_infos_extract <- function(nc, var_name, lon_name, lat_name, bbox, ...){
    if(is.na(lon_name) || is.na(lat_name)){
        xy_name <- get_latlon_dimnames(nc)
        lon_name <- xy_name$lon
        lat_name <- xy_name$lat
    }
    if(is.na(var_name)){
        var_name <- nc$var[[1]]$name
    }

    dim_dots <- get_dim_ellipsis(...)
    if(!is.null(dim_dots)){
        for(j in seq_along(dim_dots)){
            tmp <- nc$dim[[dim_dots[[j]]$name]]$vals
            if(!any(tmp == dim_dots[[j]]$value)){
                stop(paste('Wrong value for the dimension',
                     paste0('"', dim_dots[[j]]$name, '".'),
                     'Select one of the following values:',
                     paste0(tmp, collapse = ', '))
                    )
            }
        }
    }

    pos <- get_dim_order(nc, var_name, lon_name, lat_name, dim_dots)
    index <- lapply(seq_along(pos), function(i) 1)

    lon <- nc$dim[[lon_name]]$vals
    lon <- format_longitude_180(lon)
    lat <- nc$dim[[lat_name]]$vals

    if(!is.null(bbox)){
        bbox[1] <- format_longitude_180(bbox[1])
        bbox[3] <- format_longitude_180(bbox[3])
        plon <- abs(diff(sort(lon)[1:2]))
        if(abs(bbox[3] - bbox[1]) < plon){
            clon <- (bbox[1] + bbox[3])/2
            index[[pos$ilon]] <- find_intervals(clon, lon, plon)
        }else{
            index[[pos$ilon]] <- which(lon >= bbox[1] & lon <= bbox[3])
        }

        plat <- abs(diff(sort(lat)[1:2]))
        if(abs(bbox[4] - bbox[2]) < plat){
            clat <- (bbox[2] + bbox[4])/2
            index[[pos$ilat]] <- find_intervals(clat, lat, plat)
        }else{
            index[[pos$ilat]] <- which(lat >= bbox[2] & lat <= bbox[4])
        }
    }else{
        index[[pos$ilon]] <- seq_along(lon)
        index[[pos$ilat]] <- seq_along(lat)
        rx <- range(lon)
        ry <- range(lat)
        bbox <- c(rx[1], ry[1], rx[2], ry[2])
    }

    if(!is.null(dim_dots)){
        for(j in seq_along(dim_dots)){
            nv <- dim_dots[[j]]$name
            tmp <- nc$dim[[nv]]$vals
            np <- pos[[nv]]
            index[[np]] <- which(tmp == dim_dots[[j]]$value)
        }
    }

    lon <- lon[index[[pos$ilon]]]
    lat <- lat[index[[pos$ilat]]]
    ox <- order(lon)
    oy <- order(lat)
    lon <- lon[ox]
    lat <- lat[oy]
    nlon <- length(lon)
    nlat <- length(lat)

    list(varid = var_name, lon = lon, lat = lat,
         index = index, pos = pos, ox = ox, oy = oy,
         nx = nlon, ny = nlat, bbox = bbox)
}

get_dim_infos_aggregate <- function(nc, var_name, lon_name, lat_name, extra_dim){
    if(is.na(lon_name) || is.na(lat_name)){
        xy_name <- get_latlon_dimnames(nc)
        lon_name <- xy_name$lon
        lat_name <- xy_name$lat
    }
    if(is.na(var_name)){
        var_name <- nc$var[[1]]$name
    }

    dim_dots <- get_dim_list(extra_dim)
    if(!is.null(dim_dots)){
        for(j in seq_along(dim_dots)){
            tmp <- nc$dim[[dim_dots[[j]]$name]]$vals
            if(!any(tmp == dim_dots[[j]]$value)){
                stop(paste('Wrong value for the dimension',
                     paste0('"', dim_dots[[j]]$name, '".'),
                     'Select one of the following values:',
                     paste0(tmp, collapse = ', '))
                    )
            }
        }
    }

    pos <- get_dim_order(nc, var_name, lon_name, lat_name, dim_dots)
    index <- lapply(seq_along(pos), function(i) 1)

    lon <- nc$dim[[lon_name]]$vals
    lon <- format_longitude_180(lon)
    lat <- nc$dim[[lat_name]]$vals

    index[[pos$ilon]] <- seq_along(lon)
    index[[pos$ilat]] <- seq_along(lat)
    if(!is.null(dim_dots)){
        for(j in seq_along(dim_dots)){
            nv <- dim_dots[[j]]$name
            tmp <- nc$dim[[nv]]$vals
            np <- pos[[nv]]
            index[[np]] <- which(tmp == dim_dots[[j]]$value)
        }
    }

    ox <- order(lon)
    oy <- order(lat)
    lon <- lon[ox]
    lat <- lat[oy]
    nlon <- length(lon)
    nlat <- length(lat)

    list(varid = var_name, lon = lon, lat = lat,
         pos = pos, index = index, ox = ox,
         oy = oy, nx = nlon, ny = nlat)
}

extract_nc_data <- function(nc_array, info){
    args <- c(list(x = nc_array), info$index,
              list(drop = TRUE))
    tmp <- do.call(`[`, args)

    if(info$pos$ilon < info$pos$ilat){
        dim(tmp) <- c(info$nx, info$ny)
        tmp <- tmp[info$ox, info$oy, drop = FALSE]
    }else{
        dim(tmp) <- c(info$ny, info$nx)
        tmp <- tmp[info$oy, info$ox, drop = FALSE]
        tmp <- t(tmp)
    }

    return(tmp)
}

find_intervals <- function(x, vec, px = NA){
    nl <- length(vec)
    ov <- order(vec)
    vec <- vec[ov]
    if(is.na(px)) px <- abs(diff(vec[1:2]))
    tmp <- c(vec[1] - px/2, vec + px/2)
    ix <- findInterval(x, tmp,
                       rightmost.closed = TRUE,
                       left.open = TRUE)
    l <- ix == 0
    if(any(l)) ix[l] <- 1
    r <- ix > nl
    if(any(r)) ix[r] <- nl
    ov[ix]
}

format_longitude_180 <- function(lon){
    ((lon + 180) %% 360) - 180
}

format_longitude_360 <- function(lon){
    (lon + 360) %% 360
}

format_pattern <- function(x){
    # exclude %
    symbols <- "[-?.,;:'_+=()!@#$^&*|~`{}]"
    if(!grepl(symbols, x)) return(x)

    symbols <- strsplit(symbols, '')[[1]]
    for(s in symbols){
        x <- gsub(paste0('\\', s), paste0('\\\\', s), x)
    }
    return(x)
}

double_backslash_non_alnum <- function(strings){
    for(i in seq_along(strings)){
        expr <- gregexpr("[^[:alnum:]]", strings[i])
        ex <- expr[[1]]
        if(ex[1] == -1) next

        chr <- rep('', length(ex))
        for(j in seq_along(chr)){
            chr[j] <- substr(strings[i], ex[j], ex[j])
        }
        chr <- chr[!duplicated(chr)]
        for(v in chr){
            pt0 <- paste0('\\', v)
            pt1 <- paste0('\\', pt0)
            strings[i] <- gsub(pt0, pt1, strings[i])
        }
    }

    return(strings)
}

extract_filename_dates <- function(filenames, fileformat){
    expr <- gregexpr('%', fileformat)[[1]]
    len <- rep(2, length(expr))
    ret <- NULL
    if(expr[1] != -1){
        re <- FALSE
        ss <- 1
        se <- nchar(fileformat)
        nl <- length(expr)
        for(i in 1:nl){
            re <- c(re, TRUE, FALSE)
            ss <- c(ss, expr[i], expr[i] + len[i])
            j <- nl - i + 1
            se <- c(expr[j] - 1, expr[j] + len[j] - 1, se)
        }

        res <- lapply(seq_along(re), function(i){
            v <- substr(fileformat, ss[i], se[i])
            if(v == '') v <- NULL
            if(re[i]) v <- NULL
            v
        })

        inul <- sapply(res, is.null)
        if(!all(inul)){
            res <- do.call(c, res[!inul])
            res <- res[!duplicated(res)]
            res <- double_backslash_non_alnum(res)
            pattern <- paste0(res, collapse = '|')
            ret <- gsub(pattern, '', filenames)
        }
    }
    check <- grepl('[^[:digit:]]', ret)
    if(any(check))
        stop('Unambiguous netCDF file names format.')

    return(ret)
}

nbdays_of_months <- function(dates){
    # dates: class Date
    x <- format(dates, '%Y%m')
    x <- lapply(x, function(m){
        em <- 28:31
        y <- paste0(m, em)
        y <- as.Date(y, '%Y%m%d')
        y <- em[!is.na(y)]
        rev(y)[1]
    })
    do.call(c, x)
}

nbdays_of_years <- function(dates){
    # dates: class Date
    y1 <- as.Date(format(dates, '%Y-01-01'))
    y2 <- as.Date(format(dates, '%Y-12-31'))
    n <- difftime(y2, y1, units = 'days')
    as.numeric(n) + 1
}

dates_add_months <- function(dates, length){
    # dates: class Date
    yr <- as.numeric(format(dates, '%Y'))
    mo <- as.numeric(format(dates, '%m'))
    year <- floor(yr + (mo + length - 1) / 12)
    mon <- (mo + length - 1) %% 12 + 1
    day <- format(dates, '%d')
    tmp <- paste(year, mon, '01', sep = '-')
    nd <- nbdays_of_months(as.Date(tmp))
    id <- as.numeric(day) > nd
    day[id] <- nd[id]
    as.Date(paste(year, mon, day, sep = '-'))
}

season_define_dates <- function(start_date, end_date, season_length){
    # start_date, end_date: class Date
    year1 <- as.numeric(format(start_date, '%Y'))
    year2 <- as.numeric(format(end_date, '%Y'))
    seasons <- lapply(year1:year2, function(y){
        s1 <- as.Date(paste0(y, '-', 1:12, '-01'))
        s2 <- dates_add_months(s1, season_length - 1)
        n <- nbdays_of_months(s2)
        s2 <- as.Date(paste0(format(s2, '%Y-%m-'), n))
        data.frame(start = s1, end = s2)
    })
    do.call(rbind, seasons)
}

season_format_dates <- function(seas_dates){
    # seas_dates: yyyymmyyymm or yyyy-mm_yyy-mm
    if(nchar(seas_dates[1]) == 12){
        seas1 <- paste0(substr(seas_dates, 1, 6), '01')
        seas1 <- as.Date(seas1, '%Y%m%d')
        seas2 <- paste0(substr(seas_dates, 7, 12), '01')
        seas2 <- as.Date(seas2, '%Y%m%d')
    }else{
        seas <- strsplit(seas_dates, '_')
        seas <- lapply(seas, function(x) as.Date(paste0(x, '-01')))
        seas1 <- do.call(c, lapply(seas, '[[', 1))
        seas2 <- do.call(c, lapply(seas, '[[', 2))
    }
    nbd <- nbdays_of_months(seas2)
    seas2 <- as.Date(paste0(format(seas2, '%Y-%m-'), nbd))

    list(start = seas1, end = seas2)
}

season_filename_dates <- function(seas_dates){
    # seas_dates: yyyymmyyymm
    seas <- season_format_dates(seas_dates)
    start <- format(seas$start, '%Y-%m')
    end <- format(seas$end, '%Y-%m')

    paste0(start, '_', end)
}

season_middle_dates <- function(seas_dates){
    # seas_dates: yyyymmyyymm or yyyy-mm_yyy-mm
    seas <- season_format_dates(seas_dates)
    ndays <- difftime(seas$end + 1, seas$start, units = 'days')

    seas$start + as.numeric(ndays)/2
}

assign_class <- function(obj, class_name){
    class(obj) <- append(class(obj), class_name)
    obj
}

color_Brewer_Fun <- function(n, color_name = 'RdBu', reverse = TRUE, nc = NULL){
    if(is.null(nc)){
        namelist <- c('BrBG', 'PiYG', 'PRGn', 'PuOr',
                      'RdBu','RdGy', 'RdYlBu', 'RdYlGn',
                      'Spectral','Accent', 'Dark2', 'Paired',
                      'Pastel1', 'Pastel2', 'Set1', 'Set2',
                      'Set3', 'Blues', 'BuGn', 'BuPu',
                      'GnBu', 'Greens', 'Greys', 'Oranges',
                      'OrRd', 'PuBu', 'PuBuGn', 'PuRd',
                      'Purples', 'RdPu', 'Reds', 'YlGn',
                      'YlGnBu', 'YlOrBr', 'YlOrRd')
        maxcolors <- c(11, 11, 11, 11,
                       11, 11, 11, 11,
                       11, 8, 8, 12,
                       9, 8, 9, 8,
                       12, 9, 9, 9,
                       9, 9, 9, 9,
                       9, 9, 9, 9,
                       9, 9, 9, 9,
                       9, 9, 9)
        nc <- maxcolors[which(namelist == color_name)]
    }

    kolor <- RColorBrewer::brewer.pal(nc, color_name)
    if(reverse) kolor <- rev(kolor)
    kolor <- c(kolor[1:5], kolor[6], "#FFFFFF", kolor[6], kolor[7:11])
    foo <- grDevices::colorRampPalette(kolor)
    foo(n)
}

format_x_axis_labels <- function(x){
    ax <- ifelse(x < 0 & x > -180, 'W', ifelse(x > 0 & x < 180, 'E', ''))
    ax <- paste('paste(', paste(abs(x), '*degree'), ',', ax, ')', collapse = ',')
    eval(parse(text = paste0('expression(', ax, ')')))
}

format_y_axis_labels <- function(y){
    ax <- ifelse(y < 0, 'S', ifelse(y > 0, 'N', ''))
    ax <- paste('paste(', paste(abs(y), '*degree'), ',', ax, ')', collapse = ',')
    eval(parse(text = paste0('expression(', ax, ')')))
}

init_default_list_args <- function(inupt_args, init_pars){
    pars_name <- names(init_pars)
    data_name <- names(inupt_args)
    inm <- data_name %in% pars_name
    if(any(inm)){
        for(n in data_name[inm])
            init_pars[[n]] <- inupt_args[[n]]
    }

    return(init_pars)
}
