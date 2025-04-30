get_latlon_dimnames <- function(nc){
    dim_names <- sapply(nc$dim, '[[', 'name')
    dims <- tolower(dim_names)
    lo <- grep('lon', dims)
    if(length(lo) == 0) lo <- which(dims == 'x')
    if(length(lo) == 0)
        stop('Unknown longitude dimension name.')
    lon <- dim_names[lo[1]]

    la <- grep('lat', dims)
    if(length(la) == 0) la <- which(dims == 'y')
    if(length(la) == 0)
        stop('Unknown latitude dimension name.')
    lat <- dim_names[la[1]]

    list(lon = lon, lat = lat)
}

get_latlon_order <- function(nc, var_name, lon_name, lat_name){
    dim_var <- sapply(nc$var[[var_name]]$dim, '[[', 'name')
    ilon <- which(dim_var == lon_name)
    ilat <- which(dim_var == lat_name)

    list(ilon = ilon, ilat = ilat)
}

format_longitude <- function(lon){
    ((lon + 180) %% 360) - 180
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
    if(any(check)) stop('Unambiguous netCDF file names format.')

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
