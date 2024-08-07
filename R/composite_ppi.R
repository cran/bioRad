#' Create a composite of multiple plan position indicators (`ppi`)
#'
#' Combines multiple plan position indicators (`ppi`) into a single
#' `ppi`. Can be used to make a composite of `ppi`'s from multiple
#' radars.
#'
#' @inheritParams integrate_to_ppi
#' @param x A list of `ppi` objects.
#' @param param Character (vector). One or more parameter name(s) to composite.
#'   To composite all available scan parameters use `all` (default).
#' @param method Character (vector). Compositing method(s), either `mean`,
#'   `min`, `max` or `idw`. To apply different methods for each of the
#'   parameters, provide a vector with the same length as `param`.
#' @param idw_max_distance Numeric. Maximum distance from the radar to consider
#'   in inverse distance weighting. Measurements beyond this distance will have
#'   a weighting factor of zero.
#' @param idp Numeric. Inverse distance weighting power.
#' @param coverage Logical. When `TRUE`, adds an additional `coverage` parameter
#'   to the `ppi` indicating the number of `ppi`s covering a single composite
#'   pixel.
#'
#' @return A `ppi` object.
#'
#' @export
#'
#' @details
#' The function can combine multiple `ppi`s of different scan elevations of the
#' same radar or `ppi`s of different radars. The coordinates of the returned
#' `ppi` object are in the WGS84 datum, unless a different `crs` is provided. If
#' only `res` is provided, but no `crs` is set, `res` is in meters and the
#' origin of the composite `ppi` is set to the mean(lat, lon) location.
#'
#' The `method` parameter determines how values of different `ppi`s at the same
#' geographic location are combined:
#' * `mean`: Compute the average value.
#' * `max`: Compute the maximum value. If `ppi`s are of the same radar and the
#' same polar volume, this computes a max product, showing the maximum detected
#' signal at that geographic location.
#' * `min`: Compute the minimum value.
#' * `idw`: This option is useful primarily when compositing `ppi`s of multiple
#' radars. Performs an inverse distance weighting, where values are weighted
#' according to 1/(distance from the radar)^`idp`.
#'
#' Argument `method` determines how values of different ppi's at the same
#' geographic location are combined.
#'
#' * `mean`: Compute the average value
#' * `max`: Compute the maximum value. If ppi's are of the same radar and the
#'   same polar volume, this computes a max product, showing the maximum
#'   detected signal at that geographic location.
#' * `min`: Compute the minimum value
#' * `idw`: This option is useful primarily when compositing ppi's of multiple
#'   radars. Performs an inverse distance weighting, where values are weighted
#'   according to 1/(distance from the radar)^`idp`
#'
#' The coordinates system of the returned `ppi` is a WGS84
#' (lat, lon) datum, unless a different `crs` is provided. If only
#' `res` is provided, but no `crs` is set, `res` is in
#' meter units and the origin of the composite `ppi` is set to the
#' mean (lat, lon) location.
#'
#' This function is a prototype and under active development
#'
#' @examples
#' # Locate and read the polar volume example file
#' \donttest{
#' pvolfile <- system.file("extdata", "volume.h5", package = "bioRad")
#' pvol <- read_pvolfile(pvolfile)
#'
#' # Calculate a ppi for each elevation scan, 1000m grid resolution
#' ppis <- lapply(pvol$scans, project_as_ppi, grid_size=1000)
#'
#' # Overlay the ppis, calculating the maximum value observed
#' # across the available scans at each geographic location
#' composite <- composite_ppi(ppis, method = "max", res=1000)
#'
#' # Plot the calculated max product on the basemap
#' if (all(sapply(c("ggspatial","prettymapr", "rosm"), requireNamespace, quietly = TRUE))) {
#' map(composite)
#' }
#' }
composite_ppi <-
  function(x,
           param = "all",
           nx = 100,
           ny = 100,
           xlim,
           ylim,
           res,
           crs,
           raster = NA,
           method = "max",
           idp = 2,
           idw_max_distance = NA,
           coverage = FALSE) {
  if (FALSE %in% sapply(x, is.ppi)) {
    stop("'x' should be an object of class ppi, or a list of objects of class ppi")
  }
  if (!assertthat::is.count(nx) && missing(res)) stop("'nx' should be an integer")
  if (!assertthat::is.count(ny) && missing(res)) stop("'ny' should be an integer")
  if (!missing(xlim)) {
    if (length(xlim) != 2 ||
        !is.numeric(xlim)) {
      stop("'xlim' should be a numeric vector of length two")
    }
    if (is.na(xlim[1]) |
        is.na(xlim[2]) |
        xlim[1] > xlim[2]) {
      stop(paste("'xlim' should be a vector with two numeric values",
                 "for lower and upper bound respectively"))
    }
  }
  if (!missing(ylim)) {
    if (length(ylim) != 2 ||
        !is.numeric(ylim)) {
      stop("'ylim' should be a numeric vector of length two")
    }
    if (is.na(ylim[1]) |
        is.na(ylim[2]) |
        ylim[1] > ylim[2]) {
      stop(
        paste(
          "'ylim' should be a vector with two numeric values for",
          "lower and upper bound respectively"
        )
      )
    }
  }
  if (!missing(res)) {
    assertthat::assert_that(is.numeric(res))
    assertthat::assert_that(length(res) <= 2)
    t_res <- res
  } else {
    t_res <- NULL
  }
  # check crs argument as in raster::raster()
  if (!missing(crs)) {
    t_crs <- sp::CRS(as.character(raster::projection(crs)))
  } else {
    t_crs <- NULL
  }

  if (!all(method %in% c("max", "min", "mean", "idw"))) stop("'method' should be one or multiple of 'max', 'mean', 'min' or 'idw'")
  if (length(method) != length(param) & length(method) != 1) stop("'method' should be of length 1 or length(param)")
  assertthat::assert_that(assertthat::is.flag(coverage))

  if (length(param) == 1 && param == "all") {
    param <- names(x[[1]]$data)
  }
  ppis <- lapply(x, `[.ppi`, i = param)

  lons <- sapply(ppis, function(x) x$geo$bbox["lon", ])
  lats <- sapply(ppis, function(x) x$geo$bbox["lat", ])
  if (!missing(xlim)) {
    lons <- xlim
  }
  if (!missing(ylim)) {
    lats <- ylim
  }
  lons.radar <- sapply(ppis, function(x) x$geo$lon)
  lats.radar <- sapply(ppis, function(x) x$geo$lat)
  elangles <- sapply(ppis, function(x) x$geo$elangle)
  bbox <- matrix(c(min(lons), min(lats), max(lons), max(lats)),
    nrow = 2,
    ncol = 2, dimnames = dimnames(ppis[[1]]$geo$bbox)
  )

  if (!assertthat::are_equal(raster, NA)) {
    r <- raster::raster(raster)
  } else {
    d_crs <- sp::CRS("+proj=longlat +datum=WGS84")
    if (!is.null(t_res) && !is.null(t_crs)) {

      r <- raster::raster(ext = raster::extent(c(min(lons), max(lons), min(lats), max(lats))), crs = t_crs, resolution = t_res)
    } else if (!is.null(t_crs) && is.null(t_res)) {
      r <- raster::raster(ncols = nx, nrows = ny, ext = raster::extent(c(min(lons), max(lons), min(lats), max(lats))), crs = t_crs)
    } else if (is.null(t_crs) && !is.null(t_res)) {
      r <- raster::raster(ext = raster::extent(c(min(lons), max(lons), min(lats), max(lats))), crs = d_crs)
      t_crs <- sp::CRS(paste0("+proj=aeqd +units=m +ellps=WGS84 +lat_0=", mean(lats), " +lon_0=", mean(lons)))
      r <- raster::projectExtent(r, t_crs)
      raster::res(r) <- t_res
    } else {
      r <- raster::raster(ncols = nx, nrows = ny, ext = raster::extent(c(min(lons), max(lons), min(lats), max(lats))), crs = d_crs)
    }
  }

  # initialize all values of the grid to NA
  suppressWarnings(r <- raster::setValues(r, NA))
  assertthat::assert_that(dim(r)[1] > 1, msg = "output raster has only a single pixel in x dimension, increase raster resolution")
  assertthat::assert_that(dim(r)[2] > 1, msg = "output raster has only a single pixel in y dimension, increase raster resolution")

  spGrid = methods::as(r, 'SpatialGridDataFrame')

  names(spGrid@data) <- names(ppis[[1]]$data)[1]

  if (coverage) {
    ppis <-
      lapply(ppis, function(x) {
        x$data$coverage <- 1
        return(x)
      })
    param <- c(param, "coverage")
  }

  # merge
  projs <- sapply(ppis,
    function(x) {
      sp::over(
        suppressWarnings(
          sp::spTransform(
            spGrid,
            sp::CRS(sp::proj4string(x$data))
          )
        ),
        x$data[param]
      )
    }
  )

  for (p in param) {
    if (p == "coverage") next()
    if (length(param) > 1) {
      merged <- projs[p, ]
    } else {
      merged <- projs
    }

    if (length(method) > 1) {
      param_method <- method[match(p, param)]
    } else {
      param_method <- method
    }

    if (param_method == "max") {
      spGrid@data[, p] <-
        do.call(function(...)
          pmax(..., na.rm = TRUE), merged)
    }
    if (param_method == "min") {
      spGrid@data[, p] <-
        do.call(function(...)
          pmin(..., na.rm = TRUE), merged)
    }
    if (param_method == "mean") {
      as.data.frame(merged) %>% rowMeans(na.rm = TRUE) ->
        spGrid@data[, p]
    }
    if (param_method == "idw") {
      brick_data <-
        suppressWarnings(
          raster::brick(raster::brick(spGrid), nl = length(merged))
          )
      brick_weights <- brick_data

      #weights<-raster::pointDistance(as.matrix(data.frame(x=lons.radar,y=lats.radar)), sp::coordinates(raster::raster(spGrid)),lonlat=TRUE)
      for(i in 1:length(merged)){
        brick_data <- raster::setValues(brick_data, merged[[i]], layer=i)
        latlon.radar <- unique(data.frame(lat = c(lats.radar), lon = c(lons.radar)))
        if (is.null(t_res)) {
          weights <- suppressWarnings(raster::pointDistance(as.matrix(data.frame(x = latlon.radar$lon, y = latlon.radar$lat))[i, ], sp::coordinates(raster::raster(spGrid)), lonlat = TRUE))
        } else {
          d <- data.frame(lon = latlon.radar$lon, lat = latlon.radar$lat)
          sp::coordinates(d) <- c("lon", "lat")
          sp::proj4string(d) <- d_crs
          proj.radar <- as.data.frame(sp::spTransform(d, t_crs))
          weights <- suppressWarnings(raster::pointDistance(as.matrix(data.frame(x = proj.radar[,1], y = proj.radar[,2]))[i, ], sp::coordinates(raster::raster(spGrid)), lonlat = FALSE))

        }
        if (!is.na(idw_max_distance)) weights[weights > idw_max_distance] <- NA
        weights <- 1 / (weights^idp)
        brick_weights <- raster::setValues(brick_weights, weights, layer = i)
      }
      spGrid@data[, p] <-
        as.vector(
          raster::weighted.mean(brick_data, brick_weights, na.rm = TRUE)
        )
    }
  }

  if (coverage) {
    cov <- !is.na(do.call("cbind", projs["coverage", ]))
    spGrid@data$coverage <- rowSums(cov)
  }

  ppi.out <- list(data = spGrid, geo = list(
    lat = lats.radar, lon = lons.radar,
    elangle = elangles, bbox = bbox,
    merged = TRUE
  ))
  class(ppi.out) <- "ppi"
  ppi.out
}
