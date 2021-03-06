% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/composite_ppi.R
\name{composite_ppi}
\alias{composite_ppi}
\title{Create a composite of multiple plan position indicators (\code{ppi})}
\usage{
composite_ppi(
  x,
  param = "DBZH",
  nx = 100,
  ny = 100,
  xlim,
  ylim,
  res,
  crs,
  raster = NA,
  method = "max",
  idp = 2,
  idw_max_distance = NA
)
}
\arguments{
\item{x}{A list of \code{ppi} objects.}

\item{param}{Scan parameter to composite.}

\item{nx}{number of raster pixels in the x (longitude) dimension}

\item{ny}{number of raster pixels in the y (latitude) dimension}

\item{xlim}{x (longitude) range}

\item{ylim}{y (latitude) range}

\item{res}{numeric vector of length 1 or 2 to set the resolution of the raster (see \link[raster]{res}).
If this argument is used, arguments \code{nx} and \code{ny} are ignored. Unit is identical to \code{xlim} and \code{ylim}.}

\item{crs}{character or object of class CRS. PROJ.4 type description of a Coordinate Reference System (map projection).
When 'NA' (default), an azimuthal equidistant projection with origin at the radar location is used.
To use a WSG84 (lat,lon) projection, use crs="+proj=longlat +datum=WGS84"}

\item{raster}{(optional) RasterLayer with a CRS. When specified this raster topology is used for the output, and nx, ny, res
arguments are ignored.}

\item{method}{string. Compositing method, one of "mean", "min", "max" or "idw"}

\item{idp}{numeric. inverse distance weighting power}

\item{idw_max_distance}{numeric. Maximum distance from the radar to consider in
inverse distance weighting. Measuruments beyond this distance will have a
weighting factor of zero.}
}
\value{
A \code{\link[=summary.ppi]{ppi}}.
}
\description{
Combines multiple plan position indicators (\code{ppi}) into a single
\code{ppi}. Can be used to make a composite of \code{ppi}'s from multiple
radars.
}
\details{
This function composites multiple ppi objects into a ppi object that
combines all data.

Either multiple ppi's of different scan elevation of the same radar may be combined,
or ppi's of different radars can be composited.

Argument \code{method} determines how values of different ppi's at the same
geographic location are combined.
\describe{
\item{\code{"mean"}}{Compute the average value}
\item{\code{"max"}}{Compute the maximum value. If ppi's are of the same radar
and the same polar volume, this computes a max product, showing the maximum
detected signal at that geographic location.}
\item{\code{"min"}}{Compute the minimum value}
\item{\code{"idw"}}{This option is useful primarily when compositing ppi's of
multiple radars. Performs an inverse distance weighting, where values are
weighted according to 1/(distance from the radar)^\code{idp}}
}

The coordinates system of the returned \code{ppi} is a WGS84
(lat,lon) datum.

This function is a prototype and under active development
}
\examples{
# locate example volume file:
pvolfile <- system.file("extdata", "volume.h5", package = "bioRad")

# load the file:
example_pvol <- read_pvolfile(pvolfile)

# calculate a ppi for each elevation scan
my_ppis <- lapply(example_pvol$scans, project_as_ppi)

# overlay the ppi's, calculating the maximum value observed
# across the available scans at each geographic location
my_composite <- composite_ppi(my_ppis, method="max")

\dontrun{
# download basemap
bm <- download_basemap(my_composite)

# plot the calculated max product on the basemap
map(my_composite, bm)
}
}
