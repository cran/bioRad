% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvolfile.R
\name{get_odim_object_type}
\alias{get_odim_object_type}
\title{Check the \code{data} type of an ODIM HDF5 file}
\usage{
get_odim_object_type(file)
}
\arguments{
\item{file}{Character. Path of the file to check.}
}
\value{
Character. \code{PVOL} for polar volume, \code{VP} for vertical profile,
otherwise \code{NA}.
}
\description{
Checks what type of \code{data} object is contained in an ODIM HDF5 file. See
\href{https://github.com/adokter/vol2bird/blob/master/doc/OPERA2014_O4_ODIM_H5-v2.2.pdf}{ODIM specification},
Table 2 for a full list of existing ODIM file object types.
}
\examples{
# Locate the polar volume example file
pvolfile <- system.file("extdata", "volume.h5", package = "bioRad")

# Check the data type
get_odim_object_type(pvolfile)
}
\seealso{
\itemize{
\item \code{\link[=is.pvolfile]{is.pvolfile()}}
\item \code{\link[=is.vpfile]{is.vpfile()}}
}
}
