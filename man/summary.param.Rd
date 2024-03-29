% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param.R
\name{summary.param}
\alias{summary.param}
\alias{is.param}
\title{Inspect a parameter (\code{param})}
\usage{
\method{summary}{param}(object, ...)

is.param(x)
}
\arguments{
\item{object}{A \code{param} object.}

\item{...}{Additional arguments affecting the summary produced.}

\item{x}{A \code{param} object.}
}
\value{
For \code{\link[=is.param]{is.param()}}: \code{TRUE} for an object of class \code{param}, otherwise
\code{FALSE}.
}
\description{
R base functions for inspecting a parameter (\code{param}) object.
}
\details{
A parameter is a quantity/variable measured by the radar during a scan (or
sweep). These are organized along radar range (bins) and azimuth (rays). Scan
parameters are named according to the OPERA data information model (ODIM),
see Table 16 in the \href{https://github.com/adokter/vol2bird/blob/master/doc/OPERA2014_O4_ODIM_H5-v2.2.pdf}{ODIM specification}.

Commonly available parameters are:
\itemize{
\item \code{DBZH}, \code{DBZ}: (Logged) reflectivity factor in dBZ.
\item \code{TH}, \code{T}: (Logged) uncorrected reflectivity factor in dBZ.
\item \code{VRADH}, \code{VRAD}: Radial velocity in m/s. Radial velocities towards the
radar are negative, while radial velocities away from the radar are positive.
\item \code{RHOHV}: Correlation coefficient (unitless). Correlation between the
vertically and horizontally polarized reflectivity factor.
\item \code{PHIDP}: Differential phase in degrees.
\item \code{ZDR}: (Logged) differential reflectivity in dB.
}
}
\examples{
# Extract the DBZH parameter from a scan
param <- get_param(example_scan, "DBZH")

# Check if it is an object of class param
is.param(param)

# Get summary info for this parameter
param # Same as summary(param) or print(param)
}
\seealso{
\itemize{
\item \code{\link[=get_param]{get_param()}}
}
}
