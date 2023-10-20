## -----------------------------------------------------------------------------
options("sp_evolution_status" = 2L)
library(bioRad)

## -----------------------------------------------------------------------------
# define a range grid from 0 to 200 km:
my_range <- seq(0, 200000, 1000)
# plot the beam height for each range, for a 0.5 degree elevation beam:
plot(my_range, beam_height(my_range, elev = .5), xlab = "range [m]", ylab = "beam height [m]")

## -----------------------------------------------------------------------------
# plot the beam width, for a beam opening angle of 1 degree (typical for most weather radars):
plot(my_range, beam_width(my_range, beam_angle = 1), xlab = "range [m]", ylab = "beam width [m]")

## -----------------------------------------------------------------------------
# plot the beam profile, for a 0.5 degree elevation beam at 50 km distance from the radar:
plot(beam_profile(height = 0:4000, 50000, 0.5), 0:4000, xlab = "normalized radiated energy", ylab = "height [m]", main = "beam elevation: 0.5 deg, distance=50km")
# plot the beam profile, for a 2 degree elevation beam at 50 km distance from the radar:
plot(beam_profile(height = 0:4000, 50000, 2), 0:4000, xlab = "normalized radiated energy", ylab = "height [m]", main = "beam elevation: 2 deg, distance=50km")

## -----------------------------------------------------------------------------
# plot the combined beam profile for a 0.5 and 2.0 degree elevation beam at 50 km distance from the radar:
plot(beam_profile(height = 0:4000, 50000, c(0.5, 2)), 0:4000, xlab = "normalized radiated energy", ylab = "height [m]", main = "beam elevations: 0.5,2 deg, distance=50km")

## -----------------------------------------------------------------------------
# let's load an example polar volume:
pvolfile <- system.file("extdata", "volume.h5", package = "bioRad")
example_pvol <- read_pvolfile(file = pvolfile)
# a vertical profile can also be calculated from the polar volume directly, using
# calculate_vp(pvolfile)
# but for now we will use bioRad's example vertical profile already calculated:
example_vp

## -----------------------------------------------------------------------------
plot(example_vp, quantity = "eta")

## -----------------------------------------------------------------------------
dbz_to_eta(5, wavelength = 5.3)

## -----------------------------------------------------------------------------
# extract the first scan from the polar volume:
my_scan <- example_pvol$scans[[1]]
# project it as a PPI on the ground:
my_ppi <- project_as_ppi(my_scan, range_max = 100000)
# plot it
plot(my_ppi)

## -----------------------------------------------------------------------------
# let's use a 500 metre spatial grid (res), and restrict to 100x100 km area
my_corrected_ppi <- integrate_to_ppi(example_pvol, example_vp, res = 500, xlim = c(-100000, 100000), ylim = c(-100000, 100000))
my_corrected_ppi

## -----------------------------------------------------------------------------
# plot the adjustment factor R:
plot(my_corrected_ppi, param = "R")

## -----------------------------------------------------------------------------
plot(my_corrected_ppi, param = "VIR")

## -----------------------------------------------------------------------------
bm <- "osm"
map(my_corrected_ppi, map=bm, param = "VIR", alpha = .5)

## -----------------------------------------------------------------------------
# calculate overlap between vertical profile of birds
# and the vertical radiation profile emitted by the radar:
bpo <- beam_profile_overlap(example_vp, get_elevation_angles(example_pvol), seq(0, 100000, 1000), quantity = "eta")
# plot the calculated overlap:
plot(bpo)

## -----------------------------------------------------------------------------
plot(my_corrected_ppi, param = "overlap")

