## ----setup, echo=FALSE, message=FALSE-----------------------------------------
SHOW_ANSWERS <- FALSE
if (Sys.info()["sysname"] == "Linux") prefix <- "/home/adriaan" else prefix <- "/Users/amd427"
if (SHOW_ANSWERS) knitr::opts_knit$set(root.dir = normalizePath(paste(prefix, "/Dropbox/RadAero19/bioRad practical/data", sep = "")))
# knitr::opts_chunk$set(eval=FALSE)
Sys.setenv(TZ = "UTC")
library(bioRad)

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # make sure you start with a fresh R session
#  # load the bioRad package
#  library(bioRad)
#  # check the package version
#  packageVersion("bioRad")

## ---- eval=FALSE--------------------------------------------------------------
#  library(devtools)
#  install_github("adokter/bioRad")

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # bring up the package general help page:
#  ?bioRad

## ---- eval=FALSE--------------------------------------------------------------
#  # make a new local directory on your machine where to download data for this practical
#  # replace the string below with the path of that directory:
#  HOME <- "your/personal/working/directory/"
#  # check that the directory exists. If the next statement evaluates to FALSE, something went wrong: the directory does not exist or you didn't specify its path correctly
#  file.exists(HOME)
#  # we will make HOME our work directory, the default folder in which R will look
#  # for files, and where it will output newly generated files.
#  setwd(HOME)
#  # Finally, we set the local time zone to UTC, so all plotted time axes will be in UTC
#  Sys.setenv(TZ = "UTC")

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # start your local Docker installation
#  # we first test whether R can communicate with Docker:
#  check_docker()

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # the bioRad package comes with an example radar volume file, that we will inspect first
#  # first locate this example file on our computer:
#  my_filename <- system.file("extdata", "volume.h5", package = "bioRad")
#  # print the local path of the volume file:
#  my_filename
#  # load the file into R:
#  my_pvol <- read_pvolfile(my_filename)
#  ## print some information about the polar volume
#  my_pvol
#  # print information about the polar scans contained in this polar volume:
#  my_pvol$scans

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS------------------------------------
#  my_pvol <- read_pvolfile("example_pvol.h5")

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS------------------------------------
#  # The default summary information of a `pvol` object contains information
#  # on the scans (sweeps) and their moments:
#  my_pvol

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS------------------------------------
#  # We can also extract the elevation angles from the polar volume as follows:
#  get_elevation_angles(my_pvol)

## ---- eval=SHOW_ANSWERS, warning=FALSE----------------------------------------
#  # (if you haven't done so already) load the polar volume data from the example_pvol.h5 file you just downloaded
#  my_pvol <- read_pvolfile("example_pvol.h5")
#  # let's extract the third scan, which was collected at 1.5 degree elevation:
#  my_scan <- my_pvol$scans[[3]]
#  # print some information about this scan:
#  my_scan
#  # let's plot the reflectivity factor parameter of the scan in a range - azimuth coordinate system:
#  plot(my_scan, param = "DBZH")

## ---- eval=SHOW_ANSWERS, warning=FALSE----------------------------------------
#  # before we can plot the scan, we need to project it on a Cartesian grid,
#  # i.e. we need to make a Plan Position Indicator (PPI)
#  my_ppi <- project_as_ppi(my_scan)
#  # print some information about this ppi:
#  my_ppi
#  # you can see we projected it on a 500 meter grid
#  # (check the manual of the project_as_ppi function to see how you can
#  # change the grid size (argument grid_size) and the maximum distance
#  # from the radar up to where to plot data (argument range_max))
#  #
#  # Now we are ready to plot the ppi, for example let's plot reflectivity factor DBZH:
#  plot(my_ppi, param = "DBZH")

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # Plot the correlation coefficient (RHOHV):
#  plot(my_ppi, param = "RHOHV")

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # Plot the radial velocity (VRADH):
#  plot(my_ppi, param = "VRADH")

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # Answer:
#  #
#  # The texture (spatial variability) of the radial velocity is considerably smoother
#  # in areas with precipitation than in areas with biology. Note: we see this especially
#  # at C-band radars (as in this example), at S-band the difference can be less clear.

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # Answer:
#  # The radial velocity (VRADH) PPI shows that biological scatterers have
#  # a higher speed than the precipitation. This indicates the
#  # biological scatterers must have a high self-propelled speed, which is
#  # typical for birds, not for insects.

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # It is often informative to plot radar data on a base layer.
#  # first download the background image:
#  basemap <- download_basemap(my_ppi)
#  # then overlay the PPI on the satellite image, restricting the color scale from -20 to 15 dBZ:
#  map(my_ppi, map = basemap, param = "DBZH", zlim = c(-20, 15))

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # Usually we would load processed vertical profiles (vp files) by:
#  # my_vplist <- read_vpfiles("./your/directory/with/processed/profiles/goes/here")
#  # my_vplist contains after running the command a list of vertical profile (vp) objects
#  # To save time, we load these data directly from file
#  load("KBRO20170514.RData")
#  # print the length of the vplist object. It should contain 95 profiles
#  length(my_vplist)

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # let's extract a profile from the list, in this example the 41st profile:
#  my_vp <- my_vplist[[41]]
#  # print some info for this profile to the console
#  my_vp
#  # test whether this profile was collected at day time:
#  check_night(my_vp)
#  # plot the vertical profile, in terms of reflectivity factor
#  plot(my_vp, quantity = "dbz")
#  # plot the vertical profile, in terms of (linear) reflectivity
#  plot(my_vp, quantity = "eta")

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # let'splot the vertical profile, in terms of bird density
#  plot(my_vp, quantity = "dens")
#  # print the currently assumed radar cross section (RCS) per bird:
#  rcs(my_vp)

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # Answer:
#  # All bird densities will be a factor 10 times lower.

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # let's change the RCS to 110 cm^2
#  rcs(my_vp) <- 110

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # After changing the RCS above, we simply plot the vertical profile again
#  plot(my_vp)
#  # Indeed the densities are scaled down by a factor 10

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # convert the list of vertical profiles into a time series:
#  my_vpts <- bind_into_vpts(my_vplist)
#  # print summary information
#  my_vpts
#  # time series objects can be subsetted, just as you may be used to with vectors
#  # here we subset the first 50 timesteps:
#  my_vpts[1:50]
#  # here we extract a single timestep, which gives you back a vertical profile class object:
#  my_vpts[100]
#  # to extract all the dates of all profiles in the time series:
#  my_vpts$datetime
#  # to plot the full time series:
#  plot(my_vpts)
#  # check the help file for the plotting function of profile time series
#  # Because profile timeseries are of class 'vpts', it's associated plotting function
#  # is called plot.vpts:
#  ?plot.vpts

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # make a subselection for night time only
#  index_night <- check_night(my_vpts)
#  # index_night is a logical vector that specifies each profile whether it occurred at night or not:
#  index_night
#  # now subset our vpts using this selection:
#  my_vpts_night <- my_vpts[index_night]
#  # plot this smaller time series:
#  plot(my_vpts_night)

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # Answer:
#  #
#  # At 1500 meter 6 UTC the wind barbs have 2 full flags and one half flag.
#  # Therefore the ground speed is approximately 2x5 + 2.5 = 12.5 m/s.

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # First extract the profile at 6 UTC:
#  vp_6UTC <- filter_vpts(my_vpts_night, nearest = "2017-05-14 06:00")
#  plot(vp_6UTC, quantity = "ff")

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # Plot the ground speed (ff):
#  plot(vp_6UTC, quantity = "ff")
#  # Speeds at 1500 metre is approximately 12 m/s, very close to our earlier reading above of 12.5 m/s.

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # Let's continue with the vpts object created in the previous example.
#  # The vertically integrated quantities are calculated as follows:
#  my_vpi <- integrate_profile(my_vpts)
#  # The my_vpi object you created is a vpi class object, which is an acronym for "vertical profile integrated". It has its own plot method, which by default plots migration traffic rate (MTR):
#  plot(my_vpi)
#  # you can also plot vertically integrated densities (VID):
#  plot(my_vpi, quantity = "vid")
#  # the gray and white shading indicates day and night, which is calculated
#  # from the date and the radar position. You can also turn this off:
#  plot(my_vpi, night_shade = FALSE)
#  # plot the cumulative number of birds passing the radar, i.e. migration traffic (mt):
#  plot(my_vpi, quantity = "mt")
#  # execute `?plot.vpi` to open the help page listing all the options.
#  ?plot.vpi

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # Answer:
#  #
#  # VID = (200 birds / km^3) * (1 km) + (100 birds / km^3) * (0.5 km)
#  #     = 250 birds / km^2

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # Answer:
#  #
#  # MTR = (200 birds / km^3) * (50 km / hour) * (1 km) + (100 birds / km^3) * (100 km / hour) * (0.5 km)
#  #     = 10000 birds / km / hour + 5000 birds / km / hour
#  #     = 15000 birds / km / hour

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # Answer:
#  #
#  # MT = MTR * (3 hour) = (15000 birds / km / hour) * (3 hour)
#  #    = 45000 birds / km

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # instead of vertically integrated density (VID), you can use vertically integrated reflectivity (VIR):
#  plot(my_vpi, quantity = "vir")
#  # instead of migration traffic rate (MTR), you can use the reflectivity traffic rate (RTR):
#  plot(my_vpi, quantity = "rtr")
#  # instead of migration traffic (MT), you can use the reflectivity traffic (RT):
#  plot(my_vpi, quantity = "rt")

## ---- eval=SHOW_ANSWERS-------------------------------------------------------
#  # load a time series for the KBGM radar in Binghamton, NY
#  load("KBGM20170527-20170602.RData")
#  # print the loaded vpts time series for this radar:
#  my_vpts
#  # plot the bird density over time:
#  plot(my_vpts, quantity = "dens")

## ---- echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE---------------------
#  # also show meteorological signals:
#  plot(my_vpts, quantity = "DBZH")
#  
#  # Periods with high reflectivities extending to high altitudes indicate precipitation,
#  # i.e. second half of the second night, and on and off during the fourth night.

## ---- eval=FALSE--------------------------------------------------------------
#  # start your local Docker installation
#  # we first test whether R can communicate with Docker:
#  check_docker()

## ---- eval=FALSE--------------------------------------------------------------
#  # download a polar volume file you want to process and put it in your home directory.
#  my_file <- "write_your_file_to_be_processed_here"
#  # Alternatively, continue with the polar volume that comes with the package:
#  my_file <- system.file("extdata", "volume.h5", package = "bioRad")
#  # run vol2bird
#  # we set autoconf to TRUE, to let vol2bird figure out the optimal settings by itself
#  my_vp <- calculate_vp(my_file, autoconf = TRUE)
#  # vp is now a 'vp' profile object, that you can examine as in the previous exercises
#  # alternatively, you may also store the profile as a hdf5 file, which is what we will do next:
#  calculate_vp(my_file, "my_vpfile.h5", autoconf = TRUE)
#  # your work directory should now contain a new file 'my_vpfile.h5'
#  # check that we can read this file, and retrieve the vertical profile from it:
#  vp <- read_vpfiles("my_vpfile.h5")

## ---- eval=FALSE--------------------------------------------------------------
#  # read the filenames of the polar volumes you want to process
#  my_files <- list.files("your/directory/with/volumes/goes/here/", full.names = TRUE)
#  # print the filenames
#  my_files
#  # create output directory for processed profiles
#  outputdir <- "~/processed_data"
#  dir.create(outputdir)
#  # let's loop over the files and generate profiles
#  for (file_in in my_files) {
#    # generate the output filename for the input file
#    file_out <- paste(outputdir, "/", basename(file_in), "_vp.h5", sep = "")
#    # generate the profile, and write it as a hdf5 file (for later reference)
#    # we turn autoconfiguration on, so vol2bird chooses the optimal settings for the file automatically
#    vp <- calculate_vp(file_in, file_out, autoconf = TRUE)
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  # we assume outputdir contains the path to the directory with processed profiles
#  my_vpfiles <- list.files(outputdir, full.names = TRUE)
#  # print them
#  my_vpfiles
#  # read them
#  my_vplist <- read_vpfiles(my_vpfiles)

## ---- eval=FALSE--------------------------------------------------------------
#  # make a time series of profiles:
#  my_vpts <- bind_into_vpts(my_vplist)
#  # plot them between 0 - 3 km altitude:
#  plot(my_vpts, ylim = c(0, 3000))

