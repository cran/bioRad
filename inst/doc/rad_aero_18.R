## ---- eval=FALSE---------------------------------------------------------
#  # make sure you start with a fresh R session
#  # load the bioRad package
#  library(bioRad)
#  # check the package version
#  packageVersion("bioRad")
#  # make sure you have the latest version (0.3.0). If you have an older version, download and reinstall as follows:
#  library(devtools)
#  install_github("adokter/bioRad")

## ---- eval=FALSE---------------------------------------------------------
#  # bring up the package general help page:
#  ?bioRad

## ---- eval=FALSE---------------------------------------------------------
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

## ---- eval=FALSE---------------------------------------------------------
#  # start your local Docker installation
#  # we first test whether R can communicate with Docker:
#  check_docker()

## ---- eval=FALSE---------------------------------------------------------
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

## ---- eval=FALSE---------------------------------------------------------
#  # (if you haven't done so already) load the polar volume data from the volume.h5 file you just downloaded
#  my_pvol <- read_pvolfile("example_pvol.h5")
#  # let's extract the third scan, which was collected at 1.5 degree elevation:
#  my_scan <- my_pvol$scans[[3]]
#  # print some information about this scan:
#  my_scan
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

## ---- eval=FALSE---------------------------------------------------------
#  # It is often informative to plot radar data on a base layer, such as google earth maps.
#  # first download the background image:
#  satellite_basemap <- download_basemap(my_ppi, maptype = "satellite")
#  # then overlay the PPI on the satellite image:
#  map(my_ppi, map = satellite_basemap, param = "DBZH", zlim = c(-20, 15))

## ---- eval=FALSE---------------------------------------------------------
#  # Usually we would load processed vertical profiles (vp files) by:
#  # my_vplist <- read_vpfiles("./your/directory/with/processed/profiles/goes/here")
#  # my_vplist contains after running the command a list of vertical profile (vp) objects
#  # To save time, we load these data directly from file
#  load("KBRO20170514.RData")
#  # print some information on the vplist object. It should contain 95 profiles
#  my_vplist

## ---- eval=FALSE---------------------------------------------------------
#  # let's extract a profile from the list, in this example the 41st profile:
#  my_vp <- my_vplist[[41]]
#  # print some info for this profile to the console
#  my_vp
#  # test whether this profile was collected at day time:
#  check_night(my_vp)
#  # plot the vertical profile, in terms of reflectivity factor
#  plot(my_vp, quantity = "dbz")
#  # plot the vertical profile, in terms of reflectivity
#  plot(my_vp, quantity = "eta")

## ---- eval=FALSE---------------------------------------------------------
#  # plot the vertical profile, in terms of bird density
#  plot(my_vp, quantity = "dens")
#  # print the currently assumed radar cross section (RCS) per bird:
#  rcs(my_vp)

## ---- eval=FALSE---------------------------------------------------------
#  # let's change the RCS to 110 cm^2
#  rcs(my_vp) <- 110

## ---- eval=FALSE---------------------------------------------------------
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

## ---- eval=FALSE---------------------------------------------------------
#  # make a subselection for night time only
#  index_night <- check_night(my_vpts)
#  # index_night is a logical vector that specifies each profile whether it occurred at night or not:
#  index_night
#  # now subset our vpts using this selection:
#  my_vpts_night <- my_vpts[index_night]
#  # plot this smaller time series:
#  plot(my_vpts_night)

## ---- eval=FALSE---------------------------------------------------------
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

## ---- eval=FALSE---------------------------------------------------------
#  # print the currently assumed radar cross section:
#  rcs(my_vpi)
#  # instead of VID, you can use vertically integrated reflectivity (VIR):
#  plot(my_vpi, quantity = "vir")
#  # instead of MTR, you can use the reflectivity traffic rate (RTR):
#  plot(my_vpi, quantity = "rtr")
#  # instead of MT, you can use the reflectivity traffic (RT):
#  plot(my_vpi, quantity = "rt")

## ---- eval=FALSE---------------------------------------------------------
#  # load a time series for the KBGM radar in Binghamton, NY
#  load("KBGM20170527-20170602.RData")
#  # print the loaded vpts time series for this radar:
#  my_vpts
#  # plot the bird density over time:
#  plot(my_vpts, quantity = "dens")
#  # also show meteorological signals:
#  plot(my_vpts, quantity = "DBZH")

## ---- eval=FALSE---------------------------------------------------------
#  # start your local Docker installation
#  # we first test whether R can communicate with Docker:
#  check_docker()

## ---- eval=FALSE---------------------------------------------------------
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

## ---- eval=FALSE---------------------------------------------------------
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

## ---- eval=FALSE---------------------------------------------------------
#  # we assume outputdir contains the path to the directory with processed profiles
#  my_vpfiles <- list.files(outputdir, full.names = TRUE)
#  # print them
#  my_vpfiles
#  # read them
#  my_vplist <- read_vpfiles(my_vpfiles)

## ---- eval=FALSE---------------------------------------------------------
#  # make a time series of profiles:
#  my_vpts <- bind_into_vpts(my_vplist)
#  # plot them between 0 - 3 km altitude:
#  plot(my_vpts, ylim = c(0, 3000))

