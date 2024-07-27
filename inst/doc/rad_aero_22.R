## ----setup, echo=FALSE, message=FALSE-----------------------------------------
SHOW_ANSWERS <- FALSE
if (Sys.info()["sysname"] == "Linux") prefix <- "/home/adriaan/" else prefix <- "/Users/amd427/"
if (SHOW_ANSWERS) knitr::opts_knit$set(root.dir = normalizePath(paste(prefix, "git/training/colorado2022/", sep = "")))
# knitr::opts_chunk$set(eval=FALSE)
Sys.setenv(TZ = "UTC")
library(bioRad)

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # make sure you start with a fresh R session
#  # load the bioRad package
#  library(bioRad)
#  # check the package version
#  packageVersion("bioRad")

## ----eval=FALSE---------------------------------------------------------------
#  # bring up the package general help page:
#  ?bioRad

## ----eval=FALSE---------------------------------------------------------------
#  # make a new local directory on your machine for this practical
#  # replace the string below with the path of that directory:
#  HOME <- "your/personal/working/directory/"
#  # check that the directory exists. If the next statement evaluates to FALSE, something went wrong: the directory does not exist or you didn't specify its path correctly
#  file.exists(HOME)
#  # we will make HOME our work directory, the default folder in which R will look
#  # for files, and where it will output newly generated files.
#  setwd(HOME)
#  # next, we will create two directories where we will be storing data:
#  dir.create("./data_vpts") # here we will store vertical profile time series (vpts data)
#  dir.create("./data_pvol") # here we will store polar volumes (pvol data)
#  # Finally, we set the local time zone to UTC, so all plotted time axes will be in UTC
#  Sys.setenv(TZ = "UTC")

## ----eval=SHOW_ANSWERS, results='hide'----------------------------------------
#  # Let's first download the NEXRAD polar volume files for the KHGX radar (Houston)
#  # for a 15 minute period in 2017:
#  download_pvolfiles(date_min=as.POSIXct("2017-05-04 01:25:00"), date_max=as.POSIXct("2017-05-04 01:40:00"), radar="KHGX", directory="./data_pvol")
#  # store the filenames in my_pvolfiles
#  my_pvolfiles <- list.files("./data_pvol", recursive = TRUE, full.names = TRUE, pattern="KHGX")
#  # print to console our files:
#  my_pvolfiles
#  # let's load the first of our downloaded files:
#  my_pvol <- read_pvolfile(my_pvolfiles[1])

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS-------------------------------------
#  # The default summary information of a `pvol` object contains information
#  # on the scans (sweeps) and their moments:
#  my_pvol

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS-------------------------------------
#  # We can also extract the elevation angles from the polar volume as follows:
#  get_elevation_angles(my_pvol)

## ----eval=SHOW_ANSWERS, warning=FALSE-----------------------------------------
#  # let's extract the scan collected at 1.5 degree elevation from our polar volume:
#  my_scan <- get_scan(my_pvol, 0.5)
#  # print some information about this scan:
#  my_scan
#  # let's plot the reflectivity factor parameter of the scan in a range - azimuth coordinate system:
#  plot(my_scan, param = "DBZH")

## ----eval=SHOW_ANSWERS, warning=FALSE-----------------------------------------
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

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # Plot the correlation coefficient (RHOHV):
#  plot(my_ppi, param = "RHOHV")

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # Plot the radial velocity (VRADH):
#  plot(my_ppi, param = "VRADH")

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # Answer:
#  # * Precipitation areas are characterized by high correlation coefficients (typically > 0.95),
#  #   i.e. the top left corner of the image is precipitation
#  # * The blue radial velocity of the precipitation indicates it is moving towards the radar.
#  # * The radial velocity field of the biology shows areas south-west of the radar moving
#  #   towards the radar (blue), and areas north-east of the radar moving away from it (red).
#  #   The biology is therefore moving towards the north-east

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # Answer:
#  # * Precipitation moves with the wind field. Since the biology moves into a different direction
#  #   than the wind and a very different speed, we can be sure these are birds.
#  #   The biological scatterers have a high self-propelled speed, which is
#  #   typical for birds, not for insects.
#  # * Note: This is an S-band radar. In C-band radars you will typically see that the texture
#  #   (spatial variability) of the radial velocity is considerably smoother in areas with
#  #   precipitation than in areas with biology.

## ----eval=SHOW_ANSWERS, warning=FALSE-----------------------------------------
#  # It is often informative to plot radar data on a base layer.
#  # First choose a base layer from the list of rosm::osm.types()
#  basemap = "osm"
#  # then overlay the PPI on the basemap, restricting the color scale from -20 to 40 dBZ:
#  map(my_ppi, map = basemap, param = "DBZH", zlim = c(-20, 40))

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # Screen out the reflectivity areas with RHOHV < 0.95
#  my_ppi_clean <- calculate_param(my_ppi, DBZH = ifelse(RHOHV > 0.95, NA, DBZH))
#  # plot the original and cleaned up reflectivity:
#  map(my_ppi, map = basemap, param = "DBZH", zlim = c(-20, 40))
#  map(my_ppi_clean, map = basemap, param = "DBZH", zlim = c(-20, 40))

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # apply the MistNet model to the polar volume file and load it as a polar volume (pvol):
#  my_pvol <- apply_mistnet(my_pvolfiles[1])
#  # mistnet will add additional parameters to the
#  # elevation scans at 0.5, 1.5, 2.5, 3.5 and 4.5 degrees
#  # let's extract the scan closest to 0.5 degrees:
#  my_scan <- get_scan(my_pvol, 0.5)
#  # plot some summary info about the scan to the console:
#  my_scan

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # let's add depolarization ratio (DR) as a parameter (following Kilambi 2018):
#  my_ppi <- calculate_param(my_ppi, DR = 10 * log10((1+ ZDR - 2 * (ZDR^0.5) * RHOHV) /
#    (1 + ZDR+ 2 * (ZDR^0.5) * RHOHV)))

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#    # plot the depolarization ratio, using a viridis color palette:
#  map(my_ppi, map = basemap, param = "DR", zlim=c(-25,-5), palette = viridis::viridis(100))

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS-------------------------------------
#  # Now let us screen out the reflectivity areas with DR < -12 dB:
#  my_ppi_clean <- calculate_param(my_ppi, DBZH = ifelse(DR < -12, NA, DBZH))
#  # plot the original and cleaned up reflectivity:
#  map(my_ppi, map = basemap, param = "DBZH", zlim = c(-20, 40))
#  map(my_ppi_clean, map = basemap, param = "DBZH", zlim = c(-20, 40))

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # as before, project the scan as a ppi:
#  my_ppi <- project_as_ppi(my_scan)
#  # plot the probability for the WEATHER class
#  plot(my_ppi, param = 'WEATHER')
#  # plot the final segmentation result:
#  # plot the probability for the WEATHER class
#  plot(my_ppi, param = 'CELL')
#  # let's remove the identified precipitation area (and additional fringe) from the ppi, and plot it:
#  my_ppi_clean <- calculate_param(my_ppi, DBZH = ifelse(CELL >= 1, NA, DBZH))
#  map(my_ppi_clean, map=basemap, param = 'DBZH')

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS-------------------------------------
#  # To remove the additional fringe around the rain segmentation by MistNet
#  # we want to screen out CELL values equal to 2 only.
#  my_ppi_clean <- calculate_param(my_ppi, VRADH = ifelse(CELL >= 2, NA, VRADH))
#  map(my_ppi_clean, map=basemap, param = 'VRADH')

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # Usually we would load processed vertical profiles (vp files) by:
#  # my_vplist <- read_vpfiles("./your/directory/with/processed/profiles/goes/here")
#  # my_vplist contains after running the command a list of vertical profile (vp) objects
#  # To save time, we load these data directly from file
#  my_vplist <- readRDS("data_vpts/KBRO20170514.rds")
#  # print the length of the vplist object. It should contain 95 profiles
#  length(my_vplist)

## ----eval=SHOW_ANSWERS--------------------------------------------------------
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

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # let's plot the vertical profile, in terms of bird density
#  plot(my_vp, quantity = "dens")
#  # print the currently assumed radar cross section (RCS) per bird:
#  rcs(my_vp)

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # Answer:
#  # All bird densities will be a factor 10 times lower.

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # let's change the RCS to 110 cm^2
#  rcs(my_vp) <- 110

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # After changing the RCS above, we simply plot the vertical profile again
#  plot(my_vp)
#  # Indeed the densities are scaled down by a factor 10

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # convert the list of vertical profiles into a time series:
#  my_vpts <- bind_into_vpts(my_vplist)
#  # time series objects can be subsetted, just as you may be used to with vectors
#  # here we subset the first 50 timesteps:
#  my_vpts[1:50]
#  # here we extract a single timestep, which gives you back a vertical profile class object:
#  my_vpts[100]
#  # to plot the full time series:
#  plot(my_vpts)
#  # check the help file for the plotting function of profile time series
#  # Because profile timeseries are of class 'vpts', it's associated plotting function
#  # is called plot.vpts:
#  ?plot.vpts

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # filter our vpts for night time
#  my_vpts_night <- filter_vpts(my_vpts, night=TRUE)
#  # plot this smaller time series:
#  plot(my_vpts_night)

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # Answer:
#  #
#  # At 1500 meter 6 UTC the wind barbs have 2 full flags and one half flag.
#  # Therefore the ground speed is approximately 2x5 + 2.5 = 12.5 m/s.

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # First extract the profile at 6 UTC:
#  vp_6UTC <- filter_vpts(my_vpts_night, nearest = "2017-05-14 06:00")
#  plot(vp_6UTC, quantity = "ff")

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # Plot the ground speed (ff):
#  plot(vp_6UTC, quantity = "ff")
#  # Speeds at 1500 metre is approximately 12 m/s, confirming our earlier reading above of 12.5 m/s.

## ----eval=SHOW_ANSWERS--------------------------------------------------------
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

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # Answer:
#  #
#  # VID = (200 birds / km^3) * (1 km) + (100 birds / km^3) * (0.5 km)
#  #     = 250 birds / km^2

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # Answer:
#  #
#  # MTR = (200 birds / km^3) * (50 km / hour) * (1 km) + (100 birds / km^3) * (100 km / hour) * (0.5 km)
#  #     = 10000 birds / km / hour + 5000 birds / km / hour
#  #     = 15000 birds / km / hour

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # Answer:
#  #
#  # MT = MTR * (3 hour) = (15000 birds / km / hour) * (3 hour)
#  #    = 45000 birds / km
#  # for a 10 kilometer transect:
#  # Number_of_birds = 10 km * 45000 birds / km = 450000 birds

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # instead of vertically integrated density (VID), you can use vertically integrated reflectivity (VIR):
#  plot(my_vpi, quantity = "vir")
#  # instead of migration traffic rate (MTR), you can use the reflectivity traffic rate (RTR):
#  plot(my_vpi, quantity = "rtr")
#  # instead of migration traffic (MT), you can use the reflectivity traffic (RT):
#  plot(my_vpi, quantity = "rt")

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # load a time series for the KBGM radar in Binghamton, NY
#  my_vpts <- readRDS("data_vpts/KBGM20170527-20170602.rds")
#  # print the loaded vpts time series for this radar:
#  my_vpts
#  # plot the bird density over time:
#  plot(my_vpts, quantity = "dens")

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # also show meteorological signals:
#  plot(my_vpts, quantity = "DBZH")
#  
#  # Periods with high reflectivities extending to high altitudes indicate precipitation,
#  # i.e. second half of the second night, and on and off during the fourth night.

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # define ranges from 0 to 2500000 meter (250 km), in steps of 100 m:
#  range <- seq(0, 250000, 100)
#  
#  # plot the beam height of the 0.5 degree elevation beam:
#  plot(range, beam_height(range, 0.5), ylab = "beam height [m]", xlab = "range [m]", type='l')
#  
#  # let's add the lower and upper altitude of the beam, as determined by the beam width:
#  points(range, beam_height(range, 0.5)-beam_width(range)/2, type='l',lty=3)
#  points(range, beam_height(range, 0.5)+beam_width(range)/2, type='l',lty=3)

## ----eval=SHOW_ANSWERS, warning=FALSE, results='hide'-------------------------
#  # download a polar volume for the KBRO radar in Brownsville, TX
#  download_pvolfiles(date_min=as.POSIXct("2017-05-14 05:50:00"), date_max=as.POSIXct("2017-05-14 06:00:00"), radar="KBRO", directory="./data_pvol")
#  # Load all the polar volume filenames downloaded so far for the KBRO radar:
#  my_pvolfiles <- list.files("./data_pvol", recursive = TRUE, full.names = TRUE, pattern="KBRO")
#  # we will process the first one into a vp:
#  my_pvolfile <- my_pvolfiles[1]
#  # calculate the profile, using MistNet to remove precipitation:
#  # we calculate 60 layers of 50 meter width, so up to 30*100=3000 m.
#  my_vp <- calculate_vp(my_pvolfile, n_layer=60, h_layer=50, sd_vvp_threshold = 1)

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS, warning=FALSE----------------------
#  # plot the profile
#  plot(my_vp)
#  
#  # Answer:
#  # * The vertical profile shows altitude bands of about 300 meter width
#  #   We therefore expect the radar to start having difficulty to resolve this profile
#  #   when the beam shape becomes broader, so less than 50 km
#  #   (this is why we typically use 35 km as the maximum range to use in vertical profile estimation.)
#  # * the birds fly up to about 2km in this profile.
#  #   At about 200000 m (200 km) the lower end of the beam no longer overlaps with the migration layer
#  #   Therefore, if we birds are flying according to the same altitude distribution everywhere,
#  #   we expect the radar to become fully blind for birds beyond 200 km

## ----eval=SHOW_ANSWERS,warning=FALSE------------------------------------------
#  # We will use the piping operator %>% of magrittr package to
#  # execute multiple operations in one statement:
#  library(magrittr)
#  # first, load the polar volume:
#  my_pvolfile %>% read_pvolfile() -> my_pvol
#  # Next, let's calculate a PPI for the 1.5 degree elevation scan
#  # Finally, we calculated the vertically integrated PPI
#  my_ppi_integrated <- integrate_to_ppi(pvol=my_pvol,vp=my_vp, res=1000)

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS,warning=FALSE-----------------------
#  #
#  my_pvol %>%
#    get_scan(elev = 1.0) %>%
#    project_as_ppi(raster=raster::raster(my_ppi_integrated$data)) ->
#    my_ppi
#  # let's take the logarithm of the vertically integrated density, so
#  # we can compare it more directly to DBZH, which is also a logarithmic quantity:
#  my_ppi_integrated <- calculate_param(my_ppi_integrated, logVID=log(VID))
#  # plot both ppi's:
#  plot(my_ppi)
#  plot(my_ppi_integrated, param="logVID", zlim=c(-10,10))

## ----echo=SHOW_ANSWERS, eval=SHOW_ANSWERS-------------------------------------
#  # The altitude distribution in this case shows two migration layers.
#  # In the 1.5 PPI scan these show up as two concentric rings of density.
#  # Since in a conventional PPI a larger distance from the radar also
#  # implies a higher altitude, we see the rings showing up at the
#  # ranges where the beam intersects each of the migration layers.
#  #
#  # The vertically integrated PPI is more straightforward to interpret
#  # Here we correct for aforementioned beam effects and integrate the
#  # density over altitude, giving a more realistic reconstruction of the
#  # spatial distribution of migrants.

## ----eval=SHOW_ANSWERS,results='hide',warning=FALSE---------------------------
#  # First we download more data, for a total of one additional hour for the same radar:
#  download_pvolfiles(date_min=as.POSIXct("2017-05-04 01:40:00"), date_max=as.POSIXct("2017-05-04 02:40:00"), radar="KHGX", directory="./data_pvol")
#  # We will process all the polar volume files downloaded so far:
#  my_pvolfiles <- list.files("./data_pvol", recursive = TRUE, full.names = TRUE, pattern="KHGX")
#  my_pvolfiles
#  # create output directory for processed profiles
#  outputdir <- "./data_vp"
#  dir.create(outputdir)

## ----eval=FALSE, echo=FALSE---------------------------------------------------
#  # let's loop over the files and generate profiles
#  start=Sys.time()
#  for (file_in in my_pvolfiles) {
#    # generate the output filename for the input file
#    file_out <- paste(outputdir, "/", basename(file_in), "_vp.h5", sep = "")
#    # generate the profile, and write it as a hdf5 file (for later reference)
#    # we enclose the calculate_vp statement in tryCatch to keep going after errors with specific files
#    vp <- tryCatch(calculate_vp(file_in, file_out, mistnet = TRUE, h_layer=50, n_layer=60), error = function(e) NULL)
#  }
#  end=Sys.time()
#  # calculate the processing time that has passed:
#  end-start

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # we assume outputdir contains the path to the directory with processed profiles
#  my_vpfiles <- list.files(outputdir, full.names = TRUE, pattern="KHGX")
#  # print them
#  my_vpfiles
#  # read them
#  my_vplist <- read_vpfiles(my_vpfiles)

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  # make a time series of profiles:
#  my_vpts <- bind_into_vpts(my_vplist)
#  # plot them between 0 - 3 km altitude:
#  plot(my_vpts, ylim = c(0, 3000))
#  # because of the rain moving in, our ability to estimate bird profiles slowly drops:
#  # let's visualize rain by plotting all aerial reflectivity:
#  plot(my_vpts, ylim = c(0, 3000), quantity="DBZH")
#  

## ----eval=SHOW_ANSWERS--------------------------------------------------------
#  process_file <- function(file_in){
#    # construct output filename from input filename
#    file_out <- paste(outputdir, "/", basename(file_in), "_vp.h5", sep = "")
#    # run calculate_vp()
#    vp <- tryCatch(calculate_vp(file_in, file_out, mistnet = FALSE, h_layer=50, n_layer=60), error = function(e) NULL)
#    if(is.vp(vp)){
#      # return TRUE if we calculated a profile
#      return(TRUE)
#    } else{
#      # return FALSE if we did not
#      return(FALSE)
#    }
#  }
#  
#  # To process a file, we can now run
#  process_file(my_pvolfiles[1])

## ----eval=FALSE---------------------------------------------------------------
#  # load the parallel library
#  library(parallel)
#  # detect how many cores we can use. We will keep 2 cores for other tasks and use the rest for processing.
#  number_of_cores = detectCores() - 2
#  # Available cores:
#  number_of_cores
#  # let's loop over the files and generate profiles
#  start=Sys.time()
#  mclapply(my_pvolfiles, process_file, mc.cores = number_of_cores)
#  end=Sys.time()
#  # calculate the processing time that has passed:
#  end-start

