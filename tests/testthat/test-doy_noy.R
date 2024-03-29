vp <- example_vp
vpts <- example_vpts
vpi <- integrate_profile(example_vpts)
pvolfile <- system.file("extdata", "volume.h5", package = "bioRad")
pvol <- read_pvolfile(pvolfile)

test_that("doy_noy() returns error on incorrect parameters", {
  expect_error(doy("not_a_object"),
               regexp = 'argument "lat" is missing, with no default',
               fixed = TRUE)
  expect_error(doy("not_a_object", lat = 23),
               regexp = 'argument "lon" is missing, with no default',
               fixed = TRUE)
  expect_error(doy("not_a_object", lon = 23,lat = 4000),
               regexp = "invalid coordinates",
               fixed = TRUE)
  expect_error(doy("not_a_object", lon = 23,lat = 40),
               regexp = "character string is not in a standard unambiguous format",
               fixed = TRUE)
  expect_error(noy("not_a_object"),
               regexp = 'argument "lat" is missing, with no default',
               fixed = TRUE)
  expect_error(noy("not_a_object", lat = 23),
               regexp = 'argument "lon" is missing, with no default',
               fixed = TRUE)
  expect_error(noy("not_a_object", lon = 23,lat = 4000),
               regexp = "invalid coordinates",
               fixed = TRUE)
  expect_error(noy("not_a_object", lon = 23,lat = 40),
               regexp = "non-numeric argument to binary operator",
               fixed = TRUE)
})

test_that("doy_noy() returns a number", {
  expect_type(doy(as.POSIXct("2015-01-01 18:00:00", tz="UTC"), 12.8517, 56.3675), "integer")
  expect_type(doy(vp), "integer")
  expect_type(doy(vpts), "integer")
  expect_type(doy(vpi), "integer")
  expect_type(doy(pvol), "integer")

  expect_type(noy(as.POSIXct("2015-01-01 18:00:00", tz="UTC"), 12.8517, 56.3675), "integer")
  expect_type(noy(vp), "integer")
  expect_type(noy(vpts), "integer")
  expect_type(noy(vpi), "integer")
  expect_type(noy(pvol), "integer")
})

test_that("doy_noy() returns correct doy and noy", {
  expect_equal(doy(as.POSIXct("2015-01-01 18:00:00", tz="UTC"), 12.8517, 56.3675), 1)
  # noy starts with New Years eve night as first night
  expect_equal(noy(as.POSIXct("2015-01-01 00:00:00", tz="UTC"), 12.8517, 56.3675), 1)
  expect_equal(noy(as.POSIXct("2015-01-01 18:00:00", tz="UTC"), 12.8517, 56.3675), 2)
})
