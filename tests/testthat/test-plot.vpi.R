example_vpi <- integrate_profile(example_vpts)

test_that("plot.vpi() returns error on incorrect parameters", {
  # use plot.vpts() to avoid defaulting to base plot()
  expect_error(plot.vpi("not_a_vpi"),
               regexp = 'inherits(x, "vpi") is not TRUE',
               fixed = TRUE)
  expect_error(plot(example_vpi, quantity = "not_a_quantity"),
               regexp = "quantity `not_a_quantity` not found in vpi object.",
               fixed = TRUE)

  # Test error on "param" instead of "quantity"
  expect_error(plot(example_vpi, param = "ff"),
               regexp = "unknown function argument 'param`. Did you mean `quantity`?",
               fixed = TRUE)
  # Return error when lon or lat is not an numeric
  expect_error(
    plot.vpi(example_vpi, lon = 'a'),
    regexp = "No latitude/longitude found in attribute data, please provide lat and lon arguments when night_shade=TRUE.",
    fixed = TRUE
  )
  expect_error(
    plot.vpi(example_vpi, lat = 'a'),
    regexp = "No latitude/longitude found in attribute data, please provide lat and lon arguments when night_shade=TRUE.",
    fixed = TRUE
  )
  expect_error(
    plot.vpi(example_vpi, lon = NA),
    regexp = "No latitude/longitude found in attribute data, please provide lat and lon arguments when night_shade=TRUE.",
    fixed = TRUE
  )
})

test_that("plot.vpi() produces plots", {
  expect_s3_class(recordPlot(plot(example_vpi)), "recordedplot")
  expect_s3_class(recordPlot(plot(example_vpi, quantity = "vir")), "recordedplot")
  expect_s3_class(recordPlot(plot(example_vpi, quantity = "mtr")), "recordedplot")
  expect_s3_class(recordPlot(plot(example_vpi, quantity = "dd")), "recordedplot")
})
