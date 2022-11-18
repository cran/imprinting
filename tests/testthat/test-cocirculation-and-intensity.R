
test_that("Observation year less than 13 years from birth year gives known-good probabilities.", {
  this_epi_data <- get_country_cocirculation_data("United States", 2010)
  INTENSITY_DATA <- get_country_intensity_data("United States", 2010, min_specimens = 50)

  # Birth Year 2007
  expect_equal(get_p_infection_year(2007, 2010, INTENSITY_DATA, 2010, 0.28), c(0.20244365, 0.18857213, 0.42628896, 0.03914203), tolerance = 0.000001)
  # Birth Year 2010
  expect_equal(get_p_infection_year(2010, 2010, INTENSITY_DATA, 2010, 0.28), c(0.2142477), tolerance = 0.000001)
})


test_that("Observation year greater than 13 years from birth year gives known-good probabilities.", {
  this_epi_data <- get_country_cocirculation_data("United States", 2010)
  INTENSITY_DATA <- get_country_intensity_data("United States", 2010, min_specimens = 50)

  # Birth Year 1995
  expect_equal(get_p_infection_year(1995, 2010, INTENSITY_DATA, 2010, 0.28), c(0.22467298, 0.24317721, 0.07694539, 0.03014844, 0.11330399, 0.05630441, 0.03586139, 0.04841809, 0.06950272, 0.01227025, 0.02001303, 0.01554628, 0.01089872), tolerance = 0.000001)

  # Birth Year 1998
  expect_equal(get_p_infection_year(1998, 2010, INTENSITY_DATA, 2010, 0.28), c(0.066230545, 0.248907939, 0.123690388, 0.078780853, 0.106365603, 0.152684635, 0.026955481, 0.043964926, 0.034152297, 0.023942471, 0.022301922, 0.050416058, 0.004629224), tolerance = 0.000001)
})

test_that("Number of probabilities are same as number of years in input.", {
  this_epi_data <- get_country_cocirculation_data("United States", 2010)
  INTENSITY_DATA <- get_country_intensity_data("United States", 2010, min_specimens = 50)

  # Long span
  birth_year <- 1998
  obs_year <- 2010

  expect_equal(length(get_p_infection_year(birth_year, obs_year, INTENSITY_DATA, obs_year, 0.28)), obs_year - birth_year + 1)

  # Short span
  birth_year <- 2000
  obs_year <- 2003

  expect_equal(length(get_p_infection_year(birth_year, obs_year, INTENSITY_DATA, obs_year, 0.28)), obs_year - birth_year + 1)

  # One year
  birth_year <- 2008
  obs_year <- 2008

  expect_equal(length(get_p_infection_year(birth_year, obs_year, INTENSITY_DATA, obs_year, 0.28)), obs_year - birth_year + 1)
})


test_that("Misspelled country name gives an error.", {
  expect_error(get_country_cocirculation_data("United States?", 2010))
  expect_error(get_country_intensity_data("UnitedStates", 2010, min_specimens = 50))
  expect_error(get_imprinting_probabilities(observation_years = obs_year, countries = c("UnitedStates")))
})

test_that("Invalid year gives an error.", {
  expect_error(get_country_cocirculation_data("United States", -1))
  expect_error(get_country_intensity_data("United States", 1900, min_specimens = 50))
})

test_that("Numeric (non-NA) probabilities are returned for post-2017 observation years.", {
  obs_year <- 2022
  min_year <- 1918

  probs <- get_imprinting_probabilities(observation_years = obs_year, countries = c("United States"))

  expect_false(any(is.na(probs$imprinting_prob)))
})

test_that("Countries with low-quality intensity data return appropriate intensity values.", {
  intensities_Germany <- get_country_intensity_data(
    country = c("Germany"),
    max_year = 2022,
    min_specimens = 50
  )
  intensities_Iraq <- get_country_intensity_data(
    country = c("Iraq"),
    max_year = 2022,
    min_specimens = 50
  )

  expect_false(any(is.na(intensities_Germany$intensity)))
  expect_false(any(is.na(intensities_Iraq$intensity)))
  expect_true(all(intensities_Germany$intensity >= -2.5 & intensities_Germany$intensity <= 2.5))
  expect_true(all(intensities_Iraq$intensity >= -2.5 & intensities_Iraq$intensity <= 2.5))
})


test_that("Cocirculation returns the type from the output format parameter.", {
  cocirc_matrix <- get_country_cocirculation_data("United States", "2019", output_format = "matrix")
  cocirc_tibble <- get_country_cocirculation_data("United States", "2019", output_format = "tibble")

  expect_true(all(class(cocirc_matrix) == c("matrix", "array")))
  expect_true(all(class(cocirc_tibble) == c("tbl_df", "tbl", "data.frame")))
})
