
test_that("Observation years greater than current year raises error.", {
  obs_year <- as.numeric(format(Sys.Date(), "%Y")) + 1
  expect_error(get_imprinting_probabilities(observation_years = obs_year, countries = c("United States")))
})

test_that("Probabilities are not returned for birth cohorts not yet alive in the year of observation.", {
  probs <- get_imprinting_probabilities(observation_years = c(2000, 2022), countries = "Brazil")
  expect_true(all(probs$year >= probs$birth_year))
})

test_that("Range of years returned equals range of years passed.", {
  obs_year <- 2018
  min_year <- 1918

  probs <- get_imprinting_probabilities(observation_years = obs_year, countries = c("United States"))

  expect_equal(max(probs$birth_year), obs_year)
  expect_equal(min(probs$birth_year), min_year)
})

test_that("Observation year prior to 1996 returns known, valid probabilities", {
  probs <- get_imprinting_probabilities(observation_years = 1919, countries = "Aruba", df_format = "long")

  expect_equal(probs$imprinting_prob, c(0.70, 0.91, 0.00, 0.00, 0.00, 0.00, 0.30, 0.09), tolerance = 0.01)
})


test_that("Known, valid imprinting probabilities are returned", {
  probs <- get_imprinting_probabilities(observation_years = 2020, countries = "Canada", df_format = "long") %>%
    dplyr::filter(subtype == "A/H1N1") %>%
    pull(imprinting_prob)

  expect_equal(probs[1:10], c(0.0617228, 0.2165782, 0.3095202, 0.1402262, 0.4345298, 0.2954388, 0.4644255, 0.4274897, 0.3563988, 0.2927145), tolerance = 0.01)
})
