INTENSITY_DATA <- readRDS(system.file("extdata", "INTENSITY_DATA.rds", package = "imprinting"))

#' Calculate the probability imprinting occurs n years after birth
#'
#' Given an individual's birth year, the year of observation, and pre-calculated influenza circulation intensities, calculate the probability that the first influenza infection occurs exactly 0, 1, 2, ... 12 years after birth.
#'
#' @details The probability of primary influenza infection n years after birth is calculated based on a modified [geometric distribution](https://en.wikipedia.org/wiki/Geometric_distribution): let p be the average annual probability of a primary influenza infection. Then the probability that primary infection occurs n=0,1,2,... years after birth is \eqn{p*(1-p)^{n}}.
#'
#' This function modifies the geometric model above to account for changes in annual circulation intensity, so that annual probabilities of primary infection \eqn{p_i} are scaled by the intensity in calendar year i. Details are given in \doi{https://doi.org/10.1126/science.aag1322}{Gostic et al. Science, (2016)}.
#'
#' @param birth_year year of birth (numeric). Must be between 1918 and the current calendar year.
#' @param observation_year year of observation, which affects the birth cohort's age.
#' @param intensity_df data frame of annual intensities, output by [get_country_intensity_data()].
#' @param max_year maximum year for which to output probabilities. Must be greater than or equal to observation_year. (If in doubt, set equal to observation year.)
#' @param baseline_annual_p_infection average annual probability of primary infection. The default, 0.28, was estimated using age-seroprevalence data in \doi{https://doi.org/10.1126/science.aag1322}{Gostic et al. Science, (2016)}.
#'
#' @return a vector whose entries show the probability that a person born
#' in year 0 was first infected by influenza in year 0, 1, 2, 3, ...12
#' We only consider the first 13 probabilities (i.e. we assume everyone
#' imprints before age 13. These outputs are not normalized, so the
#' vector sum asymptotically approaches one, but is not exactly equal
#'  to one. For cohorts born <13 years prior to the year of observation,
#'  the output vector will have less than 13 entries.
#'
#' @examples
#' # For a cohort under 12 years old and born in 2000, return the
#' # probabilities of primary infection in 2000, 2001, ... 2012:
#' get_p_infection_year(
#'   birth_year = 2000,
#'   observation_year = 2022,
#'   intensity_df = get_country_intensity_data("Canada", 2022),
#'   max_year = 2022
#' )
#'
#' # If the cohort is still under age 12 at the time of observation, return
#' # a truncated vector of probabilities:
#' get_p_infection_year(
#'   birth_year = 2020,
#'   observation_year = 2022,
#'   intensity_df = get_country_intensity_data("Mexico", 2022),
#'   max_year = 2022
#' )
#'
#' @export
get_p_infection_year <- function(birth_year,
                                 observation_year,
                                 intensity_df,
                                 max_year,
                                 baseline_annual_p_infection = 0.28) {
  ## INPUTS
  ##    - year in which an individual was born (birth.year)
  ##    - year in which the individual became infected with bird flu (infection year)
  ## OUTPUTS
  ##    - vector of 13 probabilities, the first representing the probability of first flu infection in the first year of life (age 0), the second representing the probability of first flu infection in the second year of life (age 1), and so on up to the 13th year of life (age 12)
  stopifnot(observation_year <= max_year)
  stopifnot(birth_year <= observation_year)
  # Weighted attack rate = annual prob infection weighted by circulation intensity
  weighted.attack.rate <- baseline_annual_p_infection * (intensity_df$intensity)
  names(weighted.attack.rate) <- intensity_df$year
  ################# Calculations ---------------
  possible_imprinting_years <- birth_year:min(birth_year + 12, observation_year) # Calendar years of first infection (ages 0-12)
  nn <- length(possible_imprinting_years) # How many possible years of first infection? (should be 13)
  valid_attack_rates <- weighted.attack.rate[as.character(possible_imprinting_years)] # Get weighted attack rates corresponding to possible years of first infection
  attack_rate_complements <- matrix(rep(1 - valid_attack_rates, nn), nn, nn, byrow = T)
  ## Create matrices of 0s and 1s, which will be used below to vectorize and speed calculations
  infection_year <- not_infection_year <- matrix(0, nn, nn)
  diag(infection_year) <- 1 # Fill in diagonal of one with 1s for years of first infection
  not_infection_year[lower.tri(not_infection_year)] <- 1 # Fill in sub-diagonal for all the years since birth in which the individual escaped infection.
  # Exact probability of escaping infection in the previous (x-1) years, and becoming infected in year x
  prod.mat <- (valid_attack_rates * infection_year) + (attack_rate_complements * not_infection_year)
  # Fill in upper triangle with 1s to make multiplication possible
  prod.mat[upper.tri(prod.mat)] <- 1
  # Take product across rows
  p_ij <- apply(prod.mat, 1, prod)
  p_ij # Output probability of first infection in year i given birth year
}


to_long_df <- function(outlist) {
  # bind column name variables to function to avoid nonstandard evaluation issues in CRAN
  year_country <- year <- country <- birth_year <- NULL

  ## Reformat the list of matrix outputs into a long data frame
  reformat_one_list_element <- function(ll) {
    ## ll is a matrix whose columns represent birth years, and rows represent unique countries and years of observation
    mat_rownames <- rownames(ll) ## Extract the country-year rownames
    as_tibble(ll) %>% ## cast to tibble
      mutate(year_country = mat_rownames) %>% ## Make the country-year rownames into a column
      extract(year_country, into = c("year", "country"), regex = "(\\d{4})(\\w.+)", convert = T) %>%
      pivot_longer(-c(year, country), values_to = "imprinting_prob", names_to = "birth_year") %>%
      mutate(birth_year = as.integer(birth_year))
  }
  ## Apply the reformatting function to all list elements
  list_names <- names(outlist)
  subtypes <- gsub("(.+)_probs", "\\1", list_names) ## extract the imprinted subtype from list names
  list_of_dfs <- lapply(outlist, reformat_one_list_element) ## reformat each matrix into a data frame
  names(list_of_dfs) <- subtypes ## Get the list of the subtype represented by each matrix in the list
  bind_rows(list_of_dfs, .id = "subtype") ## Bind all subtypes into a single long data frame and return
}


#' Calculate imprinting probabilities
#'
#' For each country and year of observation, calculate the probability that cohorts born in each year from 1918 through the year of observation imprinted to a specific influenza A virus subtype (H1N1, H2N2, or H3N2), or group (group 1 contains H1N1 and H2N2; group 2 contains H3N2).
#'
#' @param observation_years year(s) of observation in which to output imprinting probabilities. The observation year, together with the birth year, determines the birth cohort's age when calculating imprinting probabilities. Cohorts <=12 years old at the time of observation have some probability of being naive to influenza.
#' @param countries a vector of countries for which to calculate imprinting probabilities. Run `show_available_countries()` for a list of valid inputs, and proper spellings.
#' @param annual_frequencies an optional input allowing users to specify custom circulation frequencies for arbitrary types of imprinting in order to study, e.g. imprinting to specific strains, clades, or imprinting by vaccination. If nothing is input, the default is to calculate subtype-specific probabilities (possible imprinting types are A/H1N1, A/H2N2, A/H3N2, or naive). See Details.
#' @param df_format must be either 'long' (default) or 'wide'. Controls whether the output data frame is in long format (with a single column for calculated probabilities and a second column for imprinting subtype), or wide format (with four columns, H1N1, H2N2, H3N2, and naive) showing the probability of each imprinting status.
#'
#' @details Imprinting probabilities are calculated following \doi{https://doi.org/10.1126/science.aag1322}{Gostic et al. Science, (2016)}. Briefly, the model first calculates the probability that an individual's first influenza infection occurs 0, 1, 2, ... 12 years after birth using a modified geometric waiting time model. The annual circulation intensities output by [get_country_intensity_data()] scale the probability of primary infection in each calendar year.
#'
#' Then, after calculating the probability of imprinting 0, 1, 2, ... calendar years after birth, the model uses data on which subtypes circulated in each calendar year (from [get_country_cocirculation_data()]) to estimate that probability that a first infection was caused by each subtype. See [get_country_cocirculation_data()] for details about the underlying data sources.
#'
#' To calculate other kinds of imprinting probabilities (e.g. for specific clades, strains, or to include pediatric vaccination), users can specify custom circulation frequencies as a list, `annual_frequencies`. This list must contain one named element for each country in the `countries` input vector. Each list element must be a data frame or tibble whose first column is named "year" and contains numeric years from 1918:max(`observation_years`). Columns 2:N of the data frame must contain circulation frequencies that sum to 1 across each row, and each column must have a unique name indicating the exposure kind. E.g. column names could be {"year", "H1N1", "H2N2", "H3N2", "vaccinated"} to include probabilities of imprinting by vaccine, or {"year", "3C.3A", "not_3C.3A"} to calculate clade-specific probabilities.  Do not include a naive column. Any number of imprinting types is allowed, but the code is not optimized to run efficiently when the number of categories is very large. Frequencies within the column must be supplied by the user. See [Vieira et al. 2021](https://www.nature.com/articles/s41467-021-24566-y) for methods to estimate circulation frequencies from sequence databases like [GISAID](https://gisaid.org/) or the [NCBI Sequence Database](https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database).
#'
#' See `vignette("custom-imprinting-types")` for use of a custom `annual_frequencies` input.
#'
#' @return
#' * If `format=long` (the default), a long tibble with columns showing the imprinting subtype (H1N1, H2N2, H3N2, or naive), the year of observation, the country, the birth year, and the imprinting probability.
#' * If `format=wide`, a wide tibble with each row representing a country, observation year, and birth year, and with a column for each influenza A subtype (H1N1, H2N2, and H3N2), or the probability that someone born in that year remains naive to influenza and has not yet imprinted.
#' For cohorts >12 years old in the year of observation, the probability of remaining naive is 0, and the subtype-specific probabilities are normalized to sum to 1. For cohorts <=12 years old in the year of observation, the probability of remaining naive is non-zero. For cohorts not yet born at the time of observation, all output probabilities are 0.
#'
#' @examples
#' # ===========================================================
#' # Get imprinting probabilities for one country and year
#' get_imprinting_probabilities(2022, "United States")
#' # ===========================================================
#' # Return the same outputs in wide format
#' get_imprinting_probabilities(2022,
#'   "United States",
#'   df_format = "wide"
#' )
#'
#' @export
get_imprinting_probabilities <- function(observation_years,
                                         countries,
                                         annual_frequencies = NULL,
                                         df_format = "long") {
  # bind column name variables to function to avoid nonstandard evaluation issues in CRAN
  country <- year <- birth_year <- subtype <- NULL

  ## Input checks
  current_year <- as.numeric(format(Sys.Date(), "%Y"))
  if (!all(observation_years >= 1918 & observation_years <= current_year)) {
    stop("observation_years must be a numeric vector with values between 1918 and the current calendar year.")
  }
  if (!all(countries %in% pull(show_available_countries(), country))) {
    problem_inputs <- countries[!(countries %in% pull(show_available_countries()))]
    stop(sprintf("You input the following country names, which are invalid: \n\n%s\n\nRun `show_available_countries()` to see a list of valid countries.", paste(problem_inputs, collapse = ", ")))
  }
  max_year <- max(observation_years)
  stopifnot(max_year <= as.numeric(format(Sys.Date(), "%Y")))
  birth_years <- 1918:max_year
  infection_years <- birth_years
  nn_birth_years <- length(birth_years)
  ## annual frequencies
  if (length(annual_frequencies) == 0) { # If no custom inputs given
    annual_frequencies <- lapply(countries, function(cc) {
      get_country_cocirculation_data(cc, max_year) %>%
        select(1:4)
    })
    names(annual_frequencies) <- countries
  }
  # ## Check annual frequencies list
  # Must be a list
  if (!is.list(annual_frequencies)) {
    stop("annual_frequencies must be a named list. See Details in ?get_imprinting_probabilities.")
  }
  # List must be named, names must match `countries`
  if (!all(countries %in% names(annual_frequencies))) {
    stop("annual_frequencies must be a named list whose names match the countries input vector. See Details in ?get_imprinting_probabilities.")
  }
  # All data frames within the list must contain the year column
  if (!all(sapply(annual_frequencies, function(ll) {
    # All data frames contain a column, year
    names(ll)[1] == "year"
  }))) {
    stop("The first column of all data frames in annual_frequencies must be named `year`.")
  }
  if (!all(sapply(annual_frequencies, function(ll) {
    # All year columns contain years 1918:max_year
    all(1918:max_year %in% ll$year)
  }))) {
    stop("The `year` column of all data frames in annual_frequencies must contain numeric values from 1918:max(observation_years)")
  }
  if (!all(sapply(annual_frequencies, function(ll) {
    # All type frequencies sum to 1 within each row (year)
    all(abs(rowSums(ll[, -1]) - 1) <= 1e-7)
  }))) {
    stop("Each row in the annual_frequencies data frames (not including the year column) must sum to 1.")
  }

  ## For each country, get imprinting probabilities
  # for (this_country in countries) {
  imprinting_probs <- lapply(countries, function(this_country) {
    who_region <- get_WHO_region(this_country)
    # get country-specific circulation intensity
    this_intensity_data <- get_country_intensity_data(this_country, max_year, min_specimens = 50)
    stopifnot(!any(is.na(this_intensity_data$intensity)))
    # get circulation fractions
    these_annual_frequencies <- annual_frequencies[[this_country]]
    # Calculate country-specific imprinting probs
    lapply(1:length(observation_years), function(jj) {
      ## Loop across birth years
      lapply(1918:observation_years[jj], FUN = function(bb) {
        get_probs_one_birth_year(
          this_birth_year = bb,
          this_observation_year = observation_years[jj],
          max_year = max_year,
          this_intensity_data = this_intensity_data,
          these_annual_frequencies = these_annual_frequencies
        )
      }) %>%
        bind_rows()
    }) %>% # end loop across observation years
      bind_rows() %>%
      mutate(country = this_country)
  }) %>% # end loop across countries
    bind_rows() %>%
    select(year, country, birth_year, !c(year, country, birth_year))


  if (df_format == "wide") {
    return(imprinting_probs)
  } else {
    stopifnot(df_format == "long")
    return(imprinting_probs %>%
      pivot_longer(-c(1:3),
        names_to = "subtype",
        values_to = "imprinting_prob"
      ) %>%
      arrange(subtype, dplyr::desc(birth_year)))
  }
}


get_probs_one_birth_year <- function(this_birth_year,
                                     this_observation_year,
                                     max_year,
                                     this_intensity_data,
                                     these_annual_frequencies) {
  # bind column name variables to function to avoid nonstandard evaluation issues in CRAN
  year <- NULL
  ## If the cohort has not yet been born, return NA
  if (this_birth_year > this_observation_year) {
    imprinting_probs <- rep(NA, ncol(these_annual_frequencies))
    names(imprinting_probs) <- c(names(these_annual_frequencies)[-1], "naive")
  } else {
    ## Loop across birth years
    # get possible years of first infection for this birth year
    # first infections can occur up to age 11, or up until the current year, whichever comes first
    n_infection_years <- min(12, this_observation_year - this_birth_year)
    valid_infection_years <- this_birth_year + (0:n_infection_years)
    # get year-specific probabilities of primary infection
    inf.probs <- get_p_infection_year(
      birth_year = this_birth_year,
      observation_year = this_observation_year,
      baseline_annual_p_infection = 0.28,
      max_year = max_year,
      intensity_df = this_intensity_data
    )
    # If all 13 possible years of infection have passed, normalize so that the probability of imprinting from age 0-12 sums to 1
    if (length(inf.probs) == 13) inf.probs <- inf.probs / sum(inf.probs)
    # Else, don't normalize and extract the probability of remaiing naive below.
    # Combine primary infection probabilities with data on what strains circulated each year
    freq_mat <- these_annual_frequencies %>%
      dplyr::filter(year %in% valid_infection_years) %>% # pull out relevant years of primary infection
      select(-1) %>%
      as.matrix()
    imprinting_probs <- colSums(freq_mat * inf.probs) # Get type-specific probabilities
    imprinting_probs <- c(imprinting_probs, "naive" = 1 - sum(imprinting_probs)) # Add the naive probability
  }
  # Return
  c(
    year = this_observation_year,
    birth_year = this_birth_year,
    imprinting_probs
  )
}
