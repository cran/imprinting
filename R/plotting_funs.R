## Functions to plot the outputs of `get_imprinting_probabilities`

#' Plot imprinting probabilities for a single country and year
#'
#' Generate a stacked barplot, where each bar represents a birth cohort, and the colors within the bar show the probabilities that someone born in that cohort has a particular imprinting status. If the data frame contains more than one country or observation year, the first country-year is plotted by default. Specify other countries and years using the country and year options.
#'
#' @param imprinting_df A long data frame of imprinted probabilities output by [get_imprinting_probabilities()]. If the data frame contains more than one country and year, on the first will be plotted.
#' @param country An optional country name to plot. The input country name must exist in the imprinting_df.
#' @param year Similar to country, and optional input specifying the year for which to plot.
#'
#' @return No return value. Opens a plot of the data frame.
#' @examples
#' # Generate imprinting probabilities for one country and year
#' imprinting_df <- get_imprinting_probabilities(
#'   observation_years = 1920,
#'   countries = "Aruba"
#' )
#' plot_one_country_year(imprinting_df)
#'
#' # If we generate probabilities for more than one country and year,
#' imprinting_df <- get_imprinting_probabilities(
#'   observation_years = c(1922, 1925),
#'   countries = c(
#'     "Algeria",
#'     "South Africa"
#'   )
#' )
#'
#' # The default is to plot the first country year in the outputs
#' plot_one_country_year(imprinting_df)
#'
#' # Or, specify a country and year of interest (both must exist in the
#' # imprinting_df).
#' plot_one_country_year(imprinting_df,
#'   country = "South Africa",
#'   year = 1925
#' )
#'
#' @export
plot_one_country_year <- function(imprinting_df,
                                  country = NULL,
                                  year = NULL) {
  # bind column name variables to function to avoid nonstandard evaluation issues in CRAN
  birth_year <- imprinting_prob <- subtype <- NULL

  # if countries and years are specificied, check that they exist in the data frame
  if (length(country) > 0) {
    stopifnot(country %in% imprinting_df$country)
    stopifnot(length(country) == 1)
  }
  if (length(year) > 0) {
    stopifnot(year %in% imprinting_df$year)
    stopifnot(length(year) == 1)
  }

  if (ncol(imprinting_df) > 5) {
    stop("imprinting_df must be in long format, output from get_imprinting_probabilities()")
  }
  ## This function plots imprinting patterns for a single country-year
  ## If the data frame contains more than one country-year, it plots the first listed
  countries <- unique(imprinting_df$country)
  years <- unique(imprinting_df$year)
  obs_year <- ifelse(length(year) == 0, years[1], year)
  this_country <- ifelse(length(country) == 0, countries[1], country)
  axis_ticks <- seq(1920, obs_year, by = 10)
  replace_these <- which(axis_ticks %in% c(1960, 1970, 1980))
  axis_ticks[replace_these] <- c(1957, 1968, 1977) ## Replace 3 axis ticks with pandemic years
  axis_tick_labs <- sapply(axis_ticks, function(yr) {
    sprintf("%i\n%i", yr, obs_year - yr)
  }) # Label each tick with the birth year/current age of that cohort
  x_axis_text <- sprintf("Birth year\nAge in %s", obs_year) ## Axis label is birth year[newline]age now
  colors <- c("dodgerblue1", "lightblue", "firebrick2", "gray")
  imprinting_df %>%
    dplyr::filter(country == this_country) %>%
    dplyr::filter(year == obs_year) %>%
    dplyr::filter(birth_year <= obs_year) %>%
    ggplot() +
    geom_bar(aes(x = birth_year, y = imprinting_prob, fill = subtype), stat = "identity", width = 1) +
    xlab(x_axis_text) +
    ylab("Imprinting fraction") +
    scale_color_manual(values = colors, aesthetics = c("color", "fill")) +
    scale_x_continuous(breaks = axis_ticks, labels = axis_tick_labs) +
    ggtitle(sprintf("Probabilities for %s in %i", this_country, obs_year))
}


#' Plot imprinting probabilities for up to five country-years
#'
#' For each country and year, generate two plots:
#' * A stacked barplot, where each bar represents a birth cohort, and
#' the colors within the bar show the probabilities that someone born
#' in that cohort has a particular imprinting status, for the first
#' observation year.
#' * A lineplot showing the age-specific probability of imprinting to
#' H3N2 in the first and last observation year. When the data contain
#' more than one observation year, this plot shows how cohorts age over
#' time.
#'
#' @param imprinting_df A long data frame of imprinted probabilities output by [get_imprinting_probabilities()]. Up to five countries and an arbitrary span of years can be plotted.
#' @return No return value. Opens a plot of the data frame.
#'
#' @examples
#' imprinting_df <- get_imprinting_probabilities(
#'   observation_years = c(1920, 1921),
#'   countries = c("Oman", "Indonesia")
#' )
#' plot_many_country_years(imprinting_df)
#' @export
plot_many_country_years <- function(imprinting_df) {
  # bind column name variables to function to avoid nonstandard evaluation issues in CRAN
  country <- plot_number <- year <- birth_year <- imprinting_prob <- subtype <- age_at_observation <- pandemic_1968 <- NULL


  if (ncol(imprinting_df) > 5) {
    stop("imprinting_df must be in long format, output from get_imprinting_probabilities()")
  }
  countries <- unique(imprinting_df$country)
  if (length(countries) > 5) {
    warning("Plotting only the first 5 countires.")
  }
  years <- unique(imprinting_df$year)
  max_obs_year <- max(years)
  min_obs_year <- min(years)
  axis_ticks <- seq(1920, max_obs_year, by = 10)
  # replace_these <- which(axis_ticks %in% c(1960, 1970, 1980))
  # axis_ticks[replace_these] <- c(1957, 1968, 1977) ## Replace 3 axis ticks with pandemic years
  axis_tick_labs <- sapply(axis_ticks, function(yr) {
    sprintf("(%i) %i", max_obs_year - yr, yr)
  }) # Label each tick with the birth year/current age of that cohort
  x_axis_text <- sprintf("(age in %s) birth year", max_obs_year) ## Axis label is birth year[newline]age now
  colors <- c("dodgerblue1", "lightblue", "firebrick2", "gray")
  ## Print a max of 5 countries
  imprinting_df <- mutate(imprinting_df, plot_number = ceiling(as.numeric(as.factor(country)) / 5))
  ## For each country, plot barplots in the max year of observation
  bar_plots <- imprinting_df %>%
    dplyr::filter(plot_number == 1) %>% ## Only include up to 5 countries
    dplyr::filter(year == max_obs_year) %>%
    ggplot() +
    geom_bar(aes(x = birth_year, y = imprinting_prob, fill = subtype, color = subtype), stat = "identity", width = 1, linetype = 1) +
    xlab(x_axis_text) +
    ylab("imprinting fraction") +
    scale_color_manual(values = colors, aesthetics = c("color", "fill")) +
    scale_x_continuous(breaks = axis_ticks, labels = axis_tick_labs) +
    ggtitle(sprintf("Probabilities in %i", max_obs_year)) +
    facet_grid(country ~ .) +
    theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2))
  ## For each country, plot lineplots of the H1N1 imprinting probs, in the min and max years of observation
  line_plots <- imprinting_df %>%
    dplyr::filter(plot_number == 1) %>%
    dplyr::filter(year %in% c(max_obs_year, min_obs_year)) %>%
    mutate(
      age_at_observation = year - birth_year,
      pandemic_1968 = year - 1968
    ) %>%
    dplyr::filter(age_at_observation >= 0) %>%
    ggplot() +
    geom_vline(aes(xintercept = pandemic_1968, lty = as.factor(year)), show.legend = F) +
    geom_line(aes(x = age_at_observation, y = imprinting_prob, lty = as.factor(year)), color = "firebrick2") +
    xlab("Age at time of observation") +
    ylab("H3N2\nImprinting fraction") +
    geom_label(aes(x = min_obs_year - 1968, y = 0.1, label = "born 1968"), angle = 90, size = 10 / .pt, hjust = .5) +
    scale_linetype(name = "Year of observation\n  \n \n \n ") +
    ggtitle("Aging of cohort") +
    facet_grid(country ~ .) +
    theme(legend.position = "bottom") +
    guides(lty = guide_legend(nrow = 2))
  if (length(years) > 1) {
    line_plots <- line_plots + geom_segment(aes(x = min_obs_year - 1968, xend = max_obs_year - 1968, y = .25, yend = 0.25), arrow = arrow(length = unit(.05, "in")))
  }
  ## return an image with barplots on the left and lineplots on the right
  cowplot::plot_grid(bar_plots, line_plots, ncol = 2)
}
