% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_funs.R
\name{plot_many_country_years}
\alias{plot_many_country_years}
\title{Plot imprinting probabilities for up to five country-years}
\usage{
plot_many_country_years(imprinting_df)
}
\arguments{
\item{imprinting_df}{A long data frame of imprinted probabilities output by \code{\link[=get_imprinting_probabilities]{get_imprinting_probabilities()}}. Up to five countries and an arbitrary span of years can be plotted.}
}
\value{
No return value. Opens a plot of the data frame.
}
\description{
For each country and year, generate two plots:
\itemize{
\item A stacked barplot, where each bar represents a birth cohort, and
the colors within the bar show the probabilities that someone born
in that cohort has a particular imprinting status, for the first
observation year.
\item A lineplot showing the age-specific probability of imprinting to
H3N2 in the first and last observation year. When the data contain
more than one observation year, this plot shows how cohorts age over
time.
}
}
\examples{
imprinting_df <- get_imprinting_probabilities(
  observation_years = c(1920, 1921),
  countries = c("Oman", "Indonesia")
)
plot_many_country_years(imprinting_df)
}
