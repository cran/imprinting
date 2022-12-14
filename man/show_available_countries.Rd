% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_import_funs.R
\name{show_available_countries}
\alias{show_available_countries}
\title{Show a list of all available countries}
\usage{
show_available_countries()
}
\value{
A data frame of valid country names.
}
\description{
Lists all available countries, with valid spelling and formatting. Each country in the list matches or can be mapped to a country with data in \href{https://apps.who.int/flumart/Default?ReportNo=12}{WHO Flu Mart}. (Note: for convenience, this package sometimes uses different country names or spellings than Flu Mart.)
}
\examples{
show_available_countries()
}
