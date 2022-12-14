% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_import_funs.R
\name{get_template_data}
\alias{get_template_data}
\title{Get country-independent flu circulation data for 1918-1996}
\usage{
get_template_data()
}
\value{
A tibble with the following columns:
\itemize{
\item year
\item \code{A/H1N1}, \code{A/H2N2}, and \code{A/H3N2} show the fraction of influenza cases caused by each subtype.
\item \code{A} = \code{A/H1N1} + \code{A/H2N2} + \code{A/H3N2}
\item \code{B} is a placeholder for future calculate of influenza B imprinting probabilities, which currently contains \code{NA}.
\item \code{group1} and \code{group2} show the fraction of cases caused by group 1 subtypes (H1N1 and H2N2), or group 2 (H3N2).
\item \code{data_from} notes the data source.
}
}
\description{
\code{get_template_data()} returns a tibble showing the fraction of influenza A cases caused by subtype H1N1, H2N2, or H3N2 in each year from 1918-1996. These data are country-independent. Country-specific data are only available from 1997 on.
\itemize{
\item For years 1918-1976 only one influenza A subtype circulated, so all fractions are 0 or 1.
\item From 1977-1996 H1N1 and H3N2 both circulated. \code{get_template_data()} reports the fraction of influenza A-positive specimens of each subtype observed in US flu surveillance. See \href{https://jamanetwork.com/journals/jama/fullarticle/195750}{Thompson et al. JAMA, 2003}, Table 1.
\item Country-specific data from \href{https://apps.who.int/flumart/Default?ReportNo=12}{WHO Flu Mart} will be appended to this template in later steps.
}
}
\seealso{
\doi{https://doi.org/10.1126/science.aag1322}{Gostic et al. Science, (2016)} for detailed methods.
}
