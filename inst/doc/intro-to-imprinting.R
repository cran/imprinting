## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)
library(dplyr)
library(ggplot2)
library(tidyr)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("devtools", repos='http://cran.us.r-project.org')

## -----------------------------------------------------------------------------
library(devtools) # Load the package

## ----eval=FALSE---------------------------------------------------------------
#  devtools::install_github("cobeylab/imprinting")

## -----------------------------------------------------------------------------
library(imprinting) 

## ----include=FALSE------------------------------------------------------------
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
get_imprinting_probabilities(observation_years = 2022, countries = "United States")

## -----------------------------------------------------------------------------
get_imprinting_probabilities(observation_years = 2022, 
                             countries = "United States", 
                             df_format = 'wide')

## -----------------------------------------------------------------------------
get_imprinting_probabilities(observation_years = c(2005, 2011, 2012, 2022), 
                             countries = "United States", 
                             df_format = 'wide') %>%
  dplyr::filter(birth_year == 2000) %>%
  mutate(age_at_observation = year-birth_year) %>%
  select(c(1,2,3,8,4:7))

## -----------------------------------------------------------------------------
show_available_countries() %>%
  print(n = 200)

## -----------------------------------------------------------------------------
many_probabilities = get_imprinting_probabilities(observation_years = c(2000, 2019, 2022),
                                                  countries = c('Brazil', 'Afghanistan', 'Estonia', 'Finland')) 
## Store the outputs in a variable called many_probabilities
many_probabilities ## View the outputs in the console

## ----eval=FALSE---------------------------------------------------------------
#  # View the outputs in a separate window.
#  View(many_probabilities)
#  
#  # Save the outputs as a .csv file in your current working directory.
#  write_csv(many_probabilities, 'many_probabilities.csv')

## ----fig.width = 7------------------------------------------------------------
head(many_probabilities)
plot_one_country_year(many_probabilities) 

## ----fig.width = 7------------------------------------------------------------
plot_one_country_year(many_probabilities %>% 
                        dplyr::filter(country == 'Estonia', year == 2019))

## ----fig.width = 7, fig.height = 5--------------------------------------------
plot_many_country_years(many_probabilities)

## -----------------------------------------------------------------------------
get_country_cocirculation_data('United States', 2022)

## -----------------------------------------------------------------------------
get_country_intensity_data(country = 'China', max_year = 2022)

## -----------------------------------------------------------------------------
probs = get_p_infection_year(birth_year = 2000,
                     observation_year = 2022,
                     intensity_df = get_country_intensity_data('Mexico', 2022),
                     max_year = 2022)
names(probs) = as.character(2000+(0:12))
probs
sum(probs) ## Raw probabilities are not yet normalized.

norm_probs = probs/sum(probs) ## Normalize
sum(norm_probs)

