## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)
library(dplyr)

## -----------------------------------------------------------------------------
library(imprinting)
## Start with subtype-specific fractions for H1N1, H2N2, H3N2
US_frequencies = get_country_cocirculation_data(country = 'United States', max_year = 2022) %>%
  select(1:4)
head(US_frequencies)

## -----------------------------------------------------------------------------
## Add a vaccination column
US_frequencies <- US_frequencies %>%
  mutate(vaccination = c(rep(0, 77), seq(.5, .75, length = 26), .75, .75), # Add a vaccination column
         `A/H1N1` = `A/H1N1`*(1-vaccination), # Assume only non-vaccinated children have primary
         `A/H2N2` = `A/H2N2`*(1-vaccination), # infections; multiply the subtype-specific circulation
         `A/H3N2` = `A/H3N2`*(1-vaccination)) # fractions by one minus the year's vaccination probability.
tail(US_frequencies, n = 30)

## -----------------------------------------------------------------------------
Germany_frequencies <- get_country_cocirculation_data(country = 'Germany',
                                                      max_year = 2022) %>%
  select(1:4) %>%
  mutate(vaccination = c(rep(0, 87), seq(.05, .75, length = 16), .75, .75),
         `A/H1N1` = `A/H1N1`*(1-vaccination), # Assume only non-vaccinated children have primary
         `A/H2N2` = `A/H2N2`*(1-vaccination), # infections; multiply the subtype-specific circulation
         `A/H3N2` = `A/H3N2`*(1-vaccination))
tail(Germany_frequencies, 20)

## -----------------------------------------------------------------------------
## Check that all frequencies sum to 1
rowSums(US_frequencies[,2:5])
rowSums(Germany_frequencies[,2:5])

## -----------------------------------------------------------------------------
# Wrap the country-specific frequencies into a named list
input_list = list("United States" = US_frequencies,
                  "Germany" = Germany_frequencies)

## Calculate probabilities
get_imprinting_probabilities(observation_years = 2022, 
                             countries = c("United States", "Germany"), 
                             annual_frequencies = input_list, 
                             df_format = "wide")

