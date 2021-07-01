## ------------- ##
## Table of Data ##
## --------------##

# read data and merge
moddat <- merge(fies_subnat,
                covars,
                all.x=T, all.y=F) %>%
  na.omit %>%
  data.frame


## summary statistics  ---------------------------------------------------------
sumstat <- moddat %>% 
  dplyr::select(-c(YEAR, GDLCODE, ISO3, region, population, rural_perc)) %>% 
  dplyr::rename("Moderate Food Insecurity" = fies.mod,
                "Severe Food Insecurity" = fies.sev,
                "Stunting" = stunting,
                "Wasting" = wasting,
                "Mean Years of Schooling" = school_mean,
                "Topographic Ruggedness" = ruggedness,
                "GDP Per Capita" = gdp_percap,
                "Gini Coefficient" = gini,
                "Malaria Mortality Rate" = mal_falciparum,
                "Mean Annual Precipitation" = precip,
                "Poverty Headcount Index" = hci,
                "Mean Temperature" = tave,
                "Urban Percentage" = urban_perc,
                "Water Scarcity" = ws_share)

stargazer(sumstat, digits = 2, digits.extra = 0, summary.stat = c("n", "mean", "median", "sd", "min", "max"))


## share of observations by income  --------------------------------------------
cty <- ne_countries(returnclass='sf') %>%
  st_drop_geometry() %>% 
  filter(region_wb != 'Antarctica') %>% 
  select(ISO3 = iso_a3, region_wb = income_grp) %>%
  mutate(region_wb = case_when(region_wb == "1. High income: OECD" ~ "High Income",
                               region_wb == "2. High income: nonOECD" ~ "High Income",
                               region_wb == "3. Upper middle income" ~ "Upper-Middle Income",
                               region_wb == "4. Lower middle income" ~ "Lower-Middle Income",
                               region_wb == "5. Low income" ~ "Low Income"))


sharebyincome <- left_join(moddat, cty, by = "ISO3") %>% 
  group_by(region_wb) %>% 
  summarise(n = n(),
            share = n()/nrow(moddat)*100)


## list of all countries and years  --------------------------------------------
cy <- moddat %>% 
  dplyr::select(ISO3, YEAR) %>%
  dplyr::distinct() %>% 
  mutate(i = "X") %>% 
  pivot_wider(names_from = YEAR, values_from = i) %>% 
  replace(is.na(.), "") %>% 
  mutate(Country = countrycode(.$ISO3, origin = "iso3c", destination = "country.name")) %>% 
  dplyr::select(Country, ISO3, "2014", "2015", "2016", "2017", "2018")


stargazer(cy, summary = F, rownames = F)