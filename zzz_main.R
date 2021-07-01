
## set working directory -------------------------------------------------------

wd <- "C:/Users/bmuel/Desktop/GitHub/wu-fies_model"
setwd(wd)

path_cache <- file.path(wd, "cache")
path_data <- file.path(wd, "data")
path_lib <- file.path(wd, "lib")
path_code <- file.path(wd, "code")
path_viz <- file.path(wd, "viz")


## load libraries --------------------------------------------------------------

library(pacman)
p_load(tidyverse, countrycode, rnaturalearth, sf, zoo, broom,
       reshape2, lubridate, glmnet, data.table, BMS, mnormt,
       randomForest, randomForestSRC, imputeTS, telegram,
       cowplot, gridGraphics, rgeos, viridis, stargazer)


## rji telegram bot ------------------------------------------------------------
rji <- TGBot$new(token = bot_token('rji'))
rji$set_default_chat_id(user_id('me'))


## load data -------------------------------------------------------------------

load(file.path(path_cache, "fies_subnat.RData"))
load(file.path(path_cache, "covars.RData"))
load(file.path(path_cache, "gdl.RData"))
load(file.path(path_cache, "u5.population.RData"))


## source files ----------------------------------------------------------------

sapply(list.files(path_lib, full.names = TRUE), source)
