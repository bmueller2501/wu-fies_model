## ---------------------------- ##
## Table with Accuracy Measures ##
## -----------------------------## 

approach <- c("lasso", "bayes_ssvs", "randomforest")

# set up accuracy data frame
acc <- data.frame(expand.grid(list(model = c("moderat", "severe"), approach = approach, 
                                   mae = NA, rmse = NA, r2 = NA,
                                   mae_lmi = NA, rmse_lmi = NA, r2_lmi = NA,
                                   mae_li = NA, rmse_li = NA, r2_li = NA
                                   ))) %>% 
  mutate(approach = as.character(approach),
         model = as.character(model))

## merge income groups ---------------------------------------------------------
a <- "randomforest"
moddat <- read.csv(file.path(path_data, a, paste0("moddat-", a,".csv")))

cty <- ne_countries(returnclass='sf') %>%
  filter(region_wb != 'Antarctica')

region_cols <- c("High Income"=magma(9)[2],
                 "Upper-Middle Income"=magma(9)[4],
                 "Lower-Middle Income"=magma(9)[6],
                 "Low Income"=magma(9)[8])

regdat <- cty %>% 
  select(ISO3 = iso_a3, region_wb = income_grp) %>% #use world bank regions by income
  filter(ISO3 %in% unique(moddat$ISO3)) %>% 
  # rename and change order from high income to low income
  mutate(region_wb = case_when(region_wb == "1. High income: OECD" ~ "High Income",
                               region_wb == "2. High income: nonOECD" ~ "High Income",
                               region_wb == "3. Upper middle income" ~ "Upper-Middle Income",
                               region_wb == "4. Lower middle income" ~ "Lower-Middle Income",
                               region_wb == "5. Low income" ~ "Low Income")) %>% 
  mutate(region_wb = factor(region_wb),
         region_wb = ordered(region_wb, levels = c("High Income", "Upper-Middle Income",
                                                   "Lower-Middle Income", "Low Income")))


## calculate MAE, RMSE and R2 --------------------------------------------------
for(a in approach) {
  
  moddat <- read.csv(file.path(path_data, a, paste0("moddat-", a,".csv")))
  
  moddat <- merge(moddat, st_drop_geometry(regdat)) %>%
    select(-region) %>%
    rename(region = "region_wb") %>%
    mutate(region = factor(region, levels=rev(c("High Income", "Upper-Middle Income",
                                                "Lower-Middle Income", "Low Income"))))
  
  
  # prediction errors
  moddat$res.mod <- moddat$fies.mod - moddat$fies.mod.pred.cv
  moddat$res.sev <- moddat$fies.sev - moddat$fies.sev.pred.cv
  
  # mean absolute error
  acc$mae[acc$approach==a&acc$model=="moderat"] <- fmae(moddat$fies.mod, moddat$fies.mod.pred.cv)
  acc$mae[acc$approach==a&acc$model=="severe"] <- fmae(moddat$fies.sev, moddat$fies.sev.pred.cv)
  
  # mean absolute error
  acc$rmse[acc$approach==a&acc$model=="moderat"] <- frmse(moddat$fies.mod, moddat$fies.mod.pred.cv)
  acc$rmse[acc$approach==a&acc$model=="severe"] <- frmse(moddat$fies.sev, moddat$fies.sev.pred.cv)
  
  # R^2
  acc$r2[acc$approach==a&acc$model=="moderat"] <- fr2(moddat$fies.mod, moddat$fies.mod.pred.cv)
  acc$r2[acc$approach==a&acc$model=="severe"] <- fr2(moddat$fies.sev, moddat$fies.sev.pred.cv)
  
  # MAE, RMSE and R^2 for lower-middle income and lower income
  m <- moddat %>% filter(region == "Lower-Middle Income")
  acc$mae_lmi[acc$approach==a&acc$model=="moderat"] <- fmae(m$fies.mod, m$fies.mod.pred.cv)
  acc$mae_lmi[acc$approach==a&acc$model=="severe"] <- fmae(m$fies.sev, m$fies.sev.pred.cv)
  acc$rmse_lmi[acc$approach==a&acc$model=="moderat"] <- frmse(m$fies.mod, m$fies.mod.pred.cv)
  acc$rmse_lmi[acc$approach==a&acc$model=="severe"] <- frmse(m$fies.sev, m$fies.sev.pred.cv)
  acc$r2_lmi[acc$approach==a&acc$model=="moderat"] <- fr2(m$fies.mod, m$fies.mod.pred.cv)
  acc$r2_lmi[acc$approach==a&acc$model=="severe"] <- fr2(m$fies.sev, m$fies.sev.pred.cv)

  
  m <- moddat %>% filter(region == "Low Income")
  acc$mae_li[acc$approach==a&acc$model=="moderat"] <- fmae(m$fies.mod, m$fies.mod.pred.cv)
  acc$mae_li[acc$approach==a&acc$model=="severe"] <- fmae(m$fies.sev, m$fies.sev.pred.cv)
  acc$rmse_li[acc$approach==a&acc$model=="moderat"] <- frmse(m$fies.mod, m$fies.mod.pred.cv)
  acc$rmse_li[acc$approach==a&acc$model=="severe"] <- frmse(m$fies.sev, m$fies.sev.pred.cv)
  acc$r2_li[acc$approach==a&acc$model=="moderat"] <- fr2(m$fies.mod, m$fies.mod.pred.cv)
  acc$r2_li[acc$approach==a&acc$model=="severe"] <- fr2(m$fies.sev, m$fies.sev.pred.cv)
  
}

acc <- acc %>% dplyr::arrange(model)


## create table ----------------------------------------------------------------
t <- acc %>% 
  dplyr::select(Model = model, Approach = approach, 
                MAE = mae, RMSE = rmse, R2 = r2,
                # "MAE - Lower-Middle Income" = mae_lmi,
                # "RMSE - Lower-Middle Income" = rmse_lmi,
                # "R2 - Lower-Middle Income" = r2_lmi,
                "MAE - Low Income" = mae_li,
                "RMSE - Low Income" = rmse_li,
                "R2 - Low Income" = r2_li)
t$Model[t$Model == "moderat"] <- "Moderat"
t$Model[t$Model == "severe"] <- "Severe"
t$Approach[t$Approach == "lasso"] <- "LASSO"
t$Approach[t$Approach == "bayes_ssvs"] <- "BSSVS"
t$Approach[t$Approach == "randomforest"] <- "RF"


stargazer(t, summary = F, digits = 4, rownames = F)