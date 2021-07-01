## --------------------------- ##
## Visualize Model Predictions ##
## ----------------------------##

source("C:/Users/bmuel/Desktop/GitHub/wu-fies_model/zzz_main.R")
options(scipen=100)
s <- TRUE #save plots?


## world bank income map  ------------------------------------------------------
a <- "randomforest"
moddat <- read.csv(file.path(path_data, a, paste0("moddat-", a,".csv")))

cty <- ne_countries(returnclass='sf') %>%
  filter(region_wb != 'Antarctica')

region_cols <- c("High Income"=magma(9)[8],
                 "Upper-Middle Income"=magma(9)[6],
                 "Lower-Middle Income"=magma(9)[4],
                 "Low Income"=magma(9)[2])

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
p <- ggplot(regdat) + 
  geom_sf(aes(fill=as.factor(region_wb)), color=NA) + 
  coord_sf(crs='+proj=robin') + 
  scale_fill_manual(values=region_cols) + 
  guides(fill=guide_legend(ncol=2)) + 
  geom_sf(data=cty, color='#000000', fill=NA, size=0.15) +
  labs(fill='') +
  theme_void() +
  theme(legend.position = "bottom"); p
leg <- get_legend(p)
if(s) ggsave(p, filename = file.path(path_viz, paste0("map-worldbank_income.pdf")), width = 8, height = 4)


## out-of-sample all modeling approaches----------------------------------------
approach <- c("lasso", "bayes_ssvs", "randomforest")

modvsobs <- list()
resvsmod <- list()


for(a in approach) {
  
  moddat <- read.csv(file.path(path_data, a, paste0("moddat-", a,".csv")))
  preddat <- read.csv(file.path(path_data, a, paste0("preddat-", a,".csv")))
  
  moddat <- merge(moddat, st_drop_geometry(regdat)) %>%
    select(-region) %>%
    rename(region = "region_wb") %>%
    mutate(region = factor(region, levels=rev(c("High Income", "Upper-Middle Income",
                                                "Lower-Middle Income", "Low Income"))))
  
  # modeled vs observed rates
  modvsobs[[paste0(a, ".mod")]] <- ggplot(moddat) + 
    geom_point(aes(x=fies.mod, y=fies.mod.pred.cv, color = region), alpha = 0.6) + 
    scale_color_manual(values=region_cols) + 
    geom_abline(intercept = 0, slope = 1, color = "black", size = 0.5, linetype = "dashed") +
    xlim(0, 1) + ylim(0, 1) +
    theme_bw() +
    theme(legend.position = "none") +
    # labs(caption=paste0("Mean Average Error: ",  round(mae.mod, 4),
    #                     "\nR-Squared: ", round(r2.mod, 4)),
    #      x="Observed Moderate Food Insecurity",
    #      y="Modeled Moderate Food Insecurity")
    labs(x="Observed Moderate Food Insecurity",
         y="Modeled Moderate Food Insecurity")
  
  modvsobs[[paste0(a, ".sev")]] <- ggplot(moddat) + 
    geom_point(aes(x=fies.sev, y=fies.sev.pred.cv, color = region), alpha = 0.6) + 
    scale_color_manual(values=region_cols) + 
    geom_abline(intercept = 0, slope = 1, color = "black", size = 0.5, linetype = "dashed") +
    xlim(0, 1) + ylim(0, 1) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x="Observed Severe Food Insecurity",
         y="Modeled Severe Food Insecurity")
  
  
  # error vs modeled rates
  moddat <- moddat %>% 
    mutate(rescv.mod = fies.mod - fies.mod.pred.cv,
           rescv.sev = fies.sev - fies.sev.pred.cv)
  
  resvsmod[[paste0(a, ".mod")]] <- ggplot(moddat) + 
    geom_point(aes(x=fies.mod.pred.cv, y=rescv.mod, color = region), alpha = 0.6) + 
    scale_color_manual(values=region_cols) +
    geom_abline(intercept = 0, slope = 0, color = "black", size = 0.5, linetype = "dashed") +
    xlim(0, 1) + 
    ylim(-0.4, 0.4) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x="Modeled Moderate Food Insecurity",
         y="Prediction Error")
  
  resvsmod[[paste0(a, ".sev")]] <- ggplot(moddat) + 
    geom_point(aes(x=fies.sev.pred.cv, y=rescv.sev, color = region), alpha = 0.6) + 
    scale_color_manual(values=region_cols) +
    geom_abline(intercept = 0, slope = 0, color = "black", size = 0.5, linetype = "dashed") +
    xlim(0, 1) + 
    ylim(-0.4, 0.4) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x="Modeled Severe Food Insecurity",
         y="Prediction Error")
  
}


plot_grid(plot_grid(ggdraw()+draw_label("LASSO", fontface = "bold"),
          plot_grid(modvsobs[["lasso.mod"]], modvsobs[["lasso.sev"]], nrow = 1, labels = "AUTO"),
          ncol = 1, rel_heights=c(0.1, 0.9)),
          
          plot_grid(ggdraw()+draw_label("BSSVS", fontface = "bold"),
          plot_grid(modvsobs[["bayes_ssvs.mod"]], modvsobs[["bayes_ssvs.sev"]], nrow = 1, labels = "AUTO"),
          ncol = 1, rel_heights=c(0.1, 0.9)),    
          
          plot_grid(ggdraw()+draw_label("RF", fontface = "bold"),
          plot_grid(modvsobs[["randomforest.mod"]], modvsobs[["randomforest.sev"]], nrow = 1, labels = "AUTO"),
          ncol = 1, rel_heights=c(0.1, 0.9)),
          
          leg,
          
          rel_heights=c(0.3, 0.3, 0.3, 0.1),
          ncol = 1)
if(s) ggsave(filename = file.path(path_viz, paste0("out_sample-modvsobs.pdf")), width = 10, height = 12)


plot_grid(plot_grid(ggdraw()+draw_label("LASSO", fontface = "bold"),
                    plot_grid(resvsmod[["lasso.mod"]], resvsmod[["lasso.sev"]], nrow = 1, labels = "AUTO"),
                    ncol = 1, rel_heights=c(0.1, 0.9)),
          
          plot_grid(ggdraw()+draw_label("BSSVS", fontface = "bold"),
                    plot_grid(resvsmod[["bayes_ssvs.mod"]], resvsmod[["bayes_ssvs.sev"]], nrow = 1, labels = "AUTO"),
                    ncol = 1, rel_heights=c(0.1, 0.9)),    
          
          plot_grid(ggdraw()+draw_label("RF", fontface = "bold"),
                    plot_grid(resvsmod[["randomforest.mod"]], resvsmod[["randomforest.sev"]], nrow = 1, labels = "AUTO"),
                    ncol = 1, rel_heights=c(0.1, 0.9)),
          
          leg,
          
          rel_heights=c(0.3, 0.3, 0.3, 0.1),
          ncol = 1)
if(s) ggsave(filename = file.path(path_viz, paste0("out_sample-errorvsmod.pdf")), width = 10, height = 12)