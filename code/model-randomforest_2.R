## ------------------ ##
## Random Forest (RF) ##
## -------------------##

## setup -----------------------------------------------------------------------
source("C:/Users/bmuel/Desktop/GitHub/wu-fies_model/zzz_main.R")

print("RF: Start with model-randomforest_2.R!")
rji$sendMessage("RF: Start with model-randomforest_2.R!")

# read data and merge
moddat <- merge(fies_subnat,
                covars,
                all.x=T, all.y=F) %>%
  na.omit %>%
  data.frame

# preddat <- covars %>%
#   data.frame

# set up model
vars <- names(moddat)[!names(moddat) %in% c("ISO3", "GDLCODE", "fies.mod.rur",
                                            "fies.sev.rur", "fies.mod.urb", "fies.sev.urb",
                                            "urban", "rural", "fies.sev", "fies.mod",
                                            "population", "YEAR", "rural_perc", "region")]

# apply logit transformation
moddat <- moddat %>%
  mutate(fies.mod.logit = logit(fies.mod),
         fies.sev.logit = logit(fies.sev))


## re-run LOOCV with optimal parameters ----------------------------------------
# from model_randomforest_1_wigeo.R
prm <- read.csv(file.path(path_data, "randomforest/prm-randomforest.csv"))

# what Matt got:
# prm <- data.frame(node = c(16, 14),
#                   mtry = 2, 
#                   depth = -1, 
#                   ix = c(1,2), 
#                   mod.rsq = c(1,0), 
#                   sev.rsq = c(0,1))
# mod: mtry = 2, nodesize = 16, depth = -1
# sev: mtry = 2, nodesize = 14, depth = -1

mod.prm <- prm[which.max(prm$mod.r2), ] 
sev.prm <- prm[which.max(prm$sev.r2), ] 

rji$sendMessage("RF: Start with re-run CV for best hyperparameters!")

moddat <- moddat %>% mutate(fies.mod.pred.cv = NA, fies.sev.pred.cv = NA)

i <- mod.prm$ix
for (iso3 in unique(moddat$ISO3)){
  print(paste0("RF mod: ", round(sum(!is.na(moddat$fies.mod.pred.cv))/nrow(moddat)*100, 2), "%"))
  rji$sendMessage(paste0("RF mod: ", round(sum(!is.na(moddat$fies.mod.pred.cv))/nrow(moddat)*100, 2), "%"))
  
  ix <- moddat$ISO3 != iso3
    
  mtry <- prm$mtry[i]
  node <- prm$node[i]
  if(prm$depth[i] < 0){
    depth <- prm$depth[i]
  } else{
    depth <- NULL
  }
  
  # moderate model
  rf.mod <- rfsrc(formula = as.formula(paste("fies.mod.logit", 
                                             paste(vars, collapse = "+"), sep= "~")),
                  data = moddat[ix, ],
                  ntree = 5000, 
                  mtry = mtry,
                  nodesize = node,
                  depth = depth)
  
  # project and apply inverse logit transformation
  moddat$fies.mod.pred.cv[!ix] <- inv.logit(predict(rf.mod, moddat[!ix,])$predicted)
}

i <- sev.prm$ix
for (iso3 in unique(moddat$ISO3)){
  print(paste0("RF sev: ", round(sum(!is.na(moddat$fies.sev.pred.cv))/nrow(moddat)*100, 2), "%"))
  rji$sendMessage(paste0("RF sev: ", round(sum(!is.na(moddat$fies.sev.pred.cv))/nrow(moddat)*100, 2), "%"))
  
  ix <- moddat$ISO3 != iso3
  
  mtry <- prm$mtry[i]
  node <- prm$node[i]
  if(prm$depth[i] < 0){
    depth <- prm$depth[i]
  } else{
    depth <- NULL
  }
  
  # severe model
  rf.sev <- rfsrc(formula = as.formula(paste("fies.sev.logit", 
                                             paste(vars, collapse = "+"), sep= "~")),
                  data = moddat[ix, ],
                  ntree = 5000, 
                  mtry = mtry,
                  nodesize = node,
                  depth = depth)
  
  # project and apply inverse logit transformation
  moddat$fies.sev.pred.cv[!ix] <- inv.logit(predict(rf.sev, moddat[!ix,])$predicted)
}

print("RF: Done!")
rji$sendMessage("RF: Done!")


## run full models -------------------------------------------------------------

# moderate model
rf.mod <- rfsrc(formula = as.formula(paste("fies.mod.logit", 
                                           paste(vars, collapse = "+"), sep= "~")),
                data = moddat,
                ntree = 10000, 
                mtry = mod.prm$mtry,
                nodesize = mod.prm$node,
                depth = mod.prm$depth,
                block.size = 1)

# severe model
rf.sev <- rfsrc(formula = as.formula(paste("fies.sev.logit", 
                                           paste(vars, collapse = "+"), sep= "~")),
                data = moddat,
                ntree = 10000, 
                mtry = sev.prm$mtry,
                nodesize = sev.prm$node,
                depth = sev.prm$depth,
                block.size = 1)

moddat$fies.mod.logit <- NULL
moddat$fies.sev.logit <- NULL

# project and apply inverse logit transformation
moddat$fies.sev.pred <- inv.logit(as.numeric(predict(rf.sev, moddat)$predicted))
# preddat$fies.sev.pred <- inv.logit(as.numeric(predict(rf.sev, preddat)$predicted))
moddat$fies.mod.pred <- inv.logit(as.numeric(predict(rf.mod, moddat)$predicted))
# preddat$fies.mod.pred <- inv.logit(as.numeric(predict(rf.mod, preddat)$predicted))

print("RF: Done!")
rji$sendMessage("RF: Done!")

# write.csv(preddat, file.path(path_data, "randomforest/preddat-randomforest.csv"), row.names = F)
# write.csv(moddat, file.path(path_data, "randomforest/moddat-randomforest.csv"), row.names = F)
save(list=c("rf.mod", "rf.sev"), file = file.path(path_data, "randomforest/models-randomforest.RData"))

print("RF: Done with model-randomforest_2.R!")
rji$sendMessage("RF: Done with model-randomforest_2.R!")
