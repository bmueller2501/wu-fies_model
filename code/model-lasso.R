## ------------------------------------------------------- ##
## Least Absolute Shrinkage and Selection Operator (LASSO) ##
## --------------------------------------------------------## 

## setup -----------------------------------------------------------------------
source("C:/Users/bmuel/Desktop/GitHub/wu-fies_model/zzz_main.R")

print("LAS: Start with model-lasso.R!")
try(rji$sendMessage("LAS: Start with model-lasso.R!"))

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


## find optimal tuning parameters ----------------------------------------------
try(rji$sendMessage("LAS: Start with CV at the country level!"))


x <- model.matrix(as.formula(paste0('fies.mod.logit ~ ', paste0(vars, collapse=' + '))), 
                  data=moddat)

# check first at which lambda value we have more than intercept
lam <- seq(0, 1, by = 0.001)
lam <- rev(lam)

n <- 2
for(l in lam) {
  mod <- glmnet(x, moddat$fies.mod.logit, alpha=1, lambda = l)
  
  tmp_coeffs <- coef(mod)
  mdf <- data.frame(term = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], estimate =tmp_coeffs@x,
                    stringsAsFactors=F)
  
  if(nrow(mdf) > n) {tlam <- l; n <- Inf}
}

tlam <- floor(tlam*100)/100 #round to lower full value and hard code for mod and sev

# set up parameter combinations
prm <- expand.grid(list(lam=seq(0, tlam, by = tlam/10000)))
prm$ix <- 1:nrow(prm)

prm$mod.r2 <- NA
prm$sev.r2 <- NA

# try all combinations and calculate R2
for(i in sample(prm$ix[is.na(prm$sev.r2)])) {
  print(paste0("LAS: ", round(sum(!is.na(prm$sev.r2))/nrow(prm)*100, 2), "%"))
  rji$sendMessage(paste0("LAS: ", round(sum(!is.na(prm$sev.r2))/nrow(prm)*100, 2), "%"))
  for (iso3 in unique(moddat$ISO3)){

    
    ix <- moddat$ISO3 != iso3
    dat <- moddat[ix, ]
    
    l <- prm$lam[i]
    
    # moderate model
    x <- model.matrix(as.formula(paste0('fies.mod.logit ~ ', paste0(vars, collapse=' + '))), 
                      data=dat)
    mod <- glmnet(x, dat$fies.mod.logit, alpha=1, lambda = l)
    
    tmp_coeffs <- coef(mod)
    mdf <- data.frame(term = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], estimate =tmp_coeffs@x,
                      stringsAsFactors=F)
    
    
    # severe model
    x <- model.matrix(as.formula(paste0('fies.sev.logit ~ ', paste0(vars, collapse=' + '))), 
                      data=dat)
    mod <- glmnet(x, dat$fies.sev.logit, alpha=1, lambda = l)
    
    tmp_coeffs <- coef(mod)
    sdf <- data.frame(term = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], estimate = tmp_coeffs@x,
                      stringsAsFactors=F)
    
    
    # project and apply inverse logit transformation
    moddat$fies.mod.pred.cv[!ix] <- mdf$estimate[mdf$term == '(Intercept)']
    moddat$fies.sev.pred.cv[!ix] <- sdf$estimate[sdf$term == '(Intercept)']
    for (j in 2:nrow(mdf)) moddat$fies.mod.pred.cv[!ix] <- 
      moddat$fies.mod.pred.cv[!ix] + moddat[!ix, mdf$term[j]]*mdf$estimate[j]
    for (j in 2:nrow(sdf)) moddat$fies.sev.pred.cv[!ix] <- 
      moddat$fies.sev.pred.cv[!ix] + moddat[!ix, sdf$term[j]]*sdf$estimate[j]
    moddat$fies.mod.pred.cv[!ix] <- inv.logit(moddat$fies.mod.pred.cv[!ix])
    moddat$fies.sev.pred.cv[!ix] <- inv.logit(moddat$fies.sev.pred.cv[!ix])
    
  }

  # calculate R2
  prm$mod.r2[i] <- fr2(moddat$fies.mod, moddat$fies.mod.pred.cv)
  prm$sev.r2[i] <- fr2(moddat$fies.sev, moddat$fies.sev.pred.cv)
  
}


# write.csv(prm, file.path(path_data, "lasso/prm-lasso.csv"), row.names = F)

mod.prm <- prm[which.max(prm$mod.r2), ]
sev.prm <- prm[which.max(prm$sev.r2), ]

try(rji$sendMessage("LAS: Done!"))


## re-run LOOCV with optimal parameters ----------------------------------------
try(rji$sendMessage("LAS: Start with re-run CV for best hyperparameters!"))

moddat <- moddat %>% mutate(fies.mod.pred.cv = NA, fies.sev.pred.cv = NA)

i <- mod.prm$ix
for (iso3 in unique(moddat$ISO3)){
  
  print(paste0("LAS mod: ", round(sum(!is.na(moddat$fies.mod.pred.cv))/nrow(moddat)*100, 2), "%"))
  try(rji$sendMessage(paste0("LAS mod: ", round(sum(!is.na(moddat$fies.mod.pred.cv))/nrow(moddat)*100, 2), "%")))
  
  ix <- moddat$ISO3 != iso3
  dat <- moddat[ix, ]
  
  l <- prm$lam[i]

  # moderate model
  x <- model.matrix(as.formula(paste0('fies.mod.logit ~ ', paste0(vars, collapse=' + '))), 
                    data=dat)
  mod <- glmnet(x, dat$fies.mod.logit, alpha=1, lambda = l)
  
  tmp_coeffs <- coef(mod)
  mdf <- data.frame(term = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], estimate =tmp_coeffs@x,
                    stringsAsFactors=F)
  
  # project and apply inverse logit transformation
  moddat$fies.mod.pred.cv[!ix] <- mdf$estimate[mdf$term == '(Intercept)']
  for (j in 2:nrow(mdf)) moddat$fies.mod.pred.cv[!ix] <- moddat$fies.mod.pred.cv[!ix] + moddat[!ix, mdf$term[j]]*mdf$estimate[j]
  moddat$fies.mod.pred.cv[!ix] <- inv.logit(moddat$fies.mod.pred.cv[!ix])
}

i <- sev.prm$ix
for (iso3 in unique(moddat$ISO3)){
  
  print(paste0("BAY sev: ", round(sum(!is.na(moddat$fies.sev.pred.cv))/nrow(moddat)*100, 2), "%"))
  try(rji$sendMessage(paste0("BAY sev: ", round(sum(!is.na(moddat$fies.sev.pred.cv))/nrow(moddat)*100, 2), "%")))
  
  ix <- moddat$ISO3 != iso3
  dat <- moddat[ix, ]
  
  l <- prm$lam[i]
  
  # severe model
  x <- model.matrix(as.formula(paste0('fies.sev.logit ~ ', paste0(vars, collapse=' + '))), 
                    data=dat)
  mod <- glmnet(x, dat$fies.sev.logit, alpha=1, lambda = l)
  
  tmp_coeffs <- coef(mod)
  sdf <- data.frame(term = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], estimate = tmp_coeffs@x,
                    stringsAsFactors=F)
  
  # project and apply inverse logit transformation
  moddat$fies.sev.pred.cv[!ix] <- sdf$estimate[sdf$term == '(Intercept)']
  for (j in 2:nrow(sdf)) moddat$fies.sev.pred.cv[!ix] <- moddat$fies.sev.pred.cv[!ix] + moddat[!ix, sdf$term[j]]*sdf$estimate[j]
  moddat$fies.sev.pred.cv[!ix] <- inv.logit(moddat$fies.sev.pred.cv[!ix])
}

print("LAS: Done!")
try(rji$sendMessage("LAS: Done!"))


## run full models -------------------------------------------------------------
print("LAS: Start with run full models with best parameters!")
try(rji$sendMessage("LAS: Start with run full models with best parameters!"))

# moderate model
l <- mod.prm$l

x <- model.matrix(as.formula(paste0('fies.mod.logit ~ ', paste0(vars, collapse=' + '))), 
                  data=moddat)
mod <- glmnet(x, moddat$fies.mod.logit, alpha=1, lambda = l)

tmp_coeffs <- coef(mod); tmp_coeffs
mdf <- data.frame(term = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], estimate =tmp_coeffs@x,
                  stringsAsFactors=F)

# severe model
l <- sev.prm$l

x <- model.matrix(as.formula(paste0('fies.sev.logit ~ ', paste0(vars, collapse=' + '))), 
                  data=moddat)
mod <- glmnet(x, moddat$fies.sev.logit, alpha=1, lambda = l)

tmp_coeffs <- coef(mod); tmp_coeffs
sdf <- data.frame(term = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], estimate = tmp_coeffs@x,
                  stringsAsFactors=F)

moddat$fies.mod.logit <- NULL
moddat$fies.sev.logit <- NULL

# project and apply inverse logit transformation
moddat$fies.mod.pred <- mdf$estimate[mdf$term == '(Intercept)']
moddat$fies.sev.pred <- sdf$estimate[sdf$term == '(Intercept)']
# preddat$fies.mod.pred <- mdf$estimate[mdf$term == '(Intercept)']
# preddat$fies.sev.pred <- sdf$estimate[sdf$term == '(Intercept)']
for (i in 2:nrow(mdf)) {
  moddat$fies.mod.pred <- moddat$fies.mod.pred + moddat[ , mdf$term[i]]*mdf$estimate[i]
  # project and apply inverse logit transformation
}
for (i in 2:nrow(sdf)) {
  moddat$fies.sev.pred <- moddat$fies.sev.pred + moddat[ , mdf$term[i]]*mdf$estimate[i]
  # project and apply inverse logit transformation
}
moddat$fies.mod.pred <- inv.logit(moddat$fies.mod.pred)
moddat$fies.sev.pred <- inv.logit(moddat$fies.sev.pred)
# preddat$fies.mod.pred <- inv.logit(preddat$fies.mod.pred)
# preddat$fies.sev.pred <- inv.logit(preddat$fies.sev.pred)


print("LAS: Done!")
try(rji$sendMessage("LAS: Done!"))

# write.csv(moddat, file.path(path_data, "lasso/moddat-lasso.csv"), row.names = F)
# write.csv(preddat, file.path(path_data, "lasso/preddat-lasso.csv"), row.names = F)

print("LAS: Done with model-lasso.R!")
try(rji$sendMessage("LAS: Done with model-lasso.R!"))
