## ------------------ ##
## Random Forest (RF) ##
## -------------------##


# on wigeo server

## setup -----------------------------------------------------------------------
# load packages
library(pacman)
p_load(tidyverse, countrycode, rnaturalearth, sf, zoo, broom,
       reshape2, lubridate, glmnet, data.table, BMS, mnormt,
       randomForest, randomForestSRC, imputeTS, telegram,
       cowplot, gridGraphics, rgeos, viridis, stargazer)

# functions
logit <- function(x){log(x/(1-x))}
inv.logit <- function(x){1/(1 + exp(-x))}
fr2 <- function(y, yhat) {(1-(sum((y-yhat)^2)/sum((y-mean(yhat))^2)))}

# telegram bot
rji <- TGBot$new(token = bot_token('rji'))
rji$set_default_chat_id(user_id('me'))

try(rji$sendMessage("RF: Start with model-randomforest_1_wigeo.R!"))

# read data and merge
load("./input/fies_subnat.RData")
load("./input/covars.RData")

moddat <- merge(fies_subnat,
                covars,
                all.x=T, all.y=F) %>%
  na.omit %>%
  data.frame

# set up model
vars <- names(moddat)[!names(moddat) %in% c('ISO3', 'GDLCODE', 'fies.mod.rur',
                                            'fies.sev.rur', 'fies.mod.urb', 'fies.sev.urb',
                                            'urban', 'rural', 'fies.sev', 'fies.mod',
                                            'population', 'YEAR', 'rural_perc', 'region')]

# apply logit transformation
moddat <- moddat %>%
  mutate(fies.mod.logit = logit(fies.mod),
         fies.sev.logit = logit(fies.sev))


## find optimal tuning parameters ----------------------------------------------

# set up parameter combinations
prm <- expand.grid(list(node=c(1:10, 12, 14, 16, 18, 20),
                        mtry=c(12:1),
                        depth=c(1:12, -1)))
prm$ix <- 1:nrow(prm)

prm$mod.r2 <- NA
prm$sev.r2 <- NA

s <- 234*2
prml <- split(prm, rep(1:ceiling(nrow(prm)/s), each=s, length.out=nrow(prm)))


# try all combinations and calculate R2
for(p in 1:length(prml)) {
  
  prm <- prml[[p]]
  
  try(rji$sendMessage(paste0(p, "RF: Start with CV at the country level!")))
  
  # create a list (in a very ugly way)
  ls <- rep(list(moddat %>% dplyr::select(ISO3,
                                          fies.mod, fies.sev,
                                          fies.mod.logit, fies.sev.logit)), 
            length(unique(moddat$ISO3)))
  ls <- rep(ls, length(prm$ix))
  x <- 1
  for(i in prm$ix) {
    for(iso3 in unique(moddat$ISO3)) {
      ls[[x]]$i <- i
      ls[[x]]$ix <- ifelse(ls[[x]]$ISO3 == iso3, FALSE, TRUE)
      
      names(ls)[x] <- paste0(i, "_", iso3)
      x <- x + 1
    }
  }
  v <- moddat %>% dplyr::select(all_of(vars))
  
  # use multiple cores
  cl <- makeSOCKcluster(25)
  registerDoSNOW(cl)
  progress <- function(n) {
    print(paste0(p, "RF: ", round(n/length(ls)*100, 2), "%"))
    try(rji$sendMessage(paste0(p, "RF: ", round(n/length(ls)*100, 2), "%")))
  }
  opts <- list(progress=progress)
  
  dat <- NA
  dat <- foreach(l=1:length(ls), .options.snow=opts, .combine = "rbind", .packages="randomForestSRC") %dopar% {
    d <- ls[[l]]
    d <- cbind(d, v)
    
    mtry <- prm$mtry[prm$ix == d$i[1]]
    node <- prm$node[prm$ix == d$i[1]]
    if(prm$depth[prm$ix == d$i[1]] < 0){
      depth <- prm$depth[prm$ix == d$i[1]]
    } else{
      depth <- NULL
    }
    
    # moderate model
    rf.mod <- randomForestSRC::rfsrc(formula = as.formula(paste("fies.mod.logit", paste(vars, collapse = "+"), sep= "~")),
                                     # forest = TRUE
                                     data = d[d$ix, ],
                                     ntree = 1000,
                                     mtry = mtry,
                                     nodesize = node,
                                     depth = depth)
    
    # severe model
    rf.sev <- randomForestSRC::rfsrc(formula = as.formula(paste("fies.sev.logit", paste(vars, collapse = "+"), sep= "~")),
                                     # forest = TRUE
                                     data = d[d$ix, ],
                                     ntree = 1000,
                                     mtry = mtry,
                                     nodesize = node,
                                     depth = depth)
    
    # project and apply inverse logit transformation
    d$fies.mod.pred.cv[!d$ix] <- inv.logit(predict(rf.mod, d[!d$ix,])$predicted)
    d$fies.sev.pred.cv[!d$ix] <- inv.logit(predict(rf.sev, d[!d$ix,])$predicted)
    
    return(d[!d$ix, c("i", "fies.mod", "fies.mod.pred.cv", "fies.sev", "fies.sev.pred.cv")])
  }
  
  stopCluster(cl)
  saveRDS(dat, paste0("dat/dat", p ,".rds"))
  
  # calculate R2
  for(i in prm$ix) {
    prm$mod.r2[prm$ix == i] <- fr2(dat[dat$i == i,]$fies.mod, dat[dat$i == i,]$fies.mod.pred.cv)
    prm$sev.r2[prm$ix == i] <- fr2(dat[dat$i == i,]$fies.sev, dat[dat$i == i,]$fies.sev.pred.cv)
  }
  
  saveRDS(prm, paste0("prm/prm", p ,".rds"))
  
}


# combine all prm and write csv
prml <- lapply(list.files(path = "./prm", pattern="*.rds", full.names = T), readRDS)
prm <- do.call(rbind, prml)

# write.csv(prm, "prm-randomforest.csv", row.names = F)

try(rji$sendMessage("RF: Done!"))
try(rji$sendMessage("RF: Done with model-randomforest_1_wigeo.R!"))