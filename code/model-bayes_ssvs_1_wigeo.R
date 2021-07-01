## --------------------------------------------------------- ##
## Bayesian Stochastic Search and Variable Selection (BSSVS) ##
## ----------------------------------------------------------## 


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

try(rji$sendMessage("BAY: Start with model-bayes_ssvs_1_wigeo.R!"))

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


## BSSVS function --------------------------------------------------------------
# combination of "BaysianOLS.R" and "SSVS.R"

bayes.ssvs <- function(y, bigX, nburn, nsave, spike, slab, fig) {
  
  # ssvs <- 1 #scale the distributions of the spike and slab prior, something between 1-10 maybe, 1 is standard (from course)
  # 
  # # bring data into corrrect format
  # data <- moddat %>%
  #   select(fies.mod, fies.sev, vars)
  # summary(data)
  # 
  # # model moderate food insecurity
  # y <- as.matrix(data %>% select(fies.mod))
  # bigX <- cbind(matrix(1, nrow = nrow(data), ncol = 1), as.matrix(data %>% select(vars)))
  # colnames(bigX) <- c("(Intercept)", vars)
  # nburn <- 10000         #number of burn ins
  # nsave <- 10000         #number of saved draws
  # 
  # fig <- T
  
  y <- y
  X <- bigX
  N <- nrow(X)
  
  # ------------------ #
  # --- Check data --- #
  # ------------------ #
  
  # let"s look at OLS estimates
  # coefficient estimates
  # (X"X)^(-1) X"y
  beta.ols   <- solve(t(X) %*% X, tol = rcond(t(X) %*% X)*0.99999) %*% t(X) %*% y
  
  # variance estimate
  errors.ols <- y - X %*% beta.ols
  sig2.ols   <- t(errors.ols) %*% errors.ols / (N-ncol(X))
  
  # you can get the same via Rs linear regression function
  # note the -1 is necessary as we already included an intercept in X
  # reg <- lm(y~X-1)
  # summary(reg)
  
  
  # --------------------------- #
  # --- GIBBS SAMPLER SETUP --- #
  # --------------------------- #
  
  nsave <- nsave         #number of saved draws
  nburn <- nburn         #number of burn ins
  ntot  <- nsave + nburn #number of total iterations
  P     <- ncol(X)       #number of explanatories
  N     <- nrow(X)       #number of observations
  
  
  # ------------------- #
  # --- PRIOR SETUP --- #
  # ------------------- #
  
  # ssvs prior on coefficients
  sig2h <- slab        #ssvs slab  (high) variance
  sig2l <- spike       #ssvs spike (low)  variance
  
  b0    <- matrix(0, nrow=P)   #prior mean of beta
  B0    <- diag(sig2h,  P)     #prior var-cov of beta
  B0inv <- diag(1/sig2h, P)    #inverse of prior var-cov of beta
  
  # how does our normal - normal spike and slab prior look like?
  #par(mfrow=c(1,1),mar=c(1,1,1,1))
  if(fig == T) {
    # png("figures/bayes_ssvs/spike_slab_prior_bayes_ssvs.png")
    plot(density(rnorm(10000, 0, sqrt(sig2h))),ylim=c(0,2)) #slab
    lines(density(rnorm(10000, 0, sqrt(sig2l))), col="red")   #spike
    # dev.off()
  }
  
  # inverse gamma prior on the variance
  c0 <- 2 #prior shape of sigma2
  C0 <- 1 #prior rate  of sigma2
  
  # how does this prior distribution look like?
  #png("figures/bayes_ssvs/variance_prior_bayes_ssvs_wt.png")
  par(mfrow=c(1,1), mar=c(2,2,2,2))
  if(fig == T) plot(density(1/rgamma(20000, c0, C0)))
  #dev.off()
  
  
  # note that we implicitly set prior inclusion probability for all variables to 0.5
  # could also implement a version with other prior inclusion probabilities easily
  
  
  # ----------------------- #
  # --- STARTING VALUES --- #
  # ----------------------- #
  
  # set all inclusion indicators to 1
  delta.draw <- matrix(1, P, 1)
  
  # set variance to 1
  sig2.draw  <- 1
  
  # set coefficients to 0
  beta.draw  <- matrix(0, P, 1)
  
  # --------------- #
  # --- STORAGE --- #
  # --------------- #
  
  beta.store  <- matrix(NA, nsave, P)
  sig2.store  <- matrix(NA, nsave, 1)
  delta.store <- matrix(NA, nsave, P)
  
  # -------------------------------- #
  # --- GIBBS SAMPLING ALGORITHM --- #
  # -------------------------------- #
  
  # compute sufficient statistics
  XX <- t(X)%*%X
  XY <- t(X)%*%y
  
  # progess bar
  # pb <- txtProgressBar(min = 0, max = ntot, style = 3)
  
  
  for (irep in 1:ntot){ # MCMC LOOP START
    
    # STEP 1: SAMPLE BETA GIVEN SIGMA & DATA
    
    # compute posterior quantities
    Bn_inv <- B0inv + XX / sig2.draw
    Bn     <- solve(Bn_inv)
    bn     <- Bn %*% (B0inv %*% b0 + XY / sig2.draw)
    
    # draw from a multivariate normal distribution
    # (there are several ways to do that, here I use a package for simplicity)
    beta.draw <- mnormt::rmnorm(1, bn, Bn)
    
    
    # STEP 2: SAMPLE SIGMA2 GIVEN BETA & DATA
    
    # compute errors
    e  <- y - X %*% beta.draw
    # sum of squared residuals (equivalent to sum((e)^2))
    ee <- t(e) %*% e
    
    # compute posterior quantities
    ck <- c0 +  N/2
    Ck <- C0 + ee/2
    
    # sample sigma2 from inverse gamma posterior
    sig2.draw <- 1/rgamma(1, ck, Ck)
    
    
    
    # STEP 3: STOCHASTIC SEARCH VARIABLE SELECTION
    
    # this will be WAY more numerically stable if you do everything in log space
    # however, for didactic purposes i stay away from the logs
    # it"s easier to see what happens in that way
    
    # likelihood that coefficients come from spike component (=exclusion)
    ll.spike    <- dnorm(beta.draw, 0, sqrt(sig2l))
    
    # likelihood that coefficients come from slab component (=inclusion)
    ll.slab     <- dnorm(beta.draw, 0, sqrt(sig2h))
    
    # compute inclusion probability
    # to get probabilities, we have to normalize!
    # here, we could in principal add our prior and multiply it with the likelihood
    # implicitly, prior inclusion probability = 0.5 and it drops out of the equation!
    pip         <- ll.slab / (ll.spike + ll.slab)
    
    # now we sample our delta inclusion indicators by flipping a coin with
    # success probability being equal to the posterior inclusion probability
    delta.draw  <- ifelse(pip > runif(P), 1, 0)
    
    # finally, use the inclusion indicator delta to update your prior var-cov 
    # coefficients that are included get the large variance
    # coefficients that are excluded get the low variance
    scaling     <- delta.draw * sig2h + (1-delta.draw) * sig2l
    B0inv       <- diag(1 / scaling, P)
    
    
    # STEP 4: STORE POSTERIOR SAMPLES
    
    #only after burn-in period!
    if(irep > nburn){
      
      beta.store[irep-nburn,]  <- beta.draw
      sig2.store[irep-nburn,]  <- sig2.draw
      delta.store[irep-nburn,] <- delta.draw 
      
    }
    
    
    # STEP 4: PROGRESS
    # (in case you want to know how many iterations we already did)
    # if(irep%%(ntot/100)==0){  
    #   #Sys.sleep(0.0001)
    #   # update progress bar
    #   setTxtProgressBar(pb, irep)}
    
  } # MCMC LOOP END
  
  
  # -------------------------- #
  # --- MCMC CONVERGENCE ----- #
  # -------------------------- #
  
  if(fig == T) {
    #coefficient 1, true value and OLS estimate
    plot.ts(beta.store[,1])
    abline(a=beta.ols[1],b=0, col="green")
    
    #coefficient 2, true value and OLS estimate
    plot.ts(beta.store[,ncol(beta.store)])
    abline(a=beta.ols[ncol(beta.store)],b=0, col="green")
  }
  
  # -------------------------- #
  # --- PIP ANALYSIS --------- #
  # -------------------------- #
  
  # we can now look at the posterior inclusion probabilities for all coefficients
  # these are simply the average over all draws for delta
  pip.post  <- apply(delta.store, 2, mean)
  beta.post <- apply(beta.store, 2, mean)
  
  cbind(round(beta.post,2), round(pip.post,2))
  #View(cbind(colnames(X), round(beta.post,2), round(pip.post,2))[order(-pip.post),])
  
  
  # -------------------------- #
  # --- POSTERIOR ANALYSIS --- #
  # -------------------------- #
  
  # usually, the reader is not interested in you plotting full posterior distributions
  # we can use a few statistics to summarize the posterior distribution of any given parameter:
  
  # ... plot all posterior distributions
  if(fig == T) {
    # png("figures/bayes_ssvs/posterior_bayes_ssvs.png")
    par(mfrow=c(ceiling(P/2),2), mar=c(2,2,2,2))
    for(i in 1:P) {
      plot(density(beta.store[,i]), main = colnames(X)[i], 
           xlim = c(median(beta.store[,i])-3*sd(beta.store[,i]),median(beta.store[,i])+3*sd(beta.store[,i])))
      abline(v=mean(beta.store[,i]), col="red")
      abline(v=median(beta.store[,i]), col="green")
      abline(v=quantile(beta.store[,i], prob = 0.05), col="blue", lty=2)
      abline(v=quantile(beta.store[,i], prob = 0.95), col="blue", lty=2)
    }
    # dev.off()
  }
  
  if(fig == T) {
    # png("figures/bayes_ssvs/covar_box_bayes_ssvs.png")
    par(mfrow=c(ceiling(P/2),2),mar=c(2,2,2,2))
    for(i in 1:P) {
      boxplot(beta.store[,i], main = colnames(X)[i], horizontal = T)
    }
    # dev.off()
  }
  
  # # ... posterior means
  # apply(beta.store, 2, mean) #how does that compare to ols estimates?
  # 
  # # ... posterior medians
  # apply(beta.store, 2, median)
  # 
  # # ... posterior standard deviations
  # apply(beta.store, 2, sd)
  # 
  # # ... highest posterior density intervals / credible intervals (nice interpretation!)
  # apply(beta.store, 2, quantile, prob = c(0.05,0.95))
  
  
  # summerize output
  df <- data.frame(
    term = colnames(X), 
    mean = apply(beta.store, 2, mean), #beta mean
    pip = apply(delta.store,2, mean), #posterior inclusion probability
    median = apply(beta.store, 2, median), #beta median
    sd = apply(beta.store, 2, sd), #standard deviation
    
    cred0.05 = apply(beta.store, 2, quantile, prob = 0.05), #posterior density intervals, 0.05 and 0.95
    cred0.95 = apply(beta.store, 2, quantile, prob = 0.95),
    cred0.25 = apply(beta.store, 2, quantile, prob = 0.25), #posterior density intervals, 0.25 and 0.75
    cred0.75 = apply(beta.store, 2, quantile, prob = 0.75)
  )
  
  return(df)
}


## find optimal tuning parameters ----------------------------------------------
try(rji$sendMessage("BAY: Start with CV at the country level!"))

# set up parameter combinations
prm <- expand.grid(list(spike = seq(5, 15, by = 1)/1000,
                        slab = seq(5000, 15000, by = 1000)/1000))
prm$ix <- 1:nrow(prm)

prm$mod.r2 <- NA
prm$sev.r2 <- NA

# try all combinations and calculate R2
for(i in prm$ix) {
  print(paste0("BAY: ", round(sum(!is.na(prm$sev.r2))/nrow(prm)*100, 2), "%"))
  rji$sendMessage(paste0("BAY: ", round(sum(!is.na(prm$sev.r2))/nrow(prm)*100, 2), "%"))
  for (iso3 in unique(moddat$ISO3)){
    
    ix <- moddat$ISO3 != iso3
    dat <- moddat[ix, ]
    
    spike <- prm$spike[i]
    slab <- prm$slab[i]
    
    # moderate model
    y <- as.matrix(dat %>% select(fies.mod.logit))
    bigX <- cbind(matrix(1, nrow = nrow(dat), ncol = 1), as.matrix(dat %>% select(vars)))
    colnames(bigX) <- c("(Intercept)", vars)
    nburn <- 1000         #number of burn ins
    nsave <- 9000         #number of saved draws
    fig <- F
    mdf <- bayes.ssvs(y, bigX, nburn, nsave, spike, slab, fig)
    
    
    # severe model
    y <- as.matrix(dat %>% select(fies.sev.logit))
    bigX <- cbind(matrix(1, nrow = nrow(dat), ncol = 1), as.matrix(dat %>% select(vars)))
    colnames(bigX) <- c("(Intercept)", vars)
    nburn <- 1000        #number of burn ins
    nsave <- 9000        #number of saved draws
    fig <- F
    sdf <- bayes.ssvs(y, bigX, nburn, nsave, spike, slab, fig)
    
    # project and apply inverse logit transformation
    moddat$fies.mod.pred.cv[!ix] <- mdf[["mean"]][mdf$term == "(Intercept)"]
    moddat$fies.sev.pred.cv[!ix] <- sdf[["mean"]][sdf$term == "(Intercept)"]
    for (j in 2:nrow(mdf)){
      moddat$fies.mod.pred.cv[!ix] <- moddat$fies.mod.pred.cv[!ix] + moddat[!ix, mdf$term[j]]*mdf[["mean"]][j]
      moddat$fies.sev.pred.cv[!ix] <- moddat$fies.sev.pred.cv[!ix] + moddat[!ix, sdf$term[j]]*sdf[["mean"]][j]
    }
    moddat$fies.mod.pred.cv[!ix] <- inv.logit(moddat$fies.mod.pred.cv[!ix])
    moddat$fies.sev.pred.cv[!ix] <- inv.logit(moddat$fies.sev.pred.cv[!ix])
    
  }
  
  # calculate R2
  prm$mod.r2[i] <- fr2(moddat$fies.mod, moddat$fies.mod.pred.cv)
  prm$sev.r2[i] <- fr2(moddat$fies.sev, moddat$fies.sev.pred.cv)
  
}

# write.csv(prm, "prm-bayes_ssvs.csv", row.names = F)

try(rji$sendMessage("BAY: Done!"))
try(rji$sendMessage("BAY: Done with model-bayes_ssvs_1_wigeo.R!"))
