## ------------------------------------------------- ##
## Visualize Variable Importance and Number of Trees ##
## --------------------------------------------------##

s <- T #save plots?

## variable importance ---------------------------------------------------------

# change variable names
vars_full <- c("Stunting", "Wasting", "Mean Years of Schooling", "Topographic Ruggedness",
               "GDP Per Capita", "Gini Coefficient", "Malaria Mortality Rate", "Mean Annual Precipitation",
               "Poverty Headcount Index", "Mean Temperature", "Urban Percentage", "Water Scarcity") #has to be in the same order as "vars"

v.mod <- vimp(rf.mod, importance='permute')
v.sev <- vimp(rf.sev, importance='permute')

df <- data.frame(var=c(names(v.mod$importance), names(v.sev$importance)),
                 val=c(v.mod$importance, v.sev$importance),
                 mod=rep(c('Moderate', 'Severe'), each=12),
                 lab=rep(vars_full, 2))

df$lab <- factor(df$lab, levels=vars_full[order(v.sev$importance)])

convertToOutcomeScale <- function(x){
  # determine how an error of x on a logit transformed scale
  # affects the mean absolute error after being transformed back to [0,1]
  mean(abs(moddat$fies.mod - inv.logit(logit(moddat$fies.mod) + x)))
}

df$val <- sapply(df$val, convertToOutcomeScale)

ggplot(df) + 
  geom_bar(aes(x=lab, y=val), stat='identity') +
  coord_flip() + 
  facet_grid(. ~ mod) + 
  theme_bw() + 
  labs(x='', y='Increase in Model Error When Variable is Permutated')
if(s) ggsave(filename = file.path(path_viz, paste0("rf-vimp.pdf")), width=7, height=3.75)


## error rate and number of trees  ------------------------------------------------
erdf <- data.frame(var=c(v.mod$err.rate, v.sev$err.rate),
                   trees=c(1:10000, 1:10000),
                   mod=rep(c('Moderate', 'Severe'), each=10000)) %>%
  na.omit

ggplot(erdf) + 
  geom_line(aes(x=trees, y=var)) + 
  facet_wrap(. ~ mod, scales='free_y') + 
  theme_bw() +
  labs(x='Number of Trees', y='Error Rate')
if(s) ggsave(filename = file.path(path_viz, paste0("rf-errortrees.pdf")), width=7, height=3.75)






