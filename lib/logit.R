# logit

logit <- function(x){
  log(x/(1-x))
}

inv.logit <- function(x){
  1/(1 + exp(-x))
}
