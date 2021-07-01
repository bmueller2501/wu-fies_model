# functions for model validation

fmae <- function(y, yhat) {mean(abs(y - yhat))}
frmse <- function(y, yhat) {sqrt((1/length(y))*sum((y - yhat)^2))}
fr2 <- function(y, yhat) {(1-(sum((y-yhat)^2)/sum((y-mean(yhat))^2)))}
