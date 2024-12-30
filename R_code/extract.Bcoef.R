extract.Bcoef <- function(fit,data){
  
  stopifnot(class(fit)[[1]] =="gam")
  model.formula <- as.formula(fit$formula[[3]])
  sm <- smoothCon(model.formula, data = data)[[1]]
  
  qrc <- qr(t(sm$C))
  K <- length(fit$coefficients)
  beta <- c(0,as.vector(fit$coefficients)[2:K])
  beta <- qr.qty(qrc, beta)
  delta <- beta + fit$coefficients[1] 
  
  return(delta)
  
}