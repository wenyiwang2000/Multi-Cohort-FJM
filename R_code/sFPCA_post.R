sFPCA_post <- function(data, S0_fit, S1_fit, 
                       mfaces_fit, Beta, pve = 0.99){
  
  p <- length(Beta)
  S0 = S0_fit$S
  S1 = S1_fit$S
  
  object = mfaces_fit
  tnew = object$argvals.new
  Bnew <- vector("list", p)
  B <- vector("list", p)
  for(k in 1:p){
    Bnew[[k]] <- object$fit[[k]]$Bnew %*% object$G_invhalf
    B[[k]] <- object$fit[[k]]$B %*% object$G_invhalf
  }
  Bbasis = object$fit[[k]]$Bnew %*% object$G_invhalf
  Bn = Bnew[[1]] / nrow(object$var.error.new)
  Theta0 = t(Bn) %*% S0 %*% Bn
  Theta1 = t(Bn) %*% S1 %*% Bn
  
  Theta_all <- do.call("rbind", lapply(1:p, function(j){
    do.call("cbind",lapply(1:p, function(k){ Beta[j]*Beta[k]*(Theta0 + (j==k)*Theta1) }))
  }))
  Bnew <- do.call(bdiag, Bnew) 
  Chat <- Bnew%*%Theta_all%*%t(Bnew)
  
  
  Theta_W <- bdiag(lapply(1:p, function(j){ Beta[j]^2*Theta1 }))
  Theta_U = Theta_all - Theta_W
  
  Chat_diag = as.vector(diag(Chat))  
  Cor = diag(1/sqrt(Chat_diag))%*%Chat%*%diag(1/sqrt(Chat_diag))
  
  B <- do.call(bdiag, B)
  #Chat.diag.pred = diag(as.matrix(tcrossprod(B%*%Matrix(Theta_all),B)))
  
  Eig <- eigen(Theta_all)
  npc <- which.max(cumsum(Eig$values)/sum(Eig$values)>pve)[1]
  
  eigenfunctions = matrix(Bnew%*%Eig$vectors[,1:min(npc, p*length(tnew))],
                          ncol=min(npc, p*length(tnew)))
  eigenvalues = Eig$values[1:min(npc, p*length(tnew))]
  
  Eig0 <- eigen(Theta0)
  npc0 <- which.max(cumsum(Eig0$values)/sum(Eig0$values)>pve)[1]
  
  eigenfunctions0 = matrix(Bbasis%*%Eig0$vectors[,1:min(npc0, length(tnew))],
                           ncol=min(npc0, length(tnew)))
  eigenvalues0 = Eig0$values[1:min(npc0, length(tnew))]
  
  Eig1 <- eigen(Theta1)
  npc1 <- which.max(cumsum(Eig1$values)/sum(Eig1$values)>pve)[1]
  
  eigenfunctions1 = matrix(Bbasis%*%Eig1$vectors[,1:min(npc1, length(tnew))],
                           ncol=min(npc1, length(tnew)))
  eigenvalues1 = Eig1$values[1:min(npc1, length(tnew))]
  
  sigma2 <- rep(NaN, p)
  for(j in 1:p){
    sigma2[j] = object$fit[[j]]$sigma2 
  }
  var.error.hat <- object$var.error.hat
  var.error.new <- object$var.error.new
  
  res <- list(fit = object$fit, Theta = Theta_all, 
              Theta0 = Theta0, Theta1 = Theta1,
              Theta_U = Theta_U, Theta_W = Theta_W,
              Chat.new = Chat, Cor.new = Cor, 
              npc = npc, eigenfunctions = eigenfunctions, 
              eigenvalues = eigenvalues, U = Eig$vectors[,1:npc],
              npc0 = npc0, eigenfunctions0 = eigenfunctions0, 
              eigenvalues0 = eigenvalues0, U0 = Eig0$vectors[,1:npc0],
              npc1 = npc1, eigenfunctions1 = eigenfunctions1, 
              eigenvalues1 = eigenvalues1, U1 = Eig1$vectors[,1:npc1],
              var.error.hat = var.error.hat, 
              var.error.new = var.error.new, pve=pve, 
              G_invhalf = object$G_invhalf, sigma2 = sigma2, Beta = Beta)
  class(res) <- "mface.sparse"
  return(res)
}