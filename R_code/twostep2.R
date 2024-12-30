
#library(mfaces)
#source("mface.sparse2.R")
#source("sFPCA_post.R")
#source("predict.MMFPCA2.R")
#source("predict.face.sparse.inner_pre.R")


solveB <- function(x,y){
  sum(x*y)/sum(x^2)
}

## post-process covariance
sFPCA_fit <- function(S_fit, pve = 0.99){
  eigS = eigen(S_fit)
  sel = (eigS$values > 0)
  eigvalues = eigS$values[sel]
  eigvectors = eigS$vectors[,sel]
  npc = which.max(cumsum(eigvalues) / sum(eigvalues) > pve)
  eigvalues = eigvalues[1:npc] # <---
  eigvectors = eigvectors[,1:npc] # <---
  if (npc==1){
    S_fit = eigvectors %*% t(eigvectors) * eigvalues
  }else{
    S_fit = eigvectors %*% diag(eigvalues) %*% t(eigvectors)
  }
  return(list(S = S_fit, eigvalues = eigvalues, eigvectors = eigvectors, 
              pve = pve, npc = npc))
}


twostep2 <- function(longdata, knots, knots.option="equally-spaced", argvals.new){
  
  
  J <- length(longdata)
  m <- length(argvals.new)
  n <- length(unique(longdata[[1]]$subj))
  
  
  fit <- mface.sparse2(data=longdata, argvals.new=argvals.new, knots=knots, 
                       knots.option=knots.option, newdata=longdata)
  
  # Get covariance matrix
  sel <- 1:m
  C <- as.matrix(fit$Chat.new)
  Cov <- lapply(1:J, function(j){
    lapply(1:J, function(k){C[sel+(j-1)*m,sel+(k-1)*m]})
  }) 
  
  # Estimate scaling parameters
  gamma <- rep(1,J)
  for(j in 2:J){
    idx <- (1:J)[c(-1,-j)]
    x <- do.call("rbind",Cov[[1]][idx])
    y <- do.call("rbind",Cov[[j]][idx])
    gamma[j] <-  solveB(x=x,y=y)
  }
  B <- unlist(lapply(1:(J-1),function(j){ gamma[j]*gamma[(j+1):J] })) 
  
  # Solve S0
  S0 <- matrix(NaN,m,m)
  S1 <- S0
  for(j in 1:m){
    for(k in 1:m){
      y <- unlist(lapply(1:(J-1), function(jj){
        sapply((jj+1):J, function(kk){(Cov[[jj]][[kk]])[j,k]}) })) 
      S0[j,k] <- solveB(x=B, y=y)
    }
  }
  S0 = forceSymmetric(S0)
  S0_fit = sFPCA_fit(S0)
  S0 = S0_fit$S
  # Solve S1 
  x <- gamma^2
  for(j in 1:m){
    for(k in 1:m){
      y <- sapply(1:J, function(jj){ (Cov[[jj]][[jj]])[j,k] - gamma[jj]^2*S0[j,k] })
      S1[j,k] <- solveB(x=x, y=y)
    }
  }
  S1_fit = sFPCA_fit(S1)
  S1 = S1_fit$S
  
  
  ## The post-process MFMM estimates
  face_fit <- sFPCA_post(data = longdata, S0_fit = S0_fit, S1_fit = S1_fit, 
                         mfaces_fit = fit, Beta = gamma, pve = 0.99)
  
  
  ## The score prediction by conditional expectation method
  xi <- matrix(NA, n, (face_fit$npc0+J*face_fit$npc1))
  xi0 <- matrix(NA, n, face_fit$npc0)
  xi_cov <- list()
  xi0_cov <- list()
  for(i in 1:n){
    
    dat_i_pred <- lapply(1:J, function(x){
      sel <- which(longdata[[x]]$subj==i)
      dat_i <- longdata[[x]][sel,]
      dat_i_pre <- data.frame(subj=rep(dat_i$subj[1], nrow(dat_i)+length(argvals.new)),
                               argvals = c(rep(NA,nrow(dat_i)),argvals.new), 
                               y=rep(NA, nr=nrow(dat_i)+length(argvals.new)))
      dat_i_pre[1:nrow(dat_i), ] <- dat_i
      dat_i_pre
    })
    pred <- predict.MMFPCA2(object=face_fit, newdata=dat_i_pred)
    
    xi[i,] <- pred$scores$scores
    xi0[i,] <- pred$scores0$scores
    
    xi_cov[[i]] <- pred$scores$Sig[[1]]
    xi0_cov[[i]] <- pred$scores0$Sig[[1]]
  }
  
  face_fit$scores$scores <- xi
  face_fit$scores0$scores <- xi0
  face_fit$scores$Sig <- xi_cov
  face_fit$scores0$Sig <- xi0_cov
  
  return(face_fit)
}









