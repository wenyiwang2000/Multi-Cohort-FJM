predict.MMFPCA2 <- function(object,newdata,calculate.scores=T, ...){
  
  ## check inputs
  if(class(object)!="mface.sparse") stop("'fit' has to be a mface.sparse object")
  p <- length(newdata) 
  Bnew <- vector("list", p)
  var.error.pred <- vector("list", p)
  Bi <- vector("list", p)
  Bi_pred <- vector("list", p)
  B <- vector("list", p)
  mu.pred <- vector("list", p)
  uy.pred <- vector("list", p)
  uscores.pred <- vector("list", p)
  ucov.pred <- vector("list", p)
  for(k in 1:p){
    temp <- newdata[[k]]
    pred <- predict.face.sparse.inner(object$fit[[k]], temp)
    Bnew[[k]] <- object$fit[[k]]$Bnew %*% object$G_invhalf
    Bi[[k]] <- lapply(1:length(pred$Bi),function(x) pred$Bi[[x]]%*%object$G_invhalf)
    Bi_pred[[k]] <- lapply(1:length(pred$Bi_pred),function(x) pred$Bi_pred[[x]]%*%object$G_invhalf)
    B[[k]] <- pred$B %*% object$G_invhalf
    
    var.error.pred[[k]] <- pred$var.error.pred
    mu.pred[[k]] <- pred$mu.pred
    uy.pred[[k]] <- pred$y.pred
    uscores.pred[[k]] <- pred$scores$scores
    ucov.pred[[k]] <- pred$cov.pred
  }
  
  subj.pred = lapply(1:p, function(x){newdata[[x]]$subj})
  subj_unique.pred = unique(subj.pred[[1]])
  y.pred = lapply(1:p, function(x){newdata[[x]][, c(-1,-2)]})
  u.pred = lapply(1:p, function(x){newdata[[x]][, c(-1,-2)]})  # <--- u
  w.pred = lapply(1:p, function(x){newdata[[x]][, c(-1,-2)]})   # <--- w
  se.pred = lapply(1:p, function(x){0*y.pred[[x]]})
  
  Theta = object$Theta
  Theta_U = object$Theta_U
  Theta_W = object$Theta_W
  Theta0 = object$Theta0
  Theta1 = object$Theta1
  npc = object$npc
  npc0 = object$npc0
  npc1 = object$npc1
  b = object$Beta   # <------
  Bbasis = unlist(Bi)
  
  
  Bnew <- do.call(bdiag, Bnew)
  Bi <- unlist(Bi)
  Bi_pred <- unlist(Bi_pred)
  B <- do.call(bdiag, B) 
  cov.pred <- matrix(0,length(unlist(mu.pred)),length(unlist(mu.pred)))
  
  scores = list(subj=subj_unique.pred,
                scores = matrix(NA,nrow=length(subj_unique.pred),ncol=(npc0+p*npc1)),
                u = matrix(NA,nrow=length(subj_unique.pred),ncol=(nrow(Theta0)+p*nrow(Theta1))),
                Sig = list())
  
  scores0 = list(subj=subj_unique.pred,
                 scores = matrix(NA,nrow=length(subj_unique.pred),ncol=npc0),
                 u = matrix(NA,nrow=length(subj_unique.pred),ncol=nrow(Theta0)),
                 Sig = list())
  
  
  for(i in 1:length(subj_unique.pred)){
    
    sel.pred = lapply(1:p, function(x){which(subj.pred[[x]]==subj_unique.pred[i])})
    pred.points <- lapply(1:p, function(x){newdata[[x]]$argvals[sel.pred[[x]]]})
    mu.predi <- lapply(1:p, function(x){mu.pred[[x]][sel.pred[[x]]]})
    var.error.predi <-  lapply(1:p, function(x){var.error.pred[[x]][sel.pred[[x]]]})
    
    y.predi = lapply(1:p, function(x){y.pred[[x]][sel.pred[[x]]] - mu.predi[[x]]})
    sel.pred.obs = lapply(1:p, function(x){which(!is.na(y.predi[[x]]))})
    obs.points <- lapply(1:p, function(x){pred.points[[x]][sel.pred.obs[[x]]]})
    
    if(!is.null(obs.points)){
      
      B3i.pred = Bi_pred[[i]]
      B3i = Bi[[i]]
      for(j in 1:(p-1)){
        B3i.pred = bdiag(B3i.pred, Bi_pred[[i+j*length(subj_unique.pred)]])
        B3i = bdiag(B3i, Bi[[i+j*length(subj_unique.pred)]])
      }
      Chati = tcrossprod(B3i%*%Theta,B3i)
      
    
      temp = unlist(lapply(1:p, function(x){var.error.predi[[x]][sel.pred.obs[[x]]]}))
      Ri = diag(temp)
      Vi.inv = as.matrix(solve(Chati + Ri))
      Vi.pred = as.matrix(tcrossprod(B3i.pred%*%Theta,B3i.pred))
      
      Hi = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i)%*%Vi.inv)
      Hi_u = as.matrix(B3i.pred%*%tcrossprod(Theta_U,B3i)%*%Vi.inv) # <---
      Hi_w = as.matrix(B3i.pred%*%tcrossprod(Theta_W,B3i)%*%Vi.inv) # <---
      
      
      Theta0Bi = lapply(1:p, function(x){tcrossprod(Theta0,Bbasis[[x]])})
      Theta0Bimat =  as.matrix(do.call("cbind",lapply(1:p, function(x){b[x]*Theta0Bi[[x]]})))
      
      vy.predi = unlist(lapply(1:p, function(x){y.predi[[x]][sel.pred.obs[[x]]]}))
      ui0 = Theta0Bimat%*%Vi.inv %*% vy.predi
      Si0 = Theta0Bimat%*%Vi.inv%*%t(Theta0Bimat)
      
      Theta1Bi = lapply(1:p, function(x){tcrossprod(Theta1,Bbasis[[x]])})
      l0 = c(0,cumsum(mapply(function(sel.pred.obs){length(sel.pred.obs)}, sel.pred.obs=sel.pred.obs)))
      Theta1Bimat = do.call("rbind",lapply(1:p, function(j){
        temp = matrix(0, nrow=nrow(Theta1Bi[[j]]), ncol=l0[p+1])
        idx = (l0[j]+1):l0[j+1]
        temp[ ,idx] = as.matrix(b[j]*Theta1Bi[[j]])
        temp }) )
      ThetaBmat = rbind(Theta0Bimat, Theta1Bimat)
      
      
      ui = ThetaBmat%*%Vi.inv %*% vy.predi
      Si = ThetaBmat%*%Vi.inv%*%t(ThetaBmat)
      scores$u[i,] = as.vector(ui)
      
      
      ll = c(0,cumsum( mapply(function(sel.pred){length(sel.pred)}, sel.pred=sel.pred)))
      y.temp = as.numeric(Hi%*%vy.predi) + unlist(mu.predi)
      u.temp = as.numeric(Hi_u%*%vy.predi) 
      w.temp = as.numeric(Hi_w%*%vy.predi) 
      # only for x_ij
      temp = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i))
      cov.pred = Vi.pred - temp%*%Vi.inv%*%t(temp)
      se.temp = sqrt(diag(cov.pred))
      for (k in 1:p) {
        y.pred[[k]][sel.pred[[k]]] = y.temp[(ll[k]+1):ll[k+1]]
        u.pred[[k]][sel.pred[[k]]] = u.temp[(ll[k]+1):ll[k+1]]
        w.pred[[k]][sel.pred[[k]]] = w.temp[(ll[k]+1):ll[k+1]]
        se.pred[[k]][sel.pred[[k]]] = se.temp[(ll[k]+1):ll[k+1]]
      }
        
      
      ## predict scores
      if(calculate.scores==TRUE){ 
        U = bdiag(matrix(t(object$U0),nrow=npc0),
                  bdiag(lapply(1:p, function(j){matrix(t(object$U1),nrow=npc1)})))
        temp = U%*%ui
        temp = as.matrix(temp)
        scores$scores[i,1:(npc0+p*npc1)] = temp[,1]
        scores$Sig[[i]] = diag(c(object$eigenvalues0, rep(object$eigenvalues1,p))) - 
          U %*% as.matrix(Si) %*% t(U)
        
        temp0 = matrix(t(object$U0),nrow=npc0)%*%ui0
        temp0 = as.matrix(temp0)
        scores0$scores[i,1:npc0] = temp0[,1]
        scores0$Sig[[i]] = diag(object$eigenvalues0,nrow=npc0) - 
          (matrix(t(object$U0),nrow=npc0)) %*% as.matrix(Si0) %*% matrix(object$U0,ncol=npc0)
      }
    }
  }
  
  Chat.pred = as.matrix(tcrossprod(B%*%Matrix(Theta),B))
  
  
  return(list(object=object,newdata=newdata,B=B,
              y.pred = y.pred,mu.pred=mu.pred,
              u.pred = u.pred, w.pred = w.pred,
              var.error.pred=var.error.pred,
              scores = scores, scores0 = scores0, 
              cov.pred = cov.pred,
              se.pred = se.pred, uscores.pred = uscores.pred,
              Chat.diag.pred = diag(Chat.pred),
              ucov.pred = ucov.pred, uy.pred = uy.pred))  
}
