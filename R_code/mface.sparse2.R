
# 1
mface.sparse2 <- function(data, newdata = NULL, center = TRUE, argvals.new = NULL, 
                         knots = 7, knots.option="equally-spaced", 
                         p = 3, m = 2, lambda = NULL, lambda_mean = NULL, 
                         lambda_bps = NULL, 
                         search.length = 14, lower = -3, upper = 10, 
                         calculate.scores = FALSE, pve = 0.99)
{
  require(face)
  require(matrixcalc)
  require(Matrix)
  require(splines)
  require(mgcv)
  
  #########################
  ####step 0: read in data
  #########################
  
  tnew <- argvals.new
  p.m = p ## this is splines degree
  p <- length(data) ## dimension of the multivariate functional data
  
  object <- vector("list", p)
  Theta0 <- vector("list", p)
  Bnew <- vector("list", p)
  C <- vector("integer", p)
  # sigma2 <- vector("numeric", p)
  var.error.pred <- vector("list", p)
  var.error.hat <- vector("list", p)
  var.error.new <- vector("list", p)
  mu.pred <- vector("list", p)
  Bi <- vector("list", p)
  Bi_pred <- vector("list", p)
  B <- vector("list", p)
  
  W <- data
  
  #############################################
  ####step 1: univariate fpca for each response
  #############################################
  
  for(k in 1:p){
    temp <- data[[k]]
    temp_new = NULL
    if(!is.null(newdata)) {
      temp_new <- newdata[[k]]
    }
    object[[k]] <- face.sparse.inner(data = temp, newdata = temp_new, center = center, 
                                     argvals.new = tnew, knots = knots, 
                                     knots.option = knots.option, p = p.m, m = m, 
                                     lambda = lambda, lambda_mean = lambda_mean, 
                                     search.length = search.length, lower = lower, upper = upper, 
                                     calculate.scores = calculate.scores, pve = pve)
    Theta0[[k]] <- object[[k]]$Theta0
    if(k==1){
      G <- crossprod(object[[k]]$Bnew) / nrow(object[[k]]$Bnew)
      eig_G <- eigen(G, symmetric = T)
      G_half <- eig_G$vectors %*% diag(sqrt(eig_G$values)) %*% t(eig_G$vectors)
      G_invhalf <- eig_G$vectors %*% diag(1/sqrt(eig_G$values)) %*% t(eig_G$vectors)
    }
    Bnew[[k]] <- object[[k]]$Bnew %*% G_invhalf
    C[k] <- ncol(object[[k]]$Bnew)
    # sigma2[k] <- object[[k]]$sigma2
    var.error.hat[[k]] <- object[[k]]$var.error.hat
    var.error.new[[k]] <- object[[k]]$var.error.new
    if(!is.null(newdata)) {
      mu.pred[[k]] <- object[[k]]$mu.pred
      var.error.pred[[k]] <- object[[k]]$var.error.pred
      Bi[[k]] <- lapply(1:length(object[[k]]$Bi),function(x) object[[k]]$Bi[[x]]%*%G_invhalf)
      Bi_pred[[k]] <- lapply(1:length(object[[k]]$Bi_pred),function(x) object[[k]]$Bi_pred[[x]]%*%G_invhalf)
      B[[k]] <- object[[k]]$B %*% G_invhalf
    }
    if(center){
      W[[k]][, "y"] <- W[[k]][, "y"] - object[[k]]$mu.hat # demean each response
    }
  }
  names(object) <- names(data)
  
  ###############################################
  ####step 2: compute empirical cross-covariances
  ###############################################
  
  xprod <- vector("list", p*(p-1)/2)
  N2 <- vector("list", p*(p-1)/2)
  usubj <- unique(W[[1]]$subj)
  for(i in 1:length(usubj)){
    crossCov <- vector("list", p*(p-1)/2)
    ind <- 1
    for(k in 1:(p-1))
      for(j in (k+1):p){
        temp_k <- W[[k]][which(W[[k]]$subj == usubj[i]),]
        temp_j <- W[[j]][which(W[[j]]$subj == usubj[i]),]
        N2[[ind]][i] <- nrow(temp_k) * nrow(temp_j)
        crossCov[[ind]] <- as.vector(temp_k[,"y"]%*%t(temp_j[,"y"]))
        xprod[[ind]][[i]] <- data.frame(arg1 = rep(temp_k$argvals, nrow(temp_j)),
                                        arg2 = rep(temp_j$argvals, each = nrow(temp_k)),
                                        crossCov = crossCov[[ind]])  
        ind <- ind + 1
      }
  }
  
  ####################################
  ####step 3: smooth cross-covariance
  ####################################
  
  # smooth.cov <- vector("list",p*(p-1)/2)
  Theta <- vector("list", p*(p-1)/2)
  bps.lambda <- matrix(NaN, p*(p-1)/2,2)
  for(k in 1:(p*(p-1)/2)) {
    temp <- do.call(rbind, xprod[[k]])
    fit <-  bps(temp$crossCov, temp$arg1, temp$arg2, knots=knots, p = p.m, m = m, 
                lower = lower, upper = upper, search.length = search.length,
                lambda = lambda_bps, knots.option=knots.option, N2 = N2[[k]])
    Theta[[k]] <- fit$Theta
    bps.lambda[k,] <- fit$lambda
  }
  
  Theta_all <- do.call(bdiag, Theta0)
  C_idx <- c(0, cumsum(C))
  ind <- 1
  for(k in 1:(p-1)){
    for(j in (k+1):p){
      sel1 <- C_idx[k] + 1:C[k]
      sel2 <- C_idx[j] + 1:C[j]
      Theta_all[sel1,sel2] <- Theta[[ind]]
      Theta_all[sel2,sel1] <- t(Theta[[ind]])
      ind <- ind + 1
    }
  }
  Theta_all <- as.matrix(Theta_all)
  
  
  G_halfmat <-do.call(bdiag, lapply(1:p, function(x) G_half))
  Theta_all <- G_halfmat %*% Theta_all %*% G_halfmat
  
  ####################################################
  ####step 4: calculate estimated covariance function
  ####################################################
  
  ##make covariance matrix psd
  # Eig <- eigen(cov_all)
  Eig <- eigen(Theta_all)
  sel <- which(Eig$values>0)
  if(length(sel)>1) Theta_all <- matrix.multiply(Eig$vectors[,sel], Eig$values[sel])%*%t(Eig$vector[,sel])
  if(length(sel)==1) Theta_all <- Eig$vectors[,sel]%*%t(Eig$vector[,sel])*Eig$values[sel]
  
  Bnew <- do.call(bdiag, Bnew)#%*%G_invhalfmat
  Chat <- Bnew%*%Theta_all%*%t(Bnew)
  
  Chat_diag = as.vector(diag(Chat))  
  Cor = diag(1/sqrt(Chat_diag))%*%Chat%*%diag(1/sqrt(Chat_diag))
  
  npc <- which.max(cumsum(Eig$values[sel])/sum(Eig$values[sel])>pve)[1]
  
  eigenfunctions = matrix(Bnew%*%Eig$vectors[,1:min(npc, p*length(tnew))],
                          ncol=min(npc, p*length(tnew))) #/ sqrt((max(tnew) - min(tnew))/(length(tnew) - 1))
  eigenvalues = Eig$values[1:min(npc, p*length(tnew))] #* ((max(tnew) - min(tnew))/(length(tnew) - 1))
  
  ##############################
  ####step 5: calculate variance
  ##############################
  
  # var.error.hat <- do.call(cbind, var.error.hat)
  var.error.new <- do.call(cbind, var.error.new)
  
  Chat.raw.new = as.matrix(tcrossprod(Bnew%*%Matrix(Theta_all),Bnew)) + diag(var.error.new) 
  Chat.raw.diag.new = as.vector(diag(Chat.raw.new)) 
  Cor.raw.new = diag(1/sqrt(Chat.raw.diag.new))%*%Chat.raw.new%*%diag(1/sqrt(Chat.raw.diag.new))
  
  #######################
  ####step 6: prediction
  #######################
  
  if(!is.null(newdata)){
    
    subj.pred = lapply(newdata, function(x) {x$subj})
    subj_unique.pred = unique(subj.pred[[1]])
    y.pred = lapply(newdata, function(x) {x[, "y"]})
    se.pred = lapply(y.pred, function(x) {0*x})
    
    Bi <- unlist(Bi)
    Bi_pred <- unlist(Bi_pred)
    B <- do.call(bdiag, B) 
    
    scores = list(subj=subj_unique.pred,
                  scores = matrix(NA,nrow=length(subj_unique.pred),ncol=npc),
                  u = matrix(NA,nrow=length(subj_unique.pred),ncol=nrow(Theta_all)),
                  Sig = list())
    
    for(i in 1:length(subj_unique.pred)){
      sel.pred = lapply(subj.pred, function(x){
        which(x==subj_unique.pred[i])})
      pred.points <- mapply(function(newdata, sel.pred) {
        newdata$argvals[sel.pred]}, newdata = newdata, sel.pred = sel.pred, 
        SIMPLIFY = F)
      mu.predi <- mapply(function(mu.pred, sel.pred) {
        mu.pred[sel.pred]}, mu.pred = mu.pred, sel.pred = sel.pred, 
        SIMPLIFY = F)
      var.error.predi <- mapply(function(var.error.pred, sel.pred) {
        var.error.pred[sel.pred]}, var.error.pred = var.error.pred,
        sel.pred = sel.pred, SIMPLIFY = F)
      
      y.predi = mapply(function(y.pred, sel.pred, mu.predi) {
        y.pred[sel.pred] - mu.predi}, y.pred = y.pred, 
        sel.pred = sel.pred, mu.predi = mu.predi, SIMPLIFY = F)
      sel.pred.obs = lapply(y.predi, function(x){
        unname(which(!is.na(x)))})
      obs.points <- mapply(function(pred.points, sel.pred.obs) {
        pred.points[sel.pred.obs]}, pred.points = pred.points, 
        sel.pred.obs = sel.pred.obs, SIMPLIFY = F)
      
      if(!any(unlist(lapply(obs.points, function(x){is.null(x)}))==T)){
        var_sel <- mapply(function(var.error.predi, sel.pred.obs) {
          var.error.predi[sel.pred.obs]
        }, var.error.predi = var.error.predi, 
        sel.pred.obs = sel.pred.obs, SIMPLIFY = F)
        var = unlist(lapply(var_sel, function(x) {mean(x)}))
        len_sel.pred.obs <- unlist(lapply(sel.pred.obs, function(x){
          length(x)}))
        if(sum(var)==0&sum(len_sel.pred.obs) < npc)
          stop("Measurement error estimated to be zero and there are fewer observed points thans PCs; scores
              cannot be estimated.")
        
        B3i.pred = Bi_pred[[i]]
        B3i = Bi[[i]]
        for(j in 1:(p-1)){
          B3i.pred = bdiag(B3i.pred, Bi_pred[[i+j*length(subj_unique.pred)]])
          B3i = bdiag(B3i, Bi[[i+j*length(subj_unique.pred)]])
        }
        Chati = tcrossprod(B3i%*%Theta_all,B3i)
        
        if(sum(len_sel.pred.obs)==1) Ri = unlist(var_sel)
        if(sum(len_sel.pred.obs)>1)  Ri = diag(unlist(var_sel))
        Vi.inv = as.matrix(solve(Chati + Ri))
        Vi.pred = as.matrix(tcrossprod(B3i.pred%*%Theta_all,B3i.pred))
        Hi = as.matrix(B3i.pred%*%tcrossprod(Theta_all,B3i)%*%Vi.inv)
        y.predi_sel = unname(unlist(mapply(function(y.predi, sel.pred.obs){
          y.predi[sel.pred.obs]}, y.predi = y.predi, 
          sel.pred.obs = sel.pred.obs, SIMPLIFY = F)))
        ui =tcrossprod(Theta_all,B3i)%*%Vi.inv %*% y.predi_sel
        Si = tcrossprod(Theta_all,B3i)%*%Vi.inv%*%t(tcrossprod(Theta_all,B3i))
        scores$u[i,] = as.vector(ui)
        temp0 = as.numeric(Hi%*%y.predi_sel) + unlist(mu.predi)
        len_sel.pred = unlist(lapply(sel.pred, function(x){length(x)}))
        ind_sel.pred = unname(cumsum(len_sel.pred))
        ind_temp <- vector("list", p)
        for(j in 1:p){
          if(j==1){
            ind_temp[[j]] = 1:ind_sel.pred[1]
          }else{
            ind_temp[[j]] = (1+ind_sel.pred[j-1]):ind_sel.pred[j]
          }
        }
        temp = as.matrix(B3i.pred%*%tcrossprod(Theta_all,B3i))
        for(j in 1:p){
          y.pred[[j]][sel.pred[[j]]] <- temp0[ind_temp[[j]]]
          if(sum(len_sel.pred.obs) >1){
            se.pred[[j]][sel.pred[[j]]] = sqrt(diag(Vi.pred - temp%*%Vi.inv%*%t(temp)))[ind_temp[[j]]]
          }
          if(sum(len_sel.pred.obs) ==1){
            se.pred[[j]][sel.pred[[j]]] = sqrt(Vi.pred[1,1] - Vi.inv[1,1]*temp%*%t(temp))
          }
        }
        
        ## predict scores
        if(calculate.scores==TRUE){ 
          temp = matrix(t(Eig$vectors[,1:npc]),nrow=npc)%*%ui
          temp = as.matrix(temp)
          scores$scores[i,1:npc] = temp[,1]
          scores$Sig[[i]] = diag(eigenvalues) - 
            (matrix(t(Eig$vectors[,1:npc]),nrow=npc)) %*% as.matrix(Si) %*% matrix(Eig$vectors[,1:npc],ncol=npc)
        }
      }
    }
    
    #Chat.diag.pred = diag(as.matrix(tcrossprod(B%*%Matrix(Theta_all),B)))
    Chat.diag.pred = NULL
  }
  if(is.null(newdata)){
    y.pred=NULL
    mu.pred = NULL
    var.error.pred = NULL
    Chat.diag.pred = NULL
    # cov.pred = NULL
    se.pred = NULL
    scores=NULL
    
  }
  
  res <- list(fit = object, Theta = Theta_all, 
              Chat.new = Chat, Cor.new = Cor, 
              npc = npc, eigenfunctions = eigenfunctions, 
              eigenvalues = eigenvalues, 
              var.error.hat = var.error.hat, 
              var.error.new = var.error.new, Chat.raw.diag.new = Chat.raw.diag.new, 
              Cor.raw.new = Cor.raw.new,
              y.pred = y.pred, mu.pred = mu.pred, var.error.pred = var.error.pred, 
              Chat.diag.pred = Chat.diag.pred, 
              se.pred = se.pred, scores = scores,
              G_invhalf = G_invhalf, bps.lambda = bps.lambda, 
              U = Eig$vectors[,1:npc], 
              argvals.new = argvals.new,
              center=center, knots=knots, knots.option = knots.option,
              p = p.m, m = m,
              lower=lower,upper=upper,search.length=search.length,
              calculate.scores=calculate.scores,pve=pve)
  class(res) <- "mface.sparse"
  return(res)
  
}





# 2
face.sparse.inner <- function(data, newdata = NULL, W = NULL,
                              center=TRUE,argvals.new=NULL,
                              knots=7, knots.option="equally-spaced",
                              p=3,m=2,lambda=NULL,lambda_mean=NULL,
                              search.length=14,
                              lower=-3,upper=10, 
                              calculate.scores=FALSE,pve=0.99){
  
  #########################
  ####step 0: read in data
  #########################
  check.data(data)
  if(!is.null(newdata)){ check.data(newdata,type="predict")}
  
  y <- data$y
  t <- data$argvals
  subj <- data$subj
  tnew <- argvals.new
  if(is.null(tnew)) tnew <- seq(min(t),max(t),length=100)
  
  fit_mean <- NULL
  
  knots.initial <- knots
  #########################
  ####step 1: demean
  #########################
  r <- y
  mu.new <- rep(0,length(tnew))
  if(center){
    fit_mean <- face::pspline(data,argvals.new=tnew,knots=knots.initial,lambda=lambda_mean)
    mu.new <- fit_mean$mu.new
    r <- y - fit_mean$fitted.values 
  }
  #########################
  ####step 2:raw estimates
  #########################
  indW <- F # whether identity W
  if(is.null(W)) indW <- T
  
  raw <- raw.construct(data.frame("argvals" = t, "subj" = subj, "y" = as.vector(r)))
  C <- raw$C
  st <- raw$st
  N <- raw$st
  N2 <- raw$N2
  if(indW) W <- raw$W  
  n0 <- raw$n0
  
  delta <- Matrix((st[,1]==st[,2]) * 1) # sparse
  
  #########################
  ####step 3: smooth
  #########################
  knots <- construct.knots(t,knots,knots.option,p)
  
  List <- pspline.setting(st[,1],knots=knots,p,m,type="simple",knots.option=knots.option)
  B1 <- List$B
  B1 <- Matrix(B1)
  DtD <- List$P
  
  B2 = spline.des(knots=knots, x=st[,2], ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  c = dim(B1)[2]
  c2 = c*(c+1)/2
  B = Matrix(t(KhatriRao(Matrix(t(B2)),Matrix(t(B1)))))
  G = Matrix(duplication.matrix(c))
  
  BtWB = matrix(0,nrow=c^2,ncol=c^2)
  Wdelta = c()
  WC = c()
  for(i in 1:n0){
    seq = (sum(N2[1:i])-N2[i]+1):(sum(N2[1:i]))
    B3 = Matrix(matrix(B[seq,],nrow=length(seq)))
    W3 = W[[i]] # don't form a large W
    BtWB = BtWB + crossprod(B3, W3%*%B3)
    Wdelta <- c(Wdelta,as.matrix(W3 %*% delta[seq]))
    WC <- c(WC,as.matrix(W3 %*% C[seq]))
  }
  
  GtBtWBG = crossprod(G,BtWB%*%G)
  
  BG = B%*%G # sparse
  detWde <- crossprod(delta,Wdelta) # detWde = sum(delta)
  GtBtWdelta <- crossprod(BG,Wdelta)
  XtWX <- rbind(cbind(GtBtWBG,GtBtWdelta), cbind(t(GtBtWdelta),detWde))
  
  eSig = eigen(XtWX,symmetric=TRUE)
  V = eSig$vectors
  E = eSig$values
  E = E + 0.000001*max(E)
  Sigi_sqrt = matrix.multiply(V,1/sqrt(E))%*%t(V)
  
  P = crossprod(G,Matrix(suppressMessages(kronecker(diag(c),DtD))))%*%G
  
  Q = bdiag(P,0)
  tUQU = crossprod(Sigi_sqrt,(Q%*%Sigi_sqrt))
  Esig = eigen(tUQU,symmetric=TRUE)
  
  U = Esig$vectors
  s = Esig$values
  A0 <- Sigi_sqrt%*%U
  X <- cbind(BG,delta)
  A = as.matrix(X%*%A0) # F=XA dense
  
  AtA = crossprod(A) # diff
  f = crossprod(A,C) # diff
  ftilde = crossprod(A,WC) # diff
  
  c2 <- c2 + 1
  g <- rep(0, c2)
  G1 <- matrix(0,c2,c2)
  mat_list <- list()
  
  for(i in 1:n0){
    seq = (sum(N2[1:i])-N2[i]+1):(sum(N2[1:i]))
    Ai = matrix(A[seq,],nrow=length(seq))
    AitAi = crossprod(Ai) #t(Ai)%*%Ai
    Wi = W[[i]]
    
    fi = crossprod(Ai,C[seq]) # t(Fi)Ci
    Ji = crossprod(Ai,Wi%*%C[seq])
    Li = crossprod(Ai,Wi%*%Ai)
    g = g + Ji*fi
    G1 = G1 + AitAi*(Ji%*%t(ftilde))
    
    LList <- list()
    LList[[1]] = AitAi
    LList[[2]] = Li
    mat_list[[i]] = LList
    
  }
  
  
  Lambda <- seq(lower,upper,length=search.length)
  Gcv <- 0*Lambda
  gcv <- function(x){
    lambda <- exp(x)
    d <- 1/(1+lambda*s)
    ftilde_d <- ftilde*d
    cv0 <- -2*sum(ftilde_d*f)
    cv1 <-  sum(ftilde_d*(AtA%*%ftilde_d))
    cv2 <-  2*sum(d*g)
    cv3 <-  -4*sum(d*(G1%*%d))
    cv4 <- sum(unlist(sapply(mat_list,function(x){
      a <- x[[1]]%*%ftilde_d
      b <- x[[2]]%*%ftilde_d
      2*sum(a*b*d)
    })))
    cv <- cv0 + cv1 + cv2 + cv3 + cv4
    return(cv)
  }
  if(is.null(lambda)){
    Lambda <- seq(lower,upper,length=search.length)
    Length <- length(Lambda)
    Gcv <- rep(0,Length)
    for(i in 1:Length) 
      Gcv[i] <- gcv(Lambda[i])
    i0 <- which.min(Gcv)
    lambda <- exp(Lambda[i0])
  }
  
  alpha <- matrix.multiply(A0,1/(1+lambda*s))%*%ftilde
  Theta <- G %*% alpha[1:c2-1]
  Theta <- matrix(Theta,c,c)         # parameter estimated (sym)
  Theta0 <- Theta
  sigma2 <- alpha[c2]
  if(sigma2 <= 0.000001) {                                               
    warning("error variance cannot be non-positive, reset to 1e-6!")    
    sigma2 <- 0.000001                                                  
  }
  
  Eigen <- eigen(Theta,symmetric=TRUE)
  Eigen$values[Eigen$values<0] <- 0
  npc <- sum(Eigen$values>0) #which.max(cumsum(Eigen$values)/sum(Eigen$values)>pve)[1]
  if(npc >1){
    Theta <- matrix.multiply(Eigen$vectors[,1:npc],Eigen$values[1:npc])%*%t(Eigen$vectors[,1:npc])
    # Theta_half <- matrix.multiply(Eigen$vectors[,1:npc],sqrt(Eigen$values[1:npc]))
  }
  if(npc==1){
    Theta <- Eigen$values[1]*suppressMessages(kronecker(Eigen$vectors[,1],t(Eigen$vectors[,1])))
    # Theta_half <- sqrt(Eigen$values[1])*Eigen$vectors[,1]
  }
  Eigen <- eigen(Theta,symmetric=TRUE)
  
  #########################
  ####step 4: calculate estimated covariance function
  #########################
  Bnew = spline.des(knots=knots, x=tnew, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  Chat.new = as.matrix(tcrossprod(Bnew%*%Matrix(Theta),Bnew)) 
  Chat.diag.new = as.vector(diag(Chat.new))  
  Cor.new = diag(1/sqrt(Chat.diag.new))%*%Chat.new%*%diag(1/sqrt(Chat.diag.new))
  Eigen.new = eigen(Chat.new,symmetric=TRUE)
  npc = which.max(cumsum(Eigen$values)/sum(Eigen$values)>pve)[1] #determine number of PCs
  eigenfunctions = matrix(Eigen.new$vectors[,1:min(npc,length(tnew))],ncol=min(npc,length(tnew)))
  eigenvalues = Eigen.new$values[1:min(npc,length(tnew))]
  eigenfunctions = eigenfunctions*sqrt(length(tnew))/sqrt(max(tnew)-min(tnew))
  eigenvalues = eigenvalues/length(tnew)*(max(tnew)-min(tnew))
  
  
  #########################
  ####step 5: calculate variance
  #########################
  var.error.hat <- rep(sigma2,length(t))
  var.error.new <- rep(sigma2,length(tnew))
  
  
  
  Chat.raw.new = as.matrix(tcrossprod(Bnew%*%Matrix(Theta),Bnew)) + diag(var.error.new) 
  Chat.raw.diag.new = as.vector(diag(Chat.raw.new)) 
  Cor.raw.new = diag(1/sqrt(Chat.raw.diag.new))%*%Chat.raw.new%*%diag(1/sqrt(Chat.raw.diag.new))
  #########################
  ####step 6: prediction
  #########################
  
  
  if(!is.null(newdata)){
    
    mu.pred <- rep(0,length(newdata$argvals))
    var.error.pred <- rep(sigma2,length(newdata$argvals))
    if(center){
      mu.pred <- predict(fit_mean,newdata$argvals)
    }
    
    subj.pred = newdata$subj
    subj_unique.pred = unique(subj.pred)
    y.pred = newdata$y
    # Chat.diag.pred = 0*y.pred
    se.pred = 0*y.pred
    
    B = spline.des(knots=knots, x=newdata$argvals, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
    Chat.pred = as.matrix(tcrossprod(B%*%Matrix(Theta),B))
    Chat.diag.pred = diag(Chat.pred)
    
    scores = list(subj=subj_unique.pred,
                  scores = matrix(NA,nrow=length(subj_unique.pred),ncol=npc),
                  u = matrix(NA,nrow=length(subj_unique.pred),ncol=nrow(Theta))
    )
    
    Bi_pred <- vector("list", length(subj_unique.pred))
    Bi <- vector("list", length(subj_unique.pred))
    
    for(i in 1:length(subj_unique.pred)){
      sel.pred = which(subj.pred==subj_unique.pred[i])
      lengthi = length(sel.pred)
      
      pred.points <- newdata$argvals[sel.pred]
      mu.predi <- mu.pred[sel.pred]
      var.error.predi <- var.error.pred[sel.pred]
      
      y.predi = y.pred[sel.pred] - mu.predi
      sel.pred.obs = which(!is.na(y.predi))
      obs.points <- pred.points[sel.pred.obs]
      if(!is.null(obs.points)){
        var <- mean(var.error.predi[sel.pred.obs])
        if(var==0&length(sel.pred.obs) < npc)
          stop("Measurement error estimated to be zero and there are fewer observed points thans PCs; scores
               cannot be estimated.")
        B3i.pred = spline.des(knots=knots, x=pred.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
        B3i = spline.des(knots=knots, x=obs.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
        Bi_pred[[i]] <- B3i.pred
        Bi[[i]] <- B3i
        Chati = tcrossprod(B3i%*%Theta,B3i)
        # Chat.diag.pred[sel.pred] = diag(Chati)
        if(length(sel.pred.obs)==1) Ri = var.error.predi[sel.pred.obs]
        if(length(sel.pred.obs)>1) Ri = diag(var.error.predi[sel.pred.obs])
        Vi.inv = as.matrix(solve(Chati + Ri))
        Vi.pred = tcrossprod(B3i.pred%*%Theta,B3i.pred)
        Hi = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i)%*%Vi.inv)
        ui =tcrossprod(Theta,B3i)%*%Vi.inv %*%y.predi[sel.pred.obs]
        scores$u[i,] = as.vector(ui)
        y.pred[sel.pred] = as.numeric(Hi%*%y.predi[sel.pred.obs]) + mu.predi
        temp = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i))
        if(length(sel.pred.obs) >1){
          se.pred[sel.pred] = sqrt(diag(Vi.pred - temp%*%Vi.inv%*%t(temp)))
        }
        if(length(sel.pred.obs) ==1){
          se.pred[sel.pred] = sqrt(Vi.pred[1,1] - Vi.inv[1,1]*temp%*%t(temp))
        }
        
        ## predict scores
        if(calculate.scores==TRUE){ 
          temp = matrix(t(eigenfunctions),nrow=npc)%*%(as.matrix(Bnew)%*%ui)/sum(eigenfunctions[,1]^2)
          temp = as.matrix(temp)
          scores$scores[i,1:npc] = temp[,1]
        }
      }
    }
  }## if(is.null(newdata))
  if(is.null(newdata)){
    y.pred=NULL
    mu.pred = NULL
    var.error.pred = NULL
    Chat.diag.pred = NULL
    se.pred = NULL
    scores=NULL
    Bi <- NULL
    Bi_pred <- NULL
    B <- NULL
    
  }
  
  
  
  res <- list(newdata=newdata, W = W, y.pred = y.pred, Theta=Theta,argvals.new=tnew, 
              mu.new = mu.new, Chat.new=Chat.new, var.error.new = var.error.new,
              Cor.new = Cor.new, eigenfunctions = eigenfunctions, eigenvalues = eigenvalues,
              Cor.raw.new = Cor.raw.new, Chat.raw.diag.new = Chat.raw.diag.new,
              scores = scores, calculate.scores=calculate.scores,
              mu.hat = fit_mean$fitted.values,var.error.hat = var.error.hat,
              mu.pred = mu.pred, var.error.pred = var.error.pred, Chat.diag.pred = Chat.diag.pred,
              se.pred = se.pred,
              fit_mean = fit_mean, lambda_mean=fit_mean$lambda,
              lambda=lambda,Gcv=Gcv,Lambda=Lambda,knots=knots,knots.option=knots.option,s=s,npc=npc, p = p, m=m,
              center=center,pve=pve,sigma2=sigma2, r = r, DtD = DtD,
              Theta0 = Theta0, Bnew = Bnew, Bi = Bi, Bi_pred = Bi_pred, B = B)
  
  class(res) <- "face.sparse"
  return(res)
}


# 3
check.data <- function(data,type="fit"){
    if(type=="fit") {
      if(!is.data.frame(data)|is.null(data$y)|is.null(data$subj)|is.null(data$argvals))
        stop("'data' should be a data frame with three variables:argvals,subj and y")
      
      if(sum(is.na(data))>0) stop("No NA values are allowed in the data")
      
    }
    
    if(type=="predict"){
      if(!is.data.frame(data)|is.null(data$y)|is.null(data$subj)|is.null(data$argvals))
        stop("'newdata' should be a data frame with three variables:argvals,subj and y") 
    }
    return(0)
}



# 4
raw.construct <- function(data,include.diag=TRUE){
    
    y <- data$y
    t <- data$argvals
    subj <- data$subj
    
    subj_unique <- unique(subj)
    n <- length(subj_unique)
    C <- c()
    st <- matrix(NA,ncol=2,nrow=0)
    N <- c()
    N2 <- c()
    n0 <- 0
    W <- list(length=n)
    for(i in 1:n){
      
      r1 <- y[subj==subj_unique[i]]
      t1 <- t[subj==subj_unique[i]]
      m1 <- length(t1)
      n0 <- n0 + 1
      if(m1>1){
        if(include.diag) {
          N2 <-c(N2,m1*(m1+1)/2)  # <------
          sel = 1:N2[n0]
        }
        
        if(!include.diag) {
          N2 <-c(N2,m1*(m1-1)/2)  # <------
          sel = setdiff(1:(m1*(m1+1)/2), c(1,1 + cumsum(m1:1)[1:(m1-1)]))
        }
        
        st <- rbind(st,cbind(vech(kronecker(t1,t(rep(1,m1)))),
                             vech(kronecker(rep(1,m1),t(t1))))[sel,])
        C <- c(C,vech(kronecker(r1,t(r1)))[sel]) 
        
        
        N <- c(N,m1)
        # N2 <-c(N2,m1^2)
        
        W[[i]] <- sparseMatrix(1:N2[n0],1:N2[n0],x=rep(1,N2[n0]))# <----
        #if(include.diag) diag(W[[i]])[c(1,1 + cumsum(m1:1)[1:(m1-1)])] <- 1/2
      }## for if(m1>1)
      if(m1==1){
        if(include.diag){
          N2 <- c(N2,1)
          st <- rbind(st,c(t1,t1))
          C <- c(C,r1^2)
          N <- c(N,1)
          W[[i]] <- matrix(1,1,1)
        }
        if(!include.diag){
          N2 <- c(N2,0)
          N <- c(N,1)
          W[[i]] <- NULL
        }
      }
    }##for i
    
    res <- list("C" = C,
                "st" = st,
                "N" = N,
                "N2" = N2,
                "W" = W,
                "n0" = n0)
    return(res)
}


# 5
construct.knots <- function(argvals,knots,knots.option,p){
    
    if(length(knots)==1){
      allknots <- select.knots(argvals,knots,p=p,option=knots.option)
    }
    
    if(length(knots)>1){
      K = length(knots)-1 
      knots_left <- 2*knots[1]-knots[p:1+1]
      knots_right <- 2*knots[K] - knots[K-(1:p)]
      if(p>0) allknots <- c(knots_left,knots,knots_right)
      if(p==0) allknots <- knots
    }
    
    return(allknots)
    
}



# 6
pspline.setting <- function(x,knots=select.knots(x,35),
           p=3,m=2,weight=NULL,type="full",
           knots.option="equally-spaced"){
    
    # x: the marginal data points
    # knots: the list of interior knots or the numbers of interior knots
    # p: degrees for B-splines, with defaults values 3
    # m: orders of difference penalty, with default values 2
    # knots.option: type of knots placement, with default values "equally-spaced"
    
    #require(splines)
    #require(Matrix)
    
    ### design matrix 
    K = length(knots)-2*p-1
    B = spline.des(knots=knots, x=x, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
    
    bs = "ps"
    s.object = s(x=x, bs=bs, k=K+p,m=c(p-1,2), sp=NULL)
    
    if(knots.option == "quantile"){
      bs = "bs"
      s.object = s(x=x, bs=bs, k=K+p,m=c(p,2), sp=NULL)
    }
    
    
    object  = smooth.construct(s.object,data = data.frame(x=x),knots=list(x=knots))
    P =  object$S[[1]]
    if(knots.option == "quantile") P = P / max(abs(P))*10 # rescaling
    
    if(is.null(weight)) weight <- rep(1,length(x))
    
    if(type=="full"){
      
      Sig = crossprod(matrix.multiply(B,weight,option=2),B)
      eSig = eigen(Sig)
      V = eSig$vectors
      E = eSig$values
      if(min(E)<=0.0000001) {#cat("Warning! t(B)%*%B is singular!\n");
        #cat("A small identity matrix is added!\n");
        E <- E + 0.000001;
        
      }
      Sigi_sqrt = matrix.multiply(V,1/sqrt(E))%*%t(V)
      
      tUPU = Sigi_sqrt%*%(P%*%Sigi_sqrt)
      Esig = eigen(tUPU,symmetric=TRUE)
      U = Esig$vectors
      s = Esig$values
      s[(K+p-m+1):(K+p)]=0
      A = B%*%(Sigi_sqrt%*%U)
    }
    
    if(type=="simple"){
      
      A = NULL
      s = NULL
      Sigi_sqrt = NULL
      U = NULL
    }
    List = list(
      "A" = A,
      "B" = B,
      "s" = s,
      "Sigi.sqrt" = Sigi_sqrt,
      "U" = U,
      "P" = P)
    
    return(List)
}


# 7
matrix.multiply <- function(A,s,option=1){
    if(option==2)
      return(A*(s%*%t(rep(1,dim(A)[2]))))
    if(option==1)
      return(A*(rep(1,dim(A)[1])%*%t(s)))
}



# 8
kr <- function (A, B, w, byrow = TRUE) 
{
  if (byrow) {
    if (nrow(A) != nrow(B)) 
      stop("Dimensions of the matrices do not match.")
    if (missing(w)) 
      w <- rep(1, nrow(A))
    if (nrow(A) != length(w)) 
      stop("Length of the weight does not match with the dimension of the matrices.")
    cola <- ncol(A)
    colb <- ncol(B)
    colab <- cola * colb
    
    expr <- paste("rbind(", paste(rep("A", colb), collapse = ","), 
                  ")", sep = "")
    A <- eval(parse(text = expr))
    A <- matrix(c(A), nrow(B), ncol = colab)
    A <- w * A
    expr2 <- paste("cbind(", paste(rep("B", cola), collapse = ","), 
                   ")", sep = "")
    B <- eval(parse(text = expr2))
  }
  else {
    if (ncol(A) != ncol(B)) 
      stop("Dimensions of the matrices do not match.")
    if (missing(w)) 
      w <- rep(1, ncol(A))
    if (ncol(A) != length(w)) 
      stop("Length of the weight does not match with the dimension of the matrices.")
    rowa <- nrow(A)
    rowb <- nrow(B)
    rowab <- rowa * rowb
    A <- matrix(rep(A, each = rowb), rowab, )
    A <- A * matrix(w, nrow(A), ncol(A), byrow = TRUE)
    expr <- paste("rbind(", paste(rep("B", rowa), collapse = ","), 
                  ")", sep = "")
    B <- eval(parse(text = expr))
  }
  return(A * B)
}


# 9
igcv.criteria <- function(B,Y,mat_list,f,P1,P2,N2,Rho,w,c2){
  S <- w * P1 + (1 - w) * P2
  EigenS <- eigen(S)
  
  s <- EigenS$values
  s[s<0.00000001] <- 0
  U = EigenS$vectors
  tU = t(EigenS$vectors)
  
  ftilde = tU %*% f
  tftilde = t(ftilde)
  
  g <- rep(0, c2)
  G <- matrix(0,c2,c2)
  Li_list <- list()
  
  # Time <- proc.time()
  for(i in 1:length(N2)){
    
    fi = tU %*% mat_list[[i]][[1]]
    # fi = mat_list[[i]][[1]]
    Li = as.matrix(tU %*% mat_list[[i]][[2]]) %*% U
    # Li = mat_list[[i]][[2]]
    
    g = g + fi*fi
    G = G + Li*(fi%*%tftilde)
    
    Li_list[[i]] = Li
    
  }
  
  igcv <- function(x){
    d <- 1/(1+x*s)
    ftilde_d <- ftilde*d
    cv0 <- -2*sum(ftilde_d*ftilde)
    cv1 <-  sum(ftilde_d^2)
    cv2 <-  2*sum(d*g)
    cv3 <-  -4*sum(d*(G%*%d))
    cv4 <- sum(unlist(sapply(Li_list,function(x){
      a <- x%*%ftilde_d
      2*sum(a^2*d)
    })))
    cv <- cv0 + cv1 + cv2 + cv3 + cv4
    return(cv)
  }
  
  Length <- length(Rho)
  iGcv <- rep(0,Length)
  # Time <- proc.time()
  for(i in 1:Length) 
    iGcv[i] <- igcv(Rho[i])
  # print(proc.time()-Time)
  index <- which.min(iGcv)
  rho.min <- Rho[index]
  res <- list("rho"=rho.min, "igcv" = iGcv[index])
  return(res)
}



# 10
igcv.wrapper <- function(Y, X, P, N2, lower, upper, search.length, 
                         search.length0 = search.length){
  
  Rho <- exp(seq(lower,upper,length=search.length))
  List <- seq(0.01,0.99,length=search.length0) 
  fit_list <- list(length=length(List))
  # Time <- proc.time()
  BtB <- crossprod(X)
  EigenB <- eigen(BtB)
  EigenB$values <- sapply(EigenB$values, function(x) max(x,0.00000001))
  B_inv_half <- EigenB$vectors%*%diag(sqrt(1/EigenB$values))%*%t(EigenB$vectors)
  
  P1 <- B_inv_half%*%P[[1]]%*%B_inv_half
  P2 <- B_inv_half%*%P[[2]]%*%B_inv_half
  
  B <- as.matrix(X %*% B_inv_half)
  f = crossprod(B,Y) 
  c2 = ncol(BtB)
  
  mat_list <- list()
  for(i in 1:length(N2)){
    
    seq = (sum(N2[1:i])-N2[i]+1):(sum(N2[1:i]))
    Bi = matrix(B[seq,], nrow = length(seq))  
    fi = crossprod(Bi,Y[seq])
    Li = crossprod(Bi)
    
    LList <- list()
    LList[[1]] = fi
    LList[[2]] = Li
    mat_list[[i]] = LList
    
  }
  # print(proc.time()-Time)
  # Time <- proc.time()
  for(i in 1:length(List)){
    w = List[i]
    fit_list[[i]] <- igcv.criteria(B,Y,mat_list,f,P1,P2,N2,Rho,w,c2)
  }
  # print(proc.time()-Time)
  index <- which.min(sapply(fit_list,function(x) x$igcv))
  res <- fit_list[[index]]
  res$lambda <- c(List[index]*res$rho, (1-List[index])*res$rho)
  
  return(res)
}



# 11
bps <- function(y, x, z, N2, knots=10, p=3, m=2, lower=-3, upper=10, 
           search.length=14,
           knots.option="equally-spaced",lambda=NULL){
    
    
    data <- data.frame("x" = x, "z" = z)
    
    # bs <- rep(bs,2)[1:2]
    knots <- rep(knots,2)[1:2]
    p <- rep(p,2)[1:2]
    m <- rep(m,2)[1:2]
    knots.option <- rep(knots.option,2)[1:2]
    
    ## x-axis
    x.knots <- construct.knots(data$x,knots[1],knots.option[1],p[1])
    x.List <- pspline.setting(data$x,knots=x.knots,p[1],m[1],
                              type="simple",knots.option=knots.option[1])
    
    # x.object <- s(x,bs=bs[1],k=k[1])
    # x.knots <- NULL
    # if(knots.option[1] == "quantile") x.knots <- select.knots(x,knots = k[1] -3)
    
    # object <- smooth.construct(x.object, 
    # data = data.frame(x = data$x),knots=list(x=x.knots))
    
    Bx <- x.List$B
    Px <- x.List$P
    cx <- nrow(Px)
    
    x.object = list(x.List=x.List,
                    knots = x.knots,
                    knots.option = knots.option[1])
    
    
    ## z-axis
    # z.object <- s(z,bs=bs[2],k=k[2])
    # z.knots <- NULL
    # if(knots.option[2] == "quantile") z.knots <- select.knots(z,knots = k[2] -3)
    
    # object <- smooth.construct(z.object,
    # data = data.frame(z = data$z),knots=list(z=z.knots))
    
    z.knots <- construct.knots(data$z,knots[2],knots.option[2],p[2])
    z.List <- pspline.setting(data$z,knots=z.knots,p[2],m[2],
                              type="simple",knots.option=knots.option[2])
    
    Bz <- z.List$B
    Pz <- z.List$P
    cz <- nrow(Pz)
    
    z.object = list(z.List=z.List,
                    knots = z.knots,
                    knots.option = knots.option[2])
    
    B <- Matrix(kr(as.matrix(Bz),as.matrix(Bx)))
    P <- list()
    P[[1]] <- Matrix(kronecker(diag(cz),Px))
    P[[2]] <- Matrix(kronecker(Pz,diag(cx)))
    
    if(is.null(lambda)){
      fit <- igcv.wrapper(Y=y,X=B,P=P,N2=N2,
                          lower=lower,upper=upper,search.length=search.length)
      lambda <- fit$lambda
    }
    
    if(!is.null(lambda)) lambda <- rep(lambda,2)[1:2]
    
    Ptotal <- P[[1]]*lambda[1] + P[[2]]*lambda[2] 
    # print(lambda)
    temp <- tcrossprod(solve(crossprod(B) + Ptotal),B)
    theta <- temp%*%y
    Theta <- matrix(theta,nrow=ncol(Bx),ncol=ncol(Bz))
    
    
    res <- list(fitted.values = as.vector(B%*%theta), B = B,
                theta = theta,lambda = lambda,
                Theta=Theta,x.object = x.object,z.object=z.object, p = p)
    class(res) <- "bps"
    return(res)
}







