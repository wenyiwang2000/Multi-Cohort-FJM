FJM <- function(longdata, time, event, X0 = NULL, argvals.new, 
                  npc0 = NULL, npc1 = NULL, face_fit = NULL,
                  nMC0 = 100, nMC1 = 5000, knots = 9,
                  Maxiter = 1000, Miniter = 20, rho = 0.5, 
                  tol0 = 1e-3, tol1 = 5e-3, tol2 = 1e-3, tol3 = 1e-7, 
                  seed = 723, trace = F, calculate.scores = T){
  
  require(survival)
  require(MASS)
  require(Matrix)
  require(splines)
  
  #######################
  ####step 0: Data Input
  #######################
  # argvals.new = tnew; pve = 0.99; npc = length(Lambda0); # <---
  # nMC=1000; Maxiter = 1000; Miniter = 20; trace = T;
  # tol0 = 1e-3; tol1 = 5e-3; tol2 = 1e-3
  
  ####################################################
  ####step 1: Initialization thru Two-step Procedure
  ####################################################
  set.seed(seed)
  p <- length(longdata) 
  n <- length(time)
  nMC <- nMC0
  
  if(is.null(npc0) || is.null(npc1)){
    npc0 <- face_fit$npc0
    npc1 <- face_fit$npc1
  }
  
  G <- crossprod(face_fit$fit$y1$Bnew) / nrow(face_fit$fit$y1$Bnew)
  eig_G <- eigen(G, symmetric = T)
  G_half <- eig_G$vectors %*% diag(sqrt(eig_G$values)) %*% t(eig_G$vectors)
  face_fit$G_half <- G_half
  
  if((npc0 <= face_fit$npc0) & (npc1 <= face_fit$npc1)){
    scr_idx <- c(1:npc0, c(sapply(1:p, function(x){face_fit$npc0+(x-1)*face_fit$npc1+1:npc1})))
  }else if((npc0 > face_fit$npc0) & (npc1 <= face_fit$npc1)){
    if(npc0-face_fit$npc0==2){
      scr_idx <- c(1:face_fit$npc0,face_fit$npc0+1,face_fit$npc0+2, c(sapply(1:p, function(x){face_fit$npc0+(x-1)*face_fit$npc1+1:npc1})))
    }else{
      scr_idx <- c(1:face_fit$npc0,rep(face_fit$npc0+1,npc0-face_fit$npc0), c(sapply(1:p, function(x){face_fit$npc0+(x-1)*face_fit$npc1+1:npc1})))
    }
  }else if((npc0 <= face_fit$npc0) & (npc1 > face_fit$npc1)){
    scr_idx <- c(1:npc0, c(sapply(1:p, function(x){face_fit$npc0+(x-1)*face_fit$npc1+c(1:face_fit$npc1,rep(face_fit$npc1,npc1-face_fit$npc1))})))
  }else{
    if(npc0-face_fit$npc0==2){
      scr_idx <- c(1:face_fit$npc0,face_fit$npc0+1,face_fit$npc0+2, c(sapply(1:p, function(x){face_fit$npc0+(x-1)*face_fit$npc1+c(1:face_fit$npc1,rep(face_fit$npc1,npc1-face_fit$npc1))})))
    }else{
      scr_idx <- c(1:face_fit$npc0,rep(face_fit$npc0+1,npc0-face_fit$npc0), c(sapply(1:p, function(x){face_fit$npc0+(x-1)*face_fit$npc1+c(1:face_fit$npc1,rep(face_fit$npc1,npc1-face_fit$npc1))})))
    }
    #scr_idx <- c(1:face_fit$npc0,rep(face_fit$npc0+1,npc0-face_fit$npc0), c(sapply(1:p, function(x){face_fit$npc0+(x-1)*face_fit$npc1+c(1:face_fit$npc1,rep(face_fit$npc1,npc1-face_fit$npc1))})))
  }
  #scr_idx <- c(1:npc0, c(sapply(1:p, function(x){face_fit$npc0+(x-1)*face_fit$npc1+1:npc1})))
  scores_mu <- face_fit$scores$scores[,scr_idx]
  scores_Sig <- face_fit$scores$Sig
  scores_Sig <- mapply(function(scores_Sig) {scores_Sig[scr_idx, scr_idx]},
                       scores_Sig = scores_Sig, SIMPLIFY = F)
  
  D0_hat <- face_fit$eigenvalues0[sapply(1:npc0,function(x){min(x,face_fit$npc0)})] 
  D1_hat <- face_fit$eigenvalues1[sapply(1:npc1,function(x){min(x,face_fit$npc1)})]
  
  u0_hat <- as.matrix(as.matrix(face_fit$U0)[,sapply(1:npc0,function(x){min(x,face_fit$npc0)})])
  u1_hat <- as.matrix(as.matrix(face_fit$U1)[,sapply(1:npc1,function(x){min(x,face_fit$npc1)})])
  
  #D0_hat <- face_fit$eigenvalues0[1:npc0] 
  #D1_hat <- face_fit$eigenvalues1[1:npc1] 
  
  #u0_hat <- as.matrix(as.matrix(face_fit$U0)[,1:npc0])
  #u1_hat <- as.matrix(as.matrix(face_fit$U1)[,1:npc1])
  
  Bi_bar <- lapply(1:p, function(x){
    sapply(1:n, function(i){matrix((face_fit$fit[[x]])$fit_mean$B[which(longdata[[x]]$subj==i),], 
                                   ncol=length((face_fit$fit[[x]])$fit_mean$theta))%*%face_fit$G_invhalf})
  })
  alpha_hat <- sapply(1:p, function(x){as.matrix(face_fit$fit[[x]]$fit_mean$theta)})
  sigma2_hat <- sapply(1:p, function(x){face_fit$fit[[x]]$sigma2})
  b_hat <- sapply(1:p, function(x){face_fit$Beta[x]})
  sub_ind <- lapply(1:p, function(x){sapply(1:n, function(i) {sum(longdata[[x]]$subj == i)})})
  
  
  if(!is.null(X0)){
    q = ncol(X0)
    cox_fit <- coxph(formula=Surv(time,event)~cbind(X0, scores_mu[,1:npc0])) ################
  }else{
    q = 0
    cox_fit <- coxph(formula=Surv(time,event)~scores_mu[,1:npc0]) ################
  }
  
  h_hat <- coxph.detail(cox_fit)$hazard # caution
  nevent <- coxph.detail(cox_fit)$nevent
  beta_hat <- unname(cox_fit$coefficients)
  #beta_hat[is.na(beta_hat)] <- 2.65
  Tj <- coxph.detail(cox_fit)$time
  
  ######################################################
  ####step 2: Joint Estimation thru EM Algorithm (E-step)
  ######################################################
  ## main loop of EM
  scores_mu <- split(t(scores_mu), rep(1:nrow(scores_mu),
                                       each = ncol(scores_mu)))
  
  ## grid of event time
  #Tj <- sort(unique(time[event != 0]))
  K = length(Tj)
  Ti_ind <- sapply(1:n, function(i) {max(1, sum(Tj <= time[i]))})
  
  Theta <- list("D0" = D0_hat, "D1" = D1_hat, 
                "sigma2" = sigma2_hat, "alpha" = alpha_hat,
                "h" = h_hat, "beta" = beta_hat, 
                "u0" = u0_hat, "u1" = u1_hat, "b" = b_hat)
  
  d1 <- c()
  d2 <- c()
  D0_d <- c()
  D1_d <- c()
  h_d <- c()
  beta_d <- c()
  sigma2_d <- c()
  alpha_d <- c()
  u0_d <- c()
  u1_d <- c()
  b_d <- c()
  ll <- c()
  ll_rel <- c()
  
  for(iter in 1:Maxiter){
    print(iter)
    
    if(iter >= Miniter) nMC = nMC1
    if(npc0==1){
      D0_hat <- as.matrix(diag(as.matrix(Theta$D0)))
    }else{
      D0_hat <- diag(Theta$D0)
    }
    if(npc1==1){
      D1_hat <- as.matrix(diag(as.matrix(Theta$D1)))
    }else{
      D1_hat <- diag(Theta$D1)
    }
    h_hat <- Theta$h
    beta_hat <- Theta$beta
    sigma2_hat <- Theta$sigma2
    alpha_hat <- Theta$alpha
    u0_hat <- Theta$u0
    u1_hat <- Theta$u1
    b_hat <- Theta$b
    
    phi0_hat = lapply(1:p, function(x){
      mapply(function(Bi_bar0){Bi_bar0%*%u0_hat}, Bi_bar0=Bi_bar[[x]], SIMPLIFY = F)}) 
    phi1_hat = lapply(1:p, function(x){
      mapply(function(Bi_bar0){Bi_bar0%*%u1_hat}, Bi_bar0=Bi_bar[[x]], SIMPLIFY = F)}) 
    
    
    ## update multivariate normal distribution
    if(iter > 1){ # should be > 1 (was >=1)
      Lam0 = rbind(D0_hat, matrix(0, p*npc1, npc0))
      Lam1_bd = kronecker(diag(p), D1_hat)
      Lam1 = rbind(matrix(0, npc0, p*npc1), Lam1_bd)
      Lam = bdiag(D0_hat, Lam1_bd)
      
      cor_hat <- lapply(1:n, function(i) 
      {Reduce(cbind, sapply(1:p, function(x){tcrossprod(b_hat[x]*Lam0, phi0_hat[[x]][[i]])}, simplify = F)) + 
          Reduce(cbind, sapply(1:p, function(x){tcrossprod(b_hat[x]*Lam1[, 1:npc1+(x-1)*npc1], 
                                                           phi1_hat[[x]][[i]])}, simplify = F))})
      
      phi0mat <- lapply(1:n, function(i) {
        as.matrix(bdiag(lapply(1:p, function(x){phi0_hat[[x]][[i]]})))})
      phi1mat <- lapply(1:n, function(i) {
        as.matrix(bdiag(lapply(1:p, function(x){phi1_hat[[x]][[i]]})))})
      D0mat <- kronecker(tcrossprod(b_hat), D0_hat)
      D1mat <- kronecker(diag(b_hat^2), D1_hat)
      
      cov_hat <- lapply(1:n, function(i) 
      {phi0mat[[i]]%*%D0mat%*%t(phi0mat[[i]]) + 
          phi1mat[[i]]%*%D1mat%*%t(phi1mat[[i]])})
      
      R_hat <- lapply(1:n, function(i){
        as.matrix(bdiag(lapply(1:p, function(x){
          if(sub_ind[[x]][i]==1){
            sigma2_hat[x]
          }else{
            sigma2_hat[x]*diag(sub_ind[[x]][i])}
        })))  })
      
      Vinv <- mapply(function(cov_hat, R_hat) {solve(cov_hat + R_hat)}, # Vinv = inv(cov(yi))
                     cov_hat = cov_hat, R_hat = R_hat, SIMPLIFY = F)
      scores_mu <- mapply(function(cor_hat, Vinv, W){ # scores_mu = E(eta_i|y_i)
        mu0 = cor_hat %*% Vinv %*% W
        if((mu0[1]*beta_hat[q+1])>8) mu0[1] = 8/(beta_hat[q+1])
        return(mu0)}, cor_hat = cor_hat, Vinv = Vinv, W = W, SIMPLIFY = F)
      scores_Sig <- mapply(function(cor_hat, Vinv) {  # scores_Sig = cov(eta_i|y_i)
        Lam - tcrossprod(cor_hat %*% Vinv, cor_hat)}, 
        cor_hat = cor_hat, Vinv = Vinv, SIMPLIFY = F)
    }
    ## MC samples of scores
    xi = mapply(function(mu,Sig){mvrnorm(nMC, matrix(mu), Sig)}, 
                mu = scores_mu, Sig = scores_Sig, SIMPLIFY = F)
    
    if(!is.null(X0)){ # expX = exp(z^T_i*gamma_z + eta^T_i*gamma_eta)
      expX0 = exp(X0 %*% beta_hat[1:q])
      expX = lapply(xi, function(xi){exp(as.matrix(xi[,1:npc0]) %*% beta_hat[-(1:q)])})  ################
      expX = mapply(function(expX0, expX){expX0*expX}, 
                    expX0 = expX0, expX = expX, SIMPLIFY = F)
    }else{
      expX = lapply(xi, function(xi){exp(xi[,1:npc0] %*% beta_hat)}) ################
    }
    
    H = sapply(1:n, function(i) {sum(h_hat[1:Ti_ind[i]])}, 
               simplify = F)
    hti = as.list(h_hat[Ti_ind])
    
    ## f(T_i, delta_i|.)
    fti = mapply(function(expX, H, hti, event){
      S = exp(- (H * expX))
      if(event==1){
        fti0 = hti * expX * S
      }else{
        fti0 = S
      }
      if(sum(fti0)==0){
        fti0 = matrix(tol0, nr=nMC, nc=1)
      }
      fti0
    }, expX = expX, H = H, hti = hti, event = event, SIMPLIFY = F)
    
    ## \sum{f(T_i, delta_i|.)}/nMC : den is the denominator in E_{g(eta_i)}
    den <- lapply(fti, mean)
    
    ## f(T_i, delta_i|.) / den
    pxi <- mapply(function(f, d) {c(f / d)}, f = fti, d = den, SIMPLIFY = F) ################
    
    ## E[xi]: Exi = E_i{g(eta_i)}
    Exi <- mapply(function(xi, pxi) {colMeans(xi*pxi)},
                  xi = xi, pxi = pxi, SIMPLIFY = F)
    
    ## E[X]: EX = E_i{(z_i,eta_i)}
    if(!is.null(X0)){
      EX <- lapply(1:n, function(i){c(X0[i, ], Exi[[i]][1:npc0])}) ####################
    }else{
      EX <- lapply(1:n, function(i){Exi[[i]][1:npc0]}) ####################
    }
    
    ## E[xi^2]: Exi2 = E_i{eta_i*eta^T_i}  #################### may change it because of zeta_ij
    Exi2 <- mapply(function(xi, pxi) {crossprod(xi, xi*pxi)/ nMC}, 
                   xi = xi, pxi = pxi, SIMPLIFY = F)
    
    ## E[expX]: E_i{exp(z^T_i*gamma_z + eta^T_i*gamma_eta)}
    EexpX <- mapply(function(pxi, expX) {mean(expX*pxi)},
                    pxi = pxi, expX = expX, SIMPLIFY = F)
    
    ## E[xi expX]: E_i{(eta_i) * exp(z^T_i*gamma_z + eta^T_i*gamma_eta)}
    ExiexpX <- mapply(function(xi, pxi, expX) {colMeans(as.matrix(xi[,1:npc0])*c(expX)*pxi)}, ###################
                      xi = xi, pxi = pxi, expX = expX, SIMPLIFY = F)
    
    ## E[X expX]: E_i{(z_i, eta_i) * exp(z^T_i*gamma_z + eta^T_i*gamma_eta)}
    if(!is.null(X0)){
      EXexpX <- lapply(1:n, function(i){c(X0[i, ] * EexpX[[i]], 
                                          ExiexpX[[i]])})
    }else{
      EXexpX <- ExiexpX
    }
    
    ######################################################
    ####step 3: Joint Estimation thru EM Algorithm (M-step)
    ######################################################
    ## D
    D_hat <- diag(diag(Reduce("+", Exi2)) / n)
    D0_hat <- D_hat[1:npc0, 1:npc0]
    D1_hat <- Reduce("+", sapply(1:p, function(x){
      D_hat[npc0+(x-1)*npc1+1:npc1, npc0+(x-1)*npc1+1:npc1]}, 
      simplify = F)) / p
    
    ##1.estimate alpha
    PExi <- lapply(1:p, function(x) {
      Exi0 = do.call("rbind",Exi)[,c(1:npc0, npc0+(x-1)*npc1+1:npc1)]
      unlist(lapply(1:n, function(i){
        PExi0 = phi0_hat[[x]][[i]] %*% Exi0[i,1:npc0] 
        PExi1 = phi1_hat[[x]][[i]] %*% Exi0[i,npc0+1:npc1]
        PExi0 + PExi1 }) ) })
    
    alpha_hat <- as.matrix( do.call("cbind", lapply(1:p, function(x){
      data_temp = data.frame("subj"=longdata[[x]][,1], "argvals"=longdata[[x]][,2],"y"=longdata[[x]][,3]-b_hat[x]*PExi[[x]])
      fit.mean <- face::pspline(data=data_temp, argvals.new=tnew, knots=knots)
      fit.mean$theta }) )) 
 
    ##2.estimate sigma2
    W_ind = sapply(1:p, function(x){
      as.matrix(cbind(longdata[[x]][,1],
                      longdata[[x]][,3] - face_fit$fit[[x]]$fit_mean$B %*% alpha_hat[,x])) }, simplify = F)
    W_p <- lapply(1:n, function(i) {
      sapply(1:p, function(x) {W_ind[[x]][W_ind[[x]][,1] == i, -1]}, simplify = F)})
    W <- lapply(W_p, function(W_p) {Reduce(c, W_p)})
    
    WtW = unname(Reduce(rbind, lapply(W_p, function(W_p) {
      sapply(1:p, function(x){crossprod(W_p[[x]])})})))
    
    WtP0 = lapply(1:n, function(j) {t(Reduce(rbind, sapply(1:p, function(x){
      crossprod(W_p[[j]][[x]], phi0_hat[[x]][[j]])}, simplify = F)))} )
    WtP1 = lapply(1:n, function(j) {t(Reduce(rbind, sapply(1:p, function(x){
      crossprod(W_p[[j]][[x]], phi1_hat[[x]][[j]])}, simplify = F)))} )
    
    P0tP0 = lapply(1:p, function(x){
      lapply(1:n, function(j) {crossprod(phi0_hat[[x]][[j]])}) })
    P1tP1 = lapply(1:p, function(x){
      lapply(1:n, function(j) {crossprod(phi1_hat[[x]][[j]])}) })
    P0tP1 = lapply(1:p, function(x) {
      mapply(function(phi0_hat, phi1_hat) {crossprod(phi0_hat, phi1_hat)}, 
             phi0_hat = phi0_hat[[x]], phi1_hat = phi1_hat[[x]], SIMPLIFY = F) })
    
    ## sigma_k^2
    b_xxt =  lapply(1:n, function(j){
      sapply(1:p, function(x) {
        sum(diag(P0tP0[[x]][[j]] %*% Exi2[[j]][1:npc0,1:npc0])) +
          sum(diag(P1tP1[[x]][[j]] %*% Exi2[[j]][npc0+npc1*(x-1)+1:npc1,npc0+npc1*(x-1)+1:npc1])) +
          2*sum(diag(P0tP1[[x]][[j]] %*% Exi2[[j]][npc0+npc1*(x-1)+1:npc1,1:npc0]))
      }) })
    
    temp = t(mapply(function(Exi,WtP0,WtP1,b_xxt)
    {sapply(1:p, function(x) {
      -2*b_hat[x]*WtP0[,x]%*%Exi[1:npc0] - 
        2*b_hat[x]*WtP1[,x]%*%Exi[npc0+npc1*(x-1)+1:npc1] +
        b_hat[x]^2*b_xxt[x]})},
    Exi = Exi, WtP0 = WtP0, WtP1 = WtP1, b_xxt = b_xxt))
    
    sigma2_hat <- (colSums(WtW) + colSums(temp)) / unlist(lapply(1:p,function(x){nrow(longdata[[x]])}))
    
    ##3.estimate b
    b_hat[1] <- 1
    b_x = lapply(1:n, function(j) {
      sapply(2:p, function(x) {
        phi0_hat[[x]][[j]] %*% Exi[[j]][1:npc0] + phi1_hat[[x]][[j]] %*% Exi[[j]][npc0+npc1*(x-1)+1:npc1]
      }, simplify = F)})
    b_res = mapply(function(b_x, W) {sapply(2:p, function(x) {
      crossprod(b_x[[x-1]], W[[x]])})}, b_x = b_x, W = W_p)
    
    if(p > 2){
      b_hat[-1] <- rowSums(b_res) / rowSums(Reduce(cbind,b_xxt))[-1]
    }else{
      b_hat[-1] <- sum(b_res) / rowSums(Reduce(cbind,b_xxt))[-1]
    }
    
    ##4.estimate u0, phi
    P1Exi2 = lapply(1:p, function(j){ # Psi*E(xi*zeta)
      do.call("rbind", mapply(function(phi1_temp, Exi20){
        phi1_temp%*%t(matrix(Exi20[1:npc0, npc0+npc1*(j-1)+1:npc1],nrow=npc0))}, 
        phi1_temp=phi1_hat[[j]], Exi20=Exi2) ) }) 
    Exi0_long <- lapply(1:p, function(j){
      Exi_m <- do.call("rbind",Exi)
      do.call("cbind",lapply(1:npc0, function(xx){rep(Exi_m[,xx], times=sub_ind[[j]])}))  }) 
    Exi20_long <- lapply(1:p, function(j){
      Exi2_m <- do.call("rbind",mapply(function(Exi20){diag(Exi20)}, Exi20=Exi2, SIMPLIFY=F))
      do.call("cbind",lapply(1:npc0, function(xx){rep(Exi2_m[,xx], times=sub_ind[[j]])}))  })  
    Exi1_long <- lapply(1:p,function(j){
      Exi_m <- do.call("rbind",Exi)
      do.call("cbind",lapply(1:npc1, function(xx){rep(Exi_m[,npc0+npc1*(j-1)+xx], times=sub_ind[[j]])}))
    }) 
    Exi21_long <- lapply(1:p,function(j){
      Exi2_m <- do.call("rbind",mapply(function(Exi20){diag(Exi20)}, Exi20=Exi2, SIMPLIFY=F))
      do.call("cbind",lapply(1:npc1, function(xx){rep(Exi2_m[,npc0+npc1*(j-1)+xx], times=sub_ind[[j]])}))
    })
    
    for(jj in 1:10){
      u0_pre = u0_hat
      for(k in 1:npc0){
        if(npc0 == 1){
          data_temp <- do.call("rbind", lapply(1:p, function(j){
            b_sqExi2_P1Exi2 <- (b_hat[j]/sqrt(Exi20_long[[j]][,k]))*P1Exi2[[j]][,k] #(beta_j/sq(Exi2))*Psi*Exi2
            Exi_sqExi2_Wij <- (Exi0_long[[j]][,k]/sqrt(Exi20_long[[j]][,k]))*W_ind[[j]][,2] #(Exi/sq(Exi2))*(y_ijk-u_ijk)
            y_temp <- Exi_sqExi2_Wij - b_sqExi2_P1Exi2 
            bsqExi2 <- b_hat[j] * sqrt(Exi20_long[[j]][,k]) #bj*sqrt(Exi2)
            argvals <- longdata[[j]]$argvals
            cbind(argvals, y_temp, bsqExi2)  }))
        } else {
          data_temp <- do.call("rbind", lapply(1:p, function(j){
            phi0_temp <- mapply(function(Bi_bar0){Bi_bar0%*%u0_hat}, Bi_bar0=Bi_bar[[j]], SIMPLIFY=F)
            P0Exi2 <- do.call("rbind", mapply(function(phi0_temp, Exi20){t(t(phi0_temp)*Exi20[1:npc0,1:npc0][k,])}, 
                                              phi0_temp=phi0_temp, Exi20=Exi2)) #Phi*Exi2
            b_sqExi2_P0Exi2 <- (b_hat[j]/sqrt(Exi20_long[[j]][,k]))*(rowSums(P0Exi2)-P0Exi2[,k]) #(beta_j/sq(Exi2))*Phi*Exi2
            b_sqExi2_P1Exi2 <- (b_hat[j]/sqrt(Exi20_long[[j]][,k]))*P1Exi2[[j]][,k] #(beta_j/sq(Exi2))*Psi*Exi2
            Exi_sqExi2_Wij <- (Exi0_long[[j]][,k]/sqrt(Exi20_long[[j]][,k]))*W_ind[[j]][,2] #(Exi/sq(Exi2))*(y_ijk-u_ijk)
            y_temp <- Exi_sqExi2_Wij - b_sqExi2_P1Exi2 - b_sqExi2_P0Exi2
            bsqExi2 <- b_hat[j] * sqrt(Exi20_long[[j]][,k]) #bj*sqrt(Exi2)
            argvals <- longdata[[j]]$argvals
            cbind(argvals, y_temp, bsqExi2)  }))
        }
        data_temp <- data.frame(data_temp)
        fit <- gam(y_temp ~ s(argvals, by=bsqExi2, bs='ps', k=knots+3, m=c(2,2)), drop.intercept=T, data=data_temp)
        u0_hat[,k] <- face_fit$G_half %*% fit$coefficients
      }
      if(max(abs(u0_pre - u0_hat)) < 1e-6) break
    }
    
    temp <- tcrossprod(u0_hat %*% D0_hat, u0_hat)
    temp <- 0.5*(temp+t(temp))
    eig0 <- eigen(temp)
    #eig0 <- eigen(tcrossprod(u0_hat %*% D0_hat, u0_hat))
    if(npc0 == 1){
      #D0_hat <- diag(as.matrix(eig0$values[1:npc0]))
      D0_hat <- diag(as.matrix(max(eig0$values[1:npc0],(sum(D0_hat)+sum(D1_hat))*0.0005)))
    }else{
      #D0_hat <- diag(eig0$values[1:npc0]) # D0_hat = D0
      D0_hat <- diag(sapply(eig0$values[1:npc0],function(x){max(x,(sum(D0_hat)+sum(D1_hat))*0.0005)}))
    }
    
    for(k in 1:npc0){
      if(crossprod(u0_hat[,k], eig0$vectors[,k]) < 0){
        u0_hat[,k] <- - eig0$vectors[,k]
      }else{
        u0_hat[,k] <- eig0$vectors[,k] 
      }
    }
    
    ##5.estimate u1, psi
    phi0_hat = lapply(1:p, function(j){
      mapply(function(Bi_bar0){Bi_bar0%*%u0_hat}, Bi_bar0=Bi_bar[[j]], SIMPLIFY = F)}) 
    P0Exi2 = lapply(1:p, function(j){ # Phi*E(xi*zeta)
        do.call("rbind", mapply(function(phi0_temp, Exi20){
          phi0_temp%*%t(matrix(Exi20[npc0+npc1*(j-1)+1:npc1,1:npc0],ncol=npc0))}, 
          phi0_temp=phi0_hat[[j]], Exi20=Exi2) ) })  
    
    for(jj in 1:10){
      u1_pre = u1_hat
      for(k in 1:npc1){
        if(npc1 == 1){
          data_temp <- do.call("rbind", lapply(1:p, function(j){
              b_sqExi2_P0Exi2 <- (b_hat[j]/sqrt(Exi21_long[[j]]))*P0Exi2[[j]] #(beta_j/sq(Exi2))*Phi*Exi2
              Exi_sqExi2_Wij <- (Exi1_long[[j]][,k]/sqrt(Exi21_long[[j]]))*W_ind[[j]][,2] #(Exi/sq(Exi2))*(y_ijk-u_ijk)
              y_temp <- Exi_sqExi2_Wij - b_sqExi2_P0Exi2 
              bsqExi2 <- b_hat[j] * sqrt(Exi21_long[[j]]) #bj*sqrt(Exi2)
              argvals <- longdata[[j]][,2]
              cbind(argvals, y_temp, bsqExi2)  })) 
        }else{
          data_temp <- do.call("rbind", lapply(1:p, function(j){
              phi1_temp <- mapply(function(Bi_bar0){Bi_bar0%*%u1_hat}, Bi_bar0=Bi_bar[[j]], SIMPLIFY=F)
              P1Exi2 <- do.call("rbind", mapply(function(phi1_temp, Exi20){t(t(phi1_temp)*Exi20[npc0+npc1*(j-1)+1:npc1,npc0+npc1*(j-1)+1:npc1][k,])}, 
                                                phi1_temp=phi1_temp, Exi20=Exi2)) #Psi*Exi2
              b_sqExi2_P1Exi2 <- (b_hat[j]/sqrt(Exi21_long[[j]][,k]))*(rowSums(P1Exi2)-P1Exi2[,k]) #(beta_j/sq(Exi2))*Psi*Exi2
              b_sqExi2_P0Exi2 <- (b_hat[j]/sqrt(Exi21_long[[j]][,k]))*P0Exi2[[j]][,k] #(beta_j/sq(Exi2))*Phi*Exi2
              Exi_sqExi2_Wij <- (Exi1_long[[j]][,k]/sqrt(Exi21_long[[j]][,k]))*W_ind[[j]][,2] #(Exi/sq(Exi2))*(y_ijk-u_ijk)
              y_temp <- Exi_sqExi2_Wij - b_sqExi2_P0Exi2 - b_sqExi2_P1Exi2
              bsqExi2 <- b_hat[j] * sqrt(Exi21_long[[j]][,k]) #bj*sqrt(Exi2)
              argvals <- longdata[[j]][,2]
              cbind(argvals, y_temp, bsqExi2)  }))
        }
        data_temp <- data.frame(data_temp)
        colnames(data_temp) <- c("argvals","y_temp","bsqExi2")
        fit <- gam(y_temp ~ s(argvals, by=bsqExi2, bs='ps', k=knots+3, m=c(2,2)), drop.intercept=T, data=data_temp)
        u1_hat[,k] <- face_fit$G_half %*% fit$coefficients
      }
      if(max(abs(u1_pre - u1_hat)) < 1e-6) break
    }
    
    temp <- tcrossprod(u1_hat %*% D1_hat, u1_hat)
    temp <- 0.5*(temp+t(temp))
    eig1 <- eigen(temp)
    #eig1 <- eigen(tcrossprod(u1_hat %*% D1_hat, u1_hat))
    if(npc1 == 1){
      D1_hat <- diag(as.matrix(max(eig1$values[1:npc1],(sum(D0_hat)+sum(D1_hat))*0.0005)))
    }else{
      D1_hat <- diag(sapply(eig1$values[1:npc1],function(x){max(x,(sum(D0_hat)+sum(D1_hat))*0.0005)})) # D1_hat = D1
    }
    
    for(k in 1:npc1){
      if(crossprod(u1_hat[,k], eig1$vectors[,k]) < 0){
        u1_hat[,k] <- - eig1$vectors[,k]
      }else{
        u1_hat[,k] <- eig1$vectors[,k] 
      }
    }
    
    ##6.estimate h0
    h_de <- sapply(Tj, function(x){sum(unlist(EexpX)*(time>=x))}) 
    h_hat <- nevent / h_de
    
    
    ## score of beta
    Si <- mapply(function(event, EX, H, EXexpX) {
      event*EX -  H*EXexpX}, event = event, EX = EX,
      H = H, EXexpX = EXexpX, SIMPLIFY = F)
    S <- Reduce("+", Si)
    
    ## information matrix
    I <- Reduce("+", lapply(Si, tcrossprod)) - tcrossprod(S) / n
    
    ##7.estimate beta
    gBeta = rho * solve(I, S,tol=1e-80)
    beta_hat <- beta_hat + gBeta
    
    ## collect updates 
    if(npc0==1){
      diag_D0 = diag(as.matrix(D0_hat))
    }else{
      diag_D0 = diag(D0_hat)
    }
    
    if(npc1==1){
      diag_D1 = diag(as.matrix(D1_hat))
    }else{
      diag_D1 = diag(D1_hat)
    }
    
    Theta_new <- list("D0" = diag_D0, "D1" = diag_D1, 
                      "sigma2" = sigma2_hat, 
                      "alpha" = alpha_hat, 
                      "h" = h_hat, "beta" = beta_hat, 
                      "u0" = u0_hat, "u1" = u1_hat, 
                      "b" = b_hat)
    
    # print(iter) # <------------------------
    
    # marginal (observed data) likelihood 
    if(iter > 1){
      loglik <- mapply(function(sub_ind, W, cov_hat, R_hat, Vinv) {
        - 0.5 * (length(W) * log(2*pi) + 
                   as.numeric(determinant(cov_hat+R_hat, logarithm = TRUE)$modulus) +
                   (t(W) %*% Vinv %*% W)[1,1])
      },
      W = W, cov_hat = cov_hat, R_hat = R_hat, Vinv = Vinv)
      loglik <-  sum(loglik + log(unlist(den)))
      ll <- c(ll, loglik)
      ll_rel <- c(ll_rel, abs(loglik-ll[length(ll)-1]) / 
                    (abs(ll[length(ll)-1]) + tol0))
    }
    
    
    ## convergence control (maximum absoulte relative change)
    d_rel <- mapply(function(Theta, Theta_new) {
      abs(Theta - Theta_new) / (abs(Theta) + tol0)}, 
      Theta = Theta, Theta_new = Theta_new)
    d_abs <- mapply(function(Theta, Theta_new) {
      abs(Theta - Theta_new)}, 
      Theta = Theta, Theta_new = Theta_new)
    max_rel <- max(unlist(d_rel))
    max_abs <- max(unlist(d_abs))
    if(iter > 2){
      flag <- (max_rel < tol1) || (max_abs < tol2) || (ll_rel[iter-2] < tol3)
    }else{
      flag <- (max_rel < tol1) || (max_abs < tol2) 
    }
    
    d1 <- c(d1, max_rel)
    d2 <- c(d2, max_abs)
    
    D0_d <- c(D0_d, c(max(d_rel$D0), max(d_abs$D0)))
    D1_d <- c(D1_d, c(max(d_rel$D1), max(d_abs$D1)))
    sigma2_d <- c(sigma2_d, c(max(d_rel$sigma2), max(d_abs$sigma2)))
    alpha_d <- c(alpha_d, c(max(d_rel$alpha), max(d_abs$alpha)))
    h_d <- c(h_d, c(max(d_rel$h), max(d_abs$h)))
    beta_d <- c(beta_d, c(max(d_rel$beta), max(d_abs$beta)))
    u0_d <- c(u0_d, c(max(d_rel$u0), max(d_abs$u0)))
    u1_d <- c(u1_d, c(max(d_rel$u1), max(d_abs$u1)))
    b_d <- c(b_d, c(max(d_rel$b), max(d_abs$b)))
    
    
    if(iter > Miniter && flag){
      Theta <- Theta_new
      break
    }else{
      Theta <- Theta_new
    }
  }  # end iter
  
  D0_d <- matrix(D0_d, nrow=2)
  D1_d <- matrix(D1_d, nrow=2)
  sigma2_d <- matrix(sigma2_d, nrow=2)
  alpha_d <- matrix(alpha_d, nrow=2)
  h_d <- matrix(h_d, nrow=2)
  beta_d <- matrix(beta_d, nrow=2)
  u0_d <- matrix(u0_d, nrow=2)
  u1_d <- matrix(u1_d, nrow=2)
  b_d <- matrix(b_d, nrow=2)
  
  d_details <- list("D0" = D0_d, "D1" = D1_d, 
                    "sigma2" = sigma2_d,
                    "alpha" = alpha_d, 
                    "h" = h_d, "beta" = beta_d, 
                    "u0" = u0_d, "u1" = u1_d, "b" = b_d)
  
  ## degree of freedom
  dof0 <- prod(dim(Theta$alpha)) + (npc0 + npc1) * (nrow(Theta$u0) + 1) + 
    npc0 + 2 * p - 1
  if(!is.null(X0)){
    dof0 <- dof0 + q
  }
  
  AIC0 <- -2 * loglik + 2 * dof0
  BIC0 <- -2 * loglik + log(n) * dof0
  
  dof <- dof0 - npc0 * (npc0 + 1) / 2 - npc1 * (npc1 + 1) / 2 
  AIC <- -2 * loglik + 2 * dof
  BIC <- -2 * loglik + log(n) * dof
  
  ## information matrix of D0 D1
  if(npc0 == 1){
    Si_D0 = -1/2 * 1 / Theta$D0 + 1/2 * 1 / (Theta$D0)^2 * 
      sapply(1:n, function(i){Exi2[[i]][1:npc0, 1:npc0]})
    S_D0 = sum(Si_D0)
    I_D0 <- crossprod(Si_D0) - S_D0^2 / n
  }else{
    Si_D0 = -1/2 * 1 / Theta$D0 + 1/2 * 1 / (Theta$D0)^2 * 
      sapply(1:n, function(i){diag(Exi2[[i]][1:npc0, 1:npc0])})
    S_D0 = rowSums(Si_D0)
    I_D0 <- tcrossprod(Si_D0) - tcrossprod(S_D0) / n
  }
  
  if(npc1 == 1){
    Si_D1 = -1/2 * p / Theta$D1 + 1/2 * 1 / (Theta$D1)^2 * 
      sapply(1:n, function(i){
        Reduce("+", sapply(1:p, function(x){
          Exi2[[i]][npc0+(x-1)*npc1+1:npc1, npc0+(x-1)*npc1+1:npc1]}, 
          simplify = F))})
    S_D1 = sum(Si_D1)
    I_D1 <- crossprod(Si_D1) - S_D1^2 / n
  }else{
    Si_D1 = -1/2 * p / Theta$D1 + 1/2 * 1 / (Theta$D1)^2 * 
      sapply(1:n, function(i){
        Reduce("+", sapply(1:p, function(x){
          diag(Exi2[[i]][npc0+(x-1)*npc1+1:npc1, npc0+(x-1)*npc1+1:npc1])}, 
          simplify = F))})
    S_D1 = rowSums(Si_D1)
    I_D1 <- tcrossprod(Si_D1) - tcrossprod(S_D1) / n
  }
  
  
  if(trace==T){
    par(mfrow=c(1,2),mar=c(5,4.8,3,1))
    
    plot(d1, main = bquote(Theta), ylab = "relative change", 
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol1, col="red", lwd=2)
    plot(d2, main = bquote(Theta), ylab = "absolute change", 
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol2, col="red", lwd=2)
    
    plot(d_details$D0[1,], main = bquote(lambda[0]), ylab = "relative change", 
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol1, col="red", lwd=2)
    plot(d_details$D0[2,], main = bquote(lambda[0]), ylab = "absolute change", 
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol2, col="red", lwd=2)
    
    plot(d_details$D1[1,], main = bquote(lambda[1]), ylab = "relative change", 
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol1, col="red", lwd=2)
    plot(d_details$D1[2,], main = bquote(lambda[1]), ylab = "absolute change", 
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol2, col="red", lwd=2)
    
    plot(d_details$sigma2[1,], main = bquote(sigma^2),
         ylab = "relative change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol1, col="red", lwd=2)
    plot(d_details$sigma2[2,], main = bquote(sigma^2),
         ylab = "absolute change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol2, col="red", lwd=2)
    
    plot(d_details$alpha[1,], main = bquote(alpha), ylab = "relative change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol1, col="red", lwd=2)
    plot(d_details$alpha[2,], main = bquote(alpha), ylab = "absolute change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol2, col="red", lwd=2)
    
    plot(d_details$h[1,], main = bquote(h), ylab = "relative change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol1, col="red", lwd=2)
    plot(d_details$h[2,], main = bquote(h), ylab = "absolute change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol2, col="red", lwd=2)
    
    plot(d_details$beta[1,], main = bquote(beta), ylab = "relative change", 
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol1, col="red", lwd=2)
    plot(d_details$beta[2,], main = bquote(beta), ylab = "absolute change", 
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol2, col="red", lwd=2)
    
    plot(d_details$u0[1,], main = bquote(u[0]), ylab = "relative change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol1, col="red", lwd=2)
    plot(d_details$u0[2,], main = bquote(u[0]), ylab = "absolute change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol2, col="red", lwd=2)
    
    plot(d_details$u1[1,], main = bquote(u[1]), ylab = "relative change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol1, col="red", lwd=2)
    plot(d_details$u1[2,], main = bquote(u[1]), ylab = "absolute change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol2, col="red", lwd=2)
    
    plot(d_details$b[1,], main = bquote(b), ylab = "relative change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol1, col="red", lwd=2)
    plot(d_details$b[2,], main = bquote(b), ylab = "absolute change",
         xlab = "iter", cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
    abline(h=tol2, col="red", lwd=2)
    
    plot(ll, ylab = "log likelihood", xlab = "iter", 
         cex.lab=2.0,cex.axis = 2.0,cex.main = 2.0)
  } 
  
  if(calculate.scores==T){
    phi0_hat = lapply(1:p, function(x){lapply(1:n, function(i){Bi_bar[[x]][[i]]%*%u0_hat})})
    phi1_hat = lapply(1:p, function(x){lapply(1:n, function(i){Bi_bar[[x]][[i]]%*%u1_hat})})
    
    Lam0 = rbind(D0_hat, matrix(0, p*npc1, npc0))
    Lam1_bd = kronecker(diag(p), D1_hat)
    Lam1 = rbind(matrix(0, npc0, p*npc1), Lam1_bd)
    Lam = bdiag(D0_hat, Lam1_bd)
    
    cor_hat <- lapply(1:n, function(i) 
    {Reduce(cbind, sapply(1:p, function(x){tcrossprod(b_hat[x]*Lam0, phi0_hat[[x]][[i]])}, simplify = F)) +
        Reduce(cbind, sapply(1:p, function(x){tcrossprod(b_hat[x]*Lam1[, 1:npc1+(x-1)*npc1], 
                                                         phi1_hat[[x]][[i]])}, simplify = F))})
    
    
    phi0mat <- lapply(1:n, function(i) {
      as.matrix(bdiag(lapply(1:p, function(x){phi0_hat[[x]][[i]]})))})
    phi1mat <- lapply(1:n, function(i) {
      as.matrix(bdiag(lapply(1:p, function(x){phi1_hat[[x]][[i]]})))})
    D0mat <- kronecker(tcrossprod(b_hat), D0_hat)
    D1mat <- kronecker(diag(b_hat^2), D1_hat)
    cov_hat <- lapply(1:n, function(i) 
    {phi0mat[[i]]%*%D0mat%*%t(phi0mat[[i]]) + 
        phi1mat[[i]]%*%D1mat%*%t(phi1mat[[i]])})
    
    R_hat <- lapply(1:n, function(i){
      as.matrix(bdiag(lapply(1:p, function(x){
        if(sub_ind[[x]][i]==1){
          sigma2_hat[x]
        }else{
          sigma2_hat[x]*diag(sub_ind[[x]][i])}
      })))  })
    
    Vinv <- mapply(function(cov_hat, R_hat) {solve(cov_hat + R_hat)}, 
                   cov_hat = cov_hat, R_hat = R_hat, SIMPLIFY = F)
    scores_mu <-  mapply(function(cor_hat, Vinv, W){cor_hat %*% Vinv %*% W}, 
                         cor_hat = cor_hat, Vinv = Vinv, W = W, SIMPLIFY = T)
    scores_mu = t(scores_mu)
    scores_Sig <- mapply(function(cor_hat, Vinv) {
      Lam - tcrossprod(cor_hat %*% Vinv, cor_hat)}, 
      cor_hat = cor_hat, Vinv = Vinv, SIMPLIFY = F)
  }else{
    scores_mu <- NULL
    scores_Sig <- NULL
  }
  
  return(list(Theta = Theta, npc0 = npc0, npc1 = npc1, 
              face_fit = face_fit, scr_idx = scr_idx, 
              cox_beta = unname(cox_fit$coefficients),
              iter = iter, d = rbind(d1, d2), d_details = d_details, 
              loglik = loglik, dof0 = dof0, AIC0 = AIC0, BIC0 = BIC0, 
              dof = dof, AIC = AIC, BIC = BIC, 
              ll = ll, ll_rel = ll_rel, 
              I = I, I_D0 = I_D0, I_D1 = I_D1,
              scores = scores_mu, scores_Sig = scores_Sig))
  
}
