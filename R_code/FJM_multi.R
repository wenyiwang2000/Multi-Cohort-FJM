
FJM_multi <- function(longdata, time, event, X0 = NULL, w= NULL, 
                       variable, argvals.new, knots = 9,
                       npc0 = NULL, npc1 = NULL, init = NULL,
                       nMC0 = 100, nMC1 = 5000, 
                       Maxiter = 1000, Miniter = 20, rho = 0.5, 
                       tol0 = 1e-3, tol1 = 5e-3, tol2 = 1e-3, tol3 = 1e-7, 
                       seed = 721, calculate.scores = T){
  
  require(survival)
  require(MASS)
  require(Matrix)
  require(splines)
  
  ################################
  ####step 1: Initialization 
  ################################
  tick  <- proc.time()
  set.seed(seed)
  # The number of cohort
  c <- length(longdata)
  p <- mapply(function(longdata){length(longdata)},longdata=longdata)    
  n <- mapply(function(time){length(time)},time=time)
  if(is.null(w)){
    w <- rep(1,c)
  }
  #w <- n[1]/n
  p_all <- length(variable)
  n_all <- sum(n)
  N_obs_v <- rep(0,c)
  for(i in 1:c){
    for(j in 1:p[i]){
      N_obs_v[i] <- N_obs_v[i]+sum(length(longdata[[i]][[j]]$y))
    }
  }
  N_obs <- sum(N_obs_v)
  nMC <- nMC0
  tnew <- argvals.new
  if(is.null(npc0) || is.null(npc1)){
    npc0 <- init$cohort1$npc0
    npc1 <- init$cohort1$npc1
  }
  
  
  
  scores_mu <- mapply(function(init){init$scores}, init=init, SIMPLIFY = F)
  scores_Sig <- lapply(1:c,function(x){init[[x]]$scores_Sig})
  
  #eigenvalue
  D0_hat <-  rowMeans( matrix(mapply(function(init){init$Theta$D0[1:npc0]},init=init),ncol=3 ))
  D1_hat <-  rowMeans( matrix(mapply(function(init){init$Theta$D1[1:npc1]},init=init),ncol=3 ) )
  
  #theta
  u0_hat <- as.matrix(as.matrix(init[[1]]$Theta$u0)[,1:npc0]) #######
  u1_hat <- as.matrix(as.matrix(init[[1]]$Theta$u1)[,1:npc1]) #######
  
  
  Bi_bar <- lapply(1:c, function(x){ lapply(1:p[x], function(xx){ sapply(1:n[x], function(i){
    matrix((init[[x]]$B[[xx]])[which(longdata[[x]][[xx]]$subj==i),], 
           ncol=nrow(init[[x]]$G_invhalf))%*%init[[x]]$G_invhalf}) }) })
  
  alpha_hat <- lapply(1:p_all, function(x){ 
    coho <- 1:c
    p_idx <- variable[[x]][coho]
    mapply(function(coho,p_idx){
      as.matrix(init[[coho]]$Theta$alpha[,p_idx])},coho=coho,p_idx=p_idx)
  })
  
  yij = lapply(1:p_all, function(x){ #for updating alpha_hat
    coho <- 1:c
    p_idx <- variable[[x]][coho]
    lapply(1:length(coho), function(i){longdata[[coho[i]]][[p_idx[i]]] }) 
  })
  
  
  
  sigma2_hat <- sapply(1:p_all, function(x){
    coho <- which(!is.na(variable[[x]]))
    mean(unlist(lapply(coho, function(xx){init[[xx]]$Theta$sigma2[variable[[x]][xx]]} )))  })
  #beta_j
  b_hat <- sapply(1:p_all, function(x){
    coho <- which(!is.na(variable[[x]]))
    mean(unlist(lapply(coho, function(xx){init[[xx]]$Theta$b[(variable[[x]][xx])]} ))) })
  #number of observation of each subj for each variable in each cohort
  sub_ind <- lapply(1:c, function(x){
    lapply(1:p[x], function(xx){sapply(1:n[x], function(i) {sum(longdata[[x]][[xx]]$subj == i)})})
  })
  
  
  # hazards model
  scores_com <- do.call("rbind",  lapply(1:c, function(x){matrix(scores_mu[[x]][,1:npc0],ncol=npc0)}))
  time_com <- unlist(time)
  event_com <- unlist(event)
  # For all subject
  if(!is.null(X0)){
    #number of covariates
    q <- ncol(X0[[1]])
    X0_com <-  do.call("rbind", X0)
    cox_fit <- coxph(formula=Surv(time_com,event_com)~cbind(X0_com, scores_com)) ################
  }else{
    q <- 0
    cox_fit <- coxph(formula=Surv(time_com,event_com)~scores_com) ################
  }
  h_hat <- coxph.detail(cox_fit)$hazard # caution
  nevent <- coxph.detail(cox_fit)$nevent
  Tj <- coxph.detail(cox_fit)$time
  
  # For cohort to get beta_hat(Gamma)
  beta_hat <- lapply(1:c, function(x) {
    if(!is.null(X0)) {
      cox_fit <- coxph(formula=Surv(time[[x]],event[[x]])~cbind(X0[[x]], scores_mu[[x]][,1:npc0])) 
    } else{
      cox_fit <- coxph(formula=Surv(time[[x]],event[[x]])~scores_mu[[x]][,1:npc0]) 
    }
    unname(cox_fit$coefficients) })
  
  
  ######################################################
  ####step 2: Joint Estimation thru EM Algorithm (E-step)
  ######################################################
  ## main loop of EM
  scores_mu <- mapply(function(scores_mu){
    split(t(scores_mu), rep(1:nrow(scores_mu), each = ncol(scores_mu)))
  }, scores_mu=scores_mu)
  
  ## grid of event time  
  #Tj <- sort(unique(time_com[event_com != 0]))
  K = length(Tj)
  Ti_ind <- lapply(1:c, function(x){mapply(function(time){max(1, sum(Tj <= time))},time=time[[x]])})
  
  
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
    
    
    phi0_hat = lapply(1:c, function(x){ lapply(1:p[x], function(xx){
      mapply(function(Bi_bar0){Bi_bar0%*%u0_hat}, Bi_bar0=Bi_bar[[x]][[xx]], SIMPLIFY = F)})
    })
    phi1_hat = lapply(1:c, function(x){ lapply(1:p[x], function(xx){
      mapply(function(Bi_bar0){Bi_bar0%*%u1_hat}, Bi_bar0=Bi_bar[[x]][[xx]], SIMPLIFY = F)}) 
    })
    
    
    ## update multivariate normal distribution
    if(iter > 1){ # should be > 1 (was >=1)
      Lam0 = lapply(1:c, function(x){rbind(D0_hat, matrix(0, p[x]*npc1, npc0))})
      Lam1_bd = lapply(1:c, function(x){kronecker(diag(p[x]), D1_hat)})
      Lam1 = lapply(1:c, function(x){rbind(matrix(0, npc0, p[x]*npc1), Lam1_bd[[x]]) })
      Lam = lapply(1:c, function(x){bdiag(D0_hat, Lam1_bd[[x]])})
      
      cor_hat <- lapply(1:c, function(x){ lapply(1:n[x], function(i) { 
        Reduce(cbind, sapply(1:p[x], function(xx){
          p_idx <- unname(which(do.call("rbind",variable)[,x]==xx))
          tcrossprod(b_hat[p_idx]*Lam0[[x]], phi0_hat[[x]][[xx]][[i]])}, simplify = F)) + 
          Reduce(cbind, sapply(1:p[x], function(xx){
            p_idx <- unname(which(do.call("rbind",variable)[,x]==xx))
            tcrossprod(b_hat[p_idx]*Lam1[[x]][,1:npc1+(xx-1)*npc1], phi1_hat[[x]][[xx]][[i]])}, simplify = F))}) })
      
      
      phi0mat <- lapply(1:c, function(x){ lapply(1:n[x], function(i) {
        as.matrix(bdiag(lapply(1:p[x], function(xx){phi0_hat[[x]][[xx]][[i]]}))) }) })
      phi1mat <- lapply(1:c, function(x){ lapply(1:n[x], function(i) {
        as.matrix(bdiag(lapply(1:p[x], function(xx){phi1_hat[[x]][[xx]][[i]]}))) }) })
      
      temp = do.call("rbind",variable)
      b_hat_list <- lapply(1:c, function(x){b_hat[which(!is.na(temp[,x]))]})
      D0mat <- mapply(function(b_hat){kronecker(tcrossprod(b_hat), D0_hat)},
                      b_hat=b_hat_list, SIMPLIFY = F)
      D1mat <- mapply(function(b_hat){kronecker(diag(b_hat^2), D1_hat)},
                      b_hat=b_hat_list, SIMPLIFY = F)
      
      cov_hat <- lapply(1:c, function(x){ lapply(1:n[x], function(i) 
      {phi0mat[[x]][[i]]%*%D0mat[[x]]%*%t(phi0mat[[x]][[i]]) + 
          phi1mat[[x]][[i]]%*%D1mat[[x]]%*%t(phi1mat[[x]][[i]])}) })
      
      
      sigma2_hat_list <- lapply(1:c, function(x){sigma2_hat[which(!is.na(temp[,x]))]})
      R_hat <- lapply(1:c, function(x){ lapply(1:n[x], function(i){
        as.matrix(bdiag(lapply(1:p[x], function(xx){
          if(sub_ind[[x]][[xx]][i]==1){
            sigma2_hat_list[[x]][xx]
          }else{
            sigma2_hat_list[[x]][xx]*diag(sub_ind[[x]][[xx]][i])
          }  })))  }) })
      
      
      Vinv <- lapply(1:c, function(x){ # Vinv = inv(cov(yi))
        mapply(function(cov_hat, R_hat) {solve(cov_hat + R_hat)}, 
               cov_hat = cov_hat[[x]], R_hat = R_hat[[x]], SIMPLIFY = F) })
      scores_mu <- lapply(1:c, function(x){ # scores_mu = E(eta_i|y_i)
        lapply(1:n[x], function(i){
          mu0 = cor_hat[[x]][[i]] %*% Vinv[[x]][[i]] %*% W[[x]][[i]]
          if((mu0[1]*beta_hat[[x]][q+1])>8) mu0[1] = 8/(beta_hat[[x]][q+1])
          return(mu0) }) })
      scores_Sig <- lapply(1:c, function(x){ # scores_Sig = cov(eta_i|y_i)
        mapply(function(cor_hat, Vinv) { 
          Lam[[x]] - tcrossprod(cor_hat %*% Vinv, cor_hat)}, 
          cor_hat = cor_hat[[x]], Vinv = Vinv[[x]], SIMPLIFY = F) })
    }
    ## MC samples of scores
    xi = lapply(1:c, function(x){ mapply(function(mu,Sig){mvrnorm(nMC, matrix(mu), Sig)}, 
                                         mu = scores_mu[[x]], Sig = scores_Sig[[x]], SIMPLIFY = F) }) 
    
    if(!is.null(X0)){ # expX = exp(z^T_i*gamma_z + eta^T_i*gamma_eta)
      expX0 = lapply(1:c, function(x){exp(X0[[x]] %*% beta_hat[[x]][1:q])})
      expX1 = lapply(1:c, function(x){
        lapply(xi[[x]], function(xii){exp(matrix(xii[,1:npc0],ncol=npc0) %*% beta_hat[[x]][-(1:q)])}) 
      }) ################
      expX = lapply(1:c, function(x){ 
        mapply(function(expX0, expX1){expX0*expX1}, 
               expX0 = expX0[[x]], expX1 = expX1[[x]], SIMPLIFY = F) }) ################
    }else{
      expX = lapply(1:c, function(x){ lapply(xi[[x]], function(xii){exp(matrix(xii[,1:npc0],ncol=npc0) %*% beta_hat[[x]])}) }) ################
    }
    
    H = lapply(1:c, function(x) {lapply(1:n[x],function(i){sum(h_hat[1:Ti_ind[[x]][i]])}) })
    hti = lapply(1:c, function(x) {lapply(1:n[x],function(i){h_hat[Ti_ind[[x]][i]]}) })
    
    ## f(T_i, delta_i|.)
    fti = lapply(1:c, function(x) {
      mapply(function(expX, H, hti, event){
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
      }, expX = expX[[x]], H = H[[x]], hti = hti[[x]], event = event[[x]], SIMPLIFY = F) })
    
    
    ## \sum{f(T_i, delta_i|.)}/nMC : den is the denominator in E_{g(eta_i)}
    den <- lapply(1:c,function(x){lapply(fti[[x]], mean)})
    
    ## f(T_i, delta_i|.) / den
    pxi <- lapply(1:c, function(x){mapply(function(f, d) {c(f / d)}, 
                                          f = fti[[x]], d = den[[x]], SIMPLIFY = F) }) ################
    
    ## E[xi]: Exi = E_i{g(eta_i)}
    Exi <-lapply(1:c, function(x){mapply(function(xi, pxi) {colMeans(xi*pxi)},
                                         xi = xi[[x]], pxi = pxi[[x]], SIMPLIFY = F) })
    
    ## E[X]: EX = E_i{(z_i,eta_i)}
    if(!is.null(X0)){
      EX <- lapply(1:c, function(x){
        lapply(1:n[x], function(i){c(X0[[x]][i, ], Exi[[x]][[i]][1:npc0])}) }) ####################
    }else{
      EX <- lapply(1:c, function(x){
        lapply(1:n[x], function(i){Exi[[x]][[i]][1:npc0]}) })####################
    }
    
    ## E[xi^2]: Exi2 = E_i{eta_i*eta^T_i}  #################### may change it because of zeta_ij
    Exi2 <- lapply(1:c, function(x) { mapply(function(xi, pxi) 
    {crossprod(xi, xi*pxi)/nMC}, xi = xi[[x]], pxi = pxi[[x]], SIMPLIFY = F) })
    
    ## E[expX]: E_i{exp(z^T_i*gamma_z + eta^T_i*gamma_eta)}
    EexpX <- lapply(1:c, function(x) { mapply(function(pxi, expX) {mean(expX*pxi)},
                                              pxi = pxi[[x]], expX = expX[[x]], SIMPLIFY = F) })
    ## E[xi expX]: E_i{(eta_i) * exp(z^T_i*gamma_z + eta^T_i*gamma_eta)}
    ExiexpX <- lapply(1:c, function(x) { 
      mapply(function(xi, pxi, expX) {colMeans(matrix(xi[,1:npc0],ncol=npc0)*c(expX)*pxi)}, ###################
             xi = xi[[x]], pxi = pxi[[x]], expX = expX[[x]], SIMPLIFY = F) })
    
    ## E[X expX]: E_i{(z_i, eta_i) * exp(z^T_i*gamma_z + eta^T_i*gamma_eta)}
    if(!is.null(X0)){
      EXexpX <- lapply(1:c, function(x){ 
        lapply(1:n[x], function(i){c(X0[[x]][i, ]*EexpX[[x]][[i]], ExiexpX[[x]][[i]]) }) })
    }else{
      EXexpX <- ExiexpX
    }
    
    ######################################################
    ####step 3: Joint Estimation thru EM Algorithm (M-step)
    ######################################################
    ## D
    w_all <- rep(w,n)
    D0_hat <- ( do.call("cbind", lapply(1:c, function(x){ 
      matrix(mapply(function(Exi2){diag(Exi2)[1:npc0]}, Exi2=Exi2[[x]]),nrow=npc0 )})) %*% matrix(w_all,ncol = 1))/sum(w_all)
    D0_hat <- as.vector(D0_hat)
    if(npc0>1) {D0_hat<-diag(D0_hat)}
    D1_hat <- (matrix(w_all,nrow = 1) %*% do.call("rbind", lapply(1:c, function(x){ 
      do.call("rbind", mapply(function(Exi2){ colMeans(matrix(diag(Exi2)[-(1:npc0)],nc=npc1,byrow=T)) }, 
                              Exi2=Exi2[[x]], SIMPLIFY = F)) }) ) )/sum(w_all)
    D1_hat <- as.vector(D1_hat)
    if(npc1>1) {D1_hat<-diag(D1_hat)}
    

    ##1.estimate alpha
    
    PExi <- lapply(1:p_all, function(x) {
      coho <- 1:c
      p_idx <- variable[[x]][coho]
      (lapply(1:length(coho), function(xx){
        if(is.na(p_idx[xx])){NA}else{
        Exi0 = do.call("rbind",Exi[[coho[xx]]])[,c(1:npc0, npc0+(p_idx[xx]-1)*npc1+1:npc1)]
        unlist(lapply(1:n[coho[xx]], function(i){
          PExi0 = phi0_hat[[coho[xx]]][[p_idx[xx]]][[i]] %*% Exi0[i,1:npc0] 
          PExi1 = phi1_hat[[coho[xx]]][[p_idx[xx]]][[i]] %*% Exi0[i,npc0+1:npc1]
          PExi0 + PExi1 }))  }}))  })
    
    nm <- knots+3
    alpha_hat <- lapply(1:p_all, function(x){
      coho <- 1:c
      p_idx <- variable[[x]][coho]
      do.call("cbind",lapply(1:length(coho),function(xx){
        if(is.na(p_idx[xx])){NA}else{
        data_temp <- data.frame("subj"=yij[[x]][[xx]][,1], "argvals"=yij[[x]][[xx]][,2],"y"=yij[[x]][[xx]][,3]-b_hat[x]*PExi[[x]][[xx]])
        #fit.mean <- face::pspline(data=data_temp, argvals.new=tnew, knots=knots)
        #fit.mean$theta
        fitm <- gam(y ~ s(argvals, bs='ps', k=nm, m=c(2,2)), data=data_temp)
        c(pen.edf(fitm),extract.Bcoef(fitm,data_temp))}
      }))
    })
    
    df_mean <- unlist(lapply(1:p_all,function(x){alpha_hat[[x]][1,]}))
    alpha_hat <- lapply(1:p_all,function(x){alpha_hat[[x]][-1,]})
    
    
    ##2.estimate sigma2
  
    W_ind = lapply(1:c, function(x){ sapply(1:p[x], function(xx){ # yij - uij  
      p_idx <- unname(which(do.call("rbind",variable)[,x]==xx))
      as.matrix(cbind(longdata[[x]][[xx]][,1], 
                     longdata[[x]][[xx]][,3] - init[[x]]$B[[xx]]%*%alpha_hat[[p_idx]][,x]))
    }, simplify = F)  })
    
    
    W_p <- lapply(1:c, function(x){ lapply(1:n[x], function(i) { # yij - uij  
      sapply(1:p[x], function(xx) {W_ind[[x]][[xx]][W_ind[[x]][[xx]][,1] == i, -1]}, 
             simplify = F)})  })
    W <- lapply(1:c, function(x){ lapply(1:n[x], function(i){Reduce(c, W_p[[x]][[i]])}) })
    
    
    WtW = lapply(1:c, function(x) { # for sigma_j: (yij - uij)^2
      W_p0 = W_p[[x]]
      unname(Reduce(rbind, lapply(W_p0, function(W_p0) { 
        sapply(1:p[x], function(xx){crossprod(W_p0[[xx]])}) }) )) })
    WtP0 = lapply(1:c, function(x){ # for sigma_j: (yij - uij)^T phi
      lapply(1:p[x], function(xx) {Reduce(rbind, sapply(1:n[x], function(j){
        crossprod(W_p[[x]][[j]][[xx]], phi0_hat[[x]][[xx]][[j]])}, simplify = F)) })  })
    WtP1 = lapply(1:c, function(x){ # for sigma_j: (yij - uij)^T psi
      lapply(1:p[x], function(xx) {Reduce(rbind, sapply(1:n[x], function(j){
        crossprod(W_p[[x]][[j]][[xx]], phi1_hat[[x]][[xx]][[j]])}, simplify = F)) })  })
    P0tP0 = lapply(1:c, function(x){ # for sigma_j
      lapply(1:p[x], function(xx){ lapply(1:n[x], function(j) 
      {crossprod(phi0_hat[[x]][[xx]][[j]])}) }) })
    P1tP1 = lapply(1:c, function(x){ # for sigma_j
      lapply(1:p[x], function(xx){ lapply(1:n[x], function(j) 
      {crossprod(phi1_hat[[x]][[xx]][[j]])}) }) })
    P0tP1 = lapply(1:c, function(x){ # for sigma_j
      lapply(1:p[x], function(xx) {
        mapply(function(phi0_hat, phi1_hat) {crossprod(phi0_hat, phi1_hat)}, 
               phi0_hat = phi0_hat[[x]][[xx]], phi1_hat = phi1_hat[[x]][[xx]], SIMPLIFY = F) }) })
    
    
    ## sigma_j^2 & beta_j
    tr = lapply(1:c, function(x){
      t(sapply(1:n[x], function(j){ sapply(1:p[x], function(xx) {
        sum(diag(P0tP0[[x]][[xx]][[j]] %*% Exi2[[x]][[j]][1:npc0,1:npc0])) +
          sum(diag(P1tP1[[x]][[xx]][[j]] %*% Exi2[[x]][[j]][npc0+npc1*(xx-1)+1:npc1,npc0+npc1*(xx-1)+1:npc1])) +
          2*sum(diag(P0tP1[[x]][[xx]][[j]] %*% Exi2[[x]][[j]][npc0+npc1*(xx-1)+1:npc1,1:npc0])) }) })) })
    part2 = lapply(1:c, function(x){ 
      sapply(1:p[x], function(xx){ 
        p_idx <- unname(which(do.call("rbind",variable)[,x]==xx))
        sapply(1:n[x], function(j){  
          - 2*b_hat[p_idx]*WtP0[[x]][[xx]][j,]%*%Exi[[x]][[j]][1:npc0] - 
            2*b_hat[p_idx]*WtP1[[x]][[xx]][j,]%*%Exi[[x]][[j]][npc0+npc1*(xx-1)+1:npc1] + 
            b_hat[p_idx]^2*tr[[x]][j,xx] }) })  })
    
    sigma2_hat = sapply(1:p_all, function(x){
      coho <- which(!is.na(variable[[x]]))
      p_idx <- variable[[x]][coho]
      WtW_temp <- c(); part2_temp <- c(); mi=0
      for(i in 1:length(coho)){
        WtW_temp <- c(WtW_temp,w[coho[i]]*WtW[[coho[i]]][,p_idx[i]])
        part2_temp <- c(part2_temp,w[coho[i]]*part2[[coho[i]]][,p_idx[i]])
        mi <- mi + nrow(longdata[[coho[i]]] [[p_idx[i]]])*w[coho[i]]
      }
      (sum(WtW_temp) + sum(part2_temp)) / mi
    })
    
    
    ##3.estimate b
    b_x = lapply(1:c, function(x){ lapply(1:n[x], function(j) {
      sapply(1:p[x], function(xx) {
        phi0_hat[[x]][[xx]][[j]] %*% Exi[[x]][[j]][1:npc0] + phi1_hat[[x]][[xx]][[j]] %*% Exi[[x]][[j]][npc0+npc1*(xx-1)+1:npc1]
      }, simplify = F)})  })
    b_res = lapply(1:c, function(x){
      mapply(function(b_x, W) {sapply(1:p[x], function(xx) {
        crossprod(b_x[[xx]], W[[xx]])})}, b_x = b_x[[x]], W = W_p[[x]])
    })
    b_hat <- sapply(1:p_all, function(x){
      coho <- which(!is.na(variable[[x]]))
      p_idx <- variable[[x]][coho]
      b_res_temp <- c(); tr_temp <- c()
      for(i in 1:length(coho)){
        b_res_temp <- c(b_res_temp,w[coho[i]]*b_res[[coho[i]]][p_idx[i],])
        tr_temp <- c(tr_temp,w[coho[i]]*tr[[coho[i]]][,p_idx[i]])
      } 
      sum(b_res_temp) / sum(tr_temp) })
    b_hat[1] <- 1
    
    
    ##4.estimate u0, phi   
    P1Exi2 = lapply(1:c, function(x){ # Psi*E(xi*zeta)
      lapply(1:p[x], function(j){
        do.call("rbind", mapply(function(phi1_temp, Exi20){
          phi1_temp%*%t(matrix(Exi20[1:npc0, npc0+npc1*(j-1)+1:npc1],nrow=npc0))}, 
          phi1_temp=phi1_hat[[x]][[j]], Exi20=Exi2[[x]]) ) })  })
    Exi0_long <- lapply(1:c, function(x){lapply(1:p[x], function(j){
      Exi_m <- do.call("rbind",Exi[[x]])
      do.call("cbind",lapply(1:npc0, function(xx){rep(Exi_m[,xx], times=sub_ind[[x]][[j]])}))  })  })
    Exi20_long <- lapply(1:c, function(x){lapply(1:p[x], function(j){
      Exi2_m <- do.call("rbind",mapply(function(Exi20){diag(Exi20)}, Exi20=Exi2[[x]], SIMPLIFY=F))
      do.call("cbind",lapply(1:npc0, function(xx){rep(Exi2_m[,xx], times=sub_ind[[x]][[j]])}))  })  })
    Exi1_long <- lapply(1:c, function(x){lapply(1:p[x],function(j){
      Exi_m <- do.call("rbind",Exi[[x]])
      do.call("cbind",lapply(1:npc1, function(xx){rep(Exi_m[,npc0+npc1*(j-1)+xx], times=sub_ind[[x]][[j]])}))
      })  })
    Exi21_long <- lapply(1:c, function(x){lapply(1:p[x],function(j){
      Exi2_m <- do.call("rbind",mapply(function(Exi20){diag(Exi20)}, Exi20=Exi2[[x]], SIMPLIFY=F))
      do.call("cbind",lapply(1:npc1, function(xx){rep(Exi2_m[,npc0+npc1*(j-1)+xx], times=sub_ind[[x]][[j]])}))
    })  })
    
    df_theta0 <- rep(0,npc0)
    for(jj in 1:10){
      u0_pre = u0_hat
      for(k in 1:npc0){
        if(npc0 == 1){
          data_temp <- do.call("rbind", lapply(1:c, function(x){
            do.call("rbind", lapply(1:p[x], function(j){
              p_idx <- unname(which(do.call("rbind",variable)[,x]==j))
              b_sqExi2_P1Exi2 <- (b_hat[p_idx]/sqrt(Exi20_long[[x]][[j]]))*P1Exi2[[x]][[j]] #(beta_j/sq(Exi2))*Psi*Exi2
              Exi_sqExi2_Wij <- (Exi0_long[[x]][[j]]/sqrt(Exi20_long[[x]][[j]]))*W_ind[[x]][[j]][,2] #(Exi/sq(Exi2))*(y_ijk-u_ijk)
              
              y_temp <- Exi_sqExi2_Wij - b_sqExi2_P1Exi2
              bsqExi2 <- b_hat[p_idx] * sqrt(Exi20_long[[x]][[j]]) #bj*sqrt(Exi2)
              argvals <- longdata[[x]][[j]][,2]
              weight <- rep(w[x],length(argvals))
              cbind(argvals, y_temp, bsqExi2,weight)  })) }) )
          } else {
            data_temp <- do.call("rbind", lapply(1:c, function(x){
              do.call("rbind", lapply(1:p[x], function(j){
                p_idx <- unname(which(do.call("rbind",variable)[,x]==j))
                phi0_temp <- mapply(function(Bi_bar0){Bi_bar0%*%u0_hat}, Bi_bar0=Bi_bar[[x]][[j]], SIMPLIFY=F)
                P0Exi2 <- do.call("rbind", mapply(function(phi0_temp, Exi20){t(t(phi0_temp)*Exi20[1:npc0,1:npc0][k,])}, 
                                                 phi0_temp=phi0_temp, Exi20=Exi2[[x]])) #Phi*Exi2
                b_sqExi2_P0Exi2 <- (b_hat[p_idx]/sqrt(Exi20_long[[x]][[j]][,k]))*(rowSums(P0Exi2)-P0Exi2[,k]) #(beta_j/sq(Exi2))*Phi*Exi2
                b_sqExi2_P1Exi2 <- (b_hat[p_idx]/sqrt(Exi20_long[[x]][[j]][,k]))*P1Exi2[[x]][[j]][,k] #(beta_j/sq(Exi2))*Psi*Exi2
                Exi_sqExi2_Wij <- (Exi0_long[[x]][[j]][,k]/sqrt(Exi20_long[[x]][[j]][,k]))*W_ind[[x]][[j]][,2] #(Exi/sq(Exi2))*(y_ijk-u_ijk)
                
                y_temp <- Exi_sqExi2_Wij - b_sqExi2_P1Exi2 - b_sqExi2_P0Exi2
                bsqExi2 <- b_hat[p_idx] * sqrt(Exi20_long[[x]][[j]][,k]) #bj*sqrt(Exi2)
                argvals <- longdata[[x]][[j]][,2]
                weight <- rep(w[x],length(argvals))
                cbind(argvals, y_temp, bsqExi2,weight)  })) }) )
          }
        
        data_temp <- data.frame(data_temp)
        colnames(data_temp) <- c("argvals","y_temp","bsqExi2","weight")
        fit <- gam(y_temp ~ s(argvals, by=bsqExi2, bs='ps', k=knots+3, m=c(2,2)), drop.intercept=T, data=data_temp,weights = weight)
        #fit <- gam(y_temp ~ s(argvals, by=bsqExi2, bs='ps', k=knots+3, m=c(2,2)), sp=0, drop.intercept=T, data=data_temp)
        u0_hat[,k] <- init[[1]]$G_half %*% fit$coefficients
        df_theta0[k] <- pen.edf(fit)
      }
      if(max(abs(u0_pre - u0_hat)) < 1e-6) break
    }
    
    
    eig0 <- eigen(tcrossprod(u0_hat %*% D0_hat, u0_hat))
    if(npc0 == 1){
      D0_hat <- diag(as.matrix(eig0$values[1:npc0]))
    }else{
      D0_hat <- diag(eig0$values[1:npc0]) # D0_hat = D0
    }
    
    for(k in 1:npc0){
      if(crossprod(u0_hat[,k], eig0$vectors[,k]) < 0){
        u0_hat[,k] <- - eig0$vectors[,k]
      }else{
        u0_hat[,k] <- eig0$vectors[,k] 
      }
    }
    
    ##5.estimate u1, psi
    phi0_hat = lapply(1:c, function(x){ lapply(1:p[x], function(xx){
      mapply(function(Bi_bar0){Bi_bar0%*%u0_hat}, Bi_bar0=Bi_bar[[x]][[xx]], SIMPLIFY = F)}) })
    P0Exi2 = lapply(1:c, function(x){ # Phi*E(xi*zeta)
      lapply(1:p[x], function(j){
        do.call("rbind", mapply(function(phi0_temp, Exi20){
          phi0_temp%*%t(matrix(Exi20[npc0+npc1*(j-1)+1:npc1,1:npc0],ncol=npc0))}, 
          phi0_temp=phi0_hat[[x]][[j]], Exi20=Exi2[[x]]) ) })  })
    
    df_theta1 <- rep(0,npc1)
    for(jj in 1:10){
      u1_pre = u1_hat
      for(k in 1:npc1){
        if(npc1 == 1){
          data_temp <- do.call("rbind", lapply(1:c, function(x){
            do.call("rbind", lapply(1:p[x], function(j){
              p_idx <- unname(which(do.call("rbind",variable)[,x]==j))
              b_sqExi2_P0Exi2 <- (b_hat[p_idx]/sqrt(Exi21_long[[x]][[j]]))*P0Exi2[[x]][[j]] #(beta_j/sq(Exi2))*Phi*Exi2
              Exi_sqExi2_Wij <- (Exi1_long[[x]][[j]][,k]/sqrt(Exi21_long[[x]][[j]]))*W_ind[[x]][[j]][,2] #(Exi/sq(Exi2))*(y_ijk-u_ijk)
              
              y_temp <- Exi_sqExi2_Wij - b_sqExi2_P0Exi2 
              bsqExi2 <- b_hat[p_idx] * sqrt(Exi21_long[[x]][[j]]) #bj*sqrt(Exi2)
              argvals <- longdata[[x]][[j]][,2]
              weight <- rep(w[x],length(argvals))
              cbind(argvals, y_temp, bsqExi2,weight)  })) }) )
        }else{
          data_temp <- do.call("rbind", lapply(1:c, function(x){
            do.call("rbind", lapply(1:p[x], function(j){
              p_idx <- unname(which(do.call("rbind",variable)[,x]==j))
              phi1_temp <- mapply(function(Bi_bar0){Bi_bar0%*%u1_hat}, Bi_bar0=Bi_bar[[x]][[j]], SIMPLIFY=F)
              P1Exi2 <- do.call("rbind", mapply(function(phi1_temp, Exi20){t(t(phi1_temp)*Exi20[npc0+npc1*(j-1)+1:npc1,npc0+npc1*(j-1)+1:npc1][k,])}, 
                                                phi1_temp=phi1_temp, Exi20=Exi2[[x]])) #Psi*Exi2
              b_sqExi2_P1Exi2 <- (b_hat[p_idx]/sqrt(Exi21_long[[x]][[j]][,k]))*(rowSums(P1Exi2)-P1Exi2[,k]) #(beta_j/sq(Exi2))*Psi*Exi2
              b_sqExi2_P0Exi2 <- (b_hat[p_idx]/sqrt(Exi21_long[[x]][[j]][,k]))*P0Exi2[[x]][[j]][,k] #(beta_j/sq(Exi2))*Phi*Exi2
              Exi_sqExi2_Wij <- (Exi1_long[[x]][[j]][,k]/sqrt(Exi21_long[[x]][[j]][,k]))*W_ind[[x]][[j]][,2] #(Exi/sq(Exi2))*(y_ijk-u_ijk)
              
              y_temp <- Exi_sqExi2_Wij - b_sqExi2_P0Exi2 - b_sqExi2_P1Exi2
              bsqExi2 <- b_hat[p_idx] * sqrt(Exi21_long[[x]][[j]][,k]) #bj*sqrt(Exi2)
              argvals <- longdata[[x]][[j]][,2]
              weight <- rep(w[x],length(argvals))
              cbind(argvals, y_temp, bsqExi2,weight)  })) }) )
        }

        data_temp <- data.frame(data_temp)
        colnames(data_temp) <- c("argvals","y_temp","bsqExi2","weight")
        fit <- gam(y_temp ~ s(argvals, by=bsqExi2, bs='ps', k=knots+3, m=c(2,2)), drop.intercept=T, data=data_temp,weights = weight)
        #fit <- gam(y_temp ~ s(argvals, by=bsqExi2, bs='ps', k=knots+3, m=c(2,2)), sp=0, drop.intercept=T, data=data_temp)
        u1_hat[,k] <- init[[1]]$G_half %*% fit$coefficients
        df_theta1[k] <- pen.edf(fit)
      }
      if(max(abs(u1_pre - u1_hat)) < 1e-6) break
    }
    
    eig1 <- eigen(tcrossprod(u1_hat %*% D1_hat, u1_hat))
    if(npc1 == 1){
      D1_hat <- diag(as.matrix(eig1$values[1:npc1]))
    }else{
      D1_hat <- diag(eig1$values[1:npc1]) # D1_hat = D1
    }

    for(k in 1:npc1){
      if(crossprod(u1_hat[,k], eig1$vectors[,k]) < 0){
        u1_hat[,k] <- - eig1$vectors[,k]
      }else{
        u1_hat[,k] <- eig1$vectors[,k] 
      }
    }
    
    
    ##6.estimate h0
    h_de <- sapply(Tj, function(x){sum(w_all*unlist(EexpX)*(time_com>=x))}) 
    h_nu <- sapply(1:length(nevent),function(x){sum(w_all*event_com*(time_com==coxph.detail(cox_fit)$time[x]))})
    h_hat <- h_nu / h_de
    
    
    ## score of beta
    Si <- lapply(1:c, function(x){
      mapply(function(event, EX, H, EXexpX) {event*EX - H*EXexpX}, 
             event=event[[x]], EX=EX[[x]], H=H[[x]], EXexpX=EXexpX[[x]], SIMPLIFY=F) })
    S <- lapply(1:c, function(x){Reduce("+", Si[[x]])})
    
    ## information matrix
    I <- lapply(1:c, function(x){
      Reduce("+", lapply(Si[[x]], tcrossprod)) - tcrossprod(S[[x]]) / n[x] })
    ##7.estimate beta
    gBeta = lapply(1:c, function(x){rho * solve(I[[x]], S[[x]])})
    beta_hat <- lapply(1:c, function(x){beta_hat[[x]] + gBeta[[x]]})
    
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
    
    #print(iter) # <------------------------
    
    # marginal (observed data) likelihood 
    if(iter > 1){
      loglik <- lapply(1:c, function(x){ mapply(function(W, cov_hat, R_hat, Vinv) {
        - 0.5 * (length(W) * log(2*pi) + 
                   as.numeric(determinant(cov_hat+R_hat, logarithm = TRUE)$modulus) +
                   (t(W) %*% Vinv %*% W)[1,1])}, 
        W = W[[x]], cov_hat = cov_hat[[x]], R_hat = R_hat[[x]], Vinv = Vinv[[x]]) })
      loglik <-  sum(unlist(loglik)*w_all + log(unlist(den))*w_all)
      ll <- c(ll, loglik)
      ll_rel <- c(ll_rel, abs(loglik-ll[length(ll)-1]) / 
                    (abs(ll[length(ll)-1]) + tol0))
    }
    
    
    ## convergence control (maximum absoulte relative change)
    Theta0 <- Theta
    Theta0$beta <- unlist(Theta$beta )
    Theta0$alpha <- unlist(Theta$alpha)
    Theta_new0 <- Theta_new
    Theta_new0$beta <- unlist(Theta_new$beta )
    Theta_new0$alpha <- unlist(Theta_new$alpha )
    d_rel <- mapply(function(Theta, Theta_new) {
      abs(Theta - Theta_new) / (abs(Theta) + tol0)}, 
      Theta = Theta0, Theta_new = Theta_new0)
    d_abs <- mapply(function(Theta, Theta_new) {
      abs(Theta - Theta_new)}, 
      Theta = Theta0, Theta_new = Theta_new0)
    max_rel <- max(unlist(d_rel),na.rm = T)
    max_abs <- max(unlist(d_abs),na.rm = T)
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
    alpha_d <- c(alpha_d, c(max(d_rel$alpha,na.rm = T), max(d_abs$alpha,na.rm = T)))
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
  dof_alpha <- sum(sapply(1:p_all,function(x){prod(dim(matrix(Theta$alpha[[x]][,colSums(is.na(Theta$alpha[[x]]))==0])))}))
  dof0 <- dof_alpha + (npc0 + npc1) * (nrow(Theta$u0) + 1) + c*npc0 + 2 * p_all - 1
  dofe <- sum(df_mean[!is.na(df_mean)]) + sum(df_theta0) + sum(df_theta1) + (npc0 + npc1) + c * npc0 + 2 * p_all - 1 
  if(!is.null(X0)){
    dof0 <- dof0 + c*q
    dofe <- dofe + c * q
  }
  
  dofnew <- rep(0,3)
  dofnew[1] <- sum(na.omit(df_mean[(1:p_all)*c-2]))+(sum(df_theta0)+sum(df_theta1)- npc0 * (npc0 - 1) / 2 - npc1 * (npc1 - 1) / 2)/c+(2 * p_all - 1)/c + npc0+q
  dofnew[2] <- sum(na.omit(df_mean[(1:p_all)*c-1]))+(sum(df_theta0)+sum(df_theta1)- npc0 * (npc0 - 1) / 2 - npc1 * (npc1 - 1) / 2)/c+(2 * p_all - 1)/c + npc0+q
  dofnew[3] <- sum(na.omit(df_mean[(1:p_all)*c]))+(sum(df_theta0)+sum(df_theta1)- npc0 * (npc0 - 1) / 2 - npc1 * (npc1 - 1) / 2)/c+(2 * p_all - 1)/c + npc0+q
  
  AIC0 <- -2 * loglik + 2 * dof0
  BIC0 <- -2 * loglik + sum(w*log(N_obs_v)) * dof0
  
  dof <- dof0 - npc0 * (npc0 + 1) / 2 - npc1 * (npc1 + 1) / 2 
  dofe <- dofe - npc0 * (npc0 + 1) / 2 - npc1 * (npc1 + 1) / 2 
  AIC <- -2 * loglik + 2 * dof
  BIC <- -2 * loglik + sum(w*log(N_obs_v)) * dof
  
  AICe <- -2 * loglik + 2 * dofe
  BICe <- -2 * loglik + sum(w*log(N_obs_v)) * dofe
  
  AICnew <- -2 * loglik + sum(w*dofnew)*2
  BICnew <- -2 * loglik + sum(w*log(N_obs_v)*dofnew)
  
  if(calculate.scores==T){
    phi0_hat = lapply(1:c, function(x){ lapply(1:p[x], function(xx){
      lapply(1:n[x], function(i){Bi_bar[[x]][[xx]][[i]]%*%u0_hat})}) })
    phi1_hat = lapply(1:c, function(x){ lapply(1:p[x], function(xx){
      lapply(1:n[x], function(i){Bi_bar[[x]][[xx]][[i]]%*%u1_hat})}) })
    
    Lam0 = lapply(1:c, function(x){rbind(D0_hat, matrix(0, p[x]*npc1, npc0))})
    Lam1_bd = lapply(1:c, function(x){kronecker(diag(p[x]), D1_hat)})
    Lam1 = lapply(1:c, function(x){rbind(matrix(0, npc0, p[x]*npc1), Lam1_bd[[x]]) })
    Lam = lapply(1:c, function(x){bdiag(D0_hat, Lam1_bd[[x]])})
    
    cor_hat <- lapply(1:c, function(x){ lapply(1:n[x], function(i) { 
      Reduce(cbind, sapply(1:p[x], function(xx){
        p_idx <- unname(which(do.call("rbind",variable)[,x]==xx))
        tcrossprod(b_hat[p_idx]*Lam0[[x]], phi0_hat[[x]][[xx]][[i]])}, simplify = F)) + 
        Reduce(cbind, sapply(1:p[x], function(xx){
          p_idx <- unname(which(do.call("rbind",variable)[,x]==xx))
          tcrossprod(b_hat[p_idx]*Lam1[[x]][,1:npc1+(xx-1)*npc1], phi1_hat[[x]][[xx]][[i]])}, simplify = F))}) })
    
    phi0mat <- lapply(1:c, function(x){ lapply(1:n[x], function(i) {
      as.matrix(bdiag(lapply(1:p[x], function(xx){phi0_hat[[x]][[xx]][[i]]}))) }) })
    phi1mat <- lapply(1:c, function(x){ lapply(1:n[x], function(i) {
      as.matrix(bdiag(lapply(1:p[x], function(xx){phi1_hat[[x]][[xx]][[i]]}))) }) })
    
    temp = do.call("rbind",variable)
    b_hat_list <- lapply(1:c, function(x){b_hat[which(!is.na(temp[,x]))]})
    D0mat <- mapply(function(b_hat){kronecker(tcrossprod(b_hat), D0_hat)},
                    b_hat=b_hat_list, SIMPLIFY = F)
    D1mat <- mapply(function(b_hat){kronecker(diag(b_hat^2), D1_hat)},
                    b_hat=b_hat_list, SIMPLIFY = F)
    cov_hat <- lapply(1:c, function(x){ lapply(1:n[x], function(i) 
    {phi0mat[[x]][[i]]%*%D0mat[[x]]%*%t(phi0mat[[x]][[i]]) + 
        phi1mat[[x]][[i]]%*%D1mat[[x]]%*%t(phi1mat[[x]][[i]])}) })
    
    sigma2_hat_list <- lapply(1:c, function(x){sigma2_hat[which(!is.na(temp[,x]))]})
    R_hat <- lapply(1:c, function(x){ lapply(1:n[x], function(i){
      as.matrix(bdiag(lapply(1:p[x], function(xx){
        if(sub_ind[[x]][[xx]][i]==1){
          sigma2_hat_list[[x]][xx]
        }else{
          sigma2_hat_list[[x]][xx]*diag(sub_ind[[x]][[xx]][i])
        }  })))  }) })
    
    
    Vinv <- lapply(1:c, function(x){ # Vinv = inv(cov(yi))
      mapply(function(cov_hat, R_hat) {solve(cov_hat + R_hat)}, 
             cov_hat = cov_hat[[x]], R_hat = R_hat[[x]], SIMPLIFY = F) })
    scores_mu <- lapply(1:c, function(x){ # scores_mu = E(eta_i|y_i)
      mapply(function(cor_hat, Vinv, W){cor_hat %*% Vinv %*% W}, 
             cor_hat = cor_hat[[x]], Vinv = Vinv[[x]], W = W[[x]], SIMPLIFY = T) })
    scores_Sig <- lapply(1:c, function(x){ # scores_Sig = cov(eta_i|y_i)
      mapply(function(cor_hat, Vinv) { 
        Lam[[x]] - tcrossprod(cor_hat %*% Vinv, cor_hat)}, 
        cor_hat = cor_hat[[x]], Vinv = Vinv[[x]], SIMPLIFY = F) })
  }else{
    scores_mu <- NULL
    scores_Sig <- NULL
  }
  tock <- proc.time()
  
  return(list(Theta = Theta, npc0 = npc0, npc1 = npc1,
              cox_beta = unname(cox_fit$coefficients),
              iter = iter, d = rbind(d1, d2), d_details = d_details, 
              loglik = loglik, dof0 = dof0, AIC0 = AIC0, BIC0 = BIC0, dofnew=dofnew,AICnew=AICnew,BICnew=BICnew,
              df_mean=df_mean, df_theta0=df_theta0, df_theta1=df_theta1, N_obs_v=N_obs_v,w=w,
              dof = dof,dofe = dofe, AIC = AIC, BIC = BIC,AICe = AICe, BICe = BICe, I = I, 
              ll = ll, ll_rel = ll_rel, init=init,
              scores = scores_mu, scores_Sig = scores_Sig,run_time=tock-tick))
  
}
