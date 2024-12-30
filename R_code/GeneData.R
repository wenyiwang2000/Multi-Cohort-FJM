data_gen <- function(n=c(715,3707,522),phi,psi,mu,beta,D0,D1,sigma2,gamma,X0,npc0,npc1,
                     alpha0=c(2.7,1.4,1.8),beta0=c(0.95,0.7,0.6),lambda0=1,rho0=20){
  require(MASS)
  # longitudinal data
  tnew <- seq(0,1,1/100)
  m = length(tnew)
  
  
  Sigma0 <- diag(D0,nrow=npc0,ncol=npc0) #var of xi
  Sigma1 <- diag(D1,nrow=npc1,ncol=npc1) #var of zeta
  mean0 <- rep(0,nrow(Sigma0))  #mean of xi
  mean1 <- rep(0,nrow(Sigma1)) #mean of zeta
  
  #xi
  xi <- lapply(1:3,function(x){mvrnorm(n[x], mean0, Sigma0)})
  U <- lapply(1:3,function(x){xi[[x]] %*% phi})
  
  #zeta
  zeta1 <- lapply(1:3,function(x){mvrnorm(n[x], mean1, Sigma1)})
  V1 <- lapply(1:3,function(x){zeta1[[x]] %*% psi})
  zeta2 <- lapply(1:3,function(x){mvrnorm(n[x], mean1, Sigma1)})
  V2 <- lapply(1:3,function(x){zeta2[[x]] %*% psi})
  zeta3 <- lapply(1:3,function(x){mvrnorm(n[x], mean1, Sigma1)})
  V3 <- lapply(1:3,function(x){zeta3[[x]] %*% psi})
  zeta4 <- lapply(1:3,function(x){mvrnorm(n[x], mean1, Sigma1)})
  V4 <- lapply(1:3,function(x){zeta4[[x]] %*% psi})
  zeta5 <- lapply(1:3,function(x){mvrnorm(n[x], mean1, Sigma1)})
  V5 <- lapply(1:3,function(x){zeta5[[x]] %*% psi})
  zeta6 <- lapply(1:3,function(x){mvrnorm(n[x], mean1, Sigma1)})
  V6 <- lapply(1:3,function(x){zeta6[[x]] %*% psi})
  zeta7 <- lapply(1:3,function(x){mvrnorm(n[x], mean1, Sigma1)})
  V7 <- lapply(1:3,function(x){zeta7[[x]] %*% psi})
  zeta8 <- lapply(1:3,function(x){mvrnorm(n[x], mean1, Sigma1)})
  V8 <- lapply(1:3,function(x){zeta8[[x]] %*% psi})
  
  #epsilon
  e1 <- lapply(1:3,function(x){matrix(rnorm(n[x]*m,0,sqrt(sigma2[1])),n[x],m)})
  e2 <- lapply(1:3,function(x){matrix(rnorm(n[x]*m,0,sqrt(sigma2[2])),n[x],m)})
  e3 <- lapply(1:3,function(x){matrix(rnorm(n[x]*m,0,sqrt(sigma2[3])),n[x],m)})
  e4 <- lapply(1:3,function(x){matrix(rnorm(n[x]*m,0,sqrt(sigma2[4])),n[x],m)})
  e5 <- lapply(1:3,function(x){matrix(rnorm(n[x]*m,0,sqrt(sigma2[5])),n[x],m)})
  e6 <- lapply(1:3,function(x){matrix(rnorm(n[x]*m,0,sqrt(sigma2[6])),n[x],m)})
  e7 <- lapply(1:3,function(x){matrix(rnorm(n[x]*m,0,sqrt(sigma2[7])),n[x],m)})
  e8 <- lapply(1:3,function(x){matrix(rnorm(n[x]*m,0,sqrt(sigma2[8])),n[x],m)})
  
  #mu
  mu1 <- lapply(1:3,function(x){matrix(rep(mu[[1]][,x],n[x]),nrow=n[x],ncol=m,byrow = T)})
  mu2 <- lapply(1:3,function(x){matrix(rep(mu[[2]][,x],n[x]),nrow=n[x],ncol=m,byrow = T)})
  mu3 <- lapply(1:3,function(x){matrix(rep(mu[[3]][,x],n[x]),nrow=n[x],ncol=m,byrow = T)})
  mu4 <- lapply(1:3,function(x){matrix(rep(mu[[4]][,x],n[x]),nrow=n[x],ncol=m,byrow = T)})
  mu5 <- lapply(1:3,function(x){matrix(rep(mu[[5]][,x],n[x]),nrow=n[x],ncol=m,byrow = T)})
  mu6 <- lapply(1:3,function(x){matrix(rep(mu[[6]][,x],n[x]),nrow=n[x],ncol=m,byrow = T)})
  mu7 <- lapply(1:3,function(x){matrix(rep(mu[[7]][,x],n[x]),nrow=n[x],ncol=m,byrow = T)})
  mu8 <- lapply(1:3,function(x){matrix(rep(mu[[8]][,x],n[x]),nrow=n[x],ncol=m,byrow = T)})
  
  # observations
  Y1 <- lapply(1:3,function(x){mu1[[x]]+beta[1]*(U[[x]]+V1[[x]])+e1[[x]]})
  Y2 <- lapply(1:3,function(x){mu2[[x]]+beta[2]*(U[[x]]+V2[[x]])+e2[[x]]})
  Y3 <- lapply(1:3,function(x){mu3[[x]]+beta[3]*(U[[x]]+V3[[x]])+e3[[x]]})
  Y4 <- lapply(1:3,function(x){mu4[[x]]+beta[4]*(U[[x]]+V4[[x]])+e4[[x]]})
  Y5 <- lapply(1:3,function(x){mu5[[x]]+beta[5]*(U[[x]]+V5[[x]])+e5[[x]]})
  Y6 <- lapply(1:3,function(x){mu6[[x]]+beta[6]*(U[[x]]+V6[[x]])+e6[[x]]})
  Y7 <- lapply(1:3,function(x){mu7[[x]]+beta[7]*(U[[x]]+V7[[x]])+e7[[x]]})
  Y8 <- lapply(1:3,function(x){mu8[[x]]+beta[8]*(U[[x]]+V8[[x]])+e8[[x]]})
  
  
  # time-to-event data
  eta <- lapply(1:3,function(x){exp(X0[[x]] %*% gamma[1:4,x]+xi[[x]] %*% gamma[5:6,x])})
  v <- lapply(1:3,function(x){runif(n[x],0,1)})
  failuretime <- lapply(1:3,function(x){(-log(v[[x]])/(lambda0*eta[[x]]))^(1/rho0)})
  censoringtime <- lapply(1:3,function(x){rbeta(n[x],alpha0[x],beta0[x])}) #beta distribution
  
  censoringtime <- lapply(1:3,function(x){
    sapply(1:n[x],function(i){min(1,censoringtime[[x]][i])})
  })
  
  event <- lapply(1:3,function(x){as.numeric(failuretime[[x]] < censoringtime[[x]])})
  rate <- lapply(1:3,function(x){mean(event[[x]])})
  time <- lapply(1:3,function(x){
    failuretime[[x]] * event[[x]] + censoringtime[[x]] * (rep(1,n[x]) - event[[x]])})
  
  ## sparse observations (for now)
  
  index = c(1,11,21,31,41,51,61,71,81,91,101)
  
  data <- lapply(1:3,function(x){
    subj <- c()
    t <- c()
    y1 <- c()
    y2 <- c()
    y3 <- c()
    y4 <- c()
    y5 <- c()
    y6 <- c()
    y7 <- c()
    y8 <- c()
    obs_time <- c()
    idx <- c()
    for(i in 1:n[x]){
      subj <- c(subj,rep(i,length(index)))
      obs_time <- c(obs_time,rep(time[[x]][i],length(index)))
      t <- c(t,tnew[index])
      y1 <- c(y1,Y1[[x]][i,index])
      y2 <- c(y2,Y2[[x]][i,index])
      y3 <- c(y3,Y3[[x]][i,index])
      y4 <- c(y4,Y4[[x]][i,index])
      y5 <- c(y5,Y5[[x]][i,index])
      y6 <- c(y6,Y6[[x]][i,index])
      y7 <- c(y7,Y7[[x]][i,index])
      y8 <- c(y8,Y8[[x]][i,index])
      idx <- c(idx, index)
    }
    data0 <- data.frame("subj"=subj, "argvals" = t, "y1" = y1, "y2" = y2,
                        "y3"=y3,"y4"=y4,"y5"=y5,"y6"=y6,"y7"=y7,"y8"=y8,
                        "obs_time"=obs_time,"idx"=idx)
    data0 <- data0[data0$argvals<=data0$obs_time,]
    colnames(data0) <- c("ID","argvals","MMSE","CDRSUM","WMS","ADAS","RAVLT",
                         "FAQ","TRAILA","SDMT","obs_time","idx")
    data0
  })
  
  data_subj <- lapply(1:3,function(x){
    data0_subj <- cbind(1:n[x],X0[[x]],event[[x]],time[[x]])
    colnames(data0_subj) <- c("ID","Age","Gender","Edu","APOE4","event","obs_time")
    data0_subj
  })
  
  
  return(list(xi=xi, U = U, zeta1=zeta1, V1=V1,zeta2=zeta2, V2=V2,
              zeta3=zeta3, V3=V3, zeta4=zeta4, V4=V4,zeta5=zeta5, V5=V5,
              zeta6=zeta6, V6=V6,zeta7=zeta7, V7=V7,zeta8=zeta8, V8=V8,
              e1=e1,e2=e2,e3=e3,e4=e4,e5=e5,e6=e6,e7=e7,e8=e8,
              alpha0=alpha0,beta0=beta0,eta=eta,rate=rate,
              tnew = tnew, m = m, data=data, data_subj=data_subj))
  
}
