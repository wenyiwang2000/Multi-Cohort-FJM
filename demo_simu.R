rm(list=ls())
require(survival)
require(Matrix)
require(splines)
require(joineRML)
require(face)
require(mgcv)
require(MASS)

# Load functions
source('./R_code/GeneData.R')
source("./R_code/sFPCA_post.R")
source("./R_code/predict.face.sparse.inner_pre.R")
source("./R_code/mface.sparse2.R")
source("./R_code/predict.MMFPCA2.R")
source("./R_code/twostep2.R")
source("./R_code/extract.Bcoef.R")
source('./R_code/FJM.R') 
source("./R_code/FJM_multi.R") # multi-cohort FJM 

# Generate data 
load("./data/initial_data_gen.RData")
set.seed(1)
out <- data_gen(n=c(715,3707,522),phi=initial_data_gen$phi,psi=initial_data_gen$psi,
                mu=initial_data_gen$mu,beta=initial_data_gen$beta,
                D0=initial_data_gen$D0,D1=initial_data_gen$D1,
                sigma2=initial_data_gen$sigma2,gamma=initial_data_gen$gamma,
                X0=initial_data_gen$X0,npc0=2,npc1=2,
                alpha0=c(2.7,1.4,1.8),beta0=c(0.95,0.7,0.6),lambda0=1,rho0=20)
data <- out$data
data_subj <- out$data_subj

tnew <- seq(0,1,1/100)
m <- length(tnew)
nm <- 7 # number of basis function

# ADNI
ADNI <- data[[1]]
ADNI_subj <- as.data.frame(data_subj[[1]])
n_ADNI <- nrow(ADNI_subj) 
time_ADNI <- ADNI_subj$obs_time
event_ADNI <- ADNI_subj$event

MMSE <- data.frame("subj"=ADNI$ID, "argvals"=ADNI$argvals, "y"=ADNI$MMSE)
CDRSB <- data.frame("subj"=ADNI$ID, "argvals"=ADNI$argvals, "y"=ADNI$CDRSUM)
WMSLM <- data.frame("subj"=ADNI$ID, "argvals"=ADNI$argvals, "y"=ADNI$WMS)
ADAS13 <- data.frame("subj"=ADNI$ID, "argvals"=ADNI$argvals, "y"=ADNI$ADAS)
RAVLT_imm <- data.frame("subj"=ADNI$ID, "argvals"=ADNI$argvals, "y"=ADNI$RAVLT)
FAQ <- data.frame("subj"=ADNI$ID, "argvals"=ADNI$argvals, "y"=ADNI$FAQ)

X0_ADNI = cbind(ADNI_subj$Age, ADNI_subj$Gender, ADNI_subj$Edu, ADNI_subj$APOE4)
longdata_ADNI <- list("y1"=MMSE, "y2"=CDRSB, "y3"=WMSLM, "y4"=ADAS13, "y5"=RAVLT_imm, "y6"=FAQ) 
n_obs_ADNI <- Reduce("+",lapply(longdata_ADNI,function(x){length(x$y)}))
face_fit_ADNI <- twostep2(longdata = longdata_ADNI, knots = 4, argvals.new = tnew)

# NACC
NACC <- data[[2]]
NACC_subj <- as.data.frame(data_subj[[2]])
n_NACC <- nrow(NACC_subj)
time_NACC <- NACC_subj$obs_time
event_NACC <- NACC_subj$event

MMSE <- data.frame("subj"=NACC$ID, "argvals"=NACC$argvals, "y"=NACC$MMSE)
CDRSB <- data.frame("subj"=NACC$ID, "argvals"=NACC$argvals, "y"=NACC$CDRSUM)
WMSLM <- data.frame("subj"=NACC$ID, "argvals"=NACC$argvals, "y"=NACC$WMS)
TRAILA <- data.frame("subj"=NACC$ID, "argvals"=NACC$argvals, "y"=NACC$TRAILA)

X0_NACC <- cbind(NACC_subj$Age, NACC_subj$Gender, NACC_subj$Edu, NACC_subj$APOE4)
longdata_NACC <- list("y1"=MMSE, "y2"=CDRSB, "y3"=WMSLM, "y7"=TRAILA)
longdata_NACC_i <- list("y1"=MMSE, "y2"=CDRSB, "y3"=WMSLM, "y4"=TRAILA)
n_obs_NACC <- Reduce("+",lapply(longdata_NACC,function(x){length(x$y)}))
face_fit_NACC <- twostep2(longdata = longdata_NACC_i, knots = 4, argvals.new = tnew)


# ROSMAP
ROSMAP <- data[[3]]
ROSMAP_subj <- as.data.frame(data_subj[[3]])
n_ROSMAP <- nrow(ROSMAP_subj)
time_ROSMAP <- as.numeric(ROSMAP_subj$obs_time)
event_ROSMAP <- ROSMAP_subj$event

MMSE <- data.frame("subj"=ROSMAP$ID, "argvals"=ROSMAP$argvals, "y"=ROSMAP$MMSE)
WMSLM <- data.frame("subj"=ROSMAP$ID, "argvals"=ROSMAP$argvals, "y"=ROSMAP$WMS)
SDMT <- data.frame("subj"=ROSMAP$ID, "argvals"=ROSMAP$argvals, "y"=ROSMAP$SDMT)

X0_ROSMAP <- cbind(ROSMAP_subj$Age, ROSMAP_subj$Gender, ROSMAP_subj$Edu, ROSMAP_subj$APOE4)
longdata_ROSMAP <- list("y1"=MMSE, "y3"=WMSLM, "y8"=SDMT)
longdata_ROSMAP_i <- list("y1"=MMSE, "y2"=WMSLM, "y3"=SDMT)
n_obs_ROSMAP <- Reduce("+",lapply(longdata_ROSMAP,function(x){length(x$y)}))
face_fit_ROSMAP <- twostep2(longdata = longdata_ROSMAP_i, knots = 4, argvals.new = tnew)



longdata = list(cohort1=longdata_ADNI,  cohort2=longdata_NACC, cohort3=longdata_ROSMAP)
time = list(cohort1=time_ADNI, cohort2=time_NACC, cohort3=time_ROSMAP)
event = list(cohort1=event_ADNI, cohort2=event_NACC, cohort3=event_ROSMAP)
X0 = list(cohort1=X0_ADNI, cohort2=X0_NACC, cohort3=X0_ROSMAP)
variable = list(y1=c(1,1,1), y2=c(2,2,NA), y3=c(3,3,2), y4=c(4,NA,NA), y5=c(5,NA,NA),
                y6=c(6,NA,NA), y7=c(NA,4,NA), y8=c(NA,NA,3))
n_obs <- c(n_obs_ADNI,n_obs_NACC,n_obs_ROSMAP)
w <- n_obs[3]/n_obs # weight for each cohort

G <- crossprod(face_fit_ADNI$fit$y1$Bnew) / nrow(face_fit_ADNI$fit$y1$Bnew)
eig_G <- eigen(G, symmetric = T)
G_half <- eig_G$vectors %*% diag(sqrt(eig_G$values)) %*% t(eig_G$vectors)
G_invhalf <- eig_G$vectors %*% diag(1/sqrt(eig_G$values)) %*% t(eig_G$vectors)


###########################################
## Get estimations by multi-cohort FJM
###########################################

# specify the number of principal components (can be selected by BIC)
npc0 <- 2 # shared
npc1 <- 2 # subject and outcome specific

# get initial value for each cohort through FJM
fit_ADNI_one <- FJM(longdata = longdata_ADNI, time = time_ADNI, event = event_ADNI, 
                    X0 = X0_ADNI, argvals.new = tnew, 
                    npc0 = npc0, npc1 = npc1,
                    face_fit = face_fit_ADNI,
                    nMC0 = 500, nMC1 = 10000, knots = 4,
                    Maxiter = 100, Miniter = 20, 
                    tol0 = 1e-3, tol1 = 5e-3, tol2 = 1e-3, tol3 = 1e-7, trace = F)
    
fit_NACC_one <- FJM(longdata = longdata_NACC_i, time = time_NACC, event = event_NACC, 
                    X0 = X0_NACC, argvals.new = tnew, 
                    npc0 = npc0, npc1 = npc1,
                    face_fit = face_fit_NACC,
                    nMC0 = 500, nMC1 = 10000, knots = 4,
                    Maxiter = 100, Miniter = 20, 
                    tol0 = 1e-3, tol1 = 5e-3, tol2 = 1e-3, tol3 = 1e-7, trace = F)
    
fit_ROSMAP_one <- FJM(longdata = longdata_ROSMAP_i, time = time_ROSMAP, event = event_ROSMAP, 
                      X0 = X0_ROSMAP, argvals.new = tnew, 
                      npc0 = npc0, npc1 = npc1,
                      face_fit = face_fit_ROSMAP,
                      nMC0 = 500, nMC1 = 10000, knots = 4,
                      Maxiter = 100, Miniter = 20, 
                      tol0 = 1e-3, tol1 = 5e-3, tol2 = 1e-3, tol3 = 1e-7, trace = F)
    
ADNI_init <- list(npc0=fit_ADNI_one$npc0, npc1=fit_ADNI_one$npc1,
                      Theta=fit_ADNI_one$Theta, G_half=G_half, G_invhalf=G_invhalf,
                      B=list(fit_ADNI_one$face_fit$fit$y1$fit_mean$B, fit_ADNI_one$face_fit$fit$y2$fit_mean$B,
                             fit_ADNI_one$face_fit$fit$y3$fit_mean$B, fit_ADNI_one$face_fit$fit$y4$fit_mean$B,
                             fit_ADNI_one$face_fit$fit$y5$fit_mean$B, fit_ADNI_one$face_fit$fit$y6$fit_mean$B),
                      scores=fit_ADNI_one$scores, scores_Sig=fit_ADNI_one$scores_Sig)
    
NACC_init <- list(npc0=fit_NACC_one$npc0, npc1=fit_NACC_one$npc1,
                      Theta=fit_NACC_one$Theta, G_half=G_half, G_invhalf=G_invhalf,
                      B=list(fit_NACC_one$face_fit$fit$y1$fit_mean$B, fit_NACC_one$face_fit$fit$y2$fit_mean$B,
                             fit_NACC_one$face_fit$fit$y3$fit_mean$B, fit_NACC_one$face_fit$fit$y4$fit_mean$B,
                             fit_NACC_one$face_fit$fit$y5$fit_mean$B),
                      scores=fit_NACC_one$scores, scores_Sig=fit_NACC_one$scores_Sig)
    
ROSMAP_init <- list(npc0=fit_ROSMAP_one$npc0, npc1=fit_ROSMAP_one$npc1, 
                        Theta=fit_ROSMAP_one$Theta, G_half=G_half, G_invhalf=G_invhalf,
                        B=list(fit_ROSMAP_one$face_fit$fit$y1$fit_mean$B, fit_ROSMAP_one$face_fit$fit$y2$fit_mean$B,
                               fit_ROSMAP_one$face_fit$fit$y3$fit_mean$B, fit_ROSMAP_one$face_fit$fit$y4$fit_mean$B,
                               fit_ROSMAP_one$face_fit$fit$y5$fit_mean$B), 
                        scores=fit_ROSMAP_one$scores, scores_Sig=fit_ROSMAP_one$scores_Sig)
    
    
init <- list(cohort1=ADNI_init, cohort2=NACC_init, cohort3=ROSMAP_init)
    
# multi-cohort FJM    
FJM_fit_multi <- FJM_multi(longdata = longdata, time = time, event = event, 
                           X0 = X0, w=w,variable=variable, argvals.new = tnew, 
                            npc0 = npc0, npc1 = npc1, knots = 4,
                            init = init, nMC0 = 500, nMC1 = 10000,
                            Maxiter = 1000, Miniter = 20, 
                            tol0 = 1e-3, tol1 = 5e-3, tol2 = 1e-3, tol3 = 1e-7)
FJM_fit_multi$BICnew # BIC











