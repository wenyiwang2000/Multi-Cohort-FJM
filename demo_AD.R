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

# Load three cohorts AD data
load("./data/ADNI.RData")
load("./data/NACC.RData")
load("./data/ROSMAP.RData")

tnew <- seq(0,1,1/100)
m <- length(tnew)
nm <- 7 # number of basis function

# ADNI
n_ADNI <- nrow(ADNI_subj) 
time_ADNI <- ADNI_subj$obs_time
event_ADNI <- ADNI_subj$event

idx1 <- which(is.na(ADNI$MMSE))
MMSE <- data.frame("subj"=ADNI$ID[-idx1], "argvals"=ADNI$argvals[-idx1], "y"=ADNI$MMSE[-idx1])
idx2 <- which(is.na(ADNI$CDRSB))
CDRSB <- data.frame("subj"=ADNI$ID[-idx2], "argvals"=ADNI$argvals[-idx2], "y"=ADNI$CDRSB[-idx2])
idx3 <- which(is.na(ADNI$WMSLM))
WMSLM <- data.frame("subj"=ADNI$ID[-idx3], "argvals"=ADNI$argvals[-idx3], "y"=ADNI$WMSLM[-idx3])
idx4 <- which(is.na(ADNI$ADAS13))
ADAS13 <- data.frame("subj"=ADNI$ID[-idx4], "argvals"=ADNI$argvals[-idx4], "y"=ADNI$ADAS13[-idx4])
idx5 <- which(is.na(ADNI$RAVLT_imm))
RAVLT_imm <- data.frame("subj"=ADNI$ID[-idx5], "argvals"=ADNI$argvals[-idx5], "y"=ADNI$RAVLT_imm[-idx5])
idx6 <- which(is.na(ADNI$FAQ))
FAQ <- data.frame("subj"=ADNI$ID[-idx6], "argvals"=ADNI$argvals[-idx6], "y"=ADNI$FAQ[-idx6])

X0_ADNI = cbind(ADNI_subj$AGE, ADNI_subj$PTGENDER, ADNI_subj$PTEDUCAT, ADNI_subj$APOE4)
longdata_ADNI <- list("y1"=MMSE, "y2"=CDRSB, "y3"=WMSLM, "y4"=ADAS13, "y5"=RAVLT_imm, "y6"=FAQ) 
n_obs_ADNI <- Reduce("+",lapply(longdata_ADNI,function(x){length(x$y)}))
face_fit_ADNI <- twostep2(longdata = longdata_ADNI, knots = 4, argvals.new = tnew)

# NACC
n_NACC <- nrow(NACC_subj)
time_NACC <- NACC_subj$obs_time
event_NACC <- NACC_subj$event

idx1 <- which(is.na(NACC$NACCMMSE))
MMSE <- data.frame("subj"=NACC$ID[-idx1], "argvals"=NACC$argvals[-idx1], "y"=NACC$NACCMMSE[-idx1])
CDRSB <- data.frame("subj"=NACC$ID, "argvals"=NACC$argvals, "y"=NACC$CDRSUM)
idx3 <- which(is.na(NACC$WMSLM))
WMSLM <- data.frame("subj"=NACC$ID[-idx3], "argvals"=NACC$argvals[-idx3], "y"=NACC$WMSLM[-idx3])
idx4 <- which(is.na(NACC$TRAILA))
TRAILA <- data.frame("subj"=NACC$ID[-idx4], "argvals"=NACC$argvals[-idx4], "y"=NACC$TRAILA[-idx4])
mylambda <- 0.5
TRAILA$y <- ((TRAILA$y)^mylambda-1)/mylambda

X0_NACC <- cbind(NACC_subj$AGE, NACC_subj$SEX, NACC_subj$EDUC, NACC_subj$NACCNE4S)
longdata_NACC <- list("y1"=MMSE, "y2"=CDRSB, "y3"=WMSLM, "y7"=TRAILA)
n_obs_NACC <- Reduce("+",lapply(longdata_NACC,function(x){length(x$y)}))
face_fit_NACC <- twostep2(longdata = longdata_NACC_i, knots = 4, argvals.new = tnew)


# ROSMAP
n_ROSMAP <- nrow(ROSMAP_subj)
time_ROSMAP <- as.numeric(ROSMAP_subj$obs_time)
event_ROSMAP <- ROSMAP_subj$event

idx1 <- which(is.na(ROSMAP$MMSE))
MMSE <- data.frame("subj"=ROSMAP$ID[-idx1], "argvals"=ROSMAP$argvals[-idx1], "y"=ROSMAP$MMSE[-idx1])
idx2 <- which(is.na(ROSMAP$WMSLM))
WMSLM <- data.frame("subj"=ROSMAP$ID[-idx2], "argvals"=ROSMAP$argvals[-idx2], "y"=ROSMAP$WMSLM[-idx2])
idx3 <- which(is.na(ROSMAP$SDMT))
SDMT <- data.frame("subj"=ROSMAP$ID[-idx3], "argvals"=ROSMAP$argvals[-idx3], "y"=ROSMAP$SDMT[-idx3])

X0_ROSMAP <- cbind(ROSMAP_subj$AGE, ROSMAP_subj$SEX, ROSMAP_subj$EDUC, ROSMAP_subj$APOE4)
longdata_ROSMAP <- list("y1"=MMSE, "y3"=WMSLM, "y8"=SDMT)
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
npc0 <- 4 # shared
npc1 <- 3 # subject and outcome specific

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
# BIC
FJM_fit_multi$BICnew 

# Estimated scaling parameters beta_j
round(FJM_fit_multi$Theta$b,4)

# Estimated coefficients gamma in survival model
gamma_est = mapply(function(x){round(x,4)},x=FJM_fit_multi$Theta$beta)
gamma_SE = mapply(function(I){round(sqrt(diag(solve(I))),4)}, I=FJM_fit_multi$I)
pvalues <- round(2*pnorm(-abs(gamma_est / gamma_SE)),4)
res <- cbind(gamma_est[,1], pvalues[,1], gamma_est[,2], pvalues[,2], gamma_est[,3], pvalues[,3])
colnames(res) <- c("ADNI", "pvalue", "NACC", "pvalue", "ROSMAP", "pvalue")
rownames(res) <- c("Age","Gender","Edu","APOE4","Shared1","Shared2","Shared3","Shared4")
res

# Estimated mean functions
Bnew <- fit_ADNI_one$face_fit$fit$y1$Bnew
mu <- lapply(1:8,function(x){
  Bnew %*% FJM_fit_multi$Theta$alpha[[x]]
})

# Estimated eigenvalues
FJM_fit_multi$Theta$D0 
FJM_fit_multi$Theta$D1 

# Estimated eigenfunctions and scores
Bbar <- Bnew %*% fit_ADNI_one$face_fit$G_invhalf
phi <-  as.matrix(Bbar %*% FJM_fit_multi$Theta$u0)/sqrt(96)
psi <-  as.matrix(Bbar %*% FJM_fit_multi$Theta$u1)/sqrt(96)
scores <- FJM_fit_multi$scores



