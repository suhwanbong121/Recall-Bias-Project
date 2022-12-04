########################################################
### Application of the Proposed Methods 
########################################################
source("basic_functions.R")
library(latex2exp)

## You need to create "wls_anger.csv" first. 
## Use WLS_to_dataset.R to create this csv file. 
dataset = read.csv("wls_anger.csv")

########################################################
### Data Construction
########################################################

## Outcome variable: anger score
## We use the 75th percentile (10 points) to dichotomize the anger score
dataset$anger75 = as.integer(dataset$anger >= 10)

## Exposure variable: childhood abuse
## We use a binary abuse exposure variable 
## Abused if there was abuse by either father or mother ("some" or "a lot" response) 
dataset$abuse = as.integer(dataset$abusefa >= 2 | dataset$abusemo >= 2)

## There are 7 covariates - (to adjust for confounders)
## (1) Sex, (2) age at the time of the interview, (3) father's education, 
## (4) mother's education, (5) parental income, (6) farm background, 
## (7) an indicator of parents' marital problems or single parent. 

## Make (1), (6), (7) binary variables
dataset$female = as.integer(dataset$sex == 2)
dataset$marital.prob = as.integer(dataset$problems <= 1)
dataset$rural = as.integer(dataset$farmback == 1)

## Use a log transformation on (5)
dataset$log.income = log(dataset$income)

## Now, we finally choose the variables to use.
## 1st column: outcome
## 2nd column: exposure
## 3rd to the last column: covariates
newdata = dataset[, c("anger75", "abuse", "female", "age", "edufa", "edumo", "marital.prob", "log.income", "rural")]

########################################################
### Data Analysis
########################################################

## We will find the estimates for each eta.
eta0.seq = seq(0, 0.5, by = 0.01)
eta1.seq = seq(0, 0.5, by = 0.01)
strata = 20
bsize = 50

## Exposure is under-reported
## ML method
ml.mat = matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    res = ML.logistic.under(data = newdata, eta = c(eta0.seq[i], eta1.seq[j]))
    ml.mat[i,j] = res$risk.diff
  }
}

## Propensity score stratification
prop.mat = matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
prop.score.model = glm(abuse ~ female + marital.prob + age + edufa + edumo + log.income + rural, data = newdata, family = "binomial", x = TRUE)
prop.score.star = prop.score.model$fitted.values
count.mat.prop = SP.set(newdata, strata + 1, prop.score.star)
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    prop.mat[i,j] = SP.inference.under(count.mat = count.mat.prop, eta = c(eta0.seq[i], eta1.seq[j]))$risk.diff
  }
}

## Prognostic score stratification
prog.mat = matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
newdata_ts = newdata[newdata$abuse==1,]
prog.score.model = glm(anger75 ~ abuse + female + marital.prob + age + edufa + edumo + log.income + rural, data = newdata_ts, family = "binomial", x = TRUE)
gamma_x = prog.score.model$coefficients[-c(1,2)]
prog.score.star = as.matrix(newdata[,-c(1,2)]) %*% as.vector(gamma_x)
count.mat.prog = SP.set(newdata, strata + 1, prog.score.star)
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    prog.mat[i,j] = SP.inference.under(count.mat = count.mat.prog, eta = c(eta0.seq[i], eta1.seq[j]))$risk.diff
  }
}

## Blocking method
block.mat = matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
set.seed(0)
block.mat = Block.set(newdata, bsize = bsize)
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    block.mat[i,j] = Block.inference.under(block.mat = block.mat, eta = c(eta0.seq[i], eta1.seq[j]))$risk.diff
  }
}

########################################################
### Figure 1

## Bounds for the ATE - Stratification
## Define a function to calculate a bound
bound = function(count.mat, delta){
  a.star = count.mat[,1] # exposed & case
  b.star = count.mat[,2] # exposed & control
  c.star = count.mat[,3] # unexposed & case
  d.star = count.mat[,4] # unexposed & control
  n.vec = a.star + b.star + c.star + d.star
  
  p1U = a.star/(a.star + b.star*(1-delta))
  p1L = a.star/(a.star + b.star/(1-delta))
  p0L = (c.star-delta*a.star/(1-delta))/(c.star + d.star - delta*a.star/(1-delta))
  p0U = (c.star)/(c.star + d.star - delta*b.star/(1-delta))
  
  pU = p1U-p0L
  pL = p1L-p0U
  
  UB = sum(n.vec*pU)/sum(n.vec)
  LB = sum(n.vec*pL)/sum(n.vec)
  
  return (c(LB,UB))
}

## Select delta sequence and calculate the bounds
delta.seq = seq(0,0.5, by = 0.01)
l = length(delta.seq)
bound.prop.U = bound.prog.U = bound.block.U = rep(NA, length(delta.seq))
bound.prop.L = bound.prog.L = bound.block.L = rep(NA, length(delta.seq))
for(i in 1:length(delta.seq)){
  bound.prop.L[i] = bound(count.mat.prop, delta.seq[i])[1]
  bound.prop.U[i] = bound(count.mat.prop, delta.seq[i])[2]
  bound.prog.L[i] = bound(count.mat.prog, delta.seq[i])[1]
  bound.prog.U[i] = bound(count.mat.prog, delta.seq[i])[2]
  bound.block.L[i] = bound(count.mat.block, delta.seq[i])[1]
  bound.block.U[i] = bound(count.mat.block, delta.seq[i])[2]
}

## Bounds for the ATE - ML
## Lower bound at eta0 = delta and eta1 = 0
## Upper bound at eta0 = 0 and eta1 = delta
bound.ml.L = ml.mat[,1]
bound.ml.U = ml.mat[1,]

## Plot
jpeg("Bound.jpeg",width=6, height = 4.5, units='in',res=1200)
plot(delta.seq,bound.ml.L,type="l",col=1, ylim = c(-0.15, 0.35), xlab = TeX(r'($\delta$)'), ylab = "Average Treatment Effect (ATE)")
lines(delta.seq,bound.ml.U,type="l",col=1)
lines(delta.seq,bound.prop.L,type="l",lty = 2)
lines(delta.seq,bound.prop.U,type="l",lty = 2)
lines(delta.seq,bound.prog.L,type="l",lty = 3)
lines(delta.seq,bound.prog.U,type="l",lty = 3)
lines(delta.seq,bound.block.L,type="l",lty = 4)
lines(delta.seq,bound.block.U,type="l",lty = 4)
lines(delta.seq,rep(0,length(delta.seq)),lwd = 2, lty =1, col="red")
legend("topright",c("ML","Prop","Prog","Block"), lwd = 1, lty = c(1,2,3,4),cex = 0.5)
dev.off()
########################################################

########################################################
### Figure 2

## Point estimates of the ATE
jpeg("ATE.jpeg",width=6, height = 4.5, units='in',res=1200)
plot(eta0.seq, diag(ml.mat), type = "l", ylim = c(0.05,0.12), xlab = TeX(r'($\eta_0 = \eta_1$)'),ylab = "Average Treatment Effect (ATE)", lty = 1)
lines(eta0.seq, diag(prop.mat), type = "l", lty = 2)
lines(eta0.seq, diag(prog.mat), type = "l", lty = 3)
lines(eta0.seq, diag(block.mat), type = "l", lty = 4)
legend("topright",c("ML","Prop","Prog","Block"),lwd = 1, lty=c(1,2,3,4),cex = 0.5)
dev.off()
########################################################


########################################################
### Table 3

## Calculating bootstrap confidence intervals
library(doParallel)
eta.seq = seq(0, 0.5, by = 0.1)
nsim = 500
strata = 20
bsize = 50

cl <- parallel::makeCluster(10) # need to specify the number of cores to use. 
doParallel::registerDoParallel(cl)

## Calculating bootstrap standard deviations
boot.mat = foreach(k = 1:nsim, .combine = "rbind") %dopar% {
  N = dim(newdata)[1]
  sample.index = sample(1:N, N, replace = T)
  boot.data = newdata[sample.index, ]
  
  ## Bootstrap propensity stratification
  prop.score.model = glm(anger75 ~ abuse + female + marital.prob + age + edufa + edumo + log.income + rural, data = boot.data, family = "binomial", x = TRUE)
  prop.score.star = prop.score.model$fitted.values
  boot.prop.mat = SP.set(boot.data, strata + 1, prop.score.star)
  
  ## Bootstrap prognostic stratification
  newdata_ts = boot.data[boot.data$abuse==1,]
  prog.score.model = glm(anger75 ~ abuse + female + marital.prob + age + edufa + edumo + log.income + rural, data = newdata_ts, family = "binomial", x = TRUE)
  gamma_x = prog.score.model$coefficients[-c(1,2)]
  prog.score.star = as.matrix(boot.data[,-c(1,2)]) %*% as.vector(gamma_x)
  boot.prog.mat = SP.set(boot.data, strata + 1, prog.score.star)
  
  ## Bootstrap blocking
  boot.block.mat = Block.set(boot.data, bsize = bsize)
  
  ## Bootstrap estimators
  boot.ML.est = boot.Prop.est = boot.Prog.est = boot.Block.est = rep(NA, length(eta.seq))
  for(j in 1:length(eta.seq)){
    eta.val = eta.seq[j]
    boot.ML.est[j] = ML.logistic.under(data = boot.data, eta = c(eta.val, eta.val))$risk.diff
    boot.Prop.est[j] = SP.inference.under(count.mat = boot.prop.mat, eta = c(eta.val, eta.val))$risk.diff
    boot.Prog.est[j] = SP.inference.under(count.mat = boot.prog.mat, eta = c(eta.val, eta.val))$risk.diff
    boot.Block.est[j] = Block.inference.under(count.mat = boot.block.mat, eta = c(eta.val, eta.val))$risk.diff
  }
  c(boot.ML.est, boot.Prop.est, boot.Prog.est, boot.Block.est)
}
parallel::stopCluster(cl)

boot.mat.ML = boot.mat[,1:length(eta.seq)]
boot.mat.Prop = boot.mat[,(length(eta.seq)+1):(2*length(eta.seq))]
boot.mat.Prog = boot.mat[,(2*length(eta.seq)+1):(3*length(eta.seq))]
boot.mat.Block = boot.mat[,(3*length(eta.seq)+1):(4*length(eta.seq))]

## Bootstrap standard deviations
ML.sd = apply((boot.mat.ML), 2, sd)
Prop.sd = apply((boot.mat.Prop), 2, sd)
Prog.sd = apply((boot.mat.Prog), 2, sd)
Block.sd = apply((boot.mat.Block), 2, sd)

## Estimators
ML.est = diag(ml.mat)[seq(from = 1, by = 10, length = 6)]
Prop.est = diag(prop.mat)[seq(from = 1, by = 10, length = 6)]
Prog.est = diag(prog.mat)[seq(from = 1, by = 10, length = 6)]
Block.est = diag(block.mat)[seq(from = 1, by = 10, length = 6)]
########################################################


########################################################
### Figure 3

## Contour plots for the ATE
jpeg("ATE_ML.jpeg",width=6, height = 4.5, units='in',res=1200)
contour(eta0.seq, eta1.seq, ml.mat, xlab = TeX(r'($\eta_0$)'), ylab = TeX(r'($\eta_1$)'), nlevels=20)
dev.off()
jpeg("ATE_Prop.jpeg",width=6, height = 4.5, units='in',res=1200)
contour(eta0.seq, eta1.seq, prop.mat, xlab = TeX(r'($\eta_0$)'), ylab = TeX(r'($\eta_1$)'), nlevels=20)
dev.off()
jpeg("ATE_Prog.jpeg",width=6, height = 4.5, units='in',res=1200)
contour(eta0.seq, eta1.seq, prog.mat, xlab = TeX(r'($\eta_0$)'), ylab = TeX(r'($\eta_1$)'), nlevels=20)
dev.off()
jpeg("ATE_Block.jpeg",width=6, height = 4.5, units='in',res=1200)
contour(eta0.seq, eta1.seq, block.mat, xlab = TeX(r'($\eta_0$)'), ylab = TeX(r'($\eta_1$)'), nlevels=20)
dev.off()
########################################################



