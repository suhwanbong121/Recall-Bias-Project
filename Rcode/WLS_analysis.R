########################################################
### Application of the Proposed Methods 
########################################################
source("basic_functions.R")
library(latex2exp)

## You need to create "wls_anger.csv" first. 
## Use WLS_to_dataset.R to create this csv file. 
dataset <- read.csv("wls_anger.csv")

########################################################
### Data Construction
########################################################

## Outcome variable: anger score
## We use the 75th percentile (10 points) to dichotomize the anger score
dataset$anger75 <- as.integer(dataset$anger >= 10)

## Exposure variable: childhood abuse
## We use a binary abuse exposure variable 
## Abused if there was abuse by either father or mother ("some" or "a lot" response) 
dataset$abuse <- as.integer(dataset$abusefa >= 2 | dataset$abusemo >= 2)

## There are 7 covariates - (to adjust for confounders)
## (1) Sex, (2) age at the time of the interview, (3) father's education, 
## (4) mother's education, (5) parental income, (6) farm background, 
## (7) an indicator of parents' marital problems or single parent. 

## Make (1), (6), (7) binary variables
dataset$female <- as.integer(dataset$sex == 2)
dataset$marital.prob <- as.integer(dataset$problems <= 1)
dataset$rural <- as.integer(dataset$farmback == 1)

## Use a log transformation on (5)
dataset$log.income <- log(dataset$income)

## Now, we finally choose the variables to use.
## 1st column: outcome
## 2nd column: exposure
## 3rd to the last column: covariates
newdata <- dataset[, c("anger75", "abuse", "female", "age", "edufa", "edumo", "marital.prob", "log.income", "rural")]

########################################################
### Data Analysis
########################################################

## We will find the estimates for each eta.
eta0.seq <- seq(0, 0.5, by = 0.01)
eta1.seq <- seq(0, 0.5, by = 0.01)
strata <- 20
bsize <- 50

## Exposure is under-reported
## ML method
ml.mat <- as.matrix(read.csv("ATE_ML.csv"))
ml.mat <- matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    ml.mat[i,j] <- ML.logistic.under(data = newdata, eta = c(eta0.seq[i], eta1.seq[j]))$risk.diff
  }
}

## Propensity score stratification
prop.mat <- as.matrix(read.csv("ATE_Prop.csv"))
prop.mat <-matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
prop.score.model <-glm(abuse ~ female + marital.prob + age + edufa + edumo + log.income + rural, data = newdata, family = "binomial", x = TRUE)
prop.score.star <-prop.score.model$fitted.values
count.mat.prop <-SP.set(newdata, strata + 1, prop.score.star)
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    prop.mat[i,j] <-SP.inference.under(count.mat = count.mat.prop, eta = c(eta0.seq[i], eta1.seq[j]))$risk.diff
  }
}

## Prognostic score stratification
prog.mat <- as.matrix(read.csv("ATE_Prog.csv"))
prog.mat <- matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
newdata_ts <- newdata[newdata$abuse==1,]
prog.score.model <- glm(anger75 ~ abuse + female + marital.prob + age + edufa + edumo + log.income + rural, data = newdata_ts, family = "binomial", x = TRUE)
gamma_x <- prog.score.model$coefficients[-c(1,2)]
prog.score.star <- as.matrix(newdata[,-c(1,2)]) %*% as.vector(gamma_x)
count.mat.prog <- SP.set(newdata, strata + 1, prog.score.star)
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    prog.mat[i,j] <- SP.inference.under(count.mat = count.mat.prog, eta = c(eta0.seq[i], eta1.seq[j]))$risk.diff
  }
}

## Blocking method
block.mat <- as.matrix(read.csv("ATE_Block.csv"))
block.mat <-matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
set.seed(2022)
mat <-Block.set(newdata, bsize = bsize)
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    block.mat[i,j] <- Block.inference.under(mat, eta = c(eta0.seq[i], eta1.seq[j]))$risk.diff
  }
}

########################################################
### Table 3

## Calculating bootstrap confidence intervals
library(doParallel)
eta.seq <- seq(0, 0.5, by = 0.1)
nsim <- 500
strata <- 20
bsize <- 50

cl <- parallel::makeCluster(10) # need to specify the number of cores to use. 
doParallel::registerDoParallel(cl)

## Calculating bootstrap standard deviations
boot.mat <- foreach(k = 1:nsim, .combine = "rbind") %dopar% {
  N <- dim(newdata)[1]
  sample.index <- sample(1:N, N, replace = T)
  boot.data <- newdata[sample.index, ]
  
  ## Bootstrap propensity stratification
  prop.score.model <- glm(anger75 ~ abuse + female + marital.prob + age + edufa + edumo + log.income + rural, data = boot.data, family = "binomial", x = TRUE)
  prop.score.star <- prop.score.model$fitted.values
  boot.prop.mat <- SP.set(boot.data, strata + 1, prop.score.star)
  
  ## Bootstrap prognostic stratification
  newdata_ts <- boot.data[boot.data$abuse==1,]
  prog.score.model <- glm(anger75 ~ abuse + female + marital.prob + age + edufa + edumo + log.income + rural, data = newdata_ts, family = "binomial", x = TRUE)
  gamma_x <- prog.score.model$coefficients[-c(1,2)]
  prog.score.star <- as.matrix(boot.data[,-c(1,2)]) %*% as.vector(gamma_x)
  boot.prog.mat <- SP.set(boot.data, strata + 1, prog.score.star)
  
  ## Bootstrap blocking
  boot.block.mat <- Block.set(boot.data, bsize = bsize)
  
  ## Bootstrap estimators
  boot.ML.est = boot.Prop.est = boot.Prog.est = boot.Block.est = rep(NA, length(eta.seq))
  for(j in 1:length(eta.seq)){
    eta.val <- eta.seq[j]
    boot.ML.est[j] <- ML.logistic.under(data = boot.data, eta = c(eta.val, eta.val))$risk.diff
    boot.Prop.est[j] <- SP.inference.under(count.mat = boot.prop.mat, eta = c(eta.val, eta.val))$risk.diff
    boot.Prog.est[j] <- SP.inference.under(count.mat = boot.prog.mat, eta = c(eta.val, eta.val))$risk.diff
    boot.Block.est[j] <- Block.inference.under(count.mat = boot.block.mat, eta = c(eta.val, eta.val))$risk.diff
  }
  c(boot.ML.est, boot.Prop.est, boot.Prog.est, boot.Block.est)
}
parallel::stopCluster(cl)

boot.mat.ML <- boot.mat[,1:length(eta.seq)]
boot.mat.Prop <- boot.mat[,(length(eta.seq)+1):(2*length(eta.seq))]
boot.mat.Prog <- boot.mat[,(2*length(eta.seq)+1):(3*length(eta.seq))]
boot.mat.Block <- boot.mat[,(3*length(eta.seq)+1):(4*length(eta.seq))]

## Bootstrap standard deviations
ML.sd <- apply((boot.mat.ML), 2, sd)
Prop.sd <- apply((boot.mat.Prop), 2, sd)
Prog.sd <- apply((boot.mat.Prog), 2, sd)
Block.sd <- apply((boot.mat.Block), 2, sd)

## Estimators
ML.est <- diag(ml.mat)[seq(from = 1, by = 10, length = 6)]
Prop.est <- diag(prop.mat)[seq(from = 1, by = 10, length = 6)]
Prog.est <- diag(prog.mat)[seq(from = 1, by = 10, length = 6)]
Block.est <- diag(block.mat)[seq(from = 1, by = 10, length = 6)]
########################################################

library(ggplot2)
library(cowplot)

########################################################
### Figure 1

CType = c("ML" = 1, "Block" = 4)
LType = c("ML" = "solid", "Prop" = "twodash", "Prog" = "dotted", "Block" = "dashed")
LType2 = c("ML" = "solid", "Block" = "dashed")

## Bounds for the ATE
delta.seq <- seq(0, 0.5, by = 0.01)

## Lower bound at eta0 = delta and eta1 = 0
## Upper bound at eta0 = 0 and eta1 = delta
bound.ml.L <- ml.mat[,1]
bound.ml.U <- ml.mat[1,]
bound.prop.L <- prop.mat[,1]
bound.prop.U <- prop.mat[1,]
bound.prog.L <- prog.mat[,1]
bound.prog.U <- prog.mat[1,]
bound.block.L <- block.mat[,1]
bound.block.U <- block.mat[1,]

## Bounds for the ATE with eta0 <= eta1
bound.ml.L2 <- cummin(diag(ml.mat))
bound.prop.L2 <- cummin(diag(prop.mat))
bound.prog.L2 <- cummin(diag(prog.mat))
bound.block.L2 <- cummin(diag(block.mat))

data1 <- data.frame(delta.seq, bound.ml.L, bound.ml.L2, bound.ml.U, ml = diag(ml.mat),
                    bound.prop.L, bound.prop.L2, bound.prop.U, prop = diag(prop.mat),
                    bound.prog.L, bound.prog.L2, bound.prog.U, prog = diag(prog.mat),
                    bound.block.L, bound.block.L2, bound.block.U, block = diag(block.mat))

## Bounds for the ATE - Figure 1(a)
plot1 <- ggplot(data = data1, mapping = aes(x = delta.seq)) +
  labs(title = "(a)", x = TeX(r'($\delta$)'), y = "") +
  geom_line(aes(y = bound.ml.L, linetype = "ML")) +
  geom_line(aes(y = bound.ml.U, linetype = "ML")) +
  geom_line(aes(y = bound.prop.L, linetype = "Prop")) +
  geom_line(aes(y = bound.prop.U, linetype = "Prop")) +
  geom_line(aes(y = bound.prog.L, linetype = "Prog")) +
  geom_line(aes(y = bound.prog.U, linetype = "Prog")) +
  geom_line(aes(y = bound.block.L, linetype = "Block")) +
  geom_line(aes(y = bound.block.U, linetype = "Block")) +
  geom_hline(yintercept = 0, col = "red") +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## Bounds for the ATE - Figure 1(b)
plot2 <- ggplot(data = data1, mapping = aes(x = delta.seq)) +
  labs(title = "(b)", x = TeX(r'($\delta$)'), y = "") +
  geom_line(aes(y = bound.ml.L2, linetype = "ML")) +
  geom_line(aes(y = bound.ml.U, linetype = "ML")) +
  geom_line(aes(y = bound.prop.L2, linetype = "Prop")) +
  geom_line(aes(y = bound.prop.U, linetype = "Prop")) +
  geom_line(aes(y = bound.prog.L2, linetype = "Prog")) +
  geom_line(aes(y = bound.prog.U, linetype = "Prog")) +
  geom_line(aes(y = bound.block.L2, linetype = "Block")) +
  geom_line(aes(y = bound.block.U, linetype = "Block")) +
  geom_hline(yintercept = 0, col = "red") +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## ATE when eta0 = eta1 - Figure 1(c)
plot3 <- ggplot(data = data1, mapping = aes(x = delta.seq)) +
  labs(title = "(c)", x = TeX(r'($\eta_0 = \eta_1$)'), y = "") +
  ylim(0.06,0.12) +
  geom_line(aes(y = ml, linetype = "ML")) +
  geom_line(aes(y = prop, linetype = "Prop")) +
  geom_line(aes(y = prog, linetype = "Prog")) +
  geom_line(aes(y = block, linetype = "Block")) +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## Figure 1(d)
data2 <- data.frame(eta.seq, ML.est, ML.sd, Block.est, Block.sd)
plot4 <- ggplot(data = data2, aes(x = eta.seq)) +
  labs(title = "(d)", x = TeX(r'($\eta_0 = \eta_1$)'), y = "") +
  ylim(-0.02,0.25) +
  geom_line(aes(y = Block.est, linetype = "Block", col = "Block")) +
  geom_point(aes(y = Block.est, col = "Block")) +
  geom_errorbar(width = .02, aes(ymin = Block.est - 1.96 * Block.sd, ymax = Block.est + 1.96 * Block.sd, col = "Block", linetype = "Block")) +
  geom_line(aes(y = ML.est, linetype = "ML", col = "ML")) +
  geom_point(aes(y = ML.est, col = "ML")) +
  geom_errorbar(width = .02, aes(ymin = ML.est - 1.96 * ML.sd, ymax = ML.est + 1.96 * ML.sd, col = "ML", linetype = "ML")) +
  geom_hline(yintercept = 0, col = "red") +
  scale_color_manual(name = "", values = CType) +
  scale_linetype_manual(name = "", values = LType2) +
  theme_test() +
  theme(legend.position = c(0.15, 0.9), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## Figure 1
plot_grid(plot1, plot2, plot3, plot4, ncol = 2, hjust = -1, byrow = TRUE)
########################################################

########################################################
### Figure 2

## Contour plots for the ATE
par(mfrow = c(2, 2), mar = c(3, 3, 2, 2),  mgp=c(2, 1, 0))
contour(eta0.seq, eta1.seq, ml.mat, font.main = 1, main = "(a) Maximum Likelihood", ylab = TeX(r'($\eta_1$)'), nlevels=20, cex.main=1)
contour(eta0.seq, eta1.seq, prop.mat, font.main = 1, main = "(b) Propensity Score Stratification", nlevels=20, cex.main=1)
contour(eta0.seq, eta1.seq, prog.mat, font.main = 1, main = "(c) Prognostic Score Stratification", xlab = TeX(r'($\eta_0$)'), ylab = TeX(r'($\eta_1$)'), nlevels=20, cex.main=1)
contour(eta0.seq, eta1.seq, block.mat, font.main = 1, main = "(d) Blocking", xlab = TeX(r'($\eta_0$)'), nlevels=20, cex.main=1)
########################################################
