########################################################
### Codes for Web Appendix
########################################################

# Use the same dataset
source("basic_functions.R")
library(latex2exp)
library(tidyverse)
library(cowplot)
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
### Table 1-3 - Baseline covariate balance adjustment

## Table 1
block.mat = Block.set(newdata, bsize = 20)
table1 <- Check.covbal(newdata, block.mat, eta = c(0, 0))
table2 <- Check.covbal(newdata, block.mat, eta = c(0, 0.5))
table3 <- Check.covbal(newdata, block.mat, eta = c(0.5, 0.5))
########################################################


## We will find the estimates -- ATT, RR, OR -- for each eta (0 <= eta <= 0.5)
eta0.seq <- seq(0, 0.5, by = 0.1)
eta1.seq <- seq(0, 0.5, by = 0.1)

## Set the number of strata
strata <- 20

## Set the block size
bsize <- 50

## Exposure is under-reported
## ML method
ml.rr <- matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
ml.or <- matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    ml.rr[i,j] <- ML.logistic.under(data = newdata, eta = c(eta0.seq[i], eta1.seq[j]))$est.RR
    ml.or[i,j] <- ML.logistic.under(data = newdata, eta = c(eta0.seq[i], eta1.seq[j]))$est.OR
  }
}

## Propensity score stratification
prop.att <-matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
prop.rr <-matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
prop.or <-matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
prop.score.model <-glm(abuse ~ female + marital.prob + age + edufa + edumo + log.income + rural, data = newdata, family = "binomial", x = TRUE)
prop.score.star <-prop.score.model$fitted.values
count.mat.prop <-SP.set(newdata, strata + 1, prop.score.star)
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    prop.att[i,j] <-SP.inference.under(count.mat = count.mat.prop, eta = c(eta0.seq[i], eta1.seq[j]))$ATT
    prop.rr[i,j] <-SP.inference.under(count.mat = count.mat.prop, eta = c(eta0.seq[i], eta1.seq[j]))$est.RR
    prop.or[i,j] <-SP.inference.under(count.mat = count.mat.prop, eta = c(eta0.seq[i], eta1.seq[j]))$est.OR
  }
}

## Prognostic score stratification
prog.att <- matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
prog.rr <- matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
prog.or <- matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
newdata_ts <- newdata[newdata$abuse==1,]
prog.score.model <- glm(anger75 ~ abuse + female + marital.prob + age + edufa + edumo + log.income + rural, data = newdata_ts, family = "binomial", x = TRUE)
gamma_x <- prog.score.model$coefficients[-c(1,2)]
prog.score.star <- as.matrix(newdata[,-c(1,2)]) %*% as.vector(gamma_x)
count.mat.prog <- SP.set(newdata, strata + 1, prog.score.star)
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    prog.att[i,j] <- SP.inference.under(count.mat = count.mat.prog, eta = c(eta0.seq[i], eta1.seq[j]))$ATT
    prog.rr[i,j] <- SP.inference.under(count.mat = count.mat.prog, eta = c(eta0.seq[i], eta1.seq[j]))$est.RR
    prog.or[i,j] <- SP.inference.under(count.mat = count.mat.prog, eta = c(eta0.seq[i], eta1.seq[j]))$est.OR
  }
}

## Blocking method
block.att <-matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
block.rr <-matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
block.or <-matrix(NA, nrow = length(eta0.seq), ncol = length(eta1.seq))
set.seed(2022)
mat <-Block.set(newdata, bsize = bsize)
for(i in 1:length(eta0.seq)){
  for(j in 1:length(eta1.seq)){
    block.att[i,j] <- Block.inference.under(mat, eta = c(eta0.seq[i], eta1.seq[j]))$ATT
    block.rr[i,j] <- Block.inference.under(mat, eta = c(eta0.seq[i], eta1.seq[j]))$est.RR
    block.or[i,j] <- Block.inference.under(mat, eta = c(eta0.seq[i], eta1.seq[j]))$est.OR
  }
}

########################################################
### Figure 1 -- ATT
LType = c("Prop" = "twodash", "Prog" = "dotted", "Block" = "dashed")

## Bounds for the ATT
delta.seq <- seq(0, 0.5, by = 0.1)

## Lower bound at eta0 = delta and eta1 = 0
## Upper bound at eta0 = 0 and eta1 = delta
bound.prop.L.ATT <- prop.att[,1]
bound.prop.U.ATT <- prop.att[1,]
bound.prog.L.ATT <- prog.att[,1]
bound.prog.U.ATT <- prog.att[1,]
bound.block.L.ATT <- block.att[,1]
bound.block.U.ATT <- block.att[1,]

## Bounds for the ATT with eta0 <= eta1
bound.prop.L2.ATT <- cummin(diag(prop.att))
bound.prog.L2.ATT <- cummin(diag(prog.att))
bound.block.L2.ATT <- cummin(diag(block.att))

data2 <- data.frame(delta.seq, bound.prop.L.ATT, bound.prop.L2.ATT, bound.prop.U.ATT, prop = diag(prop.att),
                    bound.prog.L.ATT, bound.prog.L2.ATT, bound.prog.U.ATT, prog = diag(prog.att),
                    bound.block.L.ATT, bound.block.L2.ATT, bound.block.U.ATT, block = diag(block.att))

## Bounds for the ATT - Figure 1(a)
plot1 <- ggplot(data = data2, mapping = aes(x = delta.seq)) +
  labs(title = "(a)", x = TeX(r'($\delta$)'), y = "") +
  geom_line(aes(y = bound.prop.L.ATT, linetype = "Prop")) +
  geom_line(aes(y = bound.prop.U.ATT, linetype = "Prop")) +
  geom_line(aes(y = bound.prog.L.ATT, linetype = "Prog")) +
  geom_line(aes(y = bound.prog.U.ATT, linetype = "Prog")) +
  geom_line(aes(y = bound.block.L.ATT, linetype = "Block")) +
  geom_line(aes(y = bound.block.U.ATT, linetype = "Block")) +
  geom_hline(yintercept = 0, col = "red") +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## Bounds for the ATT - Figure 1(b)
plot2 <- ggplot(data = data2, mapping = aes(x = delta.seq)) +
  labs(title = "(b)", x = TeX(r'($\delta$)'), y = "") +
  geom_line(aes(y = bound.prop.L2.ATT, linetype = "Prop")) +
  geom_line(aes(y = bound.prop.U.ATT, linetype = "Prop")) +
  geom_line(aes(y = bound.prog.L2.ATT, linetype = "Prog")) +
  geom_line(aes(y = bound.prog.U.ATT, linetype = "Prog")) +
  geom_line(aes(y = bound.block.L2.ATT, linetype = "Block")) +
  geom_line(aes(y = bound.block.U.ATT, linetype = "Block")) +
  geom_hline(yintercept = 0, col = "red") +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## ATT when eta0 = eta1 - Figure 1(c)
plot3 <- ggplot(data = data2, mapping = aes(x = delta.seq)) +
  labs(title = "(c)", x = TeX(r'($\eta_0 = \eta_1$)'), y = "") +
  geom_line(aes(y = prop, linetype = "Prop")) +
  geom_line(aes(y = prog, linetype = "Prog")) +
  geom_line(aes(y = block, linetype = "Block")) +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## Figure 1
plot_grid(plot1, plot2, plot3,  ncol = 2, hjust = -1, byrow = TRUE)
########################################################

########################################################
### Figure 2 - Risk Ratio
LType = c("ML" = "solid", "Prop" = "twodash", "Prog" = "dotted", "Block" = "dashed")

## Bounds for the RR
delta.seq <- seq(0, 0.5, by = 0.1)

## Lower bound at eta0 = delta and eta1 = 0
## Upper bound at eta0 = 0 and eta1 = delta
bound.ml.L.RR <- ml.rr[,1]
bound.ml.U.RR <- ml.rr[1,]
bound.prop.L.RR <- prop.rr[,1]
bound.prop.U.RR <- prop.rr[1,]
bound.prog.L.RR <- prog.rr[,1]
bound.prog.U.RR <- prog.rr[1,]
bound.block.L.RR <- block.rr[,1]
bound.block.U.RR <- block.rr[1,]

## Bounds for the RR with eta0 <= eta1
bound.ml.L2.RR <- cummin(diag(ml.rr))
bound.prop.L2.RR <- cummin(diag(prop.rr))
bound.prog.L2.RR <- cummin(diag(prog.rr))
bound.block.L2.RR <- cummin(diag(block.rr))

data3 <- data.frame(delta.seq, bound.ml.L.RR, bound.ml.L2.RR, bound.ml.U.RR, ml = diag(ml.rr),
                    bound.prop.L.RR, bound.prop.L2.RR, bound.prop.U.RR, prop = diag(prop.rr),
                    bound.prog.L.RR, bound.prog.L2.RR, bound.prog.U.RR, prog = diag(prog.rr),
                    bound.block.L.RR, bound.block.L2.RR, bound.block.U.RR, block = diag(block.rr))

## Bounds for the RR - Figure 2(a)
plot1 <- ggplot(data = data3, mapping = aes(x = delta.seq)) +
  labs(title = "(a)", x = TeX(r'($\delta$)'), y = "") +
  geom_line(aes(y = bound.ml.L.RR, linetype = "ML")) +
  geom_line(aes(y = bound.ml.U.RR, linetype = "ML")) +
  geom_line(aes(y = bound.prop.L.RR, linetype = "Prop")) +
  geom_line(aes(y = bound.prop.U.RR, linetype = "Prop")) +
  geom_line(aes(y = bound.prog.L.RR, linetype = "Prog")) +
  geom_line(aes(y = bound.prog.U.RR, linetype = "Prog")) +
  geom_line(aes(y = bound.block.L.RR, linetype = "Block")) +
  geom_line(aes(y = bound.block.U.RR, linetype = "Block")) +
  geom_hline(yintercept = 1, col = "red") +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## Bounds for the RR - Figure 2(b)
plot2 <- ggplot(data = data3, mapping = aes(x = delta.seq)) +
  labs(title = "(b)", x = TeX(r'($\delta$)'), y = "") +
  geom_line(aes(y = bound.ml.L2.RR, linetype = "ML")) +
  geom_line(aes(y = bound.ml.U.RR, linetype = "ML")) +
  geom_line(aes(y = bound.prop.L2.RR, linetype = "Prop")) +
  geom_line(aes(y = bound.prop.U.RR, linetype = "Prop")) +
  geom_line(aes(y = bound.prog.L2.RR, linetype = "Prog")) +
  geom_line(aes(y = bound.prog.U.RR, linetype = "Prog")) +
  geom_line(aes(y = bound.block.L2.RR, linetype = "Block")) +
  geom_line(aes(y = bound.block.U.RR, linetype = "Block")) +
  geom_hline(yintercept = 1, col = "red") +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## RR when eta0 = eta1 - Figure 2(c)
plot3 <- ggplot(data = data3, mapping = aes(x = delta.seq)) +
  labs(title = "(c)", x = TeX(r'($\eta_0 = \eta_1$)'), y = "") +
  geom_line(aes(y = ml, linetype = "ML")) +
  geom_line(aes(y = prop, linetype = "Prop")) +
  geom_line(aes(y = prog, linetype = "Prog")) +
  geom_line(aes(y = block, linetype = "Block")) +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## Figure 2
plot_grid(plot1, plot2, plot3,  ncol = 2, hjust = -1, byrow = TRUE)
########################################################

########################################################
### Figure 3 - Odds Ratio
LType = c("ML" = "solid", "Prop" = "twodash", "Prog" = "dotted", "Block" = "dashed")

## Bounds for the OR
delta.seq <- seq(0, 0.5, by = 0.1)

## Lower bound at eta0 = delta and eta1 = 0
## Upper bound at eta0 = 0 and eta1 = delta
bound.ml.L.OR <- ml.or[,1]
bound.ml.U.OR <- ml.or[1,]
bound.prop.L.OR <- prop.or[,1]
bound.prop.U.OR <- prop.or[1,]
bound.prog.L.OR <- prog.or[,1]
bound.prog.U.OR <- prog.or[1,]
bound.block.L.OR <- block.or[,1]
bound.block.U.OR <- block.or[1,]

## Bounds for the RR with eta0 <= eta1
bound.ml.L2.OR <- cummin(diag(ml.or))
bound.prop.L2.OR <- cummin(diag(prop.or))
bound.prog.L2.OR <- cummin(diag(prog.or))
bound.block.L2.OR <- cummin(diag(block.or))

data4 <- data.frame(delta.seq, bound.ml.L.OR, bound.ml.L2.OR, bound.ml.U.OR, ml = diag(ml.or),
                    bound.prop.L.OR, bound.prop.L2.OR, bound.prop.U.OR, prop = diag(prop.or),
                    bound.prog.L.OR, bound.prog.L2.OR, bound.prog.U.OR, prog = diag(prog.or),
                    bound.block.L.OR, bound.block.L2.OR, bound.block.U.OR, block = diag(block.or))

## Bounds for the OR - Figure 3(a)
plot1 <- ggplot(data = data4, mapping = aes(x = delta.seq)) +
  labs(title = "(a)", x = TeX(r'($\delta$)'), y = "") +
  geom_line(aes(y = bound.ml.L.OR, linetype = "ML")) +
  geom_line(aes(y = bound.ml.U.OR, linetype = "ML")) +
  geom_line(aes(y = bound.prop.L.OR, linetype = "Prop")) +
  geom_line(aes(y = bound.prop.U.OR, linetype = "Prop")) +
  geom_line(aes(y = bound.prog.L.OR, linetype = "Prog")) +
  geom_line(aes(y = bound.prog.U.OR, linetype = "Prog")) +
  geom_line(aes(y = bound.block.L.OR, linetype = "Block")) +
  geom_line(aes(y = bound.block.U.OR, linetype = "Block")) +
  geom_hline(yintercept = 1, col = "red") +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## Bounds for the OR - Figure 3(b)
plot2 <- ggplot(data = data4, mapping = aes(x = delta.seq)) +
  labs(title = "(b)", x = TeX(r'($\delta$)'), y = "") +
  geom_line(aes(y = bound.ml.L2.OR, linetype = "ML")) +
  geom_line(aes(y = bound.ml.U.OR, linetype = "ML")) +
  geom_line(aes(y = bound.prop.L2.OR, linetype = "Prop")) +
  geom_line(aes(y = bound.prop.U.OR, linetype = "Prop")) +
  geom_line(aes(y = bound.prog.L2.OR, linetype = "Prog")) +
  geom_line(aes(y = bound.prog.U.OR, linetype = "Prog")) +
  geom_line(aes(y = bound.block.L2.OR, linetype = "Block")) +
  geom_line(aes(y = bound.block.U.OR, linetype = "Block")) +
  geom_hline(yintercept = 1, col = "red") +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## OR when eta0 = eta1 - Figure 3(c)
plot3 <- ggplot(data = data4, mapping = aes(x = delta.seq)) +
  labs(title = "(c)", x = TeX(r'($\eta_0 = \eta_1$)'), y = "") +
  geom_line(aes(y = ml, linetype = "ML")) +
  geom_line(aes(y = prop, linetype = "Prop")) +
  geom_line(aes(y = prog, linetype = "Prog")) +
  geom_line(aes(y = block, linetype = "Block")) +
  scale_linetype_manual(name = "", values = LType) +
  theme_test() +
  theme(legend.position = c(0.15, 0.85), legend.title = element_blank(), legend.text = element_text(size=7), legend.key.size = unit(0.4, 'cm'), legend.key.width = unit(0.9, "cm"), plot.title = element_text(hjust = 0.5, size = 8)) 

## Figure 3
plot_grid(plot1, plot2, plot3, ncol = 2, hjust = -1, byrow = TRUE)
########################################################