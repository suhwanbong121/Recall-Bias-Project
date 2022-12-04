########################################################
### Creating A Dataset From the Wisconsin Longitudinal Study
########################################################
## You will need to visit https://www.ssc.wisc.edu/wlsresearch/data/
## to download "wls_bl_14_01.dta" data. 

library(foreign)
library(tidyverse)
library(haven)
library(labelled)
dataset = read_dta("wls_bl_14_01.dta")


sub.dataset = dataset[, c("familypub", "nw036rer", "nw037rer", "nua34rec", "z_sexrsp", "z_ra029re", "edfa57q", "edmo57q", "bmpin1", "nw001rer", "rlur57")]

grad.only.data = sub.dataset[, c("familypub", "edfa57q", "edmo57q", "bmpin1", "rlur57")]
paired.data =sub.dataset[, c("familypub", "nw036rer", "nw037rer", "nua34rec", "z_sexrsp", "z_ra029re", "nw001rer")]


## Select data without missing
grad.only.data.na.rm = remove_missing(grad.only.data)
paired.data.na.rm = remove_missing(paired.data)

combined.data = merge(paired.data.na.rm, grad.only.data.na.rm, by = "familypub")
combined.data.na.rm = remove_missing(combined.data)

new.dataset = subset(combined.data.na.rm, combined.data.na.rm$nw036rer >= 0)
new.dataset = subset(new.dataset, new.dataset$nw037rer >= 0)
new.dataset = subset(new.dataset, new.dataset$nw001rer >= 0)
new.dataset = subset(new.dataset, new.dataset$nua34rec >= 0)
new.dataset = subset(new.dataset, new.dataset$z_sexrsp >= 0)
new.dataset = subset(new.dataset, new.dataset$z_ra029re >= 0)
new.dataset = subset(new.dataset, new.dataset$edfa57q >= 0)
new.dataset = subset(new.dataset, new.dataset$edmo57q >= 0)
new.dataset = subset(new.dataset, new.dataset$bmpin1 >= 0)
new.dataset = subset(new.dataset, new.dataset$rlur57 >= 0)


################
## Save the dataset
colnames(new.dataset) = c("familypub", "abusefa", "abusemo", "anger", "sex", "age", "problems", "edufa", "edumo", "income", "farmback")

write.csv(new.dataset, file = "wls_anger.csv", row.names = F)
