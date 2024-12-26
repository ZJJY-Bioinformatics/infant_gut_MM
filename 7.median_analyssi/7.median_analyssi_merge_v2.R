xbox::chdir("../7.median_analyssi")
rm(list = ls())
source("../addin.r")

set.seed(123)
rm(list = ls())
require(tidyverse)

c2_m = read.csv("../8.antibio_days_or/Antibiotics_fix_days.csv",row.names = 1) %>% distinct()

# 开始分析
m = c2_m
m$Sepsis_style = factor(m$Sepsis_style,levels = c("N","LOS"))
m$Sepsis_style = as.numeric(m$Sepsis_style)
m$type = factor(m$type,levels = c("fast","slow"))
m$type = as.numeric(m$type)


#中介分析主体
require(tidyverse)
library(mediation)
library(tidyverse)
library(lmerTest)
library(broom.mixed)

# 从这里开始循环-------
medi_result = data.frame()

iv_fct_name = "Antibiotics_Duration_2w_fix"
mv_fct_name = "slope"
dv_fct_name = "Sepsis_style"

cor_iv_mv = cor.test(m[[iv_fct_name]],m[[mv_fct_name]],method = "spearman")[["estimate"]][["rho"]]
p_iv_mv = cor.test(m[[iv_fct_name]],m[[mv_fct_name]],method = "spearman")[["p.value"]]
cor_mv_dv = cor.test(m[[mv_fct_name]],m[[dv_fct_name]],method = "spearman")[["estimate"]][["rho"]]
p_mv_dv = cor.test(m[[mv_fct_name]],m[[dv_fct_name]],method = "spearman")[["p.value"]]

# 计算正向中介分析-----  
mData <- m[,c(iv_fct_name,mv_fct_name,dv_fct_name)] 

colnames(mData)[1]=iv_fct_name
colnames(mData)[2]=mv_fct_name
colnames(mData)[3]=dv_fct_name
colnames(mData) = gsub(" ","_",colnames(mData))
colnames(mData) = gsub("^\\d-","",colnames(mData))
colnames(mData) = make.names(colnames(mData))

id = colnames(mData)

iv = id[1]
mv = id[2]
dv = id[3]

set.seed("20240131")

fitM = lm(as.formula(paste(mv,"~",iv)), data = mData)
fitY <- lm(as.formula(paste(dv,"~",iv,"+",mv)), data=mData)


fitMedBoot <- mediate(fitM, fitY, boot=TRUE, sims=999, treat=iv, mediator=mv)
summary(fitMedBoot)-> summ_fit
summ_fit

fit1 = lm(as.formula(paste(mv,"~",dv)), data = mData)
fit2 = lm(as.formula(paste(iv,"~",mv)), data = mData)
fit3 = lm(as.formula(paste(iv,"~",dv)), data = mData)

i = 1
medi_result[i,1] = summ_fit[["d.avg.p"]]
medi_result[i,2] = summ_fit[["d.avg"]] 
medi_result[i,3] = summ_fit[["d.avg.ci"]][1]
medi_result[i,4] = summ_fit[["d.avg.ci"]][2]
medi_result[i,5] = summ_fit[["z.avg.p"]]
medi_result[i,6] = summ_fit[["z.avg"]]
medi_result[i,7] = summ_fit[["z.avg.ci"]][1]
medi_result[i,8] = summ_fit[["z.avg.ci"]][2]
medi_result[i,9] = iv
medi_result[i,10] = mv
medi_result[i,11] = dv
medi_result[i,12] = cor_iv_mv
medi_result[i,13] = p_iv_mv
medi_result[i,14] = cor_mv_dv
medi_result[i,15] = p_mv_dv

colnames(medi_result) = c("ACME_p",
                          "ACME_beta",
                          "ACME_beta_CI_lower",
                          "ACME_beta_CI_upper",
                          "ADE_p","ADE_beta",
                          "ADE_beta_CI_lower",
                          "ADE_beta_CI_upper",
                          "iv",
                          "mv",
                          "dv",
                          "cor_iv_mv",
                          "p_iv_mv",
                          "cor_mv_dv",
                          "p_mv_dv")

medi_result_cis = medi_result

write_csv(medi_result_cis,paste0("medi_result_trans_",iv,"_-",mv,"_-",dv,".csv"))