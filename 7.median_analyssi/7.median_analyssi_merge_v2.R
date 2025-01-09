xbox::chdir("./7.median_analyssi")
rm(list = ls())

# 先整理数据
# 合并队列 1 和队列 2 的数据---------------
# 队列2的数据-----------------
set.seed(123)
rm(list = ls())
require(tidyverse)
require(ggstatsplot)

data_s = read.csv("../18.caoshi_diff/caoshi_antibio_data.csv",row.names = 1)
m = data_s

m$anti_time = factor(m$anti_time,levels = c("short","long"))
m$anti_time = as.numeric(m$anti_time)
m$anti_sepsis = factor(m$anti_sepsis,levels  = c("N","LOS"))
m$anti_sepsis = as.numeric(m$anti_sepsis)
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

iv_fct_name = "anti_time"
mv_fct_name = "type"
dv_fct_name = "anti_sepsis"

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

set.seed(2024)
fitMedBoot <- mediate(fitM, fitY, boot=TRUE, sims=999, treat=iv, mediator=mv)
summary(fitMedBoot)-> summ_fit
summ_fit
saveRDS(summ_fit,"mediate_result_2024.rds")
