xbox::chdir("../10.maaslin2_diff")
rm(list  = ls())
#-------
library(tidyverse)
library(mlr3verse)
library(mlr3extralearners)
library(mlr3tuning)
library(data.table)
library(pROC)
library(mlr3viz)
library(future)
library(ggstatsplot)

set.seed(123)

library(phyloseq)
require(microbiomeutilities)
require(ggpubr)
require(microbiome)
require(vegan)

#-------------
getwd()
rm(list = ls())

c1_fs = read.csv("../4.bac_slow_fast/矫正前每个队列随机13_allsample_peakvalue///Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv",row.names = 1)
all_meta = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "ID"))

c1_fs = read.csv("../4.bac_slow_fast/矫正前每个队列随机13_allsample_peakvalue//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv",row.names = 1)
all_meta2 = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "SampleID"))

all_meta %>% rename(
  "Weight" = "Birth_weight",
  "Gender" = "sex",
  "RDS" = "NRDS",
  "Antibiotics_Duration_all" = "Antibiotic_duration",
  "Intrauterine_infection" = "Infection") %>% 
  mutate(Feeding = "C") %>%
  mutate(Natural_delivery = 
           case_when(
             Delivery == "S" ~ "N",
             Delivery == "E" ~ "Y")) %>%
  mutate(Weight = Weight*1000)-> all_meta

nid = intersect(colnames(all_meta),colnames(all_meta2))

all_meta = rbind(all_meta[,nid],all_meta2[,nid])

all_meta = all_meta %>% mutate(
  Feeding = case_when(
    Feeding == "A" ~ "A",
    Feeding != "A" ~ "B"
  )
) 

bac_i = read.csv("../11.bac_index/bac_index_risk.csv",row.names = 1)

all_meta = all_meta %>% filter(is.na(Onset_Day) | Age_DOL < Onset_Day) 


all_meta = all_meta %>% 
  left_join(bac_i %>% select(SampleID,slope,score,risk),by= c("Sample_ID" = "SampleID"))


all_meta = all_meta %>% column_to_rownames("Sample_ID")

all_meta = all_meta[!is.na(all_meta$slope),]

all_meta$type <- ifelse(all_meta$type == "fast", 0, 1)


input_metadata = all_meta[,c(1,303:326)]
input_data = all_meta[,c(1:302,324)]

input_data$slope = (input_data$slope - min(input_data$slope,na.rm = T))/(max(input_data$slope,na.rm = T)-min(input_data$slope,na.rm = T))


library(Maaslin2)
set.seed(20231214)
fit_data_sepsis = Maaslin2(
  input_data, input_metadata,
  min_prevalence = 0.001,
  normalization = "NONE",
  transform = "NONE",
  output = 'maaslin2_diff_mix_Sepsis_style', 
  fixed_effects = c("Sepsis_style","GA_daily", 
                    "Weight", "Feeding", "Gender", 
                    "Age_DOL","Intrauterine_infection"),
  random_effects = c("PatientID"),
  reference = c('Sepsis_style,N')
  )