setwd("../18.caoshi_diff/")
require(tidyverse)

rm(list = ls())
set.seed("20240724")
# caoshi cohort ------------------
c2_fs = read.csv("../4.bac_slow_fast/Sepsis_add_trend_cluster.meta.csv",row.names = 1)

c2_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv")
c2_meta = c2_meta[,c(2,499:552)]


all_meta = c2_fs %>% dplyr::select(Sample_ID,type,slope) %>% 
  right_join(c2_meta,by = c("Sample_ID" = "SampleID"))

c2_m = all_meta %>% dplyr::select(PatientID,
                                  Antibiotics_Duration_1w,
                                  Antibiotics_Duration_2w,
                                  Antibiotics_Duration_all,
                                  Age_DOL,
                                  slope,
                                  Onset_Day,
                                  Intrauterine_infection,
                                  Antibiotic_type,
                                  type,
                                  Sepsis_style,
                                  GA_daily ,Natural_delivery ,Feeding ,Gender , Weight)

c2_m = c2_m %>% filter(is.na(Onset_Day) | (Antibiotics_Duration_2w < Onset_Day))

c1_fs = read.csv("../4.bac_slow_fast/Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv")

c1_meta = c1_meta[,c(2,452:542)]

all_meta2 = c1_fs %>% select(Sample_ID,type,slope) %>% inner_join(c1_meta,by = c("Sample_ID" = "SampleID"))

all_meta2 %>% rename(
  "Weight" = "Birth_weight",
  "Gender" = "sex",
  "RDS" = "NRDS",
  "Antibiotics_Duration_all" = "Antibiotic_duration",
  "Intrauterine_infection" = "Infection") %>% 
  mutate(Feeding = "A") %>%
  mutate(Natural_delivery = 
           case_when(
             Delivery == "S" ~ "N",
             Delivery == "E" ~ "Y")) %>%
  mutate(Weight = Weight*1000) %>%
  mutate(Antibiotic_type = NA)-> all_meta

c1_m = all_meta %>% dplyr::select(PatientID,
                                  Antibiotics_Duration_1w,
                                  Antibiotics_Duration_2w,
                                  Antibiotics_Duration_all,
                                  Age_DOL,
                                  slope,
                                  Onset_Day,
                                  Intrauterine_infection,
                                  Antibiotic_type,
                                  type,
                                  Sepsis_style,
                                  GA_daily ,Natural_delivery ,Feeding ,Gender , Weight)

c1_m = c1_m %>% filter(is.na(Onset_Day) | (Antibiotics_Duration_2w < Onset_Day))

#------
c_m = rbind(c2_m,c1_m)

c2_m_p = c_m %>% 
  filter(!is.na(type)) %>%
  group_by(PatientID) %>% 
  slice_head(n = 1) %>% 
  mutate(Feeding = case_when(
  Feeding == "A" ~ "A",
  Feeding != "A" ~ "B"
))

all_meta2 = c2_m_p 
# 构建巢氏队列------
library("MatchIt")

mode_feeding <- as.character(names(sort(table(all_meta2$Feeding), decreasing = TRUE)[1]))
all_meta2$Feeding[is.na(all_meta2$Feeding)] <- mode_feeding

mode_Weight <- mean(all_meta2$Weight, na.rm = TRUE)
all_meta2$Weight[is.na(all_meta2$Weight)] <- mode_Weight

colnames(all_meta2)
f = "type"

# 循环
all_meta3 = all_meta2[!is.na(all_meta2[f]),]

all_meta3$temp <- ifelse(all_meta3[[f]]  == "fast",1,0)

# 缺失值
all_meta3$Intrauterine_infection[is.na(all_meta3$Intrauterine_infection)]  = "Y"
all_meta3 = all_meta3 %>% filter(!is.na(type))

m.out <- matchit(temp ~ GA_daily + 
                   Natural_delivery + 
                   Feeding + 
                   Gender +
                   Weight + 
                   Intrauterine_infection, 
                 data = all_meta3, method = "nearest", 
                 caliper = 0.2,replace = F)
# HAIYOU0.1
re = (summary(m.out))

data.frame(re[["sum.matched"]]) -> pms_data
write.csv(pms_data,'pms_data.csv')

library(cobalt)
love.plot(m.out, 
          stats = c("mean.diffs", "var.ratios"),
          abs = TRUE, 
          thresholds = c(m = 0.2, v = 2),
          var.order = "unadjusted")
ggsave("pms_data.pdf",width = 5,height = 3.5)

m.data <- match.data(m.out)

#----------
data = m.data

data %>% mutate(anti_use = Antibiotic_type) %>% mutate(anti_time = case_when(
  Antibiotics_Duration_2w > 9 ~ "long",
  Antibiotics_Duration_2w <= 9 ~ "short",
)) %>% mutate(anti_sepsis = Sepsis_style) -> data

anti_c = data %>% select(slope,type,
                         PatientID,Antibiotics_Duration_2w,
                         starts_with("anti_")) %>%
  distinct() 

write.csv(anti_c,"caoshi_antibio_data.csv")

anti_c %>% pivot_longer(
  starts_with("anti_"),names_to = "Group",values_to = "Class"
) -> anti_plot

require(ggstatsplot)

grouped_ggpiestats(
  data         = anti_plot,
  x            = Class,
  y            = type,
  grouping.var = Group,
  legend.title = "Class"
)->p1

p1

ggsave("caoshi_antibio_data.pdf",width = 20,height = 6)


saveRDS(p1,"caoshi_antibio_data_p.rds")


