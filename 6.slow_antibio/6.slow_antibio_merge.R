xbox::chdir("6.slow_antibio")
rm(list = ls())

require(tidyverse)
c2_fs = read.csv("../4.bac_slow_fast/Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c2_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv")
c2_meta = c2_meta[,c(2,499:552)]


all_meta = c2_fs %>% dplyr::select(Sample_ID,type,slope) %>% right_join(c2_meta,by = c("Sample_ID" = "SampleID"))

c2_m = all_meta %>% dplyr::select(PatientID,
                                  Antibiotics_Duration_1w,
                                  Antibiotics_Duration_2w,
                                  Antibiotics_Duration_all,
                                  slope,
                                  Onset_Day,
                                  type,
                                  Sepsis_style)

c2_m = c2_m %>% filter(is.na(Onset_Day) | (Onset_Day >= Antibiotics_Duration_2w))

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
  mutate(Weight = Weight*1000)-> all_meta

c1_m = all_meta %>% dplyr::select(PatientID,
                                  Antibiotics_Duration_1w,
                                  Antibiotics_Duration_2w,
                                  Antibiotics_Duration_all,
                                  slope,
                                  Onset_Day,
                                  type,
                                  Sepsis_style)

c1_m = c1_m %>% filter(is.na(Onset_Day) | (Onset_Day >= Antibiotics_Duration_2w))

c2_m = rbind(c2_m,c1_m)

#---------
antib_plot = c2_m %>% group_by(PatientID) %>% slice_head(n = 1) %>%
  pivot_longer(starts_with("Antibiotics"),
               names_to = "Antibiotics", 
               values_to = "Days")

ggplot(antib_plot,aes(x = Sepsis_style, y = Days))+
  geom_boxplot(width = 0.3,outlier.shape = NA, coef = 0)+
  geom_jitter(aes(color = Sepsis_style),width = 0.2, alpha = 0.8)+
  scale_color_manual(values = c("#ff6bda", "#53c6d6"))+
  facet_wrap(".~Antibiotics",scales = "free_y")+
  ggthemes::theme_few() -> p
require(ggpubr)
p + stat_compare_means(
  comparisons = list( c("LOS", "N")),
  hide.ns = F,
  method = "wilcox.test", label = "p.signif")

ggsave("Sepsis_Antibio_Antibiotics_Duration_diff.pdf",width = 8,height = 4)

write.csv(as.data.frame(table(antib_plot$Antibiotics,antib_plot$Sepsis_style)),
          "Sepsis_Antibio_Antibiotics_Duration_diff.csv")


ggplot(antib_plot %>% filter(!is.na(type)),aes(x = type, y = Days))+
  geom_boxplot(width = 0.3,outlier.shape = NA, coef = 0)+
  geom_jitter(aes(color = type),width = 0.2, alpha = 0.8)+
  scale_color_manual(values = c("#ff6bda", "#53c6d6"))+
  facet_wrap(".~Antibiotics",scales = "free_y")+
  ggthemes::theme_few() -> p

require(ggpubr)
p + stat_compare_means(
  comparisons = list( c("slow", "fast")),
  hide.ns = F,
  method = "wilcox.test")

ggsave("slowfast_Antibio_Antibiotics_Duration_diff.pdf",width = 8,height = 4)

write.csv(as.data.frame(table(antib_plot$Antibiotics,antib_plot$Sepsis_style)),
          "slowfast_Antibio_Antibiotics_Duration_diff.csv")


