xbox::chdir("../8.antibio_days_or")
rm(list = ls())

require(tidyverse)
require(ggstatsplot)
c2_fs = read.csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv",row.names = 1)

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
                                  type,
                                  Sepsis_style)

c2_m = c2_m %>% filter(is.na(Onset_Day) | Age_DOL > Onset_Day)


c1_fs = read.csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv")

c1_meta = c1_meta[,c(2,452:542)]

all_meta2 = c1_fs %>% dplyr::select(Sample_ID,type,slope) %>% inner_join(c1_meta,by = c("Sample_ID" = "SampleID"))

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
                                  Age_DOL,
                                  slope,
                                  Onset_Day,
                                  type,
                                  Sepsis_style)

c1_m = c1_m %>% filter(is.na(Onset_Day) | Age_DOL > Onset_Day)

c_m = rbind(c2_m,c1_m)

colpal = readRDS("../colpal.rds")

c2_m_p = c_m %>% group_by(PatientID) %>% slice_head(n = 1) 

write.csv(c2_m_p,"Antibiotics_fix_days.csv")

antib_plot = c2_m_fix2 %>%
  pivot_longer(starts_with("Antibiotics"),
               names_to = "Antibiotics", 
               values_to = "Days")

ggplot(antib_plot,aes(x = Sepsis_style, y = Days))+
  geom_boxplot(width = 0.3,outlier.shape = NA, coef = 0)+
  geom_jitter(aes(color = Sepsis_style),width = 0.2, alpha = 0.8)+
  scale_color_manual(values = colpal)+
  facet_wrap(".~Antibiotics",scales = "free_y")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())->p

require(ggpubr)
p + stat_compare_means(
  comparisons = list( c("LOS", "N")),
  hide.ns = F,
  method = "wilcox.test", label = "p.signif")

ggsave("Sepsis_Antibio_Antibiotics_Duration_diff_fix.pdf",width = 8,height = 4)

write.csv(as.data.frame(table(antib_plot$Antibiotics,antib_plot$Sepsis_style)),
          "Sepsis_Antibio_Antibiotics_Duration_diff_fix.csv")

ggplot(antib_plot %>% filter(!is.na(type)) %>% filter(!grepl("fix",Antibiotics)),aes(x = type, y = Days))+
  geom_boxplot(aes(fill = type),width = 0.3,outlier.shape = NA, coef = 0)+
  scale_fill_manual(values = c("#ff8bbf","#a0a5e1"))+
  facet_wrap(".~Antibiotics",scales = "free_y")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())->p

require(ggpubr)
p + stat_compare_means(
  comparisons = list( c("slow", "fast")),
  hide.ns = F,
  method = "wilcox.test")

ggsave("slowfast_Antibio_Antibiotics_Duration_diff.pdf",width = 6,height = 3)

write.csv(as.data.frame(table(antib_plot$Antibiotics,antib_plot$Sepsis_style)),
          "slowfast_Antibio_Antibiotics_Duration_diff.csv")


data = c_m %>% group_by(PatientID) %>% slice_head(n = 1)


or_sample = data %>% 
  dplyr::select(PatientID,Antibiotics_Duration_2w,type,Sepsis_style) %>% distinct()

or_data1 = as.data.frame(table(or_sample$Antibiotics_Duration_2w,or_sample$type)) %>%
  pivot_wider(names_from = "Var2", values_from = "Freq") %>%
  mutate(RR = slow/(slow+fast))

or_data2 = as.data.frame(table(or_sample$Antibiotics_Duration_2w,or_sample$Sepsis_style)) %>%
  pivot_wider(names_from = "Var2", values_from = "Freq") %>%
  mutate(RR = LOS/(LOS+N))

or_data_all = cbind(or_data1[,c(1,4)],or_data2[,c(4)])
colnames(or_data_all) = c("Day","Type","Sepsis")
write.csv(or_data_all,"total_Sepsis_SLOW_merge_RR.csv")

#
or_data_all = or_data_all %>% hablar::retype() %>%
  pivot_longer(c("Type","Sepsis"),names_to = "Class", values_to = "RR")
or_data_all$RR = round(or_data_all$RR ,2)
require(ggpmisc)
ggplot(data = or_data_all %>% filter(Day >= 6),
       aes(x = Day,y = RR, fill = Class, color = Class))+
  geom_point(size = 3,shape = 23)+
  geom_smooth(aes(x = Day,y = RR, fill = Class, color = Class),
              method = "lm",se = F)+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1)+
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey")+
  geom_text(aes(label = RR),vjust = -0.7)+
  scale_color_manual(values = c("#ff6bda", "#53c6d6")) +
  scale_fill_manual(values = c("#ff6bda", "#53c6d6")) +
  scale_x_continuous(breaks = seq(1,14,1),labels = seq(1,14,1))+
  coord_cartesian(ylim = c(0,1))+
  ggthemes::theme_few()+
  labs(y = "Risk Ratio (RR)" ,x = "Antibiotics Duration (2 Days)")

ggsave("Sepsis_SLOW_merge_RR_1d.pdf",width = 6,height = 4)


ggplot(data = or_data_all %>% filter(Class == "Type") ,aes(x = Day,y = RR, fill = Class, color = Class))+
  geom_point(size = 3,shape = 23)+
  geom_smooth(data = or_data_all %>% filter(Class == "Type"),
              aes(x = Day,y = RR, fill = Class, color = Class),
              method = "lm")+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1)+
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey")+
  geom_text(aes(label = RR),vjust = -0.7)+
  scale_color_manual(values = c("#ff6bda", "#53c6d6")) +
  scale_fill_manual(values = c("#ff6bda", "#53c6d6")) +
  scale_x_continuous(breaks = seq(1,14,1),labels = seq(1,14,1))+
  coord_cartesian(ylim = c(0,1))+
  ggthemes::theme_few()+
  labs(y = "Risk Ratio (RR)" ,x = "Antibiotics Duration (Days)")

ggsave("Sepsis_SLOW_merge_RR_0_14d.pdf",width = 5.5,height = 4)


ggplot(data = or_data_all %>% filter(Class == "Type") ,aes(x = Day,y = RR, fill = Class, color = Class))+
  geom_point(size = 3,shape = 23)+
  geom_smooth(data = or_data_all %>% filter(Class == "Type")%>% filter(Day >= 7),
              aes(x = Day,y = RR, fill = Class, color = Class),
              method = "lm")+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1)+
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey")+
  geom_text(aes(label = RR),vjust = -0.7)+
  scale_color_manual(values = c("#ff6bda", "#53c6d6")) +
  scale_fill_manual(values = c("#ff6bda", "#53c6d6")) +
  scale_x_continuous(breaks = seq(1,14,1),labels = seq(1,14,1))+
  coord_cartesian(ylim = c(0,1))+
  ggthemes::theme_few()+
  labs(y = "Risk Ratio (RR)" ,x = "Antibiotics Duration (Days)")

ggsave("Sepsis_SLOW_merge_RR_7_14d.pdf",width = 5.5,height = 4)


# 换一种展示方式
data = c_m %>% group_by(PatientID) %>% slice_head(n = 1)
or_sample = data %>% 
  dplyr::select(PatientID,Antibiotics_Duration_2w,type,Sepsis_style) %>% distinct()

or_data1 = as.data.frame(table(or_sample$Antibiotics_Duration_2w,or_sample$type)) %>%
  pivot_wider(names_from = "Var2", values_from = "Freq") %>%
  mutate(RR = slow/(slow+fast)) %>%
  mutate(n = slow+fast)%>%
  pivot_longer(-c(Var1,RR,n),names_to = "group", values_to = "num") %>%
  mutate(type = "Type")

or_data2 = as.data.frame(table(or_sample$Antibiotics_Duration_2w,or_sample$Sepsis_style)) %>%
  pivot_wider(names_from = "Var2", values_from = "Freq") %>%
  mutate(RR = LOS/(LOS+N)) %>%
  mutate(n = LOS+N) %>%
  pivot_longer(-c(Var1,RR,n),names_to = "group", values_to = "num") %>%
  mutate(type = "Sepsis")

all_or_data = rbind(or_data1,or_data2)
all_or_data$group = factor(all_or_data$group, levels = c("N","LOS","fast","slow"))
ggplot(all_or_data)+
  geom_bar(aes(x = Var1,y = num,fill = group),position = "fill",stat = "identity")+
  geom_point(aes(x = Var1,y = RR),color = "red")+
  ggalt::geom_xspline(aes(x = Var1,y = RR,group = 1),spline_shape=1,size = 3,color= "white",alpha = 0.8)+
  geom_text(aes(y = 1.1, x = Var1, label = paste0("n=",n)),size = 2.5)+
  geom_text(aes(x = Var1,y = RR,label = round(RR,2)),vjust = 1,size = 4)+
  facet_wrap(".~type")+
  scale_fill_manual(values = c("#53c6d6","#ff6bda","#53c6d6","#ff6bda")) +
  labs(x ="Antibiotics Duration",y = "Patient Ratio")+
  ggthemes::theme_few()

ggsave("Sepsis_SLOW_merge_RR_stack.pdf",width = 8,height = 4)
dev.off()
dev.off()

data = c_m %>% group_by(PatientID) %>% slice_head(n = 1) %>% filter(!is.na(Antibiotics_Duration_2w))
result <- data.frame(matrix(ncol=15, nrow=nrow(data)))
colnames(result) <- c('PatientID', paste0(1:14))

for (i in 1:nrow(data)) {
  patient_id <- data$PatientID[i]
  antibiotics_duration_1w <- data$Antibiotics_Duration_1w[i]
  antibiotics_duration_2w <- data$Antibiotics_Duration_2w[i]

  daily_antibiotics <- rep(0, 14)
  
  if (antibiotics_duration_1w > 0) {
    daily_antibiotics[1:min(7, antibiotics_duration_1w)] <- 1
    daily_antibiotics[8:min(14, antibiotics_duration_2w)] <- 1
  }
  
  daily_antibiotics <- cumsum(daily_antibiotics)
  
  result[i, ] <- c(patient_id, daily_antibiotics)
}

antibio_day = left_join(data %>% dplyr::select(PatientID,type,Sepsis_style),result,by = "PatientID")
write.csv(antibio_day,"Antibiotics_Duration_patient_sum.csv")


antibio_day_plot = antibio_day %>% 
  pivot_longer(-c(PatientID,type,Sepsis_style),names_to = "Days",values_to = "Antibiotics_Duration(Days)")

antibio_day_plot = antibio_day_plot %>% hablar::retype()

antibio_day_plot_diff = antibio_day_plot %>% 
  group_by(Days) %>%
  rstatix::wilcox_test(`Antibiotics_Duration(Days)`~Sepsis_style) %>%
  rstatix::add_significance(p.col = "p")

write.csv(antibio_day_plot_diff,"Antibiotics_Duration_diff_between_Sepsis.csv")

ggplot()+
  geom_smooth(data = antibio_day_plot, aes(x = Days, y = `Antibiotics_Duration(Days)`, fill = Sepsis_style,color = Sepsis_style))+
  geom_text(data = antibio_day_plot_diff,aes(x= Days,y = 12,label = p),color = "red")+
  geom_vline(xintercept = 9, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 7, linetype = "dashed", color = "grey")+
  scale_color_manual(values = c("#ff6bda", "#53c6d6")) +
  scale_fill_manual(values = c("#ff6bda", "#53c6d6")) +
  scale_x_continuous(breaks = seq(1,14,1),labels = seq(1,14,1))+
  scale_y_continuous(breaks = c(7,14),labels = c(7,14))+
  coord_cartesian(ylim = c(0,14))+
  labs(title = "Daily Accumulated Antibiotic Courses for Patients")+
  ggthemes::theme_few()+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("Accumulated_Antibiotics_Duration_diff_between_Sepsis.pdf",
       width = 6,height = 4)

antibio_day_plot_diff = antibio_day_plot %>% 
  group_by(Days) %>%
  rstatix::wilcox_test(`Antibiotics_Duration(Days)`~type) %>%
  rstatix::add_significance(p.col = "p")

write.csv(antibio_day_plot_diff,"Antibiotics_Duration_diff_between_type.csv")

ggplot()+
  geom_smooth(data = antibio_day_plot %>% filter(!is.na(type)), 
              aes(x = Days, y = `Antibiotics_Duration(Days)`, 
                  fill = type,color = type))+
  geom_text(data = antibio_day_plot_diff,aes(x= Days,y = 12,label = p),color = "red")+
  geom_vline(xintercept = 9, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 7, linetype = "dashed", color = "grey")+
  scale_color_manual(values = c("#ff6bda", "#53c6d6")) +
  scale_fill_manual(values = c("#ff6bda", "#53c6d6")) +
  scale_x_continuous(breaks = seq(1,14,1),labels = seq(1,14,1))+
  scale_y_continuous(breaks = c(7,14),labels = c(7,14))+
  coord_cartesian(ylim = c(0,14))+
  labs(title = "Daily Accumulated Antibiotic Courses for Patients")+
  ggthemes::theme_few()+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("Accumulated_Antibiotics_Duration_diff_between_type.pdf",
       width = 6,height = 4)