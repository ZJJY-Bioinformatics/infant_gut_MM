xbox::chdir("./3.bac_exp")
rm(list = ls())
library(tidyverse)
library(lmerTest)
library(broom.mixed)

rm(list = ls())
c1_fs = read.csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv",row.names = 1)
all_meta = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "ID"))

c1_fs = read.csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv",row.names = 1)
all_meta2 = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "SampleID"))

all_meta %>% rename(
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
  mutate(Weight = Weight*1000)-> all_meta

nid = intersect(colnames(all_meta),colnames(all_meta2))

ndata = read.csv("../1.data_tidy/all_cohort.sum_7cohort_final.csv",row.names = 1)
nid = c(nid[1:2],colnames(ndata)[2:48],nid[304:324])

all_meta = rbind(all_meta[,nid],all_meta2[,nid])
data = all_meta %>% filter(Age_DOL <=42)

write.csv(data,"merge_c1_c2_bacexp_meta.csv")

fit = lmer(Klebsiella ~ type * Age_DOL + GA_daily + Weight + Gender + Feeding + (1|PatientID), data = data)
res = summary(fit)
res = tidy(fit)

label = paste0("Type: P=", round(res$p.value[2], 3), "\n Age_DOL: P=", round(res$p.value[3], 3), "\nType * Age_DOL: P=", round(res$p.value[7], 3))
colpal = readRDS("../colpal.rds")
ggplot(data, aes(Age_DOL, Klebsiella, color = type, fill = type)) + 
  geom_smooth(alpha = 0.2,size = 2) + 
  annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
  ggthemes::theme_few()+
  scale_color_manual(values = colpal) + 
  scale_fill_manual(values = colpal) + 
  xlab("Days")

lmm_table = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ type * Age_DOL + GA_daily + Weight + Gender + Feeding +(1|PatientID), data = data2)
  res = tidy(fit) %>% 
    filter(term %in% c("typeslow", "Age_DOL", "typeslow:Age_DOL")) %>% 
    mutate(y = {{taxa}}) %>% 
    select(y, everything())
  return(res)
}

#taxa_id = colnames(data)[3:303][apply(data[,3:303],2,sum) != 0]
taxa_id = colnames(data)[3:49][apply(data[,3:49],2,sum) != 0]

all = map_dfr(taxa_id, ~lmm_table(.x))

all %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr")) -> all_fdr
write_tsv(all_fdr, "LMM-fastslow + Age_DOL_adjust_confounder_result_FDR.tsv")

lmm_func = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ type * Age_DOL + GA_daily + Weight + Gender + Feeding + (1|PatientID), data = data2)
  res = tidy(fit)
  
  label = paste0("Type: P=", round(res$p.value[2], 3), "\n Age_DOL: P=", round(res$p.value[3], 3), "\nType * Age_DOL: P=", round(res$p.value[9], 3))
  
  p =   ggplot(data2, aes(Age_DOL, feature, color = type, fill = type)) + 
    geom_smooth(alpha = 0.2) + 
    annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
    ggthemes::theme_few()+
    scale_color_manual(values = c("#1B8574","#bf48a8")) + 
    scale_fill_manual(values = c("#1B8574","#bf48a8")) + 
    theme(axis.text = element_text(color = "black")) +
    labs(x = "Day", y = taxa)
  
  return(p)
}


Ps = map(taxa_id , ~lmm_func(.x))
pdf("LMM-fastslow + Age_DOL_adjust_confounder_result_FDR.pdf", width = 5, height = 3)
walk(Ps, print)
dev.off()


taxt_data = as.data.frame(apply(data[,3:49],2,mean)) %>% rownames_to_column("bac")
colnames(taxt_data)[2] = "exp"
colnames(all_fdr)[1] = "bac"
all_fdr %>% select(bac,term,p.value,FDR) %>%
  left_join(taxt_data, by = "bac") -> plot_fdr

plot_fdr %>%
  filter(term == "typeslow") %>%
  mutate(sign = case_when(
  FDR <= 0.05 ~ "Significant",
  FDR > 0.05 ~ "Insignificant",
)) %>%
  mutate(type = "fastslow") -> plot_fdr1

ggplot()+
  geom_point(data = plot_fdr1,aes(x = exp,size = exp, color =  sign, y = -log(FDR)))+
  geom_text(data = plot_fdr1 %>% filter(exp > 0.01),aes(x = exp,y = -log(FDR),label = bac))+
  geom_hline(yintercept = -log(0.05),linetype = "dashed",color = "#bcbcbc")+
  scale_color_manual(values = c("#bcbcbc","#cc0000"))+
  labs(x = "Average relative abundance of bacteria")+
  ggthemes::theme_few()
ggsave("LMM-fastslow_fdr_signif.pdf",width = 6,height = 4)
write.csv(plot_fdr1,"LMM-fastslow_fdr_signif.csv")

# sepsis之间----------
rm(list = ls())
c1_fs = read.csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv",row.names = 1)
all_meta = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "ID"))

c1_fs = read.csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv",row.names = 1)
all_meta2 = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "SampleID"))

all_meta %>% rename(
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
  mutate(Weight = Weight*1000)-> all_meta

nid = intersect(colnames(all_meta),colnames(all_meta2))

all_meta = rbind(all_meta[,nid],all_meta2[,nid])
data = all_meta %>% filter(Age_DOL <=42)


fit = lmer(Klebsiella ~ Sepsis_style * Age_DOL + GA_daily + Weight + Gender +Feeding + (1|PatientID), data = data)
res = summary(fit)
res = tidy(fit)

label = paste0("Sepsis_style: P=", round(res$p.value[2], 3), "\n Age_DOL: P=", round(res$p.value[3], 3), "\nSepsis_style * Age_DOL: P=", round(res$p.value[7], 3))

ggplot(data, aes(Age_DOL, Klebsiella, color = Sepsis_style, fill = Sepsis_style)) + 
  geom_smooth(alpha = 0.2,size = 2) + 
  annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
  theme_minimal() + 
  scale_color_manual(values = c("#bf48a8","#1B8574")) + 
  scale_fill_manual(values = c("#bf48a8","#1B8574")) + 
  xlab("Days")

lmm_table = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ Sepsis_style * Age_DOL + GA_daily + Weight + Gender +Feeding+ (1|PatientID), data = data2)
  res = tidy(fit) %>% 
    filter(term %in% c("Sepsis_styleN", "Age_DOL", "Sepsis_styleN:Age_DOL")) %>% 
    mutate(y = {{taxa}}) %>% 
    select(y, everything())
  return(res)
}
taxa_id = colnames(data)[3:49][apply(data[,3:49],2,sum) != 0]
all = map_dfr(taxa_id, ~lmm_table(.x))

all %>% 
  filter(grepl("Sepsis",term)) %>%
  mutate(FDR = p.adjust(p.value, method = "fdr")) -> all_fdr
write_tsv(all_fdr, "LMM-Sepsis_style + Age_DOL_adjust_confounder_result_FDR.xls")


lmm_func = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ Sepsis_style * Age_DOL + GA_daily + Weight + Gender+ Feeding + (1|PatientID), data = data2)
  res = tidy(fit)
  
  label = paste0("Sepsis_style: P=", round(res$p.value[2], 3), "\n Age_DOL: P=", round(res$p.value[3], 3), "\nSepsis_style * Age_DOL: P=", round(res$p.value[9], 3))
  
  p =   ggplot(data2, aes(Age_DOL, feature, color = Sepsis_style, fill = Sepsis_style)) + 
    geom_smooth(alpha = 0.2) + 
    annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
    ggthemes::theme_few()+
    scale_color_manual(values = c("#bf48a8","#1B8574")) + 
    scale_fill_manual(values = c("#bf48a8","#1B8574")) + 
    labs(x = "Day", y = taxa)
  
  return(p)
}

Ps = map(taxa_id, ~lmm_func(.x))
pdf("LMM-Sepsis_style_interaction_scatterplot_AgeDOL_Sepsis_first.pdf", width = 5, height = 3)
walk(Ps, print)
dev.off() 


taxt_data = as.data.frame(apply(data[,3:49],2,mean)) %>% rownames_to_column("bac")

colnames(taxt_data)[2] = "exp"
colnames(all_fdr)[1] = "bac"
all_fdr %>% select(bac,term,p.value,FDR) %>%
  left_join(taxt_data, by = "bac") -> plot_fdr

plot_fdr %>%
  filter(term == "Sepsis_styleN:Age_DOL") %>%
  mutate(sign = case_when(
    FDR <= 0.05 ~ "Significant",
    FDR > 0.05 ~ "Insignificant",
  )) %>% mutate(type = "Sepsis_style")-> plot_fdr2



ggplot()+
  geom_point(data = plot_fdr2,aes(x = exp,size = exp, color =  sign, y = -log(FDR)))+
  geom_text(data = plot_fdr2 %>% filter(exp > 0.01),aes(x = exp,y = -log(FDR),label = bac))+
  geom_hline(yintercept = -log(0.05),linetype = "dashed",color = "#bcbcbc")+
  scale_color_manual(values = c("#bcbcbc","#cc0000"))+
  labs(x = "Average relative abundance of bacteria")+
  ggthemes::theme_few()
ggsave("LMM-sepsis_fdr_signif.pdf",width = 6,height = 4)
write.csv(plot_fdr,"LMM-sepsis_fdr_signif.csv")


lmm_func = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  
  fit = lmer(feature ~ Sepsis_style * Age_DOL + GA_daily + Weight + Gender + Feeding + (1|PatientID), data = data2)
  res = tidy(fit)
  
  label = paste0("Sepsis_style: P=", round(res$p.value[2], 3), 
                 "\nAge_DOL: P=", round(res$p.value[3], 3), 
                 "\nSepsis_style * Age_DOL: P=", round(res$p.value[9], 3))
  
  data2 = data2 %>% mutate(taxa = taxa)
  

  return(list(data = data2, label = label))
}


plot_data_list = map(taxa_id, lmm_func)


plot_data = map_dfr(plot_data_list, "data")

labels = map_chr(plot_data_list, "label")
labels_df = data.frame(taxa = taxa_id, label = labels)
plot_data$taxa = factor(plot_data$taxa,levels = rev(names(sort(apply(data[,3:49],2,sum)))))
labels_df$taxa = factor(labels_df$taxa,levels = rev(names(sort(apply(data[,3:49],2,sum)))))

p = ggplot(plot_data, aes(Age_DOL, feature, color = Sepsis_style, fill = Sepsis_style)) + 
  geom_smooth(alpha = 0.2) + 
  facet_wrap(~ taxa, ncol = 5,scales = "free_y") +  # 按taxa进行分面
  ggthemes::theme_few() +
  scale_color_manual(values = c("#bf48a8","#1B8574")) + 
  scale_fill_manual(values = c("#bf48a8","#1B8574")) + 
  xlab("Days") +
  geom_text(data = labels_df, aes(x = -Inf, y = Inf, label = label), hjust = -0.1, vjust = 1.1, inherit.aes = FALSE)  # 显示每个taxa的注释

ggsave("LMM-Sepsis_style_interaction_facetplot_AgeDOL_Sepsis.pdf", plot = p, width = 20, height = 25)


xbox::chdir("3.bac_exp")
rm(list = ls())
library(tidyverse)
library(lmerTest)
library(broom.mixed)

c1_fs = read.csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv",row.names = 1)
all_meta = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "ID"))

c1_fs = read.csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv",row.names = 1)
all_meta2 = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "SampleID"))

all_meta %>% rename(
  "Weight" = "Birth_weight",
  "Gender" = "sex"
) -> all_meta

nid = intersect(colnames(all_meta),colnames(all_meta2))

ndata = read.csv("../1.data_tidy/all_cohort.sum_7cohort_final.csv",row.names = 1)

nid = c(nid[1:2],colnames(ndata)[2:48],nid[304:316])


all_meta = rbind(all_meta[,nid],all_meta2[,nid])
data = all_meta %>% filter(Age_DOL <=42)

# bac_index
bac_index = read.csv("../11.bac_index/bac_index_risk.csv",row.names = 1)

data  = data %>% left_join(bac_index %>% select(SampleID,score,risk), by = c("Sample_ID" = "SampleID"))

fit = lmer(Klebsiella ~ risk * Age_DOL + GA_daily + Weight + Gender + (1|PatientID), data = data)
res = summary(fit)
res = tidy(fit)

label = paste0("Risk: P=", round(res$p.value[2], 3), "\n Age_DOL: P=", round(res$p.value[3], 3), "\nRisk * Age_DOL: P=", round(res$p.value[7], 3))
table(data$risk)
ggplot(data, aes(Age_DOL, Klebsiella, color = risk, fill = risk)) + 
  geom_smooth(alpha = 0.2,size = 2) + 
  annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
  theme_minimal() + 
  scale_color_manual(values = c("#bf48a8","#1B8574")) + 
  scale_fill_manual(values = c("#bf48a8","#1B8574")) + 
  xlab("Days")

lmm_table = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ risk * Age_DOL + GA_daily + Weight + Gender + (1|PatientID), data = data2)
  res = tidy(fit) %>% 
    filter(term %in% c("risklow", "Age_DOL", "risklow:Age_DOL")) %>% 
    mutate(y = {{taxa}}) %>% 
    select(y, everything())
  return(res)
}
taxa_id = colnames(data)[3:49][apply(data[,3:49],2,sum) != 0]
all = map_dfr(taxa_id, ~lmm_table(.x))

all %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr")) -> all_fdr
write_tsv(all_fdr, "LMM-risk + Age_DOL_adjust_confounder_result_FDR.tsv")


lmm_func = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ risk * Age_DOL + GA_daily + Weight + Gender + (1|PatientID), data = data2)
  res = tidy(fit)
  
  label = paste0("Risk: P=", round(res$p.value[2], 3), "\n Age_DOL: P=", round(res$p.value[3], 3), "\nRisk * Age_DOL: P=", round(res$p.value[7], 3))
  
  p =   ggplot(data2, aes(Age_DOL, feature, color = risk, fill = risk)) + 
    geom_smooth(alpha = 0.2) + 
    annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
    theme_minimal() + 
    scale_color_manual(values = c("#bf48a8","#1B8574")) + 
    scale_fill_manual(values = c("#bf48a8","#1B8574")) + 
    theme(axis.text = element_text(color = "black")) +
    labs(x = "Day", y = taxa)
  
  return(p)
}

Ps = map(taxa_id, ~lmm_func(.x))
pdf("LMM-risk_style_interaction_scatterplot_AgeDOL_Sepsis.pdf", width = 6, height = 4)
walk(Ps, print)
dev.off()

taxt_data = as.data.frame(apply(data[,3:49],2,mean)) %>% rownames_to_column("bac")

colnames(taxt_data)[2] = "exp"
colnames(all_fdr)[1] = "bac"
all_fdr %>% select(bac,term,p.value,FDR) %>%
  left_join(taxt_data, by = "bac") -> plot_fdr

plot_fdr %>%
  filter(term == "risklow:Age_DOL") %>%
  mutate(sign = case_when(
    FDR <= 0.05 ~ "Significant",
    FDR > 0.05 ~ "Insignificant",
  )) %>% mutate(type = "Sepsis_style")-> plot_fdr2

ggplot()+
  geom_point(data = plot_fdr2,aes(x = exp,size = exp, color =  sign, y = -log(FDR)))+
  geom_text(data = plot_fdr2 %>% filter(exp > 0.01),aes(x = exp,y = -log(FDR),label = bac))+
  geom_hline(yintercept = -log(0.05),linetype = "dashed",color = "#bcbcbc")+
  scale_color_manual(values = c("#bcbcbc","#cc0000"))+
  labs(x = "Average relative abundance of bacteria")+
  ggthemes::theme_few()
ggsave("LMM-risk_fdr_signif.pdf",width = 6,height = 4)
write.csv(plot_fdr,"LMM-risk_fdr_signif.csv")