xbox::chdir("./22.robert_vaild")
rm(list = ls())
#wwowowowoowowoowowowowoowowoowowowoowowoowoowowoowowowo

require(tidyverse)
library(mlr3verse)
library(mlr3extralearners)
library(mlr3tuning)
library(data.table)
library(pROC)
library(mlr3viz)
library(future)
merge_sum = read.csv("../1.data_tidy/all_cohort_sample_num_8cohort_0905.csv",row.names = 1)
table(merge_sum$cohort)
merge_sum = merge_sum %>% filter(Age_DOL <= 42)

# dir.create("bac_profile")

top10 = readRDS("../19.pseudotime/top10.rds")

meta_name = colnames(merge_sum)[49:ncol(merge_sum)]
top10_data = merge_sum[,c("Sample_ID",top10,meta_name)]

top10_data$Other = 1-apply(top10_data[2:11],1, sum) 
top10_data$Other[top10_data$Other < 0] = 0

colnames(top10_data)[ncol(top10_data)] = "Others"

dodge_pdata = top10_data %>% 
  pivot_longer(-c(Sample_ID,Age_DOL,cohort,
                  PatientID,Sepsis_style,NEC), values_to = "Abundance", names_to = "Bac")

top10_data %>%  
  arrange(Age_DOL) -> temp1

dodge_pdata$Sample_ID = factor(dodge_pdata$Sample_ID ,
                               levels = temp1$Sample_ID)

cohorts <- c("CALM2005_1_A", 
             "CALM2005_2_A", 
             "Rao_2021_A",
             "Lauren_2022_M",
             "Olm_2019_M",
             "Brooks_2017_M",
             "Raveh_2015_M",
             "Robert_2024_M")

dodge_pdata$Bac = factor(dodge_pdata$Bac, levels = c(top10,"Others"))
dodge_pdata$cohort = factor(dodge_pdata$cohort,levels = cohorts)

palette =  readRDS("../../colpal.rds")


dodge_pdata_p <- dodge_pdata %>%
  group_by(Age_DOL, Bac) %>%
  summarise(Abundance = mean(Abundance))

ggplot(dodge_pdata_p %>% filter(Age_DOL <= 42 & Age_DOL != 0) ,aes(x = Age_DOL,y  = Abundance, fill =  Bac))+
  geom_bar(stat = "identity", position = "fill")+
  scale_fill_manual(values = palette) +
  labs(y = "relative abundance", x = "Time (Days)")+
  coord_cartesian(xlim = c(0,42)) +
  scale_x_continuous(breaks = c(seq(1,42,3),42),
                     labels = c(seq(1,42,3),42))+
  ggthemes::theme_clean()+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("bac_profile/8cohort_sum_dodge.pdf",width = 8,height = 4)


dodge_pdata_p <- dodge_pdata %>%
  group_by(Age_DOL, cohort,Bac) %>%
  summarise(Abundance = mean(Abundance))

ggplot(dodge_pdata_p %>% filter(Age_DOL <= 42 & Age_DOL != 0) ,aes(x = Age_DOL,y  = Abundance, fill =  Bac))+
  geom_bar(stat = "identity", position = "fill")+
  scale_fill_manual(values = palette) +
  labs(y = "relative abundance", x = "Time (Days)")+
  coord_cartesian(xlim = c(0,42)) +
  scale_x_continuous(breaks = c(seq(1,42,3),42),
                     labels = c(seq(1,42,3),42))+
  scale_y_continuous(breaks = c(0,1),
                     labels = c(0,1))+
  facet_wrap(".~cohort",ncol = 1,strip.position = "right")+
  ggthemes::theme_few()
ggsave("bac_profile/8cohort_sum_dodge_split.png",width = 7,height = 12)
ggsave("bac_profile/8cohort_sum_dodge_split.pdf",width = 7,height = 12)


source("../addin.r")
require(vegan)
count_data = round(merge_sum[,2:48]*10000)
rownames(count_data) = merge_sum$Sample_ID
alpha <- alpha_diversity(count_data)
cbind(merge_sum[,49:53],alpha) -> plot_a_data

mod_top10_data = top10_data
mod_top10_data$Others = as.numeric(plot_a_data$Simpson)
dodge_pdata = mod_top10_data %>% 
  pivot_longer(-c(Sample_ID,Age_DOL,cohort,
                  PatientID,Sepsis_style,NEC), values_to = "Abundance", names_to = "Bac")
top10_data %>%  
  arrange(Age_DOL) -> temp1

dodge_pdata$Sample_ID = factor(dodge_pdata$Sample_ID ,levels = temp1$Sample_ID)

dodge_pdata$Bac = factor(dodge_pdata$Bac, levels = c(top10,"Others"))
dodge_pdata$cohort = factor(dodge_pdata$cohort,levels = cohorts)

ggplot(dodge_pdata %>% filter(Age_DOL <= 42) %>% filter(Bac != "Others"), aes(x = Age_DOL,y  = Abundance,color =  Bac,fill = Bac))+
  geom_smooth(method = "loess",span = 0.9)+
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(y = "Relative abundance", x = "Age DOL")+
  facet_wrap("cohort~Bac",ncol = 10,scales= "free_y")+
  coord_cartesian(xlim = c(0,43))+
  ggthemes::theme_few()+
  theme(axis.text = element_text(size = 7))+
  theme(strip.text = element_blank())

ggsave("bac_profile/8cohort_sum_smooth_facet.png",width = 13,height = 5.5)
ggsave("bac_profile/8cohort_sum_smooth_facet.pdf",width = 13,height = 5.5)


dodge_p = dodge_pdata %>% 
  filter(Age_DOL <= 42) %>% 
  filter(Bac != "Others") %>%
  select(Sample_ID,PatientID,Age_DOL,cohort,Bac,Abundance)

library(dplyr)
library(tidyr)
library(psych)

dodge_p %>% 
  select(Sample_ID,Age_DOL,Bac,cohort,Abundance) %>%
  pivot_wider(id_cols = c(Bac,Sample_ID,Age_DOL),names_from = cohort,values_from = Abundance) %>%
  select(-Sample_ID) %>%
  group_by(Age_DOL,Bac) %>%
  mutate(across(.cols = everything(), .fns = ~sum(., na.rm = TRUE))) %>%
  group_by(Age_DOL) %>%
  slice_head(n = 10) -> re

re %>%
  group_by(Bac) %>%  
  do({
    icc_data <- column_to_rownames(.,"Age_DOL")  
    icc_value <- ICC(icc_data[,-1])$results  
    tibble(
      Bac = unique(.$Bac),
      ICC_value = icc_value$ICC[6],  
      p_value = icc_value$p[6]  
    )
  }) -> p_values

write.csv(p_values,"top10_bac_smoothplot_ICC_p.csv")


ggplot(dodge_pdata %>% filter(Age_DOL <= 42) %>% filter(Bac != "Others"), aes(x = Age_DOL,y  = Abundance,color =  Bac,fill = Bac))+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(y = "Relative abundance", x = "Age DOL")+
  facet_wrap("cohort~Bac",ncol = 10,scales= "free_y")+
  coord_cartesian(xlim = c(0,43))+
  ggthemes::theme_few()+
  theme(axis.text = element_text(size = 7))+
  theme(strip.text = element_blank())

ggsave("bac_profile/8cohort_sum_line_facet.png",width = 14,height = 5)
ggsave("bac_profile/8cohort_sum_line_facet.pdf",width = 14,height = 5)

learner = readRDS("../4.bac_slow_fast/model_pridict_learner.rds")

all_data = merge_sum %>% filter(cohort == "Robert_2024_M")
feature = all_data[,-c(49,51,52,53)] %>% filter(Age_DOL <=  42 )

c_task = as_task_regr(feature, target = "Age_DOL", id = "regression")

c_task$set_col_roles("Sample_ID", roles = "name")

c_prediction = learner$predict(c_task)

c_prediction$score(msr("regr.mse"))
c_prediction$score(msr("regr.rmse"))
c_prediction$score(msr("regr.srho"))
c_prediction$score(msr("regr.rsq"))
c_prediction$score(msr("regr.mae"))
data.frame(regr.mse = c_prediction$score(msr("regr.mse")),
           regr.rmse = c_prediction$score(msr("regr.rmse")),
           regr.srho = c_prediction$score(msr("regr.srho")),
           regr.mae = c_prediction$score(msr("regr.mae")),
           regr.rsq = c_prediction$score(msr("regr.rsq"))) -> regr.result
write.csv(regr.result,"model_vaild/regr.result_test.csv")

best_tuned_pred <- data.frame(c_prediction$print())
age_prd = best_tuned_pred

best_tuned_pred$PatientID = all_data$PatientID
best_tuned_pred$Sample_ID = all_data$Sample_ID
best_tuned_pred$Age_DOL = all_data$Age_DOL
best_tuned_pred$Sepsis_style = all_data$Sepsis_style
best_tuned_pred$NEC = all_data$NEC
best_tuned_pred$cohort = all_data$cohort

colpal = readRDS("../colpal.rds")
require(ggpmisc)
ggplot(best_tuned_pred,aes(x = truth, y = response))+
  geom_jitter(fill = "#E0B794",alpha = 0.5,size = 3,color = "white",shape = 21)+
  geom_smooth(method = 'lm',formula = y~x,color = "#ff8bbf",se = T)+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1)+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Actual day of life", y = "Predicted day of life")

ggsave("model_vaild/mix_cohort_predict_model.pdf", width = 3.5,height = 3)

pre_age_pdata = best_tuned_pred
colnames(pre_age_pdata)

data  = best_tuned_pred
patients <- unique(data$PatientID)


lm_results <- data.frame(patient = character(),
                         slope = numeric(),
                         intercept = numeric(),
                         stringsAsFactors = FALSE)

for (p in patients) {

  data_subset <- data[data$PatientID == p, ]

  linear_model <- lm(response ~ truth, data = data_subset)
  
  slope <- coef(linear_model)[2]
  intercept <- coef(linear_model)[1]

  lm_results <- rbind(lm_results, data.frame(patient = p, slope = slope, intercept = intercept))
}


print(lm_results)

slop_data = data %>% left_join(lm_results,by = c("PatientID" = "patient"))

slop_data_unique = slop_data %>% select(slope,PatientID) %>% distinct()

sv = slop_data_unique$slope[!is.na(slop_data_unique$slope)]

slope_th  = median(round(sv,2))

slop_data %>% mutate(type = case_when(
  slope >= slope_th  ~ "fast",
  slope < slope_th ~ "slow",
)) -> slop_data

slop_data = slop_data %>% mutate(group = NEC)
slop_data_unique = slop_data %>% select(group,slope,intercept,type,PatientID) %>% distinct()

contingency_table <- table(slop_data_unique$type,slop_data_unique$group)
contingency_table_str <- paste(capture.output(print(contingency_table)), collapse = "\n")

fisher_test <- fisher.test(contingency_table)
fisher_test
odds_vaule = fisher_test$estimate["odds ratio"]
fisher_p =  fisher_test[["p.value"]]
fisher_p 

data.frame(fisher_p = fisher_p,
           odds = odds_vaule,
           data = contingency_table_str) -> fisher.result
write.csv(fisher.result,"slow_fast/fast_nec_fisher.result.csv")


ggplot(data = slop_data %>% filter(is.na(type)), aes(x = truth, y = response,color = type,fill = type)) +
  geom_point(size = 3, alpha = .1) +
  stat_smooth(method = "lm", formula = y ~ x, 
              linewidth = 0.5, aes(group = PatientID),se = F,
              alpha = 0.1,color = "#E4E4E3") +
  stat_smooth(method = "lm", formula = y ~ x, linewidth = 1.5) +
  labs(x = "Actual day of life", y = "Predicted day of life")+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1) + 
  scale_color_manual(values = colpal) +
  scale_fill_manual(values = colpal) +
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank()) -> p2


ggplot(data = slop_data, aes(x = slope, y = intercept, color = type)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_vline(xintercept = slope_th, linetype = "dashed", color = "blue") +
  annotate(geom = "text",label = paste0("(",round(slope_th,4),")"),
           x = slope_th, y = 8.08)+
  scale_color_manual(values = colpal) +
  scale_fill_manual(values = colpal) +
  labs(x = "Slope", y = "Intercept") +
  coord_cartesian(xlim = c(-2.5,2.5))+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())-> p3

require(patchwork)

p3+p2+plot_layout(ncol = 2)+ plot_annotation(tag_levels = 'A')

ggsave("slow_fast//ident_fast_split_group.pdf",width = 10,height = 4)

require(ggstatsplot)
ggpiestats(
  data         = slop_data_unique,
  y            = group,
  x            = type,
  legend.title = "Sepsis",
)+ scale_fill_manual(values = colpal)

ggsave("slow_fast/fast_nec_fisher.result.pdf",width = 8,height = 6)

write.csv(slop_data,"slow_fast/Sepsis_add_trend_cluster.meta.csv")

slop_data_unique = slop_data %>% select(slope,intercept,type,PatientID) %>% distinct()

slop_data_unique = slop_data_unique %>% arrange(slope)
slop_data_unique$PatientID = factor(slop_data_unique$PatientID,levels = slop_data_unique$PatientID)
cplot = slop_data_unique %>% filter(!is.na(slope))
cplot$slope[cplot$slope >= 3] = 3
cplot$slope[cplot$slope <= -3] = -3
ggplot(data = cplot) +
  geom_linerange(aes(x = PatientID,ymin = 0,ymax = slope,color = type))+
  scale_color_manual(values =  c("#ff8bbf","#a0a5e1"))+
  ggthemes::theme_clean()+
  theme(axis.text.x = element_blank())+
  coord_polar()

ggsave("slow_fast/ident_fast_split_group_circ.pdf",width = 6,height = 4)

ggplot(data = slop_data %>% filter(slope > -3), aes(x = truth, y = response)) +
  geom_point(aes(fill = slope),size = 3, alpha = 0.2,color = "grey",shape = 21,show.legend = T) +
  stat_smooth(method = "lm", formula = y ~ x, 
              linewidth = 0.8, aes(group = PatientID),se = F,
              alpha = 0.6,color = "#E4E4E3") +
  stat_smooth(method = "lm", formula = y ~ x, 
              linewidth = 2,se = T,
              alpha = 0.6,color = "#D95050",fill = "#fbc7cc") +
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1) + 
  scale_fill_gradient2(low = "#b56bff",mid = "#d6c6ed",high = "#ff5a8c",midpoint = 0.36) +
  labs(x = "Age DOL", y = "Predicted Age")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank()) 

ggsave("slow_fast/ident_fast_split_group_line.pdf",width = 5,height = 4)


ggplot(data = slop_data_unique %>% filter(slope > -3)) +
  geom_boxplot(aes(y = slope,x = type,fill = type))+
  scale_fill_manual(values =  c("#ff8bbf","#a0a5e1"))+
  ggthemes::theme_clean()+
  theme(axis.text.x = element_blank(),
        panel.grid.major.y = element_blank())

ggsave("slow_fast/ident_fast_split_group_boxplot.pdf",width = 6,height = 4)


merge_sum = read.csv("../1.data_tidy/all_cohort_sample_num_8cohort_0905.csv",row.names = 1)
merge_sum = merge_sum %>% filter(Age_DOL <= 42)
data = merge_sum %>% filter(cohort == "Robert_2024_M")
slop_data = read.csv("slow_fast/Sepsis_add_trend_cluster.meta.csv",row.names = 1)

data = data %>% left_join(slop_data %>% select(Sample_ID,type),by = "Sample_ID")


sid = c("Enterococcus", "Klebsiella", "Staphylococcus", "Streptococcus", "Escherichia", "Bifidobacterium", "Clostridium", "Acinetobacter", "Lactobacillus", "Lactococcus")
sdata = data %>% select(c("Sample_ID",sid,"Age_DOL","type"))

sdata_line = sdata %>% pivot_longer(sid,names_to = "bac",values_to = "exp")
sdata_line$bac = factor(sdata_line$bac,levels = sid)
ggplot(sdata_line %>% filter(!is.na(type)), aes(Age_DOL, y = exp, color = type, fill = type)) + 
  geom_smooth(alpha = 0.2,method = "loess", span = 0.95) + 
  scale_color_manual(values = c("#1B8574","#bf48a8")) + 
  scale_fill_manual(values =  c("#1B8574","#bf48a8")) + 
  facet_wrap(".~bac",scales = "free",ncol = 5)+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Day", y = "Relative abundance")

ggsave("bac_diff/top9_shape_bac_smoothplot.pdf",width = 12,height = 4)

sdata_bar = sdata %>% pivot_longer(sid,names_to = "bac",values_to = "exp")
sdata_bar$bac = factor(sdata_bar$bac,levels = sid)

sdata_bar %>% rstatix::group_by(bac) %>%
  rstatix::wilcox_test(exp ~ type) %>%
  mutate(p = round(p,3))-> df_temp

sdata_bar$bac = factor(sdata_bar$bac,levels = rev(c(sid,"Others")))
require(ggradar)
sdata_rad = 
  sdata %>% pivot_longer(sid,names_to = "bac",values_to = "exp") %>%
  group_by(type,bac) %>%
  summarise(exp = mean(exp)) %>%
  pivot_wider(type,names_from = "bac",values_from = "exp")
sdata_rad = sdata_rad[,c("type",sid)]
require(ggradar)
ggradar(sdata_rad,grid.mid = 0.15,
        grid.max = 0.3,values.radar = c("0%", "15%", "30%"),)
write.csv(df_temp,"bac_diff//AgeGroup_top10_bac_radar_42.csv")
ggsave("bac_diff//AgeGroup_top10_bac_radar_42.pdf",width = 8,height = 6)


dir.create("DLpeptiase")
mnod = read_tsv("robert_NOD2.shortbred.count.txt")
mnod = mnod[apply(mnod[,-1],1,sum) != 0,]
ndata = data.frame(apply(mnod[,-1],2,sum)) 
ndata = ndata %>% rownames_to_column("SampleID")
colnames(ndata) = c("SampleID","nod2")

slop_data = read.csv("slow_fast/Sepsis_add_trend_cluster.meta.csv",row.names = 1)
m_meta = read_csv("robert_metegenome_meta.txt")
m_meta = m_meta %>% select(Run,`Sample Name`)

ndata1 = ndata %>% inner_join(m_meta,by = c("SampleID" = "Run")) %>%
  select(-SampleID)
colnames(ndata1) = c("nod2","Sample_ID")


data = ndata1 %>% 
  inner_join( slop_data %>% select(Sample_ID,PatientID,Age_DOL,type),by = "Sample_ID")

write.csv(data,"./DLpeptiase/DLpepiase_metagenome_sample.csv")

fit = lmer(nod2 ~ type * Age_DOL +(1|PatientID), data = data)
res = tidy(fit)
label = paste0("type: P=", round(res$p.value[2], 3), "\n Age_DOL: P=", round(res$p.value[3], 3), "\ntype * Age_DOL: P=", round(res$p.value[4], 3))

ggplot(data, aes(Age_DOL, nod2, color = type, fill = type)) +
  geom_smooth(alpha = 0.2,method = "loess", span = 0.95) +
  annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) +
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  scale_color_manual(values = c("#bf48a8","#1B8574")) +
  scale_fill_manual(values = c("#bf48a8","#1B8574")) +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "Day", y = "DL−endopeptidase expression(Count)")

ggsave("./DLpeptiase/NOD2_exp_smooth.pdf",width = 4,height = 3)

data %>% mutate(
  Age_Group = case_when(
    Age_DOL <= 7 ~ "0d-7d",
    Age_DOL <= 14 ~ "8d-14d",
    Age_DOL <= 21 ~ "15d-21d",
    Age_DOL <= 28 ~ "22d-28d",
    Age_DOL >= 29 ~ "29d-42d",
  )
) -> ndata1


ndata1$Age_Group = factor(ndata1$Age_Group,
                          levels = c( "0d-7d","8d-14d","15d-21d" ,"22d-28d", 
                                      "29d-42d"))
df_p_temp <- ndata1%>% 
  rstatix::group_by(Age_Group) %>% 
  rstatix::t_test(nod2 ~ type) %>% 
  rstatix::add_xy_position() %>%
  rstatix::add_significance(p.col = "p") %>%
  mutate(p = trunc(p*10000)/10000)

write_csv(df_p_temp,"./DLpeptiase/Age_group_score_diff_boxplot.csv")

ggplot(ndata1 %>% filter(!is.na(type)),aes(x = Age_Group, y = nod2))+
  geom_boxplot(aes(fill = type),alpha = 0.5,width = 0.5,outlier.shape = NA,outlier.stroke = NA) +
  geom_smooth(aes(group = type,color = type), se = FALSE)+
  geom_text(data = df_p_temp , inherit.aes = F,aes(x = Age_Group,
                                                   y = y.position/8, 
                                                   label = p),color = "red")+
  coord_cartesian(ylim = c(0,200))+
  scale_color_manual(values = c("#ff8bbf","#a0a5e1"))+
  scale_fill_manual(values = c("#ff8bbf","#a0a5e1"))+
  labs(y = "DL−endopeptidase expression(Count)",
       x = "DOL Group")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("./DLpeptiase/Age_group_NOD2_diff_boxplot.pdf",width = 5,height = 3)

mnod = read_tsv("robert_NOD2.shortbred.count.txt")

mnod = mnod[apply(mnod[,-1],1,sum) != 0,]

nod_info = read_csv("../16.NOD2_analysis//DLendo_class_new.csv")

mnod_sum = data.frame(mnod[,1],apply(mnod[,-1],1,sum))
colnames(mnod_sum) = c("ID","sum")

mnod_sum = mnod_sum %>% left_join(nod_info,by = "ID")

mnod_sum = mnod_sum %>% arrange(sum)

mnod_sum$ID = factor(mnod_sum$ID,levels = mnod_sum$ID)
require(patchwork)

colpal = readRDS("../colpal.rds")

mnod_sum %>% filter(!is.na(genus)) %>% 
  mutate(r = sum/sum(mnod_sum$sum)) %>% 
  write.csv("./DLpeptiase/NOD2_source_bac.csv")

vdata = read.csv("./DLpeptiase/NOD2_source_bac.csv",row.names = 1)

colnames(vdata)

vdata_plot = vdata %>% select(ID,genus,r) %>%
  mutate(lab = paste0(genus,".sp...","(",round(r*100,1),"%)","-",ID)) %>%
  arrange(desc(r))
vdata_plot$lab = factor(vdata_plot$lab,levels = vdata_plot$lab)
ggplot(vdata_plot,
       aes(x = "", y = r, fill = lab)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +  # 将坐标系设置为极坐标
  scale_fill_manual(values = xbox::need_colors(32)) +  # 设置颜色
  theme_classic() +
  theme(axis.line = element_blank(),  # 隐藏坐标轴线
        axis.text = element_blank(),  # 隐藏坐标轴标签
        axis.ticks = element_blank(),  # 隐藏坐标轴刻度
        axis.title = element_blank()) +  # 隐藏坐标轴标题
  labs(fill = "DL-endopeptidase Source")  # 设置图例标题
ggsave("./DLpeptiase/NOD2_source_bac.pdf",width = 10,height = 10)


library(ggalluvial)
cp = c("#889CCF", "#D95050", "#A0DCEC", "#7DCCC5", "#E0B794", "#DECFE6", "#F37E72", "#80B1D3", "#FBF6B5", "#F6B0E3", "#D9D9D8", "#C1A6EE", "#BAD5F7", "#A7D1CD", "#FDCEBA", "#7DCAE6", "#F4C573", "#F2AECB", "#ff8bbf", "#85dab1", "#fbc7cc", "#49c8cb", "#a0a5e1", "#ff8bbf", "#a0a5e1", "#ff8bbf")

ggplot(vdata_plot, aes(axis1 = "Start", axis3 = genus,axis2 = ID, y = r)) +
  geom_alluvium(fill = "#d9d9d8",color = "black") +
  geom_stratum(aes(fill = genus),width = 0.6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = cp)+
  theme_minimal() +
  labs(title = "Sankey Diagram using ggalluvial", x = "", y = "Normalized Weight")

ggsave("./DLpeptiase/NOD2_source_bac_sankey.pdf",width = 8,height = 4)
