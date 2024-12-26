xbox::chdir("../26.fastslow_vaild")

#以所有队列的的数据建模-------------
library(tidyverse)
library(mlr3verse)
library(mlr3extralearners)
library(mlr3tuning)
library(data.table)
library(pROC)
library(mlr3viz)
library(future)

rm(list = ls())
set.seed(202410)
colpal = readRDS("../colpal.rds")

c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv",row.names = 1)
c2_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv",row.names = 1)

c1_meta %>% rename(
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
  mutate(Weight = Weight*1000)-> c1_meta

nid = intersect(colnames(c2_meta),colnames(c1_meta))

# share 47 genus--------
alldata =  read.csv("../1.data_tidy/all_cohort.sum_7cohort_final.csv",row.names = 1)
alldata = alldata %>% mutate(PatientID = paste0(cohort,"_",PatientID))
sid = colnames(alldata)[2:48]
nid = c("SampleID",sid,nid[303:323])

all_meta = rbind(c1_meta[,nid],c2_meta[,nid]) %>%
  mutate(PatientID = paste0(cohort,"_",PatientID))

data = all_meta %>% filter(Age_DOL <= 42)
write.csv(data,"c1c2_merge_47g_42d.csv")


  filter(Age_DOL <= Onset_Day | is.na(Onset_Day)) %>%
  rename("Sample_ID" = "SampleID")

data_onset1 = data_onset %>% 
  filter(Age_DOL <= 10) %>% 
  filter(!Onset_Day < 10| is.na(Onset_Day))

all_data = data_onset1[,colnames(alldata)]

# load model-----
learner = readRDS("./model_pridict_learner.rds")

test_c = all_data[,-c(49,51,52,53)]

c_task = as_task_regr(test_c, target = "Age_DOL", id = "regression")

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
write.csv(regr.result,"regr.result_test.csv")

best_tuned_pred <- data.frame(c_prediction$print())
age_prd = best_tuned_pred

best_tuned_pred$PatientID = all_data$PatientID
best_tuned_pred$Sample_ID =  all_data$Sample_ID
best_tuned_pred$Age_DOL =  all_data$Age_DOL
best_tuned_pred$Sepsis_style =  all_data$Sepsis_style
best_tuned_pred$NEC =  all_data$NEC
best_tuned_pred$cohort =  all_data$cohort

pre_age_pdata = best_tuned_pred %>% filter(Sepsis_style != "EOS")

pre_age_pdata = best_tuned_pred %>% filter(Sepsis_style != "EOS")
write.csv(pre_age_pdata,"pre_age_pdata_10d.csv")


temp = as.data.frame(unique(pre_age_pdata[,c("PatientID","Sepsis_style")]))
annot_data = cbind(as.data.frame(table(temp['Sepsis_style'])),as.data.frame(table(pre_age_pdata$Sepsis_style)))
colnames(annot_data) =c("group","patient","g","sample")

ggplot(data = pre_age_pdata , aes(x = truth, y = response,color = Sepsis_style)) +
  geom_point(alpha = 0.03)+
  stat_smooth(method = "lm", formula = y ~ x, linewidth = 1.5, aes(color = Sepsis_style,)) +
  labs(x = "Actual day of life", y = "Predicted day of life")+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1) +
  geom_text(data = annot_data, aes(x = Inf, y = -Inf,color = group,label = paste0(g,"S",sample, "P", patient)),hjust = 1,vjust = -1) +
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("Sepsis_fastslow_by_cohort_nopoint.png",width = 4,height = 2.8)
ggsave("Sepsis_fastslow_by_cohort_nopoint.pdf",width = 4,height = 2.8)


re = data %>% 
  select(PatientID,Onset_Day) %>%
  distinct()


temp = as.data.frame(unique(pre_age_pdata[,c("PatientID","NEC")]))
annot_data = cbind(as.data.frame(table(temp['NEC'])),as.data.frame(table(pre_age_pdata$NEC)))
colnames(annot_data) =c("group","patient","g","sample")

ggplot(data = pre_age_pdata %>% filter(!is.na(NEC)), aes(x = truth, y = response,color = NEC)) +
  geom_point(alpha = 0.03)+
  stat_smooth(method = "lm", formula = y ~ x, linewidth = 1.5, aes(color = NEC)) +
  labs(x = "Actual day of life", y = "Predicted day of life")+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1) +
  geom_text(data = annot_data, aes(x = Inf, y = -Inf,color = group,label = paste0(g,"S",sample, "P", patient)),hjust = 1,vjust = -1) +
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("NEC_fastslow_by_cohort_nopoint.png",width = 4,height = 2.8)
ggsave("NEC_fastslow_by_cohort_nopoint.pdf",width = 4,height = 2.8)

pre_age_pdata %>% mutate(
  status = case_when(
    is.na(Sepsis_style) & is.na(NEC) ~ "Unspecified Health Status",
    Sepsis_style == "N" & NEC == "N" ~ "Apparent Health",
    Sepsis_style != "N" | NEC != "N" ~ "Diseased"
  )
) -> pre_age_pdata

temp = as.data.frame(unique(pre_age_pdata[,c("PatientID","status")]))
annot_data = cbind(as.data.frame(table(temp['status'])),as.data.frame(table(pre_age_pdata$status)))
colnames(annot_data) =c("group","patient","g","sample")

ggplot(data = pre_age_pdata %>% filter(!is.na(status)), aes(x = truth, y = response,color = status)) +
  geom_point(alpha = 0.06)+
  stat_smooth(method = "lm", formula = y ~ x, linewidth = 1.5, aes(color = status)) +
  labs(x = "Actual day of life", y = "Predicted day of life")+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1) +
  geom_text(data = annot_data, aes(x = Inf, y = -Inf,color = group,label = paste0(g,"S",sample, "P", patient)),hjust = 1,vjust = -1) +
  scale_fill_manual(values = c("#a0a5e1","#ff8bbf"))+
  scale_color_manual(values = c("#a0a5e1","#ff8bbf"))+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("heathy_fastslow_by_cohort_nopoint.png",width = 5,height = 2.8)
ggsave("heathy_fastslow_by_cohort_nopoint.pdf",width = 5,height = 2.8)

colnames(pre_age_pdata)
source("../addin.r")
diff_slope_data = data.frame()
table(pre_age_pdata$status)
for(i in c("Sepsis_style","NEC","status")){
  data  = pre_age_pdata %>% rename(group = all_of(i))
  group_ls = unique(data$group)
  group_ls = group_ls[!is.na(group_ls)]
  if(i == "status"){
    group_ls= c("Apparent Health", "Diseased")
  }
  
  healthy_df <- subset(data,group==group_ls[2])
  dis_df <- subset(data,group==group_ls[1])
  
  model_healthy <- lm(response ~ truth,data = healthy_df)
  model_disease <- lm(response ~ truth,data = dis_df)
  
  resid_healthy <- residuals(model_healthy)
  resid_disease <- residuals(model_disease)
  
  healthy_df$residuals <- resid_healthy
  dis_df$residuals <- resid_disease
  
  test <- wilcox.test(healthy_df$residuals,dis_df$residuals)
  stat_df <- data.frame("wilcox.test",
                        statistic = test$statistic,
                        p.value = test$p.value,
                        i)
  diff_k <- diffslope(x1 = healthy_df$truth,y1 = healthy_df$response,
                      x2 = dis_df$truth,y2 = dis_df$response,permutations = 1000)
  dif_df <- data.frame(group_name = "slope_diff",
                       diff = diff_k$slope.diff,
                       p.value = diff_k$signif,
                       i)
  diff_slope_data = bind_rows(diff_slope_data,as.data.frame(dif_df),
                              as.data.frame(stat_df))
}

write.csv(diff_slope_data,"Sepsis_nec_heathy_slope_diff.csv")

colnames(pre_age_pdata)

ggplot(data = pre_age_pdata %>% filter(!is.na(NEC)), aes(x = truth, y = response,
                                                         color = status,fill = status)) +
  geom_point(aes(color = status),size = 0.7, alpha = .2) +
  stat_smooth(method = "lm", formula = y ~ x, 
              linewidth = 0.5, aes(group = PatientID),se = F,
              alpha = 0.1,color = "#E4E4E3") +
  stat_smooth(method = "lm", formula = y ~ x, 
              linewidth = 1.5, aes(color = status),se = T) +
  theme_bw() +
  labs(x = "Actual day of life", y = "Predicted day of life")+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1) + 
  coord_cartesian(ylim= c(3,28))+
  scale_fill_manual(values = c("#a0a5e1","#ff8bbf"))+
  scale_color_manual(values = c("#a0a5e1","#ff8bbf"))+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("Heathy_every_patient_sepsis_line.pdf",width = 4.5,height = 2.8)

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
slope_th = 1.35
slop_data %>% mutate(type = case_when(
  slope >= slope_th  ~ "fast",
  slope < slope_th ~ "slow",
)) -> slop_data
slop_data = slop_data %>% rename("group" = "Sepsis_style")

slop_data_unique = slop_data %>% filter(group != "EOS" | is.na(group))%>% select(group,slope,intercept,type,PatientID) %>% distinct()

contingency_table <- table(slop_data_unique$type,slop_data_unique$group)
contingency_table 
contingency_table_str <- paste(capture.output(print(contingency_table)), collapse = "\n")

fisher_test <- fisher.test(contingency_table)
chi_test = chisq.test(contingency_table)

odds_vaule = fisher_test$estimate["odds ratio"]
fisher_p =  chi_test[["p.value"]]
fisher_p 

data.frame(fisher_p = fisher_p,
           odds = odds_vaule,
           data = contingency_table_str) -> fisher.result
write.csv(fisher.result,"Sepsis_fisher.result.csv")
ggplot(data = slop_data, aes(x = truth, y = response,color = type,fill = type)) +
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

ggsave("ident_fast_split_group.pdf",width = 10,height = 4)

require(ggstatsplot)
ggpiestats(
  data         = slop_data_unique,
  y            = group,
  x            = type,
  legend.title = "Sepsis",
)+ scale_fill_manual(values = colpal)

ggsave("Sepsis_ratio_sepsis_patient.pdf",width = 8,height = 6)

save.image("Sepsis_env_var_240410.rdata")

otu_data = alldata %>% filter(Age_DOL <= 42) %>%
  left_join(slop_data_unique  %>% select(PatientID,type),by = "PatientID")

write.csv(otu_data,"allsmaple_Sepsis_add_trend_cluster.meta.csv")

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

ggsave("ident_fast_split_group_circ.pdf",width = 6,height = 4)

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

ggsave("ident_fast_split_group_line.pdf",width = 5,height = 4)


ggplot(data = slop_data_unique %>% filter(slope > -3),
       aes(y = slope,x = type,fill = type)) +
  geom_boxplot()+
  ggpubr::stat_compare_means()+
  scale_fill_manual(values =  c("#ff8bbf","#a0a5e1"))+
  ggthemes::theme_clean()+
  theme(axis.text.x = element_blank(),
        panel.grid.major.y = element_blank())

ggsave("ident_fast_split_group_boxplot.pdf",width = 3,height = 4)

#OR分析-----------
rm(list = ls())
c1_fs = read.csv("./Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv",row.names = 1)
all_meta = c1_fs %>% select(SampleID,type) %>% inner_join(c1_meta,by = c("SampleID" = "ID"))

c1_fs = read.csv("./Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv",row.names = 1)
all_meta2 = c1_fs %>% select(SampleID,type) %>% inner_join(c1_meta,by = c("SampleID" = "SampleID"))

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

data = all_meta %>% filter(Age_DOL <=42) %>%
  mutate(trend_cluster = type) 

all_meta = all_meta %>% filter(is.na(Onset_Day) | (Onset_Day < Age_DOL))


data = all_meta %>% filter(Age_DOL <=42) %>%
  mutate(trend_cluster = type) 

table(data$Feeding)

meta = data[304:325]

colnames(meta)

library(tidyverse)

meta %>% group_by(PatientID) %>%
  slice_head(n = 1) %>%
  ungroup() -> meta

meta$Sepsis_style[is.na(meta$Sepsis_style)] = "N"
df = meta %>% filter(Sepsis_style %in% c("LOS","N"))
df = meta
df <- hablar::retype(df)

df %>% mutate(
  Sepsis_style = case_when(
    Sepsis_style == "N" ~ "N",
    Sepsis_style == "LOS" ~ "Y",
  )
) -> df

require(ggstatsplot)

df$trend_cluster <- ifelse(df$trend_cluster == "fast", 0, 1)

results <- data.frame(Factor = character(), Test = character(), Statistic = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

cols_to_convert <- colnames(df)[c(4,6,7,8,9,10,11,12,13,
                                  21)]
convert_to_binary <- function(x) {
  ifelse(x == "Y" | x == "F", 1, 0)  # "B" 和 "Y" 转换为 1，其他转换为 0
}
df <- df %>%
  mutate(across(all_of(cols_to_convert), convert_to_binary))

cols_to_scale <- colnames(df)[c(14,15,16,18,19)]
convert_to_scale <- function(x){
  z = (x - min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))
}
df <- df %>%
  mutate(across(all_of(cols_to_scale), convert_to_scale))


df %>% mutate(Feeding = case_when(
  Feeding == "A" ~ 0,
  Feeding == "B" ~ 1,
  Feeding == "C" ~ 0.5,
)) -> df

library(lme4)
library(tidyverse)
library(lmerTest)
library(broom.mixed)

results <- data.frame(Factor = character(), Test = character(), Statistic = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

for(factor in colnames(df)[-c(1,2,3,5,17,22)]){
  test_data = data.frame("trend" = df$trend_cluster, "factor" = df[[factor]])
  
  glm_model <- glm(trend ~ factor, data = test_data, family = binomial)
  summary_glm <- summary(glm_model)
  
  if (any(grepl("factor",rownames(summary_glm$coefficients)))) {
    id = grep("factor",rownames(summary_glm$coefficients),value = T)
    coeff <- summary_glm$coefficients[id, "Estimate"]
    odds_value <- exp(coeff)
    conf_int <- 0
    lower_conf <- 0
    upper_conf <- 0
  } else {
    odds_value <- NA
    lower_conf <- NA
    upper_conf <- NA
  }
  
  # 将结果添加到结果数据框中
  results <- rbind(results, data.frame(Factor = factor, 
                                       Test = "Logistic Regression",
                                       Statistic = coeff,
                                       lower_conf = lower_conf,
                                       upper_conf = upper_conf,
                                       P_value = summary_glm$coefficients[id, "Pr(>|z|)"], 
                                       Mean1 = NA, 
                                       Mean2 = NA, 
                                       Median1 = NA, 
                                       Median2 = NA,
                                       Ratio = NA, 
                                       OddsRatio = odds_value, 
                                       ContingencyTable = NA, 
                                       stringsAsFactors = FALSE))
}

write.csv(results,"mergec1c2_fast_vs_slow_meta_Diff_glm_OR.csv")

# 绘图
results %>% 
  mutate(ratio = case_when(
    !is.na(OddsRatio) ~ as.numeric(OddsRatio),
    is.na(OddsRatio) ~ Ratio
  )
  ) -> data_plot


data_plot %>% arrange(Test,ratio) -> data_plot
colnames(data_plot)
data_plot$Factor

did = c("BPD","WMI", "NEC","Sepsis_style")

data_plot1 = data_plot %>% filter(Factor %in% did)
data_plot2 = data_plot %>% 
  filter(!Factor %in% did) %>% 
  filter(!is.na(Factor)) %>%
  filter(Factor != "cohort")

data_plot1$Factor = factor(data_plot1$Factor,levels =  did)

ggplot()+
  geom_col(data= data_plot1 ,aes(y = Factor, x = ratio),fill = "#D95050",
           color = "#d9d9d8")+
  geom_text(data= data_plot1,
            aes(y = Factor, x = ratio,label = round(P_value,2)),
            hjust = -0.5,
            color = "#D95050")+
  geom_vline(xintercept = 1,linetype = "dashed",color = "grey")+
  labs(x = "Odds ratio",y = "")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank()) -> p1

did = c("Natural_delivery",
        "Intrauterine_infection",
        "PE",
        "Feeding","Gender","GA_daily","Weight",
        "Antibiotics_Duration_all",
        "Antibiotics_Duration_2w")
data_plot$Ratio
data_plot2 = data_plot %>% 
  filter(Factor %in% did) %>% 
  filter(!is.na(Factor)) %>%
  filter(Factor != "cohort") %>%
  arrange(ratio)

data_plot2$Factor =  factor(data_plot2$Factor,levels = unique(data_plot2$Factor))
require(patchwork)

data_plot2 %>% mutate(ratio1 = log(ratio))-> data_plot2


ggplot()+
  geom_col(data= data_plot2 ,aes(y = Factor, x = ratio1),width = 0.9)+
  geom_col(data= data_plot2 %>% filter(P_value <= 0.05),
           aes(y = Factor, x = ratio1),fill = "#85dab1",width = 0.9)+
  geom_text(data= data_plot2 %>% filter(P_value <= 0.05),
            aes(y = Factor, x = ratio1,label = "*"),
            color = "red",size = 10)+
  geom_vline(xintercept = 0,color = "black")+
  labs(x = "Log(Ratio)",y = "")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())-> p2

p1/p2+plot_layout(heights = c(1, 2))


ggsave("mergec1c2_fast_vs_slow_meta_Diff_glm_or.pdf",width = 4,height = 4)

require(tidyverse)
library(survival)
library(ggsurvfit)
library(survminer)
require(survminer)
library(gridExtra)
rm(list = ls())

c1_fs = read.csv("./Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv",row.names = 1)
all_meta = c1_fs %>% select(SampleID,type) %>% inner_join(c1_meta,by = c("SampleID" = "ID"))

c1_fs = read.csv("./Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv",row.names = 1)
all_meta2 = c1_fs %>% select(SampleID,type) %>% inner_join(c1_meta,by = c("SampleID" = "SampleID"))

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

ndata = read.csv("../1.data_tidy/all_cohort.sum_7cohort_final.csv",row.names = 1)
nid = c(nid[1:2],colnames(ndata)[2:48],nid[304:324])

all_meta = rbind(all_meta[,nid],all_meta2[,nid])

all_meta %>% 
  mutate(
    OS.time.y = case_when(
      !is.na(Onset_Day) ~ Onset_Day,
      TRUE ~ 42
    ) 
  ) %>% mutate(group = type,OS = case_when(
    Sepsis_style == "LOS" ~ 1,
    Sepsis_style == "N" ~ 0,
  )) -> sur_data2

table(sur_data2$Feeding)

# bread or not bread
sur_data = sur_data2 %>% mutate(
  Feeding = case_when(
    Feeding == "A" ~ "A",
    Feeding != "A" ~ "B"
  )
)

sur_data %>% group_by(PatientID) %>%
  #slice_max(order_by =  OS.time.y, n = 1) %>%
  slice_head(n=1)-> sur_data



#sur_data = sur_data2

surv_data <- with(sur_data, Surv(OS.time.y, OS))
fit <- survfit2(surv_data ~ group, data = sur_data)

ggsurvplot(fit, data = sur_data,
           fun = "event",
           pval = TRUE, conf.int = FALSE,
           conf.int.style="step",
           surv.scale = "percent",
           palette = c("#ff8bbf","#a0a5e1"),
           censor = F,
           censor.size = 10,
           ylim = c(0,0.3),
           risk.table = TRUE,
           risk.table.title = "Number at risk",
           risk.table.height = 0.45) -> p1

cox_model <- coxph(Surv(OS.time.y, OS) ~ group, 
                   data = sur_data)
model_summary <- summary(cox_model)
hazard_ratios <- exp(coef(cox_model))
conf_intervals <- confint(cox_model,level = 0.95)
p_values <- model_summary$coefficients[, "Pr(>|z|)"]
results_df <- data.frame(
  Variable = rownames(model_summary$coefficients),
  Hazard_Ratio = round(hazard_ratios, 2),
  Lower_CI = hazard_ratios + round(conf_intervals[,1], 2),
  Upper_CI = hazard_ratios +  round(conf_intervals[,2], 2),
  P_Value = round(p_values, 4),
  Class = "OS"
)

write_csv(results_df,"cohort1_cox_p_HR_result_single.csv")

g1 <- ggplotGrob(p1$plot)
g1
ggsave("cohort1_cox_curve.pdf",g1 ,width = 4, height = 3.5, units = "in")
p1
ggsave("cohort1_cox_curve_table.pdf" ,width = 4, height = 2, units = "in")

library(tidyverse)
library(lmerTest)
library(broom.mixed)

rm(list = ls())
ndata = read.csv("../1.data_tidy/all_cohort.sum_7cohort_final.csv",row.names = 1)
data = ndata %>% filter(Age_DOL <= 42)

fastslow = read.csv("./allsample_Sepsis_add_trend_cluster.meta.csv",row.names = 1)

data = data %>% 
  inner_join(fastslow %>%
               select(Sample_ID,type),
             by = c("Sample_ID" = "Sample_ID"))


write.csv(data,"plot_data_c1c2.csv")

sid = readRDS("../12.shap_value/top10.rds")

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

ggsave("top9_shape_bac_smoothplot_all_type_all_10d_1.35.pdf",width = 12,height = 4)

rm(list = ls())
mnod = read_tsv("../16.NOD2_analysis/NOD2.shortbred.count.txt")
mnod = mnod[apply(mnod[,-1],1,sum) != 0,]

ndata = data.frame(apply(mnod[,-1],2,sum)) 
ndata = ndata %>% rownames_to_column("SampleID")
colnames(ndata) = c("SampleID","nod2")

slop_data = read.csv("./Sepsis_add_trend_cluster.meta.csv",row.names = 1)

data = ndata %>% 
  inner_join(slop_data %>% select(SampleID,PatientID,Age_DOL,type),by = "SampleID")

write.csv(data,"./DLpepiase_metagenome_sample.csv")
fit = lmer(nod2 ~ type * Age_DOL +(1|PatientID), data = data)
res = tidy(fit)
label = paste0("type: P=", round(res$p.value[2], 3), "\n Age_DOL: P=", round(res$p.value[3], 3), "\ntype * Age_DOL: P=", round(res$p.value[4], 3))

ggplot(data %>% filter(!is.na(type)), aes(Age_DOL, nod2, color = type, fill = type)) +
  geom_smooth(alpha = 0.2,method = "loess", span = 0.95) +
  annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) +
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  scale_color_manual(values = c("#bf48a8","#1B8574")) +
  scale_fill_manual(values = c("#bf48a8","#1B8574")) +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "Day", y = "DL−endopeptidase expression(Count)")

ggsave("./NOD2_exp_smooth.pdf",width = 4,height = 3)

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

write_csv(df_p_temp,"./Age_group_score_diff_boxplot.csv")

ggplot(ndata1 %>% filter(!is.na(type)),aes(x = Age_Group, y = nod2))+
  geom_boxplot(aes(fill = type),alpha = 0.5,width = 0.5,outlier.shape = NA,outlier.stroke = NA) +
  geom_smooth(aes(group = type,color = type), se = FALSE)+
  geom_text(data = df_p_temp , inherit.aes = F,aes(x = Age_Group,
                                                   y = y.position/8, 
                                                   label = p),color = "red")+
  coord_cartesian(ylim = c(0,300))+
  scale_color_manual(values = c("#ff8bbf","#a0a5e1"))+
  scale_fill_manual(values = c("#ff8bbf","#a0a5e1"))+
  labs(y = "DL−endopeptidase expression(Count)",
       x = "DOL Group")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("./Age_group_NOD2_diff_boxplot.pdf",width = 5,height = 3)



rm(list = ls())
c1_fs = read.csv("./Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv",row.names = 1)
all_meta = c1_fs %>% select(SampleID,type) %>% inner_join(c1_meta,by = c("SampleID" = "ID"))

c1_fs = read.csv("./Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv",row.names = 1)
all_meta2 = c1_fs %>% select(SampleID,type) %>% inner_join(c1_meta,by = c("SampleID" = "SampleID"))

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

data = all_meta %>% filter(Age_DOL <=42) %>%
  mutate(trend_cluster = type) 

all_meta = all_meta %>% filter(is.na(Onset_Day) | (Onset_Day < Age_DOL))


data = all_meta %>% filter(Age_DOL <=42) %>%
  mutate(trend_cluster = type) 

pre_age_pdata = read.csv("pre_age_pdata_10d.csv")
data = data %>% left_join(pre_age_pdata %>% select(truth,response,Sample_ID),by = c("SampleID" = "Sample_ID"))

data %>% mutate(
  status = case_when(
    is.na(Sepsis_style) & is.na(NEC) & is.na(BPD) & is.na(WMI) ~ "Unspecified Health Status",
    Sepsis_style == "N" & NEC == "N" & BPD== "N" & WMI == "N"~ "Apparent Health",
    Sepsis_style != "N" | NEC != "N" | BPD!= "N"| WMI != "N"~ "Diseased"
  )
) -> pre_age_pdata

temp = as.data.frame(unique(pre_age_pdata[,c("PatientID","status")]))
annot_data = cbind(as.data.frame(table(temp['status'])),as.data.frame(table(pre_age_pdata$status)))
colnames(annot_data) =c("group","patient","g","sample")

ggplot(data = pre_age_pdata %>% filter(!is.na(status)), aes(x = truth, y = response,color = status)) +
  geom_point(alpha = 0.06)+
  stat_smooth(method = "lm", formula = y ~ x, linewidth = 1.5, aes(color = status)) +
  labs(x = "Actual day of life", y = "Predicted day of life")+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1) +
  geom_text(data = annot_data, aes(x = Inf, y = -Inf,color = group,label = paste0(g,"S",sample, "P", patient)),hjust = 1,vjust = -1) +
  scale_fill_manual(values = c("#a0a5e1","#ff8bbf"))+
  scale_color_manual(values = c("#a0a5e1","#ff8bbf"))+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("4disease_fastslow_by_cohort_nopoint.pdf",width = 4,height = 2.8)

source("../addin.r")
diff_slope_data = data.frame()
table(pre_age_pdata$status)
for(i in c("status")){
  data  = pre_age_pdata %>% rename(group = all_of(i))
  group_ls = unique(data$group)
  group_ls = group_ls[!is.na(group_ls)]
  if(i == "status"){
    group_ls= c("Apparent Health", "Diseased")
  }
  
  healthy_df <- subset(data,group==group_ls[2])
  dis_df <- subset(data,group==group_ls[1])
  
  model_healthy <- lm(response ~ truth,data = healthy_df)
  model_disease <- lm(response ~ truth,data = dis_df)
  
  resid_healthy <- residuals(model_healthy)
  resid_disease <- residuals(model_disease)
  
  healthy_df$residuals <- resid_healthy
  dis_df$residuals <- resid_disease
  
  test <- wilcox.test(resid_healthy,dis_df$residuals)
  stat_df <- data.frame("wilcox.test",
                        statistic = test$statistic,
                        p.value = test$p.value,
                        i)
  healthy_df_clean <- na.omit(healthy_df[, c("truth", "response")])
  dis_df_clean <- na.omit(dis_df[, c("truth", "response")])
  
  diff_k <- diffslope(
    x1 = healthy_df_clean$truth,
    y1 = healthy_df_clean$response,
    x2 = dis_df_clean$truth,
    y2 = dis_df_clean$response,
    permutations = 1000
  )
  
  dif_df <- data.frame(group_name = "slope_diff",
                       diff = diff_k$slope.diff,
                       p.value = diff_k$signif,
                       i)
  diff_slope_data = bind_rows(diff_slope_data,as.data.frame(dif_df),
                              as.data.frame(stat_df))
}

write.csv(diff_slope_data,"4disease_slope_diff.csv")
