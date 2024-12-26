xbox::chdir("../5.factor")
require(tidyverse)

rm(list = ls())
set.seed("123")
c1_fs = read.csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv",row.names = 1)
all_meta = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "ID"))

c1_fs = read.csv("../4.bac_slow_fast/Sepsis_add_trend_cluster.meta.csv",row.names = 1)
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

glm_model <- glm(trend_cluster ~ 
                   Antibiotics_Duration_2w+
                   Natural_delivery + 
                   Gender + 
                   Feeding+
                   PE+
                   Weight+
                   GA_daily+
                   Intrauterine_infection, data = df,family = "binomial")

summary_glm <- summary(glm_model)
summary_glm = tidy(glm_model)
summary_glm = summary_glm %>% mutate(OR = exp(estimate))

write.csv(summary_glm,"glm_mix_mutivar_stat_addALL.csv")

glm_model <- glm(trend_cluster ~ 
                   Antibiotics_Duration_2w+
                   Natural_delivery + 
                   Gender + 
                   Feeding+
                   PE+
                   GA_daily+
                   Intrauterine_infection, data = df,family = "binomial")

summary_glm <- summary(glm_model)
summary_glm = tidy(glm_model)
summary_glm = summary_glm %>% mutate(OR = exp(estimate))

write.csv(summary_glm,"glm_mix_mutivar_stat_addGA.csv")

glm_model <- glm(trend_cluster ~ 
                   Antibiotics_Duration_2w+
                   Natural_delivery + 
                   Gender + 
                   Feeding+
                   PE+
                   Weight+
                   Intrauterine_infection, data = df,family = "binomial")

summary_glm <- summary(glm_model)
summary_glm = tidy(glm_model)
summary_glm = summary_glm %>% mutate(OR = exp(estimate))

write.csv(summary_glm,"glm_mix_mutivar_stat_addWeight.csv")

require(ggcor)
quickcor(df[,c(4,6,7,12,14,15,16,20,21,22)],cor.test = TRUE, use = "complete.obs") +
  geom_square(data = get_data(type = "lower", show.diag = FALSE)) +
  geom_mark(data = get_data(type = "upper", show.diag = FALSE), size = 2.5) +
  geom_abline(slope = -1, intercept = 11)+
  scale_fill_gradient2(low = "#41ab5d",mid = "#f1fb67",high = "#e4007c")

ggsave("mutivar_cor_spearman.pdf",width = 8,height = 6)

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

results =  read.csv("mergec1c2_fast_vs_slow_meta_Diff_glm_OR.csv",row.names = 1)
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
  geom_col(data= data_plot1 ,aes(y = Factor, x = ratio))+
  geom_col(data= data_plot1 %>% filter(P_value <= 0.05),
           aes(y = Factor, x = ratio),fill = "#D95050")+
  geom_vline(xintercept = 1,linetype = "dashed",color = "grey")+
  labs(x = "Odds ratio",y = "")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank()) -> p1

did = c("Natural_delivery",
        "Intrauterine_infection",
        "PE",
        "Feeding","Gender","GA_daily","Weight",
        "Antibiotics_Duration_all",
        "Antibiotics_Duration_2w",
        "GDM")
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
ggsave("mergec1c2_fast_vs_slow_meta_Diff_glm_or.png",width = 4,height = 4)


xbox::heatpoint(df$Weight,df$GA_daily) -> re
xbox::xplot(re)+
  labs(x = "Weight",
       y = "GA_daily")

df = meta %>% filter(Sepsis_style %in% c("LOS","N"))
df[, c("Sepsis_style","NEC","trend_cluster", "Gender", "Weight", "Feeding", 
       "GA_daily", "Natural_delivery", "Intrauterine_infection", 
       "PE", "GDM")] -> df_m

df_m_numeric <- model.matrix(~.-1, data = df_m)

cor_matrix <- cor(df_m_numeric, use = "complete.obs", method = "spearman")
write.csv(cor_matrix,"mutivar_cor_spearman.csv")

require(car)
fit <- lm(Sepsis_styleLOS ~ ., data = as.data.frame(df_m_numeric))
vif_values <- vif(fit)
write.csv(vif_values,"mutivar_vif_matrix.csv")

setwd("../5.factor/")
library(phyloseq)
require(microbiomeutilities)
require(ggpubr)
require(microbiome)
require(vegan)

otu_data2 = data[,c(3:303)]
phe = data[,c(304:325)]

otu_exp1 = otu_data2

set.seed(666)
source("../addin.r")

meio.mdf <- as.matrix.data.frame(otu_exp1)
rownames(meio.mdf) <- rownames(otu_exp1)
meio.mdf_clean <- meio.mdf[rowSums(meio.mdf) != 0, ]
arg_distance <- vegdist(meio.mdf_clean, method = "bray")

s_id = attributes(arg_distance)$Labels
results=NULL

for(i in colnames(phe)){
  adonis_re <- tryCatch(
    {
      adonis2(arg_distance ~ phe[s_id, ][[i]], permutations = 1000, na.action = na.omit)
    },
    error = function(e) {
      message(paste("Error occurred in adonis2 function for column", i, ":", e))
      return(NULL)
    }
  )
  
  if (!is.null(adonis_re)) {
    summary(adonis_re)
    rownames(adonis_re)[1] = i
    results <- rbind(results, adonis_re[1, ])
  }
}

write.csv(results,"merge_adonis_clin_data_R2_result.csv")

adonis_data = results

adonis_data$R2_2 = round(adonis_data$R2, digits = 4)

adonis_data$R2_2 = paste0(adonis_data$R2_2*100,"%")

adonis_data %>% rownames_to_column("feature") %>% arrange(R2) -> adonis_data

adonis_data$feature = factor(adonis_data$feature , levels = unique(adonis_data$feature))

adonis_data %>%
  mutate(pif = case_when(
    `Pr(>F)`>0.05 ~ "ns",
    `Pr(>F)`<0.001 ~ "***",
    `Pr(>F)`<0.01 ~ "**",
    `Pr(>F)`<= 0.05 ~ "*",
  )) -> adonis_data

ggplot(adonis_data , aes(x = R2, y = feature, fill = R2))+
  geom_col()+
  geom_text(aes(label = pif,hjust =0),size = 3,color = "red")+
  scale_fill_gradientn(colors = c('navyblue', 'darkmagenta', 'darkorange1')) +
  theme_classic()
ggsave("merge_adonis_clin_data_R2_result.pdf", width = 6, height = 10)

ggplot(adonis_data %>% 
         filter(!feature %in% c("PatientID","Group_Time","Onset_Day",
                                "trend_cluster","cohort")) , aes(x = R2, y = feature, fill = R2))+
  geom_col()+
  geom_text(aes(label = pif,hjust =0),size = 3,color = "red")+
  scale_fill_gradientn(colors = c('navyblue', 'darkmagenta', 'darkorange1')) +
  theme_classic()
ggsave("merge_adonis_clin_data_R2_result.pdf", width = 5, height = 3)

setwd("./5.factor/")
require(tidyverse)
results = read.csv("merge_adonis_clin_data_R2_result.csv",row.names = 1)

adonis_data = results

adonis_data$R2_2 = round(adonis_data$R2, digits = 4)

adonis_data$R2_2 = paste0(adonis_data$R2_2*100,"%")

adonis_data %>% rownames_to_column("feature") %>% arrange(R2) -> adonis_data

nid = c("Age_DOL",
        "Feeding",
        "Antibiotics_Duration_all",
        "Natural_delivery",
        "Weight",
        "GDM",
        "PE",
        "GA_daily",
        "Gender",
        "Intrauterine_infection")

adonis_data  = adonis_data %>% filter(feature %in% nid)

adonis_data$R2_2 = round(adonis_data$R2, digits = 4)

adonis_data$R2_2 = paste0(adonis_data$R2_2*100,"%")

adonis_data %>% arrange(R2) -> adonis_data
colnames(adonis_data)[1] = "feature"
adonis_data$feature = factor(adonis_data$feature , levels = unique(adonis_data$feature))

adonis_data %>%
  mutate(pif = case_when(
    `Pr..F.`>0.05 ~ "ns",
    `Pr..F.`<0.001 ~ "***",
    `Pr..F.`<0.01 ~ "**",
    `Pr..F.`<= 0.05 ~ "*",
  )) %>%
  mutate(test = case_when(
    `Pr..F.`>0.05 ~ "ns",
    TRUE ~ "signif"
  )) -> adonis_data1

adonis_data1 = adonis_data1 %>% 
  arrange(R2)
adonis_data1$feature = factor(adonis_data1$feature,levels = 
                                adonis_data1$feature)
ggplot(adonis_data1 , aes(x = R2, y = feature))+
  geom_col(aes(fill = test))+
  geom_text(aes(label = pif,hjust =0),size = 3,color = "red")+
  scale_fill_manual(values = c("ns" = "#e5e5d8", "signif" = "#00d4c6")) +
  coord_cartesian(xlim = c(0,0.08))+
  labs(y = "CALM2005_Merge")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y  = element_blank()) -> p3

require(tidyverse)
xbox::chdir("../5.factor")
nid = c("Age_DOL",
        "Feeding",
        "Intrauterine_infection",
        "Antibiotics_Duration_all",
        "PE",
        "Weight",
        "RDS",
        "GA_daily",
        "BPD",
        "NEC",
        "Gender",
        "Sepsis_style",
        "GDM")

nid = c("Age_DOL",
        "Feeding",
        "Antibiotics_Duration_all",
        "Delivery",
        "Weight",
        "GDM",
        "PE",
        "GA_daily",
        "Gender",
        "Intrauterine_infection")


ad1 = read.csv("./cohort1_adonis_clin_data_R2_result.csv")

adonis_data  = ad1 %>% filter(X %in% nid)

adonis_data$R2_2 = round(adonis_data$R2, digits = 4)

adonis_data$R2_2 = paste0(adonis_data$R2_2*100,"%")

adonis_data %>% arrange(R2) -> adonis_data
colnames(adonis_data)[1] = "feature"
adonis_data$feature = factor(adonis_data$feature , levels = unique(adonis_data$feature))

adonis_data %>%
  mutate(pif = case_when(
    `Pr..F.`>0.05 ~ "ns",
    `Pr..F.`<0.001 ~ "***",
    `Pr..F.`<0.01 ~ "**",
    `Pr..F.`<= 0.05 ~ "*",
  )) %>%
  mutate(test = case_when(
    `Pr..F.`>0.05 ~ "ns",
    TRUE ~ "signif"
  )) -> adonis_data1

adonis_data1 = adonis_data1 %>% 
  add_row(feature = "Feeding",R2 = 0.0001,test = "ns") %>%
  arrange(R2)
adonis_data1$feature = factor(adonis_data1$feature,levels = adonis_data1$feature)
ggplot(adonis_data1 , aes(x = R2, y = feature))+
  geom_col(aes(fill = test))+
  geom_text(aes(label = pif,hjust =0),size = 3,color = "red")+
  scale_fill_manual(values = c("ns" = "#e5e5d8", "signif" = "#d4af37")) +
  coord_cartesian(xlim = c(0,0.08))+
  labs(y = "CALM2005_A")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y  = element_blank()) -> p1

# c2 的绘图

ad2 = read.csv("./cohort2_adonis_clin_data_R2_result.csv")

nid = c("Age_DOL",
        "Feeding",
        "Antibiotics_Duration_all",
        "Natural_delivery",
        "Weight",
        "GDM",
        "PE",
        "GA_daily",
        "Gender",
        "Intrauterine_infection")


adonis_data  = ad2 %>% filter(X %in% nid)

adonis_data$R2_2 = round(adonis_data$R2, digits = 4)

adonis_data$R2_2 = paste0(adonis_data$R2_2*100,"%")

adonis_data %>% arrange(R2) -> adonis_data
colnames(adonis_data)[1] = "feature"
adonis_data$feature = factor(adonis_data$feature , levels = unique(adonis_data$feature))

adonis_data %>%
  mutate(pif = case_when(
    `Pr..F.`>0.05 ~ "ns",
    `Pr..F.`<0.001 ~ "***",
    `Pr..F.`<0.01 ~ "**",
    `Pr..F.`<= 0.05 ~ "*",
  )) %>%
  mutate(test = case_when(
    `Pr..F.`>0.05 ~ "ns",
    TRUE ~ "signif"
  )) -> adonis_data2

ggplot(adonis_data2 , aes(x = R2, y = feature))+
  geom_col(aes(fill = test))+
  geom_text(aes(label = pif,hjust =0),size = 3,color = "red")+
  scale_fill_manual(values = c("ns" = "#e5e5d8", "signif" = "#00d4c6")) +
  coord_cartesian(xlim = c(0,0.08))+
  labs(y = "CALM2005_B")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y  = element_blank()) -> p2
require(patchwork)
p2/p1/p3

ggsave("c1_c2_adonis_clin_data_R2_result.pdf", width = 6, height = 8)

ggplot(adonis_data2 , aes(x = R2, y = feature))+
  geom_col(aes(fill = test))+
  geom_text(aes(label = pif,hjust =0),size = 5,color = "red")+
  scale_fill_manual(values = c("ns" = "#e5e5d8", "signif" = "#85dab1")) +
  coord_cartesian(xlim = c(0,0.08))+
  labs(y = "CALM2005_B")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y  = element_blank())

ggsave("c2_adonis_clin_data_R2_result.pdf", width = 5, height = 3)
ggsave("c2_adonis_clin_data_R2_result.pdf", width = 5, height = 3)
