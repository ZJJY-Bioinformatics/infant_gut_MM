xbox::chdir("../9.slow_fast_sepsis_ab")
rm(list = ls())
library(tidyverse)
library(vegan)

all_data =  read.csv("../1.data_tidy/all_cohort.sum_7cohort_final.csv",row.names = 1)
all_data = all_data %>% filter(Age_DOL <=  42 )

sf_meta = read.csv("../4.bac_slow_fast////Sepsis_add_trend_cluster.meta.csv",row.names = 1)
sf_meta = sf_meta %>% select(Sample_ID,type,response)

all_data = all_data %>% left_join(sf_meta,by = c("Sample_ID"))

bac_index = read.csv("../11.bac_index/bac_index_risk.csv",row.names = 1)

all_data  = all_data %>% left_join(bac_index %>% select(SampleID,score,risk), by = c("Sample_ID" = "SampleID"))

# a value 

normal_data =  read.csv("../1.data_tidy/vmt_infant_stools_A.sum.csv",row.names = 1)
normal_data$Age_DOL= as.numeric(sub(".*\\.BSD","",normal_data$SampleID))

normal_data = normal_data %>%
  rename("Sample_ID" = "SampleID")

normal_data$Sepsis_style = "Control"
normal_data$NEC = "Control"
normal_data$type = "Control"
normal_data$response = 0
normal_data$score = 0
normal_data$risk = "Control"

normal_data1 = normal_data[,colnames(all_data)]

all_data = bind_rows(all_data,normal_data1)

all_data %>% mutate(
  Age_Group = case_when(
    Age_DOL <= 3 ~ "0d-3d",
    Age_DOL <= 7 ~ "4d-7d",
    Age_DOL <= 14 ~ "8d-14d",
    Age_DOL <= 21 ~ "15d-21d",
    Age_DOL <= 28 ~ "22d-28d",
    Age_DOL <= 35 ~ "29d-35d",
    Age_DOL >= 36 ~ "36d-42d",
  )
) -> all_data

write.csv(all_data,"all_cohort.sum_7cohort_final_Add_VMT.csv")

source("../addin.r")
count_data = round(all_data[,2:48]*10000)
rownames(count_data) = all_data$Sample_ID
alpha <- alpha_diversity(count_data)
all_data = cbind(all_data,alpha)
meio.mdf <- as.matrix.data.frame(all_data[,2:48])
rownames(meio.mdf) <- all_data$Sample_ID
meio.mdf[meio.mdf < 0.00001] = 0.00001

# BC----------
require(reshape2)
meio.bray <- vegdist(meio.mdf, method = "bray")
meio_bray_b = melt(as.matrix(meio.bray))
meio_bray_b =  meio_bray_b %>% filter(Var1!=Var2)
meio_bray_b = meio_bray_b %>% left_join(all_data %>% select("Sample_ID","type","Age_Group"),by = c("Var1" = "Sample_ID"))
meio_bray_b = meio_bray_b %>% left_join(all_data %>% select("Sample_ID","type","Age_Group"),by = c("Var2" = "Sample_ID"))

meio_bray_b %>% 
  filter(Var1 != Var2) %>%
  filter(Age_Group.x == Age_Group.y) %>%
  filter(type.x == "Control"& type.y %in% c("fast","slow")) %>%
  filter(!is.na(value)) %>%
  select(type.y,value,Age_Group.y) %>%
  rename("target"="type.y","dis" = "value","time"="Age_Group.y") %>%
  mutate(method = "Bray–Curtis") -> mdplot_b


# kulczynski
meio.bray_k <- vegdist(meio.mdf, method = "kulczynski")
meio_bray_k = melt(as.matrix(meio.bray_k))
meio_bray_k =  meio_bray_k %>% filter(Var1!=Var2)
meio_bray_k = meio_bray_k %>% left_join(all_data %>% select("Sample_ID","type","Age_Group"),by = c("Var1" = "Sample_ID"))
meio_bray_k = meio_bray_k %>% left_join(all_data %>% select("Sample_ID","type","Age_Group"),by = c("Var2" = "Sample_ID"))

meio_bray_k %>% 
  filter(Var1 != Var2) %>%
  filter(Age_Group.x == Age_Group.y) %>%
  filter(type.x == "Control"& type.y %in% c("fast","slow")) %>%
  filter(!is.na(value)) %>%
  select(type.y,value,Age_Group.y) %>%
  rename("target"="type.y","dis" = "value","time"="Age_Group.y") %>%
  mutate(method = "kulczynski") -> mdplot_k

# weight unifrac--------
library(phyloseq)
library(reshape2)
physeq = readRDS("../1.data_tidy/all_cohort.sum_7cohort_final_phylose.rds")
W.unifrac = readRDS("weight.unifrac_distance.rds")
meio_bray_wu = melt(as.matrix(W.unifrac))
meio_bray_wu =  meio_bray_wu %>% filter(Var1!=Var2)
meio_bray_wu = meio_bray_wu %>% left_join(all_data %>% select("Sample_ID","type","Age_Group"),by = c("Var1" = "Sample_ID"))
meio_bray_wu = meio_bray_wu %>% left_join(all_data %>% select("Sample_ID","type","Age_Group"),by = c("Var2" = "Sample_ID"))

meio_bray_wu %>% 
  filter(Var1 != Var2) %>%
  filter(Age_Group.x == Age_Group.y) %>%
  filter(type.x == "Control"& type.y %in% c("fast","slow")) %>%
  select(type.y,value,Age_Group.y) %>%
  rename("target"="type.y","dis" = "value","time"="Age_Group.y") %>%
  mutate(method = "weight Unifrac") -> mdplot_w

# unweight unifrac------
U.unifrac = readRDS("unweight.unifrac_distance.rds")
meio_bray_uu = melt(as.matrix(U.unifrac))
meio_bray_uu =  meio_bray_uu %>% filter(Var1!=Var2)
meio_bray_uu = meio_bray_uu %>% left_join(all_data %>% select("Sample_ID","type","Age_Group"),by = c("Var1" = "Sample_ID"))
meio_bray_uu = meio_bray_uu %>% left_join(all_data %>% select("Sample_ID","type","Age_Group"),by = c("Var2" = "Sample_ID"))

meio_bray_uu %>% 
  filter(Var1 != Var2) %>%
  filter(Age_Group.x == Age_Group.y) %>%
  filter(type.x == "Control"& type.y %in% c("fast","slow")) %>%
  filter(!is.na(value)) %>%
  select(type.y,value,Age_Group.y) %>%
  rename("target"="type.y","dis" = "value","time"="Age_Group.y") %>%
  mutate(method = "Unweight Unifrac") -> mdplot_u

mdplot  = bind_rows(mdplot_b,mdplot_k,mdplot_w,mdplot_u)

unique(mdplot$time)
mdplot$time = factor(mdplot$time,levels = c("0d-3d","4d-7d","29d-35d","36d-42d"))

# combine plot 
library(ggsignif)
colpal = readRDS("../colpal.rds")
require(gghalves)
ggplot() +
  geom_half_violin(
    data =mdplot, aes(x = time, y = 1-dis, fill = target,split = target),
    position = "identity",width = 0.8
  )+
  geom_boxplot(data =mdplot, aes(x = time, y = 1-dis, color = target),
               width = 0.25,fill = "white",outlier.shape = NA,coef = 0)+
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             label = "p.signif",
                             size = 3,
                             label.y = 0.1) +
  labs(x = "Method", y = "1-Distance")+
  scale_fill_manual(values = c("#ff8bbf","#a0a5e1"))+
  scale_color_manual(values = c("#ff8bbf","#a0a5e1"))+
  facet_wrap(".~method",ncol = 1,scales = "free")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("beta_distance_4method_bk.pdf",width = 5,height = 8)

colnames(mdplot)

median_values <- mdplot %>%
  group_by(time,method) %>%
  summarise(median_dis = median(1 - dis,na.rm = T))


ggplot(mdplot,aes(x = time,y = 1-dis))+
  geom_boxplot(aes(fill = target),width = 0.25,outlier.shape = NA,coef = 0)+
  geom_line(data = median_values, aes(x = factor(time), y = median_dis, group = 1), color = "blue", linetype = "dashed") +
  geom_point(data = median_values, aes(x = factor(time), y = median_dis), color = "#ff8bbf")+
  scale_fill_manual(values = c("#FDCEBA","#FBF6B5","#F2AECB","#49c8cb"))+
  ggpubr::stat_compare_means(
    comparisons = list(c("0d-3d","4d-7d"),
                       c("4d-7d","29d-35d"),
                       c("29d-35d","36d-42d")),
    method = "wilcox.test", 
    label = "p.value",
    size = 3,
    label.y = 0.1)+
  facet_wrap(".~method",ncol = 1,scales = "free")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("beta_distance_4method_group.pdf",width = 5,height = 8)

ggplot(mdplot,aes(x = time,y = 1-dis))+
  geom_boxplot(aes(fill = time),width = 0.25,outlier.shape = NA,coef = 0)+
  geom_line(data = median_values, aes(x = factor(time), y = median_dis, group = 1), color = "blue", linetype = "dashed") +
  geom_point(data = median_values, aes(x = factor(time), y = median_dis), color = "#ff8bbf")+
  scale_fill_manual(values = c("#FDCEBA","#FBF6B5","#F2AECB","#49c8cb"))+
  ggpubr::stat_compare_means(
    comparisons = list(c("0d-3d","4d-7d"),
                       c("4d-7d","29d-35d"),
                       c("29d-35d","36d-42d")),
    method = "wilcox.test", 
     label = "p.value",
     size = 3,
     label.y = 0.1)+
  facet_wrap(".~method",ncol = 1,scales = "free")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("beta_distance_4method_nogroup.pdf",width = 5,height = 8)

ggplot(mdplot, aes(x = time, y = 1-dis, fill = target)) +
  geom_boxplot(outlier.shape = NA,width = 0.5) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3,label.y = 0.5)+
  facet_wrap(".~method",ncol = 1)

ggsave("beta_distance_4method_add_signif.pdf",width = 5,height = 6)
saveRDS(mdplot,"beta_distance_4method.rds")

ggplot() +
  geom_half_violin(
    data =mdplot %>% filter(time%in% c("29d-35d","36d-42d")), 
    aes(x = method, y = 1-dis, fill = target,split = target),
    position = "identity",width = 0.8
  )+
  geom_boxplot(data =mdplot %>% filter(time%in% c("29d-35d","36d-42d")), 
               aes(x = method, y = 1-dis, color = target),
               width = 0.25,fill = "white",outlier.shape = NA,coef = 0)+
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             label = "p.signif",
                             size = 3,
                             label.y = 0.1) +
  labs(x = "Method", y = "1-Distance")+
  scale_fill_manual(values = c("#ff8bbf","#a0a5e1"))+
  scale_color_manual(values = c("#ff8bbf","#a0a5e1"))+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("beta_distance_4method_merge.pdf",width = 5,height = 3)


ggplot(data =mdplot %>% filter(time%in% c("29d-35d","36d-42d")), 
       aes(x = method, y = 1-dis, color = target)) +
  geom_boxplot(outlier.shape = NA,width = 0.5) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p", size = 3,label.y = 0.5)

ggsave("beta_distance_4method_merge_add_signif.pdf",width = 5,height = 3)
