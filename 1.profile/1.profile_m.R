require(patchwork)
library(vegan)
require(tidyverse)
xbox::chdir("./1.profile")
rm(list = ls())

# read data------------------------
merge_sum = read.csv("../1.data_tidy/all_cohort_sample_num_8cohort_0905.csv",row.names = 1)
# merge_sum = read.csv("../1.data_tidy/all_cohort.sum_7cohort_final.csv",row.names = 1)
# main----------
xbox::chdir("final")

write.csv(data.frame(table(merge_sum$cohort)),"cohort_sample_num.csv")
merge_sum = merge_sum %>% filter(Age_DOL <= 42)
write.csv(data.frame(table(merge_sum$cohort)),"cohort_sample_num_42d.csv")
m = merge_sum %>% dplyr::select(PatientID,cohort) %>% distinct()
write.csv(data.frame(table(m$cohort)),"cohort_patient_num.csv")
merge_sum$PatientID = as.character(merge_sum$PatientID)

# add Dol group----
merge_sum %>% mutate(
  Age_Group = case_when(
    Age_DOL <= 3 ~ "0d-3d",
    Age_DOL <= 7 ~ "4d-7d",
    Age_DOL <= 14 ~ "8d-14d",
    Age_DOL <= 21 ~ "15d-21d",
    Age_DOL <= 28 ~ "22d-28d",
    Age_DOL >= 29 ~ "29d-42d",
  )
) -> merge_sum

# add 3 days a period
merge_sum$Day_Group <- cut(merge_sum$Age_DOL, 
                           breaks = seq(0, 105, by = 3), 
                           labels = FALSE, include.lowest = TRUE)

write.csv(merge_sum, "1.profile_plot_data.csv")

# 物种聚类
top10 <- names(rev(sort(apply(merge_sum[,2:48],2,sum))))[1:10]
meta_name = colnames(merge_sum)[49:ncol(merge_sum)]
top10_data = merge_sum[,c("Sample_ID",top10,meta_name)]

top10_data$Other = 1-apply(top10_data[2:11],1, sum) # 有些加起来大于1？
top10_data$Other[top10_data$Other < 0] = 0

colnames(top10_data)[ncol(top10_data)] = "Others"

# 展开数据
dodge_pdata = top10_data %>% 
  pivot_longer(-c(Sample_ID,Age_DOL,cohort,Age_Group,
                  PatientID,Sepsis_style,NEC,
                  Day_Group), values_to = "Abundance", names_to = "Bac")

# 绘制堆叠图
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

# 赋值颜色
palette =  readRDS("../../colpal.rds")

#saveRDS(c(palette,c("Robert_2024_M" = "#E0B794")),"../../colpal.rds")

# 均值汇总----
dodge_pdata_p <- dodge_pdata %>%
  group_by(Age_DOL, Bac) %>%
  summarise(Abundance = mean(Abundance))
# 汇总的结果
ggplot(dodge_pdata_p %>% filter(Age_DOL <= 42 & Age_DOL != 0) ,aes(x = Age_DOL,y  = Abundance, fill =  Bac))+
  geom_bar(stat = "identity", position = "fill")+
  scale_fill_manual(values = palette) +
  labs(y = "relative abundance", x = "Time (Days)")+
  coord_cartesian(xlim = c(0,42)) +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42),
                     labels = c(1,7,14,21,28,35,42))+
  ggthemes::theme_clean()+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("all_cohort_top10_bac_lineplot_42.pdf",width = 8,height = 4)


dodge_pdata_pc1c2 <- dodge_pdata %>%
  filter(grepl("CALM",cohort)) %>%
  group_by(Age_DOL, Bac) %>%
  summarise(Abundance = mean(Abundance))

ggplot(dodge_pdata_pc1c2 %>% 
         filter(Age_DOL <= 42 & Age_DOL != 0),aes(x = Age_DOL,y  = Abundance, fill =  Bac))+
  geom_bar(stat = "identity", position = "fill")+
  scale_fill_manual(values = palette) +
  labs(y = "relative abundance", x = "Time (Days)")+
  coord_cartesian(xlim = c(0,42)) +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42),
                     labels = c(1,7,14,21,28,35,42))+
  ggthemes::theme_clean()+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("c1c2_cohort_top10_bac_lineplot_42.pdf",width = 8,height = 4)

# 绘制菌群的表达趋势--------------
mod_top10_data = top10_data
dodge_pdata = mod_top10_data %>% 
  pivot_longer(-c(Sample_ID,Age_DOL,cohort,Age_Group,
                  PatientID,Sepsis_style,NEC,
                  Day_Group), values_to = "Abundance", names_to = "Bac")
top10_data %>%  
  arrange(Age_DOL) -> temp1

dodge_pdata$Sample_ID = factor(dodge_pdata$Sample_ID ,levels = temp1$Sample_ID)

dodge_pdata$Bac = factor(dodge_pdata$Bac, levels = c(top10))
dodge_pdata$cohort = factor(dodge_pdata$cohort,levels = cohorts)

ggplot(dodge_pdata %>% filter(Age_DOL <= 42) %>% filter(Bac != "Others"), aes(x = Age_DOL,y  = Abundance,color =  Bac,fill = Bac))+
  geom_smooth(method = "loess", span = 0.95)+
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(y = "Relative abundance", x = "Age DOL")+
  facet_wrap("cohort~Bac",ncol = 10,scales= "free_y")+
  coord_cartesian(xlim = c(0,43))+
  ggthemes::theme_few()+
  theme(axis.text = element_text(size = 7))+
  theme(strip.text = element_blank())

ggsave("top10_bac_smoothplot_42DOL_bac_color_adiversity.pdf",width = 13,height = 5.5)

# 计算P值

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
  group_by(Bac) %>%  # 按菌种分组
  do({
    icc_data <- column_to_rownames(.,"Age_DOL")  # 选择队列列
    icc_value <- ICC(icc_data[,-1])$results  # 计算 ICC
    tibble(
      Bac = unique(.$Bac),
      ICC_value = icc_value$ICC[6],  # 提取 ICC(1,1) 的值
      p_value = icc_value$p[6]  # 提取 ICC(1,1) 对应的 P 值
    )
  }) %>% mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
  ))-> p_values

write.csv(p_values,"top10_bac_smoothplot_ICC_p.csv")

# a diversity
source("../../addin.r")
file_up<-list.files("../../1.data_tidy/",pattern = "*.sum.csv")

for (i in 1:length(file_up)) {
  a<-read.csv(file.path("../../1.data_tidy/",file_up[i]))
  a = a[,-1]
  colnames(a)[1] = "Sample_ID"
  nn = unique(a$cohort)[1]
  assign(nn,a)
  rm(a)
}

colnames(CALM2005_1_A) = sub("\\[","_",colnames(CALM2005_1_A))
colnames(CALM2005_1_A) = sub("\\]","_",colnames(CALM2005_1_A))

colnames(Rao_2021_A) = sub("\\[","_",colnames(Rao_2021_A))
colnames(Rao_2021_A) = sub("\\]","_",colnames(Rao_2021_A))

colnames(CALM2005_2_A) = sub("\\[","_",colnames(CALM2005_2_A))
colnames(CALM2005_2_A) = sub("\\]","_",colnames(CALM2005_2_A))

colnames(Olm_2019_M) = sub("\\[","_",colnames(Olm_2019_M))
colnames(Olm_2019_M) = sub("\\]","_",colnames(Olm_2019_M))

colnames(Lauren_2022_M) = sub("\\[","_",colnames(Lauren_2022_M))
colnames(Lauren_2022_M) = sub("\\]","_",colnames(Lauren_2022_M))

colnames(Raveh_2015_M) = sub("\\[","_",colnames(Raveh_2015_M))
colnames(Raveh_2015_M) = sub("\\]","_",colnames(Raveh_2015_M))

colnames(Brooks_2017_M) = sub("\\[","_",colnames(Brooks_2017_M))
colnames(Brooks_2017_M) = sub("\\]","_",colnames(Brooks_2017_M))

colnames(Robert_2024_M) = sub("\\[","_",colnames(Robert_2024_M))
colnames(Robert_2024_M) = sub("\\]","_",colnames(Robert_2024_M))

mnid = colnames(merge_sum)[2:48]


data_exp=CALM2005_1_A[,1:450]
count_data = round(data_exp[,-1]*10000)
rownames(count_data) = data_exp$Sample_ID
alpha <- alpha_diversity(count_data)
alpha %>% rownames_to_column("Sample_ID") -> plot_a_data_temp1
top_v = data.frame(apply(count_data,2,sum)/sum(apply(count_data,2,sum)))[mnid,]
data.frame(bac = mnid,ratio = top_v) -> top10_ratio_1
top10_ratio_1$cohort = "CALM2005_1_A"

data_exp=CALM2005_2_A[,1:497]
count_data = round(data_exp[,-1]*10000)
rownames(count_data) = data_exp$Sample_ID
alpha <- alpha_diversity(count_data)
alpha %>% rownames_to_column("Sample_ID") -> plot_a_data_temp2
top_v = data.frame(apply(count_data,2,sum)/sum(apply(count_data,2,sum)))[mnid,]
data.frame(bac = mnid,ratio = top_v) -> top10_ratio_2
top10_ratio_2$cohort = "CALM2005_2_A"

data_exp=Rao_2021_A[,1:237]
count_data = round(data_exp[,-1]*10000)
rownames(count_data) = data_exp$Sample_ID
alpha <- alpha_diversity(count_data)
alpha %>% rownames_to_column("Sample_ID") -> plot_a_data_temp3
top_v = data.frame(apply(count_data,2,sum)/sum(apply(count_data,2,sum)))[mnid,]
data.frame(bac = mnid,ratio = top_v) -> top10_ratio_3
top10_ratio_3$cohort = "Rao_2021_A"

data_exp=Lauren_2022_M[,1:134]
count_data = round(data_exp[,-1]*10000)
rownames(count_data) = data_exp$Sample_ID
alpha <- alpha_diversity(count_data)
alpha %>% rownames_to_column("Sample_ID") -> plot_a_data_temp4
top_v = data.frame(apply(count_data,2,sum)/sum(apply(count_data,2,sum)))[mnid,]
data.frame(bac = mnid,ratio = top_v) -> top10_ratio_4
top10_ratio_4$cohort = "Lauren_2022_M"

data_exp=Brooks_2017_M[,1:102]
count_data = round(data_exp[,-1]*10000)
rownames(count_data) = data_exp$Sample_ID
alpha <- alpha_diversity(count_data)
alpha %>% rownames_to_column("Sample_ID") -> plot_a_data_temp5
top_v = data.frame(apply(count_data,2,sum)/sum(apply(count_data,2,sum)))[mnid,]
data.frame(bac = mnid,ratio = top_v) -> top10_ratio_5
top10_ratio_5$cohort = "Brooks_2017_M"

data_exp=Olm_2019_M[,1:102]
count_data = round(data_exp[,-1]*10000)
rownames(count_data) = data_exp$Sample_ID
alpha <- alpha_diversity(count_data)
alpha %>% rownames_to_column("Sample_ID") -> plot_a_data_temp6
top_v = data.frame(apply(count_data,2,sum)/sum(apply(count_data,2,sum)))[mnid,]
data.frame(bac = mnid,ratio = top_v) -> top10_ratio_6
top10_ratio_6$cohort = "Olm_2019_M"

data_exp=Raveh_2015_M[,1:102]
count_data = round(data_exp[,-1]*10000)
rownames(count_data) = data_exp$Sample_ID
alpha <- alpha_diversity(count_data)
alpha %>% rownames_to_column("Sample_ID") -> plot_a_data_temp7
top_v = data.frame(apply(count_data,2,sum)/sum(apply(count_data,2,sum)))[mnid,]
data.frame(bac = mnid,ratio = top_v) -> top10_ratio_7
top10_ratio_7$cohort = "Raveh_2015_M"

data_exp=Robert_2024_M[,1:132]
count_data = round(data_exp[,-1]*10000)
rownames(count_data) = data_exp$Sample_ID
alpha <- alpha_diversity(count_data)
alpha %>% rownames_to_column("Sample_ID") -> plot_a_data_temp8
top_v = data.frame(apply(count_data,2,sum)/sum(apply(count_data,2,sum)))[mnid,]
data.frame(bac = mnid,ratio = top_v) -> top10_ratio_8
top10_ratio_8$cohort = "Robert_2024_M"


top10_ratio = rbind(top10_ratio_1,top10_ratio_2,top10_ratio_3,top10_ratio_4,
                    top10_ratio_5,top10_ratio_6,top10_ratio_7,top10_ratio_8)

top10_ratio %>% 
  pivot_wider(cohort,values_from = ratio,names_from = bac) %>%
  write.csv("top10_ratio_cohort.csv")


a_data = rbind(plot_a_data_temp1,plot_a_data_temp2,plot_a_data_temp3,plot_a_data_temp4,
               plot_a_data_temp5,plot_a_data_temp6,plot_a_data_temp7,plot_a_data_temp8)
write.csv(a_data,"a_diversity_value.csv")

alpha = a_data %>% column_to_rownames("Sample_ID")
cbind(merge_sum[,c(1,49:55)],alpha[merge_sum$Sample_ID,]) -> plot_a_data


adata = plot_a_data %>% 
  hablar::retype() %>%
  pivot_longer(-c(Sample_ID,Age_DOL,cohort,Age_Group,
                  PatientID,Sepsis_style,NEC,
                  Day_Group), values_to = "Value", names_to = "Diversity")

adata$cohort = factor(adata$cohort,levels = cohorts)


ggplot(adata %>% filter(Age_DOL <= 42 ) %>% filter(Diversity != "goods_Coverage"),
       aes(x = Age_DOL, y = Value,color = cohort))+
  geom_smooth(method = "loess", span = 0.95, linewidth = 1,
              aes(fill = cohort),alpha = 0.2) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(y = "a Diversity", x = "Time (Days)")+
  facet_wrap(".~Diversity",scales = "free_y",nrow = 1)+
  ggthemes::theme_few()

ggsave("cohort_age_adiversity_smoothplot.pdf",width = 10,height = 3)


ggplot(adata %>% filter(Age_DOL <= 42 ) %>% filter(Diversity != "goods_Coverage"),
       aes(x = Age_DOL, y = Value))+
  geom_smooth(method = "loess",
              color = "#D95050",fill = "#F2AECB",se = T,
              span = 0.95, linewidth = 1,alpha = 0.2) +
  labs(y = "a Diversity", x = "Time (Days)")+
  facet_wrap(".~Diversity",scales = "free_y",nrow = 1)+
  ggthemes::theme_few()

ggsave("cohort_age_adiversity_smoothplot_combine.pdf",width = 10,height = 3)



adata %>% 
  select(Sample_ID,Age_DOL,Diversity,cohort,Value) %>%
  pivot_wider(id_cols = c(Diversity,Sample_ID,Age_DOL),names_from = cohort,values_from = Value) %>%
  select(-Sample_ID) %>%
  group_by(Age_DOL,Diversity) %>%
  mutate(across(.cols = everything(), .fns = ~sum(., na.rm = TRUE))) %>%
  group_by(Age_DOL) %>%
  slice_head(n = 6) -> re 

re %>%
  group_by(Diversity) %>%  # 按菌种分组
  do({
    icc_data <- column_to_rownames(.,"Age_DOL")  # 选择队列列
    icc_value <- ICC(icc_data[,-1])$results  # 计算 ICC
    tibble(
      Diversity = unique(.$Diversity),
      ICC_value = icc_value$ICC[6],  # 提取 ICC(1,1) 的值
      p_value = icc_value$p[6]  # 提取 ICC(1,1) 对应的 P 值
    )
  }) %>% mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
  ))-> p_values

write.csv(p_values,"a_diversity_ICC_p.csv")

# b diversity-----------
meio.mdf <- as.matrix.data.frame(merge_sum[,2:48])
rownames(meio.mdf) <- merge_sum$Sample_ID

meio.mdf[meio.mdf < 0.00001] = 0.00001
require(vegan)
meio.bray <- vegdist(meio.mdf, method = "bray")
pcoa.meio.bray <- cmdscale(meio.bray, k = 2, eig = T)

pcoa.meio.bray.plotting <- as.data.frame(pcoa.meio.bray$points)
colnames(pcoa.meio.bray.plotting) <- c("axis_1", "axis_2")

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
pc1_r = pcoa.meio.bray$eig[1]/(sum(pcoa.meio.bray$eig))
pc2_r = pcoa.meio.bray$eig[2]/(sum(pcoa.meio.bray$eig))
# Age_DOL
pcoa.meio.bray.plotting$Sample_ID <- merge_sum$Sample_ID
pcoa.meio.bray.plotting$Age_DOL <- merge_sum$Age_DOL
pcoa.meio.bray.plotting$Age_Group <- merge_sum$Age_Group
pcoa.meio.bray.plotting$cohort <- merge_sum$cohort
pcoa.meio.bray.plotting$NEC <- merge_sum$NEC
pcoa.meio.bray.plotting$Sepsis_style <- merge_sum$Sepsis_style
pcoa.meio.bray.plotting$Day_Group <- merge_sum$Day_Group
pcoa.meio.bray.plotting$Sample_ID <- merge_sum$Sample_ID

# 查找代表属
top10_p = top10_data

top10_p$rep_genus = colnames(top10_p)[c(2:5)][apply(top10_p[,top10[1:4]], 1, which.max)]
top10_p$Shannon = as.numeric(plot_a_data$Shannon)
top10_p$cohort = plot_a_data$cohort

top10_p = top10_p %>% mutate(rep_genus =
                               case_when(
                                 Shannon > 2.7 ~ "Others",
                                 TRUE ~ rep_genus
                               )) %>% mutate(
                                 rep_genus = replace_na(rep_genus,"Others")
                               )

# top10_p <- top10_p %>%
#   group_by(cohort) %>%
#   mutate(
#     threshold_05 = quantile(Shannon, probs = 0.9, na.rm = TRUE), # 每组的 25% 分位数
#     rep_genus = case_when(
#       Shannon >= threshold_05 ~ "Others", # 如果 Simpson 值小于等于分位数，标记为 Others
#       TRUE ~ rep_genus                   # 否则保留原来的值
#     )
#   ) %>%
#   ungroup() %>% # 取消分组
#   mutate(rep_genus = replace_na(rep_genus, "Others"))


table(top10_p$rep_genus)


pcoa.meio.bray.plotting$rep_genus <- top10_p$rep_genus
pcoa.meio.bray.plotting$PatientID <- top10_p$PatientID

write.csv(pcoa.meio.bray.plotting,"pcao_meio,bray.data.csv")
# anova p值检验

source("../../addin.r")
anova_data = data.frame()

for(i in c("Age_DOL","Age_Group","cohort","rep_genus")){
  anova_result = adonis3(as.formula(paste("meio.bray","~",i)),
                         data = pcoa.meio.bray.plotting,
                         permutations = 999,
                         na.action=na.omit)
  anova_temp = data.frame(factor = i , Pvalue = anova_result[1,5],R2 = round(anova_result[1,3], 3))
  anova_data = rbind(anova_data,anova_temp)
}
write.csv(anova_data,"bdiversity_adonis.csv")

write.csv(pcoa.meio.bray.plotting,"bdiversity_data.csv",row.names = T)

colpal = readRDS("../../colpal.rds")
# 直接分组-----------------
pcoa.meio.bray.plotting %>% group_by(rep_genus) %>%
  summarise(pc1 = mean(axis_1),
            pc2 = mean(axis_2),
            p90_age_dol = quantile(Age_DOL, na.rm = TRUE,0.1),
            p10_age_dol = quantile(Age_DOL, na.rm = TRUE,0.9),
            median_age_dol = median(Age_DOL, na.rm = TRUE),
            sd_age_dol = sd(Age_DOL, na.rm = TRUE), 
            se_age_dol = sd(Age_DOL, na.rm = TRUE) / sqrt(n())) -> g_p

ggplot(pcoa.meio.bray.plotting %>% filter(Age_DOL <= 42)) +
  geom_point(aes(x = axis_1, y = axis_2, colour = rep_genus),size = 3) +
  geom_point(data = g_p,aes(x =  pc1, y = pc2),size = 2)+
  geom_label(data = g_p,aes(x =  pc1, y = pc2,
                            label = paste0(round(median_age_dol,1),"(","SE = ",round(se_age_dol,2),")")))+
  ggthemes::theme_few() +
  scale_color_manual(values = colpal)+
  xlab(paste0("PCoA 1 (",round(pc1_r,3)*100,")%")) +
  ylab(paste0("PCoA 2 (",round(pc2_r,3)*100,")%")) +
  #annotate(geom = 'text', label = paste0('ANOVA P:',"0.001", "\n R2:",anova_data[4,3]*100,"%"), x = -Inf, y = -Inf, hjust = -0.1,vjust = -1)+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("rep_genus_pcoa.pdf",width = 8.5,height = 6)

# Age_Group 做
pcoa.meio.bray.plotting$Age_Group = factor(pcoa.meio.bray.plotting$Age_Group,
                                           levels = c("0d-3d" ,"4d-7d", "8d-14d","15d-21d","22d-28d","29d-42d"))
ggplot(pcoa.meio.bray.plotting %>% filter(Age_DOL <= 42), 
       aes(x = axis_1, y = axis_2, colour = Age_Group)) +
  geom_point(size = 1) +
  ggthemes::theme_few(base_size = 10) +
  scale_color_manual(values = c("#eef2ea","#c7eab8","#a9de92","#a3b899","#386823","#050a03"))+
  xlab(paste0("PCoA 1 (",round(pc1_r,3)*100,")%")) +
  ylab(paste0("PCoA 2 (",round(pc2_r,3)*100,")%")) +
  annotate(geom = 'text', label = paste0('ANOVA P:',"0.001", "\n R2:",anova_data[2,3]*100,"%"), x = -Inf, y = -Inf, hjust = -0.1,vjust = -1)+
  ggtitle(paste("Age_Group"))->p2
# cohort
pcoa.meio.bray.plotting$cohort = factor(pcoa.meio.bray.plotting$cohort,
                                        levels = cohorts)
ggplot(pcoa.meio.bray.plotting %>% filter(Age_DOL <= 42), 
       aes(x = axis_1, y = axis_2, colour = cohort),alpha = 0.2) +
  geom_point(size = 1) +
  ggthemes::theme_few(base_size = 10) +
  scale_color_manual(values = palette)+
  xlab(paste0("PCoA 1 (",round(pc1_r,3)*100,")%")) +
  ylab(paste0("PCoA 2 (",round(pc2_r,3)*100,")%")) +
  annotate(geom = 'text', label = paste0('ANOVA P:',"0.001", "\n R2:",anova_data[3,3]*100,"%"), x = -Inf, y = -Inf, hjust = -0.1,vjust = -1)+
  ggtitle(paste("Cohort")) -> p3

# # 构建拟合曲线
data_filtered <- pcoa.meio.bray.plotting %>%
  arrange(Age_DOL)

fit_data <- data_filtered %>%
  select(axis_1, axis_2) %>%
  as.matrix()
principal_curve_fit <- princurve::principal_curve(fit_data)
data_filtered$trajectory_axis_1 <- principal_curve_fit$s[, 1]
data_filtered$trajectory_axis_2 <- principal_curve_fit$s[, 2]

ggplot(pcoa.meio.bray.plotting %>% filter(Age_DOL <= 42)) +
  geom_point(aes(x = axis_1, y = axis_2, colour = Age_DOL),size = 1) +
  geom_smooth(data = data_filtered,
              aes(x = trajectory_axis_1, y = trajectory_axis_2),
              colour = "black", size = 3) +
  ggthemes::theme_few() +
  scale_color_gradient2(low = "#40e0d0",midpoint = 21,high = "#f52f04")+
  xlab(paste0("PCoA 1 (",round(pc1_r,3)*100,")%")) +
  ylab(paste0("PCoA 2 (",round(pc2_r,3)*100,")%")) +
  annotate(geom = 'text', label = paste0('ANOVA P:',"0.001", "\n R2:",anova_data[1,3]*100,"%"), x = -Inf, y = -Inf, hjust = -0.1,vjust = -1)+
  ggthemes::theme_few(base_size = 10) +
  theme(panel.grid.major.y = element_blank()) +
  ggtitle(paste("DOL"))-> p5

p5+p2+p3
ggsave("mutipl_pcoa_plot.pdf",width = 18,height = 4)

# time display

pcoa.meio.bray.plotting %>% group_by(rep_genus,cohort) %>%
  summarise(pc1 = mean(axis_1),
            pc2 = mean(axis_2),
            p90_age_dol = quantile(Age_DOL, na.rm = TRUE,0.25),
            p10_age_dol = quantile(Age_DOL, na.rm = TRUE,0.75),
            median_age_dol = median(Age_DOL, na.rm = TRUE),
            sd_age_dol = sd(Age_DOL, na.rm = TRUE), 
            se_age_dol = sd(Age_DOL, na.rm = TRUE) / sqrt(n())) -> g_p

g_p = g_p %>% arrange(cohort,rep_genus)
g_p$rep_genus = factor(g_p$rep_genus,levels = c("Klebsiella","Escherichia",
                                                "Enterococcus","Staphylococcus",
                                                "Others"))

ggplot(g_p, aes(y = rep_genus, x = median_age_dol,fill = cohort)) +
  geom_point( shape = 21 , size = 5) +  # 左端的原点
  scale_fill_manual(values = colpal) +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42),
                     labels = c(1,7,14,21,28,35,42))+
  ggthemes::theme_clean() +
  theme(panel.grid.major.y = element_blank())

ggsave("mutipl_pcoa_plot_line.pdf",width = 8,height = 4)

write.csv(g_p,"mutipl_pcoa_plot_line.csv")

pcoa.meio.bray.plotting$rep_genus = factor(pcoa.meio.bray.plotting$rep_genus,levels =c("Klebsiella","Escherichia",
                                                                                       "Enterococcus","Staphylococcus",
                                                                                       "Others"))
require(gghalves)

ggplot(pcoa.meio.bray.plotting %>% 
         filter(Age_DOL <= 42),aes(x = rep_genus, y = Age_DOL))+
  geom_boxplot(aes(fill = cohort),width = 0.4,outlier.shape = NA, coef = 0)+
  scale_fill_manual(values = colpal) +
  scale_color_manual(values = colpal) +
  scale_y_continuous(breaks = c(seq(1, 42, 7), 42),
                     labels = c(seq(1, 42, 7), 42)) +
  coord_flip()+
  ggthemes::theme_clean() +
  theme(panel.grid.major.y = element_blank())


ggsave("mutipl_pcoa_plot_boxplot.pdf",width = 6,height = 4)

