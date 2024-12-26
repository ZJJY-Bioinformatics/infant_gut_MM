
xbox::chdir("./3.bac_exp_profile")
rm(list = ls())
library(tidyverse)
library(lmerTest)
library(broom.mixed)

ndata = read.csv("../1.data_tidy/all_cohort.sum_7cohort_final.csv",row.names = 1)
data = ndata %>% filter(Age_DOL <= 42)

fastslow = read_csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv")

data = data %>% 
  inner_join(fastslow %>%
               select(Sample_ID,type),
             by = c("Sample_ID" = "Sample_ID"))


write.csv(data,"plot_data_c1c2.csv")

sid = readRDS("../12.shap_value/top10.rds")

sdata = data %>% select(c("Sample_ID",sid,"Age_DOL","type"))

sdata_line = sdata %>% pivot_longer(sid,names_to = "bac",values_to = "exp")
sdata_line$bac = factor(sdata_line$bac,levels = sid)

lmm_table = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ type * Age_DOL + (1|PatientID), data = data2)
  res = tidy(fit) %>%
    filter(term %in% c("typeslow", "Age_DOL", "typeslow:Age_DOL")) %>%
    mutate(y = {{taxa}}) %>%
    select(y, everything())
  return(res)
}

all = map_dfr(sid, ~lmm_table(.x))

all %>%
  mutate(FDR = p.adjust(p.value, method = "fdr")) -> all_fdr
write_tsv(all_fdr, "LMM-risk + Age_DOL_adjust_confounder_result_FDR.tsv")


ggplot(sdata_line %>% filter(!is.na(type)), aes(Age_DOL, y = exp, color = type, fill = type)) + 
  geom_smooth(alpha = 0.2,method = "loess", span = 0.95) + 
  scale_color_manual(values = c("#ff8bbf","#a0a5e1")) + 
  scale_fill_manual(values =  c("#ff8bbf","#a0a5e1")) + 
  facet_wrap(".~bac",scales = "free",ncol = 5)+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Day", y = "Relative abundance")

ggsave("top9_shape_bac_smoothplot_all_type.pdf",width = 12,height = 4)


library(dplyr)
library(zoo)
library(entropy)

sdata_line %>% filter(!is.na(type)) %>%
  group_by(type,Age_DOL) %>%
  summarise(entropy_value = entropy(exp)) %>%
  ungroup() %>%
  rstatix::pairwise_wilcox_test(entropy_value ~ type) -> entropy_diff

write.csv(entropy_diff,"entropy_diff.csv")

sdata$Others = 1-apply(sdata[,2:11],1,sum)

# -----------
sdata = data  %>% filter(!is.na(type)) %>% select(c("Sample_ID",sid,"Age_DOL","type"))
sdata %>% mutate(
  Age_Group = case_when(
    Age_DOL <= 7 ~ "0d-7d",
    Age_DOL <= 14 ~ "8d-14d",
    Age_DOL <= 21 ~ "15d-21d",
    Age_DOL <= 28 ~ "22d-28d",
    Age_DOL >= 29 ~ "29d-42d",
  )
) -> sdata

sdata$Others = 1-apply(sdata[,2:11],1,sum)

sdata_bar = 
  sdata %>% pivot_longer(-c(Sample_ID,Age_Group,type,Age_DOL),
                         names_to = "bac",values_to = "exp") %>%
  group_by(Age_Group,type,bac) %>%
  summarise(exp = mean(exp))

sdata_bar = sdata_bar %>% 
  ungroup() %>%
  add_row(Age_Group = "F1",type = "fast",bac = "Others",exp = 1) %>%
  add_row(Age_Group = "F1",type = "slow",bac = "Others",exp = 1) %>%
  add_row(Age_Group = "F2",type = "fast",bac = "Others",exp = 1) %>%
  add_row(Age_Group = "F2",type = "slow",bac = "Others",exp = 1) 

sdata_bar$bac = factor(sdata_bar$bac,levels = rev(levels(sid)))
sdata_bar$Age_Group = factor(sdata_bar$Age_Group,levels = c("F1","0d-7d", "8d-14d","15d-21d",
                                                            "22d-28d","29d-42d","F2"))


colpal = readRDS("../colpal.rds")
colpal = c(colpal,c("Lactococcus" = "#FBF6B5",
                    "Acinetobacter" = "#F4C573",
                    "Lactobacillus" ="#80B1D3"))

ggplot(sdata_bar,aes(x = Age_Group,y  = exp, fill =  bac))+
  geom_bar(stat = "identity", position = "fill",width = 0.94) +
  scale_fill_manual(values =colpal) +
  facet_wrap(".~type")+
  labs(y = "relative abundance", x = "Time (Days)")+
  coord_polar()+
  ggthemes::theme_clean()
ggsave("AgeGroup_top10_bac_polar_42_semi_all_type.pdf",width = 15,height = 6)


# 画山峦图-------
sdata %>% pivot_longer(-c(Sample_ID,Age_Group,type,Age_DOL),
                       names_to = "bac",values_to = "exp") %>% 
  rstatix::group_by(bac) %>%
  rstatix::wilcox_test(exp ~ type) %>%
  mutate(p = round(p,3))-> df_temp

sdata_bar1 = sdata %>% pivot_longer(-c(Sample_ID,Age_Group,type,Age_DOL),
                                    names_to = "bac",values_to = "exp")

sdata_bar1$bac = factor(sdata_bar1$bac,levels = levels(sid))

require(ggridges)
ggplot(sdata_bar1 %>% filter(bac != "Others"),aes(y = bac)) +
  geom_density_ridges(
    aes(x = exp, fill =type), 
    alpha = .6, color = "white", from = 0, to = 1,scale = 2
  ) +
  labs(
    x = "relative abundance"
  ) +
  geom_text(data = df_temp%>% filter(bac != "Others"), aes(y = bac ,x= 0.8,label = paste("P =",p)),vjust  = -1)+
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("#1B8574","#bf48a8")) + 
  scale_fill_manual(values =  c("#1B8574","#bf48a8")) + 
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("ggridge_top10_all_type.pdf",width = 6,height = 4)