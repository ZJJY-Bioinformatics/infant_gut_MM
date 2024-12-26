xbox::chdir("./24.wet_valid")
require(tidyverse)
require(ggprism)
library(readxl)
rm(list = ls())
# f6a-------------
data <- read_excel("huidi.xlsx", sheet = "NOD2激活实验")
data = data[1,]
data = as.data.frame(t(data)) %>% rownames_to_column("sample") %>%
  filter(sample != "Group")
colnames(data) = c("group","od")
data$od = as.numeric(data$od)
data$group = factor(data$group,levels = data$group)
ggplot(data = data,aes(x = group, y = 1))+
  geom_tile(aes(fill  = od))+
  geom_text(aes(label = od))+
  scale_fill_gradient2(high = "#f222a9",low = "#0d98ba",mid = "white",midpoint = 0.1)+
  theme_classic()+
  labs(x = "Isolated from commercial probiotics\nIsolated from human feces",
       y = "NOD2 activation (OD value)",
       title = "E. Faecium activates NOD2 receptor with SagA")
ggsave("f6a.pdf",width = 10,height = 2)

# f6b-------------
rm(list = ls())
data <- read_excel("huidi.xlsx", sheet = "SagA")
data1 = as.data.frame(t(data)) %>% rownames_to_column("sample") %>%
  filter(sample != "Group")
colnames(data1) = c("group","saga")
data1$saga = as.numeric(data1$saga)
data1$group = factor(data1$group,levels = data1$group)
ggplot(data = data1,aes(x = group, y = 1))+
  geom_tile(aes(fill  = saga))+
  geom_text(aes(label = round(saga,2)))+
  scale_fill_gradient2(high = "#f222a9",low = "#0d98ba",mid = "white",midpoint = 30)+
  theme_classic()+
  labs(x = "Isolated from commercial probiotics\nIsolated from human feces",
       y = "SagA (ng/L)",
       title = "E. Faecium activates NOD2 receptor with SagA")
ggsave("f6b.pdf",width = 10,height = 2)

#f6c----------------
require(tidyverse)
library(survival)
library(ggsurvfit)
library(survminer)
require(survminer)
library(gridExtra)
rm(list = ls())
data <- read_excel("huidi.xlsx", sheet = "生存曲线")
data1 = data %>% mutate(ID = paste("mm",row_number()))

data1 %>% 
  pivot_longer(-c(Hours,ID),names_to = "group", values_to = "OS") %>%
  filter(!is.na(OS)) %>%
  rename("OS.time.y" = "Hours") -> sur_data

write.csv(sur_data,"f6c_table.csv")

#sur_data = read.csv("f6c_table.csv",row.names = 1)
#sur_data = bind_rows(sur_data,sur_data)

surv_data <- with(sur_data, Surv(OS.time.y, OS))
fit <- survfit2(surv_data ~ group, data = sur_data)

fit_data = tidy_survfit(fit)
colpal = readRDS("../colpal.rds")
c1 = c(
  "Ef1" = "#a0a5e1",
  "Ef1+GSK717" = "#fbc7cc",
  "Ef3" = "#A7D1CD",
  "Ef3+GSK717" = "#ff8bbf",
  "MDP" = "#49c8cb",
  "MDP+GSK717" = "#D95050",
  "PBS" = "orange"
)

# c(colpal[1:26],c1) -> colpal
# saveRDS(colpal,"../colpal.rds")

ggsurvplot(fit, data = sur_data,
           fun = "event",
           pval = TRUE, conf.int = FALSE,
           conf.int.style="step",
           surv.scale = "percent",
           palette = unname(c1),
           censor = F,
           censor.size = 10,
           ylim = c(0,0.70),
           risk.table = TRUE,
           risk.table.title = "Note the risk set sizes",
           risk.table.subtitle = "",
           risk.table.caption = "",
           risk.table.height = 0.45,
           break.time.by = 2) -> p1

pairwise_results <- pairwise_survdiff(Surv(OS.time.y, OS) ~ group, 
                                      data = sur_data, 
                                      p.adjust.method = "none")
write.csv(data.frame(pairwise_results[["p.value"]]),"f6c_pairwise_p.csv")
cox_model <- coxph(Surv(OS.time.y, OS) ~ group, data = sur_data)
model_summary <- summary(cox_model)
# 提取风险比
hazard_ratios <- exp(coef(cox_model))
# 提取风险比的置信区间
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

write_csv(results_df,"f6c.csv")

g1 <- ggplotGrob(p1$plot)
g1
ggsave("f6c.pdf",g1 ,width = 4, height = 3.5, units = "in")
p1
ggsave("f6c_table.pdf" ,width = 5, height = 2.5, units = "in")


ggforest(cox_model,data = as.data.frame(sur_data))

# f6d--------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "体重")

data1 = data %>% pivot_longer(-c(ID,group),values_to = "weight",names_to = "time")

err_data = data1 %>% select(group,time,weight)

library(rstatix)
pairwise_results <- err_data %>%
  filter(time == 12) %>% # 按时间点分组
  pairwise_wilcox_test(weight ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  filter(group2 =="PBS")


err_data_summary <- err_data %>%
  group_by(group,time) %>%
  summarise(
    mean_n = mean(weight, na.rm = TRUE),
    sem_n = sd(weight, na.rm = TRUE) / sqrt(n())
  )

pairwise_results = pairwise_results %>% left_join(err_data_summary %>% filter(time == 12) %>% select(group,mean_n),by=c("group1" = "group"))

err_data$time = as.numeric(err_data$time)
err_data_summary$time = as.numeric(err_data_summary$time)
ggplot() +
  geom_line(data = err_data_summary,
            aes(x = time, y = mean_n,color = group,group = group),
            linewidth = 1) +
  geom_point(data = err_data_summary,
             aes(x = time,y = mean_n,color = group), size = 3) +
  geom_errorbar(data = err_data_summary,
                aes(x = time,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), width = 0.1) +
  geom_text(data = pairwise_results,
            aes(x = 12.5,
                y = mean_n,
                label = stars))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  coord_cartesian(ylim = c(80,100))+
  theme_prism(base_size = 16,base_fontface = "plain",base_line_size = 1) +
  scale_y_continuous(
    breaks = seq(80, 100, 10),
    guide = "prism_offset"
  )+ 
  scale_x_continuous(
    breaks = c(0, 6, 12),
    labels = c("0h","6h" ,"12h"),
    guide = "prism_offset_minor",
  )+
  labs(y = "Body weight(% at hour 0)",
       x = "Hour")
ggsave("f6d.pdf",width = 6,height = 4)

#f3e-------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "Health score")

data1 = data %>% pivot_longer(-c(ID,group),values_to = "weight",names_to = "time")

err_data = data1 %>% select(group,time,weight)

library(rstatix)
pairwise_results <- err_data %>%
  filter(time == 12) %>% # 按时间点分组
  pairwise_wilcox_test(weight ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  filter(group2 =="PBS")


err_data_summary <- err_data %>%
  group_by(group,time) %>%
  summarise(
    mean_n = mean(weight, na.rm = TRUE),
    sem_n = sd(weight, na.rm = TRUE) / sqrt(n())
  )

pairwise_results = pairwise_results %>% left_join(err_data_summary %>% filter(time == 12) %>% select(group,mean_n),by=c("group1" = "group"))

err_data$time = as.numeric(err_data$time)
err_data_summary$time = as.numeric(err_data_summary$time)

ggplot() +
  geom_line(data = err_data_summary,
            aes(x = time, y = mean_n,color = group,group = group),
            linewidth = 1) +
  geom_point(data = err_data_summary,
             aes(x = time,y = mean_n,color = group), size = 3) +
  geom_errorbar(data = err_data_summary,
                aes(x = time,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), width = 0.1) +
  geom_text(data = pairwise_results,
            aes(x = 12.5,
                y = mean_n,
                label = stars))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 16,base_fontface = "plain",base_line_size = 1)+
  scale_y_continuous(
    breaks = seq(0, 5, 1),
    guide = "prism_offset"
  )+ 
  scale_x_continuous(
    breaks = c(0, 6, 12),
    labels = c("0h","6h" ,"12h"),
    guide = "prism_offset_minor",
  )+
  labs(y = "Sepsis health score",
       x = "Hour")
ggsave("f6e.pdf",width = 6,height = 4)

# f6f-1---------------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "肺部评分")

data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1

data1$group = factor(data1$group,levels = c("PBS","Ef1","Ef3","MDP","Ef1+GSK717","Ef3+GSK717","MDP+GSK717"))

library(rstatix)
pairwise_results <- data1 %>%
  pairwise_wilcox_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  filter(group2 =="PBS") %>%
  add_y_position()

err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,12))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 12, 2),
    guide = "prism_offset"
  )+ 
  labs(y = "Histopathology score",
       x = "",
       title = "Lung") -> p1
# f6f-2---------------------
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "肝脏评分")

data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1
data1$group = factor(data1$group,levels = c("PBS","Ef1","Ef3","MDP","Ef1+GSK717","Ef3+GSK717","MDP+GSK717"))

library(rstatix)
pairwise_results <- data1 %>%
  pairwise_wilcox_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  filter(group2 =="PBS") %>%
  add_y_position()

err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,5))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 5, 1),
    guide = "prism_offset"
  )+ 
  labs(y = "Histopathology score",
       x = "",
       title = "Liver") -> p2

# f6f-3---------------------
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "结肠评分")

data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1
data1$group = factor(data1$group,levels = c("PBS","Ef1","Ef3","MDP","Ef1+GSK717","Ef3+GSK717","MDP+GSK717"))

library(rstatix)
pairwise_results <- data1 %>%
  pairwise_wilcox_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  filter(group2 =="PBS") %>%
  add_y_position()

err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),
           color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,3))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 3, 1),
    guide = "prism_offset"
  )+ 
  labs(y = "Histopathology score",
       x = "",
       title = "Colon") -> p3 
require(patchwork)
p1+p2+p3+patchwork::plot_layout(ncol = 3,guides = "collect")+ plot_annotation(tag_levels = "a")

ggsave("f6f.pdf",width = 12.5,height = 4)


# f6g-------------

rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "IL6")

# colpal=setNames(colpal[1:7],colnames(data))
# saveRDS(colpal,"./colpal.rds")


data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1

data1$group = factor(data1$group,levels = c("PBS","Ef1","Ef3","MDP","Ef1+GSK717","Ef3+GSK717","MDP+GSK717"))



library(rstatix)
pairwise_results <- data1 %>%
  pairwise_wilcox_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  filter(group1 =="PBS") %>%
  add_y_position()

err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,1500))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 1500, 500),
    guide = "prism_offset"
  )+ 
  labs(y = "IL-6 (pg/mL)",
       x = "",
       title = "IL-6") -> p1

# f6f-2---------------------
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "IL-1B")

data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1
data1$group = factor(data1$group,levels = c("PBS","Ef1","Ef3","MDP","Ef1+GSK717","Ef3+GSK717","MDP+GSK717"))


library(rstatix)
pairwise_results <- data1 %>%
  pairwise_wilcox_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  filter(group2 =="PBS") %>%
  add_y_position()

err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,7500))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 7500, 2500),
    guide = "prism_offset"
  )+ 
  labs(y = "IL-1B (pg/mL)",
       x = "",
       title = "IL-1B") -> p2

# f6f-3---------------------
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "IFN")

data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1

data1$group = factor(data1$group,levels = c("PBS","Ef1","Ef3","MDP","Ef1+GSK717","Ef3+GSK717","MDP+GSK717"))

library(rstatix)
pairwise_results <- data1 %>%
  pairwise_wilcox_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  filter(group2 =="PBS") %>%
  add_y_position()

err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),
           color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,1000))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 1000, 200),
    guide = "prism_offset"
  )+ 
  labs(y = "IFN (pg/mL)",
       x = "",
       title = "IFN") -> p3 
# f6g-4---------------------
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "IL10")

data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1
data1$group = factor(data1$group,levels = c("PBS","Ef1","Ef3","MDP","Ef1+GSK717","Ef3+GSK717","MDP+GSK717"))

library(rstatix)
pairwise_results <- data1 %>%
  pairwise_wilcox_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  filter(group2 =="PBS") %>%
  add_y_position()

err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),
           color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,1000))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 1000, 200),
    guide = "prism_offset"
  )+ 
  labs(y = "IL10 (pg/mL)",
       x = "",
       title = "IL10") -> p4 

require(patchwork)
p1+p2+p3+p4+patchwork::plot_layout(ncol = 4,guides = "collect")+ plot_annotation(tag_levels = "a")

ggsave("f6g.pdf",width = 12.5,height = 4)

#f6h-----------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("xirao.xlsx", sheet = "nod2 activate")

data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1

library(rstatix)
pairwise_results <- data1 %>%
  pairwise_t_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  filter(group2 =="PBS") %>%
  add_y_position()

err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),
           color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,1.5))+
  scale_fill_manual(values = colpal[12:16])+
  scale_color_manual(values = colpal[12:16])+
  theme_prism(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 1.5,0.5),
    guide = "prism_offset"
  )+ 
  labs(y = "Relative absorbance",
       x = "",
       title = "")

ggsave("f6h.pdf",width = 5,height = 4)

#f6i-------------------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("xirao.xlsx", sheet = "图2WT体重")

data1 = data %>% pivot_longer(-c(ID,group),values_to = "weight",names_to = "time")

err_data = data1 %>% select(group,time,weight)

library(rstatix)
pairwise_results <- err_data %>%
  filter(time == 12) %>% # 按时间点分组
  pairwise_t_test(weight ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) 

err_data_summary <- err_data %>%
  group_by(group,time) %>%
  summarise(
    mean_n = mean(weight, na.rm = TRUE),
    sem_n = sd(weight, na.rm = TRUE) / sqrt(n())
  )

pairwise_results = pairwise_results %>% left_join(err_data_summary %>% filter(time == 12) %>% select(group,mean_n),by=c("group1" = "group"))

err_data$time = as.numeric(err_data$time)
err_data_summary$time = as.numeric(err_data_summary$time)
ggplot() +
  geom_line(data = err_data_summary,
            aes(x = time, y = mean_n,color = group,group = group),
            linewidth = 1) +
  geom_point(data = err_data_summary,
             aes(x = time,y = mean_n,color = group), size = 3) +
  geom_errorbar(data = err_data_summary,
                aes(x = time,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), width = 0.1) +
  geom_text(data = pairwise_results,
            aes(x = 12.5,
                y = mean_n,
                label = stars))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  coord_cartesian(ylim = c(96,100))+
  theme_prism(base_size = 16,base_fontface = "plain",base_line_size = 1)+
  scale_y_continuous(
    breaks = seq(96, 100, 1),
    guide = "prism_offset"
  )+ 
  scale_x_continuous(
    breaks = c(0, 6, 12),
    labels = c("0h","6h" ,"12h"),
    guide = "prism_offset_minor",
  )+
  labs(y = "Body weight(% at hour 0)",
       x = "Hour")
ggsave("f6i.pdf",width = 6,height = 4)

#f6j-----------------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("xirao.xlsx", sheet = "图3WT健康评分")

data1 = data %>% pivot_longer(-c(ID,group),values_to = "weight",names_to = "time")

err_data = data1 %>% select(group,time,weight)

library(rstatix)
pairwise_results <- err_data %>%
  filter(time == 12) %>% # 按时间点分组
  pairwise_t_test(weight ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) 


err_data_summary <- err_data %>%
  group_by(group,time) %>%
  summarise(
    mean_n = mean(weight, na.rm = TRUE),
    sem_n = sd(weight, na.rm = TRUE) / sqrt(n())
  )

pairwise_results = pairwise_results %>% left_join(err_data_summary %>% filter(time == 12) %>% select(group,mean_n),by=c("group1" = "group"))

err_data$time = as.numeric(err_data$time)
err_data_summary$time = as.numeric(err_data_summary$time)

ggplot() +
  geom_line(data = err_data_summary,
            aes(x = time, y = mean_n,color = group,group = group),
            linewidth = 1) +
  geom_point(data = err_data_summary,
             aes(x = time,y = mean_n,color = group), size = 3) +
  geom_errorbar(data = err_data_summary,
                aes(x = time,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), width = 0.1) +
  geom_text(data = pairwise_results,
            aes(x = 12.5,
                y = mean_n,
                label = stars))+
  coord_cartesian(ylim = c(0,5))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 16,base_fontface = "plain",base_line_size = 1)+
  scale_y_continuous(
    breaks = seq(0, 5, 1),
    guide = "prism_offset"
  )+ 
  scale_x_continuous(
    breaks = c(0, 6, 12),
    labels = c("0h","6h" ,"12h"),
    guide = "prism_offset_minor",
  )+
  labs(y = "Sepsis health score",
       x = "Hour")
ggsave("f6j.pdf",width = 6,height = 4)

#f6k------------------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("xirao.xlsx", sheet = "图4肺损伤")

data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1

library(rstatix)
pairwise_results <- data1 %>%
  pairwise_t_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  add_y_position()

err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,12))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 12, 2),
    guide = "prism_offset"
  )+ 
  labs(y = "Lung histopathology score",
       x = "",
       title = "") 
ggsave("f6k.pdf",width = 5,height = 4)

#f6l------------------
rm(list = ls())
colpal=readRDS("../colpal.rds")

# colpal = c(colpal,c("Nod2-/--PBS" = "#ff8bbf",
#                     "Nod2-/--L. reuteri" = "#49c8cb"))
# saveRDS(colpal,"../colpal.rds")

data <- read_excel("xirao.xlsx", sheet = "图5ko体重")

data1 = data %>% pivot_longer(-c(ID,group),values_to = "weight",names_to = "time")

err_data = data1 %>% select(group,time,weight)

library(rstatix)
pairwise_results <- err_data %>%
  filter(time == 12) %>% # 按时间点分组
  pairwise_t_test(weight ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) 

err_data_summary <- err_data %>%
  group_by(group,time) %>%
  summarise(
    mean_n = mean(weight, na.rm = TRUE),
    sem_n = sd(weight, na.rm = TRUE) / sqrt(n())
  )

pairwise_results = pairwise_results %>% left_join(err_data_summary %>% filter(time == 12) %>% select(group,mean_n),by=c("group1" = "group"))

err_data$time = as.numeric(err_data$time)
err_data_summary$time = as.numeric(err_data_summary$time)
ggplot() +
  geom_line(data = err_data_summary,
            aes(x = time, y = mean_n,color = group,group = group),
            linewidth = 1) +
  geom_point(data = err_data_summary,
             aes(x = time,y = mean_n,color = group), size = 3) +
  geom_errorbar(data = err_data_summary,
                aes(x = time,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), width = 0.1) +
  geom_text(data = pairwise_results,
            aes(x = 12.5,
                y = mean_n,
                label = stars))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  coord_cartesian(ylim = c(94,100))+
  theme_prism(base_size = 16,base_fontface = "plain",base_line_size = 1)+
  scale_y_continuous(
    breaks = seq(94, 100, 2),
    guide = "prism_offset"
  )+ 
  scale_x_continuous(
    breaks = c(0, 6, 12),
    labels = c("0h","6h" ,"12h"),
    guide = "prism_offset_minor",
  )+
  labs(y = "Body weight(% at hour 0)",
       x = "Hour")
ggsave("f6l.pdf",width = 7,height = 4)

#f6m------------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("xirao.xlsx", sheet = "图6ko健康评分")

data1 = data %>% pivot_longer(-c(ID,group),values_to = "weight",names_to = "time")

err_data = data1 %>% select(group,time,weight)

library(rstatix)
pairwise_results <- err_data %>%
  filter(time == 6) %>% # 按时间点分组
  pairwise_t_test(weight ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) 


err_data_summary <- err_data %>%
  group_by(group,time) %>%
  summarise(
    mean_n = mean(weight, na.rm = TRUE),
    sem_n = sd(weight, na.rm = TRUE) / sqrt(n())
  )

pairwise_results = pairwise_results %>% left_join(err_data_summary %>% filter(time == 12) %>% select(group,mean_n),by=c("group1" = "group"))

err_data$time = as.numeric(err_data$time)
err_data_summary$time = as.numeric(err_data_summary$time)

ggplot() +
  geom_line(data = err_data_summary,
            aes(x = time, y = mean_n,color = group,group = group),
            linewidth = 1) +
  geom_point(data = err_data_summary,
             aes(x = time,y = mean_n,color = group), size = 3) +
  geom_errorbar(data = err_data_summary,
                aes(x = time,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), width = 0.1) +
  geom_text(data = pairwise_results,
            aes(x = 12.5,
                y = mean_n,
                label = stars))+
  coord_cartesian(ylim = c(0,6))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 16,base_fontface = "plain",base_line_size = 1)+
  scale_y_continuous(
    breaks = seq(0, 6, 2),
    guide = "prism_offset"
  )+ 
  scale_x_continuous(
    breaks = c(0, 6, 12),
    labels = c("0h","6h" ,"12h"),
    guide = "prism_offset_minor",
  )+
  labs(y = "Sepsis health score",
       x = "Hour")
ggsave("f6m.pdf",width = 7,height = 4)

#f6n-----------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("xirao.xlsx", sheet = "图7肺损伤")

data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1

library(rstatix)
pairwise_results <- data1 %>%
  pairwise_t_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  add_y_position()

err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,8))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
    scale_y_continuous(
    breaks = seq(0, 2.5, 0.3),
    guide = "prism_offset"
  )+ 
  scale_x_continuous(
    breaks = seq(-10, 32, 2),
    guide = "prism_offset"
  )+ 
  labs(y = "Lung histopathology score",
       x = "",
       title = "") 
ggsave("f6n.pdf",width = 5,height = 4)

# sup1---------------
setwd("../24.wet_valid/")

rm(list =ls())

data <- read_excel("huidi.xlsx", sheet = "nod2_new")
data = data[1,]
data = as.data.frame(t(data)) %>% rownames_to_column("sample") %>%
  filter(sample != "菌株")

colnames(data) = c("group","od")
#data = data[-1,]
data$od = as.numeric(data$od)

data$group = factor(data$group,levels = c("Positive control","Negative control", 
                                              "Lactobacillus reuteri DSM 17938", 
                                              "Enterocuccus faecium 1", 
                                              "Enterocuccus faecium 3",
                                              "Lactobacillus acidophilus 1", 
                                              "Lactobacillus acidophilus 2", 
                                              "Lactobacillus fermentum 1", 
                                              "Lactobacillus fermentum 2", 
                                              "Lactobacillus fermentum 3", 
                                              "Lactobacillus fermentum 4", 
                                              "Lactobacillus gasseri 1", 
                                              "Lactobacillus gasseri 2", 
                                              "Lactobacillus gasseri 3", 
                                              "Lactobacillus gasseri 4", 
                                              "Lactobacillus salivarius", 
                                              "Escherichia coli"))

ggplot(data = data,aes(x = group, y = od))+
  geom_col(aes(fill  = od),width = 0.6)+
  #geom_text(aes(label = od))+
  scale_fill_gradient2(high = "#f222a9",
                       low = "#0d98ba",mid = "grey",midpoint = 0.1)+
  theme_prism(base_size = 16,
              base_fontface = "plain",
              base_line_size = 0.5)+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  labs(x = "",
       y = "NOD2 activation (OD value)",
       title = "Activates NOD2",
       fill = "OD values")+
  scale_y_continuous(
    breaks = c(seq(0,2,0.2)),
    guide = "prism_offset"
  )+ 
  scale_x_discrete(
    guide = "prism_offset"
  )+
  coord_cartesian(ylim = c(0,1.4))
ggsave("s6a.pdf",width = 13,height = 6)


ggplot(data = data,aes(x = group, y = od))+
  geom_linerange(aes(ymin = 0, ymax = od,color = od),linewidth = 1.5)+
  geom_point(aes(color  = od),size = 5)+
  geom_text(aes(label = od),vjust = -1)+
  scale_fill_gradient2(high = "#f222a9",
                       low = "#0d98ba",mid = "grey",midpoint = 0.1)+
  scale_color_gradient2(high = "#f222a9",
                       low = "#0d98ba",mid = "grey",midpoint = 0.1)+
  theme_prism(base_size = 16,
              base_fontface = "plain",
              base_line_size = 0.5)+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  labs(x = "",
       y = "NOD2 activation (OD value)",
       title = "Activates NOD2",
       fill = "OD values")+
  scale_y_continuous(
    breaks = c(seq(0,2,0.2)),
    guide = "prism_offset"
  )+ 
  scale_x_discrete(
    guide = "prism_offset"
  )+
  coord_cartesian(ylim = c(0,1.4))
ggsave("s6a.pdf",width = 15,height = 6)


ggplot(data = data,aes(x = group, y = 1))+
  geom_tile(aes(fill  = od))+
  geom_text(aes(label = od))+
  scale_fill_gradient2(high = "#f222a9",low = "#0d98ba",mid = "white",midpoint = 0.1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  labs(x = "",
       y = "NOD2 activation (Delta OD value)",
       title = "Activation Capacity of NOD2",
       fill = "Delta OD")

ggsave("s6a.pdf",width = 15,height = 3)

#s6b----------
require(tidyverse)
library(survival)
library(ggsurvfit)
library(survminer)
require(survminer)
library(gridExtra)
rm(list = ls())
data <- read_excel("huidi.xlsx", sheet = "new_suv")
data1 = data %>% mutate(ID = paste("mm",row_number()))

data1 %>% 
  pivot_longer(-c(Hours,ID),names_to = "group", values_to = "OS") %>%
  filter(!is.na(OS)) %>%
  rename("OS.time.y" = "Hours") -> sur_data

write.csv(sur_data,"s6b_table.csv")

#sur_data = read.csv("f6c_table.csv",row.names = 1)
#sur_data = bind_rows(sur_data,sur_data)

surv_data <- with(sur_data, Surv(OS.time.y, OS))
fit <- survfit2(surv_data ~ group, data = sur_data)

fit_data = tidy_survfit(fit)

c1 = c(
  "#a0a5e1",
  "#D95050",
  "#E0B794",
  "#49c8cb",
  "#ff8bbf",
  "orange"
)

#F6B0E3

# c(colpal[1:26],c1) -> colpal
# saveRDS(colpal,"../colpal.rds")

ggsurvplot(fit, data = sur_data,
           fun = "event",
           pval = TRUE, conf.int = FALSE,
           conf.int.style="step",
           surv.scale = "percent",
           palette = unname(c1),
           censor = F,
           censor.size = 10,
           ylim = c(0,0.70),
           risk.table = TRUE,
           risk.table.title = "Note the risk set sizes",
           risk.table.subtitle = "",
           risk.table.caption = "",
           risk.table.height = 0.45,
           break.time.by = 2) -> p1

pairwise_results <- pairwise_survdiff(Surv(OS.time.y, OS) ~ group, 
                                      data = sur_data, 
                                      p.adjust.method = "none")
write.csv(data.frame(pairwise_results[["p.value"]]),"s6b_pairwise_p.csv")
cox_model <- coxph(Surv(OS.time.y, OS) ~ group, data = sur_data)
model_summary <- summary(cox_model)
# 提取风险比
hazard_ratios <- exp(coef(cox_model))
# 提取风险比的置信区间
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

write_csv(results_df,"s6b.csv")

g1 <- ggplotGrob(p1$plot)
g1
ggsave("s6b.pdf",g1 ,width = 4, height = 3.5, units = "in")
p1
ggsave("s6b_table.pdf" ,width = 5, height = 2.5, units = "in")


ggforest(cox_model,data = as.data.frame(sur_data))

#s6c------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "liver_bac")

colpal = c(
  "DSM" = "#a0a5e1",
  "E.coli" = "#E0B794",
  "MDP" = "#49c8cb",
  "DSM+GSK" = "#D95050",
  "MDP+GSK" = "#ff8bbf",
  "PBS" = "orange"
)



data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1
data1$group = factor(data1$group,
                     levels = 
                       c("PBS","DSM","MDP",
                         "E.coli","DSM+GSK",
                         "MDP+GSK"))

library(rstatix)
pairwise_results <- data1 %>%
  pairwise_wilcox_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  add_y_position()
write_csv(as.data.frame(pairwise_results),"s6c.csv")


err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),
           color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,1000000))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 14,base_fontface = "plain",base_line_size = 1)+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 1000000, 200000),
    guide = "prism_offset"
  )+ 
  labs(y = "CFU/mL",
       x = "")
ggsave("s6c.pdf",width = 6,height = 6)
#s6d------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "mes_bac")

colpal = c(
  "DSM" = "#a0a5e1",
  "E.coli" = "#E0B794",
  "MDP" = "#49c8cb",
  "DSM+GSK" = "#D95050",
  "MDP+GSK" = "#ff8bbf",
  "PBS" = "orange"
)



data %>% 
  rownames_to_column("rowid") %>% 
  pivot_longer(-rowid,names_to = "group",values_to = "score") %>%
  filter(!is.na(score)) -> data1
data1$group = factor(data1$group,
                     levels = 
                       c("PBS","DSM","MDP",
                         "E.coli","DSM+GSK",
                         "MDP+GSK"))

library(rstatix)
pairwise_results <- data1 %>%
  pairwise_wilcox_test(score ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  add_y_position()
write_csv(as.data.frame(pairwise_results),"s6d.csv")


err_data_summary <- data1 %>%
  group_by(group) %>%
  summarise(
    mean_n = mean(score, na.rm = TRUE),
    sem_n = sd(score, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = group,y = mean_n,fill = group),
           color = "#d9d9d8",
           alpha = 0.2,width = 0.78)+
  geom_errorbar(data = err_data_summary,
                aes(x = group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = group), 
                width = 0.2,size = 1) +
  geom_jitter(data = data1,aes(x = group,y = score,color = group),
              size = 1.5,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,1500000))+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_prism(base_size = 14,base_fontface = "plain",base_line_size = 1)+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 1200000, 200000),
    guide = "prism_offset"
  )+ 
  labs(y = "CFU/mL",
       x = "")
ggsave("s6d.pdf",width = 6,height = 6)

# f6o-------------------

# # 处理原始数据
# data1 = data %>% separate(GA,into = c("w","d"),sep = "\\+") 
# data1$w = as.numeric(data1$w)
# data1$d = as.numeric(data1$d)
# 
# data1 %>% mutate(d = replace_na(d,0)) %>%
#   mutate(GA_daily = 7*w+d) %>%
#   select(-c(w,d)) %>%
#   mutate(adjust_time = Age - Onset_Day) %>%
#   mutate(adjust_time = case_when(
#     adjust_time >= 0 ~ adjust_time+1,
#     adjust_time < 0 ~ adjust_time
#   )) %>% 
#   mutate(Feeding = case_when(
#     Feeding == "A" ~ "A",
#     Feeding != "A" ~ "B"
#   )) -> data1
# 
#  
# data1 %>% group_by(Patient_ID) %>%
#   mutate(ma = max(deltaOD635),
#          mi = min(deltaOD635)) %>% 
#   ungroup() %>%
#   mutate(norm_od = (deltaOD635-mi)/(ma-mi)) %>%
#   select(-c(mi,ma)) -> data_p
# 
# write.csv(data_p,"f6o.csv")
setwd("./24.wet_valid/")
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "Sheet3")

table(data$Patient_ID,data$group)

require(rstatix)

data  %>% 
  filter(adjust_time >= -2) %>%
  mutate(on = case_when(
    adjust_time > 14 ~ "after",
    adjust_time >= 0 ~ "during",
    adjust_time < 0 ~ "before"
  )) %>%
  group_by(on) %>% 
  wilcox_test(deltaOD635 ~ group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) -> test.data

write.csv(test.data,"f6o_test.csv")


err_data_summary <- data %>%
  filter(adjust_time <= 14 & adjust_time >= -2) %>%
  group_by(group,adjust_time) %>%
  summarise(
    mean_n = mean(deltaOD635, na.rm = TRUE),
    sem_n = sd(deltaOD635, na.rm = TRUE) / sqrt(n())
  )

# ggplot() +
#   geom_line(data = err_data_summary,
#             aes(x = adjust_time, y = mean_n,color = group,group = group),
#             linewidth = 1) +
#   geom_point(data = err_data_summary,
#              aes(x = adjust_time,y = mean_n,color = group), size = 3) +
#   geom_errorbar(data = err_data_summary,
#                 aes(x = adjust_time,ymin = mean_n - sem_n, 
#                     ymax = mean_n + sem_n,color = group), width = 0.1)+
#   theme_prism(base_size = 16,base_fontface = "plain",base_line_size = 1) +
#   scale_y_continuous(
#     breaks = seq(0, 2.5, 0.3),
#     guide = "prism_offset"
#   )+ 
#   scale_x_continuous(
#     breaks = seq(-10, 20, 2),
#     guide = "prism_offset"
#   )+ 
#   labs(y = "NOD2 activation (OD value)",
#        x = "Days relative to intervention")


ggplot(data %>% filter(adjust_time <= 14 & adjust_time >= -2)) +
  # geom_rect(aes(xmin = 0, xmax = 14, ymin = 0, ymax = 2.4),
  #           fill = "#A7D1CD", alpha = 0.01) +
  # geom_rect(aes(xmin = -2, xmax = 0, ymin = 0, ymax = 2.4),
  #           fill = "grey", alpha = 0.01) +
  geom_point(aes(color = group,x = adjust_time,y = deltaOD635),size = 3)+
  # geom_line(data = err_data_summary,
  #           aes(x = adjust_time, y = mean_n,color = group,group = group),
  #           linewidth = 1) +
  # geom_point(data = err_data_summary,
  #            aes(x = adjust_time,y = mean_n,color = group), size = 3) +
  # geom_errorbar(data = err_data_summary,
  #               aes(x = adjust_time,ymin = mean_n - sem_n, 
  #                   ymax = mean_n + sem_n,color = group), width = 0.1,alpha = 0.2)+
  geom_smooth(aes(x = adjust_time,y = deltaOD635,
    fill = group,color = group),span = 0.8)+
  scale_fill_manual(values = c("con" = "#F4C573","interv" = "#49c8cb"))+
  scale_color_manual(values = c("con" = "#F4C573","interv" = "#49c8cb"))+
  theme_prism(base_size = 20,base_fontface = "plain",
              base_line_size = 0.5)+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 2.5, 0.3),
    guide = "prism_offset"
  )+ 
  scale_x_continuous(
    breaks = seq(-10, 22, 2),
    guide = "prism_offset"
  )+ 
  labs(y = "NOD2 activation (OD value)",
       x = "Days relative to intervention")+
  coord_cartesian(xlim = c(-2,14))
ggsave("f6o.pdf",width = 8.5,height = 6)

ggplot(data %>% filter(group == "interv") %>% filter(adjust_time >= -2),
       aes(x = adjust_time,y = deltaOD635)) +
  geom_rect(aes(xmin = 0, xmax = 22, ymin = 0, ymax = 2.4),
            fill = "#A7D1CD", alpha = 0.01) +
  geom_rect(aes(xmin = -2, xmax = 0, ymin = 0, ymax = 2.4),
            fill = "grey", alpha = 0.01) +
  geom_point(color = "#49c8cb")+
  geom_smooth(span = 0.8,color = "#49c8cb",fill = "#49c8cb") +
  theme_prism(base_size = 14,base_fontface = "plain",
              base_line_size = 0.5)+
  scale_y_continuous(
    breaks = seq(0, 2.5, 0.3),
    guide = "prism_offset"
  )+ 
  scale_x_continuous(
    breaks = seq(-10, 32, 2),
    guide = "prism_offset"
  )+ 
  labs(y = "NOD2 activation (OD value)",
       x = "Days relative to intervention") 
ggsave("f6o1.pdf",width = 8.5,height = 6)
  
ggplot(data %>% filter(group == "con") %>% filter(adjust_time >= -2),
       aes(x = adjust_time,y = deltaOD635)) +
  geom_rect(aes(xmin = 0, xmax = 22, ymin = 0, ymax = 2.4),
            fill = "#A7D1CD", alpha = 0.01) +
  geom_rect(aes(xmin = -2, xmax = 0, ymin = 0, ymax = 2.4),
            fill = "grey", alpha = 0.01) +
  geom_point(color = "#F4C573")+
  geom_smooth(span = 0.8,color = "#F4C573",fill = "#F4C573") +
  theme_prism(base_size = 14,base_fontface = "plain",base_line_size = 0.5)+
  scale_y_continuous(
    breaks = seq(0, 2.5, 0.3),
    guide = "prism_offset"
  )+ 
  scale_x_continuous(
    breaks = seq(-10, 32, 2),
    guide = "prism_offset"
  )+ 
  labs(y = "NOD2 activation (OD value)",
       x = "Days relative to intervention")

require(tidyverse)
require(gtsummary)

data %>% group_by(Patient_ID) %>%
  slice_head(n = 1) %>%
  ungroup() -> data_pat

data_pat$Duration = as.double(data_pat$Duration)
data_pat %>% mutate(Antibiotics_Duration = 
                      case_when(
                        Duration >= 7 ~ "long",
                        Duration < 7 ~ "short",
                      )) -> d

data_pat %>% mutate(big = 
                      case_when(
                        GA_daily >= 38*7 ~ Feeding
                      )) %>% 
  mutate(small = 
                case_when(
                  GA_daily < 38*7 ~ Feeding
                )) -> dd

theme_gtsummary_journal(journal = "jama")
theme_gtsummary_compact()
d %>%
  tbl_summary(
    include = c("Age","Nature_Delivery",
                "Gender","Feeding","Antibiotics",
                "Antibiotics_Duration","GA_daily"),
    by = group,
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c(
      "{mean} ({sd})",
      "{min}, {max}"
    ),
    missing = "no"
  ) |> 
  add_p() |>
  modify_header(label = "**Variable**") |> # update the column header
  bold_labels() ->s1

as_gt(s1) |> 
  gt::gtsave(filename = "./f7b.pdf")

ggsave("f7b.pdf")
as_hux_table(s1) %>% huxtable::quick_xlsx(file = "patient_meta_diff_raw.xlsx")

#f6k-------------
rm(list = ls())
colpal=readRDS("../colpal.rds")
data <- read_excel("huidi.xlsx", sheet = "LOS-NOD2")

data  -> data1

library(rstatix)
pairwise_results <- data1 %>%
  pairwise_t_test(OD ~ Group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", # p值 ≤ 0.001 添加 ***
    p <= 0.01  ~ "**",  # p值 ≤ 0.01 添加 **
    p <= 0.05  ~ "*",   # p值 ≤ 0.05 添加 *
    TRUE       ~ "ns"     # 其他情况不添加星号
  )) %>% 
  add_y_position()

err_data_summary <- data1 %>%
  group_by(Group) %>%
  summarise(
    mean_n = mean(OD, na.rm = TRUE),
    sem_n = sd(OD, na.rm = TRUE) / sqrt(n())
  )


ggplot()+
  geom_col(data = err_data_summary,
           aes(x = Group,y = mean_n,fill = Group),color = "#d9d9d8",
           alpha = 0.2,width = 0.65)+
  geom_errorbar(data = err_data_summary,
                aes(x = Group,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = Group), 
                width = 0.2,size = 0.75) +
  geom_jitter(data = data1,aes(x = Group,y = OD,color = Group),
              size = 3,width = 0.1)+
  geom_text(data = pairwise_results,
            aes(x = group1,
                y = y.position/3,
                label = stars))+
  coord_cartesian(ylim = c(0,1))+
  scale_fill_manual(values = c("#49c8cb","#ff8bbf"))+
  scale_color_manual(values = c("#49c8cb","#ff8bbf"))+
  theme_prism(base_size = 16,base_fontface = "plain",base_line_size = 0.5) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    guide = "prism_offset"
  )+ 
  labs(y = "NOD2 activation (Relative OD value)",
       x = "",
       title = "") 
ggsave("f6kk.pdf",width = 5,height = 6)
