
xbox::chdir("../14.cox_hr_sur_curve")
require(tidyverse)
library(survival)
library(ggsurvfit)
library(survminer)
require(survminer)
library(gridExtra)
rm(list = ls())

rm(list = ls())

c1_fs = read.csv("../4.bac_slow_fast/矫正前每个队列随机13_allsample_peakvalue///Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv",row.names = 1)
all_meta = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "ID"))

c1_fs = read.csv("../4.bac_slow_fast/矫正前每个队列随机13_allsample_peakvalue//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
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

fit_data = tidy_survfit(fit)
ggplot(fit_data,aes(x=time,y=estimate))+
  geom_point(aes(color=strata))+
  geom_line(aes(group=strata,color=strata),
            linewidth=2)+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = c(0.1,0.3))+
  scale_color_manual(values = c("#56b4e8","#d55e00"),
                     name=NULL)+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1))

#
ggsurvplot(fit, data = sur_data,
           fun = "event",
           pval = TRUE, conf.int = TRUE,
           conf.int.style="step",
           surv.scale = "percent",
           palette = c("#ff8bbf","#a0a5e1"),
           censor = F,
           censor.size = 10,
           ylim = c(0,0.3),
           risk.table = TRUE,
           risk.table.title = "Note the risk set sizes",
           risk.table.subtitle = "and remember about censoring.",
           risk.table.caption = "source code: website.com",
           risk.table.height = 0.45) -> p1

cox_model <- coxph(Surv(OS.time.y, OS) ~ group, data = sur_data)
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

sur_data = sur_data %>% mutate(Weight = Weight/1000)

# data = all_meta2 %>% filter(Age_DOL <=42)
sur_data %>% mutate(anti_use = case_when(
  Antibiotics_Duration_2w > 3 ~ "Y",
  Antibiotics_Duration_2w <= 3 ~ "N",
)) %>% mutate(anti_time = case_when(
  Antibiotics_Duration_2w > 9 ~ "blong",
  Antibiotics_Duration_2w <= 9 ~ "ashort",
)) %>% mutate(anti_sepsis = Sepsis_style) -> sur_data

cox_model = coxph(Surv(OS.time.y, OS) ~ group+
                    Gender+
                    GA_daily+
                    Natural_delivery+
                    Feeding+
                    Antibiotics_Duration_2w+
                    Intrauterine_infection+
                    PE+
                    GDM, data = as.data.frame(sur_data))

ggforest(cox_model,data = as.data.frame(sur_data))

ggsave("cox_test_muti_forest_addGA.pdf",width = 8,height = 4)

cox_model = coxph(Surv(OS.time.y, OS) ~ group+
                    Gender+
                    Weight + GA_daily+
                    Natural_delivery+
                    Feeding+
                    Intrauterine_infection+
                    Antibiotics_Duration_2w+
                    PE+
                    GDM, data = as.data.frame(sur_data))
model_summary <- summary(cox_model)
hazard_ratios <- exp(coef(cox_model))
conf_intervals <- confint(cox_model,level = 0.95)
p_values <- model_summary$coefficients[, "Pr(>|z|)"]

results_df <- data.frame(
  Variable = rownames(model_summary$coefficients),
  Hazard_Ratio = round(hazard_ratios, 2),
  Lower_CI = exp(round(conf_intervals[,1], 2)),
  Upper_CI = exp(round(conf_intervals[,2], 2)),
  P_Value = round(p_values, 4),
  Class = "OS"
)

require(tableone)
ShowRegTable(
  cox_model,                  
  exp = TRUE,               
  digits = 2,               
  pDigits = 3,             
  printToggle = TRUE,      
  quote = FALSE,           
  ciFun = confint           
) -> results_df1

sur_data

stat_data = print(CreateTableOne(vars = c("group","Gender","Weight", 
                                          "GA_daily","Natural_delivery","Feeding","Intrauterine_infection",
                                          "Antibiotics_Duration_2w","PE","GDM"),data = sur_data))



f_data = bind_cols(results_df,results_df1,stat_data[-1,])
colnames(f_data)[c(7,9)] = c("HR(95% CI)","n")

require(forestplot)
forestplot(
  f_data[, c(1, 5,9,7,8)],        
  mean = f_data[, 2],             
  lower = f_data[, 3],            
  upper = f_data[, 4],            
  zero = 1,                        
  boxsize = 0.25,                   
  graphwidth = unit(.25, "npc"),   
  xlab = "",
  xticks = c(-1, 1, 3, 5, 7,9), 
  txt_gp = fpTxtGp(                
    label = gpar(cex = 0.8),       
    ticks = gpar(cex = 1),         
    xlab = gpar(cex = 0.9),        
    title = gpar(cex = 1.2)       
  ),
  lwd.zero = 1.5,                    
  lwd.ci = 1.5,                    
  lwd.xaxis = 1.5,                  
  lty.ci = 1.5,                    
  ci.vertices = T,                 
  ci.vertices.height = 0.05,       
  clip = c(0.1, 8),                
  ineheight = unit(8, 'mm'),       
  line.margin = unit(8, 'mm'),    
  colgap = unit(6, 'mm'),          
  fn.ci_norm = "fpDrawCircleCI",  
  title = "Multivariate Cox Regression",   
  col = fpColors(                
    box = "#ff8bbf",               
    lines = "#a0a5e1",               
    zero = "black"               
  ),
  hrzl_lines = list(
    "1" = gpar(lwd = 0.5, col = "black"),
    "2" = gpar(lwd = 0.5, col = "black"),
    "3" = gpar(lwd = 0.5, col = "black"),
    "4" = gpar(lwd = 0.5, col = "black"),
    "5" = gpar(lwd = 0.5, col = "black"),
    "6" = gpar(lwd = 0.5, col = "black"),
    "7" = gpar(lwd = 0.5, col = "black"),
    "8" = gpar(lwd = 0.5, col = "black"),
    "9" = gpar(lwd = 0.5, col = "black"),
    "10" = gpar(lwd = 0.5, col = "black")
  )
) -> fig

fig
pdf("cox_test_muti_curve_diy3.pdf",width = 10,height = 4)
fig
dev.off()

write.csv(f_data,"muti_hr.csv")

cox_model = coxph(Surv(OS.time.y, OS) ~ group+
                    Gender+
                    GA_daily+
                    Natural_delivery+
                    Feeding+
                    Intrauterine_infection+
                    Antibiotics_Duration_2w+
                    PE+GDM, data = as.data.frame(sur_data))
ggforest(cox_model,data = as.data.frame(sur_data))

ggsave("cox_test_muti_forest_delweight.pdf",width = 6,height = 5.5)


# AUC --------------
require()
setwd("./14.cox_hr_sur_curve/")
rm(list = ls())
bac_i = read.csv("../11.bac_index/bac_index_risk.csv",row.names = 1)
sur_data = read.csv("sur_data.csv",row.names = 1)

sur_data = sur_data %>% filter(is.na(Onset_Day) | Age_DOL < Onset_Day) 


sur_data = sur_data %>% 
  left_join(bac_i %>% select(SampleID,slope,score,risk),by= c("Sample_ID" = "SampleID"))
require(pROC)
bid = colnames(sur_data)[3:49]

sur_data %>% mutate(
  Age_Group = case_when(
    Age_DOL <= 7 ~ "0d-7d",
    Age_DOL <= 14 ~ "8d-14d",
    Age_DOL <= 21 ~ "15d-21d",
    Age_DOL <= 28 ~ "22d-28d",
    Age_DOL >= 29 ~ "29d-42d",
  )
) -> sur_data

sur_data_unique = sur_data %>% group_by(PatientID,Age_Group) %>% slice_head(n = 1)

write.csv(as.data.frame(table(sur_data_unique$Age_Group,sur_data_unique$Sepsis_style)),
          "age_group_lOS_patient_number.csv")

sur_data[, c(3:50,53, 66,74,77)]  -> sur_datap1
aid = unique(sur_datap1$Age_Group)

for(age in aid){
  sur_datap = sur_datap1 %>% 
    filter(cohort == "CALM2005_2_A") %>%
    filter(Age_Group==age)
  
  results_list <- list()
  for(grp in c("slope",bid)){
    roc_obj <- roc(sur_datap$Sepsis_style, sur_datap[[grp]])
    coords_obj <- coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity"))
    result_df <- data.frame(
      Optimal_threshold = round(coords_obj$threshold, 1),
      Accuracy = round(coords_obj$accuracy, 2),
      Sensitivity = round(coords_obj$sensitivity, 2),
      Specificity = round(coords_obj$specificity, 2),
      AUC = round(auc(roc_obj), 2),
      Class = "ROC",
      Group = grp  
    )
    results_list[[grp]] <- result_df
  }
  
  auc_df <- do.call(rbind, results_list)
  auc_df1 = auc_df %>% group_by(Group) %>% slice_head(n = 1) %>% arrange(AUC)
  auc_df1$Age_DOL = age
  
  auc_df1$Group = factor(auc_df1$Group,levels = rev(auc_df1$Group))
  write.csv(auc_df1,paste0(age,"_roc_bac_slope_sample_c1c2.csv"))
  ggplot(auc_df1,aes(x = Group,y = AUC,fill = AUC))+
    geom_col(color = "black",width = 0.8)+
    geom_hline(yintercept = 0.5)+
    scale_fill_gradient2(low = "grey",mid = "#d6c6ed",high = "#ff5a8c",midpoint = 0.5) +
    coord_cartesian(ylim = c(0.4,0.9))+
    ggthemes::theme_clean()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1.1))+
    ggtitle(age)
  
  ggsave(paste0(age,"_roc_bac_slope_sample_c1c2.pdf"),width = 9,height = 4)
}
