xbox::chdir("../4.bac_slow_fast")

library(tidyverse)
library(mlr3verse)
library(mlr3extralearners)
library(mlr3tuning)
library(data.table)
library(pROC)
library(mlr3viz)
library(future)

rm(list = ls())
set.seed(2022)
colpal = readRDS("../colpal.rds")

# ------------
all_data =  read.csv("../1.data_tidy/all_cohort.sum_7cohort_final.csv",row.names = 1)
all_data = all_data %>% filter(Age_DOL <=  42 )

all_data$cohort = factor(all_data$cohort,levels = unique(all_data$cohort))

all_data %>% select(Sepsis_style,cohort,PatientID) -> test
test = unique(test)
write.csv(as.data.frame(table(test$Sepsis_style,test$cohort)),"Sepsis_style_cohort.patient.table.csv")

model_data = all_data %>%
  group_by(cohort) %>%
  sample_frac(1/3, replace = FALSE)

feature = model_data[,-c(49,51,52,53)] %>% filter(Age_DOL <=  42 )


task <- as_task_regr(feature, target = "Age_DOL", id = "regression")
task$set_col_roles("Sample_ID", roles = "name")
saveRDS(task,"model_pridict_task.rds")

learner <- lrn("regr.ranger", num.threads = 5, importance = "impurity")
search_space <- ps(
  mtry = p_int(lower = 10, upper = length(task$feature_names))
)
resampling <- rsmp("cv", folds = 10)
measure = msr("regr.mse")
none <- trm("none")
tuner <- tnr("grid_search", resolution=10)
at <- AutoTuner$new(
  learner = learner,
  resampling = resampling,
  measure = measure,
  search_space = search_space,
  terminator = none,
  tuner = tuner
)
set.seed(2022)

at$train(task)

at$tuning_result

learner$param_set$values = at$tuning_result$learner_param_vals[[1]]

learner$train(task)

best_nr <- at$archive$best()$batch_nr
best_pred <- at$archive$predictions(best_nr)
best_tuned_pred <- map_dfr(best_pred, ~ bind_rows(as.data.table(.x)))

model_predict_result  = best_tuned_pred
aid = model_predict_result$row_ids
model_predict_result$PatientID = model_data$PatientID[aid]
model_predict_result$Sample_ID = model_data$Sample_ID[aid]
model_predict_result$Age_DOL = model_data$Age_DOL[aid]
model_predict_result$Sepsis_style = model_data$Sepsis_style[aid]
model_predict_result$NEC = model_data$NEC[aid]
model_predict_result$cohort = model_data$cohort[aid]

source("../addin.r")
regr.result = calculate_metrics(model_predict_result$truth,model_predict_result$response)
write.csv(regr.result,"regr.result_train.csv")
saveRDS(learner,"model_pridict_learner.rds")

library(ggpmisc)

age_prd = best_tuned_pred

rid = model_data %>% select(Sample_ID) %>% pull

test_pre = all_data %>% filter(!Sample_ID %in% rid)

test_c = test_pre[,colnames(feature)]

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

best_tuned_pred$PatientID = test_pre$PatientID
best_tuned_pred$Sample_ID = test_pre$Sample_ID
best_tuned_pred$Age_DOL = test_pre$Age_DOL
best_tuned_pred$Sepsis_style = test_pre$Sepsis_style
best_tuned_pred$NEC = test_pre$NEC
best_tuned_pred$cohort = test_pre$cohort

model_predict_result$model = "Predict"
best_tuned_pred$model = "Validation"

best_tuned_pred = bind_rows(model_predict_result,best_tuned_pred)

pre_age_pdata = best_tuned_pred %>% filter(Sepsis_style != "EOS")

annot_dat = data.frame(table(best_tuned_pred$model))
colnames(annot_dat) = c("model","n")
test_d = best_tuned_pred %>% filter(model == "Validation")
pred_d = best_tuned_pred %>% filter(model == "Predict")
annot_dat$mae = c(mean(abs(pred_d$truth - pred_d$response)),mean(abs(test_d$truth - test_d$response)))
write.csv(annot_dat,"MIX_Age_DOL_predict_corplot.0d-42d.csv")


ggplot(best_tuned_pred,aes(x = truth, y = response))+
  geom_jitter(aes(fill = model),alpha = 0.5,size = 3,color = "white",shape = 21)+
  geom_smooth(method = 'lm',formula = y~x,color = "#ff8bbf",se = T)+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1)+
  scale_fill_manual(values = c("#7DCCC5","#C1A6EE"))+
  facet_wrap(".~model")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Actual day of life", y = "Predicted day of life")

ggsave("MIX_Age_DOL_predict_corplot.0d-42d.pdf", width = 8,height = 4)

write.csv(best_tuned_pred,"sum_model.data.csv")

ggplot(best_tuned_pred %>% filter(model == "Validation"),aes(x = truth, y = response))+
  geom_jitter(aes(fill = cohort),alpha = 1,size = 3,color = "white",shape = 21)+
  geom_smooth(method = 'lm',formula = y~x,color = "#5d4990",se = T)+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1)+
  scale_fill_manual(values = colpal)+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Age_DOL", y = "MA") -> p1

ggplot(best_tuned_pred %>% filter(model == "Predict"),aes(x = truth, y = response))+
   geom_jitter(aes(fill = cohort),alpha = 1,size = 3,color = "white",shape = 21)+
  geom_smooth(method = 'lm',formula = y~x,color = "#5d4990",se = T)+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1)+
  scale_fill_manual(values = colpal)+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Age_DOL", y = "MA") -> p2
require(patchwork)
require(tidyverse)
require(ggpubr)
require(ggpmisc)
p2+p1
ggsave("MIX_Age_DOL_predict_corplot.0d-42d_2.pdf", width = 8,height = 4)

pre_age_pdata = best_tuned_pred %>% filter(Sepsis_style != "EOS")

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
  coord_cartesian(ylim= c(7,28))+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("Sepsis_fastslow_by_cohort_nopoint.png",width = 4,height = 2.8)
ggsave("Sepsis_fastslow_by_cohort_nopoint.pdf",width = 4,height = 2.8)

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
  coord_cartesian(ylim= c(7,28))+
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
  coord_cartesian(ylim= c(7,28))+
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
  #slope difference test
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

density_curve <- density(sv)
peak_index <- which.max(density_curve$y)
middle_value <- density_curve$x[peak_index]

slope_th  = median(round(sv,2))

slop_data %>% mutate(type = case_when(
  slope >= slope_th  ~ "fast",
  slope < slope_th ~ "slow",
)) -> slop_data
slop_data = slop_data %>% rename("group" = "Sepsis_style")
slop_data_unique = slop_data %>% filter(group != "EOS")%>% select(group,slope,intercept,type,PatientID) %>% distinct()

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

write.csv(slop_data,"Sepsis_add_trend_cluster.meta.csv")

slop_data = read.csv("./Sepsis_add_trend_cluster.meta.csv")
slop_data_unique = slop_data %>% select(slope,intercept,type,PatientID) %>% distinct()

slop_data_unique = slop_data_unique %>% arrange(slope)
slop_data_unique$PatientID = factor(slop_data_unique$PatientID,levels = slop_data_unique$PatientID)
cplot = slop_data_unique %>% filter(!is.na(slope))
cplot$slope[cplot$slope >= 3] = 3
cplot$slope[cplot$slope <= -3] = -3
cplot$slope = (cplot$slope-min(cplot$slope))/(max(cplot$slope) - min(cplot$slope))
ggplot(data = cplot) +
  geom_linerange(aes(x = PatientID,
                     ymin = 0,
                     ymax = slope,
                     color = type),
                 linewidth = 0.35)+
  #geom_point(aes(x = PatientID,y = slope,color = type),size = 0.5)+
  scale_color_manual(values =  c("#ff8bbf","#a0a5e1"))+
  ggthemes::theme_clean()+
  theme(axis.text.x = element_blank())
ggsave("ident_fast_split_group.pdf",width = 6,height = 4)


ggplot(data = slop_data %>% 
         filter(slope > -0.2 & slope < 1), 
       aes(x = truth, y = response)) +
  geom_point(aes(color = slope),size = 1, 
             alpha = 0.2,show.legend = T) +
  stat_smooth(alpha = 0.1,method = "lm", formula = y ~ x, 
              linewidth = 1.5, aes(group = PatientID,
                                   color = slope),
              se = F) +
  stat_smooth(method = "lm", formula = y ~ x, 
              linewidth = 2,se = T,
              alpha = 0.6,color = "#D95050",fill = "#fbc7cc") +
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1) + 
  scale_fill_gradient2(low = "#a0a5e1",mid = "#e2e4f6",high = "#ff8bbf",midpoint = 0.36) +
  scale_color_gradient2(low = "#a0a5e1",mid = "#e2e4f6",high = "#ff8bbf",midpoint = 0.36) +
  labs(x = "Age DOL", y = "Predicted Age")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank()) 

ggsave("ident_fast_split_group_line.pdf",width = 5,height = 4)


ggplot(data = slop_data_unique %>% filter(slope > -3)) +
  geom_boxplot(aes(y = slope,x = type,fill = type))+
  #geom_point(aes(x = PatientID,y = slope,color = type),size = 0.5)+
  scale_fill_manual(values =  c("#ff8bbf","#a0a5e1"))+
  ggthemes::theme_clean()+
  theme(axis.text.x = element_blank(),
        panel.grid.major.y = element_blank())

ggsave("ident_fast_split_group_boxplot.pdf",width = 6,height = 4)