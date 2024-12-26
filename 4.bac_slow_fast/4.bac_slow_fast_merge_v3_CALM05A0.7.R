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

all_data =  read.csv("../1.data_tidy/all_cohort_sample_num_8cohort_0905.csv",row.names = 1)
all_data = all_data %>% filter(Age_DOL <=  42 )

all_data$cohort = factor(all_data$cohort,levels = unique(all_data$cohort))


model_data = all_data %>%
  filter(cohort == "CALM2005_2_A") %>%
  sample_frac(7/10, replace = FALSE)


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

require(ggpmisc)
ggplot(best_tuned_pred,aes(x = truth, y = response))+
  geom_jitter(aes(fill = model),alpha = 0.5,size = 3,color = "white",shape = 21)+
  geom_smooth(method = 'lm',formula = y~x,color = "#ff8bbf",se = T)+
  scale_fill_manual(values = c("#7DCCC5","#C1A6EE"))+
  facet_wrap(".~model")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Actual day of life", y = "Predicted day of life")

ggsave("MIX_Age_DOL_predict_corplot.0d-42d.pdf", width = 8,height = 4)

write.csv(best_tuned_pred,"sum_model.data.csv")

# why robert have line in predict 8 days

best_tuned_pred %>%
  filter(cohort == "Robert_2024_M") %>%
  filter(truth > 10) %>%
  filter(response < 10) %>% pull(Sample_ID) %>%
  unique() -> nnid_s

best_tuned_pred %>%
  filter(cohort == "Robert_2024_M") %>%
  filter(truth > 10) %>%
  filter(response < 10) %>% pull(PatientID) %>%
  unique() -> nnid_p

nndata = all_data %>% filter(Sample_ID %in% nnid)
nndata_p = all_data %>% filter(PatientID %in% nnid_p)

write.csv(nndata,"robert_line_patient.csv")

ggplot(best_tuned_pred %>%
         filter(PatientID %in% nnid_p),aes(x = truth, y = response))+
  geom_jitter(aes(fill = cohort),alpha = 1,size = 3,color = "white",shape = 21,show.legend = F)+
  geom_smooth(method = 'lm',formula = y~x,color = "#5d4990",se = T)+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1)+
  scale_fill_manual(values = colpal)+
  facet_wrap(".~cohort",nrow = 1,scales = "free_y")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Age_DOL", y = "MA")
ggsave("robert_line_patient.pdf",width = 4,height = 3)

colpal = readRDS("../colpal.rds")
colpal = c(c("Robert_2024_M" = "#E0B794"),colpal)

ggplot(best_tuned_pred %>%
         filter(!cohort  %in% c("CALM2005_2_A",
                                "Brooks_2017_M",
                                "Olm_2019_M",
                                "Raveh_2015_M")),aes(x = truth, y = response))+
  geom_jitter(aes(fill = cohort),alpha = 1,size = 3,color = "white",shape = 21,show.legend = F)+
  geom_smooth(method = 'lm',formula = y~x,color = "#5d4990",se = T)+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1)+
  scale_fill_manual(values = colpal)+
  facet_wrap(".~cohort",nrow = 1,scales = "free_y")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Age_DOL", y = "MA")->p1

ggplot(best_tuned_pred %>% 
         filter(cohort == "CALM2005_1_A"),
       aes(x = truth, y = response))+
  geom_jitter(aes(fill = cohort),alpha = 1,size = 3,color = "white",shape = 21,show.legend = F)+
  geom_smooth(method = 'lm',formula = y~x,color = "#5d4990",se = T)+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1)+
  scale_fill_manual(values = colpal)+
  facet_wrap(".~cohort",nrow = 1,scales = "free_y")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Age_DOL", y = "MA") -> p2

ggplot(best_tuned_pred %>% 
         filter(cohort == "CALM2005_2_A"),
       aes(x = truth, y = response,color = model))+
  geom_jitter(aes(fill = model),alpha = 1,size = 3,color = "white",shape = 21,show.legend = F)+
  geom_smooth(aes(color = model),method = 'lm',formula = y~x,se = T,show.legend = F)+
  stat_poly_eq(
    aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
    formula = y~x, parse = TRUE,position = "identity",
    label.x.npc = 0.1) +
  scale_fill_manual(values = c("#BAD5F7","#85dab1"))+
  scale_color_manual(values = c("#BAD5F7","#85dab1"))+
  facet_wrap(".~model",nrow = 1,scales = "free_y")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Age_DOL", y = "MA") -> p3
require(patchwork)
(p3)/p1

ggsave("split_cohort_predict_model.pdf",width = 9,height = 6)
p1
ggsave("split_cohort_predict_model_displayvalue.pdf",width = 20,height = 3)


best_tuned_pred %>% 
  mutate(cohort = case_when(
    model == "Predict" & cohort == "CALM2005_2_A" ~ "CALM2005_Train",
    model == "Validation" & cohort == "CALM2005_2_A" ~ "CALM2005_Test",
    TRUE ~ cohort
  )) -> best_tuned_pred

colpal = c(c("Robert_2024_M" = "#E0B794"),colpal)


colpal = c(c("Robert_2024_M" = "#E0B794","CALM2005_Test" = "#BAD5F7","CALM2005_Train" = "#85dab1"),colpal)

best_tuned_pred$cohort = factor(best_tuned_pred$cohort,
                                levels = c("CALM2005_Test", "CALM2005_Train", "CALM2005_1_A", "Rao_2021_A", "Lauren_2022_M", "Robert_2024_M", "Brooks_2017_M", "Olm_2019_M", "Raveh_2015_M")
                                
)

ggplot(best_tuned_pred %>%
         filter(!cohort  %in% c(
           "CALM2005_2_A",
           "Brooks_2017_M",
           "Olm_2019_M",
           "Raveh_2015_M")),aes(x = truth, y = response))+
  geom_point(aes(fill = cohort),alpha = 1,size = 3,color = "white",shape = 21,show.legend = F)+
  geom_smooth(method = 'lm',formula = y~x,color = "#5d4990",se = T)+
  # stat_poly_eq(
  #   aes(label =  paste(..eq.label..,..rr.label..,..p.value.label.., sep = "~~~~")),
  #   formula = y~x, parse = TRUE,position = "identity",
  #   label.x.npc = 0.1)+
  coord_cartesian(ylim = c(8,30))+
  scale_fill_manual(values = colpal)+
  facet_wrap(".~cohort",nrow = 2,scales = "free_y")+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Age DOL", y = "Predict Age(Days)")
ggsave("split_cohort_predict_model_3x2.pdf",width = 9,height = 5.6)


results <- best_tuned_pred %>%
  group_by(model, cohort) %>%
  summarise(metrics = list(calculate_metrics(truth, response)), .groups = 'drop')

results <- results %>%
  unnest(cols = c(metrics))

c("model","cohort","MAE","MSE","RMSE","Rho","R2")

write.csv(results,"split_cohort_predict_model.csv")