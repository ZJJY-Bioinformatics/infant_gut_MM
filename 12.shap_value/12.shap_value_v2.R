xbox::chdir("../12.shap_value")

library(tidyverse)
library(mlr3verse)
library(mlr3extralearners)
library(mlr3tuning)
library(data.table)
library(pROC)
#-------
rm(list = ls())

# #记载模型数据---------------
# load("../4.bac_slow_fast/矫正前每个队列随机13_allsample_peakvalue/Sepsis_env_var_240410.rdata")
# 
# task_input = task
# learner = learner
# #
# library(iml)
# #
# credit_x = task_input$data(rows = NULL,
#                             cols = task_input$feature_names)
# credit_y = task_input$data(rows = NULL,
#                             cols = task_input$target_names)
# row_names = task_input$row_names$row_name
# predictor = Predictor$new(learner, data = credit_x, y = credit_y)
# #
# shapley_results <- data.frame()
# #
# for (i in 1:nrow(credit_x)) {
#   print(i)
#   sample_row = credit_x[i,]
#   shapley = Shapley$new(predictor, x.interest = sample_row,
#                         sample.size = nrow(credit_x))
# 
#   shapley_row <- data.frame(shap_values = shapley$results$phi,
#                             row.names = shapley$results$feature)
#   gc()
#   colnames(shapley_row) = row_names[i]
#   if (i == 1) {
#     shapley_results <- shapley_row
#   } else {
#     shapley_results <- cbind(shapley_results, shapley_row)
#   }
# }
# 
# write.csv(shapley_results,"shapley_result.csv")

# shap_input_data = data.frame(credit_x)
# rownames(shap_input_data) = task_input$row_names$row_name
# write.csv(shap_input_data,"shapley_value_input_data.csv")

#输出shape值得python输入-----------------------
#xbox::chdir("12.shap_value/")
colpal = readRDS("../colpal.rds")
colpal = c(colpal,c("Stenotrophomonas" = "#FBF6B5",
                    "Acinetobacter" = "#F4C573"))


shapley_results = read.csv("shapley_result.csv",row.names = 1)
