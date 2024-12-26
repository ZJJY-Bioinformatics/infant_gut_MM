xbox::chdir("23.s1_table")
rm(list = ls())
# 加载R包----------
require(tidyverse)

c2_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv")
c2_meta = c2_meta[,c(2,499:552)]

c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv")
c1_meta = c1_meta[,c(2,452:542)]

c1_meta %>% rename(
  "Weight" = "Birth_weight",
  "Gender" = "sex",
  "RDS" = "NRDS",
  "Antibiotics_Duration_all" = "Antibiotic_duration",
  "Intrauterine_infection" = "Infection") %>% 
  mutate(Feeding = "A") %>%
  mutate(Weight = Weight*1000)-> c1_meta

write.csv(c1_meta,"c1_meta.csv")
write.csv(c2_meta,"c2_meta.csv")

c1_meta = read.csv("c1_meta.csv",row.names = 1)
c2_meta = read.csv("c2_meta.csv",row.names = 1)
require(gtsummary)

#
c1_meta$Antibiotics_Duration_1w <- as.numeric(c1_meta$Antibiotics_Duration_1w )

c1_meta %>%
tbl_summary(
  include = colnames(c1_meta)[-c(1,2,3)],
  type = all_continuous() ~ "continuous2",
  statistic = all_continuous() ~ c(
    "{N_nonmiss}",
    "{median} ({p25}, {p75})",
    "{min}, {max}"
  ),
  missing = "no"
) |> 
  add_n() |> # add column with total number of non-missing observations
  modify_header(label = "**Variable**") |> # update the column header
  bold_labels() -> s1
as_hux_table(s1) %>% huxtable::quick_xlsx(file = "CALM05_B_sample.xlsx")

c1_meta %>%
  group_by(PatientID) %>%
  slice_head(n = 1) %>%
  tbl_summary(
    include = colnames(c1_meta)[-c(1,2,3)],
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c(
      "{N_nonmiss}",
      "{median} ({p25}, {p75})",
      "{min}, {max}"
    ),
    missing = "no"
  ) |> 
  add_n() |> # add column with total number of non-missing observations
  modify_header(label = "**Variable**") |> # update the column header
  bold_labels() -> s1p
as_hux_table(s1p) %>% huxtable::quick_xlsx(file = "CALM05_B_patient.xlsx")

# 队列2 ----
c2_meta %>%
  tbl_summary(
    include = colnames(c2_meta)[-c(1,2,3)],
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c(
      "{N_nonmiss}",
      "{median} ({p25}, {p75})",
      "{min}, {max}"
    ),
    missing = "no"
  ) |> 
  add_n() |> # add column with total number of non-missing observations
  modify_header(label = "**Variable**") |> # update the column header
  bold_labels() -> s2
as_hux_table(s2) %>% huxtable::quick_xlsx(file = "CALM05_A_sample.xlsx")

c2_meta %>%
  group_by(PatientID) %>%
  slice_head(n = 1) %>%
  tbl_summary(
    include = colnames(c2_meta)[-c(1,2,3)],
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c(
      "{N_nonmiss}",
      "{median} ({p25}, {p75})",
      "{min}, {max}"
    ),
    missing = "no"
  ) |> 
  add_n() |> # add column with total number of non-missing observations
  modify_header(label = "**Variable**") |> # update the column header
  bold_labels() -> s2p
as_hux_table(s2p) %>% huxtable::quick_xlsx(file = "CALM05_A_patient.xlsx")