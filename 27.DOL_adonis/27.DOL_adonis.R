# a diversity
xbox::chdir("./27.DOL_adonis")
rm(list = ls())
require(tidyverse)
require(vegan)

source("../addin.r")
file_up<-list.files("../1.data_tidy/",pattern = "*.sum.csv")

for (i in 1:length(file_up)) {
  a<-read.csv(file.path("../1.data_tidy/",file_up[i]))
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

data_exp=CALM2005_1_A
otu_exp1=data_exp[,2:450]
dol = data_exp["Age_DOL"]
meio.mdf <- as.matrix.data.frame(otu_exp1)
rownames(meio.mdf) <- rownames(otu_exp1)
meio.mdf_clean <- meio.mdf[rowSums(meio.mdf) != 0, ]
arg_distance <- vegdist(meio.mdf_clean, method = "bray")
anova_result = adonis3(as.formula(paste("arg_distance","~","Age_DOL")),
                       data = data_exp,
                       permutations = 999,
                       na.action=na.omit)
anova_temp1 = data.frame(factor = i , Pvalue = anova_result[1,5],R2 = round(anova_result[1,3], 3))
anova_temp1$cohort = "CALM2005_1_A" 




data_exp=CALM2005_2_A
otu_exp1=data_exp[,2:497]
dol = data_exp["Age_DOL"]
meio.mdf <- as.matrix.data.frame(otu_exp1)
rownames(meio.mdf) <- rownames(otu_exp1)
meio.mdf_clean <- meio.mdf[rowSums(meio.mdf) != 0, ]
arg_distance <- vegdist(meio.mdf_clean, method = "bray")
anova_result = adonis3(as.formula(paste("arg_distance","~","Age_DOL")),
                       data = data_exp,
                       permutations = 999,
                       na.action=na.omit)
anova_temp2 = data.frame(factor = i , Pvalue = anova_result[1,5],R2 = round(anova_result[1,3], 3))
anova_temp2$cohort = "CALM2005_2_A" 


data_exp=Rao_2021_A
otu_exp1=data_exp[,2:237]
dol = data_exp["Age_DOL"]
meio.mdf <- as.matrix.data.frame(otu_exp1)
rownames(meio.mdf) <- rownames(otu_exp1)
meio.mdf_clean <- meio.mdf[rowSums(meio.mdf) != 0, ]
arg_distance <- vegdist(meio.mdf_clean, method = "bray")
anova_result = adonis3(as.formula(paste("arg_distance","~","Age_DOL")),
                       data = data_exp,
                       permutations = 999,
                       na.action=na.omit)
anova_temp3 = data.frame(factor = i , Pvalue = anova_result[1,5],R2 = round(anova_result[1,3], 3))
anova_temp3$cohort = "Rao_2021_A" 


data_exp=Lauren_2022_M
otu_exp1=data_exp[,2:134]
dol = data_exp["Age_DOL"]
meio.mdf <- as.matrix.data.frame(otu_exp1)
rownames(meio.mdf) <- rownames(otu_exp1)
meio.mdf_clean <- meio.mdf[rowSums(meio.mdf) != 0, ]
arg_distance <- vegdist(meio.mdf_clean, method = "bray")
anova_result = adonis3(as.formula(paste("arg_distance","~","Age_DOL")),
                       data = data_exp,
                       permutations = 999,
                       na.action=na.omit)
anova_temp4 = data.frame(factor = i , Pvalue = anova_result[1,5],R2 = round(anova_result[1,3], 3))
anova_temp4$cohort = "Lauren_2022_M" 


data_exp=Brooks_2017_M
otu_exp1=data_exp[,2:102]
dol = data_exp["Age_DOL"]
meio.mdf <- as.matrix.data.frame(otu_exp1)
rownames(meio.mdf) <- rownames(otu_exp1)
meio.mdf_clean <- meio.mdf[rowSums(meio.mdf) != 0, ]
arg_distance <- vegdist(meio.mdf_clean, method = "bray")
anova_result = adonis3(as.formula(paste("arg_distance","~","Age_DOL")),
                       data = data_exp,
                       permutations = 999,
                       na.action=na.omit)
anova_temp5 = data.frame(factor = i , Pvalue = anova_result[1,5],R2 = round(anova_result[1,3], 3))
anova_temp5$cohort = "Brooks_2017_M"

data_exp=Olm_2019_M
otu_exp1=data_exp[,2:102]
dol = data_exp["Age_DOL"]
meio.mdf <- as.matrix.data.frame(otu_exp1)
rownames(meio.mdf) <- rownames(otu_exp1)
meio.mdf_clean <- meio.mdf[rowSums(meio.mdf) != 0, ]
arg_distance <- vegdist(meio.mdf_clean, method = "bray")
anova_result = adonis3(as.formula(paste("arg_distance","~","Age_DOL")),
                       data = data_exp,
                       permutations = 999,
                       na.action=na.omit)
anova_temp6 = data.frame(factor = i , Pvalue = anova_result[1,5],R2 = round(anova_result[1,3], 3))
anova_temp6$cohort = "Olm_2019_M"

data_exp=Raveh_2015_M
otu_exp1=data_exp[,2:102]
dol = data_exp["Age_DOL"]
meio.mdf <- as.matrix.data.frame(otu_exp1)
rownames(meio.mdf) <- rownames(otu_exp1)
meio.mdf_clean <- meio.mdf[rowSums(meio.mdf) != 0, ]
arg_distance <- vegdist(meio.mdf_clean, method = "bray")
anova_result = adonis3(as.formula(paste("arg_distance","~","Age_DOL")),
                       data = data_exp,
                       permutations = 999,
                       na.action=na.omit)
anova_temp7 = data.frame(factor = i , Pvalue = anova_result[1,5],R2 = round(anova_result[1,3], 3))
anova_temp7$cohort = "Raveh_2015_M"



data_exp=Robert_2024_M
otu_exp1=data_exp[,2:132]
dol = data_exp["Age_DOL"]
meio.mdf <- as.matrix.data.frame(otu_exp1)
rownames(meio.mdf) <- rownames(otu_exp1)
meio.mdf_clean <- meio.mdf[rowSums(meio.mdf) != 0, ]
arg_distance <- vegdist(meio.mdf_clean, method = "bray")
anova_result = adonis3(as.formula(paste("arg_distance","~","Age_DOL")),
                       data = data_exp,
                       permutations = 999,
                       na.action=na.omit)
anova_temp8 = data.frame(factor = i , Pvalue = anova_result[1,5],R2 = round(anova_result[1,3], 3))
anova_temp8$cohort = "Robert_2024_M"


anova_data = rbind(anova_temp1,anova_temp2,anova_temp3,anova_temp4,
                   anova_temp5,anova_temp6,anova_temp7,anova_temp8)

write.csv(anova_data,"dol_adonis_cohort.csv")