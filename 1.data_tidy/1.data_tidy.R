# split every cohort for data tidy and then merge

xbox::chdir("1.data_tidy")
# library R packages ---------------
require(tidyverse)
library(phyloseq)
require(microbiomeutilities)
require(ggpubr)
require(microbiome)
require(vegan)
rm(list = ls())

# VMT-------------
amp_data = read_tsv("../raw_data/VMT/dada2.xls")
colnames(amp_data)[1] = "ID"
amp_meta =  read_tsv("../raw_data/VMT/infant-cohort_metadata-microbiome.xls")
colnames(amp_meta)[1] = "ID"

amp_meta = amp_meta %>% filter(BodySite == "Stool")  %>% filter(Species == "Infant")

amp_data = amp_data %>% column_to_rownames("ID")
amp_data = amp_data[,c(amp_meta$ID,"taxonomy")]

# out
otumat = as.matrix(amp_data[,c(1:(ncol(amp_data) - 1))])
# tax
taxmat = amp_data[,c(1,ncol(amp_data))] 
taxmat$taxonomy = gsub("[A-Za-z]__","",taxmat$taxonomy)
taxmat = taxmat %>% separate(taxonomy, sep = "; ",into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% select(-1)
taxmat[taxmat == " "] = NA
taxmat[taxmat == ""] = NA

taxmat = as.matrix(taxmat)

# meta
group_data = amp_meta
rownames(group_data) = group_data$ID
# convert
OTU = otu_table(otumat, taxa_are_rows = T)
TAX = tax_table(taxmat)
sampledata = sample_data(group_data)
rownames(sampledata) = sampledata$ID
physeq = phyloseq(OTU, TAX, sampledata)
saveRDS(physeq,"clean_data/vmt_infant_stools_16s_origin_obj.rds")

# 输出一个汇总的表格

c_levels = "Genus" #"Family"
ps.g = tax_glom(physeq, c_levels)
ps.g = transform_sample_counts(ps.g, function(OTU) OTU/sum(OTU))

otu_data = otu_tibble(ps.g) %>% 
  left_join(tax_tibble(ps.g),.,by = "FeatureID") %>%
  select(-c("FeatureID")) %>%
  select(-c("Kingdom","Phylum","Class","Order","Family","Species"))

otu_data %>% 
  dplyr::group_by(Genus) %>%
  mutate(across(.cols = everything(), .fns = ~sum(., na.rm = TRUE))) %>%
  distinct() -> otu_data

otu_data1 = otu_data %>% sjmisc::rotate_df(cn = T, rn = "SampleID")
meta_data =  sample_data(ps.g)

otu_data2 = otu_data1 %>% inner_join(meta_data,by = c("SampleID" = "ID"))
otu_data2$cohort = "vmt_infant_stools_16s"
write.csv(otu_data2,"vmt_infant_stools_A.sum.csv")

rm(list = ls())

# CALM05 B------------
amp_data = read_tsv("../raw_data/cohort1.bac.exp.xls")
colnames(amp_data)[1] = "ID"
amp_meta =  read_tsv("../raw_data/cohort1.bac.meta.xls")
colnames(amp_meta)[1] = "ID"

amp_data = amp_data %>% column_to_rownames("ID")
amp_data = amp_data[,c(amp_meta$ID,"taxonomy")]

# out
otumat = as.matrix(amp_data[,c(1:(ncol(amp_data) - 1))])
# tax
taxmat = amp_data[,c(1,ncol(amp_data))] 
taxmat$taxonomy = gsub("[A-Za-z]__","",taxmat$taxonomy)
taxmat = taxmat %>% separate(taxonomy, sep = "; ",into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% select(-1)
taxmat[taxmat == " "] = NA
taxmat[taxmat == ""] = NA

taxmat = as.matrix(taxmat)

# meta
group_data = amp_meta
rownames(group_data) = group_data$ID
# convert
OTU = otu_table(otumat, taxa_are_rows = T)
TAX = tax_table(taxmat)
sampledata = sample_data(group_data)
rownames(sampledata) = sampledata$ID
physeq = phyloseq(OTU, TAX, sampledata)
saveRDS(physeq,"cohort1_origin_obj.rds")

c_levels = "Genus" #"Family"
ps.g = tax_glom(physeq, c_levels)
ps.g = transform_sample_counts(ps.g, function(OTU) OTU/sum(OTU))

otu_data = otu_tibble(ps.g) %>% 
  left_join(tax_tibble(ps.g),.,by = "FeatureID") %>%
  select(-c("FeatureID")) %>%
  select(-c("Kingdom","Phylum","Class","Order","Family","Species"))

otu_data %>% 
  dplyr::group_by(Genus) %>%
  mutate(across(.cols = everything(), .fns = ~sum(., na.rm = TRUE))) %>%
  distinct() -> otu_data

otu_data1 = otu_data %>% sjmisc::rotate_df(cn = T, rn = "SampleID")
meta_data =  sample_data(ps.g)

otu_data2 = otu_data1 %>% inner_join(meta_data,by = "SampleID")
otu_data2$cohort = "CALM2005_1_A"
write.csv(otu_data2,"CALM2005_1_A.sum.csv")

rm(list = ls())

# CALM05 A------------
amp_data = read_tsv("../raw_data/cohort2_16s_profile.xls")
colnames(amp_data)[1] = "ID"
amp_meta =  read_tsv("../raw_data/cohort2_meta_16s_v1.xls")
amp_tax = read_tsv("../raw_data/cohort2_features_table.txt")

amp_data = amp_data %>% column_to_rownames("ID")
amp_data = amp_data[,amp_meta$ID]
colnames(amp_data) = amp_meta$SampleID
amp_data = amp_data %>% rownames_to_column("ID")
amp_data = amp_data %>% left_join(amp_tax,by = "ID")
amp_data = amp_data %>% column_to_rownames("ID")

# out
otumat = as.matrix(amp_data[,c(1:437)])
# tax
taxmat = amp_data[,c(1,438)] 
taxmat$taxonomy = gsub("[A-Za-z]__","",taxmat$taxonomy)
taxmat = taxmat %>% separate(taxonomy, sep = ";",into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% select(-1)
taxmat[taxmat == ""] = NA
taxmat = as.matrix(taxmat)
# meta
group_data = amp_meta[,-1]
rownames(group_data) = group_data$SampleID
# convert 
OTU = otu_table(otumat, taxa_are_rows = T)
TAX = tax_table(taxmat)
sampledata = sample_data(group_data)
rownames(sampledata) = sampledata$SampleID
physeq = phyloseq(OTU, TAX, sampledata)
saveRDS(physeq,"cohort2_origin_obj.rds")

c_levels = "Genus" #"Family"
ps.g = tax_glom(physeq, c_levels)
ps.g = transform_sample_counts(ps.g, function(OTU) OTU/sum(OTU))

otu_data = otu_tibble(ps.g) %>% 
  left_join(tax_tibble(ps.g),.,by = "FeatureID") %>%
  select(-c("FeatureID")) %>%
  select(-c("Kingdom","Phylum","Class","Order","Family","Species"))

otu_data %>% 
  dplyr::group_by(Genus) %>%
  mutate(across(.cols = everything(), .fns = ~sum(., na.rm = TRUE))) %>%
  distinct() -> otu_data

otu_data1 = otu_data %>% sjmisc::rotate_df(cn = T, rn = "SampleID")
meta_data =  sample_data(ps.g)

otu_data2 = otu_data1 %>% inner_join(meta_data,by = "SampleID")
otu_data2$cohort = "CALM2005_2_A"
write.csv(otu_data2,"CALM2005_2_A.sum.csv")

rm(list = ls())

# Rao_2021 data ---------------------

remove_leves = ""
remove_list = read.csv("../remove_bac.list",header = F)
remove_list = paste(remove_leves,as.vector(remove_list$V1),sep = "")
tax_data = read_csv("../raw_data/Rao_2021_A_tax.csv")
tax_data = data.frame(lapply(
  tax_data, function(x) gsub("\\(\\d+\\.\\d+\\)", "", x)
)
)
tax_data = tax_data[which(!Reduce(`|`, lapply(remove_list, grepl, tax_data$genus))),]
tax_data = data.frame(lapply(tax_data, trimws))

amp_data = read.csv("../raw_data/Rao_2021_A_exp.csv",row.names = 1,header = T)
# meta
amp_meta = amp_data[,c(811,812)]
# exp
amp_data = amp_data[,-c(811,812)]
amp_data = t(amp_data)
amp_data = as.data.frame(amp_data) %>% rownames_to_column("OTU_ID")
amp_data$OTU_ID = sub("^X","",amp_data$OTU_ID)

# meta
amp_data = amp_data %>% inner_join(tax_data,by = "OTU_ID")
# out
amp_data[,1:958] %>% column_to_rownames("OTU_ID") -> otumat
otumat = as.matrix(otumat)
# tax
taxmat = amp_data[,c(1,959:965)]
taxmat = taxmat %>% column_to_rownames("OTU_ID")
taxmat = as.matrix(taxmat)
# meta
group_data = amp_meta %>% rownames_to_column("ID")
colnames(group_data)[1] = "SampleID"
rownames(group_data) = group_data$SampleID
# convert 
OTU = otu_table(otumat, taxa_are_rows = T)
TAX = tax_table(taxmat)
sampledata = sample_data(group_data)
rownames(sampledata) = sampledata$SampleID
# tax_names
physeq = phyloseq(OTU, TAX, sampledata)
saveRDS(physeq,"Rao_2021_A_origin_obj.rds")

c_levels = "genus" #"Family"
ps.g = tax_glom(physeq, c_levels)
ps.g = transform_sample_counts(ps.g, function(OTU) OTU/sum(OTU))

otu_data = otu_tibble(ps.g) %>% 
  left_join(tax_tibble(ps.g),.,by = "FeatureID") %>%
  select(-c("FeatureID")) %>%
  select(-c("kingdom","phylum","class","order","family","species"))

otu_data %>% 
  dplyr::group_by(genus) %>%
  mutate(across(.cols = everything(), .fns = ~sum(., na.rm = TRUE))) %>%
  distinct() -> otu_data

otu_data1 = otu_data %>% sjmisc::rotate_df(cn = T, rn = "SampleID")
meta_data =  sample_data(ps.g)

otu_data2 = otu_data1 %>% inner_join(meta_data,by = c("SampleID"))
otu_data2$cohort = "Rao_2021_A"
write.csv(otu_data2,"Rao_2021_A.sum.csv")


# Lauren_2022_M---------------------
rm(list  = ls())
amp_data = read.table("../raw_data/Lauren_2022_M.csv",header = T,sep = ",")
amp_data = amp_data %>% filter(Probiotics_Summary %in% c("Before","Never"))

colnames(amp_data)
amp_data = amp_data[,1:419]
amp_meta = amp_data[,1:16]
amp_exp = amp_data[,c(1,26:419)]

amp_exp = amp_exp%>% sjmisc::rotate_df(cn = T,rn = "Genus")
amp_exp$Genus = sub("_.*","",amp_exp$Genus)

amp_exp %>% 
  dplyr::group_by(Genus) %>%
  mutate(across(.cols = everything(), .fns = ~sum(., na.rm = TRUE))) %>%
  distinct() -> otu_data
# 检查下
apply(otu_data[,-1],2,sum)

otu_data1 = otu_data %>% sjmisc::rotate_df(cn = T, rn = "SampleID")

otu_data2 = otu_data1 %>% inner_join(amp_meta,by = c("SampleID" = "SampleID"))

otu_data2$cohort = "Lauren_2022_M"
write.csv(otu_data2,"Lauren_2022_M.sum.csv")

rm(list = ls())

# Olm_2019_M---------------------
amp_data = read_csv("../raw_data/Olm_2019_M_exp.csv")
colnames(amp_data)[1] = "ID"
amp_meta =  read_csv("../raw_data/Olm_2019_M_meta.csv")

amp_data = amp_data %>% column_to_rownames("ID")
remove_leves = ""
remove_list = read.csv("../remove_bac.list",header = F)
remove_list = paste(remove_leves,as.vector(remove_list$V1),sep = "")

amp_data = amp_data[!amp_data$taxonomy %in% remove_list,]
amp_data1 = amp_data %>%
  rename("Genus" = "taxonomy")
amp_data1 %>%
  dplyr::group_by(Genus) %>%
  distinct() -> otu_data

otu_data = otu_data %>% filter(Genus != "unk")
otu_data = otu_data %>% column_to_rownames("Genus")
otu_data = t(otu_data)
otu_data1 = as.data.frame(otu_data/apply(otu_data,1,sum))

otu_data1  = otu_data1 %>% rownames_to_column("Sample_ID")

otu_data2 = otu_data1 %>% inner_join(amp_meta,by = c("Sample_ID"))

otu_data2$cohort = "Olm_2019_M"
write.csv(otu_data2,"Olm_2019_M.sum.csv")
# 分队列的研究
for(i in unique(otu_data2$campaign)){
  otu_data2 %>% filter(campaign == i) %>%
    mutate(cohort =paste0("Olm_2019_",i,"_M")) %>%
    write.csv(paste0("Olm_2019_",i,"_M.sum.csv"))
}


rm(list = ls())

# Robert_2024_M---------------------
amp_data = read_csv("../raw_data/Robert_2024_M_exp.csv")
amp_data1 =  read_csv("../raw_data/Robert_2024_M_exp_nec.csv")

id = intersect(colnames(amp_data),colnames(amp_data1))
amp_data = bind_rows(amp_data[,id],amp_data1[,id])
amp_data[!duplicated(amp_data$SampleID),] -> amp_data
colnames(amp_data)[1] = "ID"

amp_meta =  read_csv("../raw_data/Robert_2024_M_sample_meta.csv")
amp_meta1 =  read_csv("../raw_data/Robert_2024_M_sample_meta_nec.csv")
id = intersect(colnames(amp_meta),colnames(amp_meta1))
id = id[-12]
amp_meta = bind_rows(amp_meta[,id],amp_meta1[,id])
amp_meta[!duplicated(amp_meta$DNA_Sample_ID),] -> amp_meta
rm(amp_data1,amp_meta1)

amp_data = amp_data %>% column_to_rownames("ID")

amp_data1 = xbox::t_dat(amp_data)
amp_data1 = amp_data1 %>% rownames_to_column("taxonomy")
amp_data1$taxonomy = sub("_.*","",amp_data1$taxonomy)

amp_data1 = amp_data1[!amp_data1$taxonomy %in% remove_list,]

amp_data1 = amp_data1 %>%
  rename("Genus" = "taxonomy")

amp_data1 %>%
  dplyr::group_by(Genus) %>%
  mutate(across(.cols = everything(), .fns = ~sum(., na.rm = TRUE))) %>%
  distinct() -> otu_data

otu_data = otu_data %>% filter(Genus != "unk")
otu_data = otu_data %>% column_to_rownames("Genus")
otu_data = t(otu_data)
otu_data1 = as.data.frame(otu_data/apply(otu_data,1,sum))

otu_data1  = otu_data1 %>% rownames_to_column("Sample_ID")

colnames(amp_meta)
amp_meta %>% rename("Age_DOL" = "DOL",
                    "PatientID" = "Subject",
                    ) %>%
  mutate(Sepsis_style = NA,
         NEC = case_when(
           NEC_status == 1 ~ "Y",
           NEC_status == 0 ~ "N"
         )) -> amp_meta1

# 
otu_data2 = otu_data1 %>% inner_join(amp_meta1,
                                     by = c("Sample_ID" = "DNA_Sample_ID"))


otu_data2$cohort = "Robert_2024_M"

write.csv(otu_data2,"Robert_2024_M.sum.csv")

rm(list = ls())

# sum cohort -------------
rm(list = ls())
getwd()
setwd("../1.data_tidy/")

file_up<-list.files("./",pattern = "*.sum.csv")

for (i in 1:length(file_up)) {
  a<-read.csv(file_up[i])
  a = a[,-1]
  colnames(a)[1] = "Sample_ID"
  nn = unique(a$cohort)[1]
  assign(nn,a)
  rm(a)
}



# 选择对垒
# CALM2005_1_A,CALM2005_2_A,Lauren_2022_M,
# Rao_2021_A,Brooks_2017_M,Olm_2019_M,Raveh_2015_M

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

col_names_list <- lapply(list(CALM2005_1_A,
                              CALM2005_2_A,
                              Rao_2021_A,
                              Lauren_2022_M,
                              Brooks_2017_M,
                              Olm_2019_M,
                              Raveh_2015_M,
                              Robert_2024_M
                              ), colnames)

lapply(col_names_list, function(x) "Age_DOL" %in% x)

lapply(col_names_list, function(x) length(x))

common_col_names <- Reduce(intersect, col_names_list)

# 修改下ID的类型
Olm_2019_M$PatientID = as.character(Olm_2019_M$PatientID)
Brooks_2017_M$PatientID = as.character(Brooks_2017_M$PatientID)
Raveh_2015_M$PatientID = as.character(Raveh_2015_M$PatientID)
Rao_2021_A$PatientID = as.character(Rao_2021_A$PatientID)
Lauren_2022_M$PatientID = as.character(Lauren_2022_M$PatientID)
Robert_2024_M$PatientID = as.character(Robert_2024_M$PatientID)


require(tidyverse)
combined_df <- bind_rows(
  CALM2005_1_A[,common_col_names],
  CALM2005_2_A[,common_col_names],
  Rao_2021_A[,common_col_names],
  Lauren_2022_M[,common_col_names],
  Brooks_2017_M[,common_col_names],
  Olm_2019_M[,common_col_names],
  Raveh_2015_M[,common_col_names],
  Robert_2024_M[,common_col_names]
)

write.csv(as.data.frame(table(combined_df$cohort)),"all_cohort_sample_num_8cohort_0905.csv")

renorm_df = combined_df[,2:48]/apply(combined_df[,2:48],1,sum)
renorm_df[is.na(renorm_df)] = 0

renorm_df_model = cbind(combined_df$Sample_ID,renorm_df,combined_df[,49:ncol(combined_df)])

colnames(renorm_df_model) = sub("combined_df\\$","",colnames(renorm_df_model))

write.csv(combined_df,"all_cohort_sample_num_8cohort_0905.csv")