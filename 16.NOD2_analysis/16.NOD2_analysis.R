xbox::chdir("../16.NOD2_analysis/")
require(tidyverse)
rm(list = ls())
set.seed("20240724")
#c1------------------
kegg_c1 = read_tsv("CALM05_1_KEGG.Pathway/KEGG.Pathway.raw.xls")
#c2------------------
kegg_c2 = read_tsv("CALM05_2_KEGG.Pathway/KEGG.Pathway.raw.xls")

idtable  = read.table("../raw_data/cohort2_meta_16s_v1.xls",sep = "\t",header = 1,row.names = 1)
kegg_c2 = kegg_c2[,c("Pathway",rownames(idtable))]

colnames(kegg_c2)[-1] = idtable$SampleID
#merge c1 c2---------
kegg_m =  kegg_c1 %>% inner_join(kegg_c2,by = "Pathway")

colnames(kegg_m)
kegg_m[,-c(1,2)] = apply(kegg_m[,-c(1,2)], 2, function(x) x / sum(x))



library(tidyverse)
library(lmerTest)
library(broom.mixed)

c1_fs = read.csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_1_A.sum.csv",row.names = 1)
all_meta = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "ID"))

c1_fs = read.csv("../4.bac_slow_fast//Sepsis_add_trend_cluster.meta.csv",row.names = 1)
c1_meta = read.csv("../1.data_tidy/CALM2005_2_A.sum.csv",row.names = 1)
all_meta2 = c1_fs %>% select(Sample_ID,type) %>% inner_join(c1_meta,by = c("Sample_ID" = "SampleID"))

all_meta %>% rename(
  "Weight" = "Birth_weight",
  "Gender" = "sex",
  "RDS" = "NRDS",
  "Antibiotics_Duration_all" = "Antibiotic_duration",
  "Intrauterine_infection" = "Infection") %>% 
  mutate(Feeding = "A") %>%
  mutate(Natural_delivery = 
           case_when(
             Delivery == "S" ~ "N",
             Delivery == "E" ~ "Y")) %>%
  mutate(Weight = Weight*1000)-> all_meta

nid = intersect(colnames(all_meta),colnames(all_meta2))
nid = c(nid[1:2],nid[304:324])

all_meta = rbind(all_meta[,nid],all_meta2[,nid])
data = all_meta %>% filter(Age_DOL <=42)

data = data %>% filter(Age_DOL < Onset_Day | is.na(Onset_Day))

bac_i =  read.csv("../11.bac_index/bac_index_risk.csv",row.names = 1,header = T)

meta = data %>% left_join(bac_i %>% select(SampleID,score,Sepsis_sample,risk),by = c("Sample_ID" = "SampleID"))

iid = intersect(meta$Sample_ID,colnames(kegg_m))

kegg = kegg_m[,c("Pathway",iid)]

data = kegg %>% sjmisc::rotate_df(cn = T, rn = "Sample_ID") %>%
  left_join(meta,by = "Sample_ID")

kegg_test = kegg %>%
  pivot_longer(-Pathway,names_to = "SampleID",values_to = "exp") %>%
  left_join(meta,by = c("SampleID" = "Sample_ID"))


require(rstatix)

kegg_test %>%
  rstatix::group_by(Pathway) %>%
  rstatix::wilcox_test(exp ~ type) -> kegg_test_df_7d

kegg_test_df_7d = kegg_test_df_7d %>% mutate(p.adj = p.adjust(p, method = "fdr"))

fold_change <- kegg_test %>%
  group_by(Pathway,type) %>%
  summarise(median_value = mean(exp), .groups = 'drop') %>%
  drop_na() %>%
  group_by(Pathway) %>%
  summarise(fold_change = 
              (median_value[type == "fast"]) / (median_value[type == "slow"])) %>%
  mutate(log2fc= log2(fold_change))


kegg_test_a = kegg_test_df_7d %>% left_join(fold_change,by = "Pathway")

kegg_test_a$p.adj[kegg_test_a$p.adj > 0.3] = 0.3
kegg_test_a$fold_change[kegg_test_a$fold_change > 2] = 2

ggplot() +
  geom_point(data = kegg_test_a %>% 
               filter(p.adj > 0.05 | fold_change < 1.05) %>%
               filter(p.adj > 0.05 | fold_change > 0.95)  , aes(x = fold_change, y = -log10(p.adj)), color = "grey") +
  geom_point(data = kegg_test_a %>% filter(p.adj < 0.05 & fold_change > 1.05), 
             aes(x = fold_change, y = -log10(p.adj),size = -log10(p.adj),
                 alpha = -log10(p.adj)), 
             color = "#ff8bbf") +
  geom_point(data = kegg_test_a %>% filter(p.adj < 0.05 & fold_change < 0.95), 
             aes(x = fold_change, y = -log10(p.adj),size = -log10(p.adj),
                 alpha = -log10(p.adj)) , 
             color =  "#a0a5e1") +
  geom_text(size = 1,data = kegg_test_a %>% 
                             filter(p.adj < 0.01 | fold_change > 1.2) %>% 
                             filter(p.adj < 0.05), aes(x = fold_change, y = -log10(p.adj), label = Pathway)) +
  ggrepel::geom_text_repel(size = 1,data = kegg_test_a %>% 
                             filter(p.adj < 0.01 | fold_change < 0.8)%>% 
                             filter(p.adj < 0.05), aes(x = fold_change, y = -log10(p.adj), label = Pathway)) +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed") +
  geom_vline(xintercept = c(0.95, 1.05), color = "grey",linetype = "dashed") +
  coord_cartesian(ylim = c(0.5,2.2))+
  labs(x = "Fold Change(Fast/Slow)", y = "-log10(FDR)") +
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("vlcon_plot_diff_fastslow.pdf",width = 5.5,height = 4)
write.csv(kegg_test_a,"vlcon_plot_diff_fastslow.csv")

top_d = rbind(kegg_test_a %>% filter(p.adj < 0.015 & fold_change > 1.003),
              kegg_test_a %>% filter(p.adj < 0.015 & fold_change < 0.997))


kegg_test %>% mutate(
  Age_Group = case_when(
    Age_DOL <= 3 ~ "0d-3d",
    Age_DOL <= 7 ~ "4d-7d",
    Age_DOL <= 14 ~ "8d-14d",
    Age_DOL <= 21 ~ "15d-21d",
    Age_DOL <= 28 ~ "22d-28d",
    Age_DOL >= 29 ~ "29d-42d",
  )
) -> kegg_test


data2 = kegg_test %>% 
  filter(grepl("Peptidoglycan biosynthesis and degradation proteins",Pathway))
fit = lmer(exp ~ type * Age_DOL + GA_daily + Weight + Gender + Feeding + (1|PatientID), data = data2)
res = tidy(fit)

label = paste0("Type: P=", round(res$p.value[2], 3), "\n Age_DOL: P=", round(res$p.value[3], 3), "\nType * Age_DOL: P=", round(res$p.value[9], 3))

ggplot(data2, aes(Age_DOL, exp, color = type, fill = type)) +
    geom_smooth(alpha = 0.2) +
    annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) +
    ggthemes::theme_few()+
    scale_color_manual(values = c("#ff8bbf","#a0a5e1")) +
    scale_fill_manual(values = c("#ff8bbf","#a0a5e1")) +
    ggthemes::theme_clean()+
    theme(panel.grid.major.y = element_blank())+
    labs(x = "Day", y = "Peptidoglycan biosynthesis \nand degradation proteins")

ggsave("smooth_plot_diff_fastslow.pdf",width = 4.5,height = 4)

require(tidyverse)
results <- tibble::tibble(Pathway = character(),
                          Age_Group = character(),
                          fold_change = numeric(),
                          p_value = numeric(),
                          significance = character())

pathway_age_groups <- unique(kegg_test %>% select(Pathway, Age_Group))

for (i in 1:nrow(pathway_age_groups)) {
  pathway <- pathway_age_groups$Pathway[i]
  age_group <- pathway_age_groups$Age_Group[i]
  subset_data <- kegg_test %>% filter(Pathway == pathway, Age_Group == age_group)
  fast_count <- sum(subset_data$type == "fast",na.rm = T)
  slow_count <- sum(subset_data$type == "slow",na.rm = T)
  
  if (fast_count > 1 & slow_count > 1) {
    tryCatch({
      median_fast <- mean(subset_data$exp[subset_data$type == "fast"],na.rm = T)
      median_slow <- mean(subset_data$exp[subset_data$type == "slow"],na.rm = T)
      
      fold_change <- median_fast / median_slow
      
      wilcox_test_result <- wilcox.test(exp ~ type, data = subset_data)

      results <- results %>%
        add_row(Pathway = pathway,
                Age_Group = age_group,
                fold_change = fold_change,
                p_value = wilcox_test_result$p.value,
                significance = ifelse(wilcox_test_result$p.value < 0.05, "significant", "not significant"))
    }, error = function(e) {
      # 在捕获错误的情况下记录NA值
      results <- results %>%
        add_row(Pathway = pathway,
                Age_Group = age_group,
                fold_change = NA,
                p_value = NA,
                significance = NA)
    })
  } else {
    results <- results %>%
      add_row(Pathway = pathway,
              Age_Group = age_group,
              fold_change = NA,
              p_value = NA,
              significance = NA)
  }
}


results1 = results %>% mutate(p.adj = p.adjust(p_value, method = "fdr")) %>%
  mutate(p_value = round(p_value,2)) %>%
  mutate(signif = ifelse(p_value <= 0.05,"*",""))
write.csv(results1,"wilcox_diff_type_time_group.csv")

#----------
results_p = results1%>% filter(grepl("Peptidoglycan biosynthesis and degradation proteins",Pathway))
results_p$Age_Group = factor(results_p$Age_Group,
                             levels = rev(c( "0d-3d", "4d-7d","8d-14d","15d-21d" ,"22d-28d", 
                                         "29d-42d")))
results_p$fold_change = as.numeric(results_p$fold_change)

ggplot(results_p %>% filter(Age_Group != "0d-3d"), aes(y = fold_change, x = Age_Group)) +
  geom_line(aes(group = 1),color = "#a0a5e1",size = 2,alpha = 0.9)+
  geom_linerange(aes(ymin = 1, ymax = fold_change),linetype = "dashed")+
  geom_point(color = "#ff8bbf",size = 3,alpha = 0.9)+
  geom_text(aes(label = paste0("P = ",p_value)),hjust = 1.5)+
  geom_text(aes(label = signif),hjust = -1,color = "red",size = 10)+
  labs(x = "",y = "Peptidoglycan biosynthesis and degradationabundance\n(Fast/Slow)")+
  coord_flip()+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())

ggsave("peptidog_fc_timegroup.pdf",width = 3,height = 5)


# Metagenome data------------

all_gene = read_tsv("../16.NOD2_analysis/all_gene_result_out.tsv")
ko_s = read.csv("../16.NOD2_analysis/peptidog_s.csv")

sgene = all_gene %>% 
  inner_join(ko_s, by = c("KO" = "KO.Number"))
rm(all_gene)

write.csv(sgene,"peptidog_kegg_pathway.csv")

sgene$Gene = sub("[;,].*","",sgene$Gene.Name)

sgene %>% 
  dplyr::group_by(Family,Calss,taxonomy) %>%
  summarize(across(starts_with("S"),sum,.names = "{col}")) -> gene_sum

gene_pbp = gene_sum %>% filter(Family == "Class B Penicillin-Binding Protein")

nid  = intersect(meta$Sample_ID,colnames(gene_pbp))
gene_pbp$sum = apply(gene_pbp[,nid],1,sum)
gene_pbp %>% select(Family,Calss,taxonomy,sum) %>%
  group_by(taxonomy) %>%
  summarise(sum_sum = sum(sum)) %>%
  arrange(desc(sum_sum)) -> test

sgene %>% 
  dplyr::group_by(Family,Calss) %>%
  summarize(across(starts_with("S"),sum,.names = "{col}")) -> gene_sum

gene = gene_sum %>% mutate(gene = Family) %>%
  column_to_rownames("gene")
meta = kegg_test %>% select(-c(Pathway,exp)) %>% distinct()
nid  = intersect(meta$SampleID,colnames(gene))

meta_annot = meta %>% column_to_rownames("SampleID")
meta_annot = meta_annot[nid,] %>% arrange(type,Age_DOL)

gene_annot = gene[,1:2] %>% arrange(rev(Calss))

gene_exp = gene[rownames(gene_annot),rownames(meta_annot)]

redgreen = colorRampPalette(c("blue","#41ab5d","#f1fb67","#e4007c","red"),bias = 1)(200)
pheatmap::pheatmap(log10(gene_exp+0.0001),
                   annotation_col = meta_annot[,c(1,3,4)],
                   annotation_row = gene_annot,
                   cluster_cols = F,
                   cluster_rows = F,scale = "none",
                   color = redgreen,show_colnames = F)
dev.off()
write.csv(gene_exp,"peptidog_gene_allsample.heatmap.csv")


setwd("../16.NOD2_analysis/")

require(tidyverse)
gene_exp =  read.csv("peptidog_gene_allsample.heatmap.csv",row.names = 1)
require(ggcor)
quickcor(t(gene_exp),cor.test = TRUE) +
  geom_square(data = get_data(type = "lower", show.diag = FALSE)) +
  geom_mark(data = get_data(type = "upper", show.diag = FALSE), size = 2.5) +
  geom_abline(slope = -1, intercept = 14)+
  scale_fill_gradient2(low = "#41ab5d",mid = "#f1fb67",high = "#e4007c")

ggsave("gene_cor_heatmap.pdf",width = 8,height = 6)
#gene_exp = read.csv("peptidog_gene_allsample.heatmap.csv",row.names = 1)

gene_exp1 = gene_exp[apply(gene_exp > 0 ,1,sum) > 30,]

heat_data = gene_exp1 %>% rownames_to_column("gene") %>%
  pivot_longer(-gene,values_to = "exp",names_to = "SampleID") %>%
  left_join(meta_annot %>% rownames_to_column("SampleID"),by = "SampleID")
require(rstatix)
heat_data1_p = heat_data %>% 
  rstatix::group_by(gene) %>%
  rstatix::wilcox_test(exp ~ type) %>%
  add_significance()

write.csv(heat_data1_p,"peptidog_gene.heatmap_diff_p.csv")

gene_p = heat_data1_p %>% select(gene,p.signif) %>% column_to_rownames("gene")
gene_annot = cbind(gene_annot,gene_p[rownames(gene_annot),])

heat_data1 = heat_data %>% group_by(gene,type) %>%
  summarise(exp = mean(exp,na.rm =T)) %>%
  arrange(gene) %>%
  pivot_wider(names_from = gene,values_from = exp) %>%
  column_to_rownames("type")

gene_annot = gene_annot[colnames(heat_data1),] %>% arrange(Calss)

redgreen = colorRampPalette(c("blue","#41ab5d","#f1fb67","#e4007c"),bias = 0.5)(100)
redgreen = colorRampPalette(c("#41ab5d","#f1fb67"),bias = 1)(5)
pheatmap::pheatmap(scale(heat_data1[,rownames(gene_annot)]+0.0001),
                   annotation_col = gene_annot,
                   cellwidth = 15,cellheight = 30,
                   cluster_cols = F,cluster_rows = F,
                   border_color = T,
                   color = redgreen,
                   angle_col = 45,filename = "gene_ball_plot_annot.pdf")


heat_data %>% group_by(gene,type) %>%
  summarise(exp = mean(exp,na.rm =T)) -> ball_data

ball_data %>% 
  pivot_wider(names_from = type,values_from = exp) %>%
  group_by(gene) %>%
  summarise(fc = fast/slow) %>%
  mutate(type =  case_when(
    fc >= 1 ~ "fast",
    fc < 1 ~ "slow"
  ))-> ball_data

ball_data = ball_data %>% left_join(gene_annot %>% rownames_to_column("gene"),by= "gene")

ball_data %>% arrange(Calss,fc) -> ball_data

# output 

write.csv(ball_data,"ball_data.data.csv")
write.csv(heat_data1_p,"ball_data.p.csv")

setwd("../16.NOD2_analysis/")
ball_data = read.csv("ball_data.data.csv",row.names = 1)
heat_data1_p = read.csv("ball_data.p.csv",row.names = 1)

ball_data$gene = factor(ball_data$gene,levels= rev(ball_data$gene))
heat_data1_p$gene = factor(heat_data1_p$gene,levels= rev(ball_data$gene))

ggplot(ball_data %>% mutate(fc = 
  case_when(
    fc >= 1 ~ fc + 0.1,
    fc < 1 ~ fc - 0.1
  )
  ),aes(x = gene, y = fc)) + 
  geom_linerange(aes(color = type,ymin = 1,ymax= fc),
                 size = 1,alpha = 0.8) +
  geom_tile(data = heat_data1_p,aes(x = gene, y = 1,fill = p),width = 1, height = 0.2,color = "black")+
  geom_point(aes(color = type),size = 5,alpha = 1) +
  geom_text(data = heat_data1_p,aes(x = gene, y = 1,label = p.signif))+
  scale_fill_gradient2(low = "#8d0202",high = "white",midpoint = 0.1)+
  scale_y_continuous(breaks = c(seq(1,1.8,0.2)+0.1,seq(0.2,1,0.2)-0.1),
                     labels = c(seq(1,1.8,0.2),seq(0.2,1,0.2)))+
  scale_color_manual(values = c("#ff8bbf","#a0a5e1"))+
  ggthemes::theme_few()+
  theme(axis.text.x = element_text(angle = 90))

ggsave("gene_ball_plot_merge.pdf",width = 6,height = 6)


mnod = read_tsv("NOD2.shortbred.count.txt")
mnod = mnod[apply(mnod[,-1],1,sum) != 0,]
ndata = data.frame(apply(mnod[,-1],2,sum)) 
ndata = ndata %>% rownames_to_column("SampleID")
colnames(ndata) = c("SampleID","nod2")


meta = kegg_test %>% select(-c(Pathway,exp)) %>% distinct()
ndata1  = ndata %>% inner_join(meta,by = c("SampleID" = "Sample_ID"))

wirte.csv(ndata1,"DLpepiase_metagenome_sample.csv")

ndata2 = ndata1 

fit = lmer(nod2 ~ Age_DOL + GA_daily + Weight + Gender + Feeding + Antibiotics_Duration_2w+(1|PatientID), data = ndata2)
res = tidy(fit)
label = paste0("\nLMM adjust covariate\nAge_DOL: P=", round(res$p.value[2], 3))

ggplot() + 
  geom_smooth(data = ndata2 , 
              aes(Age_DOL, log10(nod2+0.1)),alpha = 0.2,
              method = "loess", 
              span = 0.95,color = "#1B8574",fill= "#1B8574") + 
  annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())+
  labs(x = "Age DOL", y = "NOD2")

ggsave("NOD2_exp_smooth.pdf",width = 4,height = 3)

ndata1 %>% mutate(
  Age_Group = case_when(
    Age_DOL <= 7 ~ "0d-7d",
    Age_DOL <= 14 ~ "8d-14d",
    Age_DOL <= 21 ~ "15d-21d",
    Age_DOL <= 28 ~ "22d-28d",
    Age_DOL >= 29 ~ "29d-42d",
  )
) -> ndata1


ndata1$Age_Group = factor(ndata1$Age_Group,
                            levels = c( "0d-7d","8d-14d","15d-21d" ,"22d-28d", 
                                        "29d-42d"))
df_p_temp <- ndata1%>% 
  rstatix::group_by(Age_Group) %>% 
  rstatix::wilcox_test(nod2 ~ type) %>% 
  rstatix::add_xy_position() %>%
  rstatix::add_significance(p.col = "p") %>%
  mutate(p = trunc(p*10000)/10000)

write_csv(df_p_temp,"Age_group_score_diff_boxplot.csv")

ggplot(ndata1,aes(x = Age_Group, y = nod2))+
  geom_boxplot(aes(fill = type),alpha = 0.5,width = 0.5,outlier.shape = NA,outlier.stroke = NA) +
  geom_smooth(aes(group = type,color = type), se = FALSE)+
  geom_text(data = df_p_temp , inherit.aes = F,aes(x = Age_Group,
                                                   y = y.position/8, 
                                                   label = p),color = "red")+
  coord_cartesian(ylim = c(0,200))+
  scale_color_manual(values = c("#ff8bbf","#a0a5e1"))+
  scale_fill_manual(values = c("#ff8bbf","#a0a5e1"))+
  ggthemes::theme_clean()+
  theme(panel.grid.major.y = element_blank())
ggsave("Age_group_NOD2_diff_boxplot.pdf",width = 5,height = 3)

setwd("../16.NOD2_analysis/")
mnod = read_tsv("NOD2.shortbred.count.txt")

mnod = mnod[apply(mnod[,-1],1,sum) != 0,]

nod_info = read_csv("DLendo_class_new.csv")

mnod_sum = data.frame(mnod[,1],apply(mnod[,-1],1,sum))
colnames(mnod_sum) = c("ID","sum")

mnod_sum = mnod_sum %>% left_join(nod_info,by = "ID")

mnod_sum = mnod_sum %>% arrange(sum)

mnod_sum$ID = factor(mnod_sum$ID,levels = mnod_sum$ID)
require(patchwork)

colpal = readRDS("../colpal.rds")
ggplot(mnod_sum %>% filter(!is.na(genus)),
       aes(x = "", y = sum, fill = genus)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +  
  scale_fill_manual(values = xbox::need_colors(16)) +  
  theme_classic() +
  theme(axis.line = element_blank(),  
        axis.text = element_blank(),  
        axis.ticks = element_blank(), 
        axis.title = element_blank()) +  
  labs(fill = "Genus")  
ggsave("NOD2_source_bac.pdf",width = 4,height = 4)

mnod_sum %>% filter(!is.na(genus)) %>% mutate(r = 
  sum/sum(mnod_sum$sum)) %>% write.csv("NOD2_source_bac.csv")

vdata = read.csv("NOD2_source_bac.csv",row.names = 1)

colnames(vdata)

vdata_plot = vdata %>% select(ID,genus,r)

library(ggalluvial)

ggplot(vdata_plot, aes(axis1 = "Start", axis3 = genus,axis2 = ID, y = r)) +
  geom_alluvium(fill = "#d9d9d8",color = "black") +
  geom_stratum(aes(fill = genus),width = 0.6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = cp)+
  theme_minimal() +
  labs(title = "Sankey Diagram using ggalluvial", x = "", y = "Normalized Weight")

ggsave("NOD2_source_bac_sankey.pdf",width = 8,height = 4)