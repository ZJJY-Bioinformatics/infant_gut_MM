require(tidyverse)
require(ggprism)
require(rstatix)
setwd("28.Poilt_metagenome/")
rm(list = ls())
meta = read_csv("./meta_poilt_metagenome.csv")
colnames(meta)

meta %>% mutate(adjust_time = 
                  case_when(
                    Age >= Onset_Day ~ Age - Onset_Day + 1,
                    Age < Onset_Day ~ Age - Onset_Day
                  )) %>% filter(Sample_ID != "SD00084983") -> meta
write.csv(meta,"meta_poilt_metagenome_adjusttime.csv")

# DL-endo profile--------------

dl_exp = read_tsv("NOD2.shortbred.count.txt")

# WP_003669787 WP_003671744 WP_003670888

dl_exp %>% 
  filter(Family %in% c("WP_003669787","WP_003671744","WP_003670888")) %>%
  summarise(across(-Family, sum)) %>% 
  xbox::t_dat() %>% as.data.frame()-> dl_exp

dl_exp$Sample_ID = rownames(dl_exp)
colnames(dl_exp) = c("DLendo","Sample_ID")
meta %>% 
  inner_join(dl_exp,by = "Sample_ID") -> dl_plot

dl_plot$DLendo = as.numeric(dl_plot$DLendo)

dl_plot  %>% 
  filter(adjust_time >= -2) %>%
  mutate(on = case_when(
    adjust_time > 14 ~ "after",
    adjust_time >= 0 ~ "during",
    adjust_time < 0 ~ "before"
  )) %>%
  group_by(on) %>% 
  wilcox_test(DLendo ~ Group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", 
    p <= 0.01  ~ "**",  
    p <= 0.05  ~ "*",  
    TRUE       ~ "ns"   
  )) -> test.data

write.csv(test.data,"DL_exp.test.csv")


ggplot(dl_plot  %>% 
         filter(adjust_time >= -2) %>%
         mutate(on = case_when(
           adjust_time > 14 ~ "after",
           adjust_time >= 0 ~ "during",
           adjust_time < 0 ~ "before"
         )))+
  geom_boxplot(aes(x = on, y = DLendo,fill = Group))+
  scale_fill_manual(values = c("Untreated Control" = "#F4C573","Probiotic Intervention" = "#49c8cb"))+
  scale_color_manual(values = c("Untreated Control" = "#F4C573","Probiotic Intervention" = "#49c8cb"))+
  theme_prism(base_size = 20,base_fontface = "plain",
              base_line_size = 0.5) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  labs(y = "RPKM of DL-endopeptidasefrom\n(from L.reuteri DSM 17938)",
       x = "Days relative to intervention")
ggsave("DL_exp_boxplot.pdf",width = 11,height = 6)


dl_exp = read_tsv("NOD2.shortbred.count.txt")

# WP_003669787 WP_003671744 WP_003670888

dl_exp %>% 
  summarise(across(-Family, sum)) %>% 
  xbox::t_dat() %>% as.data.frame()-> dl_exp

dl_exp$Sample_ID = rownames(dl_exp)
colnames(dl_exp) = c("DLendo","Sample_ID")
meta %>% 
  inner_join(dl_exp,by = "Sample_ID") -> dl_plot

dl_plot$DLendo = as.numeric(dl_plot$DLendo)

dl_plot  %>% 
  filter(adjust_time >= -2) %>%
  mutate(on = case_when(
    adjust_time > 14 ~ "after",
    adjust_time >= 0 ~ "during",
    adjust_time < 0 ~ "before"
  )) %>%
  group_by(on) %>% 
  wilcox_test(DLendo ~ Group) %>%
  mutate(stars = case_when(
    p <= 0.001 ~ "***", 
    p <= 0.01  ~ "**", 
    p <= 0.05  ~ "*",  
    TRUE       ~ "ns"   
  )) -> test.data

write.csv(test.data,"total_DL_exp.test.csv")


ggplot(dl_plot  %>% 
         filter(adjust_time >= -2) %>%
         mutate(on = case_when(
           adjust_time > 14 ~ "after",
           adjust_time >= 0 ~ "during",
           adjust_time < 0 ~ "before"
         )))+
  geom_boxplot(aes(x = on, y = DLendo,fill = Group))+
  scale_fill_manual(values = c("Untreated Control" = "#F4C573","Probiotic Intervention" = "#49c8cb"))+
  scale_color_manual(values = c("Untreated Control" = "#F4C573","Probiotic Intervention" = "#49c8cb"))+
  theme_prism(base_size = 20,base_fontface = "plain",
              base_line_size = 0.5) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  labs(y = "RPKM of DL-endopeptidasefrom\n(Total Sample)",
       x = "Days relative to intervention")
ggsave("taotal_DL_exp_boxplot.pdf",width = 11,height = 6)

# metaphlan--------

mtp_exp = read_tsv("all.sample_buglist.tsv")
colnames(mtp_exp)[1] = "taxonomy"

taxmat = mtp_exp %>% filter(grepl("k__Bacteria",taxonomy))
taxmat$taxonomy = gsub("[A-Za-z]__","",taxmat$taxonomy)
taxmat = taxmat %>% separate(taxonomy, sep = "\\|",into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Strain"))
taxmat[taxmat == " "] = NA
taxmat[taxmat == ""] = NA
taxmat[taxmat == "unclassified"] = NA

sexp = taxmat %>% filter(!is.na(Strain))

colnames(sexp) = sub("_profiled_metagenome.txt","",colnames(sexp))

sexp$tax = paste0(sexp$Species,"_",sexp$Strain)
sexp = sexp[,-c(1:8)]
sexp = sexp %>% column_to_rownames("tax")

sexp = sexp %>% sjmisc::rotate_df(rn = "Sample_ID")
sexp = sexp %>% column_to_rownames("Sample_ID")
otu_data1 = as.data.frame(sexp/apply(sexp,1,sum))

apply(otu_data1,1,sum)

otu_data1 %>% rownames_to_column("Sample_ID") %>%
  select(Sample_ID,	ends_with("SGB7095")) -> bac_exp
colnames(bac_exp)[2] = "rpkm"
meta %>% 
  left_join(bac_exp,by = c("Sample_ID")) -> bac_plot

err_data_summary <- bac_plot  %>% 
  filter(adjust_time >= -2) %>% 
  mutate(on = case_when(
    adjust_time > 14 ~ "after",
    adjust_time >= 0 ~ "during",
    adjust_time < 0 ~ "before"
  )) %>%
  filter(on != "after") %>% 
  group_by(Group,on) %>%
  summarise(
    mean_n = mean(rpkm, na.rm = TRUE),
    sem_n = sd(rpkm, na.rm = TRUE) / sqrt(n())
  )

err_data_summary$on = ifelse(err_data_summary$on=="before",0,1)

ggplot() +
  geom_line(data = err_data_summary,
            aes(x = on, y = mean_n,color = Group,group = Group),
            linewidth = 1) +
  geom_point(data = err_data_summary,
             aes(x = on,y = mean_n,color = Group), size = 3) +
  geom_errorbar(data = err_data_summary,
                aes(x = on,ymin = mean_n - sem_n, 
                    ymax = mean_n + sem_n,color = Group), width = 0.05) +
  scale_fill_manual(values = c("Untreated Control" = "#F4C573","Probiotic Intervention" = "#49c8cb"))+
  scale_color_manual(values = c("Untreated Control" = "#F4C573","Probiotic Intervention" = "#49c8cb"))+
  scale_x_continuous(breaks = c(0,1),labels = c("Pre-Intervention","During-Intervention"))+
  scale_y_continuous(breaks = c(0,0,0.01,0.02,0.03),labels = c(0,0,0.01,0.02,0.03))+
  labs(x = "",
       y = "RPKM of L.reuteri DSM 17938")+
  theme_prism(base_size = 16,base_fontface = "plain",base_line_size = 1)

ggsave("LR_exp_line_metaphlan.pdf",width = 6,height = 5)


# cor with nod2
require(tidyverse)
nod2 = read_csv("./polit.csv") 
colnames(nod2)[13] = "VS"
bac_plot %>% left_join(nod2,by = c("Patient_ID")) ->  re1

re1 %>% mutate(s = case_when(
  abs(`Age.x` - `Age.y`) <= 1 ~ "Y"
)) %>% 
  filter(s == "Y") -> re2

re2 %>%
  select(Sample_ID.x,`Delta OD Value`) ->  re
re %>% drop_na() -> re

re$Sample_ID.x

bac_otu = otu_data1[re$Sample_ID.x,]

result = data.frame()
for(i in 1:ncol(bac_otu)){
  cor_re = cor.test(re$`Delta OD Value`,bac_otu[[i]],method = "spearman")
  temp = data.frame(P = cor_re[["p.value"]], 
                    cor_r = cor_re[["estimate"]][["rho"]],
                    taxonomy = colnames(bac_otu)[i])
  result = rbind(result,temp )
}

result = result %>% mutate(padj = p.adjust(P, method = "BH"))

result = result %>% arrange(cor_r)
result$taxonomy = sub("_SGB.*","",result$taxonomy)
result[!duplicated(result$taxonomy),] ->result

result$taxonomy = factor(result$taxonomy,levels = result$taxonomy)

write.csv(result,"cor_nod2_table.csv")

result %>% 
  mutate(cclass = case_when(
    padj > 0.05& P <= 0.05 ~ "P < 0.05",
    padj <= 0.05 & P <= 0.05 ~ "FDR < 0.05",
    P > 0.05 ~ "NO"
  )) -> result1

ggplot(result1 %>% filter(cclass != "NO"))+
  geom_col(aes(y = taxonomy,x = cor_r,fill = cclass))+
  scale_fill_manual(values = c("P < 0.05" = "#fbc7cc",
                               "FDR < 0.05" = "#D95050"))+
  scale_x_continuous(breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3))+
  theme_prism(base_size = 10,base_fontface = "plain",
              base_line_size = 0.5) +
  labs(y = "",
       x = "Correlation Between NOD2 Activation")
ggsave("cor_nod2.pdf",width = 6,height = 2.5)