xbox::chdir("../5.cluster_trend")
rm(list = ls())
set.seed("2024")

library(phyloseq)
library(igraph)
library(markovchain)
library(tidyverse)

samdat.prune_prev <- function(samdat) {
  GAP_MIN <- 2
  GAP_MAX <- 42
  samdf <- data.frame(samdat)
  subjects <- unique(samdf$SubjectID)
  csub <- split(samdf, samdf$SubjectID)
  for(sub in subjects) {
    cc <- csub[[sub]]
    cc <- cc[order(cc$GDColl),]
    cc$PrevID <- c(NA, cc$SampleID[-nrow(cc)])
    del <- cc$GDColl - c(-999, cc$GDColl[-nrow(cc)])
    keep <- del>=GAP_MIN & del<=GAP_MAX
    if(sum(keep) == 0) {
      csub[[sub]] <- NULL
    } else {
      cc <- cc[keep,]
      csub[[sub]] <- cc
    }
  }
  return(do.call(rbind, csub))
}

all_data = read.csv("../1.profile/final/pcao_meio,bray.data.csv")

all_data %>% 
  dplyr::select(Sample_ID,PatientID,Age_DOL,rep_genus) %>%
  column_to_rownames("Sample_ID") -> samdf

colnames(samdf) = c("SubjectID","GDColl","type")

samdf$SampleID = rownames(samdf)

samdf$type = as.factor(samdf$type)
CSTs <- levels(samdf$type)
nstates <- nlevels(samdf$type)
samdf_prev <- samdat.prune_prev(samdf)

rownames(samdf_prev) <- samdf_prev$SampleID
samdf_prev$PrevCST <- data.frame(samdf)[samdf_prev$PrevID,"type"]
samdf_prev$CurCST <- samdf_prev$type

ttab <- table(samdf_prev$PrevCST, samdf_prev$CurCST) # prevstate=row, curstate=col
trans <- matrix(ttab, nrow=nstates)
trans <- trans/rowSums(trans)  # Normalize row sums to 1
CSTtrans <- trans
colnames(CSTtrans) <- CSTs
rownames(CSTtrans) <- CSTs
t_persist <- -1/log(diag(CSTtrans))

##plot markov chain
mcPreg <- new("markovchain", states=CSTs,
              transitionMatrix = trans, name="PregCST")
netMC <- markovchain:::.getNet(mcPreg, round = TRUE)
wts <- E(netMC)$weight/100
edgel <- get.edgelist(netMC)
edgel_with_weights <- cbind(edgel, wts)

vert.sz <- 2*sapply(states(mcPreg), 
                    function(x) nrow(unique(sample_data(samdf)[sample_data(samdf)$type==x,"SubjectID"])))
vert.sz <- log(vert.sz)*5
colnames(edgel_with_weights) = c("source","target","wts")

as.data.frame(edgel_with_weights) %>% 
  group_by(source) %>% 
  mutate(wt_s = scales::rescale(as.numeric(wts))*100) -> edgel_with_weights

write.csv(edgel_with_weights, "edge_list.csv", row.names = FALSE)
node_info <- data.frame(Node = states(mcPreg), Size = vert.sz)
write.csv(node_info, "node_sizes.csv", row.names = FALSE)

data_num = data.frame(table(all_data$rep_genus)) %>%
  mutate(num_ratio = Freq/sum(Freq))


node_info %>% left_join(data_num,by = c("Node" = "Var1")) %>%
  write.csv("node_sizes.csv", row.names = FALSE)
