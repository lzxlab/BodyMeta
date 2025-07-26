library(phyloseq)
library(tidyverse)
library(microbiome)
library(ggalluvial)
library(ggtree)
library(microbiomeMarker)
library(igraph)
library(sna)
library(ggClusterNet)
library(tidyfst)
library(SpiecEasi)
library(ggprism)
library(tibble)
library(ggpicrust2)
library(SpiecEasi)
library(microbiomeMarker)
library(ggh4x)
#BiocManager::install('ggpicrust2')
#install.packages("ggpicrust2")
setwd('~/Bodymeta/')
set.seed(123)
biom_data <- phyloseq::import_biom(BIOMfilename = "./phyloseq/taxa_table.biom", treefilename = "./phyloseq/tree.nwk")
mapping_file <- phyloseq::import_qiime_sample_data(mapfilename = "./result/metadata.txt")
ps <- merge_phyloseq(biom_data, mapping_file)
colnames(tax_table(ps))= c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")


## 1.Community diversity analysis
group = "Group"
samplesize = min(phyloseq::sample_sums(ps))
if (samplesize == 0) {
  print("0 number sequence of some samples")
  print("median number were used")
  ps11  = phyloseq::rarefy_even_depth(ps,sample.size = samplesize)
} else{
  ps11  = phyloseq::rarefy_even_depth(ps,sample.size = samplesize)
}
ps11 = phyloseq::filter_taxa(ps, function(x) sum(x ) >0 , TRUE)
mapping = phyloseq::sample_data(ps11)
mapping$Group = as.factor(mapping$Group)
metadata <- as.data.frame(mapping)

count = as.data.frame(t(ggClusterNet::vegan_otu(ps11)))
x = t(count) 
# head(x)
Shannon = vegan::diversity(x)  #
Inv_Simpson <- vegan::diversity(x, index = "invsimpson")
S <- vegan::specnumber(x)  
Pielou_evenness <- Shannon/log(S)
Simpson_evenness <- Inv_Simpson/S
est <- vegan::estimateR(x)
Richness <- est[1, ]
Chao1 <- est[2, ]
ACE <- est[4, ]
report = cbind(Shannon, Inv_Simpson, Pielou_evenness, Simpson_evenness,Richness, Chao1,ACE)
report <- as.data.frame(report)
alp = report[rownames(metadata),]
alp$Run <- rownames(metadata)
alp$Group <- metadata$Group
write.csv(alp,'./result/alp.csv',col.names = T,row.names = F)


###1.2beta diversity
group = "Group"
dist = "bray"
Micromet = "adonis"
# Obtain relative abundance
ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
metadata <- as.data.frame(mapping)
#---------PCoA#----
method ="PCoA"
if (method == "PCoA") {
  unif = phyloseq::distance(ps1_rela , method=dist, type="samples")
  pcoa = stats::cmdscale(unif, k=2, eig=T) 
  points_PCoA = as.data.frame(pcoa$points) 
  colnames(points_PCoA) = c("x", "y") 
  points_PCoA0 <- points_PCoA
  eig = pcoa$eig
}
points_PCoA <- points_PCoA[rownames(metadata),]
points_PCoA$Run <- rownames(metadata)
points_PCoA$Group <- metadata$Group
write.csv(points_PCoA,'./result/points_PCoA.csv',col.names = T,row.names = F)


#---------------NMDS#----
method = "NMDS"
if (method == "NMDS") {
  ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
  points_NMDS = ordi$points[,1:2]
  colnames(points_NMDS) = c("x", "y")
  points_NMDS <- as.data.frame(points_NMDS)
  stress = ordi$stress
  stress= paste("stress",":",round(stress,2),sep = "")
}
points_NMDS <- points_NMDS[rownames(metadata),]
points_NMDS$Run <- rownames(metadata)
points_NMDS$Group <- metadata$Group
write.csv(points_NMDS,'./result/points_NMDS.csv',col.names = T,row.names = F)



### 2.Species composition analysis
label = TRUE
sd = FALSE
label = TRUE
sd = FALSE
for (i in c("Phylum","Class","Order","Family","Genus")) {
  psdata <- ggClusterNet::tax_glom_wt(ps = ps,ranks = i)
  psdata = psdata%>%phyloseq::transform_sample_counts(function(x) {x/sum(x)} )
  Taxonomies <- psdata %>% phyloseq::psmelt()
  Taxonomies <- Taxonomies[,c('OTU','Abundance','Sample','Group')]
  file_path <- paste0('./result/taxonomy_',i,'.csv')
  write.csv(Taxonomies,file = file_path,col.names = T,row.names = F)
}

####3.Biomarker identification(LEfSe)
mm_lefse <- run_lefse(ps,group = 'Group',taxa_rank ='Genus',lda_cutoff = 2)
mm_lefse <- as.data.frame(mm_lefse@marker_table)
write.csv(mm_lefse,'./result/lefse.csv',col.names = T,row.names = F)



###4.Functional prediction
metadata <- as.data.frame(ps@sam_data)
path_abundance  <- read.delim('./picrust2/pathways_out/path_abun_unstrat.tsv.gz')
sample_name <- colnames(path_abundance)
sample_name <- sample_name[sample_name%in%rownames(metadata)]
path_abundance <- cbind(path_abundance[,1],path_abundance[,colnames(path_abundance)%in%sample_name])
colnames(path_abundance)[1] <- 'pathway'                                             
S <- nrow(metadata[!duplicated(metadata$Group),])
if(S>2){
  method='ALDEx2_Kruskal-Wallace test'
}else{
  method='ALDEx2_Wilcoxon rank test'
}

##MetacycFeature Notes
metacyc_df <- pathway_daa(abundance = path_abundance %>% column_to_rownames("pathway"),
                          metadata = metadata, group = "Group", daa_method = "ALDEx2")
metacyc_df <- metacyc_df[metacyc_df$method==method,]
metacyc_annotated_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_df, ko_to_kegg = FALSE)
write.csv(metacyc_annotated_df,'./result/function.csv',col.names = T,row.names = F)
