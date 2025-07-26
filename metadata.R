setwd('/hwdata/home/zhangdong/gutMEGA/PRJNA598080/index1/')
library(pacman)
library(stringr)
library(dplyr)
pacman::p_load(tidyverse,magrittr,stringr)
metadata <- "./rawdata/metadata_raw.txt" %>%
  read.delim(check.names = FALSE,header = T,sep = ',')

metadata <- metadata[,c('BioProject','Run','BioSample','Experiment','Sample Name')]
colnames(metadata)[5] <- c('Sample Name')
#colnames(metadata)[6] <- c('Sex')
metadata <- metadata%>%
  mutate(Group=case_when(
    str_detect(metadata$`Sample Name`,pattern="rd")~"GERD",
    str_detect(metadata$`Sample Name`,pattern="hc")~"Health"
  ))
table(metadata$Group)
metadata <- na.omit(metadata)
##添加另外的信息
metadata$`Human/Mice` <- rep(c('Human'),nrow(metadata))
metadata$`Sample location` <- rep(c('Oral'),nrow(metadata))
metadata$`Sample type` <- rep(c('saliva'),nrow(metadata))
metadata$Sex <- NA
#colnames(metadata)[5] <- 'Sex'
metadata <- metadata[,c('Run','BioProject','BioSample','Experiment','Group','Human/Mice','Sample location','Sample type','Sex')]
#write.table(metadata,file = './index17/metadata.txt',col.names = T,row.names = F,quote = F,sep = '\t')
colnames(metadata)[1] <- 'SampleID'
write.table(metadata,file = './result/metadata.txt',col.names = T,row.names = F,quote = F,sep = '\t')

