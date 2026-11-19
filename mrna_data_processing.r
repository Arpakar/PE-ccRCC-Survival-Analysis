setwd("E:\\4 2\\AKP\\Preeclampsia and IUGR\\survival analysis\\Kidney Renal Clear Cell Carcinoma (TCGA, Firehose Legacy)")
pc<-read.csv("Genes_PC.csv",header = TRUE,stringsAsFactors = FALSE)

mrna<-read.csv("data_mrna_agilent_microarray_zscores_ref_all_samples.csv",header = TRUE,stringsAsFactors = FALSE)

data<-mrna[(mrna$Hugo_Symbol %in% pc$Gene_Symbol),]
data2<-t(data)

write.table(data2, file = "data_mrna_agilent_microarray_zscores_ref_all_results.csv",sep=",",  col.names=FALSE)

