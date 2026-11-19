setwd("E:\\4 2\\AKP\\Preeclampsia and IUGR\\survival analysis\\Kidney Renal Clear Cell Carcinoma (TCGA, Firehose Legacy)")
#################################################################
##############################working with clinical data#########
#################################################################
dato<-read.csv("kirc_tcga_clinical_data.csv",header = TRUE,stringsAsFactors = FALSE)

library(dplyr)
#datos<-select(dato,Sample.ID,Race.Category,Diagnosis.Age,Tumor.Disease.Anatomic.Site,Neoplasm.Histologic.Grade,Patient.ID,Primary.Tumor.Site,Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Group.Stage,Censor.Status,Overall.survival.in.days)
datos<-select(dato,Sample.ID,Race.Category,Diagnosis.Age,Tumor.Tissue.Site,Patient.ID,Patient.Primary.Tumor.Site,Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code,Censor.Status,Overall.survival.in.days)
#colnames(dato)[c(82,83)]



#ranme the column
#datos<-dplyr::rename(datos,Patient_ID=Patient.ID,anatomic_site=Tumor.Disease.Anatomic.Site,histologic_grade=Neoplasm.Histologic.Grade,Rectime=Overall.survival.in.days,age=Diagnosis.Age,tumour_site=Primary.Tumor.Site,cancer_stage=Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Group.Stage,race=Race.Category)
datos<-dplyr::rename(datos,Patient_ID=Patient.ID,anatomic_site=Tumor.Tissue.Site,Rectime=Overall.survival.in.days,age=Diagnosis.Age,tumour_site=Patient.Primary.Tumor.Site,cancer_stage=Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code,race=Race.Category)
#replace the "" in race category
datos$race=as.character(lapply(datos$race,function(x){gsub("^$","Others",x)}))
#replace the "" in cancer_site with NA
datos$tumour_site=as.character(lapply(datos$tumour_site,function(x){gsub("^$",NA,x)}))
#replace "" in cancer_stage with NA
datos$cancer_stage=as.character(lapply(datos$cancer_stage,function(x){gsub("^$",NA,x)}))
datos$anatomic_site=as.character(lapply(datos$anatomic_site,function(x){gsub("^$",NA,x)}))
#datos$histologic_grade=as.character(lapply(datos$histologic_grade,function(x){gsub("^$",NA,x)}))
#datos$cancer_stage=as.character(lapply(datos$,function(x){gsub("^$",NA,x)}))

library(plyr)
count(datos,'cancer_stage')
count(datos,'race')
count(datos,'tumour_site')
count(datos,'anatomic_site')
#count(datos,'histologic_grade')


#################################################################
##############################working with mRNA data#########
#################################################################



#read rna expression data
expr_mr<-read.csv("data_mrna_agilent_microarray_zscores_ref_all_results.csv",header=TRUE,stringsAsFactors = FALSE)

#function for labelling each expression value
altered_test<-function(x){
  
  if(typeof(x)=="character"){
    d=x
  }
  else if (abs(x)>=1.5){ 
      d="Altered" 
    
      }
    else{
      d="Normal"
      
    }
  d
}


#altered_test(10)


#apply the function over all the colulmn to convert altered unaltered

applyfunc<-function(df,f){
  ds<-matrix(0,nrow = nrow(df),ncol=ncol(df))
  colnames(ds)<-colnames(df)
  for (i in seq(1:ncol(df))){
    ds[,i]<-(sapply(expr_mr[,i],f))
  }
  ds<-as.data.frame(ds)
}
gene_status<-applyfunc(expr_mr,altered_test)

#remove the 01 from patient iD

remove_01<-function(x){
  x<-unlist(strsplit(x,split=""))
  x<-paste(x[0:(length(x)-3)],collapse = "")
  x
}
#remove_01("TCGA-04-1332-01")
gene_status$Patient_ID<-as.character(gene_status$Patient_ID)
gene_status$Patient_ID=unlist(lapply(gene_status$Patient_ID,remove_01))
#View(gene_status)
####################################################################
#########################Merge the tables###########################
###################################################################
#gene_status$Patient_ID=as.character(gene_status$Patient_ID)
combined<-datos%>%inner_join(gene_status)
#View(combined)  
#View(datos)
#View(gene_status)

#relevel the genes as normal as reference factor
applyrevel<-function(combined){
  
  col_names<-colnames(combined)[10:ncol(combined)]
  for(i in col_names){
    #combined$i<-as.factor(combined$i)
    combined[,i]<-as.factor(combined[,i])
    combined[,i]<-relevel(combined[,i],ref="Normal")
    #combined$i<-relevel(fctocombined$i,ref="Normal")
  }
  combined
}
combined<-applyrevel(combined)
######################################################
################ Univariate analysis ################
#####################################################

# Run this right before your Univariate analysis to clean the dataset
variance_check <- sapply(combined[, 10:ncol(combined)], function(x) length(unique(x)) > 1)
genes_to_keep <- names(variance_check)[variance_check]
combined <- combined[, c(colnames(combined)[1:9], genes_to_keep)]

library(survival)
kmsurvo<-Surv(combined$Rectime,combined$Censor.Status)
applycox <- function(combined) {
  models <- list()
  col_names <- colnames(combined)[10:ncol(combined)]
  
  for(i in col_names) {
    # Check if the gene has 2 or more unique categories
    if (length(unique(combined[, i])) < 2) {
      next 
    }
    
    form <- as.formula(paste("kmsurvo ~", i))
    
    # CRITICAL FIX: Run the model and extract summary metrics immediately
    fit <- tryCatch({
      m <- coxph(form, data = combined)
      tss <- summary(m)
      # Extract coef, exp(coef), and p-value safely
      c(tss$coefficients[c(1, 2, 5)]) 
    }, error = function(e) { 
      NULL 
    })
    
    # Store only the numeric results vector, not the whole model object
    if (!is.null(fit)) {
      models[[i]] = fit
    }
  }
  
  final_mode <- as.data.frame(models) 
  final_model = t(final_mode)
  colnames(final_model) <- c("coef", "exp.coef", "p")
  
  as.data.frame(final_model)
}


fs<-applycox(combined)

#class(fs)
fs<-fs%>%mutate(gene=rownames(.))
fs

#write.csv(fs, file = "E:\\4 2\\AKP\\Preeclampsia and IUGR\\survival analysis\\Kidney Renal Clear Cell Carcinoma (TCGA, Firehose Legacy)\\univ2.csv")
######################################################
################ Multivariate Analysis################
#####################################################

fitt<-coxph(kmsurvo~.,data=combined[,10:ncol(combined)])
fitt

#write.csv(fitt, file = "E:\\4 2\\AKP\\Preeclampsia and IUGR\\survival analysis\\Kidney Renal Clear Cell Carcinoma (TCGA, Firehose Legacy)\\multv.csv")

######################################################

