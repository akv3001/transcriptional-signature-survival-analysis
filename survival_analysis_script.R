################################################
# Analysis to correlate Bone metastasis signature 
# with bone met status
################################################


## Data used ##

## Survival data used from published studies  ##
# https://www.cell.com/supplemental/S1535-6108(09)00179-2#secd11908239e374
# GEO IDs :
# GEO2603: Bone metastasis  <--- Relevant for bone metastasis in breast cancer
# GSE5327 : lung metastasis, 
# GSE2034 : Bone mets <--- Relevant for bone metastasis in breast cancer
# GSE12276 : primary breast cancer

# Install GEO Query
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")

library(GEOquery)
library(openxlsx)
library(dplyr)
library(sva) # BATCH CORRECTION using ComBat
require(survminer)
library(survival)
library('GSVA')

## Input gene signature 
gene_signature = read.xlsx('Gene Signatures.xlsx',check.names = FALSE)

#### GLP96 Microarray Reference ####

GLP96_GeneSymbolMap = read.csv('GPL96-57554.txt',skip = 16,sep="\t",strip.white = TRUE,fill=TRUE)

######################## Cohort 1 ###################
# Extract data from the GEO Query
Study_1_GEO2603 <- getGEO('GSE2603',GSEMatrix=TRUE)
show(Study_1_GEO2603)

exprs_data_study1 <- as.data.frame(Study_1_GEO2603$GSE2603_series_matrix.txt.gz@assayData$exprs)
exprs_data_study1 = merge(GLP96_GeneSymbolMap[,c("ID","Gene.Symbol")],exprs_data_study1,by.x='ID',by.y='row.names')
exprs_data_study1.summ = ddply(exprs_data_study1,'Gene.Symbol',numcolwise(mean))
rownames(exprs_data_study1.summ) = exprs_data_study1.summ$Gene.Symbol


# Extract the phenotype data
Study_1_GEO2603_metadata = pData(phenoData(Study_1_GEO2603[[1]]))

# Subset Clinical data from cell lines all the clinical data has T-
Study_1_GEO2603_metadata.clinical = Study_1_GEO2603_metadata[which( grepl(pattern = "-T",Study_1_GEO2603_metadata$title)),]


exprs_data_study1 <- as.data.frame(Study_1_GEO2603$GSE2603_series_matrix.txt.gz@assayData$exprs)
exprs_data_study1 = merge(GLP96_GeneSymbolMap[,c("ID","Gene.Symbol")],exprs_data_study1,by.x='ID',by.y='row.names')
exprs_data_study1.summ = ddply(exprs_data_study1,'Gene.Symbol',numcolwise(mean))
rownames(exprs_data_study1.summ) = exprs_data_study1.summ$Gene.Symbol

exprs_data.summ.fil = exprs_data.summ[,Study_1_GEO2603_metadata.clinical$geo_accession]

######################## Cohort 2 ###################
Study_2_GSE2034 = getGEO('GSE2034',GSEMatrix=TRUE)
show(Study_2_GSE2034)

# Extract the phenotype data
Study_2_GSE2034_metadata = pData(phenoData(Study_2_GSE2034[[1]]))
Study_2_GSE2034_metadata.clinical = read.xlsx('Study_2_GSE2034_survival_data_mets.xlsx')
Study_2_GSE2034_metadata = merge(Study_2_GSE2034_metadata,Study_2_GSE2034_metadata.clinical
                                 ,by.x='geo_accession',by.y='GEO.asscession.number')
# Subset Clinical data from cell lines all the clinical data has T-

exprs_data_study2 <- as.data.frame(Study_2_GSE2034$GSE2034_series_matrix.txt.gz@assayData$exprs)
exprs_data_study2 = merge(GLP96_GeneSymbolMap[,c("ID","Gene.Symbol")],exprs_data_study2,by.x='ID',by.y='row.names')
exprs_data_study2.summ = ddply(exprs_data_study2,'Gene.Symbol',numcolwise(mean))
rownames(exprs_data_study2.summ) = exprs_data_study2.summ$Gene.Symbol



### Combine cohorts ###



## combine clinical data
#subset_clinical_cohort1 = Study_1_GEO2603_metadata.clinical[,c('geo_accession',"bm event:ch1","bm event:ch1", "path er status:ch1","path pr status:ch1"),]
subset_clinical_cohort1 = Study_1_GEO2603_metadata.clinical[,c('geo_accession',"bm event:ch1","bmfs (yr):ch1","path er status:ch1" ),]
subset_clinical_cohort1$bmfs_months = as.numeric(subset_clinical_cohort1$`bmfs (yr):ch1`)*12
subset_clinical_cohort1 = subset_clinical_cohort1[,c("geo_accession","bmfs_months","bm event:ch1","path er status:ch1")]
colnames(subset_clinical_cohort1)  = c("geo_accession","bonemetsfree_months","bonemets_event","ER_Status")
subset_clinical_cohort1$StudyID = "GEO2603"

subset_clinical_cohort2 = Study_2_GSE2034_metadata[,c('geo_accession',"time.to.relapse.or.last.follow-up.(months)","bone relapses (1=yes, 0=no):ch1","ER.Status" ),]
colnames(subset_clinical_cohort2)  = c("geo_accession","bonemetsfree_months","bonemets_event","ER_Status")
subset_clinical_cohort2$StudyID = "GSE2034"


##### Rbind both datasets clinical metadata = 385 Patients
combined_final_clinical_metadata = rbind(subset_clinical_cohort1,subset_clinical_cohort2)


##### combine Expression
Combined_final_Expression = merge(exprs_data_study2.summ[,-1],exprs_data_study1.summ[,-1],by = 'row.names',all=TRUE)
rownames(Combined_final_Expression) = Combined_final_Expression$Row.names
Combined_final_Expression = Combined_final_Expression[,-1]


write.csv(Combined_final_Expression,'Combined_final_Expression_microarray.csv')
write.csv(combined_final_clinical_metadata,'combined_final_clinical_metadata.csv')



Combined_final_Expression.bc = sva::ComBat(Combined_final_Expression[,combined_final_clinical_metadata$geo_accession],batch = combined_final_clinical_metadata$StudyID,mean.only = TRUE)

### Score signature #####


# Input expression table is a microarray dataset - some gene probes aren't captures 
# this loop prints how many genes from the signature are lost in the analysis
# prioritizing signatures based on these values
gene_signature.list = list()
gene_signature.list_used = list() # this is the list of genes that overlaps and was used

lapply(gene_signature.list, length)
for( each_col in 1:ncol(gene_signature)){
  
  sig_name = colnames(gene_signature)[each_col]
  gene_signature.list[[sig_name]] = as.matrix(na.omit(gene_signature[,each_col]))
  
  overlap_input = intersect(rownames(Combined_final_Expression),gene_signature.list[[sig_name]])
  gene_signature.list_used[[sig_name]] = overlap_input
  print(paste0(sig_name," Total_genes=",length(gene_signature.list[[sig_name]])," Overlapping genes=",length(overlap_input)))
  
}

# check length of signature
lapply(gene_signature.list, length)
lapply(gene_signature.list_used, length)





# Use all signatures for up genes 
set.seed(50)
GSVA_signature_score = gsva(expr = as.matrix(Combined_final_Expression.bc),gset.idx.list = gene_signature.list_used)
GSVA_signature_score.up = GSVA_signature_score[grepl('Up',rownames(GSVA_signature_score)),]

GSVA_signature_score.up = t(GSVA_signature_score.up)
#merge with clinical data
GSVA_signature_score.up_clinical = merge(GSVA_signature_score.up,
                                         combined_final_clinical_metadata,
                                         by.x='row.names',by.y='geo_accession')



###### Survival analysis ########



pdf('Survival_Curves_BMF_Survival_batchcorrected.pdf',onefile = TRUE)
for (gene_signature in colnames(GSVA_signature_score.up)){
  


input_survival = GSVA_signature_score.up_clinical
subset_by = quantile(input_survival[[gene_signature]])

input_survival$GeneSignatureClass = ifelse(input_survival[gene_signature] <= subset_by[[2]],"Low", ifelse(input_survival[gene_signature] >= subset_by[[3]],"High",'-' ))
input_survival = input_survival[which(input_survival$GeneSignatureClass != "-"),]
input_survival$GeneSignatureClass = as.factor(input_survival$GeneSignatureClass)
input_survival$bonemets_event = as.numeric(input_survival$bonemets_event)
fit <- survfit(Surv(bonemetsfree_months, bonemets_event) ~ GeneSignatureClass, data = input_survival)

P=ggsurvplot(fit, data = input_survival,
             pval = TRUE,risk.table = TRUE,
             palette = 'Dark2') + # Customize colors here
  ggtitle(paste0("Signature= ",gene_signature),subtitle = paste0("Bone Free Metastasis Data: GEO2603 & GSE2034"))
print(P)

P_ER = ggboxplot(input_survival,x = )

}

dev.off()










