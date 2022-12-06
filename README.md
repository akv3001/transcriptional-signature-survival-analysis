# Transcriptomic Signature based Survival Analysis
Associating transcriptional signature with Survival Data

Question: Associate an identified bone metastasis signature with 

- Clinical cohort used was microarray data with survival stats from 2 independent datasets with shared phenotype of interest. COMBAT was used batch correct expression matrices before combining the two cohorts for downstream analysis
-  Subset gene expression data using signature of interest
-  Quantify signature score from normalized expression matric - apply signature scoring method of choice - using GSVA single sample GSEA or ssGSEA here  https://www.genepattern.org/modules/docs/ssGSEAProjection/4#gsc.tab=0  (  you can also use rank based normalization methods )
-  Stratify patients by high/low quartile signature scores based classification to prepate for surival comparison
-  Apply coxph regression to associate stratified patients with overall survival analysis
-  Multiple hypothesis correct distribution of signature associations


Broader application:
- signature scoring method can be applied to any expression matrix for downstream survival analysis

Associated publication:
*under review*


Other publication with this method documented and applied to:
https://elifesciences.org/articles/20728v1
