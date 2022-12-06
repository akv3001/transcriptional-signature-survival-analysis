# Transcriptomic Signature based Survival Analysis
Associating transcriptional signature with Survival Data

-  Subset gene expression data using signature of interest
-  Quantify signature score from normalized expression matric - apply signature scoring method of choice - using GSVA single sample GSEA or ssGSEA here  https://www.genepattern.org/modules/docs/ssGSEAProjection/4#gsc.tab=0  (  you can also use rank based normalization methods )
-  Stratify patients by high/low quartile signature scores based classification to prepate for surival comparison
-  Apply coxph regression to associate stratified patients with overall survival analysis
-  Multiple hypothesis correct distribution of signature associations
- 
