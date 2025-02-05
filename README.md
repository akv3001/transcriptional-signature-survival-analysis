# Transcriptomic Signature based Survival Analysis
Associating transcriptional signature with Survival Data. This analysis is insightful in associating model derived findings to highlight a novel gene signature's translation value by associating it to patient data

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
*published work that applied this analysis*


Other publishes work with this method also documented and applied to:

- Oncogenic BRAF disrupts thyroid morphogenesis and function via Twist expression: https://elifesciences.org/articles/20728v1 

- Histone 3 Methyltransferases Alter Melanoma Initiation and Progression Through Discrete Mechanisms:  https://pubmed.ncbi.nlm.nih.gov/35223844/

- Bone-matrix mineralization dampens integrin-mediated mechanosignalling and metastatic progression in breast cancer : https://www.nature.com/articles/s41551-023-01077-3
