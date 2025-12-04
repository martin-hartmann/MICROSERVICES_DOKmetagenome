MICROSERVICES_DOKmetagenome

R-Scripts and corresponding data files for the statistical analysis of the rhizosphere metagenome data of the drought sheltering experiment in the DOK trial, which are the base of the manuscript "How agricultural legacy shapes the drought response of the wheat rhizobiome". 

It includes input data files from universal ontologies (e.g. SEED, EC, InterPro2GO, eggNOG) and nutrient-specific ontologies (e.g. CAZy, NCycDB, PCycDB) based on raw counts (*count.tsv) and normalized counts (*RPKmaxORF.tsv) by open reading frame length.

R scripts:
1. Script for analysis of data quality based on raw counts (R-Script_Metagenome-Data-Quality.r)
2. Script for analysis of beta-diversity of universal ontologies based on rarefaction (R-Script_Metagenome-Beta-Diversity-Universal-Ontologies.r
3. Script for analysis of beta diversity from nutrient-specific ontologies (R-Script_Metagenome-Beta-Diversity-Nutrient-Ontology.r)
4. Script for univariate analysis of gene functions based on SEED (R-Script_Metagenome-Gene-Level-Analysis_SEED.r)
5. Script for indicator analysis of gene functions based on SEED (R-Script_Metagenome-Indicator-analysis.r)
