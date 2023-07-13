# MDR1 Variants

This repository contains the code used to run the analyses described in the paper "Revisiting the effects of MDR1 Variants using computational approaches" by Tal Gutman and Tamir Tuller. 

## Main Idea
MDR1 encodes for p-glycoprotein (p-gp), a transmembrane efflux pump known to be related to patient response to chemotherapy. T1236C, T2677G and T3435C are variants in the MDR1 gene that are suspected be impact the function of p-gp.
In this project we perform multiple computational analyses an order to examine the effects of these variants on the different steps of MDR1 expression, hoping to gain insight about the mechanism of action of these variants and their effect on response to chemotherapy and patient survival. 

## The Data 
Several data sources were used in this project. 
* Genomic (SNV), expression and clinical data of TCGA patients were downloaded in November 2021 (https://www.cancer.gov/tcga).
* Protein and CDS sequence of MDR1 was downloaded from Ensembl (https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000085563;r=7:87503017-87713323)
* Human CAI codon weights were calculated using protein abundance measurements from PAXdb (https://pax-db.org/dataset/9606/1502934799/) and the complete human coding sequences (CDS) downloaded from Ensembl (https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cds/). 
* Frequency per 1000 codons values were downloaded from Kazusa (https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606). 
* tAI tissue-specific weights were taken from [1] and the s weights were optimized as depicted in [2]. 
* For the co-translational folding model we downloaded orthologs for all human protein coding genes from Ensembl (https://rest.ensembl.org/documentation/info/homology_ensemblgene)

## The Code
The directory "MDR1_code" contains the code that is needed to run all the analyses depicted in the paper. 

One of the analyses regards co-translational folding, for which we use a model that is developed in our lab and is not published yet. 
The directory "co_trans_code" contains the code for the co-translational folding model, whose output is used in the notebook "Co_translational_folding_MDR1.ipynb" in the "MDR1_code" directory. 

## Requirements 

The notebooks MDR1_Code/Splicing_MDR1 and MDR1_Code/Enformer_MDR1 were ran using different environments. The environments can be created using "spliceAI_requirements.txt" and "Enformer_requirements.txt" respectively. The rest of the notebooks can be ran using the environment created with "MDR1_requirements.txt". 

## References

[1] Hernandez-Alias, X., Benisty, H., Schaefer, M. H. & Serrano, L. Translational adaptation of human viruses to the tissues they infect. Cell Rep 34, 108872 (2021)

[2] Sabi, R. & Tuller, T. Modelling the Efficiency of Codon–tRNA Interactions Based on Codon Usage Bias. DNA Research 21, 511–526 (2014).




```

