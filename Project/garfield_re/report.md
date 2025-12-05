## Abstract
Variant calling in genomic sequencing generates false positives that require 
filtering. GARFIELD (Genomic vARiant FIltering by dEep Learning moDels) uses 
deep learning to try to improve variant filtering accuracy. We successfully replicated GARFIELD-NGS results using their provided test sets, achieving AUCs within 0.3% the published values. To test the authours claim that GARFIELD was more robust at lower coverage(30x), since this data was not avialabel in the repository, ..



## Introduction


## Methods
Our study consists of two parts: (1) replication of published GARFIELD-NGS 
results using the original datasets provided in their GitHub repository, and 
(2) evaluation of model performance on independently generated variant calls 
at different coverage depths. We implemented all models using the H2O deep 
learning framework in Python.

# Replication using published dataset
We obtained the original training, validation, and test datasets from the 
GARFIELD-NGS GitHub repository (https://github.com/gedoardo83/GARFIELD-NGS). 
The datasets include:
- ILM_SNP_Pre_training.balance.txt
- ILM_SNP_Training.balance.txt  
- ILM_SNP_Validation1.txt
- ILM_SNP_Test.txt
- ILM_INDEL_Pre_training.txt
- ILM_INDEL_Training.txt
- ILM_INDEL_Validation1.txt
- ILM_INDEL_Test.txt

These datasets were derived from NA12878 Illumina exome sequencing data and 
include variants annotated with GATK quality metrics. Each variant is labeled 
as true or false based on comparison with the Genome in a Bottle (GIAB) v3.3.2 
gold-standard variant calls.


# Features

The original paper uses different feature sets for Illumina and Ion Torrent platforms due to their distinct sequencing technologies and variant calling pipelines.

**Illumina Features**
Used for both SNVs and INDELs from GATK variant calls:
1. BaseQRankSum - Rank sum test for base quality
2. ReadPosRankSum - Rank sum test for read position bias  
3. DP - Total read depth
4. FS - Fisher's exact test for strand bias
5. MQ - Root mean square mapping quality
6. MQRankSum - Rank sum test for mapping quality
7. QD - Quality by depth
8. SOR - Symmetric odds ratio for strand bias
9. QUAL - Phred-scaled quality score
10. GQ - Genotype quality for NA12878

**Ion Torrent SNV Features:**
Used for Ion SNVs from TVC (Torrent Variant Caller):
FDP, FAO, QD, FSAF, FSAR, FXX, LEN, HRUN, RBI, VARB, STB, STBP, PB, PBP, 
MLLD, SSSB, QUAL, GQ

**Ion Torrent INDEL Features :**
Used for Ion INDELs (slightly different from SNVs):
FDP, FAO, QD, FSAF, FSAR, FXX, LEN, HRUN, RBI, VARB, STB, STBP, MLLD, SSEN, 
SSEP, SSSB, QUAL, GQ






# GARFIELD-NGS Replication vs Paper Results

| Dataset        | Platform  | Variant Type | Paper Test Set AUC | Our Test Set AUC | Difference |
|----------------|----------|--------------|------------------|----------------|------------|
| Illumina SNV   | Illumina | SNV          | 0.7998           | 0.8021         | +0.0023    |
| Illumina INDEL | Illumina | INDEL        | 0.9269           | 0.9299         | +0.0030    |
| Ion SNV        | Ion      | SNV          | 0.9757           | 0.9596         | -0.0161    |
| Ion INDEL      | Ion      | INDEL        | 0.9464           | 0.9242         | -0.0222    |

