<!-- # Documentation for downloading data 
# Notice that the data files are not pushed to the github  -->
```
mkdir -p data
mkdir -p data/reference
mkdir -p data/sumstats

cd data/sumstats
```

#--------------- Download the baldness sumstats data --------------------# 

```
wget http://cnsgenomics.com/data/mpb/mpb_bolt_lmm_aut_x.tab.zip 
unzip mpb_bolt_lmm_aut_x.tab.zip 
mv mpb_bolt_lmm_aut_x.tab baldness.tab

# Metadata avaiable at: https://cnsgenomics.com/data/mpb/MPB_GWAS_summary_statistics_README.pdf
# link: https://www.ebi.ac.uk/gwas/publications/30573740
# less gwas-association-downloaded_2025-04-25-pubmedId_30573740.tsv
# mv gwas-association-downloaded_2025-04-25-pubmedId_30573740.tsv baldness_yap_et_al.tsv
```
This is saved as baldness.tab

#--------------- Five other datasets to figure out genetic correlations --------------------#

Initially I downloaded from GWAS catlog but it doesn't have allele0 and allele1, so I seached for other sources to get the data. 

#--------------- Download the income sumstats data --------------------#
# It doesn't support wget, so I manually download it to the folder
# link: https://www.ebi.ac.uk/gwas/publications/31844048
# less gwas-association-downloaded_2025-04-25-pubmedId_31844048.tsv
# mv gwas-association-downloaded_2025-04-25-pubmedId_31844048.tsv income_hill_at_al.tsv

The dataset is mannually downloaded from https://datashare.ed.ac.uk/handle/10283/3441
I chose the file Household_Income_UKBiobank.txt, because the effective sample size corresponds with what I expected. 

The dataset is saved as income.txt 

#--------------- Download the Chronotype sumstats data --------------------#
# mv gwas-association-downloaded_2025-04-25-pubmedId_30696823.tsv chronotype_Jones_et_al.tsv

wget https://pmc.ncbi.nlm.nih.gov/articles/instance/6351539/bin/41467_2018_8259_MOESM13_ESM.xlsx 
# do some manual manipulation to get the data
# data is saved as chronotype.csv  
less chronotype.csv  

#--------------- Computer leisure time --------------------#
Link: https://data.mendeley.com/datasets/mxjj6czsrd/1/files/b6fb8a18-7ffc-471a-8e11-0926c4401aeb 
This is saved as computeruse.txt 
# Overall data link for this research: https://data.mendeley.com/datasets/mxjj6czsrd/1

#----------------- AHDH --------------------# 
Link: https://figshare.com/articles/dataset/adhd2022/22564390
saved as adhd.txt 

#--------------- Cannabis dependence ---------# 
Link: https://figshare.com/articles/dataset/sud2020-cud/14842692 
I used CUD_EUR_casecontrol_public_11.14.2020 because it contains beta, which has N cases = 14,080; N controls = 343,726.

#--------------- number of children ever born and age at the first birth male ----------# 
Link: https://www.thessgac.com/papers/7
number of children ever born: sample size = 318463
first birth age for male = 103736


#--------------- Dataset Information ---------------#

# Mannually check the header for each file and then decide 

For baldness,  
SNP: SNP
CHR: chromosome
BP: base pair position
ALLELE1: effect allele
ALLELE0: other allele
BETA: Effect size
P_BOLT_LMM: Pvalue

chronotype data is downloaded for the top 10,000 SNPs only 
For chronotype, 
Locus: SNP
CHR: chromosome
POS: base pair position
Allele1: effect allele
Allele2: other allele
"Meta ln OR": Effect size
Meta P (Corrected): Pvalue

For income, 
Chr: chromosome
SNP: SNP
BPos: basepair position
Non_effect_Allele: The non-effect allele
Effect_Allele: The effect allele
Beta: The regression weight
Standard_Error_of_Beta: The standard error of the Beta
P: The P value of the association test.

For computer use, 
SNP: SNP
CHR: chromosome
BP: base pair position
ALLELE1: effect allele
ALLELE0: other allele
BETA: Effect size
P_BOLT_LMM: Pvalue

For ADHD,
SNP: SNP
CHR: chromosome
BP: base pair position
A1: effect allele
A2: other allele
OR: effect size
P: Pvalue

For cannabis use,
SNP: SNP
CHR: chromosome
BP: base pair position
A1: effect allele
A2: other allele
BETA: effect size
P: Pvalue








#--------------- Download the insomonia sumstats data --------------------#
# DOWNLOAD manually from https://www.ebi.ac.uk/gwas/publications/35835914 
# less gwas-association-downloaded_2025-04-25-pubmedId_35835914.tsv
# mv gwas-association-downloaded_2025-04-25-pubmedId_35835914.tsv insomnia_Watanabe_et_al.tsv 

I only downloaded male data, which is saved as insomonia.txt 


#--------------- Download the Age at first birth sumstats data --------------------#
# link: https://www.ebi.ac.uk/gwas/publications/34211149
less gwas-association-downloaded_2025-04-25-pubmedId_34211149.tsv 
mv gwas-association-downloaded_2025-04-25-pubmedId_34211149.tsv agefirstbirth_Mills_et_al.tsv

#--------------- Download the num of children sumstats data --------------------#
mv gwas-association-downloaded_2025-04-25-pubmedId_27798627.tsv numchildren_Barben_et_al.tsv

#--------------- Download the 1000 Genomes data --------------------# 


