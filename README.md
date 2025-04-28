# The project re-analyzes the data from Yap et al. 2018 
Citation: Yap CX, Sidorenko J, Wu Y, Kemper KE, Yang J, Wray NR, Robinson MR, Visscher PM. Dissection of genetic variation and evidence for pleiotropy in male pattern baldness. Nat Commun. 2018 Dec 20;9(1):5407. doi: 10.1038/s41467-018-07862-y. PMID: 30573740; PMCID: PMC6302097.y

# Objective: 
1) Estimate heritability --> done
2) Finish Manhattan plot + QQplot --> Done 
3) Have cell-type specific heritability --> Done
4) Compare with other traits 
5) 
6) 
7) Five traits to look at: 

    a. Insomnia: https://www.nature.com/articles/s41588-022-01124-w#Abs1 We downloaded the male dataset only 

    b. Leisure computer use: https://pmc.ncbi.nlm.nih.gov/articles/PMC7174427/

    Number of Children ever born: https://www.ncbi.nlm.nih.gov/pubmed/27798627

    c. ADHD: https://pmc.ncbi.nlm.nih.gov/articles/PMC10914347/#ABS1 
    

    Age at first birth: https://pubmed.ncbi.nlm.nih.gov/34211149/

    d. Chronotype: 

    e. Income: https://pubmed.ncbi.nlm.nih.gov/31844048/ 

# Structure 
.
├── data --> me finish this later 
|
|
├── results/                      
│   ├── plots/                   # Visualization outputs
│   │   └── enrichment_plot.png
│   ├── step1.log               # Need to re-name this 
│   ├── step1.sumstats.gz        
│   ├── step2.log
│   ├── step3.log
│   └── step3.results           
│
├── scripts/                     
│   ├── CellType.R               
│   ├── heritability.sh          
│   ├── heritibility.py          
│   ├── LDSC_EnrichmentPlot.R    
│   ├── Manhattan_qq.R          
│   ├── p_value_calculation_from_h2.py 
│
├── main.py or main.sh --> let's do this after finishing it 
├── celltype_enrichment_table.docx  
├── documentation.md               
├── get-pip.py                    
└── README.md                   


# Initializing... 
## Download GWAS summary statistics 

## Make sure ldsc is installed 
Using ldsc 2.0.1 from here: https://github.com/svdorn/ldsc-2.0.1 
Python version 3.9.

## Make sure 1000 Genomes Phase 3 reference is downloaded 
HapMap3 SNP list is used as the reference snp list. 
European 

# Re-analyzing Yap et al. 2018: 
Based on Yap et al. 2018, hair loss is a quantatitive trait. 

## Find heritability 
```
chmod 740 ./scripts/heritability.sh 
./scripts/heritability.sh 'path_to_input_file' 'path_to_ldsc_executable' 'num_sample' 'path_to_reference' 'output_directory'
```

Heritability results for balsness using p-value based on LLM-INF: 
***
Read regression weight LD Scores for 1293150 SNPs.
After merging with reference panel LD, 1139574 SNPs remain.
After merging with regression SNP LD, 1139574 SNPs remain.
Using two-step estimator with cutoff at 30.
Total Observed scale h2: 0.3212 (0.0256)
Lambda GC: 1.4926
Mean Chi^2: 2.5062
Intercept: 1.2273 (0.0198)
Ratio: 0.1509 (0.0131) 
*** 

Heritability results for balsness using p-value based on LLM: 
***
Total Observed scale h2: 0.3431 (0.0274)
Lambda GC: 1.4926
Mean Chi^2: 2.5854
Intercept: 1.2232 (0.0205)
Ratio: 0.1408 (0.0129)
Analysis finished at Sun Apr 27 12:02:33 2025
Total time elapsed: 10.51s
*** 


***
z-score: 12.546874999999998
p-value: 0.00e+00
*** 

PS: I will change the whole pipeline with one main.py and src/utility later! It's a bit messey now since it's simply preliminary results 

# Manhattan plot and QQplot


