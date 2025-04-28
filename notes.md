# Notes for Yap et al. 2018 
# Methods: 
1) MPB score is highly correlated with age. Each year age increase confers a 0.02 increase in MPB score. 
2) MPB scores were regressed on age, assessment center, self‑reported ethnicity and 40 genetic principal components (UKB European PCs), and the residuals were used as the “adjusted MPB score” (SD=1.1) 
3) After standard QC (imputation info>0.3, MAF>1e-4 for GWAS; MAF>0.01 for heritability), ~18 million autosomal SNPs and ~1.06 million X‑chr SNPs remained for association testing​. 
4) For SNP‑heritability estimation, 1.12M HapMap3 autosomal SNPs (MAF>0.01) and 262,588 X‑chr SNPs (MAF>0.01) were retained. 

# Results: 

## Heritability (Figure 2a gray bars)

1) Pedigree‐based: 

Brothers had h^2 ≈ 0.62 (SE=0.028), whereas father–son pairs had a lower h^2 ≈ 0.41 (SE=0.071); the difference was statistically significant (P=0.006), which could be explained by the effect from maternal side. Thus, brother-brother pairs explain the pedigree-based heritability better. The overall pedigree h^2 (all 1st‑degree pairs) was ≈0.59 (SE=0.026). 

Repeatability measures the consistency of a trait across repeated observations in the same individuals. In this study, MPB scores were measured at two time points in ~9,603 men (2006–2010 and 2012–2013). The Pearson correlation (r = 0.88) between the two measurements means that individuals who had more baldness at the first timepoint tended to have similar scores later on.  Repeatability is an upper bound on heritability (h²) because it includes both genetic effects and permanent environmental effects. 

The pedigree-based h² = 0.62 and SNP-based h² = 0.39–0.41, which are substantial portions of that 0.88.This means the estimated heritability (h²) from pedigree or SNP data is getting close to the repeatability value (0.88). 

2) Heritability based on unrelated men: 

SNP‑heritability (h^2_SNP) was estimated in unrelated men (rel<0.05). Jointly fitting autosomal and X‑chromosome genetic relatedness matrices yielded total h^2_SNP ≈ 0.393. The autosomal component was 0.359 (SE=0.006) and the X‑chromosome contributed 0.034 (SE=0.002). 

Thus ~8.6% of SNP-heritability resides on the X‑chromosome. (Out of total h²_SNP = 0.393: X-chromosome contributes 0.034 (from MGRM estimates). 0.034 / 0.393 ≈ 8.6%). This reflects the strong role of the androgen receptor (AR) locus on the X chromosome, which is a major known contributor to MPB. The X-chromosome is haploid in males (only one copy), but effects are strong and still quantifiable.

Figure 2a colored bars: 

Common variants (MAF, Minor Allele Frequency >0.01) carried virtually all of this signal (h^2_SNP,common = 0.41, SE=0.008. (PS: MAF is the frequency of the less common allele at a given SNP in the population.) Rare variants contributed essentially zero). Likewise, variants in low-LD regions were enriched for heritability. 

Figure 2b (right panel) shows that the 622 independent genome-wide significant SNPs (see Fig. 3) jointly explain h^2 ≈ 0.252 of the trait variance. Of this, 0.223 is from the 598 autosomal SNPs and 0.029 from the 24 X‑chr SNPs. 

PS: MGRM = Multi-GRM (Genomic Relationship Matrix) model. It allows fitting multiple GRMs simultaneously to separate contributions: One GRM from autosomal SNPs. One GRM from X-chromosome SNPs. This helps partition heritability across different parts of the genome (e.g., autosomes vs. X). 

PS: What is big K / small K analysis? This is a variance partitioning method: Big K: GRM using close relatives (higher genetic similarity). Small K: GRM using distant/unrelated individuals (lower similarity). This method simultaneously estimates pedigree-based and SNP-based heritability.

In Figure 2a colored bars: 
bKsK: Three GRMs (big-K, small-K, and X-chromosome) to estimate how much h² is explained by SNPs vs. relatives --> This account for both related and un-related individuals. 
MGRM: SNP heritability using only autosomal and X-chr GRMs (for unrelated individuals).
LDMS: Partitioned SNP heritability by MAF and LD (low vs. high).
LD: If causal variants are in low LD regions, it means they don’t have many nearby SNPs that tag them well. This can reduce power to detect them unless they're directly genotyped or imputed. Here, SNPs in low LD regions explain more heritability than expected (enriched). SNPs in high LD regions explain less heritability, even though they’re more numerous. It suggests that causal variants for MPB tend to be in low-LD regions, where SNPs are more independent and effects aren’t diluted across many correlated SNPs. 

# Genome-wide association and independent loci (Figure 3) 

The GWAS workflow diagram in Figure 3a outlines the steps (quality control, mixed-model association, COJO selection). 

A mixed-model GWAS was performed on all 205,327 men using BOLT-LMM, which accounts for population structure and relatedness. Adjusted MPB score (quantitative) was the outcome. Autosomal and X‑chromosomal SNPs were tested, with male X alleles coded as haploid (BOLT-LMM encodes male genotypes as 0/2, effectively halving effect sizes)

The Manhattan plots in Figure 3b (autosomes) and 3c (X‑chromosome) display the −log10 GWAS p-values for every SNP. Thousands of variants reached genome-wide significance (P<5×10^−8). Notably, chromosome X shows a dense cluster of highly significant SNPs (including the known AR locus). 

Figure 3d is an LD heatmap for the X‑chromosome 22.31–13.2 region containing several jointly significant SNPs from COJO (as indicated in Fig. 3c). Finally, Figure 3e (QQ-plot) shows that the observed test statistics deviate substantially from the null, consistent with a polygenic trait of very many associations

PS: 
1) BOLT-LMM is a tool for performing genome-wide association studies (GWAS) using linear mixed models (LMMs). Adjusted MPB is the outcome after dealing with covariants like age, population, etc, and BOLT-LMM estimates an effect size (beta) for each SNP: how much each SNP is associated with variation in adjusted MPB.

2) GCTA-COJO: Find independent association signals from GWAS results. It uses summary statistics and a reference panel (from genotype data) to account for LD between SNPs. Independent locus = a region of the genome where a SNP shows association not explained by LD with other SNPs. They identified 622 conditionally independent loci: 598 autosomal and 24 X-chromosomal. These are non-redundant signals of association, and help focus on biologically distinct regions involved in MPB.

3) Figure 3d: LD heatmap 
Each cell in the heatmap corresponds to LD between a pair of SNPs. Darker colors = high LD (strong correlation) Lighter = low LD (little correlation). In Figure 3d, they used the heatmap for the X-chromosome region with strong MPB association: -> Shows how several SNPs are highly correlated (in LD).


# Functional annotation and enrichment: 

They next examined the GWAS hits and heritability enrichments. Annotation of the 622 lead SNPs and all genome-wide significant SNPs (via FUMA/SNP2GENE) showed that most lie in noncoding regions (intronic or intergenic) with few coding variants​. 

Partitioned heritability using LDSC found significant enrichment of SNPs in functional genomic regions (Check Supplementary Data).

SNPs in conserved regions and regulatory marks (e.g. enhancers marked by H3K27ac, H3K4me1/3, DHS, TF binding sites) were highly enriched for MPB heritability (P < 10^−6)

Tissue-expression analysis (FUMA GENE2FUNC) highlighted skin tissues (including hair follicle relevant tissues) as enriched among MPB-associated genes. Pathway analysis (ICSNPathway) nominated reproductive hormone pathways. 

These enrichments suggest that many MPB loci affect hair growth, skin, and developmental processes – consistent with the phenotype.

# Pleiotropy and genetic correlations

They tested for genetic overlap between MPB and other traits using LDSC and LDHub. Among female-specific traits, later age at menarche was genetically correlated with less baldness (r_g ≈ −0.09, P<2e-8). Among male traits, earlier puberty indicators (younger relative age of facial hair, voice-breaking) were strongly associated with higher MPB (r_g ≈ −0.18 and −0.11, respectively)​. Phenotypically, men with higher MPB score had slightly fewer children: MPB score4 men fathered ~0.09 fewer children than score1 men, thus a selection gradient of about −0.018 per SD of MPB.

MPB shares genetic factors with puberty timing, bone biology, and metabolism, implying pleiotropic effects and potential evolutionary selection against early-onset baldness.

PS: 
1) Pleiotropy is when one gene or genetic variant affects multiple traits. They used LDSC (LD score regression) to test for genetic correlations between MPB and other traits. A significant genetic correlation implies that shared genetic variants influence both traits, which is classic pleiotropy.

MPB-associated variants were also associated with: a. Age at puberty (facial hair onset, voice breaking, age of menarche); b. Bone mineral density. c. Pancreatic beta-cell function (HOMA-B)

PS on their metadata for summary statistics: 
P_BOLT_LMM_INF (Infinitesimal Model): All SNPs have small effects (very polygenic). A polygenic background, as shown by high SNP-heritability and many small-effect SNPs. 
P_BOLT_LMM (Non-infinitesimal Model): Some SNPs have large effects (less polygenic). Some strong-effect loci (like the AR gene on the X-chromosome), which violate the infinitesimal assumption. 
































