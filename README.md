# OMR
OMR is implemented as an open source R package for two-sampling Mendelian randomization analysis using omnigenic genetic architecture. 

OMR makes an omnigenic modeling assumption on the SNP effects on the exposure variable and uses genome-wide SNPs as instrumental variables without pre-selection. OMR directly uses summary statistics and imposes a more general modeling assumption of horizontal pleiotropy than that used in the Egger regression. OMR takes LD information into account through joint modeling of all SNPs and relies on a scalable composite likelihood framework to infer causal effect. 

# Installation
OMR is implemented as an R package, which can be installed from GitHub.

####  Install from GitHub
```
library(devtools)
install_github("wanglu205/OMR")
```

# Usage
The main functions is OMR. You can find the instructions and an example by '?OMR'.

# Example
```
data(exampledata)
ny <- exampledata$ny
nx <- exampledata$nx
num.per <- exampledata$num.per
l.j <- exampledata$l.j
Z <- exampledata$Z
res=OMR(ny,nx,num.per,l.j,Z)
```

# Example: body mass index (BMI) and asthma
#### Data source
* Exposure: 2,554,637 SNPs from GIANT consortium on BMI European ancestry [Locke et al, 2015, PMID: 25673413].
* Outcome: 2,001,280 SNPs from Trans-National Asthma Genetic Consortium (TAGC) on asthma European ancestry [Demenais et al ,2018, PMID 29273806]. 
* Reference: about 84.7 million SNPs from European ancestry from the 1,000 Genomes Project Phase3 [1000 Genomes Project Consortium, 2015]. 

#### Step 1. Load package and data
The first step is to prepare the data set files. Datasets for exposure and outcome can be obtained using `wget`
```
wget https://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz data/
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/DemenaisF_29273806_GCST005212/TAGC_meta-analyses_results_for_asthma_risk.zip data/
cd data/
unzip TAGC_meta-analyses_results_for_asthma_risk.zip
```
Let’s assume you have downloaded required data set and saved in the folder data/. We will read in the data set using fread. The bim file of Reference data was loaded and the MHC region was removed. The formatting procedure is provided below; you can follow the following steps to prepare data file or chose to use the already formatted data file. The formatted data files are available below.
```
library(data.table)
library(dplyr)
X<-fread("zcat /net/fantasia/home/borang/data2/real_data/hg_19/SNP_gwas_mc_merge_nogc.tbl.uniq.gz")
Y <- fread("/net/mulan/home/yuef/realdata/asthma/raw_data/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv")
nx <- round(mean(X$N))
X <- X %>% select(SNP=SNP, beta.x = b, se.x = se, 
                   A1.x = A1, A2.x = A2) 
Y <- Y %>% select(SNP=rsid, beta.y = European_ancestry_beta_rand, se.y = European_ancestry_se_rand, 
                   A1.y = alternate_allele, A2.y = reference_allele)
ny <- 127669				   
REF<-fread("/net/fantasia/home/borang/data2/real_data/ldscore_1000_genome_less_SNP/eur_chr_all.bim")
REF <- REF %>% select(SNP=V2,CHR=V1,BP=V4,A1.ref = V5, A2.ref = V6) %>% 
       filter(!(CHR == 6 & BP>2*10^7&BP<3*10^7)) 
nref <- 503
```
#### Step 2. Harmonize data
First we will remove strand ambiguous alleles (A/T or C/G alleles) from exposure and outcome datasets. 
```
X <- X %>% filter(!((A1.x=="G"&A2.x=="C")|(A1.x=="C"&A2.x=="G")|(A1.x=="A"&A2.x=="T")|(A1.x=="T"&A2.x=="A")))
Y <- Y %>% filter(!((A1.y=="G"&A2.y=="C")|(A1.y=="C"&A2.y=="G")|(A1.y=="A"&A2.y=="T")|(A1.y=="T"&A2.y=="A")))
```
Then merge exposure and outcome data set with reference panel respectively, and find alleles that do not match for reference data set. Flip the sign of effect size if the effect allele in exposure or outcome is not the minor allele in reference panel.
```
data_x_ref = X %>% inner_join(REF,by="SNP") 
data_y_ref = Y %>% inner_join(REF,by="SNP") 

data_x_ref %>% filter(!((A1.x==A1.ref&A2.x==A2.ref) | (A1.x==A2.ref&A2.x==A1.ref)))%>% nrow()  ##0 - Alleles are the same.
data_x_ref = data_x_ref %>% mutate(beta.x = ifelse(A1.x== A1.ref,beta.x,-beta.x))

data_y_ref %>% filter(!((A1.y==A1.ref&A2.y==A2.ref) | (A1.y==A2.ref&A2.y==A1.ref)))%>% nrow() ##0 - Alleles are the same.
data_y_ref = data_y_ref %>% mutate(beta.y = ifelse(A1.y== A1.ref,beta.y,-beta.y))
```
Merge exposure and outcome data set and keep the following variables: SNP, CHR, BP, summary statistic for exposure and outcome variable.
```
data <- data_x_ref %>% inner_join(data_y_ref,by="SNP") %>%
        mutate(z.x = beta.x/se.x) %>%
		      mutate(z.y = beta.y/se.y) %>%
        select(SNP,CHR=CHR.x,BP=BP.x,z.y,z.x)###610963 SNPs
write.table(data$SNP,file="snp_list.txt",col.names = F,row.names = F,quote = F)
```
#### Step 3. Calculate LD scores
The LDSC (ref: Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.) was employed to calculate LD Scores using reference panel.
```
./plink --bfile eur_chr_all --extract snp_list.txt --make-bed --out eur_for_BMI_asthma
python ldsc.py --bfile eur_for_BMI_asthma --l2 --ld-wind-kb 10000 --out eur_for_BMI_asthma
```
#### Step 4. OMR analysis
With prepared z-scores and LD scores, we could perform OMR using:
```
Z <- cbind(data$z.y,data$z.x)
LDscore<-read.table("eur_for_BMI_asthma.l2.ldscore.gz",header = T)
l.j <- LDscore$L2
res=OMR(ny,nx,nref,l.j,Z)
```
# Results reproduced
All results from all methods used in the OMR paper can be reproduced at 

 <https://github.com/willa0205/OMRreproduce>.

Details in [here](https://github.com/willa0205/OMRreproduce)

## Our group

 <http://www.xzlab.org>.
