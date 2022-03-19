![图片名称](https://camo.githubusercontent.com/9e54064fb698af20a2b6089b4f16ec3e31f31f72b47f15a5bb215bfd2e41d1b2/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f4c6963656e73652d47504c25323076332d626c75652e737667)
# Identification of Expression QTL by QTLtools in a Rice Recombinant Inbred Line Population
A comprehensive eQTL study requires first obtaining genetic markers as well as expression profiles for each individual in the population, then taking the expression of each target gene as a trait (termed an expression trait, eTrait) and determining whether some markers are statistically associated with the eTrait by association analysis, and finally identifying candidate genes or regulatory sequences around the associated markers through various additional evidence. Usually, eQTL can be classified into two types: cis and trans. A cis-acting eQTL is an eQTL for a gene that is localized around that gene, indicating that sequence differences around that gene result in changes in expression levels. A trans-acting eQTL is an eQTL that is positioned distantly from the target gene it regulates, indicating that the expression level of the target gene is genetically regulated by distal factors (e.g., upstream transcription factors). 
>
In the following protocol, we explain how to use QTLtools to identify cis- and trans- eQTL and use qqman to visualize the results.To guide eBook authors having a better sense of the workflow layout, here we briefly introduce the specific purposes of the dir system.

1. **cache:** Here, it stores intermediate datasets or results that are generated during the preprocessing steps.
2. **graphs:** The graphs/figures produced during the analysis.
3. **input:** Here, we store the raw input data.
4. **lib:** The source code, functions, or algorithms used within the workflow.
5. **output:** The final output results of the workflow.
6. **workflow:** Step by step pipeline. 
# Installation
## Required software and installation:
### Installing [Anaconda](https://www.anaconda.com/) 
- wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/Anaconda3-4.1.0-Linux-x86.sh
- bash Anaconda3-4.1.0-Linux-x86.sh
- echo 'export PATH="~/anaconda3/bin:$PATH"' >> ~/.bashrc
- source ~/.bashrc
### conda environment
- conda create -n eQTL_Analysis_for_Rice_RIL_population
- conda activate eQTL_Analysis_for_Rice_RIL_population
### Installing [Htslib](https://github.com/samtools/htslib/)
- conda install -c bioconda htslib
### Installing [Bcftools](https://github.com/samtools/bcftools/)
- conda install -c bioconda bcftools
### Installing [Samtools](https://github.com/samtools/samtools/)
- conda install -c bioconda samtools 
### Installing [R 3.6.1](http://www.R-project.org/)
- conda install r-base=3.6.1
### Installing [QTLtools](https://qtltools.github.io/qtltools/)(download and unzip to use)
- wget https://qtltools.github.io/qtltools/binaries/QTLtools_1.2_CentOS7.8_x86_64.tar.gz
- tar xzvf QTLtools_1.2_CentOS7.8_x86_64.tar.gz
- cd QTLtools_1.2_CentOS7.8_x86_64
- ln -s QTLtools_1.2_CentOS7.8_x86_64 QTLtools
- echo 'export PATH="~/QTLtools_1.2_CentOS7.8_x86_64:$PATH"' >> ~/.bashrc
- source ~/.bashrc
# Input Data
The raw data are available from the National Center for Biotechnology Information Gene Expression Omnibus database under the accession number GSE49020.
- a.	Genotype data（VCF/BCF format）:eQTL_genotype.vcf
- b.	eTrait/phenotype data（BED format）:flag_leaf_eTrait.bed
# Major steps
### Step 1: compressed and indexed raw data
       sh Compression.sh
- Script content of Compression.sh
```ruby 
bgzip eQTL_genotype.vcf
tabix vcf eQTL_genotype.vcf.gz
bgzip flag_leaf_eTrait.bed 
tabix flag_leaf_eTrait.bed.gz
``` 
### Step 2: cis-eQTL identification with QTLtools
       sh QTLtools_cis.sh
- Script content of QTLtools_cis.sh
```ruby 
QTLtools cis --vcf eQTL_genotype.vcf.gz --bed flag_leaf_eTrait.bed.gz --permute 1000 --out flag_leaf_eTrait_cis_permutation.txt > running.log
``` 
### Note:We can also use --norminal to enforce the eTraits to be normally distributed,for example:
```ruby 
QTLtools cis --vcf eQTL_genotype.vcf.gz --bed flag_leaf_eTrait.bed.gz --norminal 0.05 --out flag_leaf_eTrait_cis_permutation.txt > running.log
``` 

### Step 3: trans-eQTL identification with QTLtools
The command 1 will generate three output files, one named “*.best.txt.gz” containing the top eQTL for each eTrait, one named “*.bins.txt.gz” containing all eQTLs with P-values below the specified threshold, and the last named “*.hits.txt.gz”, containing the details of all eQTLs with P-values above the specified threshold. The command 2 permutes all eTraits and generate three files like the command 1. The command 3 will generate the file “flag_leaf_trans_005_permutations_all.txt” which contains the data in “*.hits.txt.gz” and with an additional column that gives the estimated false discovery rate (FDR) for each eTrait by 100 permutations.
>
       sh QTLtools_trans.sh
- Script content of QTLtools_trans.sh
```ruby 
##command 1 
QTLtools trans --vcf eQTL_genotype.vcf.gz --bed flag_leaf_eTrait.bed.gz --nominal --threshold 0.05 --out flag_leaf005.trans.nominal.hits.txt.gz > running.log

##command 2
for i in {1..100};do
       QTLtools trans --vcf eQTL_genotype.vcf.gz --bed flag_leaf_eTrait.bed.gz --threshold 0.05 --permute --out flag_leaf005_trans_perm_${i} --seed ${i} > running.log
done

##command 3
zcat flag_leaf005_trans_perm_*.hits.txt.gz | gzip -c > flag_leaf005_permutations_all.txt.gz
Rscript runFDR_ftrans.R flag_leaf005.trans.nominal.hits.txt.gz flag_leaf005_permutations_all.txt.gz flag_leaf_trans_005_permutations_all.txt
``` 
### Step 4: Draw Manhattan diagram with R script
       sh get_Ehd1_eQTL_result.sh
- Script content of get_Ehd1_eQTL_result.sh
```ruby 
zcat flag_leaf005.trans.nominal.hits.txt.gz | grep OsAffx.30643.1.S1_at > Ehd1_eQTL_result.txt
``` 
        Rscript Manhattan.R
- Script content of Manhattan.R
```ruby 
library('qqman')
data <- read.table(" Ehd1_eQTL_result.txt ")
data <-data[,c(4,5,6,7)]
colnames(data)<-c('SNP','CHR','BP','P')
manhattan(data,main = "Manhattan Plot",suggestiveline =FALSE,genomewideline = FALSE,csi=1.5,cex.lab=1.5, cex.axis=1.5,cex.main=1.5,annotatePval = 5e-40, annotateTop = FALSE,xlim=c(32284500,155000000),ylim=c(0,50),cex = 1.5)
abline(h=-log10(6.17e-04), col="blue", lty=1, lwd=3) 
abline(h=-log10(6.12e-04), col="red", lty=2, lwd=3)
``` 
# Expected results
### You will get cis-eQTL results and trans-eQTL results files
- flag_leaf_eTrait_cis_permutation.txt
- flag_leaf005.trans.nominal.hits.txt.gz 
- flag_leaf_trans_005_permutations_all.txt（ Compare to the “flag_leaf005.trans.nominal.hits.txt.gz”,this file has an additional column that gives the estimated false discovery rate (FDR) for each eTrait by 100 permutations.）
### Manhattan plot of OsAffx.30643.1.S1_at
![图片名称](https://github.com/ziongfen/protocol/blob/main/graphs/Rplot.png)
# License
It is a free and open source software, licensed under (choose a license from the suggested list: [GPLv3](https://github.com/github/choosealicense.com/blob/gh-pages/_licenses/gpl-3.0.txt), [MIT](https://github.com/github/choosealicense.com/blob/gh-pages/LICENSE.md), or [CC BY 4.0](https://github.com/github/choosealicense.com/blob/gh-pages/_licenses/cc-by-4.0.txt)).
