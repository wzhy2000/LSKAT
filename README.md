# LSKAT

LSKAT is an R program that performs association testing between a set of SNPs and a longitudinal quantitative trait in population samples. The main reference for this program is

# Reference

Wang Z., Xu K., Zhang X., Wu X., and Wang Z., (2016) Longitudinal SNP-set association analysis of quantitative phenotypes. Genetic Epidemiology.

## Abstract:

To do ...

## Installation Instructions:

### Required software and packages
    
1. R (http://www.r-project.org/)
    
2. Package SKAT, mvtnorm, snpStats, snowfall.
    
3. PLINK 1.0 or 1.9 (https://www.cog-genomics.org/plink2)

Please install the required R package before you install LSKAT package. After the  installation of `SKAT`, `mvtnorm`, `snpStats` and `snowfall` package, please install the **LSKAT** as following steps.

 
### Install LSKAT on LUNIX or Mac OSX

```
git clone https://github.com/ZWang-Lab/LSKAT.git

cd LSKAT

R CMD INSTALL LSKAT

```

### Install LSKAT on Windows

1) Please download windows package from (https://github.com/wzhy2000/LSKAT/blob/master/windows/LSKAT_0.5.zip)

2) Install the package in R GUI by selecting the menu "Packages|Install package(s) from local zip files..."

##Usage instructions##

LSKAT is an R package which provides the main functions:

> 1) NULL model estimation to estimate the model parameters and the residuals.

> 2) Gene association test by LSKAT

> 3) LSKAT pipeline for PLINK data set with longitudial phenotype traits.


In general, two ways are common to do data analysis by LSKAT, one is test the genetic association between the longitudinal phenotype traits and SNP-set genotypic matrix, i.e, test LSKAT using the SNP-set matrix given the estimated NULL model which assumes no genetic effects contibute to phenotype traits.


```
## NULL model estimation
r.nodel <- longskat_est_model( phe.long.matrix, phe.cov.matrix, phe.time.matrix);

## Gene association test
r.lskat <- longskat_gene_test( r.model, snp.mat);
```

Other way is to run LSKAT on whole PLINK data set using the pipeline provided by the function `longskat_gene_plink`


```
r.lskat <- longskat_gene_plink( file.plink.bed,  file.plink.bim,  file.plink.fam,  
    file.phe.long,  file.phe.cov, NULL,  file.gene.set );
```

All functions and examples in the LSKAT are available in the manual (https://github.com/ZWang-Lab/LSKAT/blob/master/manual.pdf).
