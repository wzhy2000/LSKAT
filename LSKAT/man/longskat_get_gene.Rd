\name{longskat_get_gene}
\alias{longskat_get_gene}
\title{
Get SNP matrix of gene.
}
\description{
This function provides a interface to read all SNPs of genes from the PLINK data object. 
}
\usage{
longskat_get_gene(  gen.obj, gene.set, snp.impute="mean", verbose = FALSE )
}
\arguments{
  \item{gen.obj}{Reference class , it is a wrapper of plink data object obtained from the function \code{\link{longskat_plink_load}}. }
  \item{gene.set}{Vector of numeric or string, indicating the gene index or gene name.}
  \item{snp.impute}{String, indicating the method of SNP imputation, the default model uses the mean of each variant to replace the missing SNP data. Two optional values: 'mean' or 'random'. }
  \item{verbose}{Logical variable, indicating whether some debug information can be outputted. };
}

\details{
The parameter \code{gene.set} can be gene index (numeric) or gene names (string).
}

\value{

The returned list object could contain multiple genes. It includes 

\item{snp.mat}{List, SNP Matrices for genes in \code{gene.set}. }
\item{maf}{List, SNP MAF vectors for genes in \code{gene.set}. }
\item{nmiss}{List, SNP missing information for genes in \code{gene.set}. }
\item{gene.name}{List, gene names for genes in \code{gene.set}. }

}
\references{
Wang Z., Xu K., Zhang X., Wu X., and Wang Z., (2016) Longitudinal SNP-set association analysis of quantitative phenotypes. Genetic Epidemiology.
}

\seealso{
\code{\link{longskat_plink_load}} and \code{\link{longskat_gene_test}}
}
\examples{

## data simulation for the power test
p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-plink-load", 
      power.test=T, n.gene=20);

## Loading PLINK data set into the reference object
plk <- longskat_plink_load ( p0$file.plink.bed, p0$file.plink.bim, p0$file.plink.fam, 
      p0$file.gene.set,  verbose=TRUE);

## Retriving the SNP matrix from the reference object
snps <- longskat_get_gene( plk, 1:2, snp.impute="mean" );

## show the structures of the result.
str(snps);


}
\keyword{ PLINK }

