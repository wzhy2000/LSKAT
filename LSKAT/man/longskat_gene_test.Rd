\name{longskat_gene_test}
\alias{longskat_gene_test}
\title{
Assoictaion test using LSKAT
}
\description{
Assoictaion test for the comined effect of common and rare variants using LSKAT
}
\usage{

longskat_gene_test( r.model, 
    snp.mat, 
    weights.common=c(0.5,0.5), 
    weights.rare=c(1,25),
    rare.cutoff=NULL, 
    test.type="Joint", 
    snp.impute = "mean",
    run.cpp=F, 
    verbose=F)
}
\arguments{
  \item{r.model}{The list object obtained from the function \code{\link{longskat_est_model}}, including the estimated parameters and residuals.}
  \item{snp.mat}{Matrix with m row for individuals and n columns for the variaints(SNPs), also with the individuals' ID as the row names.}
  \item{weights.common}{a numeric vector of parameters of beta weights for common variants (default=c(0.5,0.5)).}
  \item{weights.rare}{a numeric vector of parameters of beta weights for rare variants (default=c(1,25)).}
  \item{rare.cutoff}{Numeric, a value of MAF cutoff for the rare SNPs. Only SNPs that have MAFs smaller than this are considered as rare SNP. The default criterion of rare SNP is calculated by the formula \eqn{1/\sqrt{2*sample}} }
  \item{test.type}{String, Three models can be selected, "joint", "Common.Only", "rare.Only". }
  \item{snp.impute}{String, indicating the method of SNP imputation, the default model uses the mean of each variant to replace the missing SNP data. }
  \item{run.cpp}{Logical, indicating whether C/C++ functions are used to compute LSKAT.}
  \item{verbose}{Logical variable, indicating whether some debug information can be outputted.}
}
\details{

} 
\value{
The list object is returned by this function with the following items:\cr
  \item{snp.NMISS }{Vector, the missing rate for each SNP.}
  \item{snp.MAF }{Vector, the MAF for each SNP.}
  \item{snp.total }{Numeric, the total number of variants. }
  \item{snp.rare }{Numeric, the number of rare variants.}
  \item{q.lskat }{Numeric, the statistical test of LSKAT }
  \item{p.lskat }{Numeric, the p-value of LSKAT }
  \item{q.burden }{Numeric, the statistical test of L-Burden Test }
  \item{p.burden }{Numeric, the p-value of L-Burden Test }  
}
\references{
Wang Z., Xu K., Zhang X., Wu X., and Wang Z., (2016) Longitudinal SNP-set association analysis of quantitative phenotypes. Genetic Epidemiology.

}
\note{
}
\seealso{
\code{\link{longskat_est_model}}
}

\examples{

## data simulation for the power test
p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-gene-test", 
      power.test=T, n.gene=5 );

## model estimation
r.reml <- longskat_est_model( p0$phe.long, p0$phe.cov, phe.time = NULL, g.maxiter=3,  
      method="REML", verbose=T )

print(r.reml);

## test the genes using the function 'longskat_gene_test'
for(i in 1:length(p0$snp.mat))
{
    r.lskat0 <- longskat_gene_test(r.reml, p0$snp.mat[[i]], snp.impute="mean");
    print(r.lskat0);
}

}


\keyword{ Gene-base }
\keyword{ Test }
