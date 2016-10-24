\name{longskat_snp_test}
\alias{longskat_snp_test}
\title{
None
}
\description{
None
}
\usage{
longskat_snp_test(r.model, 
    snp, 
    weights.common = c(0.5, 0.5),
    weights.rare = c(1, 25), 
    snp.impute = "mean",
    rare.cutoff = NULL, 
    run.cpp = F, 
    verbose = FALSE )
}
\arguments{
  \item{r.model}{The list object obtained from the function \code{\link{longskat_est_model}}, including the estimated parameters and residuals.}
  \item{snp}{Matrix with m row for individuals and n columns for the variaints(SNPs), also with the individuals' ID as the row names.}
  \item{weights.common}{a numeric vector of parameters of beta weights for common variants (default=c(0.5,0.5)).}
  \item{weights.rare}{a numeric vector of parameters of beta weights for rare variants (default=c(1,25)).}
  \item{rare.cutoff}{Numeric, a value of MAF cutoff for the rare SNPs. Only SNPs that have MAFs smaller than this are considered as rare SNP. The default criterion of rare SNP is calculated by the formula \eqn{1/\sqrt{2*sample}} }
  \item{run.cpp}{Logical, indicating whether C/C++ functions are used to compute LSKAT.}
  \item{verbose}{Logical variable, indicating whether some debug information can be outputted.}

}
\value{
None
}
\examples{

## data simulation for the power test
p0 <- longskat_gene_simulate( plink.format=F, power.test=T, n.gene=1);

## model estimation
r.model <- longskat_est_model( p0$phe.long, p0$phe.cov, g.maxiter=3, verbose=T);

for(i in 1:NCOL(p0$snp.mat[[1]]))
{
    ## test all SNPs in the 1st gene.
    r.lskat <- longskat_snp_test(r.model, p0$snp.mat[[1]][,i]);
    print(r.lskat);
}


}
\keyword{ SNP-base }
\keyword{ Test }
