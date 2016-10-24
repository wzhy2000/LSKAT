\name{longskat_snp_plink}
\alias{longskat_snp_plink}
\title{
Single SNP associaion test using LSKAT for PLINK data set.
}
\description{
This function provides a pipeline for the plink data set based on the LSKAT function. Except the plink data set(BED, BIM and FAM), SNP set table is same as SKAT package. 
}
\usage{
longskat_snp_plink(file.plink.bed, 
    file.plink.bim, 
    file.plink.fam,
    file.phe.long, 
    file.phe.cov,
    file.phe.time = NULL, 
    file.gene.set=NULL,
    snp.set = NULL, 
    options = list(),
    verbose=FALSE);
}
\arguments{
  \item{file.plink.bed}{File name, PLINK bed file}
  \item{file.plink.bim}{File name, PLINK bim file}
  \item{file.plink.fam}{File name, PLINK fam file}
  \item{file.phe.long}{File name, indicating phenotype file in CSV format with row name (individual ID) and header information (Measured Index). Each individual is encoded to one row data which has individual ID as row name and multiple observed values followed by row name.}
  \item{file.phe.cov}{File name, indicating covariate file in CSV format with row name  (individual ID) and header information (Covariate Index). Each individual is encoded to one row data which has individual ID as row name and covariate values followed by row name.}
  \item{file.phe.time}{File name, indicating measured times in CSV format with row name  (individual ID) and header information (Measured Time). Each individual is encoded to one row data which has individual ID as row name and covariate values followed by row name. }
  \item{file.gene.set}{File name, indicating gene set table which has two columns, 1st column is gene name, 2nd column is variant name and no header. }
  \item{snp.set}{numeric or string. indicating the gene name or gene index in the gene set table}
  \item{options}{String, see the \emph{details} }
  \item{verbose}{Logical variable, indicating whether some debug information can be outputted. }
}
\details{
None
}
\value{
A list object with 3 items is returned by this function:

\item{par}{The parameters of this function call.}
\item{mle}{The estimation of NULL Model returned by \code{\link{longskat_est_model}}, including model parameters and residuals, please refer to the structure in \code{\link{longskat_est_model}}.  }
\item{result}{Matrix containing the test result of LSKAT for all SNPs, the sub-items are described as follows.}

The structure: of \code{result} item:\cr
  \item{index}{Number, gene index }
  \item{snp.name}{String, snp name }
  \item{gene.name}{String, gene name }
  \item{chr}{Number, chromosome }
  \item{pos }{Number, position }
  \item{MAF }{Number, Minor Allel Frequency.}
  \item{NMISS }{Number, number of missing SNP.}
  \item{rare }{Logical, rare SNP or not.}
  \item{q.lskat }{Numeric, the statistical test of LSKAT }
  \item{p.lskat }{Numeric, the p-value of LSKAT }


}
\examples{

## data simulation for the power test.
p1 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-simu-snp", 
      power.test=T, n.gene=2);

## Single SNP test by LSKAT 
r.lskat <- longskat_snp_plink( p1$file.plink.bed, p1$file.plink.bim, p1$file.plink.fam, 
      p1$file.phe.long, p1$file.phe.cov, NULL, p1$file.gene.set );

## show the result
print(r.lskat);

## draw the result in the Manhattan figure.
plot(r.lskat, "tmp-simu-snp.pdf", "title")

}
\keyword{ SNP-base }
\keyword{ PLINK }
\keyword{ Test }
