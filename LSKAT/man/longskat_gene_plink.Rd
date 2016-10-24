\name{longskat_gene_plink}
\alias{longskat_gene_plink}
\title{
LSKAT test for plink data set.
}
\description{
This function provides a pipeline for the plink data set based on the LSKAT function. Except the plink data set(BED, BIM and FAM), SNP set table is same as SKAT package. 
}
\usage{
longskat_gene_plink(file.plink.bed, 
    file.plink.bim, 
    file.plink.fam, 
    file.phe.long, 
    file.phe.cov, 
    file.phe.time = NULL, 
    file.gene.set, 
    gene.set = NULL, 
    options = list(),
    verbose = FALSE)
}
\arguments{
  \item{file.plink.bed}{File name, PLINK bed file}
  \item{file.plink.bim}{File name, PLINK bim file}
  \item{file.plink.fam}{File name, PLINK fam file}
  \item{file.phe.long}{File name, indicating phenotype file in CSV format with row name (individual ID) and header information (Measured Index). Each individual is encoded to one row data which has individual ID as row name and multiple observed values followed by row name.}
  \item{file.phe.cov}{File name, indicating covariate file in CSV format with row name  (individual ID) and header information (Covariate Index). Each individual is encoded to one row data which has individual ID as row name and covariate values followed by row name.}
  \item{file.phe.time}{File name, indicating measured times in CSV format with row name  (individual ID) and header information (Measured Time). Each individual is encoded to one row data which has individual ID as row name and covariate values followed by row name. }
  \item{file.gene.set}{File name, indicating gene set table which has two columns, 1st column is gene name, 2nd column is variant name and no header. }
  \item{gene.set}{numeric or string. indicating the gene name or gene index in the gene set table.}
  \item{options}{String, see the details. }
  \item{verbose}{Logical variable, indicating whether some debug information can be outputted. }
}

\details{

The gene set file is same as the SetID file in SKAT packa which is  white-space (space or tab) seperated file with 2 columns: SetID (1st column) and SNP_ID (2nd column), Please keep in mind that there should be no header! \cr

The phenotype file, covariate file and measured time file are CSV files with header and row name(indicating the individual id), these 3 files should have same individual order. Please note that the columns in measured time file are corresponding to same columns in the phenotype file.

The \code{options} list has some important feature to tune the association test. The default values are defined in the package as the follows:

\preformatted{
 options <- list( rare.cutoff = NULL,
                time.cov      = 0,
                g.maxiter     = 20,
                weights.common= c(0.5,0.5),
                weights.rare  = c(1,25),
                run.cpp       = F,
                verbose       = F,
                n.cpu         = 1,
                snp.impute    = "mean",
                intercept     = F,
                plink.path    = NULL,
                test.type     = "Joint",
                est.method    = "REML");
}

\describe{
  \item{rare.cutoff}{Numeric, a value of MAF cutoff for the rare SNPs. Only SNPs that have MAFs smaller than this are considered as rare SNP. The default criterion of rare SNP is calculated by the formula \eqn{1/\sqrt{2*sample}} }
  \item{time.cov}{Numeric, indicating which order of time is included as extra covariates. If this value is 2, time and time square are considered as extra covariates.}
  \item{g.maxiter}{Number, indicating the number of maximum iteration is applied to MLE alogrithm. }
  \item{weights.common}{a numeric vector of parameters of beta weights for common variants (default=c(0.5,0.5)).}
  \item{weights.rare}{a numeric vector of parameters of beta weights for rare variants (default=c(1,25)).}
  \item{run.cpp}{Logical, indicating whether C/C++ functions are used to compute LSKAT.}
  \item{verbose}{Logical, indicating the computaional process output much more information than normal mode.}
  \item{snp.impute}{String, indicating the method of SNP imputation, the default model uses the mean of each variant to replace the missing SNP data. Two optional values: 'mean' or 'random'.  }
  \item{intercept}{Logical, indicating whether the intercept is considered in the NULL model.}
  \item{plink.path}{String, indicating PLINK command path, for the large PLINK data, the package will load small plink data set extracted by the plink command rather than the whole data set. }
  \item{test.type}{String, Three models can be selected, "joint", "Common.Only", "rare.Only". }
  \item{est.method}{String, indicating the estimate method for NULL model, two options, REML or ML, default is REML.}
}

}
\value{
A list object with 3 items is returned by this function:

\item{par}{The parameters of this function call.}
\item{mle}{The estimation of NULL Model returned by \code{\link{longskat_est_model}}, including model parameters and residuals, please refer to the structure in \code{\link{longskat_est_model}}.  }
\item{result}{Matrix containing the test result of LSKAT and L-burden for all genes, the sub-items are described as follows.}

The structure: of \code{result} item:\cr
  \item{index}{Number, gene index }
  \item{gene.name}{String, gene name }
  \item{chr}{Number, chromosome }
  \item{min.pos }{Number, minimum position of all variants }
  \item{snp.total }{Number, the total number of variants }
  \item{snp.rare }{Number, the number of rare variants.}
  \item{q.lskat }{Numeric, the statistical test of LSKAT }
  \item{p.lskat }{Numeric, the p-value of LSKAT }
  \item{q.burden }{Numeric, the statistical test of L-Burden test }
  \item{p.burden }{Numeric, the p-value of L-Burden test }

}
\references{
Wang Z., Xu K., Zhang X., Wu X., and Wang Z., (2016) Longitudinal SNP-set association analysis of quantitative phenotypes. Genetic Epidemiology.
}
\seealso{
\code{\link{longskat_gene_test}} and \code{\link{longskat_est_model}}
}
\examples{

## data simulation for the test of type 1 error 
p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-gene-plink", 
      power.test=F, n.gene=20);

## test all genes in the PLINK data set 
r.lskat0 <- longskat_gene_plink(p0$file.plink.bed, p0$file.plink.bim, p0$file.plink.fam, 
      p0$file.phe.long, p0$file.phe.cov, NULL, 
      p0$file.gene.set, gene.set=c(5:10), verbose=T );

## print the results
r.lskat0;

## data simulation for the power test
p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-gene-plink", 
      power.test=T, n.gene=20);

## test all genes in the PLINK data set 
r.lskat1<-longskat_snp_plink( p0$file.plink.bed, p0$file.plink.bim, p0$file.plink.fam, 
      p0$file.phe.long, p0$file.phe.cov, file.gene.set=p0$file.gene.set, 
      snp.set=c(1:199),
      options=list(time.cov=1, g.maxiter=6, n.cpu=2, test.type="Rare.only" ));

## print the results
r.lskat1;

#draw the Manhattan figure for the result
plot( r.lskat1, pdf.file="tmp-plink-lskat1.pdf", title="TEST PLINK 1", bonferroni=T);


}
\keyword{ PLINK }
\keyword{ Test }
\keyword{ Gene-base }
