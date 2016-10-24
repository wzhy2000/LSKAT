\name{longskat_plink_load}
\alias{longskat_plink_load}
\title{
Loading plink data set.
}
\description{
This function provides a pipeline for the plink data set based on the LSKAT function. Except the plink data set(BED, BIM and FAM), SNP set table is same as SKAT package. 
}
\usage{
longskat_plink_load(file.plink.bed, 
    file.plink.bim, 
    file.plink.fam,
    file.gene.set,
    plink.path=NULL,
    verbose=FALSE)
}
\arguments{
  \item{file.plink.bed}{File name, PLINK bed file, containing the packed binary SNP genotype data}
  \item{file.plink.bim}{File name, PLINK bim file, containing the SNP descriptions}
  \item{file.plink.fam}{File name, PLINK fam file, containing subject(and, possibly, family) identifiers}
  \item{file.gene.set}{File name, indicating SNP set table which has two columns, gene name in 1st column and variant name in 2nd column and no header. }
  \item{plink.path}{String, indicating PLINK command path, for the large PLINK data, the package will load partial plink data extracted by the plink command rather than the whole data set. }
  \item{verbose}{Logical Variable, indicating whether some information are outputted for debug.}
}

\details{
The SNP set file is same as the \code{SetID} file in \emph{SKAT} package which is white-space (space or tab) seperated file with 2 columns: SetID (1st column) and SNP_ID (2nd column), Please keep in mind that there should be \emph{no header!} \cr
}

\value{

This function returns a reference class ("\code{'PLINK.refer'}")for PLINK data operation which provides a interface to access the SNP information. \code{\link{longskat_get_gene }} use this object to extract the SNP matrix for genes.

}
\references{
Wang Z., Xu K., Zhang X., Wu X., and Wang Z., (2016) Longitudinal SNP-set association analysis of quantitative phenotypes. Genetic Epidemiology.
}

\seealso{
\code{\link{longskat_get_gene}}
}
\examples{
## see the example in the function longskat_get_gene
}
\keyword{ PLINK }
