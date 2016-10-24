\name{plot.LSKAT.snp.plink}
\alias{plot.LSKAT.snp.plink}
\title{
Plot Manhattan figure using SNP association test
}
\description{
This function draws manhattan figure using the p-value of LSKAT and L-Burden. 
}
\usage{
plot.LSKAT.snp.plink(r.lskat, pdf.file, title, y.max)
}
\arguments{
  \item{r.lskat}{The list object with S3 class name: LSKAT.gen.plink, obtained from \code{\link{longskat_snp_plink}}.  }
  \item{pdf.file}{String, the file name of PDF output.}
  \item{title}{String, title}
  \item{y.max}{Numeric, the maximum value of Y-axis}
  \item{bonferroni}{Boolean, indicating whether Bonferroni correction is used in the plot}
}
\seealso{
\code{\link{longskat_snp_plink}}
}
\examples{
## check the code in the function longskat_snp_plink
}
\keyword{ Plot }
\keyword{ SNP-base }
\keyword{ PLINK }
