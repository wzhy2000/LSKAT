\name{plot.LSKAT.gen.plink}
\alias{plot.LSKAT.gen.plink}
\title{ Plot Manhattan figure.
}
\description{
This function draws manhattan figure using the p-value of LSKAT and L-Burden. 
}
\usage{
plot.LSKAT.gen.plink( r.lskat, pdf.file, title, y.max, bonferroni )
}
\arguments{
  \item{r.lskat}{The list object with S3 class name: LSKAT.gen.plink, obtained from \code{\link{longskat_gene_plink}}.  }
  \item{pdf.file}{String, the file name of PDF output.}
  \item{title}{String, title}
  \item{y.max}{Numeric, the maximum value of Y-axis}
  \item{bonferroni}{Boolean, indicating whether Bonferroni correction is used in the plot}
}
\seealso{
\code{\link{longskat_gene_plink}}
}
\examples{
## check the code in the function longskat_gene_plink
}
\keyword{ Plot }
\keyword{ Gene-base }
\keyword{ PLINK }