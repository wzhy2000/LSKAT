\name{summary.LSKAT.gen.plink}
\alias{summary.LSKAT.gen.plink}
\title{
Summarizing the LSKAT results for all genes
}
\description{
Summarizing the LSKAT results for all genes in PLINK data set
}
\usage{
summary.LSKAT.gen.plink(r.lskat)
}
\arguments{
  \item{r.lskat}{The list object with S3 class name: LSKAT.gen.plink, obtained from \code{\link{longskat_gene_plink}}.  }
}
\seealso{
\code{\link{longskat_gene_plink}}
}
\examples{

## data simulation for the power test
p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-gene-plink", 
      power.test=T );

## LSKAT test
r.lskat1 <- longskat_gene_plink(p0$file.plink.bed, p0$file.plink.bim, p0$file.plink.fam, 
      p0$file.phe.long, p0$file.phe.cov, NULL, p0$file.gene.set, 
      options=list(g.maxiter=3, plink.path="plink", est.method = "ML"), verbose=T );

## sumamry information, this function
df <- summary(r.lskat1);

## check the sumamry information
show(df);

}
\keyword{ Gene-base }
\keyword{ Summary }
\keyword{ PLINK }