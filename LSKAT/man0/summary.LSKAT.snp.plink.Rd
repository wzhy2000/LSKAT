\name{summary.LSKAT.snp.plink}
\alias{summary.LSKAT.snp.plink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.LSKAT.snp.plink(r.lskat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r.lskat}{
%%     ~~Describe \code{r.lskat} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (r.lskat) 
{
    cat("  LSKAT/SNP Summary\n")
    cat("[1] Paramaters:\n")
    cat("* PHE.LOG.FILE = ", r.lskat$par$file.phe.long, "\n")
    cat("* PHE.LOG.TIME = ", r.lskat$par$file.phe.time, "\n")
    cat("* PHE.COV.FILE = ", r.lskat$par$file.phe.cov, "\n")
    cat("* Covariate Count = ", r.lskat$par$y.cov.count, "\n")
    cat("* Time Effect = ", r.lskat$par$y.cov.time, "\n")
    cat("* Weight Rare = ", r.lskat$par$w.rare, "\n")
    cat("* Weight Common = ", r.lskat$par$w.common, "\n")
    cat("[2] MLE Results:\n")
    cat("* SIGMA_A =", r.lskat$mle$par$sig_a, "\n")
    cat("* SIGMA_B =", r.lskat$mle$par$sig_b, "\n")
    cat("* SIGMA_E =", r.lskat$mle$par$sig_e, "\n")
    cat("* RHO =", r.lskat$mle$par$rho, "\n")
    cat("* MU =", r.lskat$mle$par$mu, "\n")
    cat("* Beta(Cov)=", r.lskat$mle$par$par_cov, "\n")
    cat("* Beta(Time)=", r.lskat$mle$par$par_t, "\n")
    cat("* L(min) =", r.lskat$mle$likelihood, "\n")
    cat("[3] Longitudinal SKAT Results:\n")
    cat("* SNP Count = ", NROW(r.lskat$snp), "\n")
    df.snp <- r.lskat$snp[, c(1, 2, 3, 4, 5, 10), drop = F]
    df.snp <- df.snp[order(df.snp[, 6]), , drop = F]
    df.snp.bonf <- df.snp[, 6] * NROW(df.snp)
    n.sig <- c()
    n.sig.i <- c()
    for (i in 20:4) {
        n <- length(which(df.snp.bonf < 10^(-i)))
        if (n > 0) {
            n.sig <- c(n.sig, n)
            n.sig.i <- c(n.sig.i, i)
        }
    }
    cat("  Summary of LSKAT result:\n")
    cat("* PHE.LOG.FILE = ", r.lskat$par$file.phe.long, "\n")
    cat("* PHE.COV.FILE = ", r.lskat$par$file.phe.cov, "\n")
    if (length(n.sig.i) > 0) {
        cat("* Significant Genes (Bonferroni Correct):\n")
        for (i in 1:length(n.sig.i)) cat("* <= 10^(-", n.sig.i[i], 
            ") Levels:", n.sig[i], "\n")
    }
    cat("* Top 20 SNPs (Bonferroni Correct):\n")
    df.show <- df[c(1:20), c(3, 4, 1, 2, 5, 6)]
    df.show <- cbind(df.show, df.show[, 6] * NROW(df.snp))
    colnames(df.show) <- c("SNP", "Gene", "Chr.", "Pos.", "MAF", 
        "p-value", "Bonferroni")
    show(df.show)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
