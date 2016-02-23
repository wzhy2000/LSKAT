\name{show.LSKAT.gen.ret}
\alias{show.LSKAT.gen.ret}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
show.LSKAT.gen.ret(r.lskat)
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
    cat("  LSKAT/Gene Summary\n")
    cat("[1] MLE Results:\n")
    cat("* SIGMA_A =", r.lskat$mle$par$sig_a, "\n")
    cat("* SIGMA_B =", r.lskat$mle$par$sig_b, "\n")
    cat("* SIGMA_E =", r.lskat$mle$par$sig_e, "\n")
    cat("* RHO =", r.lskat$mle$par$rho, "\n")
    cat("* MU =", r.lskat$mle$par$mu, "\n")
    cat("* Beta(Cov)=", r.lskat$mle$par$par_cov, "\n")
    cat("* Beta(Time)=", r.lskat$mle$par$par_t, "\n")
    cat("* L(min) =", r.lskat$mle$likelihood, "\n")
    cat("[2] LSKAT:\n")
    cat("* SNP total = ", r.lskat$snp.total, "\n")
    cat("* Rare SNPs = ", r.lskat$snp.rare, "\n")
    cat("* Q value = ", r.lskat$qv, "\n")
    cat("* p-value = ", r.lskat$pv, "\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
