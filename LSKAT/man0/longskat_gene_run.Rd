\name{longskat_gene_run}
\alias{longskat_gene_run}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
longskat_gene_run(r.model, snp.mat, weights.rare, weights.common, run.cpp = T, debug = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r.model}{
%%     ~~Describe \code{r.model} here~~
}
  \item{snp.mat}{
%%     ~~Describe \code{snp.mat} here~~
}
  \item{weights.rare}{
%%     ~~Describe \code{weights.rare} here~~
}
  \item{weights.common}{
%%     ~~Describe \code{weights.common} here~~
}
  \item{run.cpp}{
%%     ~~Describe \code{run.cpp} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
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
function (r.model, snp.mat, weights.rare, weights.common, run.cpp = T, 
    debug = F) 
{
    Y.delt <- r.model$y.delt
    Y.time <- r.model$y.time
    par_null <- c(r.model$par$sig_a, r.model$par$sig_b, r.model$par$sig_e, 
        r.model$par$rho)
    X <- matrix(cbind(1, r.model$y.cov))
    Z.scale <- SKAT_Scale_Genotypes(X, t(snp.mat), weights.rare = weights.rare, 
        weights.common = weights.common)
    Q <- try(est_gen_Q(Y.delt, Y.time, Z.scale$maf, Z.scale$new, 
        X, par_null, r.model$y.cov.time, run.cpp = run.cpp))
    if (class(Q) == "try-error") 
        return(NULL)
    P <- get_Q_pvale(Q$v, Q$w)
    if (debug) 
        cat(" MAF=", (Z.scale$maf), "\n")
    r.lskat <- list(mle = r.model, snp.total = length(Z.scale$maf), 
        snp.rare = Z.scale$rare, qv = Q$v, pv = P$p.value)
    class(r.lskat) <- "LSKAT.gen.ret"
    return(r.lskat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
