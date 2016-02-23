\name{longskat_gene_simulate}
\alias{longskat_gene_simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
longskat_gene_simulate(g.snp.hap1, g.snp.pos2, n.minsect = 3000, n.maxsect = 30000, n.sample = 800, par = list(a = 0.5, b = 0.5, sig_a = 0.8, sig_b = 0.8, sig_e = 0.8, beta.effect = c(0.1, 0.1, 0, 0.8), c1 = c(0.2828, 0.2828, 0, 0), times = 6, par1 = 0.7, par2 = 0.4), covariance = list(phe.dist = "mn", phe.cov = "ar1"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g.snp.hap1}{
%%     ~~Describe \code{g.snp.hap1} here~~
}
  \item{g.snp.pos2}{
%%     ~~Describe \code{g.snp.pos2} here~~
}
  \item{n.minsect}{
%%     ~~Describe \code{n.minsect} here~~
}
  \item{n.maxsect}{
%%     ~~Describe \code{n.maxsect} here~~
}
  \item{n.sample}{
%%     ~~Describe \code{n.sample} here~~
}
  \item{par}{
%%     ~~Describe \code{par} here~~
}
  \item{covariance}{
%%     ~~Describe \code{covariance} here~~
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
function (g.snp.hap1, g.snp.pos2, n.minsect = 3000, n.maxsect = 30000, 
    n.sample = 800, par = list(a = 0.5, b = 0.5, sig_a = 0.8, 
        sig_b = 0.8, sig_e = 0.8, beta.effect = c(0.1, 0.1, 0, 
            0.8), c1 = c(0.2828, 0.2828, 0, 0), times = 6, par1 = 0.7, 
        par2 = 0.4), covariance = list(phe.dist = "mn", phe.cov = "ar1")) 
{
    snp.mat <- get_snp_mat(g.snp.hap1, g.snp.pos2, n.minsect, 
        n.maxsect, n.sample, 1)
    f.simu <- f.mn.ar1
    if (covariance$phe.dist == "mn" && covariance$phe.cov == 
        "ar1") 
        f.simu <- f.mn.ar1
    if (covariance$phe.dist == "mn" && covariance$phe.cov == 
        "sad") 
        f.simu <- f.mn.sad
    if (covariance$phe.dist == "mn" && covariance$phe.cov == 
        "cm") 
        f.simu <- f.mn.cm
    if (covariance$phe.dist == "mt") 
        f.simu <- f.mt.ar1
    if (covariance$phe.dist == "msn") 
        f.simu <- f.sn.ar1
    if (covariance$phe.dist == "mmn") 
        f.simu <- f.mmn.ar1
    phe <- simu.long.phe.power(n.sample, snp.mat[[1]], par, f.simu)
    return(list(phe = phe, snp = snp.mat[[1]]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
