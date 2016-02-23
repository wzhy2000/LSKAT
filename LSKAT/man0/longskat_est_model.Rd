\name{longskat_est_model}
\alias{longskat_est_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
longskat_est_model(y.log, y.time = NULL, y.cov, y.cov.nTime = 0, g.maxiter = 20, debug = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.log}{
%%     ~~Describe \code{y.log} here~~
}
  \item{y.time}{
%%     ~~Describe \code{y.time} here~~
}
  \item{y.cov}{
%%     ~~Describe \code{y.cov} here~~
}
  \item{y.cov.nTime}{
%%     ~~Describe \code{y.cov.nTime} here~~
}
  \item{g.maxiter}{
%%     ~~Describe \code{g.maxiter} here~~
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
function (y.log, y.time = NULL, y.cov, y.cov.nTime = 0, g.maxiter = 20, 
    debug = F) 
{
    nrow <- dim(y.log)[1]
    ncol <- dim(y.log)[2]
    nCov <- dim(y.cov)[2]
    get_par <- function(par, y, y.time, y.cov) {
        par_rho <- par[1]
        if (par_rho < 0 || par_rho >= 0.99) 
            return(NaN)
        sig_a <- par[2]
        sig_b <- par[3]
        sig_e <- par[4]
        par_u <- par[5]
        par_cov <- par[c(6:(5 + nCov))]
        AR.1 <- array(0, dim = c(ncol, ncol))
        for (i in 1:ncol) for (j in 1:ncol) AR.1[i, j] <- par_rho^abs(i - 
            j)
        sigma <- diag(sig_e^2, ncol) + sig_a^2 + sig_b^2 * AR.1
        y.delt <- y - par_u
        for (i in 1:nCov) y.delt <- y.delt - y.cov[, i] * par_cov[i]
        par_t <- c()
        if (y.cov.nTime > 0) {
            par_t <- par[5 + nCov + c(1:y.cov.nTime)]
            for (i in 1:y.cov.nTime) y.delt <- y.delt - (y.time^i) * 
                par_t[i]
        }
        A <- 0
        for (i in 1:NROW(y.delt)) {
            t.sel <- which(!is.na(c(y.time[i, ])))
            sig <- sigma[t.sel, t.sel, drop = F]
            A <- A - sum(dmvnorm(y.delt[i, t.sel, drop = F], 
                rep(0, length(t.sel)), sig, log = T))
        }
        return(A)
    }
    par.init <- c(0.5, sd(y.log, na.rm = T), sd(y.log, na.rm = T), 
        sd(y.log, na.rm = T), mean(y.log, na.rm = T))
    for (i in 1:nCov) par.init <- c(par.init, 1/(mean(y.cov[, 
        i], na.rm = T)^2 + 1))
    if (y.cov.nTime > 0) {
        par.init <- c(par.init, rep(1/(mean(y.time, na.rm = T)^2 + 
            1), y.cov.nTime))
        if (is.null(y.time)) 
            y.time <- t(t(ifelse(is.na(y.log), NA, 1)) * (rep(1:NROW(y.log))))
    }
    loop.n <- 0
    min.val <- Inf
    min.par <- par.init
    while (loop.n <= g.maxiter) {
        r0 <- try(optim(par.init, get_par, y = y.log, y.time = y.time, 
            y.cov = y.cov, method = "BFGS", control = list(maxit = 500)), 
            silent = F)
        if (class(r0) == "try-error") {
            if (min.val >= 1e+08) 
                loop.n <- loop.n + 0.2
            par.init <- min.par * runif(length(min.par))
            next
        }
        loop.n <- loop.n + 1
        if (r0$convergence == 0 && r0$val < min.val) {
            if (debug) 
                cat("  LOOP =", loop.n, "/", g.maxiter, " val=", 
                  r0$val, "par=", r0$par, "\n")
            min.val <- r0$value
            min.par <- r0$par
            par.init <- min.par
        }
        par.init <- par.init * runif(length(par.init), 0.8, 1.2)
    }
    if (debug) 
        cat("  LOOP =", loop.n, " min.val=", min.val, "par=", 
            min.par, "\n")
    if (is.infinite(min.val) || is.na(min.val) || any(is.na(min.par))) 
        return(list(bSuccess = F))
    par_cov <- min.par[c(6:(5 + nCov))]
    y.delt <- y.log - min.par[5]
    for (i in 1:nCov) y.delt <- y.delt - y.cov[, i] * par_cov[i]
    par_t <- c()
    if (y.cov.nTime > 0) {
        par_t <- min.par[5 + nCov + c(1:y.cov.nTime)]
        for (i in 1:length(par_t)) y.delt <- y.delt - (y.time^i) * 
            par_t[i]
    }
    pars <- list(mu = min.par[5], rho = min.par[1], sig_a = abs(min.par[2]), 
        sig_b = abs(min.par[3]), sig_e = abs(min.par[4]), par_cov = par_cov, 
        par_t = par_t)
    r.model <- list(par = pars, likelihood = min.val, y.cov.nTime = y.cov.nTime, 
        y.delt = y.delt, y.time = y.time, y.cov = y.cov)
    class(r.model) <- "LSKAT.null.model"
    return(r.model)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
