\name{longskat_snp_plink}
\alias{longskat_snp_plink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
longskat_snp_plink(file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.cov, file.gene.set, snp.range = NULL, options = list(y.cov.count = 2, y.cov.time = 0, g.maxiter = 20, w.common = c(0.5, 0.5), w.rare = c(1, 25), run.cpp = T, debug = F, n.cpu = 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.plink.bed}{
%%     ~~Describe \code{file.plink.bed} here~~
}
  \item{file.plink.bim}{
%%     ~~Describe \code{file.plink.bim} here~~
}
  \item{file.plink.fam}{
%%     ~~Describe \code{file.plink.fam} here~~
}
  \item{file.phe.long}{
%%     ~~Describe \code{file.phe.long} here~~
}
  \item{file.phe.cov}{
%%     ~~Describe \code{file.phe.cov} here~~
}
  \item{file.gene.set}{
%%     ~~Describe \code{file.gene.set} here~~
}
  \item{snp.range}{
%%     ~~Describe \code{snp.range} here~~
}
  \item{options}{
%%     ~~Describe \code{options} here~~
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
function (file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, 
    file.phe.cov, file.gene.set, snp.range = NULL, options = list(y.cov.count = 2, 
        y.cov.time = 0, g.maxiter = 20, w.common = c(0.5, 0.5), 
        w.rare = c(1, 25), run.cpp = T, debug = F, n.cpu = 1)) 
{
    cat("[ LONGSKAT_SNP_PLINK ] Procedure\n")
    cat("Checking the optional items......\n")
    cat("* Covariate Count: ", options$y.cov.count, "\n")
    cat("* Covariate Time Effect: ", options$y.cov.time, "\n")
    if (is.null(options$n.cpu) || is.na(options$n.cpu)) 
        options$n.cpu <- 1
    cat("* Parallel Computing: ", ifelse(options$n.cpu > 1, "Yes,", 
        "No,"), options$n.cpu, "CPU(s)\n")
    if (is.null(options$debug) || is.na(options$debug)) 
        options$debug <- FALSE
    cat("* Debug Output: ", ifelse(options$debug, "Yes", "No"), 
        "\n")
    if (is.null(options$run.cpp) || is.na(options$run.cpp)) 
        options$run.cpp <- TRUE
    cat("* C/C++ Module Used Output: ", ifelse(options$run.cpp, 
        "Yes", "No"), "\n")
    if (is.null(options$w.common) || is.na(options$w.common)) 
        options$w.common <- c(0.5, 0.5)
    cat("* Beta Weights for Common SNPs: ", options$w.common[1], 
        options$w.common[2], "\n")
    if (is.null(options$w.rare) || is.na(options$w.rare)) 
        options$w.rare <- c(1, 25)
    cat("* Beta Weights for Rare SNPs: ", options$w.rare[1], 
        options$w.rare[2], "\n")
    chk.plink <- check_plink_file(file.plink.bed, file.plink.bim, 
        file.plink.fam)
    if (!chk.plink$bSuccess) 
        stop("PLINK file can not be loaded by the snpStats package.")
    chk.phe <- check_pheno_file(file.phe.long, chk.plink$family)
    if (!chk.phe$bSuccess) 
        stop("Phenotypic data file is failed to load.")
    chk.cov <- check_covariate_file(file.phe.cov, chk.plink$family, 
        options$y.cov.count)
    if (!chk.cov$bSuccess) 
        stop("Covariate data file is failed to load.")
    chk.genset <- check_geneset_file(file.gene.set)
    if (!chk.genset$bSuccess) 
        stop("Gene defintion file  is failed to load.")
    cat("Starting to load all data files......\n")
    PF <- read_gen_phe_cov(file.plink.bed, file.plink.bim, file.plink.fam, 
        file.phe.long, file.phe.time, file.phe.cov)
    PF.par <- list(file.plink.bed = file.plink.bed, file.plink.bim = file.plink.bim, 
        file.plink.fam = file.plink.fam, file.phe.long = file.phe.long, 
        file.phe.cov = file.phe.cov, file.gene.set = file.gene.set, 
        y.cov.count = options$y.cov.count, y.cov.time = options$y.cov.time, 
        w.common = options$w.common, w.rare = options$w.rare)
    Y.ncol <- NCOL(PF$phe.long)
    Y <- matrix(as.integer(as.matrix(PF$phe.log[, c(2:Y.ncol)])), 
        ncol = Y.ncol - 1)
    Y.cov <- PF$phe.cov[, 3:(2 + options$y.cov.count), drop = F]
    Y.time <- NULL
    if (!is.null(PF$phe.time)) 
        Y.time <- matrix(as.integer(as.matrix(PF$phe.time[, -1])), 
            ncol = Y.ncol - 1)
    cat("Starting to estimate the SIGMA_A, SIGMA_B, SIGMA_E and other parameters......\n")
    r.model <- longskat_est_model(Y, Y.time, Y.cov, y.cov.time = options$y.cov.time, 
        g.maxiter = options$g.maxiter, debug = options$debug)
    if (class(r.model) != "LSKAT.null.model") 
        stop("! Failed to estimate the parameters of Covariance Compoment.")
    cat("* SIGMA_A =", r.model$par$sig_a, "\n")
    cat("* SIGMA_B =", r.model$par$sig_b, "\n")
    cat("* SIGMA_E =", r.model$par$sig_e, "\n")
    cat("* RHO =", r.model$par$rho, "\n")
    cat("* MU =", r.model$par$mu, "\n")
    cat("* Beta(Cov) =", r.model$par$par_cov, "\n")
    cat("* Beta(Time) =", r.model$par$par_t, "\n")
    cat("* L(min) =", r.model$likelihood, "\n")
    snp.len <- dim(PF$snp.mat$genotypes)[2]
    if (is.null(snp.range)) 
        snp.range <- c(1:snp.len)
    if (length(which(snp.range > snp.len)) > 0) {
        warning("The snp range is out of data set.")
        snp.range <- snp.range[-which(snp.range > snp.len)]
    }
    cat("* SNP.RANGE =", min(snp.range), "-", max(snp.range), 
        "[", length(snp.range), "]\n")
    cpu.fun <- function(sect) {
        snp.range0 <- snp.range[((sect - 1) * n.percpu + 1):(sect * 
            n.percpu)]
        PF <- read_gen_phe_cov(PF.par$file.plink.bed, PF.par$file.plink.bim, 
            PF.par$file.plink.fam, PF.par$file.phe.long, PF.par$file.phe.cov)
        ret.cluster <- longskat_snp_task(r.model, snp.range0, 
            PF.par$file.gene.set, PF, weights.rare = options$w.rare, 
            weights.common = options$w.common, run.cpp = options$run.cpp, 
            debug = options$debug)
        return(ret.cluster)
    }
    lskat.ret <- c()
    tm.start <- proc.time()
    if (options$n.cpu > 1 && require(snowfall)) {
        cat("Starting parallel computing, snowfall/snow......\n")
        sfInit(parallel = TRUE, cpus = options$n.cpu, type = "SOCK")
        n.percpu <- ceiling(length(snp.range)/options$n.cpu)
        sfExport("n.percpu", "snp.range", "Y.delt", "Y.t", "X", 
            "par_null", "options", "PF.par")
        ret.cluster <- sfClusterApplyLB(1:options$n.cpu, cpu.fun)
        sfStop()
        cat("Stopping parallel computing......\n")
        lskat.ret <- do.call("rbind", ret.cluster)
    }
    else {
        cat("Starting the LSKAT estimate for each gene......\n")
        lskat.ret <- longskat_snp_task(r.model, snp.range, file.gene.set, 
            PF, weights.rare = options$w.rare, weights.common = options$w.common, 
            run.cpp = options$run.cpp, debug = options$debug)
    }
    tm <- proc.time() - tm.start
    cat("* RUNTIME =", tm[3], "seconds \n")
    ret = list(snp = lskat.ret, par = PF.par, mle = r.model)
    class(ret) <- "LSKAT.snp.plink"
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
