\name{longskat_gene_plink}
\alias{longskat_gene_plink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
longskat_gene_plink(file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.time = NULL, file.phe.cov, file.gene.set, gene.range = NULL, options = list(y.cov.count = 2, y.cov.time = 0, g.maxiter = 20, w.common = c(0.5, 0.5), w.rare = c(1, 25), run.cpp = T, debug = F, n.cpu = 1))
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
  \item{file.phe.time}{
%%     ~~Describe \code{file.phe.time} here~~
}
  \item{file.phe.cov}{
%%     ~~Describe \code{file.phe.cov} here~~
}
  \item{file.gene.set}{
%%     ~~Describe \code{file.gene.set} here~~
}
  \item{gene.range}{
%%     ~~Describe \code{gene.range} here~~
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
    file.phe.time = NULL, file.phe.cov, file.gene.set, gene.range = NULL, 
    options = list(y.cov.count = 2, y.cov.time = 0, g.maxiter = 20, 
        w.common = c(0.5, 0.5), w.rare = c(1, 25), run.cpp = T, 
        debug = F, n.cpu = 1)) 
{
    cat("[ LONGSKAT_GENE_PLINK ] Procedure.\n")
    cat("Checking the optional items......\n")
    if (missing(options)) 
        options <- get_default_options()
    else {
        options0 <- get_default_options()
        options0[names(options)] <- options
        options <- options0
    }
    cat("* Covariate Count: ", options$y.cov.count, "\n")
    cat("* Covariate Time Effect: ", options$y.cov.time, "\n")
    cat("* Parallel Computing: ", ifelse(options$n.cpu > 1, "Yes,", 
        "No,"), options$n.cpu, "CPU(s)\n")
    cat("* Debug Output: ", ifelse(options$debug, "Yes", "No"), 
        "\n")
    cat("* C/C++ Module Used Output: ", ifelse(options$run.cpp, 
        "Yes", "No"), "\n")
    cat("* Beta Weights for Common SNPs: ", options$w.common[1], 
        options$w.common[2], "\n")
    cat("* Beta Weights for Rare SNPs: ", options$w.rare[1], 
        options$w.rare[2], "\n")
    chk.plink <- check_plink_file(file.plink.bed, file.plink.bim, 
        file.plink.fam)
    if (!chk.plink$bSuccess) 
        stop("PLINK file can not be loaded by the snpStats package.")
    chk.phe <- check_pheno_file(file.phe.long, file.phe.time, 
        chk.plink$family)
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
        file.phe.time = file.phe.time, file.phe.cov = file.phe.cov, 
        file.gene.set = file.gene.set, y.cov.count = options$y.cov.count, 
        y.cov.time = options$y.cov.time, w.common = options$w.common, 
        w.rare = options$w.rare)
    Y.ncol <- NCOL(PF$phe.long)
    Y <- matrix(as.integer(as.matrix(PF$phe.long[, -1])), ncol = Y.ncol - 
        1)
    Y.cov <- PF$phe.cov[, c(3:(2 + options$y.cov.count))]
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
    cat("* L(min) =", r.model$par$likelihood, "\n")
    gen.list <- read_gen_dataset(file.gene.set)
    if (is.null(gene.range)) 
        gene.range <- c(1:gen.list$len)
    if (length(which(gene.range > gen.list$len)) > 0) {
        warning("The gene range is out of data set.")
        gene.range <- gene.range[-which(gene.range > gen.list$len)]
    }
    cat("* GENE.RANGE =", min(gene.range), "-", max(gene.range), 
        "[", length(gene.range), "]\n")
    tm.start <- proc.time()
    cpu.fun <- function(sect) {
        g.range0 <- gene.range[((sect - 1) * n.percpu + 1):(sect * 
            n.percpu)]
        PF <- read_gen_phe_cov(PF.par$file.plink.bed, PF.par$file.plink.bim, 
            PF.par$file.plink.fam, PF.par$file.phe.long, PF.par$file.phe.time, 
            PF.par$file.phe.cov)
        res.cluster <- longskat_gene_task(r.model, g.range0, 
            PF, gen.list, options$w.common, options$w.rare, options$run.cpp, 
            options$debug)
        return(res.cluster)
    }
    lskat.ret <- c()
    if (options$n.cpu > 1 && require(snowfall)) {
        cat("Starting parallel computing, snowfall/snow......\n")
        sfInit(parallel = TRUE, cpus = options$n.cpu, type = "SOCK")
        n.percpu <- ceiling(length(gene.range)/options$n.cpu)
        sfExport("n.percpu", "gene.range", "r.model", "gen.list", 
            "options", "PF.par")
        res.cluster <- sfClusterApplyLB(1:options$n.cpu, cpu.fun)
        sfStop()
        cat("Stopping parallel computing......\n")
        lskat.ret <- do.call("rbind", res.cluster)
    }
    else {
        cat("Starting the LSKAT estimate for each gene......\n")
        lskat.ret <- longskat_gene_task(r.model, gene.range, 
            PF, gen.list, options$w.common, options$w.rare, options$run.cpp, 
            options$debug)
    }
    tm <- proc.time() - tm.start
    cat("* RUNTIME =", tm[3], "seconds \n")
    ret = list(gene = lskat.ret, par = PF.par, mle = r.model)
    class(ret) <- "LSKAT.gen.plink"
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
