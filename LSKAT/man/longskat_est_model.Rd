\name{longskat_est_model}
\alias{longskat_est_model}
\title{
NULL model estimation
}
\description{
Estimating the parameters and residuals for the NULL model in LSKAT.
}
\usage{
longskat_est_model( phe.long,
    phe.cov, 
    phe.time  = NULL, 
    time.cov  = 0, 
    intercept = FALSE, 
    method    = c("REML", "ML"),
    g.maxiter = 20, 
    par.init  = list(), 
    verbose   = F)
}
\arguments{
  \item{phe.long}{Phenotype matrix with \var{m} rows denoting the individuals and \var{n} columns denoting the time points, the row name indicates the individual's ID. }
  \item{phe.cov}{Time covariate matrix with \var{m} rows denoting individuals and \var{x} columns denoting the covariate variables, the row name indicates the individual's ID. }
  \item{phe.time}{Time point matrix with \var{m} rows denoting individuals and \var{n} columns denoting the time points, the row name indicates the individual's ID. If this matrix is not specified, the default matrix is generated.}
  \item{time.cov}{Numeric, indicating whether the time exponnents are included as extra covariates, The time points are used if 1, the time points and time squares are used if 2, and so on. The default value (0) doesn't use the time covariate.}
  \item{intercept}{Logical variable, indicating whether the intercept is estimated. }
  \item{method}{String, REML or ML are available for the parameter estimation. }
  \item{g.maxiter}{Numeric, the maximum count for the iterative estimation. }
  \item{par.init}{List, the initial values for the parameter \code{rho}, \code{sig.a}, \code{sig.b}, \code{sig.e}. }
  \item{verbose}{Logical variable, indicating whether some debug information can be outputted. }
}
\details{
}

\value{
This function returns an list object with model parameters and residuals of the NULL model which assumes there is no association between genes and longitudinal phenotypes.\cr

The return object is a list with the following items:\cr
  \item{par}{List, model paramters as shown in below.}
  \item{likelihood}{Numeric, the likelihood value estimated by REML or ML.}
  \item{phe.delt}{ Residual matrix with the row name indicating the individual's ID, the structure is same as \code{phe.long}. }
  \item{phe.time}{ Time matrix, copied from the input parameter \code{phe.time}.}
  \item{phe.cov}{ Covariate matrix, copied from the input parameter \code{phe.cov}.}

The Model paramters: \code{par} has the following sub-items:\cr
  \item{intercept}{Logical variable copied from the input parameter, indicating the intercept is estimated.}
  \item{mu}{Numeric indicating the intercept value.}
  \item{cov.effect}{String indicating the coefficient of the covariates except intercept}
  \item{sig.a}{String indicating the standar deviation of individual random effects}
  \item{sig.b}{String indicating the standar deviation of individual-specific timede-pendent random effects.}
  \item{sig.e}{String indicating the standar deviation of measurement error.}
  \item{rho}{String indicating the corelation coefficient of covariance structure.}
  \item{time.cov}{Numeric, indicating whether consider times as covariate, 0 means no time effects, 1 means time effects, 2 means time effects and time square effects are included as covariates. and so on.}
  \item{time.effect}{Vector of numeric, the time coefficient of time effects. The 1st item is the coefficient for time effects, The 2nd item is the coefficient for time square effects and so on.}

After obtaining the model parameters, please use the \code{\link{longskat_gene_test}} to test the association between gene and traits.
}

\references{
Wang Z., Xu K., Zhang X., Wu X., and Wang Z., (2016) Longitudinal SNP-set association analysis of quantitative phenotypes. Genetic Epidemiology.
}
\seealso{
}
\examples{

## Data simulation using the default parameters
p0 <- longskat_gene_simulate();

## Estimating the model parameters and residuals
r.model0 <- longskat_est_model( p0$phe.long, p0$phe.cov, g.maxiter=3, verbose=T);

##print this model
print(r.model0);

}
\keyword{ NULL Model }
