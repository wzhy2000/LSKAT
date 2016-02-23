\name{plot.LSKAT.gen.plink}
\alias{plot.LSKAT.gen.plink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.LSKAT.gen.plink(r.lskat, pdf.file = NA, title = "", y.max = NA, bonferroni = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r.lskat}{
%%     ~~Describe \code{r.lskat} here~~
}
  \item{pdf.file}{
%%     ~~Describe \code{pdf.file} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{y.max}{
%%     ~~Describe \code{y.max} here~~
}
  \item{bonferroni}{
%%     ~~Describe \code{bonferroni} here~~
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
function (r.lskat, pdf.file = NA, title = "", y.max = NA, bonferroni = F) 
{
    if (is.na(pdf.file)) 
        pdf.file <- paste(r.lskat$par$file.phe.long, ".pdf", 
            sep = "")
    df.gene <- r.lskat$gene[, c(3, 4, 2, 8)]
    df.gene <- df.gene[order(df.gene[, 1], df.gene[, 2]), ]
    if (bonferroni) 
        df.gene[, 4] <- df.gene[, 4] * NROW(df.gene)
    pdf(pdf.file, width = 6, height = 3.5)
    par(mar = c(3.2, 3.5, 2.5, 1))
    par(mgp = c(1.4, 0.3, 0), tck = -0.03)
    draw_manhattan(df.gene[, c(1, 3, 4)], map.title = title, 
        ifelse(bonferroni, 1e-04, 0.05)/NROW(df.gene), 0.7, y.max = y.max)
    dev.off()
    cat("* Manhattan plot is save into ", pdf.file, "\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
