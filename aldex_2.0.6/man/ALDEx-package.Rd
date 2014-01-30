\name{ALDEx-package}
\alias{ALDEx-package}
\alias{ALDEx2}
\docType{package}
\title{
    ANOVA-Like Differential Gene Expression Analysis of Single-Organism and Meta-RNA-Seq
}
\description{
    A differential count abundance analysis for the comparison of two conditions. Uses a Dirichlet-multinomial model to infer abundance from counts, that has been optimized for three or more experimental replicates. Infers sampling variation and calculates the expected local false discovery rate given the biological and sampling variation. User can choose either Welches t-test or Wilcox rank test for significance with Benjamini-Hochberg correction for multiple hypothesis testing. NOTE: All components with 0 reads across all samples are removed prior to analysis
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\keyword{ package }
\seealso{
    \code{\link{aldex}},
    \code{\link{plot.aldex}},
    \code{\link{summary.aldex}}
}
\examples{
    # See the examples for the 'aldex', 'plot.aldex',
    # and 'summary.aldex' functions.
}