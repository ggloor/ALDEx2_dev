\name{iqlr.features}
\alias{iqlr.features}
\title{
   calculate glm and Kruskal Wallis test statistics
}
\description{
calculate the features that are to be used in the Geometric Mean calculation in clr_function.R}
\usage{
set.mode(input="default", reads, conds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
    Character argument specifying which indicies to return. 'default' returns all. 'zero' returns the nonzero count features per condition. 'iqlr' returns the features whose variance falls within the inter-quantile range of the CLR-transformed data. In cases of malformed or null queries, input defaults to 'default'.
    Additionally, input can be a numeric vector, which contains a set of indicies that represent the features to center the data against. Only for advanced users who know the invariant set of features within their data.
}
  \item{reads}{
	A data frame containing the samples and features per sample.
}
  \item{conds}{
    A vector describing which samples belong to what condition. 
}
}
\details{
    An explicit example for two conditions is shown in the `Examples' below.
}
\value{
    Outputs a list containing indicies per condition.
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\author{
    Jia Rong Wu
}

\seealso{
    \code{\link{aldex.clr}},
    \code{\link{aldex.ttest}},
    \code{\link{aldex.effect}},
    \code{\link{selex}}
}
\examples{
    # x is the output of the \code{x <- clr(data, mc.samples)} function
    # conditions is a description of the data
    # for the selex dataset, conditions <- c(rep("N", 7), rep("S", 7))
    # input can be "default", "iqlr", "zero" or numeric for advanced users
    data(selex)
    conds <- c(rep("NS", 7), rep("S", 7))
    x <- aldex.clr(selex, conds, mc.samples=2, input="default")
}
