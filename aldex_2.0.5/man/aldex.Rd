\name{aldex}
\alias{aldex}
\title{
    Compute an \code{aldex} Object
}
\description{
    Compute an ALDEx Monte Carlo analysis for the given read-counts and experimental conditions. Statistical significance by Wilcox test or Welch's t-test
}
\usage{
    aldex( reads, conditions,
           mc.samples=128, gene.ratios=FALSE, ratio.method="" )
}
\arguments{
    \item{reads}{
        an non-negative, integer-only containing \code{data.frame} that has unique names for all rows and columns, where each row is a different gene and each column represents a sequencing read-count. Rows with 0 
        reads in each sample will be deleted prior to analysis
    }
    \item{conditions}{
        a \code{vector} of length \code{ncol(reads)} that can be coerced into a \code{factor} representing the conditions of the samples; for example: `treatment' and `control', or `normal' and `disease'
    }
    \item{mc.samples}{
        the number of Monte Carlo samples to use to estimate the underlying distributions; since we are estimating central tendencies, 128 is usually sufficient
    }
	\item{gene.ratios}{
		If \code{gene.ratios=TRUE}, ALDEx will compare the expression ratio of pairs of genes between conditions, instead of expression of single genes. If \code{gene.ratios=TRUE}, \code{ratio.method} must be defined. Examples of uses for this mode include comparing sense:antisense transcription ratios for multiple genes between conditions, or comparing the transcription in riboswitches to their downstream genes and observing if that changes between conditions.
	}
	\item{ratio.method}{
		Define if \code{gene.ratios=TRUE}. If each column is a sample and gene pairs are in adjacent rows, choose \code{ratio.method="rows"}. If each sample is represented by two columns and gene pairs are in adjacent columns, choose \code{ratio.method="columns"}.
	}
}
\details{

An explicit description of the input format for the \code{reads} object is shown under `Examples', below.

Note that \code{lfdr.cutoff} is provided for convenience only; the Boolean column in the returned list is computed via these values and they can be adjusted without needing to re-run the analysis.

}
\value{
    An object of class \code{aldex} produced by the \code{aldex} function is a nested list-of-lists structure that is easy to visualize using the \code{str(x)} function. The return value is a list with the following elements:
    \item{expression}{
        a list with three elements, \code{among}, \code{within}, and \code{sample}, each of which gives the adjusted-log2 expression values among-conditions, within-conditions, and for each sample, respectively
    }
    \item{difference}{
        a list with the two elements \code{between} and \code{within} giving the expression-differences between- and within-conditions, respectively
    }
    \item{effect}{
        quantiles of the inferred effect-size calculated by scaling the difference between groups by the largest within-group difference
    }
    \item{criteria}{
        a list of elements giving the effect-size, p value, q value, and local false discovery rate for each gene
    }
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\seealso{
    \code{\link{plot.aldex}},
    \code{\link{summary.aldex}}
}
\examples{


    # The 'reads' data.frame MUST have row
    # and column names that are unique, and
    # looks like the following:
    #
    #              T1a T1b  T2  T3  N1  N2  Nx
    #   Gene_00001   0   0   2   0   0   1   0
    #   Gene_00002  20   8  12   5  19  26  14
    #   Gene_00003   3   0   2   0   0   0   1
    #   Gene_00004  75  84 241 149 271 257 188
    #   Gene_00005  10  16   4   0   4  10  10
    #   Gene_00006 129 126 451 223 243 149 209
    #       ... many more rows ...
    #
    #
    # The 'conditions' vector must be coercable
    # to a two-level factor that indicates which
    # conditions are to be compared.
    #
    # Note that we do NOT do nesting of technical
    # replicates within biological replicates; instead
    # these should be run separately to asses the
    # relative magnitudes of these effects.

    # reads <- read.delim( file="file_as_above.tsv" )
    # conditions <- c( "T", "T", "T", "T", "N", "N", "N" )
    # x   <- x <- aldex( reads, conditions, mc.samples=128, lfdr.cutoff=0.1, test="wilcox" )
    # 
    # par( pty='s', mfrow=c(1,2) )
    # plot( x, type="MA" )
    # plot( x, type="MW" )
    # 
    # sx  <- summary( x )
    # sxm <- summary( x, sort=TRUE, median.only=TRUE )

}
