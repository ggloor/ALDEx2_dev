\name{summary.aldex}
\alias{summary.aldex}
\title{
    Summarize an \code{aldex} Object
}
\description{
    Summarize, median values of the distributions. Optionally report the entire set of Monte Carlo replicates. Optionally sort based on relative effect-size.
}
\usage{
    \method{summary}{aldex}( object, ... )
}
\arguments{
    \item{object}{
        an object of class \code{aldex} produced by the \code{aldex} function
    }
    \item{...}{
        optional arguments described in `Details', below
    }
}
\details{
The following optional arguments control the format of the \code{data.frame} returned by the function.
    \describe{
        \item{digits}{
            the number of digits to round to, defaulting is three
        }
        \item{sort}{
            whether or not to sort based on decreasing effect-size, defaulting to false
        }
        \item{median.only}{
            whether or not to return medians only or multiple quantiles, defaulting to true
        }
    }
}
\value{
    The nested-list structure of the \code{aldex} object is flattened, rounded, and summarized as specified into a \code{data.frame}. This \code{data.frame} is suitable for saving via \code{write.table} as, for example, a CSV to be loaded into a spreadsheet.
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\seealso{
    \code{\link{aldex}},
    \code{\link{plot.aldex}}
}
\examples{
    # See the examples for 'aldex'.
}
