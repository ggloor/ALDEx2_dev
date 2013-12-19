\name{plot.aldex}
\alias{plot.aldex}
\title{
    Plot an \code{aldex} Object
}
\description{
    Create `MW'- or `MA'-type plots from the given \code{aldex} object.
}
\usage{
    \method{plot}{aldex}( x, ..., type=c("MW","MA","hist"),
    xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
    nbin=256, nrpoints=Inf, col=rgb(0,0,0,0.4), pch='.',
    signif.col="cyan", signif.pch=20, signif.cex=0.4,
    called.col="red", called.pch=20, called.cex=0.5,
    thres.line.col="darkgrey", thres.lwd=1.5, fdr="qval", 
    test="welches", cutoff=0.1)
}
\arguments{
    \item{x}{
        an object of class \code{aldex} produced by the \code{aldex} function
    }
    \item{...}{
        optional, unused arguments included for compatibility withe the S3 method signature
    }
    \item{type}{
        which type of plot is to be produced. MA is a Bland-Altman style plot; MW is a difference between to a variance within plot as described in the paper; hist is a histogram of p values
    }
    \item{xlab}{
        the x-label for the plot, as per the parent \code{plot} function
    }
    \item{ylab}{
        the y-label for the plot, as per the parent \code{plot} function
    }
    \item{xlim}{
        the x-limits for the plot, as per the parent \code{plot} function
    }
    \item{ylim}{
        the y-limits for the plot, as per the parent \code{plot} function
    }
    \item{nbin}{
        the number of grid points desired, as per the \code{smoothScatter} function
    }
    \item{nrpoints}{
        the number of points to be superimposed on the density, as per the \code{smoothScatter} function
    }
    \item{col}{
        the default colour of the plotted points, as per the \code{smoothScatter} function
    }
    \item{pch}{
        the default plotting symbol, as per the \code{smoothScatter} function
    }
    \item{signif.col}{
        the default colour of `significant' points with p < 0.05
    }
    \item{signif.pch}{
        the default symbol of `significant' points
    }
    \item{signif.cex}{
        the character expansion of `significant' points
    }
    \item{called.col}{
        the colour of points with local false discovery rate <= 0.1
    }
    \item{called.pch}{
        the symbol of points with local false discovery rate <= 0.1
    }
    \item{called.cex}{
        the character expansion of points with local false discovery rate <= 0.1
    }
    \item{thres.line.col}{
        the colour of the threshold line where within and between group variation is equivalent
    }
    \item{thres.lwd}{
        the width of the threshold line where within and between group variation is equivalent
    }
    \item{test}{
        the method of calculating p values, one of welches or wilcox
    }
    \item{fdr}{
        the method of characterizing fdr, one of qval or lfdr 
    }
    \item{cutoff}{
        the fdr cutoff, default 0.1 
    }
}
\details{
    This particular specialization of the \code{plot} function is relatively simple and provided for convenience. For more advanced control of the plot is is best to use the values returned by \code{summary(x)}.
}
\value{
    None.
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\seealso{
    \code{\link{aldex}},
    \code{\link{summary.aldex}},
    \code{smoothScatter}
}
\examples{
    # See the examples for 'aldex'.
}
