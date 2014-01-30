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
    nbin=256, nrpoints=Inf, col=rgb(0.3,0.3,0.3,0.4), pch='.',
    signif.col="cyan", signif.pch=20, signif.cex=0.5,
    called.col="red", called.pch=20, called.cex=0.5,
    thres.line.col=rgb(0.3,0.3,0.3,0.5), thres.lwd=1.5, fdr = "BH",
    test="welches", cutoff=0.05, rare.col=black, rare=0,
    rare.pch=20, rare.cex=0.1 )
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
    \item{test}{
        the method of calculating p values, one of welches t-test or wilcoxon
    }
    \item{fdr}{
        the method of characterizing fdr, currently only BH for Benjamini-Hochberg
    }
    \item{cutoff}{
        the Benjamini-Hochberg fdr cutoff, default 0.05 
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
        the colour of points with false discovery rate <= 0.05
    }
    \item{called.pch}{
        the symbol of points with false discovery rate <= 0.05
    }
    \item{called.cex}{
        the character expansion of points with false discovery rate <= 0.05
    }
    \item{thres.line.col}{
        the colour of the threshold line where within and between group variation is equivalent
    }
    \item{thres.lwd}{
        the width of the threshold line where within and between group variation is equivalent
    }
    \item{rare.col}{
        color for features less abundant than the mean, default darkblue
    }
    \item{rare}{
        cutoff for rare features, default 0 or the mean abundance 
    }
    \item{rare.pch}{
        the default symbol of rare points
    }
    \item{rare.cex}{
        the character expansion of rare points
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
