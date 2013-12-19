##requirements:
#reads <- read.table("ALL.data.quad.txt", header=TRUE, row.names=1)
#conditions <- c("NS", "NS","NS","NS","NS","NS","NS","S","S","S","S","S","S","S") 
#
##invocation:
##x <- aldex( reads, conditions, mc.samples=128, lfdr.cutoff=0.1, test="wilcox", fdr="qval")
#
#this version calculates both p, q and lfdr for wilcox and welches
#this version allows exploratory plotting of these values to choose appropriate cutoff values

#requires rdirichlet
"rdirichlet" <-
  function(n, alpha) {
    l <- length(alpha)
    x <- matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE)
    sm <- x%*%rep(1,l)
    return(x/as.vector(sm))
  }
  
#the cutoff, test and fdr method are used for plotting and summary purposes only
#jmd: removed test, fdr, ldfr cutoff; added gene.ratios mode and associated ratio.methods
aldex <- function( reads, conditions, mc.samples=128, gene.ratios=FALSE, ratio.method="") {


    # The 'reads' data.frame MUST have row
    # and column names that are unique, and
    # looks like the following:
    #
    #              T1a T1b  T2  T3  N1  N2
    #   Gene_00001   0   0   2   0   0   1
    #   Gene_00002  20   8  12   5  19  26
    #   Gene_00003   3   0   2   0   0   0
    #       ... many more rows ...
    #
    # The 'conditions' vector must be coercable
    # to a two-level factor that indicates which
    # conditions are to be compared. These conditions
    # will be organized alphabetically for output
    #
    # Note that we do NOT do nesting of technical
    # replicates win biological replicates; instead
    # these should be run separately to asses the
    # relative magnitudes of these effects.

    # ---------------------------------------------------------------------
    # Fully validate and coerce the data into required formats
	
	#####
	#remove all rows with reads less than the minimum
	#set by minsum 
	#####
	minsum <- 0
	
	#put into next version from yohanna
	#nozeroes <- myfunctions[-which(rowSums(myfunctions)>minsum), ]
	
	#JMD: if gene.ratios==TRUE ratio.method must be "columns" or "rows"
	if (gene.ratios==TRUE && !(ratio.method == "columns" | ratio.method == "rows")) {
		stop("When operating in gene.ratios mode, specify ratio.method as \"columns\" or \"rows\".")
		}
		
	#JMD: if gene.ratios==TRUE and ratio.methods=="columns" you need to convert the columns to interleaved rows
	if (gene.ratios ==TRUE && ratio.method=="columns") {
		#this requires an even number of columns. if there isn't, throw error
		if (ncol(reads)%%2 != 0) {
			stop("ratio.method=\"columns\" requires an even number of columns.")
			}
		#odd-numbered columns belong to group A
		aindex <- (c(1:(ncol(reads)/2))*2-1)
		#even ones belong to group B
		bindex <- (c(1:(ncol(reads)/2))*2)
		#make an empty table with twice the rows but half the columns as reads
		temp <-as.data.frame(matrix(data=NA, nrow=(nrow(reads)*2), ncol=(ncol(reads)/2)))
		#give the names of the columns of the new table the names of the columns of group A
		colnames(temp) <- colnames(reads)[aindex]
		i <- 1
		#for each row, make 2 rows in new temp table: interleave A, B, A, B etc. down the columns
		while (i<=nrow(reads)) {
			temp[(2*i-1),] <- reads[i,aindex]
			temp[(2*i),] <- reads[i,bindex]
			rownames(temp)[(2*i-1)] <- paste(rownames(reads)[i], "_A", sep="")
			rownames(temp)[(2*i)] <- paste(rownames(reads)[i], "_B", sep="")
			i <- i+1
			}
		#resulting table has half the columns and twice the rows. rows go: A,B,A,B. each column is a sample
		reads <- temp
		rm(temp)
		}
		
	#JMD: if gene.ratios==TRUE and ratio.methods=="rows", there must be an even number of rows
	if (gene.ratios==TRUE && ratio.method=="rows" && nrow(reads)%%2 != 0) {
		stop("ratio.method=\"rows\" requires an even number of rows.")
		}
	
	z <- as.numeric(apply(reads, 1, sum))
	###JMD: if not in gene ratio mode, remove any rows in which the sum of the reads < minsum
	#this bit of code was there before, I just put it in the if statement
	if (gene.ratios==TRUE) {
	    #if in gene ratio mode, you're going to have to remove any groups of TWO rows in which either one row or the other has a sum of 0
		toolow <- which(z<=minsum)
		others <- c()
		#loop through toolow
		for (i in toolow) {
			#if even add i-1 to others
			if (i %% 2==0) {
				others <- append(others, i-1)
			} else {
				#if odd add i+1 to others
				others <- append(others, i+1)
				}
			}
		#put all these in a vector
		toolow <- unique(append(toolow,others))
		#remove from reads table
		if (length(toolow)>0) {
			reads <- as.data.frame(reads[-toolow,])
			} 
	} else {
		#otherwise just remove any row in which the sum of the row is 0. this was in the aldex code before, I just moved it into the if statement.
		reads <- as.data.frame( reads[(which(z > minsum)),]  )
		}
    rm(z)
    
    if ( any( round(reads) != reads ) ) stop("not all reads are integers")
    if ( any( reads < 0 ) )             stop("one or more reads are negative")

    for ( col in names(reads) ) {
        if ( any( ! is.finite( reads[[col]] ) ) )  stop("one or more reads are not finite")
    }

    if ( length(rownames(reads)) == 0 ) stop("rownames(reads) cannot be empty")
    if ( length(colnames(reads)) == 0 ) stop("colnames(reads) cannot be empty")

    if ( length(rownames(reads)) != length(unique(rownames(reads))) ) stop ("row names are not unique")
    if ( length(colnames(reads)) != length(unique(colnames(reads))) ) stop ("col names are not unique")

    conditions <- as.factor( conditions )
    levels     <- levels( conditions )
    
    if ( length( conditions ) != ncol( reads ) ) stop("mismatch btw 'lenght(conditions)' and 'ncol(reads)'")

    if ( length( levels ) != 2 ) stop("only two condition levels are currently supported")
 
    levels <- vector( "list", length( levels ) )
    names( levels ) <- levels( conditions )
    sets <- names(levels)
    
    #generate the comparison sets from the condition levels
    setA <- which(conditions == sets[1])
    setB <- which(conditions == sets[2])

	if ( length( setA ) < 3 ) stop("require at least 3 replicates in set A")
    if ( length( setB ) < 3 ) stop("require at least 3 replicates in set B")

    for ( l in levels( conditions ) ) {
        levels[[l]] <- which( conditions == l )
        if ( length( levels[[l]] ) <3 ) stop("condition level '",l,"' has less than three replicates")
    }

    if ( mc.samples < 128 ) warning("values are unreliable when estimated with so few MC smps")

	#jmd: commented out since "test" is no longer defined by user
	#if(test != "welches" & test != "wilcox") stop("please choose wilcox or welches test for significance")	

    probs <- c(0.01,0.05,0.50,0.95,0.99)

    # ---------------------------------------------------------------------
    # Monte Carlo spl the frequencies of each spl via the Dirichlet distribution,
    # using the one-group objective reference prior of Berger and Bernardo

    nr <- nrow( reads )
    rn <- rownames( reads )

    p <- lapply( reads , function(col) { q <- t( rdirichlet( mc.samples, col + 0.5 ) ) ; rownames(q) <- rn ; q } )

    for ( i in 1:length(p) ) {
            if ( any( ! is.finite( p[[i]] ) ) ) stop("non-finite frequencies estimated")
    }

    # ---------------------------------------------------------------------
    # Take the log2 of the frequency and remove the noninformative subspace component
    # i.e., do a centered logratio transformation as per Aitchison
	
	#apply the function over elements in a list, that contains an array
    l2p <- lapply( p, function(m) {
        apply( log2(m), 2, function(col) { col - mean(col) } )
    })

    for ( i in 1:length(l2p) ) {
        if ( any( ! is.finite( l2p[[i]] ) ) ) stop("non-finite log-frequencies were unexpectedly computed")
    }

#############
##
## t test subroutine
##
#here is where t-tests are done rather than the effect size stuff below
#but would really like to keep the win-btw calculation values for plotting
#test.df <- data.frame(x1, x2, x3, y1, y2, y3)

#can make a matrix of the number of t and df values to collect data for welches (we) and wilcox (wi)

#JMD: sizes are different if gene.ratios==TRUE
if (gene.ratios==FALSE) {
	we.tfdr.matrix =  matrix(data=NA, nrow=nrow(reads), ncol=mc.samples)
	we.tqval.matrix =  matrix(data=NA, nrow=nrow(reads), ncol=mc.samples)
	we.p.matrix =  matrix(data=NA, nrow=nrow(reads), ncol=mc.samples)

	wi.tfdr.matrix =  matrix(data=NA, nrow=nrow(reads), ncol=mc.samples)
	wi.tqval.matrix =  matrix(data=NA, nrow=nrow(reads), ncol=mc.samples)
	wi.p.matrix =  matrix(data=NA, nrow=nrow(reads), ncol=mc.samples)
} else {
	we.tfdr.matrix =  matrix(data=NA, nrow=nrow(reads)/2, ncol=mc.samples)
	we.tqval.matrix =  matrix(data=NA, nrow=nrow(reads)/2, ncol=mc.samples)
	we.p.matrix =  matrix(data=NA, nrow=nrow(reads)/2, ncol=mc.samples)

	wi.tfdr.matrix =  matrix(data=NA, nrow=nrow(reads)/2, ncol=mc.samples)
	wi.tqval.matrix =  matrix(data=NA, nrow=nrow(reads)/2, ncol=mc.samples)
	wi.p.matrix =  matrix(data=NA, nrow=nrow(reads)/2, ncol=mc.samples)
	}

#mc.i is the monte carlo spl
for(mc.i in 1:mc.samples){
	
	#make the data set for each monte carlo realization
	t.input.data <- matrix(data=NA, nrow=nrow(reads), ncol=(length(l2p)))
	for(i in 1:length(l2p)){ 
		t.input.data[,i] <- l2p[[i]][,mc.i] 
	}
	
	###JMD: if gene.ratios==TRUE calculate them here
	if (gene.ratios==TRUE) {
		#this has half the rows of t input data
		diffs <- matrix(data=NA, nrow=nrow(reads)/2, ncol=(length(l2p)))
		i <- 1
		j <- 1
		while (i<nrow(t.input.data)) {
			#subtract B from A for every pair
			diffs[j,] <- t.input.data[i,]-t.input.data[i+1,]
			i <- i+2
			j <- j+1
			}
		#this is a t.input.data with the differences calculated
		t.input.data <- diffs
		}
	
	#do the t test
	#JMD: originally, wilcox data saved in we.p.matrix/we.f and welch data saved in wi.p.matrix/wi.f. carries to output: super bad, welch data saved in wilcox and vice versa. fixed :D
		x <- t(apply(t.input.data, 1, function(t.input.data){as.numeric(wilcox.test(x=t.input.data[setA],y=t.input.data[setB])[3:4])}))
		#save the values
		wi.p.matrix[,mc.i] <- as.numeric(x[,1])
		wi.f <- (fdrtool(as.numeric(x[,1]), statistic="pvalue", plot=FALSE, color.figure=FALSE, verbose=FALSE))

		x <- t(apply(t.input.data, 1, function(t.input.data){as.numeric(t.test(x=t.input.data[setA],y=t.input.data[setB])[1:3])}))
		#save the values
		we.p.matrix[,mc.i] <- x[,3]
		we.f <- (fdrtool(x[,3], statistic="pvalue", plot=FALSE, color.figure=FALSE, verbose=FALSE))	

	we.tfdr.matrix[,mc.i] <- as.numeric(we.f$lfdr)
	we.tqval.matrix[,mc.i] <- as.numeric(we.f$qval)

	wi.tfdr.matrix[,mc.i] <- as.numeric(wi.f$lfdr)
	wi.tqval.matrix[,mc.i] <- as.numeric(wi.f$qval)
	
	
}

#get the Expected values of p, q and lfdr
we.ep <- apply(we.p.matrix, 1, mean)
we.elfdr <- apply(we.tfdr.matrix,1,mean)
we.eq <- apply(we.tqval.matrix,1,mean)

wi.ep <- apply(wi.p.matrix, 1, mean)
wi.elfdr <- apply(wi.tfdr.matrix,1,mean)
wi.eq <- apply(wi.tqval.matrix,1,mean)


#############
#REMOVE p cuts peak memory
	rm(p)
	gc()
#

    # ---------------------------------------------------------------------
    # Summarize the rab win and all groups
    
    rab <- vector( "list", 3 )
    names(rab) <- c( "all", "win", "spl" )
    rab$win <- list()
    
    #this is the quantile values across all monte carlo replicates
    cl2p <- NULL
    for ( m in l2p ) cl2p <- cbind( cl2p, m )
    rab$all <- t(apply( cl2p, 1, quantile, probs=probs, type=8 ))
    rm(cl2p)
    gc()
    
    #this is the quantile values across all monte carlo replicates per level
    for ( level in levels(conditions) ) {
        cl2p <- NULL
        for ( i in levels[[level]] ) cl2p <- cbind( cl2p, l2p[[i]] )
        rab$win[[level]] <- t(apply( cl2p, 1, quantile, probs=probs, type=8 ))
        rm(cl2p)
        gc()
    }

    rab$spl <- lapply( l2p, function(m) { t(apply( m, 1, quantile, probs=probs, type=8 )) } )

    # ---------------------------------------------------------------------
    # Compute diffs btw and win groups

    l2d <- vector( "list", 2 )
    names( l2d ) <- c( "btw", "win" )
    l2d$win <- list()

    # abs( win-conditions diff ), btw smps

    for ( level in levels(conditions) ) {
        for ( l1 in sort( levels[[level]] ) ) {
            for ( l2 in sort( levels[[level]] ) ) {
                if ( l2 <= l1 ) next
                l2d$win[[level]] <- cbind( l2d$win[[level]] , abs( l2p[[l1]] - l2p[[l2]] ) )
            }
        }
    }

    # Handle the case when the groups have different spl sizes
    # get the minimum number of win spl comparisons
    ncol.wanted <- min( sapply( l2d$win, ncol ) )
    l2d$win  <- lapply( l2d$win, function(arg) { arg[,1:ncol.wanted] } )    

    # btw condition diff (signed)

    for ( l1 in levels[[1]] ) {
        for ( l2 in levels[[2]] ) {
        }}

    for ( l1 in levels[[1]] ) {
        for ( l2 in levels[[2]] ) {
            l2d$btw <- cbind( l2d$btw , ( l2p[[l2]] - l2p[[l1]] ) )
        }
    }
###last use of l2p, remove to reduce memory
	rm(l2p)
	gc()
    win.max <- matrix( 0 , nrow=nr , ncol=ncol.wanted )
    l2d$effect <- matrix( 0 , nrow=nr , ncol=ncol(l2d$btw) )
    rownames(l2d$effect) <- rn
	
###the number of elements in l2d$btw and l2d$win may leave a remainder when
  #recycling these random vectors. Warnings are suppressed because this is not an issue
  #for this calculation. In fact, any attempt to get rid of this error would
  #decrease our power as one or both vectors would need to be truncated gg 20/06/2013
  
	options(warn=-1)
    
    for ( i in 1:nr ) {
        win.max[i,] <- apply( ( rbind( l2d$win[[1]][i,] , l2d$win[[2]][i,] ) ) , 2 , max )
        l2d$effect[i,] <- l2d$btw[i,] / win.max[i,]
    }

	options(warn=0)

    rownames(win.max)   <- rn
    attr(l2d$win,"max") <- win.max
    rm(win.max)

    # ---------------------------------------------------------------------
    # Summarize diffs

    l2s <- vector( "list", 2 )
    names( l2s ) <- c( "btw", "win" )
    l2s$win <- list()

    l2s$btw <- t(apply( l2d$btw, 1, quantile, probs=probs, type=8 ))
    l2s$win  <- t(apply( attr(l2d$win,"max"), 1, quantile, probs=probs, type=8 ))
    
    effect  <- t(apply( l2d$effect, 1, quantile, probs=probs, type=8 ))


    # ---------------------------------------------------------------------
    # Done
    
    #jmd
    if (gene.ratios==TRUE) {
   		rn <- rownames(reads)[(c(1:(nrow(reads)/2))*2-1)]
   		if (ratio.method=="columns") {
   			rn <- substr(rn,1,nchar(rn)-2)
   			}
		}
    rv <- list(
        rab = rab,
        diff = l2s,
        effect = effect,
        criteria = data.frame(
            we.pval = we.ep,
            we.qval = we.eq,
            we.lfdr = we.elfdr,

            wi.pval = wi.ep,
            wi.qval = wi.eq,
            wi.lfdr = wi.elfdr,

        	row.names    = rn
        )
    )
    attr(rv,"probs") <- probs
    class(rv) <- c( "aldex", class(rv) )
    on.exit( expr=gc(), add=TRUE )
    return(rv)

}


plot.aldex <- function( x, ..., type=c("MW","MA","hist"),
    xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
    nbin=256, nrpoints=Inf, col=rgb(0,0,0,0.4), pch='.',
    signif.col="cyan", signif.pch=20, signif.cex=0.4,
    called.col="red", called.pch=20, called.cex=0.5,
    thres.line.col="darkgrey", thres.lwd=1.5, fdr = "qval",
    test="welches", cutoff=0.1
) {

    stopifnot( inherits( x, "aldex" ) )

    type <- match.arg(type)
    if (test == "welches"){
    	if ( fdr == "lfdr" ) called <- x$criteria$we.lfdr <= cutoff
   		if ( fdr == "qval" ) called <- x$criteria$we.qval <= cutoff
   	}

    if (test == "wilcox"){
    	if ( fdr == "lfdr" ) called <- x$criteria$wi.lfdr <= cutoff
   		if ( fdr == "qval" ) called <- x$criteria$wi.qval <= cutoff
   	}
    signif <- x$criteria$sig.p

    if ( is.null(col) ) col <- blues9[9]

    if ( type == "MA" ) {

        if ( is.null(xlab) ) xlab <- expression( "Median" ~~ Log[2] ~~ "Combined- Level" )
        if ( is.null(ylab) ) ylab <- expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" )

        if ( is.null(xlim) ) xlim <- range( x$rab$all[      ,"50%"] )
        if ( is.null(ylim) ) ylim <- range( x$diff$btw [      ,"50%"] )

        plot( x$rab$all[      ,"50%"], x$diff$btw[      ,"50%"], xlim=xlim, ylim=ylim, col=col, pch=pch, xlab=xlab, ylab=ylab, main=paste(c(test, fdr), sep=", ") )
        points       ( x$rab$all[signif,"50%"], x$diff$btw[signif,"50%"], pch=signif.pch, cex=signif.cex, col=signif.col )
        points       ( x$rab$all[called,"50%"], x$diff$btw[called,"50%"], pch=called.pch, cex=called.cex, col=called.col )

    } else if ( type == "MW" ) {

        if ( is.null(xlab) ) xlab <- expression( "Median" ~~ Log[2] ~~ "win-Condition diff" )
        if ( is.null(ylab) ) ylab <- expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" )

        if ( is.null(xlim) ) xlim <- range( x$diff$win [      ,"50%"] )
        if ( is.null(ylim) ) ylim <- range( x$diff$btw[      ,"50%"] )

        plot( x$diff$win[      ,"50%"], x$diff$btw[      ,"50%"], xlim=xlim, ylim=ylim, col=col, pch=pch, xlab=xlab, ylab=ylab, main=paste(c(test, fdr), sep=", ") )
        points(        x$diff$win[signif,"50%"], x$diff$btw[signif,"50%"], pch=signif.pch, cex=signif.cex, col=signif.col )
        points(        x$diff$win[called,"50%"], x$diff$btw[called,"50%"], pch=called.pch, cex=called.cex, col=called.col )
        abline( a=0, b= 1, col=thres.line.col, lwd=thres.lwd, lty=2 )
        abline( a=0, b=-1, col=thres.line.col, lwd=thres.lwd, lty=2 )

    } else if ( type == "hist" ) {
    	par(mfrow=c(1,2))
    	hist(x$criteria$we.pval, breaks=100, main="Welches t-test p values")
    	hist(x$criteria$wi.pval, breaks=100, main="Wilcox rank test p values")
    } else stop("unkown plot type")

    invisible()
}


summary.aldex <- function( object, ... ) {
    
    stopifnot( inherits( object, "aldex" ) )
    x <- object
    rm(object)

    digits <- 3
    sort   <- FALSE
    median.only <- TRUE
    minimal.only <- TRUE
    
    a <- as.list( match.call() )
    if ( ! is.null( a[["digits"]] ))        digits <- a$digits
    if ( ! is.null( a[["sort"]] ))          sort <- a$sort
    if ( ! is.null( a[["median.only"]] ))   median.only <- a$median.only
    if ( ! is.null( a[["minimal.only"]] ))   minimal.only <- a$minimal.only
    
    y <- data.frame(x)
    yn <- names(y)
    
    for ( n in yn ) if ( is.numeric( y[[n]] ) ) y[[n]] <- round( y[[n]], digits=digits )
    
    yn <- gsub( "\\.([0-9])\\.$", ".q0\\1", yn )
    yn <- gsub( "\\.([0-9][0-9])\\.$", ".q\\1", yn )
    names(y) <- yn

    if ( median.only ) {
        is.quantile <- grepl( "\\.q([0-9][0-9])$", yn )
        not.median  <- ! grepl( "\\.q50$", yn )
        wanted <- ( ! is.quantile ) | ( ! not.median )
        y <- y[ , wanted ]
    }
	
	if ( minimal.only ) {
		keep <- (! grepl("rab.spl", names(y)))
		y <- y[ , keep ]
	}
	
    if ( sort ) y <- y[ order( y$effect.q50 , decreasing=TRUE ) , ] 
    
    return(y)    
}

