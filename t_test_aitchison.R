aitchison_ttest <- function(in.table, conds){

	
	# clr transform and report t-test output for a data table
	# samples in columns, features in rows
	# conds is a list of conditions for the table

	if(min(in.table) == 0) stop("nay, nay, nay! replace 0 values with a prior")

	#clr transform
	#table.clr <- apply(selex, 2, function(x){log2(x) - mean(log2(x))})
	table.clr <- apply(in.table, 2, function(x){log2(x) - mean(log2(x))})
	
	conditions <- as.factor( conds )
 	lev <- levels(conditions)
 	
 	# t-test and standard deviation functions for the apply
	#compT <- function(x, grp){as.numeric(t.test(x[conditions == lev[1]], x[conditions == lev[2]])[3])}
	compT <- function(x, grp){t.test(x[conditions == lev[1]], x[conditions == lev[2]])$p.value}
	#compW <- function(x, grp){as.numeric(wilcox.test(x[conditions == lev[1]], x[conditions == lev[2]])[3])}

	calc.sd.max <- function(x, grp){
		max(
			c(sqrt( sum( (x[conditions == lev[1]] - mean(x[conditions == lev[1]]))^2)/length(x[conditions == lev[1]])),
			  sqrt( sum( (x[conditions == lev[2]] - mean(x[conditions == lev[2]]))^2)/length(x[conditions == lev[2]]))
			) 
		)
	}
		
	
#	var.total <- apply(table.clr, 1, function(x){sqrt( sum((x - mean(x))^2)/length(x) )})
	sd.max <- apply(table.clr, 1, calc.sd.max, conditions)
	
	dif <- apply(table.clr, 1, function(x){mean(x[conditions == lev[1]]) - mean(x[conditions == lev[2]])})

	# test functions using apply
	t <- apply(table.clr, 1, compT, conditions)
	t.bh <- p.adjust(t)
	#w <- apply(table.clr, 1, compW, conditions)

	out <- list(
		t, 
		t.bh,
		sd.max,
		dif
	)
	names(out) <- c("p.value", "p.adjust", "max_sd", "difference")
	
#x <- aitchison_ttest(selex,conds)	
	plot(sd.max, dif, pch=19, cex=0.3, col=rgb(0,0,0,0.2))
	points(sd.max[t.bh<0.5], dif[t.bh<0.5], pch=19, cex=0.3, col=rgb(1,0,0,0.5))
	abline(0,2, lty=2, col=rgb(0,0,0,0.3))
	abline(0,-2, lty=2, col=rgb(0,0,0,0.3))

	return(out)
}