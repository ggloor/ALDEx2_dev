

X <- rnorm(1000,0)
Y <- rnorm(1000,1)

diff.win <- function(x,y){
	z <- sample(x, size=length(x), replace=FALSE)
	b <- sample(y, size=length(y), replace=FALSE)
	c <- cbind(abs(z-x), abs(b-y))
	a <- apply(c,1,max)
	return(summary(a))
}

diff.btw <- function(x,y){
	z <- x-y
	return(summary(z))
}

summary.table <- matrix(data=NA, nrow=100, ncol=6)

for(i in 1:100){
	X <- rnorm(10000,0, sd=1)
	Y <- rnorm(10000,10, sd=5)
	summary.table[i,1] <- mean(Y) - mean(X)
	summary.table[i,3] <- max(c(sd(X),sd(Y)))
	summary.table[i,2] <- diff.btw(Y,X)[3]
	summary.table[i,4] <- diff.win(X,Y)[3]
	summary.table[i,5] <- summary.table[i,1]/summary.table[i,3]
	summary.table[i,6] <- summary.table[i,2]/summary.table[i,4]
}

boxplot(summary.table, names=c("mean.dif", "diff.btw", "max sd (X|Y)", "diff.win", "m/sd"
, "btw/win"))
