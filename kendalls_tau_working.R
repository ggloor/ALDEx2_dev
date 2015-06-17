# /Volumes/seq/LRGC/sailorsam/vitaminD_NZ/analysis_nzdp_nocrohns2

d <- read.table("meta_VitD_numerical_HP.txt", header=T, row.names=1, comment.char="")
otu <- read.table("Filter_Manual/OTU_table_filter1.txt", header=T, skip=1, comment.char="", sep="\t", row.names=1)

# data table of nocrohns2 set using rownames from d table
otu.d <- data.frame(otu[,rownames(d)])

# remove column
vit.d <- d$VitD_Dose
d$VitD_Dose <- NULL

# change to a real NA value
d[d=="ND"] <- NA

# make taxon variable for convenience
tax <- otu$taxonomy
 
library(ALDEx2)

# just need the clr function
x <- aldex.clr(otu.d)

# class functions to get data parameters about x
mc.instances <- numMCInstances(x)
n.features <- numFeatures(x)
n.samples <- length(getSampleIDs(x))
mean.cor <- matrix(data=0, nrow=n.features, ncol=n.features)
bh.cor <- matrix(data=0, nrow=n.features, ncol=n.features)
p.cor <- matrix(data=0, nrow=n.features, ncol=n.features)

# function to convert kendalls tau values to p values
kendallz <- function(T,n){( 3 * T * sqrt(n * (n-1)) ) / sqrt( 2*(2*(n+5)) )}
kendallp <-  function(Z){2*pnorm(-abs(Z))}

# set up the output matrices
# row contain samples, columns contain OTUs
mean.cor <- matrix(data=0, nrow=ncol(d), ncol=n.features)
bh.cor <- matrix(data=0, nrow=ncol(d), ncol=n.features)

colnames(mean.cor) <- rownames(otu.d)
colnames(bh.cor) <- rownames(otu.d)
rownames(mean.cor) <- colnames(d)
rownames(bh.cor) <- colnames(d)

# OTUs correlated with metadata using kendalls tau
for(j in 1:ncol(d)){ # each metadata
for(i in 1:mc.instances){ # each dir instance
#	print(c(j,i)) 
    # get the dir instance
	t.input <- sapply(getMonteCarloInstances(x), function(y){y[,i]})
	# calculate kendalls tau, the as.numeric() and use=  are important
	cor.rho <- apply(t.input,1 , function(x){cor(x, as.numeric(d[,j]), method="kendall", use="pairwise.complete.obs")})	
	# p and BH values
	cor.t <- kendallz(cor.rho, n.samples)
	cor.p <- kendallp(cor.t) # two-tailed test
	adj.cor <- p.adjust(cor.p) 
	
	# tabulate the values
	bh.cor[j,] <- bh.cor[j,] + ( adj.cor / mc.instances )
	mean.cor[j,] <- mean.cor[j,] + cor.rho / mc.instances
	
}
}

# turn the output into a data frame with samples in columns
t.bh <- as.data.frame(t(bh.cor))
t.cor <- as.data.frame(t(mean.cor))

# list the taxa that are differential at the BH value given
tax[which(t.bh$Gingival_Index < 0.1)] 
tax[which(t.bh$Plaque_Index < 0.1)]
tax[which(t.bh$Plaque_Accumulation < 0.1)]  
tax[which(t.bh$Age < 0.1)]
tax[which(t.bh$Frequency < 0.1)]
tax[which(t.bh$Severity < 0.1)]

# just list the line numbers (easier to compare?)
which(t.bh$Gingival_Index < 0.1) 
which(t.bh$Plaque_Index < 0.1)
which(t.bh$Plaque_Accumulation < 0.1)  
which(t.bh$Age < 0.1)
which(t.bh$Frequency < 0.1)
which(t.bh$Severity < 0.1)