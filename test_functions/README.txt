# How to run the test functions

source("clr_test_function.r")
source("ttest_test_function.r")

# run clr_test_function.r
# this is multiprocessor aware
x <- clr(data_table, mc.samples=128, verbose=TRUE)

#INPUT
# reads a data table with feature counts in rows, and samples in columns

# OUTPUT
# The output returned is a list (x) that contains Monte-Carlo instances of
# the centre log-ratio transformed values for each sample
# Access to values 
# sample IDs: names(x) 
# number of features (genes, OTUs): length(x[[1]][,1])
# number of Monte-Carlo Dirichlet instances: length(x[[1]][1,])
# feature names: rownames(x[[1]])

# run ttest_test_function.r
# input is output from clr_test_function.r
# this is not multiprocessor aware
x.tt <- clr.ttest(x,conditions)

# if using the selex dataset, conditions=c(rep("N", 7), rep("S",7))

# OUTPUT
# dataframe containing welch's expected P, welch's expected BH corrected P
# wilcoxon expected P, wilcoxom expected BH corrected P
