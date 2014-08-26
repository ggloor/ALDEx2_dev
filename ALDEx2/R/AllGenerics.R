### -------------------------------------------------------------------------
### aldex.clr
###

setGeneric("getMonteCarloInstances", function(.object, ...) standardGeneric("getMonteCarloInstances"))

setGeneric("getSampleIDs", function(.object, ...) standardGeneric("getSampleIDs"))

setGeneric("numFeatures", function(.object, ...) standardGeneric("numFeatures"))

setGeneric("numMCInstances", function(.object, ...) standardGeneric("numMCInstances"))

setGeneric("getFeatureNames", function(.object, ...) standardGeneric("getFeatureNames"))

setGeneric("getReads", function(.object, ...) standardGeneric("getReads"))

setGeneric("numConditions", function(.object, ...) standardGeneric("numConditions"))

setGeneric("getMonteCarloReplicate", function(.object, i, ...) standardGeneric("getMonteCarloReplicate"))
