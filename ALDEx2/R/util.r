aitchison.mean <- function( n, log=FALSE ) {

    # Input is a vector of non-negative integer counts.
    # Output is a probability vector of expected frequencies.
    # If log-frequencies are requested, the uninformative subspace is removed.

    n <- round( as.vector( n, mode="numeric" ) )
    if ( any( n < 0 ) ) stop("counts cannot be negative")

    a <- n + 0.5
    sa <- sum(a)

    log.p <- digamma(a) - digamma(sa)
    log.p <- log.p - mean(log.p)

    if ( log ) return(log.p)

    p <- exp( log.p - max(log.p) )
    p <- p / sum(p)
    return(p)
}


ilr.matrix <- function( d ) {

    # Returns a '(d-1)-by-d' matrix 'u' such that the Isometric
    # Log-Ratio (ILR) of 'log(p)' is 'u %*% log(p)'.

    d <- round( as.vector( d, mode="numeric" ) )[1]
    if ( d < 2 ) stop("dimension must be at least two")

    d1 <- d - 1

    u <- matrix(NA,nrow=d1,ncol=d)
    for ( i in 1 : d1 ) {
        u[i,] = sqrt(i/(i+1)) * c( rep(1,i)/i , -1 , rep(0,d1-i) )
    }

    return(u)
}


dt3 <- function( x, df, location=0, precision=1, log=FALSE ) {

    # Compute the three-parameter Student-t density. If 'mu' is the location,
    # 'lambda' is the precision (squared inverse scale), and 'nu' is the degrees
    # of freedom, then:
    #
    #   mean(t) = mu, for nu > 1, and
    #    var(t) = nu/(nu-2)/lambda, for nu > 2

    nu     <- df
    mu     <- location
    lambda <- precision

    log.d <- lgamma((nu+1)/2) - lgamma(nu/2) + 0.5 * ( log(lambda) - log(pi*nu) ) - (nu+1)/2 * log( 1 + lambda*(x-mu)^2/nu )

    if ( log ) return(log.d)
    return(exp(log.d))
}


rt3 <- function( n, df, location=0, precision=1 ) {

    # Samples the three-parameter Student-t density. If 'mu' is the location,
    # 'lambda' is the precision (squared inverse scale), and 'nu' is the degrees
    # of freedom, then:
    #
    #   mean(t) = mu, for nu > 1, and
    #    var(t) = nu/(nu-2)/lambda, for nu > 2

    nu     <- df
    mu     <- location
    lambda <- precision

    ts <- rt( n, df ) / sqrt(lambda) + location

    return(ts)
}

