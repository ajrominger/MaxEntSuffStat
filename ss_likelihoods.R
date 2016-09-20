bci <- read.csv('~/Research/datasets/stri/all/BCIS/BCIS.csv')
bci <- bci[bci$year == max(bci$year), ]

## solution for nu under SSNT
nuSSNT <- function(nu, nbar) {
    (1 - nu) / (nu * log(1/nu)) - nbar
}

## likelihood function for SSNT
ssntLogLik <- function(n, M) {
    S <- length(n)
    nbar <- mean(n)
    Mbar <- mean(M)
    
    m0 <- Mbar/nbar
    nu <- uniroot(nuSSNT, interval = c(.Machine$double.eps, 10), 
                  nbar = nbar, tol = .Machine$double.eps)$root
    
    S * log(1 / (m0*log(1/nu))) * sum(log() - M/m0)
}

