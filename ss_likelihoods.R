bci <- read.csv('~/Research/datasets/stri/all/BCIS/BCIS.csv', as.is = TRUE)
bci <- bci[bci$year == max(bci$year), ]
bciSpLevel <- aggregate(bci[, c('count', 'dbh')], list(spp = bci$spp), sum)

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
    
    S * log(1 / (m0*log(1/nu))) + 
        sum(n*log(1-nu) - lfactorial(n) + (n-1)*log(M/m0) - M/m0)
}

## likelihood function for SSME
ssmeLogLik <- function(n, M) {
    S <- length(n)
    nbar <- mean(n)
    Mbar <- mean(M)
    la1 <- log(nbar / (nbar - 1))
    la2 <- 1/Mbar
    
    S * log(la2 * (exp(la1) - 1)) - sum(la1*n + la2*M)
}


ssntLogLik(bciSpLevel$count, bciSpLevel$dbh)
ssmeLogLik(bciSpLevel$count, bciSpLevel$dbh)



## likelihood function for SSNTI
ssntiLogLik <- function(spp, n, m) {
    S <- length(unique(spp))
    nbar <- sum(n) / S
    nu <- uniroot(nuSSNT, interval = c(.Machine$double.eps, 10), 
                  nbar = nbar, tol = .Machine$double.eps)$root
    
    dat <- split(data.frame(n = n, m = m), spp)
    
    lp <- sapply(dat, function(d) {
        m0 <- sum(d$m) / nbar
        n <- sum(d$n)
        
        n*log(1 - nu) - log(n) - n*log(m0) - sum(d$m/m0)
    })
    
    return(S * log(1/log(1/nu)) + sum(lp))
}

ssntiLogLik(bci$spp, bci$count, bci$dbh)


