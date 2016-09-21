## testing ss likelihoods on BCI data

setwd('~/Dropbox/Research/MaxEntSuffStat')
source('ss_likelihoods.R')

bci <- read.csv('~/Research/datasets/stri/all/BCIS/BCIS.csv', as.is = TRUE)
bci <- bci[bci$year == max(bci$year), ]
bciSpLevel <- aggregate(bci[, c('count', 'dbh')], list(spp = bci$spp), sum)

ssntLogLik(bciSpLevel$count, bciSpLevel$dbh)
ssmeLogLik(bciSpLevel$count, bciSpLevel$dbh)
ssntiLogLik(bciSpLevel$count, bciSpLevel$dbh)
ssmeiLogLik(bciSpLevel$count, bciSpLevel$dbh)
