## ================================
## script to run SS theories on BBS
## ================================

setwd('~/Dropbox/Research/MaxEntSuffStat')
source('ss_likelihoods.R')

## load one year of BBS data
bbsYear <- 2009
bbs <- read.csv(sprintf('~/Research/datasets/bbs/db/bbs%s.csv', bbsYear), as.is = TRUE)

## load bird body size and combine with bbs data.frame
bbsBody <- read.csv('~/Research/datasets/bbs/db/speciesTableBody.csv', as.is = TRUE)
bbs$avgMass <- bbsBody$mass[match(bbs$spp, bbsBody$sppKey)]
bbs$totMass <- bbs$avgMass * bbs$abund
