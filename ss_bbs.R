## ================================
## script to run SS theories on BBS
## ================================

library(parallel)

setwd('~/Dropbox/Research/MaxEntSuffStat')
source('ss_likelihoods.R')

## load one year of BBS data
bbsYear <- 2009
bbs <- read.csv(sprintf('~/Research/datasets/bbs/db/bbs%s.csv', bbsYear), as.is = TRUE)

## load bird body size and combine with bbs data.frame
bbsBody <- read.csv('~/Research/datasets/bbs/db/speciesTableBody.csv', as.is = TRUE)
bbs$avgMass <- bbsBody$mass[match(bbs$spp, bbsBody$sppKey)]
bbs$totMass <- bbs$avgMass * bbs$abund

## write both to local directory so they're availible on github
write.csv(bbs, 'bbs2009.csv', row.names = FALSE)
write.csv(bbsBody, 'speciesTableBody.csv.csv', row.names = FALSE)

## run on all routes
bbsSSLogLik <- mclapply(split(bbs, bbs$route), mc.cores = 6, FUN = function(d) {
    ssnt <- ssntLogLik(d$abund, d$totMass)
    ssme <- ssmeLogLik(d$abund, d$totMass)
    ssnti <- ssntiLogLik(d$abund, d$totMass)
    ssmei <- ssmeiLogLik(d$abund, d$totMass)
    
    return(c(ssnt = ssnt,
             ssme = ssme,
             ssnti = ssnti,
             ssmei = ssmei,
             tot = ssnt - ssme,
             ind = ssnti - ssmei))
})

bbsSSLogLik <- data.frame(route = names(bbsSSLogLik), do.call(rbind, bbsSSLogLik))

pdf('fig_ssLogLikBBS.pdf', width = 9, height = 3)

par(mfrow = c(1, 3), mar = c(4, 4, 2, 1) + 0.1, mgp = c(2.5, 0.75, 0))
plot(-log(-bbsSSLogLik$ssme), -log(-bbsSSLogLik$ssnt), 
     xlab = expression(-log(-logLik['SSME'])), 
     ylab = expression(-log(-logLik['SSNT'])))
abline(0, 1, col = 'red')
mtext('Total mass', line = 1)

plot(-log(-bbsSSLogLik$ssmei), -log(-bbsSSLogLik$ssnti),
     xlab = expression(-log(-logLik['SSMEI'])), 
     ylab = expression(-log(-logLik['SSNTI'])))
abline(0, 1, col = 'red')
mtext('Individual mass', line = 1)

plot(bbsSSLogLik$tot, bbsSSLogLik$ind, 
     xlab = expression(logLik['SSNT'] - logLik['SSME']), 
     ylab = expression(logLik['SSNTI'] - logLik['SSMEI']))

dev.off()


## write out for james
write.csv(bbsSSLogLik[, 1:5], file = 'bbs_sizeStructured_ll.csv', row.names = FALSE)


## example figure of neutral ll - maxent ll
pdf('fig_eq9_ssInd.pdf', width = 4, height = 4)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
hist(bbsSSLogLik$ssnti - bbsSSLogLik$ssmei, 
     xlab = expression(ll[SSNTI] - ll[SSMEI]), main = '')
rect(xleft = -2, xright = 2, ybottom = par('usr')[3], ytop = par('usr')[4], 
     col = gray(0.5, 0.5), border = NA)
dev.off()

pdf('fig_eq9_ssTot.pdf', width = 4, height = 4)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
hist(bbsSSLogLik$ssnt - bbsSSLogLik$ssme, 
     xlab = expression(ll[SSNT] - ll[SSME]), main = '', 
     xlim = range(bbsSSLogLik$ssnt - bbsSSLogLik$ssme, -2, 2))
rect(xleft = -2, xright = 2, ybottom = par('usr')[3], ytop = par('usr')[4], 
     col = gray(0.5, 0.5), border = NA)
dev.off()
