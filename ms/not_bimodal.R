mat <- cbind(rnorm(1e6), rnorm(1e6))
pat <- cbind(rnorm(1e6), rnorm(1e6)) * (runif(1e6) < 0.5)
dist <- sqrt( rowSums((mat+pat)^2) )

pdf(file="not_bimodal.pdf", width=4, height=2.5, pointsize=10)
par(mar=c(4,4,1,1)+.1)
    hist(dist, breaks=50, main=expression(abs(M + P * theta)), xlab='dispersal distance')
dev.off()
