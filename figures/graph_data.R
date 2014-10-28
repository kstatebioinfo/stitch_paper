#!/usr/local/bin/Rscript

# Makes a histogram from a CSV list of gap lengths

pdf('/Users/jennifer_shelton/Desktop/stitch_paper/figures/gap_length_histogram.pdf', bg='white', width=15, height=5)
neggapvalue <- read.table("/Users/jennifer_shelton/Desktop/stitch_paper/raw_data/neg_gaps_pre_manual.tab", header=TRUE)
posgapvalue <- read.table("/Users/jennifer_shelton/Desktop/stitch_paper/raw_data/pos_gaps_pre_manual.tab", header=TRUE)

p1 <- hist(neggapvalue$negative_gaps,breaks=seq(-1680000,1300000,by=20000) )

p2 <- hist(posgapvalue$positive_gaps,breaks=seq(-1680000,1300000,by=20000) )

ballgownCol <- rgb(0,0,1,3/4)
cuffdiffCol <- rgb(1,0,0,3/4)


#plot( p1, col=c(ballgownCol),ylim=c(0,20), xlim=c(-1680000,1300000),main="Distribution of gap lengths for automated output",xlab="Gap length (bp)",ylab="Count")  # first histogram
plot( p1, col=c(ballgownCol),ylim=c(0,20), xlim=c(-1680000,1300000),main="",xlab="Gap length (bp)",ylab="Count")  # first histogram
plot( p2, col=c(cuffdiffCol),ylim=c(0,20), xlim=c(-1680000,1300000), add=TRUE)  # first histogram

abline(v=0,col="red")

legend("topright", legend=c("Negative gap lengths","Positive gap lengths"), fill=c(ballgownCol,cuffdiffCol))


dev.off()