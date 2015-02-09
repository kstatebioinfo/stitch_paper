#!/usr/local/bin/Rscript
# Makes a histogram from BioNano genome map assembly metrics

pdf('/Users/jennifer_shelton/Desktop/stitch_paper/figures/Fig_4_assembly_metrics.pdf', bg='white', width=9, height=9)
assemblymetrics <- read.csv("/Users/jennifer_shelton/Desktop/stitch_paper/figures/Supplemental_2_Assembly_parameter_tests.csv", header=TRUE)

y1 <- (assemblymetrics$Breadth_of_alignment)

y2 <- (assemblymetrics$Total_alignment_length)

y3 <- (assemblymetrics$Cumulative_length)

x1 <- factor(assemblymetrics$"Genome_map")

x1 = factor(x1,levels(x1)[c(2,1,4,3,5)])

ballgownCol <- rgb(0,0,1,3/4)

plot(x1, ylim=c(100,300),main="",xlab="Assembly name",ylab="Length (Mb)")
axis(1,tick = TRUE, labels = FALSE)

abline(h=200,col="grey",lty=2, lwd=3)

points(x1,y1, pch=18,col=c(ballgownCol),cex=3)
points(x1,y2, pch=25,col=c(ballgownCol),cex=2,lwd=3)
points(x1,y3, pch=1,cex=2, type="b",lwd=4, col="yellow2",cex.lab=1)


legend( "top",legend=c("Cumulative length","Estimated genome size","Total aligned length","Breadth of alignment"),lty=c(1,2,NA,NA),pch=c(1,NA,25,18), pt.cex=c(2,NA,2,3),col=c("yellow2","grey",ballgownCol,ballgownCol), lwd=3 )


dev.off()

