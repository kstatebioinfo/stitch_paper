#!/usr/local/bin/Rscript
# Makes a graph from flowcell dates and/or throughput

pdf('/home/irys/stitch_paper/figures/Fig_4_flowcell_date.pdf', bg='white', width=9, height=9)
flowcellmetrics <- read.csv("/home/irys/stitch_paper/figures/Supplemental_2_bnx_date_and_metrics.csv", header=TRUE)

#print $temp_bnx_lengths "date,input_bnx,total_flowcell_length_Mb,bnx_count\n";

ballgownCol <- rgb(0,0,1,3/4)

plot(as.Date(flowcellmetrics$date,"%Y-%m-%d"), flowcellmetrics$total_flowcell_length_Mb, xlab= "Month", ylab= "Cumulative length per BNX file (Mb)", pch=4,col=c(ballgownCol),cex=2, type="p",cex.lab=1,lwd=3)
abline(v=as.Date("2013-10-18"),col="grey",lty=2, lwd=3)

legend( "top",legend=c("Cumulative length (Mb)","V2 chip upgrade"),lty=c(NA,2),pch=c(4,NA), pt.cex=c(2,NA),col=c(ballgownCol,"grey"), lwd=3 )

# V2 installed on 2013-10-18

dev.off()

