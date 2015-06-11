#!/usr/local/bin/Rscript
# USAGE: Rscript histograms.R <assembly_metrics_file>
# DESCRIPTION: Runs PCA for metrics for all molecule map assemblies
library("ggplot2")
library("ggbiplot")
#pdf("/Users/jennifer_shelton/Dropbox/Project_Bionano/Additional_fig_X_test.pdf", bg='white' )
pdf("/Users/jennifer_shelton/Dropbox/Project_Bionano/Additional_fig_6_approved.pdf", bg='white' )
## Molecule map assemblies ##
#assemblies <- read.csv("/Users/jennifer_shelton/Dropbox/Project_Bionano/all_assemblies_PCA.csv", header=TRUE)
assemblies <- read.csv("/Users/jennifer_shelton/Dropbox/Project_Bionano/all_assemblies_PCA_approved.csv", header=TRUE)
# names(assemblies)
#  [1] "Type"
#  [2] "Organism"
#  [3] "Genome_FASTA_N50_Mb"
#  [4] "Genome_FASTA_Count"
#  [5] "Genome_FASTA_Length_Mb"
#  [6] "Genome_Map_N50_Mb"
#  [7] "Genome_Map_Count"
#  [8] "Genome_Map_Length_Mb"
#  [9] "Super_Scaffold_Genome_FASTA_N50_Mb"
# [10] "Super_Scaffold_Count"
# [11] "Super_Scaffold_Length_Mb"
# [12] "N50_Percent_Increase"
# [13] "Percent_in_silico_maps_aligned"
# [14] "Percent_genome_FASTA_aligned_with_molecule_maps"
# [15] "Project_ID"
# [16] "X_Molecule_Map_Coverage_Above_150_kb"
# [17] "Row_Name"
# [18] "Genus"
# [19] "Labels_per_100_kb"
## Subset draft and reference based assemblies ##
draft_and_reference_assemblies <- subset(assemblies, Type == "Draft" | Type == "Reference")
## Subset draft and reference based assemblies skipping E. coli ##
draft_and_reference_assemblies_no_e_coli <- subset(assemblies, Organism != "Escherichia coli")
subset_values <- (draft_and_reference_assemblies[, c(1,2,6,13,16,19)])
# for standard PCA
draft_and_reference_assemblies_for_pca <- subset_values[, 3:6]
# log transform
log.draft_and_reference_assemblies <- log(subset_values[, 3:6])
# Grab project id
rownames(log.draft_and_reference_assemblies) <- draft_and_reference_assemblies$Row_Name
rownames(draft_and_reference_assemblies) <- draft_and_reference_assemblies$Row_Name
# Grab genus
groups_Genus <- draft_and_reference_assemblies$Genus
####################################
# PCA log transformed (because variance increases with increasing values of X for most variables)
####################################
# apply PCA - scale. = TRUE is highly
# advisable, but default is FALSE.
log.draft_and_reference_assemblies.pca <- prcomp(log.draft_and_reference_assemblies, center = TRUE, scale. = TRUE, labels = draft_and_reference_assemblies$Row_Name)
# plot the eigenvectors
#plot(log.draft_and_reference_assemblies.pca, type = "l", main = "PCA of log transformed assembly metrics")
# print rotation and STD
print(log.draft_and_reference_assemblies.pca)
# print Variance explained
summary(log.draft_and_reference_assemblies.pca)
# biplot(log.draft_and_reference_assemblies.pca, col=c(2,3), cex=c(1/2, 3/4))
## make biplot grouped by genus
g <- ggbiplot(log.draft_and_reference_assemblies.pca, obs.scale = 1, var.scale = 1, groups = groups_Genus,labels = draft_and_reference_assemblies$Genus, circle = TRUE,labels.size = 1.5, varname.size = 2,varname.adjust = 1)
# varname.abbrev
g <- g + coord_cartesian(xlim = c(-6.0, 3.0))
g <- g + coord_cartesian(ylim = c(-4.0, 4.0))
g <- g + scale_color_discrete(name = '')
# Change title appearance
g <- g + labs(title = "Biplot for PCA of log transformed assembly metrics")
print(g)
#############################################################################
# Draft only
#############################################################################
## Subset draft and reference based assemblies ##
draft_assemblies <- subset(assemblies, Type == "Draft" )
# ## Subset draft and reference based assemblies skipping E. coli ##
# draft_and_reference_assemblies_no_e_coli <- subset(assemblies, Organism != "Escherichia coli")
subset_values <- (draft_assemblies[, c(1,2,3,5,6,12,13,16,19)])
# for standard PCA
draft_assemblies_for_pca <- subset_values[, 3:9]
# log transform
subset_values$N50_Percent_Increase <- (subset_values$N50_Percent_Increase + 1)

log.draft_assemblies <- log(subset_values[, 3:9])
# Grab project id
rownames(log.draft_assemblies) <- draft_assemblies$Row_Name
rownames(draft_assemblies) <- draft_assemblies$Row_Name
# Grab genus
groups_Genus_draft <- draft_assemblies$Genus
####################################
# PCA log transformed
####################################
# apply PCA - scale. = TRUE is highly
# advisable, but default is FALSE.
log.draft_assemblies.pca <- prcomp(log.draft_assemblies, center = TRUE, scale. = TRUE, labels = draft_assemblies$Row_Name)
# plot the eigenvectors
# plot(log.draft_assemblies.pca, type = "l", main = "PCA of log transformed assembly metrics (Draft only)")
# print rotation and STD
print(log.draft_assemblies.pca)

# print Variance explained
summary(log.draft_assemblies.pca)
# biplot(log.draft_and_reference_assemblies.pca, col=c(2,3), cex=c(1/2, 3/4))
## make biplot grouped by genus
g <- ggbiplot(log.draft_assemblies.pca, obs.scale = 1, var.scale = 1 , groups = groups_Genus_draft,varname.size = 2,varname.adjust = 1, circle = TRUE, labels = draft_assemblies$Genus, labels.size = 1.5)
# # varname.abbrev
g <- g + coord_cartesian(xlim = c(-6.0, 3.0))
g <- g + coord_cartesian(ylim = c(-4.0, 4.0))
g <- g + scale_color_discrete(name = '')
# # Change title appearance
g <- g + labs(title = "Biplot for PCA of log transformed assembly metrics (Draft only)")
print(g)
# Reset values to include 0
subset_values$N50_Percent_Increase <- (subset_values$N50_Percent_Increase - 1)
####################################
# Scatterplots (formatting from http://stackoverflow.com/questions/15271103/how-to-modify-this-correlation-matrix-plot and http://personality-project.org/r/r.graphics.html
####################################
panel.regression <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
cex = 1, col.regres = "red", ...)
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok))
    abline(stats::lm(y[ok] ~ x[ok]), col = col.regres, ...)
}
##   put correlations on the upper panels,
## with size proportional to the (absolute) correlations.
#first create a function (panel.cor)
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y)
    r.abs <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    
    test <- cor.test(x,y)
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " "))
    
    text(0.5, 0.5, txt, cex = cex * (r.abs + 0.4))
    text(.8, .8, Signif, cex=cex, col=2)
}
subset_values <- (draft_and_reference_assemblies[, c(3,5,6,13,16,19)])
names(subset_values) <- c("FASTA_N50","FASTA_Length","Genome_Map_N50","Percent_aligned","Molecule_Coverage", "Labels_per_100_kb")
# pairs(subset_values, main= "Correlation matrix for assembly metrics", lower.panel=panel.regression, upper.panel=panel.cor)
log.subset_values <- log(subset_values)
pairs(log.subset_values, main= "Correlation matrix for log transformed assembly metrics", lower.panel=panel.regression, upper.panel=panel.cor)
# names(log.draft_assemblies)
# [1] "Genome_FASTA_N50_Mb"
# [2] "Genome_FASTA_Length_Mb"
# [3] "Genome_Map_N50_Mb"
# [4] "N50_Percent_Increase"
# [5] "Percent_in_silico_maps_aligned"
# [6] "X_Molecule_Map_Coverage_Above_150_kb"
# [7] "Labels_per_100_kb"
names(log.draft_assemblies) <- c("FASTA_N50","FASTA_Length","Genome_Map_N50","N50_Percent_Increase","Percent_aligned","Molecule_Coverage", "Labels_per_100_kb")
pairs(log.draft_assemblies  , main= "Correlation matrix for log transformed assembly metrics (Draft only)", lower.panel=panel.regression, upper.panel=panel.cor)
#######################################
# Histograms
#######################################
require(gridExtra)
plot_1 <- ggplot(subset_values, aes(x = subset_values$Genome_Map_N50)) + geom_histogram(aes(y = ..density..),colour = "black",fill = "lightblue" ) + geom_density() + xlab("Genome map N50 (Mb)") + ggtitle(paste("Assembly metrics for various genomes","\n\nGenome map N50 (draft and reference based assemblies)")) + geom_vline(xintercept=1.35, colour="black", linetype = "longdash") + geom_text(aes(x=1.35, label="Tcas", y=1, size=2.5,family="Helvetica", fontface = "plain"), colour="gray40", angle=90, vjust = -0.5) + theme(legend.position="none")
plot_2 <- ggplot(draft_assemblies, aes(x = draft_assemblies$Genome_FASTA_N50_Mb)) + geom_histogram(aes(y = ..density..),colour = "black",fill = "lightblue", binwidth = 0.5 ) + geom_density() + xlab("FASTA N50 (Mb)") + ggtitle("Original FASTA N50 (draft assemblies)") + geom_vline(xintercept=1.16, colour="black", linetype = "longdash") + geom_text(aes(x=1.16, label="Tcas5.0", y=0.50, size=2.5,family="Helvetica", fontface = "plain"), colour="gray40", angle=90, vjust = 1.5) + theme(legend.position="none")
draft_assemblies$N50_Percent_Increase <- draft_assemblies$N50_Percent_Increase *100
plot_3 <- ggplot(draft_assemblies, aes(x = draft_assemblies$N50_Percent_Increase)) + geom_histogram(aes(y = ..density..),colour = "black",,fill = "lightblue" ) + geom_density() + xlab("Percent of N50 increase") + ggtitle("Percent FASTA N50 increased after super scaffolding (draft assemblies)") + geom_vline(xintercept=269, colour="black", linetype = "longdash") + geom_text(aes(x=269, label="Tcas5.0", y=0.015, size=2.5,family="Helvetica", fontface = "plain"), colour="gray40", angle=90, vjust = -0.5) + theme(legend.position="none")

g1 <- ggplotGrob( plot_1 )
g2 <- ggplotGrob( plot_2 )
g3 <- ggplotGrob( plot_3 )
grid.arrange( g1, g2, g3, nrow = 3 )

dev.off()

