# Generate pretty graphs from Codon usage test data

#Load and/or install require packages.
PkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

PkgTest(ggplot2)
PkgTest(reshape)

# Read in data frame from standard script output
GC_all <- read.table("D:/Users/Martin/Dropbox/primate_common_ancestry/primate_sequence_alignments/NEXUS/trimmed/200_codonGC_all.txt", quote="\"")
colnames(GC_all) <- c("observed","expectedMean","variance")

# Generate random draws from the normal distribution specified in the first line and plot
y = rnorm(200,GC_all[1,2], sqrt(GC_all[1,3]))
plot1 <- ggplot() + aes(y) + geom_density(fill="burlywood")
plot1 <- plot1 + xlim((GC_all[1,1]-100), (GC_all[1,2]+100))
plot1 <- plot1  + ggtitle("GC codon usage \n One taxa set (200 trials)") + theme(plot.title = element_text(size = rel(1.5), face="bold")) + theme(axis.title = element_text(size = rel(1.5)))
plot1 <- plot1 + geom_vline(xintercept=GC_all[1,1], color="red", size=1.0)
plot1

# Melt and cast the data for graphing purposes
GC_all$variance = NULL
a = melt(GC_all, measured=c("observed", "expectedMean"))
colnames(a) = c("Distribution", "Entropy")

plot2 <- ggplot(a, aes(x=Entropy, fill=Distribution)) + geom_histogram(binwidth=15)+ ggtitle("GC all codons") + theme(plot.title = element_text(lineheight=.8, face="bold")) + scale_fill_manual(values=c("purple", "aquamarine4"))
plot2
