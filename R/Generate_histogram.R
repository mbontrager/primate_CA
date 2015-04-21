# Generate pretty graphs from Codon usage test data

#Set the working directory to the directory in which the script is run.
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

#Load and/or install require packages.
PkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

PkgTest("ggplot2")
PkgTest("reshape2")

setwd("../results/")

histogramPlot = function(file){
  
  # Read in data frame from standard script output
  a <- read.table(file, quote = "\"", skip = 1)
  b <- gsub(".txt", "", file)
  colnames(a) <- c("observed","expectedMean","variance")
  
  # Generate random draws from the normal distribution specified in the first line and plot
  y = rnorm(200,a[1,2], sqrt(a[1,3]))
  plot1 <- ggplot() + aes(y) + geom_density(fill="burlywood")
  plot1 <- plot1 + xlim((a[1,1]-100), (a[1,2]+100))
  #plot1 <- plot1  + ggtitle("GC codon usage \n One taxa set (200 trials)") + theme(plot.title = element_text(size = rel(1.5), face="bold")) + theme(axis.title = element_text(size = rel(1.5)))
  plot1 <- plot1 + geom_vline(xintercept=a[1,1], color="red", size=1.0)
  plot1
  
  # Melt and cast the data for graphing purposes
  #a$variance = NULL
  c = melt(a, measured=c("observed", "expectedMean"))
  colnames(c) = c("Distribution", "Entropy")
  
  plot2 <- ggplot(c, aes(x=Entropy, fill=Distribution)) + geom_histogram(binwidth=15)+ theme(plot.title = element_text(lineheight=.8, face="bold")) + scale_fill_manual(values=c("purple", "aquamarine4"))
  plot2
}
histogramPlot("GC.All.txt")

# for (i in 1:length(files)){
#   histogramPlot(files[i])
# }