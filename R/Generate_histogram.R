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

file <- "GC.All.txt"

histogramPlot = function(file){

  # Read in data frame from standard script output
  a <- read.table(file, quote = "\"")
  b <- gsub(".txt", "", file)
  colnames(a) <- c("observed","expectedMean","variance")
  
  # Generate random draws from the normal distribution specified in the first line and plot
  y1 <- rnorm(100, a[1,2], sqrt(a[1,3]))
  sd_away <- (a[1,2] - a[1,1]) / sqrt(a[1,3])
  
  plot1 <- (ggplot() + aes(y1) + geom_density(fill="gray87", size=1.0)
  + xlim((a[1,1]-100), (a[1,2]+100))
  + theme_bw()
  + theme(axis.title = element_text(size = rel(1.3)))
  + theme(axis.text = element_text(size = rel(1.1)))
  + annotate("segment", x = a[1,1], xend = a[1,1], y = 0, yend = 0.035, color="black", size=1.0, linetype="longdash") 
  + annotate("text", x = (a[1,2]+100)/2, y = 0.05, label = sd_away)
  + annotate("errorbarh", x = 1250, xmin = a[1,1], xmax = mean(y1), y = 0.04, height = 0.004, size = 0.75)
  + xlab("Entropy") + ylab("Density")
  )
  
  #ggsave(file = "plot1.png", plot = plot1, width = 10, height = 6)
  return(plot1)
#   # Melt and cast the data for graphing purposes
#   #a$variance = NULL
#   c = melt(a, measured=c("observed", "expectedMean"))
#   colnames(c) = c("Distribution", "Entropy")
#   
#   plot2 <- ggplot(c, aes(x=Entropy, fill=Distribution)) + geom_histogram(binwidth=15)+ theme(plot.title = element_text(lineheight=.8, face="bold")) + scale_fill_manual(values=c("purple", "aquamarine4"))
#   plot2
}
histogramPlot(file)

# for (i in 1:length(files)){
#   histogramPlot(files[i])
# }