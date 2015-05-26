
# Generate pretty graphs from Codon usage test data

#Set the working directory to the directory in which the script is run.
this.dir <- dirname(parent.frame(2)$ofile)
set.seed(9373399)
setwd(this.dir)
setwd("../results/")

library("ggplot2")
library("reshape2")

files <- list.files(pattern="*.txt")

histogramPlot = function(){
    
  # Build a plot with one observed value and its matching expected distribution    
  single_plot <- (ggplot() + aes(dens) + geom_density(fill="gray87", size=0.7)
  + xlim((a[1,1]-xScale), (a[1,2]+xScale))
  + theme_bw()
  + theme(axis.title = element_text(size = rel(1.3)))
  + theme(axis.text = element_text(size = rel(1.1)))
  + xlab("Entropy") + ylab("Density")
  + ggtitle(paste(b, "Entropy Distribution"))
  + theme(plot.title = element_text(size = rel(1.6))))
  
  # Annotate the plot with relative positioning of text and bars
  yrange <- ggplot_build(single_plot)$panel$ranges[[1]]$y.range
  xrange <- ggplot_build(single_plot)$panel$ranges[[1]]$x.range
  single_plot <- single_plot +   
      annotate("segment", x = a[1,1], xend = a[1,1], y = 0, yend = yrange[2]/2.1, 
                color="black", size=1.0, linetype="longdash") +
      annotate("text", x = mean(c(a[1,1], a[1,2])), y = yrange[2]/1.7, 
               label = sd_away, size = 8) + 
      annotate("errorbarh", x = mean(c(a[1,1], a[1,2])), xmin = a[1,1], 
             xmax = mean(dens), y = yrange[2]/1.9, height = 0.004, size = 0.75)
  print(single_plot)

}

for (i in files){
    a <- read.table(i)
    b <- gsub(".txt", "", i)
    colnames(a) <- c("observed","expectedMean","variance")
    dens <- rnorm(100, a[1,2], sqrt(a[1,3]))
    maxDens <- dnorm(mean(dens), mean(dens), sqrt(a[1,3]))
    sd_away <- paste(round(((a[1,2] - a[1,1]) / sqrt(a[1,3])), 1), "SD")
    xScale <- (mean(dens) - a[1,1]) * 0.07
    histogramPlot()
}