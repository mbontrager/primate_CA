
# Generate pretty graphs from Codon usage test data

#Set the working directory to the directory in which the script is run.
this.dir <- dirname(parent.frame(2)$ofile)
set.seed(9373399)
setwd(this.dir)
setwd("../results/")

library("ggplot2")
library("reshape2")

files <- list.files(pattern="*.txt")

histogramPlot <- function(){
    
  # Build a plot with one observed value and its matching expected distribution    
  single_plot <- (ggplot() + aes(dens) + geom_density(fill="gray75", size=0.7)
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

total <- data.frame(Observed=numeric(), 
                    Expected=numeric(), 
                    Variance=numeric(), 
                    Dist=factor(),
                    Codons=factor())

for (i in files){
    a <- read.table(i)
    b <- gsub(".txt", "", i)
    d <- unlist(strsplit(b, "[.]"))
    colnames(a) <- c("Observed","Expected","Variance")
    a$Dist <- as.factor(d[1])
    a$Codons <- as.factor(d[2])
    
    # Scale plots to similar visual dimensions and annotation positions
    dens <- rnorm(100, a[1,2], sqrt(a[1,3]))
    sd_away <- paste(round(((a[1,2] - a[1,1]) / sqrt(a[1,3])), 1), "SD")
    xScale <- (mean(dens) - a[1,1]) * 0.07 # Scale the x-axis relatively
    
    # Create the plot
    histogramPlot()

    #Merge separate data frames
    total <- rbind(total, a)
}



## Begin second round of plotting
total$Variance <- NULL
total <- melt(total, measured=c("Observed", "Expected"))
colnames(total) <- c("Distribution", "Codons", "Measure", "Entropy")
plot2 <- ggplot(total, aes(x=Entropy, fill=Measure)) +
    geom_histogram(binwidth = 10) +
    scale_fill_manual(values=c("gray75", "black")) +
    theme_bw() +
    facet_grid(Distribution ~ Codons, scale = "free") +
    ggtitle("Distance Observed to Expected Entropy") +
    theme(axis.title = element_text(size = rel(1.2)), 
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          strip.text.x = element_text(size=rel(1.1)),
          strip.text.y = element_text(size=rel(1.1)),
          plot.title = element_text(size = rel(1.45), face = "bold", vjust = 1))

## Create data frame with (x, y) coordinates for annotations
#Length of the data frame
len <- length(levels(total$Distribution)) * length(levels(total$Codons))
# Declare the data frame
vars <- data.frame(expand.grid(levels(total$Distribution), 
                               levels(total$Codons)))
colnames(vars) <- c("Distribution", "Codons")
# Populate the new data frame with (x, y) coordinates:
dat <- data.frame(x = rep(800, len), y = rep(4.5, len), vars, 
                  labs = sd_away)

print(plot2)
