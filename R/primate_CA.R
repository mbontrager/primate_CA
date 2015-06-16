# Martin Bontrager
# Bret Larget
# Testing Common Ancestry at the family level in primates.

# Source script from `/R` directory

# This script takes input in the form of a list of fasta files aligned 
# in the correct reading frame, a mouse whole-genome codon usage table, the
# codon translation, and a large list of primate species (which match headers
# in the alignments). It returns a table with `repsPerGene` draws of the family
# representative. For each n in `repsPerGene` it calculates the observed 
# entropy, and draws n=trials times from three different codon usage 
# models to generate mean expected entropy and variance under the model.
# We also return the mean expected entropy stat if invariant degenerate
# positions are excluded

# User-defined Parameters------------------------------------------------------

alignedGenes <- "../data/aligned_coding_genes.txt" # Aligned gene fasta files
repsPerGene <- 50   # Number of spcecies draws to be family representatives
trials <- 100       # Draws from the expected codon usage distribution
set.seed(353204)
genSeqs <- FALSE    # If "TRUE" will output repsPerGene nuc alignment files

# Load packages and setwd------------------------------------------------------

library("ape")
library("hash")
library("RGenetics")

# Set the working directory to the directory in which the script is run.
thisDir <- dirname(parent.frame(2)$ofile)
setwd(thisDir)

# Functions -------------------------------------------------------------------

# Reads in a list of filenames with aligned nucleotide sequences 
# (starting and ending at coding positions)
primateCA <- function(geneNames){
  
  
  labels <- c("GC.All", "GC.Var", "Mouse.All", "Mouse.Var", "PerGene.All", 
              "PerGene.Var")
  entropyList <- setNames(replicate(6, matrix(0, nrow=repsPerGene, ncol=4, 
                dimnames=list(1:repsPerGene, 
                              c("Observed", "Expected", "Var", "Length"))), 
                simplify=FALSE),labels)
  
  # Generate Mouse whole-genome codon usage distribution
  mouse <- genMouseUsage()
  
  for (i in 1:length(geneNames)){
    gene <- processFile(paste("../data/aligned_coding_genes/", 
                              geneNames[i], sep=""))
    # Find the GC content of each gene
    gc <- ThirdPosGCContent(gene$dna)
    # Generate GC content per gene codon usage distribution
    codonFreqGC <- gcCodonUsage(gc)
    # Generate per-gene codon usage distribution
    codonFreqGene <- geneCodonUsage(gene)

    #Pick Which taxa to use as family representatives
    taxa <- sampleSpecies(repsPerGene, gene)

    a <- CalculateEntropy(gene, taxa, codonFreqGC, genSeqs)
    b <- CalculateEntropy(gene, taxa, mouse)
    c <- CalculateEntropy(gene, taxa, codonFreqGene)
    
    entropyList[]$GC.All <- entropyList[]$GC.All + a[[1]]
    entropyList[]$GC.Var <- entropyList[]$GC.Var + a[[2]]
    entropyList[]$Mouse.All <- entropyList[]$Mouse.All + b[[1]]
    entropyList[]$Mouse.Var <- entropyList[]$Mouse.Var + b[[2]]
    entropyList[]$PerGene.All <- entropyList[]$PerGene.All + c[[1]]
    entropyList[]$PerGene.Var <- entropyList[]$PerGene.Var + c[[2]]
  }
  return(entropyList)
}

#Calculate entropy and generate null distribution
CalculateEntropy <- function(gene, taxaList, codonFreq, generator = FALSE){

  allSites <- matrix(0, nrow = repsPerGene, ncol = 4, 
                dimnames=list(1:repsPerGene, 
                              c("Observed", "Expected", "Var", "Length")))
  varSites <- matrix(0, nrow = repsPerGene, ncol = 4, 
                dimnames=list(1:repsPerGene, 
                              c("Observed", "Expected", "Var", "Length")))
  
  for (i in 1:repsPerGene){
    # Use all degenerate sites in the alignment
    allSynSites <- AllSynonymousSites(gene, taxaList[i, ])
    # Use only variable sites in the alignment
    varSynSites <- variableSynonymousSites(gene, taxaList[i, ], i, generator)

    if (length(allSynSites) > 0){
      observedEntropyAll <- calculateObservedEntropy(gene, allSynSites, 
                                                     taxaList[i, ])
      expectedEntropyAll <- generateUsage(gene, allSynSites, 
                                          taxaList[i, ], codonFreq)
      allSites[i,1] <- observedEntropyAll
      allSites[i,2] <- mean(expectedEntropyAll)
      allSites[i,3] <- var(expectedEntropyAll)
      allSites[i,4] <- length(allSynSites)
    }

    if (length(varSynSites) > 0){
      observedEntropyVar <- calculateObservedEntropy(gene, varSynSites, 
                                                     taxaList[i, ])
      expectedEntropyVar <- generateUsage(gene, varSynSites, 
                                          taxaList[i,], codonFreq)
      varSites[i,1] <- observedEntropyVar      
      varSites[i,2] <- mean(expectedEntropyVar)
      varSites[i,3] <- var(expectedEntropyVar)
      varSites[i,4] <- length(varSynSites)
    }
  }
  a <- list(allSites, varSites)
  return(a)
}

# Convert sequence of coding dna into amino acid sequence
dnaToAA <- function(x) {
    n <- length(x)
    
    if ( n %% 3 != 0 ) {
        stop("length of sequence not a multiple of three")
    }
    m <- n %/% 3
    aa <- character(m)
  
    for (i in 1:m) {
        foo <- tolower(x[((i*3)-2):(i*3)])

        if (any(foo == '-')){
            aa[i] <- '-'
        } else if (!all(foo %in% c('a','c','g','t'))){
            aa[i] <- '-'
        } else {
          aa[i] <- codonToAAone(paste(x[(i*3)-2], x[(i*3)-1], x[i*3], sep = ""))
        } 
    }
    return(aa)
}

# Remove first word from a string: here, get genus from species name
getGenus <- function(x) {
    return( strsplit(x," ")[[1]][1] )
}

# Get family from genus using the information in primates
getFamily <- function(genus) {
    return( as.character(primates$Family[which(primates$Genus == genus)[1]]) )
}

# Read in a data file;
# Drop last row, which is duplicate of homo sapiens
processFile <- function(file) {
    gene <- read.dna(file,as.character = TRUE, format = "fasta")
    gene <- gene[-nrow(gene),]
    nspecies <- nrow(gene)
    species <- rownames(gene)
    genus <- character(nspecies)
    for ( i in 1:nspecies )
        genus[i] <- getGenus(species[i])

    family <- character(nspecies)
    for ( i in 1:nspecies )
        family[i] <- getFamily(genus[i])
    
    aa <- character(0)
    for ( i in 1:nrow(gene) ) {
        aa <- rbind(aa,dnaToAA(gene[i,]))
    }
    
    return (list(name = sub("^([^.]*).*", "\\1", file), species = species, 
                genus = genus, family = as.factor(family), dna = gene, aa = aa))
}

# calculate entropy statistic from a vector of counts
entropy <- function(x) {
  y <- x[x>0]
  p <- y/sum(y)
  return( -1 * sum(p * log(p) ) )
}

# sample codon usage from a distribution
codonSample <- function(n,p) {
  s <- sample(1:length(p), size = n, prob = p, replace = TRUE)
  x <- rep(0,n)
  for ( i in 1:n )
    x[i] <- sum(s==i)
  return(x)
}

# sample from null distribution
nullSample <- function(B,n,p) {
  out <- numeric(B)
  for ( i in 1:B ) {
    out[i] <- entropy( codonSample(n,p) )
  }
  return(out)
}
  
# Function to sample variable length vectors (deals with length=1)
samplePrimates <- function(x, ...) x[sample(length(x), ...)]

# Sample one species from each family, B times
# Return in a matrix with B rows, each row is a sample of the indicies
# Exclude families with no sequenced representatives in a particular gene
sampleSpecies <- function(B,x) {
    nfamilies <- length(levels(x$family))
    out <- matrix(0,B,nfamilies)
    for ( i in 1:nfamilies ) {
        fam <- (levels(x$family))[i]
        indices <- which( x$family == fam )
        if ( !any(is.na(indices)) )
            out[,i] <- samplePrimates(indices, size = B, replace = TRUE)
    }
    return(out)
}

# Return a vector of indices of synonymous aa sites from a group of taxa 
# Include wholly conserved sites
AllSynonymousSites <- function(x, s) {
  synSites <- vector()
  for (i in 1:length(x$aa[1,])){
    residue <- character(length(s))
    for (k in 1:length(s)){
      residue[k] <- x$aa[s[k],i]
    }

    if (all(residue[1] == residue) && (residue[1] != "W") && 
            (residue[1] != "M") && (residue[1] != "-"))
      synSites <- c(synSites, i)
  }
  return(synSites)
}

# Return a vector of indices of synonymous aa sites from a group of taxa 
# Exclude invariant degenerate sites.
variableSynonymousSites <- function(gene, taxa, iteration, generator) {
  synsites <- vector()
  for (i in 1:length(gene$aa[1,])){
    residue <- character(length(taxa))
    codon <- character(length(taxa))
    for (k in 1:length(taxa)){
      residue[k] <- gene$aa[taxa[k],i]
      startPos <- (i*3)-2
      codon[k] <- toupper(paste(gene$dna[taxa[k], startPos:(startPos+2)], 
                                sep = '', collapse = ''))
    }
    if (all(residue[1] == residue) && (residue[1] != "W") && 
            (residue[1] != "M") && (residue[1] != "-")){
      if (!(all(codon[1] == codon))){
        synsites <- c(synsites, i)
      }
    }
  }
  if (generator == TRUE){
      sequenceOutput(gene, synsites, taxa, iteration)
  }
  return(synsites)
}

# Find GC content of all third codon positions in an alignment given a matrix
# of characters `x`
ThirdPosGCContent <- function(x) {
  nsites <- ncol(x)
  if ( nsites %% 3 != 0 )
    stop("Expected number of sites to be a multiple of 3")
  third <- seq(3,nsites,3)
  dna.tab <- table(x[,third])
  gc <- sum(dna.tab[ which(dimnames(dna.tab)[[1]] %in% c('c','g')) ])
  at <- sum(dna.tab[ which(dimnames(dna.tab)[[1]] %in% c('a','t')) ])
  return( gc / (gc+at) )
}

# Generate mouse genome-wide codon usage frequency hash table
genMouseUsage <- function(){
    codons <- read.table("../data/use_mouse.txt", header = T, as.is = 1)
    codonUsageFreq <- hash()
    for (i in 1:length(codons$codon)){
        codonUsageFreq[[codons$codon[i]]] <- codons$freq[i]
    }
    return(codonUsageFreq)
}

#Calculate codon frequencies based on GC content of each gene
gcCodonUsage <- function(gc){
  a <- keys(codonTable)
  codonUsageFreq <- hash()
  
  for (i in 1:length(a)){
    x <- length(codonTable[[a[i]]])
    for (k in 1:x){
      if (x == 2){
        if (substr(codonTable[[a[i]]][k],3,3) == "G" || 
                (substr(codonTable[[a[i]]][k],3,3) == "C"))
          codonUsageFreq[[codonTable[[a[i]]][k]]] <- gc
        else {
          codonUsageFreq[[codonTable[[a[i]]][k]]] <- (1-gc)
        }
      }
      else if(x == 3){
        if (substr(codonTable[[a[i]]][k],3,3) == "G" || 
                (substr(codonTable[[a[i]]][k],3,3) == "C"))
          codonUsageFreq[[codonTable[[a[i]]][k]]] <- 1-(((1-gc)/3)*4)
        else {
          codonUsageFreq[[codonTable[[a[i]]][k]]] <- ((1-gc)/3)*2
        }
      }
      else if(x == 4){
        if (substr(codonTable[[a[i]]][k],3,3) == "G" || 
                (substr(codonTable[[a[i]]][k],3,3) == "C"))
          codonUsageFreq[[codonTable[[a[i]]][k]]] <- gc/2
        else {
          codonUsageFreq[[codonTable[[a[i]]][k]]] <- (1-gc)/2
        }
      }
      else if(x == 6){
        if (substr(codonTable[[a[i]]][k],3,3) == "G" || 
                (substr(codonTable[[a[i]]][k],3,3) == "C"))
          codonUsageFreq[[codonTable[[a[i]]][k]]] <- (gc/3)
        else {
          codonUsageFreq[[codonTable[[a[i]]][k]]] <- ((1-gc)/3)
        }
      }      
    }
  }
  return(codonUsageFreq)
}

# Calculate codon frequencies based on observed frequencies per gene
geneCodonUsage <- function(gene){
  counts <- hash()
  output <- hash()
  
  for (i in 1:length(codons$codon)){
    output[[codons$codon[i]]] <- 0
    counts[[codons$codon[i]]] <- 0
  }
  for (i in 1:length(gene$aa[1,])){
    startPos <- (i*3)-2
    for (j in 1:length(gene$species)){
      codon <- toupper(paste(gene$dna[j, startPos:(startPos+2)], 
                             sep = '', collapse = ''))
      if (has.key(codon, counts)){
        counts[[codon]] <- counts[[codon]] + 1
      }
    }  
  }
  for (i in 1:length(keys(counts))){
    c <- keys(counts)[i]
    a <- codonToAAone(c)
    if (a == "Stop"){
      a <- "X"
    }
    s <- 0
    for (j in 1:length(codonTable[[a]])){
      s <- s + counts[[codonTable[[a]][j]]]
    }
    output[[c]] <- counts[[c]] / s
  }
  return(output)
}

# Create a hash table of frequency for each codon position
cHash <- function(x){
  a <- hash()
  for (i in 1:length(x)){
    a[[x[i]]] <- 0
  }
  return(a)
}

# Generate the summed entropy from B samples of a gene given a 
# null distribution of codon usage
generateUsage <- function(gene, sites, taxaReps, codonFreq, B=trials){
  m <- matrix(nrow = B, ncol = length(sites))
  for (i in 1:length(sites)){
    aa <- gene$aa[taxaReps[1], sites[i]]
    a <- codonTable[[aa]]
    pdist <- c()
    for (j in 1:length(a)){
      pdist <- c(pdist, codonFreq[[a[j]]])
    }
    m[,i] <- nullSample(B, length(taxaReps), pdist)
  }
  n <- rowSums(m)
  return(n)
}

# Compute the entropy statistics over a vector y of synonymous sites
calculateObservedEntropy <- function(gene, sites, taxaReps){
  entropyList <- c()

  for (i in 1:length(sites)){
    startPos <- (sites[i]*3)-2
    a <- codonToAAone(paste(gene$dna[taxaReps[1], startPos:(startPos+2)], 
                            sep = '', collapse = ''))
    b <- cHash(codonTable[[a]])
    c <- c()
    
    for (j in 1:length(taxaReps)){
      c <- c(c, toupper(paste(gene$dna[taxaReps[j], startPos:(startPos+2)], 
                              sep = '', collapse = '')))
    }
    
    c <- table(c)
    d <- names(c)
    for (k in 1:length(d)){
      b[[d[k]]] <- (c[[d[k]]])
    }
    e <- c()
    for (l in 1:length(d)){
      e <- c(e, b[[d[l]]])
    }
    entropyList <- c(entropyList, entropy(e))    
  }
  return(sum(entropyList))
}

# Take a gene, list of synonymous sites, the iteration, a list of taxa, and a
# large list of matrices for output, and append further degenerate sites to 
# the matrix.
sequenceOutput <- function(gene, sites, taxa, iteration){
    
    alignLengths <- (length(sites) * 3)
    newOutput <- matrix(, 16, alignLengths)
    rownames(newOutput) <- rownames(sequences[[iteration]])
    
    for (i in rownames(newOutput)){
        if (!(i %in% levels(gene$family))){
            newOutput[i, ] <- rep("-", alignLengths)
        }
    }
    for (t in taxa){
        s <- vector()
        for (j in sites){
            startPos <- (j*3)-2
            s <- c(s, c(gene$dna[t, startPos:(startPos+2)]))
        }
        fam <- as.character(gene$family[t])
        newOutput[fam, ] <- s
    }
    sequences[[iteration]] <<- cbind(sequences[[iteration]], newOutput)
}

# Given a large list of n output alignments, where n is repsPerGene, write n 
# files with alignments in fasta format.
# Called only if genSeqs == TRUE
writeOutputToFiles <- function() {
    for (i in 1:repsPerGene){
        fileName <- paste0("../results/seq/degenAlignment_", i, ".fa")
        f <- file(fileName, "w")
        alignment <- apply(sequences[[i]], 1, paste, collapse="")
        for (j in levels(primates$Family)){
            writeLines(paste0(">", j), con=f)
            writeLines(alignment[[j]], con=f)    
        }
        close(f)
    }
}

# Read input files-------------------------------------------------------------
# Primate Family,Genus,Species table.
primates <- read.table("../data/primate-species-16fam.csv", header = T, 
                       sep = ";", quote = "")

# If sequences are to be output, set up the output structure
if (genSeqs){
    sequences <- list()
    for (i in 1:repsPerGene){
        sequences[[i]] <- matrix(, length(levels(primates$Family)), 0)
        dimnames(sequences[[i]]) <- list(levels(primates$Family))
    }
}

# List of aligned fasta files to be included
genes <- scan(alignedGenes, what = '')

# Read codon table and create a hash table with aa,[codons] key,value pairs
codons <- read.table("../data/codon_table.txt", header = T, sep = "", 
                     as.is = c(1,2))
codonTable <- hash()
for (i in 1:length(codons$codon)){
    if (all(has.key(codons$aa[i], codonTable))){
        codonTable[[codons$aa[i]]] <- append(codonTable[[codons$aa[i]]], 
                                             codons$codon[i])
    }
    else {
        codonTable[[codons$aa[i]]] <- codons$codon[i]
    } 
}

# Call main function-----------------------------------------------------------
final <- primateCA(genes)

# Output ----------------------------------------------------------------------
if (genSeqs){
    writeOutputToFiles()
}

# Write observed, mean expected, and variance expected entropy to files
write.table(final[]$GC.All, file="../results/GC.All.txt", quote=FALSE, 
            row.names=FALSE, col.names=TRUE)
write.table(final[]$GC.Var, file = "../results/GC.Var.txt", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
write.table(final[]$Mouse.All, file = "../results/Mouse.All.txt", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
write.table(final[]$Mouse.Var, file = "../results/Mouse.Var.txt", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
write.table(final[]$PerGene.All, file = "../results/PerGene.All.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(final[]$PerGene.Var, file = "../results/PerGene.Var.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
