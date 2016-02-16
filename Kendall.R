#!/usr/bin/env Rscript
# Calculates the Kendall metric on 2+ trees

# Credit to Jo Williams for showing me this metric and getting me going
# in Rscript

#options(error=traceback) 

# Libraries
library(stringr,quietly=TRUE)
library(treescape,quietly=TRUE)
library(phangorn,quietly=TRUE)
library(docopt,quietly=TRUE)

# Get pristine ARGV
#argv <- commandArgs(trailingOnly=FALSE);
#thisScript = substring(argv[grep("--file=", argv)], 8)

# Get command line parameters
doc <- "Description: This script uses the Kendall-Colijn phylogeny metric to determine the distance between two rooted trees.  See: Kendall and Colijn 2015, Arxiv

Usage: Kendall.R [--lambda=f] [options] TREE TREE... 

All results are printed to stdout.

  Options:
    -h --help        Show this screen
    --lambda=<f>     A lambda value to use in the metric. [Default: 0.5]
                     0 gives weight to a topology metric; 1 gives weight to branch lengths.
                     Must be between 0 and 1.
    --background     Create a background distribution and generate a p-value
    --rep=<i>        Number of replicates in the background distribution [Default: 1000]
                     Implies --background
    --seed=<i>       A seed for randomly generating a background distribution of trees [Default: -1]
                     Implies --background
    --plot           Generate a plot of the background distribution of Kendall vs observed value.
                     Implies --background

  Examples:
    Kendall.R --rep 1000 trees/*.dnd | column -t
    Kendall.R trees/*.dnd | sort -k8,8n | column -t
    Kendall.R --lambda 0 --lambda 1 trees/*.dnd

"

# Script options
opts <- docopt(doc)
treefiles <- opts$TREE
# Which values of lambda to calculate with? Values can be 0 to 1.
lambda <- as.double(opts$lambda)
reps    <- opts$rep
if(opts$seed >= 0){
  set.seed(opts$seed)
}
# Set background distribution if anything needs it
if(opts$plot | opts$seed > 0 | opts$rep > 0){
  opts$background=TRUE;
}

# Parameter checking
if(lambda < 0 | lambda > 1){
  stop("ERROR: lambda must be between 0 and 1")
}
if(reps < 1){
  stop("ERROR: number of reps must be > 0")
}

########################
# START Functions
logmsg <- function(msg){
  cat(paste(c("Kendall.R:",msg),sep=""),"\n", file=stderr())
}

kendallBackground <- function(treeObj, treeObj2, lambdaCoefficient, rep=1000){
  # Figure out taxa
  taxa=treeObj$tip.label
  numTaxa=length(taxa)
  # Figure out max and min branch lengths for br=
  branchLength=sort(c(treeObj$edge.length,treeObj2$edge.length))
  minLength=branchLength[1]
  maxLength=branchLength[length(branchLength)]

  # Have to make the multiphylo object.
  # The first element is the original tree.
  # The second element will be the random tree.
  # Convert to multiphylo for Kendall metric
  mytrees <- vector("list",2)
  class(mytrees) <- "multiphylo"
  mytrees[[1]]=treeObj

  # make a function for choosing random branch lengths for the 
  # rtree() function
  randBranchLength <- function(v){
    return(v[round(runif(1,1,length(v)))])
  }

  # Start recording the kendal metric in a vector.
  # It can be averaged out later.
  kendallVec=c()
  for(i in 1:rep){
    mytrees[[2]] <- rtree(numTaxa,rooted=TRUE, tip.label=taxa, br=randBranchLength(branchLength))
    mytrees <- .compressTipLabel(mytrees)

    # Find the Kendall metric between this random tree and the 
    # true tree.
    kendall <- multiDist(mytrees, lambda = lambdaCoefficient);
    kendallVec=append(kendallVec,kendall)
  }
  kendallVec
  #return(kendallVec)
}

plotBackground <- function(distribution, observed){
  library(ggplot2,quietly=TRUE)

  my_example <- rnorm(n=1000, m=24.2, sd=2.2)

  dat <- data.frame(replicate=c(1:length(distribution)), Kendall=distribution)

  binwidth <- (max(distribution)-min(distribution))/100
  my_histogram <- ggplot(data=dat, aes(x = Kendall)) +
                  geom_histogram(binwidth=binwidth) +
                  geom_vline(xintercept = observed)
  return(my_histogram)
}
## END Functions
####################

#treefiles <- c(opt$tree1, opt$tree2)
ntrees <- length(treefiles)

## Loop over files and import as multiphylo object
mytrees <- vector("list", ntrees)
class(mytrees) <- "multiphylo"
for(f in 1:ntrees) {
  ### This is actually reading the tree file, doing a midpoint root, 
  ### and storing it as a list in one slot of the vector.
  mytrees[f] <- list(midpoint(read.tree(file = treefiles[f])))
}

## Cleaning up the tip labels on the second tree because they don't match the first. 
## Idiosyncratic.
#for(t in 1:length(mytrees[[2]]$tip.label)) {
#  mytrees[[2]]$tip.label[t] <- str_split(string = mytrees[[2]]$tip.label[t], pattern = "_")[[1]][1]
#}

mytrees <- .compressTipLabel(mytrees)

# Headers
header=c("Tree1","Tree2","lambda","Kendall");
if(opts$background){
  header=append(header,c("BackgroundKendall","n","Z","p-value"));
}
cat(paste(header,sep="\t"),"\n")

# Calculating the Kendall pairwise distance
histogram=c() # saving histogram plots in case I want them later
for(t in 1:(length(mytrees)-1)){
  logmsg(c("Kendal distances for",treefiles[t]))
  
  for(u in (t+1):length(mytrees)){

    # Get the background of Kendall distributions
    if(opts$background){
      background=kendallBackground(mytrees[[t]],mytrees[[u]],lambda,rep=reps)
      backgroundMean=mean(background)
      backgroundSd = sd(background)
    }

    # List of tree files from which to calculate
    treeVector=c(mytrees[t],mytrees[u]);

    # Calculate Kendall metric
    dist=multiDist(treeVector, lambda = lambda)
    dist=dist[1]

    # Generate output
    rowVector=c(treefiles[t],treefiles[u],lambda,round(dist,digits=2))

    # background distribution
    if(opts$background){
      # Calculate Z and P
      z <- (backgroundMean - dist)/backgroundSd # Want to know whether the background is bigger than observed
      pvalue <- pnorm(z)

      # Formatting for output
      distributionString=paste(round(backgroundMean,digits=2),"±",round(backgroundSd,digits=2),sep="")
      pvalueString=round(pvalue,digits=6)
      zString=round(z,digits=2)

      # Add onto the output vector
      rowVector=append(rowVector,c(distributionString,reps,zString,pvalueString))
      
    }

    # Print the output
    cat(paste(rowVector,sep="\t"),"\n");

    if(opts$plot){
      outfile=paste(t,"_",u,".bmp",sep="")
      my_histogram=append(histogram,plotBackground(background,dist))
      logmsg(c("Printing to file",outfile))
      suppressMessages(
        ggsave(filename=outfile)
      )
    }
  }
}
#print(histogram)

