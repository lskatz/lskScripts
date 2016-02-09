#!/bin/env Rscript
# Calculates the Kendall metric on 2+ trees

# Credit to Jo Williams for showing me this metric and getting me going
# in Rscript

#options(error=traceback) 

# Libraries
library(stringr,quietly=TRUE)
library(treescape,quietly=TRUE)
library(phangorn,quietly=TRUE)
library(getopt,quietly=TRUE)

# Script options
# Which values of lambda to calculate with? Values can be 0 to 1.
lambdas <- c(0, 1)

# Import data
# Column  3: Argument mask  of  the flag , an integer
# Possible  values:  0=no argument, 1=required argument, 2=optional argument.
#spec = matrix(c(
#  'help'   , 'h', 0, "logical",   'A call for help',
#  'tree1'  , '1', 1, "character", 'The first tree',
#  'tree2'  , '2', 1, "character", 'The second tree'
#), byrow=TRUE, ncol=5);
#opt = getopt(spec);

## Get a list of all Newick files in the directory
#treefiles <- list.files(path = "LeeTesting", 
#                        all.files = FALSE,
#                        full.names = FALSE)

treefiles <- commandArgs(trailingOnly=TRUE);

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
for(t in 1:length(mytrees[[2]]$tip.label)) {
  mytrees[[2]]$tip.label[t] <- str_split(string = mytrees[[2]]$tip.label[t], pattern = "_")[[1]][1]
}
mytrees <- .compressTipLabel(mytrees)

# Headers
cat(paste("Tree1","Tree2","lambda","Kendall","\n",sep="\t"))

# Calculating the Kendall pairwise distance
for(t in 1:(length(mytrees)-1)){

  for(u in (t+1):length(mytrees)){

    # List of tree files from which to calculate
    treeVector=c(mytrees[t],mytrees[u]);

    # Calculate Kendall metric
    for(x in 1:length(lambdas)){
      dist=multiDist(treeVector, lambda = lambdas[x])
      cat(paste(treefiles[t],treefiles[u],lambdas[x],dist,"\n", sep="\t"));
    }

  }
}
