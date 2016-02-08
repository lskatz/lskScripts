#!/bin/env Rscript
# Calculates the Kendall metric on 2+ trees

# Credit to Jo Williams for showing me this metric and getting me going
# in Rscript

# Libraries
library(stringr)
library(treescape)
library(phangorn)
library(getopt)

# Import data
# Column  3: Argument mask  of  the flag , an integer
# Possible  values:  0=no argument, 1=required argument, 2=optional argument.
spec = matrix(c(
  'help'   , 'h', 0, "logical",   'A call for help',
  'tree1'  , '1', 1, "character", 'The first tree',
  'tree2'  , '2', 1, "character", 'The second tree'
), byrow=TRUE, ncol=5);
opt = getopt(spec);


## Get a list of all Newick files in the directory
#treefiles <- list.files(path = "LeeTesting", 
#                        all.files = FALSE,
#                        full.names = FALSE)
treefiles <- c(opt$tree1, opt$tree2)
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

# Calculating the Kendall pairwise distance
lambdas <- c(0, 0.5, 1.0)
for(x in 1:length(lambdas)){
  dist=multiDist(mytrees, lambda = lambdas[x])
  print(paste(c(lambdas[x],dist), sep="\t"));
}
