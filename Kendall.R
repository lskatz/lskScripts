#!/usr/bin/env Rscript
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
reps    <- 1000    # number of repetitions for background


########################
# START Functions
kendallBackground <- function(treeObj, lambdaCoefficient, rep=1000){
  # Figure out taxa
  taxa=treeObj$tip.label
  numTaxa=length(taxa)
  # Figure out max and min branch lengths for br=
  branchLength=sort(treeObj$edge.length)
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
## END Functions
####################


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
opt=c()

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
#for(t in 1:length(mytrees[[2]]$tip.label)) {
#  mytrees[[2]]$tip.label[t] <- str_split(string = mytrees[[2]]$tip.label[t], pattern = "_")[[1]][1]
#}

mytrees <- .compressTipLabel(mytrees)

# Headers
header=c("Tree1","Tree2","lambda","Kendall");
#if(match('backgroundDistribution', opt) > 0){
  header=append(header,c("BackgroundKendall","n","Z","p-value"));
#}
cat(paste(header,sep="\t"),"\n")

# Calculating the Kendall pairwise distance
for(t in 1:(length(mytrees)-1)){
  
  # Get the background of Kendall distributions with lambda==0 and lambda==1
  # The key will be 'lambda0', 'lambda0.5', etc
  background=list()
  for(x in 1:length(lambdas)){
    key <- paste("lambda",lambdas[x],sep="")
    background[[key]]=kendallBackground(mytrees[[t]],lambdas[x],rep=reps)
  }

  for(u in (t+1):length(mytrees)){

    # List of tree files from which to calculate
    treeVector=c(mytrees[t],mytrees[u]);

    # Calculate Kendall metric
    for(x in 1:length(lambdas)){
      # background distribution
      lambdakey=paste("lambda",lambdas[x],sep="")
      # Convert list to vector, and then call 'mean'
      backgroundMean=mean(unlist(background[[lambdakey]]))
      backgroundSd = sd(unlist(background[[lambdakey]]))

      # Get the actual distance
      dist=multiDist(treeVector, lambda = lambdas[x])

      # Figure out the stats
      z <- (backgroundMean - dist)/backgroundSd
      pvalue <- 2 * pnorm(z)

      # Generate output
      distributionString=paste(round(backgroundMean,digits=2),"±",round(backgroundSd,digits=2),sep="")
      pvalueString=round(pvalue,digits=6)
      zString=round(z,digits=2)
      rowVector=c(treefiles[t],treefiles[u],lambdas[x],round(dist,digits=2),distributionString,reps,zString,pvalueString);
      cat(paste(rowVector,sep="\t"),"\n");
    }

  }
}

