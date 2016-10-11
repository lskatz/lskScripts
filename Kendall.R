#!/usr/bin/env Rscript

library(docopt,quietly=TRUE)

logmsg <- function(msg){
  cat(paste(c("Kendall.R:",msg),sep=""),"\n", file=stderr())
}

# Get pristine ARGV
#argv <- commandArgs(trailingOnly=FALSE);
#thisScript = substring(argv[grep("--file=", argv)], 8)

# Get command line parameters
doc <- "Description: This script uses the Kendall-Colijn phylogeny metric to see if one or more query trees are close to a reference tree.  See: Kendall and Colijn 2015, Arxiv

Usage: Kendall.R [--lambda=f] [options] REFERENCE_TREE QUERY_TREE... 

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
    --rootnode=<i>   Perform a root at this node. 0 for midpoint root. [Default: 0]

  Examples:
    Kendall.R --rep 1000 trees/*.dnd | column -t
    Kendall.R trees/*.dnd | sort -k8,8n | column -t
    Kendall.R --lambda 0 --lambda 1 trees/*.dnd

"

# Script options
opts <- docopt(doc)

# Libraries
# These libraries are still noisy despite the quietly parameter
suppressMessages(library(stringr,quietly=TRUE))
suppressMessages(library(treescape,quietly=TRUE))
suppressMessages(library(phangorn,quietly=TRUE))
suppressMessages(library(ips,quietly=TRUE))
suppressMessages(library(ggplot2,quietly=TRUE))
suppressMessages(library(tools,quietly=TRUE))

treefiles <- c(opts$REFERENCE_TREE,opts$QUERY_TREE)
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

# The background distribution helps answer whether a query tree is
# close enough to a reference tree.  In the Lyve-SET paper, we would
# be asking for example, if kSNP3's tree (query) is close to the 
# Lyve-SET (reference) tree
# In particular, this distribution is the set of Kendall scores of
# the query tree vs random trees
#                          e.g., kSNP    e.g., Lyve-SET
kendallBackground <- function(queryTree, lambdaCoefficient, rep=1000){
  # Figure out taxa
  taxa=queryTree$tip.label
  numTaxa=length(taxa)
  # Figure out max and min branch lengths for br=
  branchLength=sort(queryTree$edge.length)
  minLength=min(branchLength)
  maxLength=max(branchLength)

  # Have to make the multiphylo object.
  # The first element is the original tree.
  # The second element will be the random tree.
  # Convert to multiphylo for Kendall metric
  mytrees <- vector("list",2)
  class(mytrees) <- "multiphylo"
  mytrees[[1]]=queryTree

  # Start recording the kendal metric in a vector.
  # It can be averaged out later.
  kendallVec=c()
  for(i in 1:rep){
    # Make a distribution of branch lengths 
    # whose min and max are defined above
    randBranchLength <- runif(length(branchLength), min=minLength, max=maxLength)

    # Generate a random tree with random taxa and branch lengths
    #mytrees[[2]] <- rtree(numTaxa,rooted=TRUE, tip.label=taxa, randBranchLength)
    mytrees[[2]] <- rcoal(numTaxa, rooted=TRUE, tip.label=taxa, br=randBranchLength)
    mytrees <- .compressTipLabel(mytrees)

    # Find the Kendall metric between this random tree and the 
    # true tree.
    kendall <- multiDist(mytrees, lambda = lambdaCoefficient);
    kendallVec=append(kendallVec,kendall)
  }
  kendallVec  # returns the kendall values
}

plotBackground <- function(distribution, observed){
  my_example <- rnorm(n=1000, m=24.2, sd=2.2)

  dat <- data.frame(replicate=c(1:length(distribution)), Kendall=distribution)

  binwidth <- (max(distribution)-min(distribution))/100
  my_histogram <- ggplot(data=dat, aes(x = Kendall)) +
                  geom_histogram(binwidth=binwidth) +
                  geom_vline(xintercept = observed)
  return(my_histogram)
}

readyTreeForComparison <- function(treefile,root.node=0){
  # Read tree
  my_tmptree <- read.tree(file = treefile)

  # Add 100% confidence where it is null.  However, is it on a 0-1
  # or a 0-100 scale first?
  node_labels <- as.numeric(my_tmptree$node.label)
  if(all(is.na(node_labels))){
    node_labels=rep(100,my_tmptree$Nnode)
  }
  max_bootstrap <- max(node_labels,na.rm=TRUE)
  node_labels[is.na(node_labels)] <- max_bootstrap
    
  node_labels=mapply(sprintf,"%0.2f",node_labels)
  my_tmptree$node.label <- as.character(node_labels)
  # collapse low-confidence nodes
  my_tmptree <- collapseUnsupportedEdges(my_tmptree, "node.label", 0.7*max_bootstrap)
  # For whatever reason, collapseUnsupportedEdges adds NA values to the end of the node.label vector.
  # Removing all NAs should be fine because the preexisting NAs have already been converted.
  #my_tmptree$node.label <- my_tmptree$node.label[!is.na(my_tmptree$node.label)]
  # midpoint root
  if(root.node==0){
    my_tmptree <- midpoint(my_tmptree)
  } else {
    # Root it according to the user's parameter. 
    # A nuance is that the root length has to be zero to 
    # be considered rooted according to the documentation.
    clade <- strsplit(root.node,',')[[1]]
    if(!is.monophyletic(my_tmptree,clade)){
      stop("ERROR: the outgroup",paste(clade,sep=",")," is not monophyletic on ",treefile)
    }
    my_tmptree <- root(my_tmptree,clade,resolve.root=TRUE)
    my_tmptree$root.edge <- 0
  }
  # Sort polytomies
  my_tmptree <- reorder(my_tmptree, order="cladewise", index.only=FALSE)
  
  # Make into a list, to make it compatible with Kendall-Colijn
  return(list(my_tmptree))
  
}

# I took the latest code from IPS so that the node
# labels were not broken
# https://github.com/heibl/ips/blob/248432094daf39da0cda373a942de41cf54ff99c/R/collapseUnsupportedEdges.R
collapseUnsupportedEdges <- function(phy, value, cutoff){
  
  if ( !inherits(phy, "phylo") ) 
    stop ("'phy' is not of class 'phylo'")
  
  if ( missing(value) ) value <- "node.label"
  
  stat <- as.numeric(phy[[value]])
  nt <- Ntip(phy)
  root.node <- nt + 1
  collapse <- which(stat < cutoff) + nt
  ## the root node cannot be collapsed:
  collapse <- setdiff(collapse, root.node)
  
  ## collapse nodes in post-order traversal!!
  ## ----------------------------------------
  for ( i in rev(collapse) ){
    
    #i <- rev(collapse)[1] # FOR DEBUGGING
    
    ## identify edges
    id <- phy$edge[, 2] == i
    id2 <- phy$edge[, 1] == i
    
    ## modify: node values
    #phy$node.label <- phy$node.label[!id]
    
    ## modify: edges
    phy$edge[id2, 1] <- phy$edge[id, 1]
    phy$edge <- phy$edge[!id, ]
    phy$edge[phy$edge > i] <- phy$edge[phy$edge > i] - 1
    
    ## modify: edge lengths
    phy$edge.length[id2] <- phy$edge.length[id2] + phy$edge.length[id]
    phy$edge.length <- phy$edge.length[!id]
    
    ## modify: node labels
    phy[["node.label"]] <- phy[["node.label"]][-(i - nt)]
    
    ## modify: number of internal nodes
    phy$Nnode <- phy$Nnode - 1
    
  }
  phy
}

## END Functions
####################

#treefiles <- c(opt$tree1, opt$tree2)
ntrees <- length(treefiles)

## Loop over files and import as multiphylo object
mytrees <- vector("list", ntrees)
class(mytrees) <- "multiphylo"
for(f in 1:ntrees) {
  basename <- file_path_sans_ext(basename(treefiles[f]))

  ### This is actually reading the tree file, doing a midpoint root, 
  ### and storing it as a list in one slot of the vector.
  logmsg(c("Reading",treefiles[f]));

  #make it compatible with Kendall-Colijn
  mytrees[f] <- readyTreeForComparison(treefiles[f],opts$rootnode)

  write.tree(mytrees[f][[1]], file=paste(basename,".flattened.dnd",sep=""))
}

mytrees <- .compressTipLabel(mytrees)



# Print headers
header=c("Reference","Query","lambda","Kendall");
if(opts$background){
  header=append(header,c("BackgroundKendall","n","Z","p-value"));
}
cat(paste(header,sep="\t"),"\n",sep="\t")

# Calculating the Kendall pairwise distance
histogram=c() # saving histogram plots in case I want them later
t <- 1
logmsg(c("Kendal distances with",treefiles[t],"as a reference"))
for(u in 1:length(mytrees)){

  # Check if taxa are comparable
  if(!identical(sort(mytrees[[t]]$tip.label) , sort(mytrees[[u]]$tip.label))){
    stop(paste(c("ERROR: tips of ",treefiles[t]," and ",treefiles[u]," are different")));
  }

  logmsg(c("Query:",treefiles[u]))

  queryTree <- mytrees[[u]]
  referenceTree <- mytrees[[t]]
  # List of tree files from which to calculate
  treeVector <- c(queryTree, referenceTree);

  # Calculate Kendall metric of the query vs reference
  dist <- multiDist(treeVector, lambda = lambda)
  dist <- dist[1]

  # Generate output
  rowVector=c(treefiles[t],treefiles[u],lambda,round(dist,digits=2))

  # Get the background of Kendall distributions
  if(opts$background){
    background=kendallBackground(queryTree, lambda,rep=reps)
    backgroundMean=mean(background)
    backgroundSd = sd(background)

    # Calculate Z and P
    z <- (dist - backgroundMean)/backgroundSd # Want to know whether the background is bigger than observed
    pvalue <- pnorm(z)

    # Formatting for output
    distributionString=paste(round(backgroundMean,digits=2),"±",round(backgroundSd,digits=2),sep="")
    pvalueString=round(pvalue,digits=4)
    zString=round(z,digits=2)

    # Add onto the output vector
    rowVector=append(rowVector,c(distributionString,reps,zString,pvalueString))
  }

  # Print the output
  cat(paste(rowVector,sep="\t"),"\n",sep="\t");

  if(opts$plot){
    outfile=paste("Kendall.",u,".bmp",sep="")
    my_histogram=append(histogram,plotBackground(background,dist))
    logmsg(c("Printing to file",outfile))
    suppressMessages(
      ggsave(filename=outfile)
    )
  }
}

