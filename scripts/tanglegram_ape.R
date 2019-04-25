#!/usr/bin/env Rscript
# Authors: Beau Bruce and Weidong Gu
# Modified by Lee Katz

library(argparse, quietly=TRUE)
parser <- ArgumentParser()
parser$add_argument("-t1", "--tree1", default=FALSE,
                        help="First tree")
parser$add_argument("-t2", "--tree2", default=FALSE,
                        help="Second tree")
parser$add_argument("-o",  "--outfile", default=FALSE,
                        help="output png file")
args <- parser$parse_args()

treefile1 <- args$tree1
treefile2 <- args$tree2
outfile   <- args$outfile

#print("Loading libraries")
#library(phytools,  quietly=TRUE)
#library(dendextend,quietly=TRUE)
#library(ape,       quietly=TRUE)
myReturn <- suppressPackageStartupMessages(c(
  library(phytools,  quietly=TRUE),
  #library(dendextend,quietly=TRUE),
  library(ape,       quietly=TRUE)
));

outbreak <- read.delim('/scicomp/home/gzu2/projects/mashtree/data/katzEtAl/Lyve-SET/outbreakStatus.tsv',
                       sep="\t", header=T, stringsAsFactors=F)
tree1 <- ladderize(midpoint.root(read.tree(treefile1)))
tree2 <- ladderize(midpoint.root(read.tree(treefile2)))

tree1 <- reorder(tree1, "postorder")
tree2 <- reorder(tree2, "postorder")

# Default minimum length
min_length <- 0.000000000000000000001
tree1$edge.length[ tree1$edge.length < min_length ] <- min_length
tree2$edge.length[ tree2$edge.length < min_length ] <- min_length

outbreakIndex    <- match(outbreak$sample[outbreak$outbreak== 1],tree1$tip.label)
nonoutbreakIndex <- match(outbreak$sample[outbreak$outbreak== 0],tree1$tip.label)
maybeoutbreakIndex <- match(outbreak$sample[outbreak$outbreak==-1],tree1$tip.label)

myColors <- c()
myColors[outbreakIndex]      <- 'red'
myColors[nonoutbreakIndex]   <- 'blue'
myColors[maybeoutbreakIndex] <- 'green'

association <- cbind(tree1$tip.label, tree1$tip.label)
png(outfile)
cophyloplot(tree1, tree2, assoc = association, space = 100, length.line=0, gap=1, show.tip.label=F, col = myColors)
dev.off()

