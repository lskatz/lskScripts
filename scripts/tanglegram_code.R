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

print("Loading libraries")
#library(phytools,  quietly=TRUE)
#library(dendextend,quietly=TRUE)
#library(ape,       quietly=TRUE)
myReturn <- suppressPackageStartupMessages(c(
  library(phytools,  quietly=TRUE),
  library(dendextend,quietly=TRUE),
  library(ape,       quietly=TRUE)
));

#eid_iso=read.delim('\\\\cdc.gov\\project\\CCID_NCZVED_DFBMD_EDEB\\Analytics\\Weidong\\LM model\\ncbi_access_dat.csv',sep=',',header=T,stringsAsFactors = F)
outbreak <- read.delim('/scicomp/home/gzu2/projects/mashtree/data/katzEtAl/Lyve-SET/outbreakStatus.tsv',
                       sep="\t", header=T, stringsAsFactors=F)
tree1 <- reorder(midpoint.root(read.tree(treefile1)), order = "cladewise")
tree2 <- reorder(midpoint.root(read.tree(treefile2)), order = "cladewise")

dend_tree1 <- force.ultrametric(tree1)
dend_tree2 <- force.ultrametric(tree2)

min_length <- 0.000000000000000000001
dend_tree1$edge.length[ dend_tree1$edge.length < min_length ] <- min_length
dend_tree2$edge.length[ dend_tree2$edge.length < min_length ] <- min_length

dend_tree1=(midpoint.root(dend_tree1))
dend_tree2=(midpoint.root(dend_tree2))

my.col=c('blue','brown','green','pink','red')

outbreakIndex    <- match(outbreak$sample[outbreak$outbreak== 1],dend_tree1$tip.label)
nonoutbreakIndex <- match(outbreak$sample[outbreak$outbreak== 0],dend_tree1$tip.label)
maybeoutbreakIndex <- match(outbreak$sample[outbreak$outbreak==-1],dend_tree1$tip.label)
myColors <- c()
myColors[outbreakIndex]      <- 'red'
myColors[nonoutbreakIndex]   <- 'blue'
myColors[maybeoutbreakIndex] <- 'green'
#outbreak$color[ outbreakIndex ] <- 'red'
#outbreak$color[nonoutbreakIndex]<- 'blue'
#col.s=my.col[as.factor(conn_l_col)]

print("untangle")
dendl <- dendextend::untangle(as.dendrogram(dend_tree1), 
                     as.dendrogram(dend_tree2), 
                     method = "step2side") 
                     #method = "labels") 
                     #method = "ladderize") 
                     #method = "random") 
                     #method = "step1side") 
                     #method = "DendSer") 

# Make the branches look nice
dendl %>% set("branches_lwd", 1) %>%
  set("labels_col", "white") -> dendl

print("entanglement...");
myEntanglement <- entanglement(dendl)
cophenetic <- cor.dendlist(dendl, method = "cophenetic")
baker      <- cor.dendlist(dendl, method = "baker")

# Start off the viz
png(outfile)
tanglegram(dendl,
           main_left='Lyve-SET',
           main_right='Mashtree',
           lab.cex=0.3,
           highlight_distinct_edges = FALSE,
           color_lines=myColors
           )
#myReturn <- text("SOMETHING", x=1, y=1)
myReturn <- dev.off();



