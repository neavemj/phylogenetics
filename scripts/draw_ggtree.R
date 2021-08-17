# script to draw ggtree in the phylogenetics module
# Matthew Neave 2021-08-04

# required libraries - need to add to conda environment
library(ggtree)
library(ape)
library(ggplot2)
library(phangorn)

library("optparse")

# use optparse to grab command line arguments for the tree

option_list <- list(
  # required args
  make_option(c("--tree"), type="character", default=NULL,
              help="tree file in newick format", metavar="character"),
  make_option(c("--metadata"), type="character", default=NULL,
              help="metadata information for tree tips", metavar="character"),
  make_option(c("--output"), type="character", default=NULL,
              help="name for the rendered output tree file", metavar="character")
)


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#print(opt$tree)
#print(opt$metadata)

# now read data for plotting the tree

tree = read.tree(opt$tree)
metadata = read.table(opt$metadata)

colnames(metadata) <- c("Sample", "Source")

# root tree at the mid-point using phangorn

mid_root <- midpoint(tree)

# add color for source
# this information comes from the metadata rule

cols <- c(
  "NEW" = "#2ecc71",
  "BACKBONE" = "#bdc3c7"
)

# now draw rooted tree using ggplot

p <- ggtree(mid_root) %<+% metadata +
  theme_tree2() +
  geom_tiplab(hjust=-0.018) +
  geom_tippoint(aes(fill=Source), shape=21, size=3, color="black") +
  scale_fill_manual(values = cols, na.value="grey30") +
  xlim(values=c(0, 1.0))

# use ggsave for the plot - required as an output in snakemake
ggsave(opt$output, width=10, height=10, limits=F)
