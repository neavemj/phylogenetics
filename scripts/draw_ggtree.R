# script to draw ggtree in the phylogenetics module
# Matthew Neave 2021-08-04

# required libraries - need to add to conda environment
library(ggtree)
library(ape)
library(ggplot2)
library(phangorn)

# the below files are specified as inputs in the snakemake rule
# if we do it this way, need to ensure filenames exactly match

tree = read.tree("04_phylogenetics/HA-sequences.treefile")
metadata = read.table("04_phylogenetics/HA-sequences.metadata.txt")

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
ggsave("04_phylogenetics/HA-sequences.tree.pdf", width=10, height=10, limits=F)
