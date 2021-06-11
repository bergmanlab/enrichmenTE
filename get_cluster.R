#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)#!/usr/bin/env Rscript

matrix_file <- args[1]
tree_file <- args[2]
pdf_file <- args[3]

# load packages
library(tidyverse)
library(phangorn)

# read data
data <-
  read_csv(matrix_file) %>%
  column_to_rownames("sample") %>%
  as.matrix()

# distance matrix
d <- dist(data, method = "euclidean")

# build NJ tree
NJ <- NJ(d)
tree <- phangorn::midpoint(NJ)
tree$tip.label <- str_remove_all(tree$tip.label, "GSF2595-")

# save tree as nwk file
write.tree(
  tree,
  file = tree_file,
  append = FALSE,
  digits = 10,
  tree.names = FALSE
)

# save tree as pdf file
pdf(
  file = pdf_file,
  width = 8,
  height = 8
)
plot(tree, use.edge.length=T,underscore=T,no.margin = T,font=1)
dev.off()
