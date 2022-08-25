
library(XML)

get_genes_from_kegg_xml <- function(kegg_file_loc){
  # parse the document
  source <- readLines(kegg_file_loc, encoding = "UTF-8")
  parsed.doc <- htmlParse(source, encoding = "UTF-8")
  # grab the nodes
  entries <- xpathSApply(parsed.doc, path='//pathway/entry[@type="gene"]/graphics/@name')
  # change to a vector
  gene_strings <- as.vector(unlist(entries))
  # remove spaces
  gene_strings <- gsub(' ', '', gene_strings)
  # split each string by the comma
  gene_lists <- strsplit(gene_strings, ',')
  # turn the list of vectors into one vector
  genes <- Reduce(c, gene_lists)
  # remove duplicates
  genes <- unique(genes)
  return(genes)
}
genes <- get_genes_from_kegg_xml('~/Downloads/hsa04657.xml')
