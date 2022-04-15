MEanalysis <- function(somMuts, pathways){
  install.packages("Rediscover", repos='http://cran.us.r-project.org')
  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repos='http://cran.us.r-project.org')
  }
  BiocManager::install("maftools", repos='http://cran.us.r-project.org')
  library("Rediscover")
  library("tidyverse")
  library("discover")
  out = data.frame(names=character(), p_value=double())
  if (nrow(pathways) > 0){
    for (row in 1:nrow(pathways)) {
      if (nchar(pathways[row, "Somatically_Mutated_Genes"])>2){
        genes <-as.vector(strsplit(substring(pathways[row, "Somatically_Mutated_Genes"], 2,nchar(pathways[row, "Somatically_Mutated_Genes"])-1), ", "))[[1]]
        A = data.matrix(somMuts)
        PM <- getPM(A)
        x= data.frame(rownames(somMuts))
        x$indices = rownames(x)
        p_val = getMutexGroup(data.matrix(A[genes,]), PM[as.numeric(c(filter(x, rownames.somMuts. %in% genes)$indices)),], "Coverage")
        out[nrow(out) + 1,] = c(pathways[row, "Pathway_Name"], p_val)
      }
    }
  }
  return(out)
}