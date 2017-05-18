annotate.mutect2.data <- function(data=NULL) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  data$ensembl_gene <- apply(X=as.matrix(data$annovar_ens_gene), MARGIN=1, FUN=remove.multiple.genes)
  data$hgnc_gene <- apply(X=as.matrix(data$annovar_gene), MARGIN=1, FUN=remove.multiple.genes)
  data$cosmic_census <- filter.for.cosmic.cancer.gene.census(data=data$hgnc_gene)
  data$aa <- apply(X=as.matrix(data$annovar_annotation), MARGIN=1, FUN=add.aminoacid.column)
  return(data)
  }
