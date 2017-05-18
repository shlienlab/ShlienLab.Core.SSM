\name{annotate.mutect2.data}
\alias{annotate.mutect2.data}
\title{
Annotate MuTect2 data.
}
\description{
After running MuTect2 on a per chromosome basis and annotating with ANNOVAR, this program adds additional annotation data for downstream analysis.
}
\usage{
annotate.mutect2.data(data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{MuTect and ANNOVAR processed variant dataframe [required]}
}
\value{
Returns the indel dataframe with appended columns including: hgnc_gene, ensembl_gene, cosmic_census and aa (amino acid).
}
\author{
Richard J. de Borja <richard.deborja@sickkids.ca>
}
\keyword{ indel }
\keyword{ gatk }
