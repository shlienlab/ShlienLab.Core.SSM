\name{create_final_filtered_dataframe}
\alias{create_final_filtered_dataframe}
\title{
Create the final filtered dataframe.
}
\description{
Creates the hard filtered final dataframe by removing PON, cFilter, soft clipped and hard clipped
filtered data from the input dataframe.
}
\usage{
create_final_filtered_dataframe(data = NULL, pon_threshold = 1, flagged_threshold = 1, soft_threshold = 0, hard_threshold = 0)
}
\arguments{
  \item{data}{dataframe containing the annotated SNV data (required)}
  \item{pon_threshold}{the number of recurrent normal samples a mutation is found in (default: 1)}
  \item{flagged_threshold}{the number of flagged filters a mutation is found in (default: 1)}
  \item{soft_threshold}{the maximum number of reads a mutation can be found in a soft clipped region (default: 0)}
  \item{hard_threshold}{the maximum number of reads a mutation can be found in a hard clipped region (default: 0)}
}
\value{
Returns a dataframe containing hard filtered SNVs.
}
\author{
Richard J. de Borja <richard.deborja@sickkids.ca>
}
\keyword{ snvs }
\keyword{ somatic }

