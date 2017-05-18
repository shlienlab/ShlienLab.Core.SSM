\name{get.mutect2.annotated.header}
\alias{get.mutect2.annotated.header}
\title{
Get the MuTect2 annotated header.
}
\description{
A method that returns a vector containing the header for the MuTect2 output data.
}
\usage{
get.mutect2.annotated.header()
}
\value{
Returns a vector containing the column header names for the annotated dataframe.  The
original header in the output file from ANNOVAR contains incorrect information and
is updated to reflect all columns in the dataframe.
}
\author{
Richard J. de Borja <richard.deborja@sickkids.ca>
}
\keyword{ ssm }
\keyword{ header }
\keyword{ mutect2 }
\keyword{ gatk }
