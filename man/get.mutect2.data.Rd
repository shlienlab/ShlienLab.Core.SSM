\name{get.mutect2.data}
\alias{get.mutect2.data}
\title{Get MuTect2 Data}
\description{
A function that imports all the MuTect2 annotated text files and returns a full dataframe
containing all the merged data.
}
\usage{
get.mutect2.data(path = ".", pattern = "hg19_multianno.txt$")
}
\arguments{
  \item{path}{Path containing the MuTect2 output text files [required]}
  \item{pattern}{The filename pattern to search [default: hg19_multianno.txt$]}
}
\value{
Returns a dataframe containing all the variant data from MuTect2 and ANNOVAR.
}
\author{
Richard J. de Borja <richard.deborja@sickkids.ca>
}
\keyword{ ssm }
\keyword{ mutect2 }
\keyword{ gatk }
