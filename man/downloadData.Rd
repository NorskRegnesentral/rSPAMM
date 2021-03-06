% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inputFunctions.R
\name{downloadData}
\alias{downloadData}
\title{Download Repository Data}
\usage{
downloadData(
  url = "https://github.com/NorskRegnesentral/HarpAndHoodSealData/archive/master.zip",
  chooseFolder = TRUE
)
}
\arguments{
\item{url}{The URL of the data repository.}

\item{chooseFolder}{Logical parameter if you want to specify to which folder to download the data.}
}
\value{
A folder containing the Data for the various populations and a Scripts folder for an example of complete analysis of the data.
}
\description{
Downloading the most updated data on harp seals in the West Ice (along the coast of Greenland) and in the East Ice (the White Sea) and the Hooded seal population in the West Ice from the Seal Data repository.
}
\examples{
downloadData(url = "")
}
\keyword{data}
\keyword{input,}
