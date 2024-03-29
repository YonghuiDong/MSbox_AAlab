% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/what2.R
\name{what2}
\alias{what2}
\title{metabolite/lipid identification}
\usage{
what2(
  xset,
  RT.DB = NULL,
  use.DB = "METABOLOMICS",
  mode = "+",
  RT.cor = TRUE,
  RT.cor.method = "L",
  frag = T,
  ppm = 10,
  rt = 25,
  ms2.rm = T,
  subgroup = NULL,
  scale_group = NULL,
  scale_factor = 1
)
}
\arguments{
\item{xset}{xcms object}

\item{RT.DB}{correction RT result from RT_correct()}

\item{use.DB}{which database is used for peak identification, choose between "Lipidomics" and "Metabolomics".}

\item{mode}{ionization mode, positive "+" or negative "-".}

\item{RT.cor}{should use corrected RT}

\item{RT.cor.method}{which corrected RT should be used? L: linear regression (default) or P: polynomial regression?}

\item{frag}{should search fragment}

\item{ppm}{mass tolerance, default value = 10}

\item{rt}{retention time search window, default value = 25}

\item{ms2.rm}{should remove MSMS data when it is included, default = TRUE}

\item{subgroup}{subset only specific groups for statistical test, default = NULL, which means statistics will be performed to all group pairs.}

\item{scale_group}{select groups needs to be scaled up or down, accoding to sample dilution or concentration.}

\item{scale_factor}{the scale factor, default = 1.}
}
\description{
metabolite/lipid identification using home database accoring to m/z and RT. MSMS fragments are also given.
}
\examples{
\dontrun{
my_compound <- what2(xset)
}
}
\author{
Yonghui Dong
}
