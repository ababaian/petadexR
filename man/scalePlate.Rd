% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scalePlate.R
\name{scalePlate}
\alias{scalePlate}
\title{Scale a measurement column between plate-level controls}
\usage{
scalePlate(
  plate.tbl,
  value.col = "colony.value",
  negative.control,
  positive.control,
  suffix = "scaled"
)
}
\arguments{
\item{plate.tbl}{A tibble returned by \code{\link[=readPlate]{readPlate()}}; must contain
\code{plateID}, \code{media}, \code{timepoint}, \code{gene}, and the column named in
\code{value.col}.}

\item{value.col}{\code{character(1)} – name of the \strong{numeric} column to
normalise.  \strong{Default:} \code{"colony.value"}.}

\item{negative.control}{\code{character(1)} – value in \code{gene} that labels
the negative-control wells (baseline \strong{0}).}

\item{positive.control}{\code{character(1)} – value in \code{gene} that labels
the positive-control wells (baseline \strong{1}).}

\item{suffix}{\code{character(1)} – suffix for the new column name.}
}
\value{
The input tibble plus one numeric column
\verb{<value.col>_<suffix>} containing the scaled values.
}
\description{
\code{scalePlate()} rescales a numeric measurement column (default
\strong{\code{colony.value}}) so that the \strong{mean} of \emph{negative-control} wells
maps to \strong{0} and the \strong{mean} of \emph{positive-control} wells maps to
\strong{1}.  All other wells are placed on the same linear scale.
}
\details{
Controls are identified via string equality in the \strong{\code{gene}} column.
Scaling is performed \strong{within each} \verb{\{plateID, media, timepoint\}}
triple so that different plates, media, and sampling times remain
independent.
}
