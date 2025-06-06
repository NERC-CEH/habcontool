% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_manip_functions.R
\name{combine_overlap_gridded}
\alias{combine_overlap_gridded}
\title{Combine Overlapping Gridded Rasters}
\usage{
combine_overlap_gridded(
  rast_folder,
  fun = "max",
  min_hab_area = NULL,
  write = TRUE,
  return_rast = TRUE,
  save_loc = NULL
)
}
\arguments{
\item{rast_folder}{Character. Path to the folder containing raster files.}

\item{fun}{Character. The function used to combine overlapping rasters (e.g., \code{"max"}, \code{"mean"}, \code{"sum"}). Default is \code{"max"}.}

\item{write}{Logical. Whether to save the resulting raster. Default is \code{TRUE}.}

\item{save_loc}{Character. Path to the folder where the output raster will be saved. Required if \code{write = TRUE}.}
}
\value{
A \code{SpatRaster} object representing the mosaicked raster.
}
\description{
DISCLAIMER: This documentation and the associated function were generated by ChatGPT.
Please verify and adapt the function to suit your specific requirements and data.
}
\details{
This function reads all raster files from a specified folder, mosaics them using a specified function, and optionally saves the output.
\itemize{
\item The function lists all raster files in \code{rast_folder} and creates a spatial collection.
\item It then mosaics the rasters using the specified function (\code{fun}).
\item If \code{write = TRUE}, the function creates the \code{save_loc} directory (if not existing) and saves the result as a \code{.tif} file.
}
}
\examples{
\dontrun{
combine_overlap_gridded("path/to/rasters", fun = "mean", write = TRUE, save_loc = "path/to/save")
}

}
