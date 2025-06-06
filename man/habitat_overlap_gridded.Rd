% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/habitat_overlap_gridded.R
\name{habitat_overlap_gridded}
\alias{habitat_overlap_gridded}
\title{Find Overlapping Habitats across large areas and Save Files}
\usage{
habitat_overlap_gridded(
  spatial_object,
  SQL_query = NULL,
  wkt_filter = character(0),
  habitat_column_name = NULL,
  buffer_distance,
  min_core_area = NULL,
  min_area = NULL,
  combine_touching_polys = TRUE,
  combine_close_polys = FALSE,
  connection_distance,
  return_rast = TRUE,
  plot_it = FALSE,
  resolution = c(10, 10),
  extent = NULL,
  save = FALSE,
  save_loc,
  save_name,
  quiet = FALSE
)
}
\arguments{
\item{spatial_object}{\code{character} An \code{sf} object or path to a spatial polygon file that can be read by \code{sf::st_read} (e.g., a geodatabase).}

\item{SQL_query}{\code{character} SQL query to select specific habitat data from the spatial dataset.}

\item{wkt_filter}{\code{character} WKT representation of a spatial filter}

\item{habitat_column_name}{\code{character} Name of the column containing habitat information.}

\item{buffer_distance}{\code{numeric} Distance to buffer habitats for analysis. Assumed to be in metres, unless otherwise defined as a "units" class.}

\item{min_core_area}{\code{numeric} Minimum area of "core" (i.e. contiguous) habitat. If not \code{NULL}, filters out polygons smaller than this value before combining polygons within \code{connect_dist} of each other. Defaults to \code{NULL}.}

\item{min_area}{\code{numeric} Minimum polygon area. Polygons below this threshold are removed prior to buffering and after connected areas have been identified}

\item{combine_touching_polys}{\code{logical} Whether to combine touching polygons into single features. Defaults to \code{TRUE}.}

\item{combine_close_polys}{\code{logical} Whether to combine polygons within \code{connection_distance} into single features. Defaults to \code{FALSE}.}

\item{connection_distance}{\code{numeric} Maximum allowable distance between habitats to consider them connected. Assumed to be in metres, unless otherwise defined as a "units" class.}

\item{plot_it}{\code{logical} Whether to generate and display plots of the habitat connectivity process. Defaults to \code{FALSE}.}

\item{extent}{\code{numeric} Optional extent for cropping the larger region, provided as (\code{xmin}, \code{ymin}, \code{xmax}, \code{ymax}). Defaults to \code{NULL}.}

\item{save}{\code{logical} Whether to save the resulting spatial objects to disk. Defaults to \code{FALSE}.}

\item{save_loc}{\code{character} Directory path where results should be saved. Required if \code{save = TRUE}.}

\item{save_name}{\code{character} File name prefix for the saved results. Required if \code{save = TRUE}.}

\item{quiet}{\code{logical} Whether to print progress messages.}

\item{resol}{\code{numeric} Resolution of the analysis raster, specified as a numeric vector (e.g., \code{c(10, 10)}).}
}
\value{
Invisibly returns a \code{sf} object representing the cropped central region with overlapping habitats.
}
\description{
Identifies overlapping habitats within a spatial dataset, applies various filters and transformations,
and optionally saves the resulting spatial objects as files.
}
\details{
This function processes spatial data to identify overlapping habitats. It applies buffering,
connection distance thresholds, and area filtering, and allows optional saving of the results.
It can also handle polygon-to-raster conversion and crop results to specified extents.
}
\note{
Ensure that the required packages (\code{tidyverse}, \code{sf}, \code{viridis}, \code{terra}, \code{patchwork}) are installed and loaded.
}
\examples{
\dontrun{
habitat_connect(
  spatial_object_loc = "path/to/spatial_data.gdb",
  SQL_query = "SELECT * FROM habitats WHERE type = 'broadleaf_woodland'",
  habitat_column_name = "habitat_type",
  buffer_distance = 100,
  connection_distance = 200,
  min_area = 500,
  combine_touching_polys = TRUE,
  combine_close_polys = TRUE,
  plot_it = TRUE,
  resol = c(10, 10),
  extent = c(-1000, 1000, -1000, 1000),
  save = TRUE,
  save_loc = "output/directory",
  save_name = "habitat_results",
  quiet = FALSE
)
}

}
