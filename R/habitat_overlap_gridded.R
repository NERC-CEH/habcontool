
#' Find Overlapping Habitats across large areas and Save Files
#'
#' Identifies overlapping habitats within a spatial dataset, applies various filters and transformations, 
#' and optionally saves the resulting spatial objects as files.
#'
#' @param spatial_object_loc `character` Path to the spatial data file (e.g., a geodatabase).
#' @param SQL_query `character` SQL query to select specific habitat data from the spatial dataset.
#' @param habitat_column_name `character` Name of the column containing habitat information.
#' @param buffer_distance `numeric` Distance to buffer habitats for analysis. Assumed to be in metres, unless otherwise defined as a "units" class.
#' @param connection_distance `numeric` Maximum allowable distance between habitats to consider them connected. Assumed to be in metres, unless otherwise defined as a "units" class.
#' @param min_hab_area `numeric` Minimum habitat area to retain in the analysis.
#' @param combine_touching_polys `logical` Whether to combine touching polygons into single features. Defaults to `TRUE`.
#' @param combine_close_polys `logical` Whether to combine polygons within `connection_distance` into single features. Defaults to `TRUE`.
#' @param plot_it `logical` Whether to generate and display plots of the habitat connectivity process. Defaults to `FALSE`.
#' @param resol `numeric` Resolution of the analysis raster, specified as a numeric vector (e.g., `c(10, 10)`).
#' @param extent_large `numeric` Optional extent for cropping the larger region, provided as `c(xmin, xmax, ymin, ymax)`. Defaults to `NULL`.
#' @param extent_central `numeric` Optional extent for cropping the central region, provided as `c(xmin, xmax, ymin, ymax)`. Defaults to `NULL`.
#' @param save `logical` Whether to save the resulting spatial objects to disk. Defaults to `FALSE`.
#' @param save_loc `character` Directory path where results should be saved. Required if `save = TRUE`.
#' @param save_name `character` File name prefix for the saved results. Required if `save = TRUE`.
#'
#' @return Invisibly returns a `sf` object representing the cropped central region with overlapping habitats.
#'
#' @details 
#' This function processes spatial data to identify overlapping habitats. It applies buffering, 
#' connection distance thresholds, and area filtering, and allows optional saving of the results. 
#' It can also handle polygon-to-raster conversion and crop results to specified extents. 
#'
#' Required helper functions must be sourced prior to using this function:
#' - `habitat_connectivity_lotus_function.R`
#' - `Hab_connect_tool_R.R`
#' - `filter_min_area.R`
#' - `combine_touching_polys.R`
#' - `poly_to_rast.R`
#' - `is_within_dist.R`
#'
#' @note Ensure that the required packages (`tidyverse`, `sf`, `viridis`, `terra`, `patchwork`) are installed and loaded.
#'
#' @examples
#' \dontrun{
#' habitat_connect(
#'   spatial_object_loc = "path/to/spatial_data.gdb",
#'   SQL_query = "SELECT * FROM habitats WHERE type = 'broadleaf_woodland'",
#'   habitat_column_name = "habitat_type",
#'   buffer_distance = 100,
#'   connection_distance = 200,
#'   min_hab_area = 500,
#'   combine_touching_polys = TRUE,
#'   combine_close_polys = TRUE,
#'   plot_it = TRUE,
#'   resol = c(10, 10),
#'   extent_large = c(-1000, 1000, -1000, 1000),
#'   extent_central = c(-500, 500, -500, 500),
#'   save = TRUE,
#'   save_loc = "output/directory",
#'   save_name = "habitat_results"
#' )
#' }
#'
#' @importFrom sf st_read st_write st_crop
#' @importFrom terra crop ext writeRaster
#' @importFrom tidyverse %>%
#' @export

habitat_overlap_gridded <- function(spatial_object_loc, SQL_query, habitat_column_name, buffer_distance,
                                    connection_distance, min_hab_area, combine_touching_polys = TRUE,
                                    combine_close_polys = TRUE, plot_it = FALSE, resol = c(10,10),
                                    extent_large = NULL, extent_central = NULL, save = FALSE, save_loc, save_name) {
  
  
  # prevent scientific notation
  options(scipen=999)
  
  # convert concatenated values to proper structures
  resol = strsplit(resol, '_')[[1]]
  if(!is.null(extent_large)) ext_large = as.numeric(strsplit(extent_large, '_')[[1]])
  if(!is.null(extent_central)) ext_central = as.numeric(strsplit(extent_central, '_')[[1]])
  
  # load in geodatabase - SQL query to get broadleaf woodland
  spatial_object <- sf::st_read(spatial_object_loc, 
                                query = SQL_query)
  
  ## run habitat connectivity function
  overlap_hab <- habitat_overlap(spatial_object = spatial_object,
                                 habitat_column_name = habitat_column_name,
                                 buffer_distance = buffer_distance,
                                 connection_distance = connection_distance,
                                 min_area = min_hab_area,
                                 combine_touching_polys = combine_touching_polys,
                                 combine_close_polys = combine_close_polys,
                                 plot_it = plot_it,
                                 resol = as.numeric(resol),
                                 extent = ext_large)
  
  # get only overlapping habitats
  overs_only <- overlap_hab$habitat_connectivity_raster
  
  # get large patches only - warning can be ignored
  large_only <- filter_min_area(overs_only, min_hab_area)
  
  # crop everything as a polygon
  print("!! Cropping to central region")
  central_only_poly <- st_crop(large_only, 
                               xmin = ext_central[1],
                               ymin = ext_central[3], 
                               xmax = ext_central[2],
                               ymax = ext_central[4])
  
  if(dim(central_only_poly)[1] == 0) stop('No polygons in central square after cropping larger extent')
  
  if(save) {
    
    dir.create(paste0(save_loc, '/central_squares/min_hab_area', min_hab_area), recursive = TRUE)
    
    sf::st_write(central_only_poly, 
                 dsn = paste0(save_loc, '/central_squares/min_hab_area', min_hab_area, '/', 
                              save_name, extent_central, '_buff', buffer_distance,
                              '_conn', connection_distance, '_habarea', min_hab_area, '.shp'),
                 append = FALSE)
    
  }
  
  
}
