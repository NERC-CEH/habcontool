
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
#' @importFrom dplyr %>%
#' @export
habitat_overlap_gridded <- function(spatial_object, 
                                    SQL_query = NULL, 
                                    habitat_column_name, 
                                    buffer_distance = 500,
                                    min_hab_area = NULL, 
                                    combine_touching_polys = TRUE,
                                    combine_close_polys = TRUE,
                                    connection_distance = 500,  
                                    plot_it = FALSE, 
                                    resolution = c(10,10),
                                    extent_large = NULL, 
                                    extent_central = NULL, 
                                    save = FALSE, 
                                    save_loc, 
                                    save_name) {
  
  
  # prevent scientific notation
  options(scipen=999)
  
  # convert concatenated values to proper structures
  if(any(grepl("_", resolution))) 
    resolution = strsplit(resolution, '_')[[1]]
  if(!is.null(extent_large) & 
     is.character(extent_large)) extent_large = as.numeric(strsplit(extent_large, '_')[[1]])
  if(!is.null(extent_central) & 
     is.character(extent_central)) extent_central = as.numeric(strsplit(extent_central, '_')[[1]])
  
  # load in geodatabase - SQL query to get broadleaf woodland
  if(is.character(spatial_object)){
    spatial_object <- sf::st_read(spatial_object, 
                                  query = SQL_query)
  } else if (!inherits(spatial_object, "sf")) {
    stop("`spatial_object` must be of class `sf`")
  } 
  
  ## from here need to lapply through the different grids if nrow of 
  ## grids >1. Need to change saving process too if want to save each square
  ## separately.
  
  ## run habitat connectivity function
  
  # if the extent is a list, lapply through all grids, if not, run once
  if(inherits(extent_large, "list")) {
    overlap_hab <- lapply(1:length(extent_large), function(x){
      
      tryCatch(
        {
          habitat_overlap(spatial_object = spatial_object,
                          habitat_column_name = habitat_column_name,
                          buffer_distance = buffer_distance,
                          connection_distance = connection_distance,
                          min_area = min_hab_area,
                          combine_touching_polys = combine_touching_polys,
                          combine_close_polys = combine_close_polys,
                          plot_it = plot_it,
                          resolution = as.numeric(resolution),
                          extent = extent_large[[x]])$habitat_connectivity_raster
          
        },
        
        error = function(cond) NULL
      )
      
    })
  } else if(inherits(extent_large, "numeric")) {
    
    tryCatch(
      {
        overlap_hab <- habitat_overlap(spatial_object = spatial_object,
                                       habitat_column_name = habitat_column_name,
                                       buffer_distance = buffer_distance,
                                       connection_distance = connection_distance,
                                       min_area = min_hab_area,
                                       combine_touching_polys = combine_touching_polys,
                                       combine_close_polys = combine_close_polys,
                                       plot_it = plot_it,
                                       resolution = as.numeric(resolution),
                                       extent = extent_large)$habitat_connectivity_raster
        
      },
      
      error = function(cond) NULL
    )
    
  } else {
    stop("`extent_large` must either be a single numeric vector with named elements 
         `xmin`, `ymin`, `xmax` and `ymax`, or a list of named vectors.")
  }
  # # get only overlapping habitats
  # overs_only <- overlap_hab$habitat_connectivity_raster
  
  central_only_poly <- lapply(1:length(overlap_hab), function(x) {
    
    tryCatch(
      {
        # check to see if there are any overlaps
        if(all(is.na(unique(terra::values(overlap_hab[[x]])))))
          stop("!! No overlaps in area")
        
        # get large patches only - warning can be ignored
        if(!is.null(min_hab_area)) {
          large_only <- filter_min_area(overlap_hab[[x]], min_hab_area)
        } else {
          large_only <- rast_to_poly(overlap_hab[[x]])
        }
        
        # crop everything as a polygon
        message("!! Cropping to central region")
        if(inherits(extent_central, "list")) {
        central_only_poly <- sf::st_crop(large_only, extent_central[[x]])
        } else if(inherits(extent_large, "numeric")) {
          central_only_poly <- sf::st_crop(large_only, extent_central)
        } else {
          stop("`extent_large` must either be a single numeric vector with named elements 
         `xmin`, `ymin`, `xmax` and `ymax`, or a list of named vectors.")
        }
        return(central_only_poly)
      },
      error = function(e) NULL
      
    )
  })
  
  # combine the grids into one
  central_only_poly <- do.call(rbind, central_only_poly)
  
  central_only_poly <- central_only_poly %>% 
    dplyr::summarise(geometry = st_union(geometry)) %>% 
    sf::st_cast("POLYGON")
  
  # if(dim(central_only_poly)[1] == 0) stop('No polygons in central square after cropping larger extent')
  
  ###### NEED TO CHECK WHETHER THIS WORKS WITH ONE GRID!!!!!
  if(save) {
    
    dir.create(paste0(save_loc, '/central_squares/min_hab_area', min_hab_area), recursive = TRUE)
    
    sf::st_write(central_only_poly, 
                 dsn = paste0(save_loc, '/central_squares/min_hab_area', min_hab_area, '/', 
                              save_name, paste(extent_central, collapse = "_"), '_buff', buffer_distance,
                              '_conn', connection_distance, '_habarea', min_hab_area, '.shp'),
                 append = FALSE)
    
  }
  
  return(central_only_poly)
  
  
}
