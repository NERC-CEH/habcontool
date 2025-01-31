
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
#' @param quiet `logical` Whether to print progress messages.
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
#'   save_name = "habitat_results",
#'   quiet = FALSE
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
                                    combine_grid = TRUE,  
                                    plot_it = FALSE, 
                                    resolution = c(10,10),
                                    extent_large = NULL, 
                                    extent_central = NULL, 
                                    save = FALSE, 
                                    save_loc, 
                                    save_name,
                                    quiet = FALSE) {
  
  
  # # prevent scientific notation
  # options(scipen=999)
  
  # convert concatenated values to proper structures
  if(any(grepl("_", resolution))) 
    resolution = strsplit(resolution, '_')[[1]]
  if(!is.null(extent_large) &
     is.character(extent_large)) extent_large = as.numeric(strsplit(extent_large, '_')[[1]])
  # if(!is.null(extent_central) & 
  #    is.character(extent_central)) extent_central = as.numeric(strsplit(extent_central, '_')[[1]])
  
  # set units for min_hab_area
  if(!is.null(min_hab_area)){
    if(class(min_hab_area) != "units") {
      if(!quiet)
        message("assuming 'min_area' is provided in metres^2")
      min_hab_area <- units::set_units(min_hab_area, 'm^2')
    }
  }
  
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
  if(inherits(extent_large, "list") || 
     inherits(extent_large, "sfc") ||
     inherits(extent_large, "sfg")) {
    
    overlap_hab <- lapply(cli::cli_progress_along(extent_large), function(x){
      
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
                          extent = extent_large[[x]],
                          quiet = quiet)$habitat_connectivity_raster
          
        }, 
        error = function(cond) warning("!! `habitat_overlap` failed with: ", cond)
      )
      
    })
  } else if(inherits(extent_large, "numeric")) {
    
    overlap_hab <- tryCatch(
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
                        extent = extent_large,
                        quiet = quiet)$habitat_connectivity_raster
        
      },
      
      error = function(cond) stop(cond) 
    )
    
  } else {
    stop("`extent_large` must be a single numeric vector with named elements 
         `xmin`, `ymin`, `xmax` and `ymax`, a list of named vectors, or an object 
         which van be passed to `sf::st_crop`.")
  }
  
  # combine to get all overlaps together
  # remove the list entries that have error messages
  if(inherits(overlap_hab, "list"))
    overlap_hab <- overlap_hab[!sapply(overlap_hab, is.character)]
  
  # combine overlaps together
  overlaps_sprc <- terra::sprc(overlap_hab)
  overlaps_mos <- terra::mosaic(overlaps_sprc, fun = "max")
  
  # remove original polygons
  orig_polys <- sum(poly_to_rast(spatial_object, field_val = 1, 
                                 resolution = resolution, 
                                 rast_extent = terra::ext(overlaps_mos), 
                                 layer_names = NULL), 
                    na.rm = TRUE)
  
  overlaps_mos <- terra::mask(overlaps_mos, orig_polys, inverse = TRUE)
  names(overlaps_mos) <- "n_overlaps"
  
  # check to see if there are any overlaps
  if(!all(terra::global(overlaps_mos, fun = "anynotNA")$anynotNA))
    stop("!! No overlaps in area")
  
  # filter the minimum area and end up converting to polygons 
  # convert it to raster!!!!!
  if(!is.null(min_hab_area)) {
    overlaps_mos <- filter_min_area(spatial_object = overlaps_mos, 
                                    min_area = min_hab_area, 
                                    combine_touching_polys = TRUE,
                                    quiet = TRUE) 
  } else {
    overlaps_mos <- rast_to_poly(overlaps_mos)
  }
  
  ## I do want to filter polygons after overlapping because too small habitat is bad
  ## need to do it AFTER combining all grids though, because it's possible area 
  # increases during the combining process.
  
  #### SHOULD I CONVERT EVERYTHING TO A RASTER OR OUTPUT AS A POLYGON?!?!?!?!
  
  if(save) {
    
    ## Save as a raster, not polygon?
    overlaps_mos
    
    message("!! Saving combined gridded polygons")
    
    dir.create(paste0(save_loc, '/combined_grids/min_hab_area', min_hab_area), recursive = TRUE)
    
    sf::st_write(combined_poly, 
                 dsn = paste0(save_loc, '/combined_grids/min_hab_area', min_hab_area, '/', 
                              save_name, '_buff', buffer_distance,
                              '_conn', connection_distance, '_habarea', min_hab_area, '.shp'),
                 append = FALSE)
    
  }
  
  # } else {
  
  if(save) {
    
    message("!! Saving individual polygon grids")
    
    dir.create(paste0(save_loc, '/central_squares/min_hab_area', min_hab_area), recursive = TRUE)
    
    sf::st_write(central_only_poly, 
                 dsn = paste0(save_loc, '/central_squares/min_hab_area', min_hab_area, '/', 
                              save_name, paste(extent_central, collapse = "_"), '_buff', buffer_distance,
                              '_conn', connection_distance, '_habarea', min_hab_area, '.shp'),
                 append = FALSE)
    
  }
  
  
  # }
  
  return(overlaps_mos)
  
}
