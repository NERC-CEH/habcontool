
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
#' @param combine_close_polys `logical` Whether to combine polygons within `connection_distance` into single features. Defaults to `FALSE`.
#' @param plot_it `logical` Whether to generate and display plots of the habitat connectivity process. Defaults to `FALSE`.
#' @param resol `numeric` Resolution of the analysis raster, specified as a numeric vector (e.g., `c(10, 10)`).
#' @param extent `numeric` Optional extent for cropping the larger region, provided as (`xmin`, `ymin`, `xmax`, `ymax`). Defaults to `NULL`.
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
#'   extent = c(-1000, 1000, -1000, 1000),
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
                                    wkt_filter = character(0),
                                    habitat_column_name = NULL, 
                                    buffer_distance,
                                    min_hab_area = NULL, 
                                    combine_touching_polys = TRUE,
                                    combine_close_polys = FALSE,
                                    connection_distance,
                                    return_rast = TRUE,
                                    plot_it = FALSE, 
                                    resolution = c(10,10),
                                    extent = NULL, 
                                    save = FALSE, 
                                    save_loc, 
                                    save_name,
                                    quiet = FALSE) {
  
  # convert concatenated values to proper structures
  if(any(grepl("_", resolution))) 
    resolution <- as.numeric(strsplit(resolution, '_')[[1]])
  
  if(!is.null(extent) &
     is.character(extent)) 
    extent <- as.numeric(strsplit(extent, '_')[[1]])
  
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
                                  query = SQL_query,
                                  wkt_filter = wkt_filter)
  } else if (!inherits(spatial_object, "sf")) {
    stop("`spatial_object` must be of class `sf`")
  } 
  
  ## from here need to lapply through the different grids if nrow of 
  ## grids >1. Need to change saving process too if want to save each square
  ## separately.
  
  ## run habitat connectivity function
  
  # if the extent is a list, lapply through all grids, if not, run once
  if(inherits(extent, "list") || 
     inherits(extent, "sfc") ||
     inherits(extent, "sfg")) {
    
    overlap_hab <- lapply(cli::cli_progress_along(extent), function(x){
      
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
                          extent = extent[[x]],
                          quiet = quiet)$habitat_connectivity_raster
          
        }, 
        error = function(cond) warning("!! `habitat_overlap` failed with: ", cond)
      )
      
    })
  } else if(inherits(extent, "numeric") | is.null(extent)) {
    
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
                        extent = extent,
                        quiet = quiet)$habitat_connectivity_raster
        
      },
      
      error = function(cond) stop(cond) 
    )
    
  } else {
    stop("`extent` must be a single numeric vector with named elements 
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
  names(overlaps_mos) <- "n_overlaps"
  
  # remove original polygons
  orig_polys <- sum(poly_to_rast(obj = spatial_object, 
                                 field_val = 1, 
                                 resolution = as.numeric(resolution), 
                                 rast_extent = terra::ext(overlaps_mos), 
                                 layer_names = NULL), 
                    na.rm = TRUE)
  
  overlaps_mos <- terra::mask(overlaps_mos, orig_polys, inverse = TRUE)
  
  # check to see if there are any overlaps
  if(!all(terra::global(overlaps_mos, fun = "anynotNA")$anynotNA)){
    message("!! No overlaps in area, returning NULL object")
    
    return(NULL)
  }
  
  # filter the minimum area and end up converting to polygons 
  if(!is.null(min_hab_area)) {
    
    overlaps_mos <- filter_min_area(spatial_object = overlaps_mos, 
                                    min_area = min_hab_area, 
                                    combine_touching_polys = TRUE,
                                    quiet = TRUE,
                                    return_rast = return_rast,
                                    combine_output_rast = TRUE)
    
  } 
  
  if(!return_rast) {
    overlaps_mos <- rast_to_poly(overlaps_mos)
  }
  
  
  ## I do want to filter polygons after overlapping because too small habitat is bad
  ## need to do it AFTER combining all grids though, because it's possible area 
  # increases during the combining process.
  
  if(save) {
    
    if(inherits(extent, "list")) 
      message("!! Saving combined outputs")
    
    if(inherits(extent, "numeric")) 
      message("!! Saving single grid output")
    
    if(is.null(extent) & !is.null(wkt_filter)) {
      message("!! Saving name based on 'wkt_filter'")
      extent <- as.numeric(sf::st_bbox(sf::st_as_sf(x = data.frame(wkt = wkt_filter), wkt = "wkt")))
    }
    
    # create save location
    save_loc_full <- paste0(save_loc, '/',
                            save_name, 
                            '_buff', buffer_distance,
                            '_conn', connection_distance, 
                            '_habarea', min_hab_area)
    
    # create the directory
    dir.create(save_loc_full, 
               recursive = TRUE)
    
    
    if(return_rast) {
      
      terra::writeRaster(overlaps_mos,
                         filename = paste0(save_loc_full, '/',
                                           save_name, 
                                           ifelse(inherits(extent, "numeric"),
                                                  paste0('_extent_', paste(extent, collapse = '_')),''),
                                           '_buff', buffer_distance,
                                           '_conn', connection_distance, '_habarea', 
                                           min_hab_area, '.tif'),
                         overwrite = TRUE)
      
    } else {
      
      sf::st_write(overlaps_mos, 
                   dsn = paste0(save_loc_full, '/', 
                                save_name, 
                                ifelse(inherits(extent, "numeric"),
                                       paste0('_extent_', paste(extent, collapse = '_')),''),
                                '_buff', buffer_distance,
                                '_conn', connection_distance, '_habarea', 
                                min_hab_area, '.shp'),
                   append = FALSE)
    }
    
  }
  
  return(overlaps_mos)
  
}
