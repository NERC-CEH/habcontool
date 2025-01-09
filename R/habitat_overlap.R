## Habitat connectivity tool for the Defra EVAST project
## tool to identify areas where new habitat could be planted to increase habitat connectivity.

# # testing
# library(terra)
# library(sf)
# library(tidyverse)
# library(viridis)
# 
# forest_data <- st_read('data/forestry/National_Forest_Inventory_Woodland_England_2019/National_Forest_Inventory_Woodland_England_2019.shp')
# forest_data
# 
# # crop to small area
# woods <- st_crop(forest_data, xmin = 392500, ymin= 380000, xmax= 395000, ymax=382000)
# plot(st_geometry(woods))
# 
# 
## geodatabase from LCM 2015
# gdb <- sf::st_read("data/landcover_map/LCM15NFI16.gdb/LCM15NFI16.gdb")
# # gdb
# 
# unique(gdb$Land_Cover)
# ## function
# ## used for debugging
# buffer_distance = 500 ## in metres
# min_area = 1000 ## metres squared
# connection_distance = 500 # number of metres to use as a measure of what's 'connected' before buffering
# spatial_object = gdb[gdb$Land_Cover == "Arable and horticulture",] # [grepl('woodland', gdb$Land_Cover),] ## object to check
# extent = data.frame(xmin = 200000, ymin= 200000, xmax= 210000, ymax=210000)
# habitat_column_name = 'FIRST_bhab'
# plot_it = TRUE
# combine_touching_polys = T
# combine_close_polys = TRUE
# resol = c(10, 10)

habitat_overlap <- function(spatial_object, habitat_column_name, extent = NULL, buffer_distance = 500,
                            connection_distance, min_area, combine_touching_polys = TRUE, combine_close_polys = TRUE, 
                            plot_it = TRUE, resol = c(10,10)) {
  
  # crop to region of interest
  if(!is.null(extent)) {
    print('!! cropping object'); object <- st_crop(spatial_object, xmin = extent[[1]], xmax = extent[[2]], 
                                                   ymin = extent[[3]], ymax = extent[[4]], crs = crs(spatial_object))
  } else if(is.null(extent)) object <- spatial_object
  
  if(dim(object)[1] == 0) stop('!! No polygons present after cropping.\nIncrease extent size or change area.')
  
  if(plot_it) {
    p1 <- ggplot() +
      geom_sf(data = object) +
      theme_bw() +
      ggtitle('Initial object')
    # print(p1)
  }
  
  # set units to metres for use in the buffering functions
  buffer_dist <- units::set_units(buffer_distance, 'm')
  min_area <- units::set_units(min_area, 'm^2')
  connection_dist <- units::set_units(connection_distance, 'm')
  
  # Combine touching polygons and those within connection_dist if combine_close == TRUE
  if(combine_touching_polys) {
    
    print('!! combining touching and/or close polygons')
    
    # run function
    comb_object <- combine_touching(comb_obj = object, variable_name = habitat_column_name, 
                                    Plot = plot_it, connect_dist = connection_dist, 
                                    combine_close = combine_close_polys)
    
  } else if (!combine_touching_polys) {
    
    comb_object <- object
    
  }
  
  # filter minimum area 
  if(!is.null(min_area)){
    obj_lrge <- comb_object %>% 
      mutate(area = st_area(st_geometry(.))) %>% 
      filter(area > min_area)
  } else { obj_lrge <- comb_object }
  
  
  if(plot_it) {
    p2 <- ggplot() +
      geom_sf(data = obj_lrge, aes(fill = poly_id)) +
      theme_bw() +
      ggtitle('Combined objects') +
      scale_fill_viridis_d() +
      theme(legend.position = 'none')
    # print(p2)
  }
  
  
  print('!! Buffering polygons')
  
  ## Buffer, calculate number of overlaps (n_overlaps), and number within double
  ## of the buffer specified.
  obj_lrge_buff <- st_buffer(obj_lrge, dist = buffer_dist) %>% 
    mutate(n_overlaps = lengths(st_intersects(.)))  
  
  
  # convert buffered polygon to raster - produces a raster stack with 1s in every non-0 cell
  # for each layer in the raster. Raster should be projected. Can also be given the name
  # of the column in the sf polygon ('field' in the rasterize() function) - to retain values
  # associated with the polygon
  buffered_object_rast <- poly_to_rast(obj = obj_lrge_buff, field_val = 1, 
                                       resolution = resol, 
                                       rast_extent = terra::ext(obj_lrge_buff)+10, 
                                       layer_names = obj_lrge_buff$variable)
  
  # check that the polygon names are in the same order
  if(!identical(obj_lrge_buff$poly_id, obj_lrge$poly_id)) stop("polys aren't in same order")
  
  # take the sum across all layers to get the overlaps (where value > 1)
  buff_obj_sum <- sum(buffered_object_rast, na.rm = TRUE)
  
  
  if(plot_it){
    
    p3 <- ggplot() +
      geom_raster(data = as.data.frame(buff_obj_sum, xy=TRUE), aes(x=x, y=y, fill = sum), alpha = 0.5) +
      geom_sf(data = obj_lrge) +
      theme_bw() +
      scale_fill_viridis(na.value = NA, name = 'Overlaps') +
      ggtitle('Initial object, buffered overlaps')
    # print(p3)
    
  }
  
  ## remove the original polygons from the raster to identify areas for habitat connectivity
  # convert original polygon into raster layer - set field to -1 to identify overlapping regions
  
  # convert original polygon into raster layer.
  # set field to -1 to distinguish them from the empty areas in the raster
  # (i.e. areas with no polygons in them)
  obj_lrge_rast <- sum(poly_to_rast(obj_lrge, field_val = -1, 
                                    resolution = resol, 
                                    rast_extent = terra::ext(obj_lrge_buff)+10, 
                                    layer_names = obj_lrge$variable), 
                       na.rm = TRUE)
  
  # set any areas that == 0 to 1
  # When multiplying this by the overlapping areas, this will turn any areas overlapping
  # with the original objects negative, but leave everything else the same
  values(obj_lrge_rast)[is.na(values(obj_lrge_rast))] <- 1
  
  # multiply new raster with buffered object with overlaps - negative values are
  # where they overlap
  overlaps_only <- buff_obj_sum*obj_lrge_rast
  
  # Set anything <=1 to NA (areas where they overlap original polygon and 
  # where there's only one buffered region)
  values(overlaps_only)[values(overlaps_only)<=1] <- NA
  
  if(plot_it){
    
    p4 <- ggplot() +
      geom_tile(data = as.data.frame(overlaps_only, xy=TRUE), aes(x=x, y=y, fill = factor(sum)), alpha = 0.5) +
      coord_quickmap() +
      theme_bw() +
      scale_fill_viridis_d(na.value = NA, name = 'Overlaps') +
      ggtitle('Overlapping regions only')
    print(p1+p2+p3+p4)
    
  }
  
  return(list(combined_object = comb_object,
              buffered_object = obj_lrge_buff,
              buffered_raster_stack = buffered_object_rast,
              buffered_raster_overlaps = buff_obj_sum,
              habitat_connectivity_raster = overlaps_only))
  
}


# ## geodatabase from LCM 2015
# gdb <- sf::st_read("data/landcover_map/LCM15NFI16.gdb/LCM15NFI16.gdb")
# # gdb
# 
# unique(gdb$Land_Cover)
# 
# bl_overlap <- habitat_overlap(spatial_object = gdb[gdb$Land_Cover == "Broadleaf woodland",],
#                               habitat_column_name = 'FIRST_bhab',
#                               buffer_distance = 200,
#                               connection_distance = 200,
#                               min_area = 1000,
#                               combine_touching_polys = TRUE,
#                               combine_close_polys = TRUE,
#                               plot_it = TRUE,
#                               resol = c(10,10),
#                               extent = data.frame(xmin = 210000, ymin= 210000, xmax= 260000, ymax=260000))
# 
# plot(bl_overlap$habitat_connectivity_raster)
