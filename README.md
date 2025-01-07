# Habitat Connectivity tool

This is an R package written for the ERAMMP and EVAST projects which is used to measure the impact of planting on the connectivity of habitats. It relies heavily on the excellent _sf_ (https://r-spatial.github.io/sf/) and _terra_ (https://rspatial.github.io/terra/) packages. 

In brief, this works by buffering polygons and identifying areas where polygons overlap. However, it also accounts for the dispersal abilities of the species/habitats of interest. For example, some birds are highly mobile, meaning that patches of suitable habitat may already be considered "connected" if they are within a certain distance of each other. Our function works by first identifying habitat patches that are already connected (as definied by a user-supplied distance) and then calculating the effect of planting/habitat managedment on these connected parcels.     

This function **should** work at any scale, from local pacthes at fine resolution to (inter)national scale habitats. Some words of warning though. These have been designed for specific use cases and as such, we cannot ensure transferability to other uses. Furthermore, the functions are almost certainly not well optimised and are, unfortunately, slow (due to time pressures and a rubbish package author yaddah yaddah). When we ran the code across a large spatial scale (entirity of Wales), we used a cluster computer to speed up processing. We would highly recommend **not** running a large spatial area on your own computer (although if you are running a coarse enough resolution it might work fine). 

