# habcontool: Identifying areas to improve habitat connectivity

This is an R package written for the ERAMMP and EVAST projects which is used to is used to determine spatial locations in which planting new trees/introducing new habitats would connect two or more habitat patches together. It relies heavily on the excellent _sf_ (https://r-spatial.github.io/sf/) and _terra_ (https://rspatial.github.io/terra/) packages. To install:

```{r install}

remotes::install_github("NERC-CEH/habcontool")
library(habcontool)

```

## Brief overview

In brief, this package buffers polygons from the _sf_ package and identifies areas where polygons overlap. This is to determine spatial locations in which introducing more of a given habitat would increase connectivity between existing habitat patches. While this code has been written in the context of tree planting under agri-environmnet schemes, it can be used for any habitat (or anything, so long as it is given in polygon format). For example, some schemes might require tree planting. This package identifies locations in which planting of trees would connect two or more habitat patches. The outputs could therefore be used to inform where planting extra trees would be most beneficial for connectivity. Furthermore, the functions also accounts for the dispersal abilities of the species/habitats of interest. For example, some birds are highly mobile, meaning that patches of suitable habitat may already be considered "connected" if they are within a certain distance of each other. Our function works by first identifying habitat patches that are already connected (as definied by a user-supplied distance) and then calculating the effect of planting/habitat management on connectivity between these patches. 

## Issues with measuring habitat connectivity

One of the issues with measuring habitat connectivity is scalability. When processing large spatial dataframes, you may end up with an extremely large computational load. Furthermore, it is difficult to split the job up and work on smaller subsets of the main data. If you cut the area down to a smaller square, you need to consider whether patches within your area of interest overlap with patches outside of your area. We address this in the `habitat_connect` function (see below), in which you can supply two extents - an inner and an outer area. The difference between these two is the maximum distance that you expect connectivity to happen at. This should determined for each use case based on your project questions/planting area or species of interest.

This function **should** work at any scale, from local patches at fine resolution to (inter)national scale habitats. Some words of warning though. These have been designed for specific use cases and as such, we cannot ensure transferability to other uses.

## Disclaimer: As slow as a sloth 

The functions are almost certainly not well optimised and are, unfortunately, slow (due to time pressures and a rubbish package author yaddah yaddah). When we ran the code across a large spatial scale (entirity of England), we used a computer cluster to speed up processing. We would highly recommend **not** running a large spatial area on your own computer (although if you are running a coarse enough resolution it might work fine).

## Package structure

### Function `habitat_overlap`

The function in this package that does the heavy lifting is `habitat_overlap`. This identifies the parcels that overlap and returns a named list comprising five different spatial aggregations generated by the function. It does this using a series of spatial helper functions. While this function can be used by itself, to identify only overlapping habitats over a large area we access it through the `habitat_connect` wrapper function.

### Function `habitat_connect`

This is a wrapper function for the `habitat_overlap` function. This should only really be used when obtaining overlaps across a large spatial scale. To do this, you provide the function with two extents, one for an outer square and one for an inner square. You should choose the relative sizes of each square based on your project questions/planting area or species of interest.  

It runs the `habitat_overlap` and crops it to a central square for ease of joining different squares together.

### Helper functions

There are also several helper functions which are used by the packages used above. We won't explain these properly here, but they do/help with various spatial manipulations which are required for the workflow above to work.

## Simple example

For a simple example using data freely available online (XXXXX), please see the vignette (yyyyy). 


