## function to treat a RasterLayer as non-geographic
## when we have coordinates we have one for every cell 
## (and propagate those to corner values for each cell in a mesh)
raster0 <- function(x, ...) {
  if (inherits(x, "BasicRaster")) {
    raster::setExtent(x, raster::extent(0, ncol(x), 0, nrow(x)))
  } else {
    x <- raster::raster(x, ...)
    raster::setExtent(x, raster::extent(0, ncol(x), 0, nrow(x)))
  }
}
f <- "PreInd_CISO.climatology.PO4.remap.nc"
library(tidync)
library(anglr)
(x <- tidync(f) )

## names (below we also need to know 'time' and 'lon' for hyper_filter)
varname <- "PO4"
xname <- "lat"
yname <- "z_t"

## steps
## 1.  get the array of data and its coordinates, these are rasters r, x, y in index extent

## this is almost a lvar/band job for raster() but we can't control the orientation
## options: tidync, stars st[,91,,,1, drop=T], angstroms::roms_yz

## a is a list of arrays, sliced as per hyper_filter
a <- x %>% hyper_filter(time = index == 1, 
lon = index == which.min(abs(lon - 180))) %>% 
  hyper_array(select_var = varname)

r <- raster0(t(a[[varname]]))

## the coordinates are on $transforms, all values are present but selected FALSE/TRUE
tr <- lapply(attr(a, "transforms"), dplyr::filter, selected)

x <- raster0(t(matrix(tr[[xname]][[xname]], nrow(tr[[xname]]), ncol = nrow(r))))
y <- raster0(t(matrix(rep(-tr[[yname]][[yname]], each = nrow(tr[[xname]])), nrow(tr[[xname]]), ncol = nrow(r))))


## 2. convert to corner-oriented mesh, we use a trick where each z value $vb[3, ] is propagated to the 
##   correct corner coordinate by as.mesh3d
##  must use triangles = FALSE else the mesh is created differently depending on the input field
mesh <- as.mesh3d(r, triangles = FALSE)
mesh$vb[1, ] <- as.mesh3d(x, triangles = FALSE)$vb[3, ]
mesh$vb[2, ] <- as.mesh3d(y, triangles = FALSE)$vb[3, ]

## alternatively we can create triangles after the fact
#mesh$it <- quadmesh:::triangulate_quads(mesh$ib)
#mesh$ib <- NULL

## 3.  convert to expanded data frame, every coordinate of every cell
fortify.mesh3d <- function(x, ...) {
  idx <- if (!is.null(x$it)) x$it else x$ib
  nc <- dim(idx)[2L]
  nr <- dim(idx)[1L]
  idx <- as.vector(idx)
  xx <- x ## workaround the tibble name-steal
  tibble::tibble(x = xx$vb[1L, idx],
                 y = xx$vb[2L, idx],
                 z = xx$vb[3L, idx],
                 group = rep(seq_len(nc), each = nr))
}
gg <- fortify.mesh3d(mesh) %>% dplyr::filter(!is.na(z))
names(gg) <- c("lat", "z_t", "PO4", "group")  ## we used generic names in the mesh conversion
library(ggplot2)
ggplot() + geom_polygon(data = gg, aes(lat, z_t/100, group = group, fill = PO4)) + 
  theme_bw() + 
  labs(x = "Latitude", y = "Depth (m)")+
  labs(title = "Pre-Industrial Phosphate",
       subtitle = "Transect at 180 Longitude")+
  theme(plot.title=element_text(hjust=0.5)) + 
  theme(plot.subtitle = element_text(hjust=0.5)) + 
  viridis::scale_fill_viridis(option = "D")







