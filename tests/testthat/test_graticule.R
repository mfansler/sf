context("sf: graticule")

test_that("graticule 1", {
  library(maps)
  world1 <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))
  world2 <- st_transform(world1, "+proj=laea +y_0=0 +lon_0=155 +lat_0=-90 +ellps=WGS84 +no_defs")
  plot(world2, graticule = st_crs(4326), axes = TRUE)
})
