---
título: "Capítulo 8. Caso práctico. Patrones de localización en el poblamiento de la Edad del hierro de Galicia (Noroeste de la Península Ibérica), perteneciente al libro Carrero-Pazos, M. (2022). *Arqueología Computacional del Territorio. Métodos y técnicas para estudiar decisiones humanas en paisajes pretéritos*. Oxford. Archaeopress."
autor: "Miguel Carrero-Pazos"
date: '2023-02-10'
---

## Cálculo de la densidad de rutas óptimas (Least Cost Paths)
  
# Cargar e instalar paquetes de trabajo específicos
# La primera vez es necesario instalar el paquete fasterRaster
# install.packages("remotes")
# remotes::install_github("adamlilith/fasterRaster")

# Establecer el directorio de datos
sourcedir <- file.path("data","raw_data")
targetdir <- file.path("data","derived_data")
if(!dir.exists(targetdir)) dir.create(targetdir)

# Instalar y cargar paquetes de trabajo específicos
required_packages <- c("rgdal", "rgeos", "sp",
                       "raster","gdistance", 
                       "leastcostpath", "sf",
                       "maptools", "spatstat",
                       "ggplot2", "RColorBrewer",
                       "fasterRaster")

if(!all(required_packages %in% installed.packages())) {
  packages_to_install <- which(!required_packages %in% installed.packages())
  install.packages(required_packages[packages_to_install],
                   dependencies = TRUE,
                   repos = "https://ftp.gwdg.de/pub/misc/cran/")
}

# Nota: Instalar el paquete "terra" desde el CRAN da error, así que mejor desde Github
# install.packages('terra', repos='https://rspatial.r-universe.dev')
library(terra)

# Cargar los paquetes requeridos
lapply(required_packages, require, character.only = TRUE)

# Pasos iniciales
# Cargar el DEM
dem <- raster(file.path(targetdir,"dem.tif"))

# Establecer la resolución del DEM en 150m
dem150m <- aggregate(dem, fact=10)

# Establecer los límites del área de estudio (polígono y línea)
r <- terra::rast(dem150m)
pol <- terra::as.polygons(r > -Inf) %>% st_as_sf
pol_sp <- sf:::as_Spatial(pol$geometry)
pol_line <- as(pol_sp, "SpatialLines")

# Cortar el DEM con el polígono del área de estudio y establecer el CRS
dem_crop <- raster::crop(dem150m, pol_sp)
raster::crs(dem_crop) <- raster::crs(dem150m)

# PARTE 1. Aproximación FETE ####
# Superficie de coste basada en la pendiente y altitud
cfs <- c("tobler", "tobler offpath", "irmischer-clarke male", "irmischer-clarke offpath male", "irmischer-clarke female", "irmischer-clarke offpath female",
         "modified tobler", "wheeled transport", "herzog", "llobera-sluckin", "campbell 2019")

# Establecer el número de vecinos para el cálculo de coste
neigh <- 4 # 4, 8, 16, 32, 48

slope_cs <- leastcostpath::create_slope_cs(dem = dem_crop, 
                                           cost_function = "tobler", 
                                           neighbours = neigh)

# Cortar la superficie de coste con el polígono del área de estudio
slope_cs <- crop_cs(slope_cs, pol_sp)

# Crear puntos regulares a lo largo de los límites del área de estudio
# Note que el código de abajo simplemente crea 100 puntos
# Este número puede ampliarse, aunque ampliará considerablemente el tiempo de procesado

Pts_ext <- spsample(pol_line, n = 100,
                    type='regular')

# Establecer las coordenadas de los puntos de muestreo
crs <- "+proj=utm +zone=29 +ellps=GRS80 +units=m +no_defs"
coords <- c("x", "y")
Pts_sf <- st_as_sf(Pts_ext, coords = coords, crs = crs)

Pts_borders <- Pts_sf %>%
  st_cast("MULTIPOINT") %>%
  st_cast("POINT")

Pts_borders <- st_coordinates(Pts_borders)
# Convertir los puntos regulares de los límites del área de estudio en puntos espaciales
coords <- data.frame(Pts_borders[,1],Pts_borders[,2])
datcoords <- cbind(1,coords)
locs <- sp::SpatialPointsDataFrame(coords,datcoords)
proj4string(locs) <- CRS("+proj=utm +zone=29 
                         +ellps=GRS80 +units=m +no_defs")

# Mapear los puntos de origen y destino sobre la superficie de costes
plot(raster(slope_cs))
plot(locs, add=T)

## FETE (From everywhere to everywhere) ----
# El código de abajo calculará las rutas de menor coste para todas las localizaciones proporcionadas,
# tomando un punto como origen y el resto como destino. El proceso se
# se repite para un segundo punto, y así sucesivamente hasta el total
fete_lcp <- leastcostpath::create_FETE_lcps(cost_surface = slope_cs, 
                                            locations = locs, 
                                            cost_distance = FALSE, 
                                            parallel = TRUE)

# Guardar el vector de las rutas óptimas como shapefile
 writeOGR(fete_lcp, 
         file.path(targetdir, "fete_sampleo_lcps.shp"),
         delete_dsn = TRUE,
         layer="fete_sampleo_lcps",
         driver = "ESRI Shapefile",
         overwrite = TRUE)

# Para cargar el shapefile de rutas, ejecutar el siguiente código:
fete_lcp <- sf::st_read(file.path(targetdir,"fete_sampleo_lcps.shp"),
                     layer = "fete_sampleo_lcps")
fete_lcp <- as(fete_lcp, "Spatial")

#Convertir los caminos de menor coste en puntos (input para el análisis de densidad)
lcp_pts <- spsample(fete_lcp, n = 1000000,
                    type='regular')

# Establecer las coordenadas de los puntos a lo largo de las rutas de menor coste
crs <- "+proj=utm +zone=29 +ellps=GRS80 +units=m +no_defs"
coords_lcp_pts <- c("x", "y")
lcp_pts_sf <- st_as_sf(lcp_pts, coords = coords_lcp_pts, 
                       crs = crs)

lcp_pts_sp <- lcp_pts_sf %>%
  st_cast("MULTIPOINT") %>%
  st_cast("POINT")

lcp_pts_sp <- st_coordinates(lcp_pts_sp)

# Definir la ventana (clase owin) para el análisis de densidad
window <- as(pol_sp, "owin")

# Crear un objeto de tipo ppp a partir de los puntos
lcp_pts_pp <- ppp(x = lcp_pts_sp[,1], y = lcp_pts_sp[,2], 
                  window = window)

# Calcular la densidad de rutas óptimas
lcp_density <- density(lcp_pts_pp, sigma = 1200, eps = 100)

## Guardar la covariable como geotiff
lcp_density_rast <- raster(lcp_density)
writeRaster(lcp_density_rast,
            file.path(targetdir, "densidad_rutas_optimas.tiff"),
            format = "GTiff", 
            overwrite = TRUE)
