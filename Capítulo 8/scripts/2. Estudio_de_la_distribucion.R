---
título: "Capítulo 8. Caso práctico. Patrones de localización en el poblamiento de la Edad del hierro de Galicia (Noroeste de la Península Ibérica), perteneciente al libro Carrero-Pazos, M. (2022). *Arqueología Computacional del Territorio. Métodos y técnicas para estudiar decisiones humanas en paisajes pretéritos*. Oxford. Archaeopress."
autor: "Miguel Carrero-Pazos"
date: '2023-02-10'
---

## Estudio de la distribución: la intensidad espacial

# Cargar los paquetes de trabajo
spatpack<-c("raster","spatstat","rgdal","maptools",
            "sf")
lapply(spatpack, require, character.only=TRUE)

# Establecer el directorio de datos
sourcedir <- file.path("data","raw_data")
targetdir <- file.path("data","derived_data")
if(!dir.exists(targetdir)) dir.create(targetdir)

# Definir crs
crs = "+proj=utm +zone=29 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# Cargar el vector de los yacimientos 
castros <- sf::st_read(file.path(sourcedir,
                                 "castros.shp"),
                       layer = "castros")
castros <- sf::st_transform(castros, 25829)

# Crear área de estudio
bb <- st_as_sfc(st_bbox(castros))
bb_sp <- as(bb, "Spatial")
# Convertir el vector del área de estudio en la ventana de análisis (para Spatstat)
bb_w <- as.owin(bb)

# Convertir el shapefile de los castros en un patrón de puntos
sppp_castros <- ppp(x = castros$xcoord, 
                    y = castros$ycoord, 
                    window = bb_w)

dem <- raster(file.path(targetdir,
                        "dem.tif"))
raster::crs(dem) <- raster::crs(crs)

# Comprobar Complete Spatial Randomness
# Función K
Kfunct_castros <- Kest(sppp_castros, correction = "Ripley")
plot(Kfunct_castros, xlim = c(0,6000), main = "Función K de Ripley")

# Densidad de castros
par(mfrow = c(1, 3))
par(mar=c(2,2,2,2))
plot(density(sppp_castros, sigma = 500), main = "")
title("Sigma = 500", line=-8)
plot(density(sppp_castros, sigma = 1000), main = "")
title("Sigma = 1000", line=-8)
plot(density(sppp_castros, sigma = 2000), main = "")
title("Sigma = 2000", line=-8)
par(mfrow = c(1, 1))

# Quadrat test
qtest <- quadrat.test(sppp_castros, 
                   method = "MonteCarlo")

qtest # para observar los resultados

# Crear el gráfico
plot(density(sppp_castros, sigma = 1000, eps=50), main = "Test de cuadrantes")
plot(qtest, cex = 0.7, col = "white", add = T)



