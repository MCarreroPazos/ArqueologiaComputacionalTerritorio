---
título: "Capítulo 8. Caso práctico. Patrones de localización en el poblamiento de la Edad del hierro de Galicia (Noroeste de la Península Ibérica), perteneciente al libro Carrero-Pazos, M. (2022). *Arqueología Computacional del Territorio. Métodos y técnicas para estudiar decisiones humanas en paisajes pretéritos*. Oxford. Archaeopress."
autor: "Miguel Carrero-Pazos"
date: '2023-02-10'
---

## Patrón monovariante

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

# Cargar las covariables raster y establecer el mismo csr
dem <- raster(file.path(targetdir,"dem.tif"))
crs(dem) <- "+proj=utm +zone=29 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
dem <- crop(dem, bb_sp)
dens_rutas <- raster(file.path(targetdir,"densidad_rutas_optimas.tif"))
crs(dens_rutas) <- "+proj=utm +zone=29 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
dens_rutas <- projectRaster(dens_rutas, dem)
i_prom_topo <- raster(file.path(targetdir,"indice_posicion_topografica.tif"))
crs(i_prom_topo) <- "+proj=utm +zone=29 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
i_prom_topo <- projectRaster(i_prom_topo, dem)
prom_visual <- raster(file.path(targetdir,"prominencia_visual.tiff"))
crs(prom_visual) <- "+proj=utm +zone=29 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
prom_visual <- projectRaster(prom_visual, dem)


# Convertir las covariables a imagen para la función rhohat
dem_im <- as.im(dem)
dens_rutas_im <- as.im(dens_rutas)
i_prom_topo_im <- as.im(i_prom_topo)
prom_visual_im <- as.im(prom_visual)

# Calcular la intensidad de los castros con respecto
# a las diferentes covariables
dem.rh <- rhohat(sppp_castros, dem_im, confidence=0.95)
dens_rutas.rh <- rhohat(sppp_castros, dens_rutas_im, confidence=0.95)
i_prom_topo.rh <- rhohat(sppp_castros, i_prom_topo_im, confidence=0.95)
prom_visual.rh <- rhohat(sppp_castros, prom_visual_im, confidence=0.95)

# Crear los gráficos de intensidad para las covariables
par(mfrow=c(2,2))
par(mar=c(4, 4, 2, 1)) # c(abajo, izquierda, derecha, arriba)
plot(dem.rh, main="", xlab="Altitud (m.s.n.m.)", ylab="",
     legend=FALSE, cex.axis=0.9)
legend("topleft", legend="a. Altitud", cex=0.9, bty='n', 
       text.font=2)
plot(dens_rutas.rh, main="", xlab="Densidad de rutas", ylab="", 
     legend=FALSE, cex.axis=0.9)
legend("topleft", legend="b. Densidad de rutas de tránsito potencial", 
       cex=0.9, bty='n', text.font=2)
plot(i_prom_topo.rh, main="", xlab="Prominencia topográfica", 
     ylab="", legend=FALSE, cex.axis=0.9)
legend("topleft", legend="c. Índice de prominencia topográfica",
       cex=0.9, bty='n', text.font=2)
plot(prom_visual.rh, main="", xlab="Visibilidad", ylab="", 
     legend=FALSE, cex.axis=0.9)
legend("topleft", legend="d. Prominencia visual", 
       cex=0.9, bty='n', text.font=2)
par(mfrow=c(1,1))

# Predecir la densidad de castros dependiendo de las covariables

par(mfrow=c(2,2))
par(mar=c(2, 2, 2, 2)) # c(abajo, izquierda, derecha, arriba)
# Altitud
pred_dem.rh <- predict(dem.rh)
plot(pred_dem.rh, las=1, main = "Altitud")
# Densidad de rutas óptimas
pred_dens_rutas.rh <- predict(dens_rutas.rh)
plot(pred_dens_rutas.rh, las=1, main = "Densidad de tránsito")
# Prominencia topográfica
pred_i_prom_topo.rh <- predict(i_prom_topo.rh)
cl   <- interp.colours(c("lightyellow", "orange" ,"red"), 100)
plot(pred_i_prom_topo.rh, col=cl, las=1, main = "Prominencia topográfica", gamma = 0.25)
# Prominencia visual
pred_prom_visual.rh <- predict(prom_visual.rh)
plot(pred_prom_visual.rh, col=cl, las=1, main = "Prominencia visual", gamma = 0.25)
par(mfrow=c(1,1))

