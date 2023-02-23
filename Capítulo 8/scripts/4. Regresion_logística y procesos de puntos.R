---
título: "Capítulo 8. Caso práctico. Patrones de localización en el poblamiento de la Edad del hierro de Galicia (Noroeste de la Península Ibérica), perteneciente al libro Carrero-Pazos, M. (2022). *Arqueología Computacional del Territorio. Métodos y técnicas para estudiar decisiones humanas en paisajes pretéritos*. Oxford. Archaeopress."
autor: "Miguel Carrero-Pazos"
date: '2023-02-10'
---

## Regresión logística y procesos de puntos

# El siguiente código está basado en Spencer y Bevan 2018,
# Carrero Pazos, Bevan y Lake, 2019 y Riris 2020.
# Links a los artículos y material supplementario:
## https://discovery.ucl.ac.uk/id/eprint/10055556/
### R code: https://github.com/cspencer905/SpencerBevan-MirabelloModelling
## https://discovery.ucl.ac.uk/id/eprint/10077897/
### R code: https://github.com/MCarreroPazos/MontePenide
## https://staffprofiles.bournemouth.ac.uk/display/journal-article/332319
## R code: attached as supplementary material on the published version

# Cargar los paquetes de trabajo
spatpack<-c("raster","spatstat","rgdal","maptools",
            "sf", "MASS", "RColorBrewer", "virtualspecies",
            "onpoint", "cowplot", "ggplot2")
lapply(spatpack, require, character.only=TRUE)

# Establecer el directorio de datos
sourcedir <- file.path("data","raw_data")
targetdir <- file.path("data","derived_data")
if(!dir.exists(targetdir)) dir.create(targetdir)

# Cargar el vector del área de estudio y los yacimientos 
castros <- sf::st_read(file.path(sourcedir,
                                 "castros.shp"),
                       layer = "castros")
castros <- as(castros, "Spatial")

# Crear área de estudio
bb <- st_as_sfc(st_bbox(castros))
bb_sp <- as(bb, "Spatial")
# Convertir el vector del área de estudio en la ventana de análisis (para Spatstat)
bb_w <- as.owin(bb_sp)

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

# Comprobar variables colineales
covs <- stack(dem, dens_rutas, i_prom_topo, prom_visual)
removeCollinearity(covs, multicollinearity.cutoff = 0.75,
                   select.variables = FALSE, sample.points = TRUE, nb.points = 100000,
                   plot = TRUE) # No hay variables colineales

# Convertir las covariables a imagen para Spatstat
dem_im <- as.im(dem)
dens_rutas_im <- as.im(dens_rutas)
i_prom_topo_im <- as.im(i_prom_topo)
prom_visual_im <- as.im(prom_visual)

# Modelo de regresión logística de primer orden
## Definir los efectos de primer orden
covar_list_mod1 <- list(i_prom_topo_im, prom_visual_im, dens_rutas_im, dem_im)
names(covar_list_mod1) <- c("i_prom_topo_im", "prom_visual_im", "dens_rutas_im", "dem_im")
mod1_efectos_primer_orden <- ~ i_prom_topo_im + prom_visual_im + dens_rutas_im + dem_im

## Crear los modelos de simulación estadística
# Modelo 1 (efectos de primer orden)
mod1 <- step(ppm(sppp_castros, trend = mod1_efectos_primer_orden, 
                 interaction = NULL, 
                 covariates = covar_list_mod1, method = "logi"))

summary(mod1)

# Superficie de intensidad (véase summary(mod1))
logodds <- raster(-1.748041e+01+(i_prom_topo_im*1.183857e+00)+(prom_visual_im*3.198531e-01))
plot(logodds, main="Superficie de intensidad\npredictiva (log-odds)", 
     col = terrain.colors(8), axes=FALSE, box=FALSE, 
     legend.args=list(text="",side=2), horizontal=TRUE)

# Modelo 2 (efectos de primer y segundo orden)
# Definir la superficie de tendencia de primer orden (eliminar la variable dem -véase summary(mod1))
covar_list_mod2 <- list(i_prom_topo_im, prom_visual_im)
names(covar_list_mod2) <- c("i_prom_topo_im", "prom_visual_im")
mod2_efectos_primer_orden <- ~ i_prom_topo_im + prom_visual_im

# Efectos de segundo orden: interacción AreaInter (r = 1000)
mod2_areaInter <- step(ppm(sppp_castros, trend = mod2_efectos_primer_orden, 
                           interaction = AreaInter(r=1000),
                           covariates = covar_list_mod2))
# Efectos de segundo orden: interacción Hardcore (hc = 200)
mod2_hardcore <- step(ppm(sppp_castros, trend = mod2_efectos_primer_orden, 
                          interaction = Hardcore(hc = 200),
                          covariates = covar_list_mod2))

# Calcular las funciones de correlación par

## Modelo con efectos de primer orden
Pcfinhom_mod1 <- envelope(mod1, fun = pcfinhom, correction = "best", nsim = 10000)
plot(Pcfinhom_mod1, xlim=c(0,5000),ylim=c(0, 5), legend=FALSE, 
     main = "Efectos de primer orden")

## Modelo con efectos de primer y segundo orden (AreaInter r = 1000)
Pcfinhom_mod2A <- envelope(mod2_areaInter, fun= pcfinhom, correction="best", nsim = 10000)
plot(Pcfinhom_mod2A, xlim=c(0,10000),ylim=c(0, 5), legend=FALSE)

## Modelo con efectos de primer y segundo orden (Hardcore hc = 198m)
Pcfinhom_mod2B <- envelope(mod2_hardcore, fun= pcfinhom, correction="best", nsim = 10000)
plot(Pcfinhom_mod2B, xlim=c(0,10000),ylim=c(0, 5), legend=FALSE)

# Guardar los resultados de la función de correlación par
# tempdir <- "data/derived_data/funciones_correlacion_par/"
# save(Pcfinhom_mod1, file = file.path(tempdir, "Pcfinhom_mod1.RData"))
# save(Pcfinhom_mod2A, file = file.path(tempdir, "Pcfinhom_mod2A.RData"))
# save(Pcfinhom_mod2B, file = file.path(tempdir, "Pcfinhom_mod2B.RData"))

# Cargar los resultados de la función de correlación par
tempdir <- "data/derived_data/funciones_correlacion_par/"
load(file=file.path(tempdir,"Pcfinhom_mod1.RData"))
load(file=file.path(tempdir,"Pcfinhom_mod2A.RData"))
load(file=file.path(tempdir,"Pcfinhom_mod2B.RData"))

# Comprobar los AICs de los diferentes modelos
AIC(mod1) # 7656.264
AIC(mod2_areaInter) # 7155.866
AIC(mod2_hardcore) # 7470.05

# Crear los gráficos para las funciones de correlación par
plot1 <- plot_quantums(Pcfinhom_mod1, 
                       xlab = "Distancia entre puntos (r)",
                       ylab = "ginhom(r)",
                       labels=c("Patrón agrupado\nsignificativo", 
                                "Patrón aleatorio", 
                                "Patrón regular\nsignificativo"),
                       quantum_size= 3,
                       quantum_position = 0,
                       quantum=T,
                       base_size = 10,
                       title="a. Efectos de primer orden") +
                       scale_x_continuous(limits = c(0, 8000))
plot2 <- plot_quantums(Pcfinhom_mod2A, 
                       xlab = "Distancia entre puntos (r)",
                       ylab = "ginhom(r)",
                       labels=c("Patrón agrupado\nsignificativo", 
                                "Patrón aleatorio", 
                                "Patrón regular\nsignificativo"),
                       quantum_size= 0.5,
                       quantum_position = 0,
                       quantum=T,
                       base_size = 10,
                       title="b. Efectos de primer y segundo orden\n Interacción AreaInter (r = 1000)",
                       legend_position = "none") +
                       scale_x_continuous(limits = c(0, 8000))
plot3 <- plot_quantums(Pcfinhom_mod2B, 
                       xlab = "Distancia entre puntos (r)",
                       ylab = "ginhom(r)",
                       labels=c("Patrón agrupado\nsignificativo", 
                                "Patrón aleatorio", 
                                "Patrón regular\nsignificativo"),
                       quantum_size= 0.5,
                       quantum_position = 0,
                       quantum=T,
                       base_size = 10,
                       title="c. Efectos de primer y segundo orden\n Interacción Hardcore (hc = 200)",
                       legend_position = "none") +
                       scale_x_continuous(limits = c(0, 8000))
ggdraw() +
  draw_plot(plot1, 0, .5, 1, .5) +
  draw_plot(plot2, 0, 0, .5, .5) +
  draw_plot(plot3, .5, 0, .5, .5)
