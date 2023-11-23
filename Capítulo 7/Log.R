---
título: "Capítulo 7 · Caso práctico. Modelos de procesos de puntos para el estudio del Megalitismo de A Serra do Barbanza (Noroeste de la Península Ibérica), perteneciente al libro Carrero-Pazos, M. (2022). *Arqueología Computacional del Territorio. Métodos y técnicas para estudiar decisiones humanas en paisajes pretéritos*. Oxford. Archaeopress."
autor: "Miguel Carrero-Pazos"
date: '2023-02-10'
---

# Cargar los paquetes de trabajo
spatpack<-c("raster","spatstat","rgdal","maptools", "MuMIn", "MASS")
lapply(spatpack, require, character.only=TRUE)

# Establecer el directorio de trabajo
setwd("C:/Users/Usuario/Desktop/BarbanzaVisualModels")

# Cargar informacion de inicio: yacimientos y vector de ?rea
tumulos <- read.table(file="BarbanzaVisualModels/csv/yacimientos.csv",header=TRUE, sep=";")
area_estudio <- readOGR(dsn="BarbanzaVisualModels/shp/area.shp", layer="area")

# Crear el patrón de puntos
area <- as(area_estudio,"owin")
sppp <- ppp(x=tumulos$UMTX, y=tumulos$UMTY, window=area)

# Cargar rasters generados con metodologias SIG
lcp_dens <- raster("BarbanzaVisualModels/grids/lcp_density100m.tiff")
total_view <- raster("BarbanzaVisualModels/grids/total_viewshed.tif")
tpi <- raster("BarbanzaVisualModels/grids/TPI.tif")
dist_ejes_cuencas <- raster("BarbanzaVisualModels/grids/dist_cuencas_hidro.tiff")
total_horizon <- raster("BarbanzaVisualModels/grids/total_horizon.tif")
visib_rutas <- raster("BarbanzaVisualModels/grids/visib_acum_desde_rutas.tiff")

# Convertir las variables a formato espacial (para spatstat)
lcp_densidad <- as.im(as(lcp_dens, "SpatialGridDataFrame"))
t_view <- as.im(as(total_view,"SpatialGridDataFrame"))
tprom <- as.im(as(tpi,"SpatialGridDataFrame"))
d_ejes_cuencas <- as.im(as(dist_ejes_cuencas,"SpatialGridDataFrame"))
t_horizon <- as.im(as(total_horizon,"SpatialGridDataFrame"))
visib_desde_rutas <- as.im(as(visib_rutas,"SpatialGridDataFrame"))

## Creacion de los modelos de procesos de puntos
# Nota inicial: los modelos seran creados con 999 simulaciones aleatorias, pero este numero puede modificarse cambiandolo en "nsim"

# Modelo 1
(model1 = ppm(sppp~lcp_densidad+t_view+tprom))

# Calcular la K residual del modelo 1
residualK_mod1 = envelope(model1,Kres,correction="best",nsim=999)

# Modelo 2
(model2 = ppm(sppp~d_ejes_cuencas+t_horizon+visib_desde_rutas))

# Calcular la K residual del modelo 2
residualK_mod2 = envelope(model2,Kres,correction="best",nsim=999)

# Modelo 3 (aleatorio)
(model3 = ppm(sppp~1)) # Modelo nulo

# Calcular la K residual del modelo 3
residualK_mod3 = envelope(model3,Kres,correction="best",nsim=999)

# Visualizar el criterio de Informacion Akaike para los diferentes modelos
AIC(model1)
AIC(model2)
AIC(model3)

# Crear un grafico con los tres modelos de K residuales
par(mar=c(3, 2, 3, 2), mfrow=c(1,3))
plot(residualK_mod1, legend = F, main="Modelo 1")
plot(residualK_mod2, legend = F,main="Modelo 2")
plot(residualK_mod3, legend = F,main="Modelo 3")
par(mfrow=c(1,1))

# Guardar la figura como archivo png en la carpeta "figuras"
png(file="BarbanzaVisualModels/figuras/Modelos_K_Residual.png",width = 600,height=350)
par(mar=c(3, 2, 3, 2), mfrow=c(1,3))
plot(residualK_mod1, legend = F, main="Modelo 1")
plot(residualK_mod2, legend = F,main="Modelo 2")
plot(residualK_mod3, legend = F,main="Modelo 3")
par(mfrow=c(1,1))
dev.off()
