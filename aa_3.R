#' ---
#' title: "Análisis de agrupamiento (cluster analysis). <br> Parte 3: Grupos (clústers), variables ambientales y mapas"
#' author: "JR"
#' date: "15 de noviembre, 2020"
#' output: github_document
#' ---

knitr::opts_chunk$set(fig.width=12, fig.height=8)

#' ## Preámbulo
#' 
#' ### Cargar paquetes
#' 
library(mapview)
library(tidyverse)
library(sf)
library(RColorBrewer)
source('biodata/funciones.R')
#' 
#' ### Cargar datos
#' 
load('biodata/Myrtaceae.Rdata')
load('biodata/matriz_ambiental.Rdata')
grupos_upgma_k2 <- readRDS('grupos_upgma_k2.RDS')
table(grupos_upgma_k2) #Importante, tener en cuenta los desiguales tamaños de los grupos
grupos_ward_k4 <- readRDS('grupos_ward_k4.RDS')
table(grupos_ward_k4)
#' 
#' ### Paletas
#' 
rojo <- colorRampPalette(brewer.pal(8, "Reds"))
rojo_inv <- colorRampPalette(rev(brewer.pal(8, "Reds")))
colores_grupos <- brewer.pal(8, "Accent")
azul <- colorRampPalette(brewer.pal(8, "Blues"))
azul_inv <- colorRampPalette(rev(brewer.pal(8, "Blues")))
verde <- colorRampPalette(brewer.pal(8, "Greens"))
verde_inv <- colorRampPalette(rev(brewer.pal(8, "Greens")))
gris <- colorRampPalette(brewer.pal(8, "Greys"))
gris_inv <- colorRampPalette(rev(brewer.pal(8, "Greys")))
naranja <- colorRampPalette(brewer.pal(8, "Oranges"))
naranja_inv <- colorRampPalette(rev(brewer.pal(8, "Oranges")))
pastel1 <- colorRampPalette(brewer.pal(8, "YlOrBr"))
pastel1_inv <- colorRampPalette(rev(brewer.pal(8, "YlOrBr")))
pastel2 <- colorRampPalette(brewer.pal(8, "YlOrRd"))
pastel2_inv <- colorRampPalette(rev(brewer.pal(8, "YlOrRd")))
PRGn <- colorRampPalette(brewer.pal(8, "PRGn"))
PRGn_inv <- colorRampPalette(rev(brewer.pal(8, "PRGn")))
PuBu <- colorRampPalette(brewer.pal(8, "PuBu"))
PuBu_inv <- colorRampPalette(rev(brewer.pal(8, "PuBu")))
#' 
#' ## Explorar efectos
#' 
#' ### Pruebas de igualdad de promedios de las variables entre 2 grupos
#' 
#' Para evaluar homogeneidad de promedios usaré las pruebas *t* (medias), basada en la distribución *t* de *Student*, y la prueba no paramétrica de la suma de rangos de Wilcoxon (medianas), usando como variable de agrupamiento los grupos establecidos en el agrupamiento UPGMA. Nota que en mi caso UPGMA clasifica los sitios en dos grupos, pero en tu caso podría ser distinto (para evaluar homogeneidad de promedios de un número mayor de grupos, ver sección siguiente).
#' 
#' Primero crearé un objeto que permita realizar tanto las pruebas como los diagramas de cajas.
#' 
(m_amb_upgma_k2 <- bci_env_grid %>%
    select_if(is.numeric) %>% select(-id) %>% 
    mutate(grupos_upgma_k2) %>%
    st_drop_geometry() %>% 
    pivot_longer(-grupos_upgma_k2, names_to = "variable", values_to = "valor"))
#' 
#' A continuación, las pruebas:
#' 
m_amb_upgma_k2 %>%
  group_by(variable) %>%
  summarise(
    p_valor_t = t.test(valor ~ grupos_upgma_k2)$p.value,
    p_valor_w = wilcox.test(valor ~ grupos_upgma_k2, exact = F)$p.value) %>%
  arrange(p_valor_t) %>%
  print(n=Inf)
#' 
#' Interesa observar las variables que obtuvieron valores de p<0.01. Reitero que, en mi caso, mis grupos resultaron muy desiguales, recordando: el grupo 1 tiene 43 sitios (43) y el grupo 2 tiene 7. Este desigual número de sitios por grupo, hace que la prueba estadística pierda potencia, porque se viola la recomendación de evitar tamaños de los tratamientos muy desiguales.
#' 
#' Por otra parte, este es un buen momento para "revisitar" tus análisis exploratorios de datos (AED), específicamente el análisis de correlación (*script* 5). Es probable que algunas de las variables ambientales que presentaron efecto entre grupos (las que obtuvieron p<0.01), te aparezca también como significativamente correlacionada con la abundancia o la riqueza en el script 5 de AED.
#' 
#' Los gráficos:
#' 
m_amb_upgma_k2 %>% 
  group_by(variable) %>% 
  ggplot() + aes(x = grupos_upgma_k2, y = valor, fill = grupos_upgma_k2) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = 'Accent') +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~ variable, scales = 'free_y')
#' 
#' Mapas:
#' 
mapa_upgma_k2 <- mapView(
  bci_env_grid %>% mutate(grupos_upgma_k2),
  layer.name = 'Grupos (2) UPGMA',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = colores_grupos[1:2],
  zcol = 'grupos_upgma_k2') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_upgma_k2
mapa_upgma_k2 %>% mapshot(
  file = 'mapa_upgma_k2.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#' 
#' Mapa de una de las variables donde se presentó efecto de su promedio (p<0.01), en este caso, Zinc (`Zn`)
#' 
mapa_zn <- mapView(
  bci_env_grid,
  layer.name = 'Zinc',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = rojo,
  zcol = 'Zn') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_zn
mapa_zn %>% mapshot(
  file = 'mapa_zinc.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#' 
#' ### Pruebas de igualdad de promedios de las variables entre 3 grupos o más
#' 
#' Objeto común:
#' 
(m_amb_ward_k4 <- bci_env_grid %>%
    select_if(is.numeric) %>% select(-id) %>% 
    mutate(grupos_ward_k4) %>%
    st_drop_geometry() %>% 
    pivot_longer(-grupos_ward_k4, names_to = "variable", values_to = "valor"))
#' 
#' Pruebas, en este caso ANOVA (evalúa homogeneidad de medias; no se cumplen muchos de los supuestos requeridos para esta prueba) y Kruskal-Wallis (evalúa homogeneidad de medianas):
#' 
m_amb_ward_k4 %>% 
  group_by(variable) %>% 
  summarise(
    p_valor_a = oneway.test(valor ~ grupos_ward_k4)$p.value,
    p_valor_k = kruskal.test(valor ~ grupos_ward_k4)$p.value) %>%
  arrange(p_valor_k) %>%
  print(n=Inf)
#' 
#' Gráficos:
#' 
m_amb_ward_k4 %>% 
  group_by(variable) %>% 
  ggplot() + aes(x = grupos_ward_k4, y = valor, fill = grupos_ward_k4) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = 'Accent') +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~ variable, scales = 'free_y')
#' 
#' Mapas:
#' 
mapa_ward_k4 <- mapView(
  bci_env_grid %>% mutate(grupos_ward_k4),
  layer.name = 'Grupos (4) Ward',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = colores_grupos[1:4],
  zcol = 'grupos_ward_k4') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_ward_k4
mapa_ward_k4 %>% mapshot(
  file = 'mapa_ward_k4.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#' 
#' Mapa de una de las variables donde se presentó efecto de su promedio (p<0.01), en este caso, Zinc (`Zn`)
#'
mapa_abundanciagl <- mapView(
  bci_env_grid,
  layer.name = 'abundancia_global',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = azul_inv,
  zcol = 'abundancia_global') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_abundanciagl
mapa_abundanciagl %>% mapshot(
  file = 'mapa_abundanciagl.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_curvatura <- mapView(
  bci_env_grid,
  layer.name = 'curvatura_perfil_media',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = rojo_inv,
  zcol = 'curvatura_perfil_media') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_curvatura
mapa_curvatura %>% mapshot(
  file = 'mapa_curvatura.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_Al <- mapView(
  bci_env_grid,
  layer.name = 'Al',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = azul_inv,
  zcol = 'Al') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_Al
mapa_Al %>% mapshot(
  file = 'mapa_Al.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_Fe <- mapView(
  bci_env_grid,
  layer.name = 'Fe',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = verde_inv,
  zcol = 'Fe') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_Fe
mapa_Fe %>% mapshot(
  file = 'mapa_ph.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_curva <- mapView(
  bci_env_grid,
  layer.name = 'curvatura_tangencial_media',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = gris_inv,
  zcol = 'curvatura_tangencial_media') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_curva
mapa_curva %>% mapshot(
  file = 'mapa_curva.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_elevacion <- mapView(
  bci_env_grid,
  layer.name = 'elevacion_media',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = naranja_inv,
  zcol = 'elevacion_media') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_elevacion
mapa_elevacion %>% mapshot(
  file = 'mapa_elevacion.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_Mn <- mapView(
  bci_env_grid,
  layer.name = 'Mn',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = pastel1_inv,
  zcol = 'Mn') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_Mn
mapa_Mn %>% mapshot(
  file = 'mapa_Mn.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_N.min. <- mapView(
  bci_env_grid,
  layer.name = 'N.min.',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = pastel2_inv,
  zcol = 'N.min.') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_N.min.
mapa_ph %>% mapshot(
  file = 'mapa_Nin.min.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_P <- mapView(
  bci_env_grid,
  layer.name = 'P',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = PuBu_inv,
  zcol = 'P') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_P
mapa_P %>% mapshot(
  file = 'mapa_P.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_heterogeneidadambiental <- mapView(
  bci_env_grid,
  layer.name = 'heterogeineidad_ambiental',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = PuBu_inv,
  zcol = 'heterogeneidad_ambiental') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_heterogeneidadambiental
mapa_heterogeneidadambiental %>% mapshot(
  file = 'mapa_heterogeneidad.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_orientacionmedia <- mapView(
  bci_env_grid,
  layer.name = 'orientacion_media',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = PuBu_inv,
  zcol = 'orientacion_media') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_orientacionmedia
mapa_orientacionmedia %>% mapshot(
  file = 'mapa_orientacion.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
mapa_riquezaglobal <- mapView(
  bci_env_grid,
  layer.name = 'riqueza_global',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = PuBu_inv,
  zcol = 'riqueza_global') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_riquezaglobal
mapa_riquezaglobal %>% mapshot(
  file = 'mapa_riquezaglobal.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_geohombrera <- mapView(
  bci_env_grid,
  layer.name = 'geomorf_hombrera_pct',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = azul_inv,
  zcol = 'geomorf_hombrera_pct') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_geohombrera
mapa_geohombrera %>% mapshot(
  file = 'mapa_geohombrera.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_Cu <- mapView(
  bci_env_grid,
  layer.name = 'Cu',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = azul_inv,
  zcol = 'Cu') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_abundanciagl
mapa_abundanciagl %>% mapshot(
  file = 'mapa_abundanciagl.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'
mapa_utm_ns <- mapView(
  bci_env_grid,
  layer.name = 'UTM.NS',
  alpha.regions = 0.6,
  map.types = 'OpenTopoMap',
  legend = T,
  col.regions = azul_inv,
  zcol = 'UTM.NS') %>%
  addStaticLabels(label = bci_env_grid$id) %>% 
  leaflet::setView(
    lng = -79.85136,
    lat = 9.15097,
    zoom = 15)
mapa_utm_ns
mapa_utm_ns %>% mapshot(
  file = 'mapa_utm_ns.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
#'