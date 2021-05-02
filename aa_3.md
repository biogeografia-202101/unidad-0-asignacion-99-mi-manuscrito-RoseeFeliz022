Análisis de agrupamiento (cluster analysis). <br> Parte 3: Grupos (clústers), variables ambientales y mapas
================
JR
15 de noviembre, 2020

``` r
knitr::opts_chunk$set(fig.width=12, fig.height=8)
```

Preámbulo
---------

### Cargar paquetes

``` r
library(mapview)
library(tidyverse)
```

    ## ── Attaching packages ──────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.2.1     ✔ purrr   0.3.3
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.3
    ## ✔ tidyr   1.0.0     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(sf)
```

    ## Linking to GEOS 3.8.0, GDAL 3.0.4, PROJ 7.0.0

``` r
library(RColorBrewer)
source('biodata/funciones.R')
```

### Cargar datos

``` r
load('biodata/Myrtaceae.Rdata')
load('biodata/matriz_ambiental.Rdata')
grupos_upgma_k2 <- readRDS('grupos_upgma_k2.RDS')
table(grupos_upgma_k2) #Importante, tener en cuenta los desiguales tamaños de los grupos
```

    ## grupos_upgma_k2
    ##  1  2 
    ## 48  2

``` r
grupos_ward_k4 <- readRDS('grupos_ward_k4.RDS')
table(grupos_ward_k4)
```

    ## grupos_ward_k4
    ##  1  2  3  4 
    ## 20 13  2 15

### Paletas

``` r
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
```

Explorar efectos
----------------

### Pruebas de igualdad de promedios de las variables entre 2 grupos

Para evaluar homogeneidad de promedios usaré las pruebas *t* (medias), basada en la distribución *t* de *Student*, y la prueba no paramétrica de la suma de rangos de Wilcoxon (medianas), usando como variable de agrupamiento los grupos establecidos en el agrupamiento UPGMA. Nota que en mi caso UPGMA clasifica los sitios en dos grupos, pero en tu caso podría ser distinto (para evaluar homogeneidad de promedios de un número mayor de grupos, ver sección siguiente).

Primero crearé un objeto que permita realizar tanto las pruebas como los diagramas de cajas.

``` r
(m_amb_upgma_k2 <- bci_env_grid %>%
    select_if(is.numeric) %>% select(-id) %>% 
    mutate(grupos_upgma_k2) %>%
    st_drop_geometry() %>% 
    pivot_longer(-grupos_upgma_k2, names_to = "variable", values_to = "valor"))
```

    ## # A tibble: 1,650 x 3
    ##    grupos_upgma_k2 variable                       valor
    ##    <fct>           <chr>                          <dbl>
    ##  1 1               heterogeneidad_ambiental       0.627
    ##  2 1               UTM.EW                    625754.   
    ##  3 1               UTM.NS                   1011569.   
    ##  4 1               geomorf_llanura_pct           10.0  
    ##  5 1               geomorf_pico_pct               0    
    ##  6 1               geomorf_interfluvio_pct        0.83 
    ##  7 1               geomorf_hombrera_pct          10.8  
    ##  8 1               geomorf_espolón/gajo_pct       7.26 
    ##  9 1               geomorf_vertiente_pct         67.1  
    ## 10 1               geomorf_vaguada_pct            3.28 
    ## # … with 1,640 more rows

A continuación, las pruebas:

``` r
m_amb_upgma_k2 %>%
  group_by(variable) %>%
  summarise(
    p_valor_t = t.test(valor ~ grupos_upgma_k2)$p.value,
    p_valor_w = wilcox.test(valor ~ grupos_upgma_k2, exact = F)$p.value) %>%
  arrange(p_valor_t) %>%
  print(n=Inf)
```

    ## # A tibble: 33 x 3
    ##    variable                      p_valor_t p_valor_w
    ##    <chr>                             <dbl>     <dbl>
    ##  1 geomorf_espolón/gajo_pct   0.0000000202    0.0475
    ##  2 geomorf_vaguada_pct        0.0000000260    0.0351
    ##  3 Zn                         0.000000139     0.119 
    ##  4 orientacion_media          0.000000148     0.0788
    ##  5 UTM.NS                     0.00000841      0.325 
    ##  6 geomorf_valle_pct          0.0000405       0.110 
    ##  7 B                          0.0000450       0.131 
    ##  8 N                          0.000151        0.310 
    ##  9 Ca                         0.000638        0.131 
    ## 10 geomorf_interfluvio_pct    0.000718        0.181 
    ## 11 elevacion_media            0.00126         0.225 
    ## 12 heterogeneidad_ambiental   0.0155          0.319 
    ## 13 pH                         0.0164          0.0354
    ## 14 pendiente_media            0.0203          0.0354
    ## 15 Mg                         0.0265          0.190 
    ## 16 riqueza_global             0.0279          0.0666
    ## 17 Al                         0.0465          0.0505
    ## 18 UTM.EW                     0.0522          0.332 
    ## 19 Cu                         0.0536          0.310 
    ## 20 Fe                         0.0766          0.674 
    ## 21 geomorf_hombrera_pct       0.149           0.220 
    ## 22 K                          0.152           0.144 
    ## 23 geomorf_pico_pct           0.159           0.764 
    ## 24 curvatura_tangencial_media 0.251           0.386 
    ## 25 Mn                         0.309           0.414 
    ## 26 geomorf_sima_pct           0.322           0.919 
    ## 27 geomorf_llanura_pct        0.349           0.0321
    ## 28 geomorf_vertiente_pct      0.396           0.0505
    ## 29 abundancia_global          0.422           0.287 
    ## 30 geomorf_piedemonte_pct     0.572           0.294 
    ## 31 P                          0.634           0.443 
    ## 32 curvatura_perfil_media     0.869           0.980 
    ## 33 N.min.                     0.897           0.901

Interesa observar las variables que obtuvieron valores de p&lt;0.01. Reitero que, en mi caso, mis grupos resultaron muy desiguales, recordando: el grupo 1 tiene 43 sitios (43) y el grupo 2 tiene 7. Este desigual número de sitios por grupo, hace que la prueba estadística pierda potencia, porque se viola la recomendación de evitar tamaños de los tratamientos muy desiguales.

Por otra parte, este es un buen momento para "revisitar" tus análisis exploratorios de datos (AED), específicamente el análisis de correlación (*script* 5). Es probable que algunas de las variables ambientales que presentaron efecto entre grupos (las que obtuvieron p&lt;0.01), te aparezca también como significativamente correlacionada con la abundancia o la riqueza en el script 5 de AED.

Los gráficos:

``` r
m_amb_upgma_k2 %>% 
  group_by(variable) %>% 
  ggplot() + aes(x = grupos_upgma_k2, y = valor, fill = grupos_upgma_k2) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = 'Accent') +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~ variable, scales = 'free_y')
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-7-1.png)

Mapas:

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
mapa_upgma_k2 %>% mapshot(
  file = 'mapa_upgma_k2.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

Mapa de una de las variables donde se presentó efecto de su promedio (p&lt;0.01), en este caso, Zinc (`Zn`)

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
mapa_zn %>% mapshot(
  file = 'mapa_zinc.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

### Pruebas de igualdad de promedios de las variables entre 3 grupos o más

Objeto común:

``` r
(m_amb_ward_k4 <- bci_env_grid %>%
    select_if(is.numeric) %>% select(-id) %>% 
    mutate(grupos_ward_k4) %>%
    st_drop_geometry() %>% 
    pivot_longer(-grupos_ward_k4, names_to = "variable", values_to = "valor"))
```

    ## # A tibble: 1,650 x 3
    ##    grupos_ward_k4 variable                       valor
    ##    <fct>          <chr>                          <dbl>
    ##  1 1              heterogeneidad_ambiental       0.627
    ##  2 1              UTM.EW                    625754.   
    ##  3 1              UTM.NS                   1011569.   
    ##  4 1              geomorf_llanura_pct           10.0  
    ##  5 1              geomorf_pico_pct               0    
    ##  6 1              geomorf_interfluvio_pct        0.83 
    ##  7 1              geomorf_hombrera_pct          10.8  
    ##  8 1              geomorf_espolón/gajo_pct       7.26 
    ##  9 1              geomorf_vertiente_pct         67.1  
    ## 10 1              geomorf_vaguada_pct            3.28 
    ## # … with 1,640 more rows

Pruebas, en este caso ANOVA (evalúa homogeneidad de medias; no se cumplen muchos de los supuestos requeridos para esta prueba) y Kruskal-Wallis (evalúa homogeneidad de medianas):

``` r
m_amb_ward_k4 %>% 
  group_by(variable) %>% 
  summarise(
    p_valor_a = oneway.test(valor ~ grupos_ward_k4)$p.value,
    p_valor_k = kruskal.test(valor ~ grupos_ward_k4)$p.value) %>%
  arrange(p_valor_k) %>%
  print(n=Inf)
```

    ## # A tibble: 33 x 3
    ##    variable                      p_valor_a p_valor_k
    ##    <chr>                             <dbl>     <dbl>
    ##  1 pH                           0.000669     0.00309
    ##  2 Cu                           0.0178       0.0144 
    ##  3 UTM.EW                       0.0147       0.0231 
    ##  4 K                            0.0865       0.0282 
    ##  5 B                            0.0000313    0.0472 
    ##  6 Zn                           0.00000343   0.0494 
    ##  7 pendiente_media              0.00163      0.0514 
    ##  8 N                            0.0000539    0.0608 
    ##  9 geomorf_vaguada_pct        NaN            0.0741 
    ## 10 N.min.                       0.403        0.0888 
    ## 11 geomorf_interfluvio_pct    NaN            0.0930 
    ## 12 orientacion_media            0.0000513    0.110  
    ## 13 riqueza_global               0.00642      0.120  
    ## 14 Ca                           0.000841     0.122  
    ## 15 Mn                           0.183        0.139  
    ## 16 geomorf_llanura_pct          0.493        0.142  
    ## 17 geomorf_vertiente_pct        0.506        0.151  
    ## 18 Mg                           0.0296       0.158  
    ## 19 elevacion_media              0.000900     0.171  
    ## 20 Al                           0.0168       0.179  
    ## 21 UTM.NS                     NaN            0.186  
    ## 22 geomorf_espolón/gajo_pct     0.0000144    0.189  
    ## 23 Fe                           0.0939       0.196  
    ## 24 abundancia_global            0.378        0.218  
    ## 25 geomorf_piedemonte_pct       0.341        0.302  
    ## 26 geomorf_valle_pct          NaN            0.315  
    ## 27 curvatura_perfil_media       0.560        0.420  
    ## 28 geomorf_hombrera_pct         0.149        0.444  
    ## 29 curvatura_tangencial_media   0.285        0.581  
    ## 30 P                            0.702        0.637  
    ## 31 geomorf_pico_pct           NaN            0.670  
    ## 32 geomorf_sima_pct           NaN            0.682  
    ## 33 heterogeneidad_ambiental     0.0142       0.753

Gráficos:

``` r
m_amb_ward_k4 %>% 
  group_by(variable) %>% 
  ggplot() + aes(x = grupos_ward_k4, y = valor, fill = grupos_ward_k4) + 
  geom_boxplot() + 
  scale_fill_brewer(palette = 'Accent') +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~ variable, scales = 'free_y')
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-12-1.png)

Mapas:

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
mapa_ward_k4 %>% mapshot(
  file = 'mapa_ward_k4.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

Mapa de una de las variables donde se presentó efecto de su promedio (p&lt;0.01), en este caso, Zinc (`Zn`)

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
mapa_abundanciagl %>% mapshot(
  file = 'mapa_abundanciagl.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
mapa_curvatura %>% mapshot(
  file = 'mapa_curvatura.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-16-1.png)

``` r
mapa_Al %>% mapshot(
  file = 'mapa_Al.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
mapa_Fe %>% mapshot(
  file = 'mapa_Fe.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-18-1.png)

``` r
mapa_curva %>% mapshot(
  file = 'mapa_curva.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
mapa_elevacion %>% mapshot(
  file = 'mapa_elevacion.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
mapa_Mn %>% mapshot(
  file = 'mapa_Mn.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-21-1.png)

``` r
mapa_N.min. %>% mapshot(
  file = 'mapa_Nin.min.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
mapa_P %>% mapshot(
  file = 'mapa_P.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-23-1.png)

``` r
mapa_heterogeneidadambiental %>% mapshot(
  file = 'mapa_heterogeneidad.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-24-1.png)

``` r
mapa_orientacionmedia %>% mapshot(
  file = 'mapa_orientacion.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-25-1.png)

``` r
mapa_riquezaglobal %>% mapshot(
  file = 'mapa_riquezaglobal.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-26-1.png)

``` r
mapa_geohombrera %>% mapshot(
  file = 'mapa_geohombrera.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
mapa_Cu
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-27-1.png)

``` r
mapa_Cu %>% mapshot(
  file = 'mapa_Cu.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```

``` r
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
```

![](aa_3_files/figure-markdown_github/unnamed-chunk-28-1.png)

``` r
mapa_utm_ns %>% mapshot(
  file = 'mapa_utm_ns.png', 
  remove_controls = c("zoomControl", "layersControl", "homeButton")
)
```
