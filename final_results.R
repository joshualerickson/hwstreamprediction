library(ggtext)
library(tidyverse)
library(ggridges)
library(wesanderson)
library(tidyverse)
source('utils.R')

# for cairo_view()

trace(grDevices::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)

gbm_topoclim_results <- read_csv('data/gbm_topoclim_results.csv')

gbm_topoclim_varimp <- read_csv('data/gbm_topoclim_varimp.csv')

gbm_topoclim_results <- gbm_topoclim_results %>% mutate(cluster_no = str_remove(cluster, 'cluster.'),
                                                        cluster_no = parse_number(cluster_no))

gbm_topoclim_varimp <- gbm_topoclim_varimp %>% mutate(cluster_no = str_remove(cluster, 'cluster.'),
                                                      cluster_no = parse_number(cluster_no))

gbm_topo_results <- read_csv('data/gbm_topo_results.csv')

gbm_topo_varimp <- read_csv('data/gbm_topo_varimp.csv')

gbm_topo_results <- gbm_topo_results %>% mutate(cluster_no = str_remove(cluster, 'cluster.'),
                                                        cluster_no = parse_number(cluster_no))

gbm_topo_varimp <- gbm_topo_varimp %>% mutate(cluster_no = str_remove(cluster, 'cluster.'),
                                                      cluster_no = parse_number(cluster_no))

nhd_results <- read_csv('data/gbm_nhd_results.csv')

nhd_results_mean <- nhd_results %>% summarise(across(is.numeric, mean))

topoclim_results_mean <- gbm_topoclim_results %>% summarise(across(is.numeric, mean)) %>% 
  mutate(model = 'topoclim') %>% 
  select(colnames(nhd_results_mean))

topo_results_mean <- gbm_topo_results %>% summarise(across(is.numeric, mean)) %>% 
  mutate(model = 'topo') %>% 
  select(colnames(nhd_results_mean))


all_mean_results <- bind_rows(topoclim_results_mean, topo_results_mean, nhd_results_mean)

# Table 2 final results ---------------------------------------------------

## Table 2 in final results
all_mean_results %>% select(c('Accuracy', 'Kappa', 'J', 'ROC', 'Specificity', 'Sensitivity')) %>% 
  kableExtra::kable() %>% kableExtra::kable_styling()

### Plotting


# Figure 5  ---------------------------------------------------------------

## figure 5 in paper  


gbm_topo_results %>% 
  bind_rows(gbm_topoclim_results) %>% 
  bind_rows(nhd_results %>% mutate(model = 'NHDPlus')) %>% 
  dplyr::select(Accuracy, J, Kappa, ROC, Sensitivity, Specificity, cluster_no, model) %>% 
  pivot_longer(1:6) %>% 
  ggplot(aes(value)) + 
  geom_density(aes(linetype = model, color = model),
               key_glyph = 'path',
               bw = 0.005,
               size = .75) +
  geom_rug(aes(color=model)) + 
  scale_linetype_manual(values = c(1,2,3),
                        labels = c('Topo-only','Topoclimatic','NHDPlus')) +
  scale_color_manual(values = c('black', '#D81B60', '#1E88E5'),
                     labels = c('Topo-only','Topoclimatic','NHDPlus')) + 
guides(
  linetype = guide_legend(keywidth = 2.5),
  color = guide_legend(keywidth = 2.5,
                      override.aes = list(
                        linetype = c(1, 2, 3),
                        size = c(.75, .75, .75),
                      reverse = T))) + 
  labs(color = 'Model Type', linetype = 'Model Type', x = '') +
  facet_wrap(~name, scales = 'free') +
  custom_theme()


# Figures 7 and 8 ---------------------------------------------------------

### figure 7 and 8
library(pdp)
library(vip)
library(gbm)
# this is a tuned GBM model being used with 4 variables and 10 cluster hold out

GBM_k10_4_model <- readRDS("data/GBM_k10_4_model.rds")

inter <- pdp::partial(GBM_k10_4_model, pred.var = c("accum30","deficitRS"), trim.outliers = TRUE, prob = TRUE)
def_npp <- pdp::partial(GBM_k10_4_model, pred.var = c("deficitRS","nppmmid30agg"), trim.outliers = TRUE, prob = TRUE) 

inter <- inter %>% dplyr::mutate(yhat = 1-yhat)
pdp_uaa_def <- plotPartial(inter,
                    levelplot = FALSE,
                    plot.pdp = F,
                    zlab = "yhat",
                    drape=T,
                    colorkey = T,
                    contour = T, 
                    main = "Topoclimatic GBM Interactions (CWD and UAA)",
                    screen = list(z = 120, x = -50, y = 0),
                    xlab = 'UAA',
                    col.regions = rev(hcl.colors(100, 'Grays')),
                    smooth = F,
                    ylab = 'CWD')
pdp_uaa_def

def_npp <- def_npp %>% dplyr::mutate(yhat = 1-yhat)

pdp_def_npp <- plotPartial(def_npp,
                    levelplot = FALSE,
                    zlab = "yhat",
                    drape=T,
                    colorkey = T,
                    contour = TRUE, 
                    main = "Topoclimatic GBM Interactions (CWD and NPP)",
                    screen = list(z = -10, x = -75, y = -30),
                    col.regions = rev(hcl.colors(100, 'Grays')),
                    xlab = 'CWD',
                    ylab = 'NPP')
pdp_def_npp

#now make the single pdp plots 

pdp_uaa <- pdp::partial(GBM_k10_4_model,
                        pred.var = "accum30",
                        trim.outliers = TRUE,
                        prob = TRUE,
                        rug = TRUE)
pdp_uaa <- pdp_uaa %>% dplyr::mutate(Variable = 'UAA') %>% 
           dplyr::rename(value = 'accum30')

pdp_npp <- pdp::partial(GBM_k10_4_model,
                        pred.var = "nppmmid30agg",
                        trim.outliers = TRUE,
                        prob = TRUE,
                        rug = TRUE)
pdp_npp <- pdp_npp %>% dplyr::mutate(Variable = 'NPP') %>% 
           dplyr::rename(value = 'nppmmid30agg')

pdp_cwd <- pdp::partial(GBM_k10_4_model,
                        pred.var = "deficitRS",
                        trim.outliers = TRUE,
                        prob = TRUE,
                        rug = TRUE)
pdp_cwd <- pdp_cwd %>% dplyr::mutate(Variable = 'CWD') %>% 
           dplyr::rename(value = 'deficitRS')

pdp_tpi <- pdp::partial(GBM_k10_4_model,
                        pred.var = "tpi30agg",
                        trim.outliers = TRUE,
                        prob = TRUE,
                        rug = TRUE)
pdp_tpi <- pdp_tpi %>% dplyr::mutate(Variable = 'TPI') %>% 
           dplyr::rename(value = 'tpi30agg')

pdp_all_together <- dplyr::bind_rows(pdp_uaa, pdp_cwd, pdp_npp, pdp_tpi) %>% 
                    dplyr::mutate(yhat = 1-yhat)

library(ggplot2)
pdp_all_together %>% 
  ggplot(aes(value, yhat)) + 
  geom_point(alpha = 0.25) +
  geom_line(aes(linetype = 'Partial Dependence'),
            key_glyph = 'timeseries') + 
  geom_rug(aes(x = NULL), size = 0.5) + 
  geom_smooth(aes(linetype = 'LOESS'), key_glyph = 'timeseries') +
  scale_linetype_manual(values = c(1,1), name = 'Model Fit') +
  guides(linetype = guide_legend(keywidth = 2.5,
                                 override.aes = list(linetype = c(1,1),
                                                     color = c('blue', 'black')))) + 
  facet_wrap(~Variable, scales = 'free') + 
  scale_x_continuous(labels = scales::comma_format()) +
  labs(x = 'Value', y = 'Probability of Occurence') +
  custom_theme()


# Spatial Maps ------------------------------------------------------------

## making the spatial maps

library(sf)
library(terra)
library(tidyverse)
proj.study <- "+proj=aea +lat_1=46 +lat_2=48 +lat_0=44 +lon_0=-109.5 +x_0=600000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

## read in the final predicted rasters
## 
topo_only <- rast('data/gbm_topo.tif')

topo_only_rc <- terra::classify(topo_only, matrix(c(0,0.5,1,
                                                 0.5,1,0), ncol = 3, byrow = T))

topoclim <- rast('data/gbm_topoclimatic.tif')

topoclim_rc <- terra::classify(topoclim, matrix(c(0,0.5,1,
                                               0.5,1,0), ncol = 3, byrow = T))

nhdStrOrdRast <- rast("data/nhdStrOrdRast.tif")

aoi1 <- mapedit::drawFeatures()
aoi2 <- mapedit::drawFeatures()
st_write(aoi1, 'hw_study.gpkg', layer = 'aoi1')
st_write(aoi2, 'hw_study.gpkg', layer = 'aoi2')
aoi <- aoi1 %>% st_transform(crs = st_crs(topo_only))

## topo only crop
topo_only_crop <- crop(topo_only_rc, vect(aoi))
output <- tempfile(fileext = '.tif')
terra::writeRaster(topo_only_crop, output, overwrite = T)

whitebox::wbt_line_thinning(output, paste0(getwd(),'/data/toc_thin.tif'))

toc_thin <- rast('data/toc_thin.tif')
toc_thin <- classify(toc_thin, matrix(c(0,NA),
                            ncol = 2,
                            byrow = T))

## topoclimatic crop
topoclim_only_crop <- crop(topoclim_rc, vect(aoi))
output <- tempfile(fileext = '.tif')
terra::writeRaster(topoclim_only_crop, output, overwrite = T)

whitebox::wbt_line_thinning(output, paste0(getwd(),'/data/tc_thin.tif'))

tc_thin <- rast('data/tc_thin.tif')
tc_thin <- classify(tc_thin, matrix(c(0,NA),
                                      ncol = 2,
                                      byrow = T))

library(stars)
nhdcrop <- crop(nhdStrOrdRast, vect(aoi)) %>% classify(matrix(c(1,4,1),
                                                              byrow = TRUE,
                                                              ncol = 3))%>% st_as_stars()

toc_stars <- st_as_stars(toc_thin) 
tc_stars <- st_as_stars(tc_thin )
nhd_stars <- st_as_stars(nhdcrop)

raster_ele <- elevatr::get_elev_raster(aoi,
                                       z = 13,
                                       clip = 'bbox')

raster_ele <- terra::terrain(rast(raster_ele), v = c('slope', 'aspect'), unit = 'radians')
raster_hill <- terra::shade(raster_ele$slope, raster_ele$aspect)

writeRaster(raster_hill, 'data/hillshade.tif', overwrite = T)
hillshade <- terra::rast("data/hillshade.tif")

hs_stars <- crop(hillshade, vect(aoi)) %>% st_as_stars()
data_sf <- read_sf("data/hw_study.gpkg", "points_spaced30m")
data_sf <- data_sf %>% dplyr::rename(stream = 'copy_TWI_1') %>% dplyr::mutate(stream = factor(stream))

data_sf_crop <- st_crop(data_sf, toc_thin)

## color pallete 'black' = topo-only, '#D81B60' = topoclimatic, '#1E88E5' = nhdplus

p1 <- ggplot() +     
  geom_stars(data = hs_stars, alpha = 0.4) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  geom_sf(data = data_sf_crop, aes(shape = stream, color = stream), size = 1, show.legend = F) + 
  scale_shape_manual(values = c(4,1),
                     labels = c("No Stream", "Stream"))+
  scale_color_manual(values =c('red', 'blue'),
                     labels = c("No Stream", "Stream")) +
    guides(shape = guide_legend(override.aes = list(alpha=1,color = c('red','blue')))) + 
  theme_void() + 
  theme(legend.position = 'bottom') + 
  labs(x = "Latitude", y = "Longitude")

p1

p2 <- ggplot() + 
  geom_stars(data = toc_stars, show.legend = FALSE)  + 
  scale_fill_gradientn(colors = 'black',na.value = "white")+ 
  ggnewscale::new_scale_fill()+ 
  geom_stars(data = hs_stars, alpha = 0.4) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  ggnewscale::new_scale_fill() +
  geom_sf(data = data_sf_crop, aes(shape = stream, color = stream),
          show.legend = F,
          size = 0, alpha = 0) +
  labs(title = 'Topo-only Model') +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

p2  


p3 <- ggplot() + 
  geom_stars(data = tc_stars, show.legend = FALSE)  + 
  scale_fill_gradientn(colors = '#D81B60',na.value = "white")+ 
  ggnewscale::new_scale_fill()+ 
  geom_stars(data = hs_stars, alpha = 0.4) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  ggnewscale::new_scale_fill() +
  geom_sf(data = data_sf_crop, aes(shape = stream, color = stream),
          show.legend = F,
          size = 0, alpha = 0) +
  labs(title = 'Topoclimatic Model') +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

p3  


p4 <- ggplot() + 
  geom_stars(data = nhd_stars, show.legend = FALSE)  + 
  scale_fill_gradientn(colors = '#1E88E5',na.value = "white")+ 
  ggnewscale::new_scale_fill()+ 
  geom_stars(data = hs_stars, alpha = 0.4) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  ggnewscale::new_scale_fill() +
  geom_sf(data = data_sf_crop,
          show.legend = F,
          size = 0, alpha = 0) +
  labs(title = 'NHDPlus') +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

p4
  
library(patchwork)
  
p2|p3|p4|p1


toc_gbm_stars <- st_as_stars(crop(topo_only,vect(aoi))) 
tc_gbm_stars <- st_as_stars(crop(topoclim, vect(aoi)))
p5 <- ggplot() + 
  geom_stars(data = 1-toc_gbm_stars, show.legend = F)  + 
  scale_fill_gradientn(colors = hcl.colors(11,'RdBu'),na.value = "white")+ 
  ggnewscale::new_scale_fill()+ 
  geom_stars(data = hs_stars, alpha = 0.4) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  ggnewscale::new_scale_fill() +
  geom_sf(data = data_sf_crop, 
          show.legend = F,
          size = 0, alpha = 0) +
  cowplot::theme_nothing() +
  theme(legend.position = 'bottom')

p5


p6 <- ggplot() + 
  geom_stars(data = 1-tc_gbm_stars, show.legend = F)  + 
  scale_fill_gradientn(colors = hcl.colors(11,'RdBu'),na.value = "white")+ 
  ggnewscale::new_scale_fill()+ 
  geom_stars(data = hs_stars, alpha = 0.4, show.legend = F) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  ggnewscale::new_scale_fill() +
  geom_sf(data = data_sf_crop,
          show.legend = F,
          size = 0, alpha = 0) +
  cowplot::theme_nothing() +
  theme(legend.position = 'bottom',
        legend.justification = "left",
        legend.key.width = unit(1.5, 'cm'))

p6

## topo only crop

aoi <- aoi2 %>% st_transform(crs = st_crs(topo_only))

topo_only_crop <- crop(topo_only_rc, vect(aoi))
output <- tempfile(fileext = '.tif')
terra::writeRaster(topo_only_crop, output, overwrite = T)

whitebox::wbt_line_thinning(output, paste0(getwd(),'/data/toc_thin.tif'))

toc_thin <- rast('data/toc_thin.tif')
toc_thin <- classify(toc_thin, matrix(c(0,NA),
                                      ncol = 2,
                                      byrow = T))

## topoclimatic crop
topoclim_only_crop <- crop(topoclim_rc, vect(aoi))
output <- tempfile(fileext = '.tif')
terra::writeRaster(topoclim_only_crop, output, overwrite = T)

whitebox::wbt_line_thinning(output, paste0(getwd(),'/data/tc_thin.tif'))

tc_thin <- rast('data/tc_thin.tif')
tc_thin <- classify(tc_thin, matrix(c(0,NA),
                                    ncol = 2,
                                    byrow = T))

library(stars)
nhdcrop <- crop(nhdStrOrdRast, vect(aoi)) %>% classify(matrix(c(1,4,1),
                                                              byrow = TRUE,
                                                              ncol = 3))%>% st_as_stars()

toc_stars <- st_as_stars(toc_thin) 
tc_stars <- st_as_stars(tc_thin )
nhd_stars <- st_as_stars(nhdcrop)

raster_ele <- elevatr::get_elev_raster(aoi,
                                       z = 13,
                                       clip = 'bbox')

raster_ele <- terra::terrain(rast(raster_ele), v = c('slope', 'aspect'), unit = 'radians')
raster_hill <- terra::shade(raster_ele$slope, raster_ele$aspect)

writeRaster(raster_hill, 'data/hillshade.tif', overwrite = T)
hillshade <- terra::rast("data/hillshade.tif")

hs_stars <- crop(hillshade, vect(aoi)) %>% st_as_stars()
data_sf <- read_sf("data/hw_study.gpkg", "points_spaced30m")
data_sf <- data_sf %>% dplyr::rename(stream = 'copy_TWI_1') %>% dplyr::mutate(stream = factor(stream))

data_sf_crop <- st_crop(data_sf, toc_thin)


p7 <- ggplot() +     
  geom_stars(data = hs_stars, alpha = 0.4) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  geom_sf(data = data_sf_crop, aes(shape = stream, color = stream), size = 1) + 
  scale_shape_manual(values = c(4,1),
                     labels = c("No Stream", "Stream"))+
  scale_color_manual(values =c('red', 'blue'),
                     labels = c("No Stream", "Stream")) +
  guides(shape = guide_legend(override.aes = list(alpha=1,color = c('red','blue')))) + 
  labs(color = '', shape = '') + 
  theme_void() + 
  theme(legend.position = 'bottom') + 
  labs(x = "Latitude", y = "Longitude")

p7

p8 <- ggplot() + 
  geom_stars(data = toc_stars, show.legend = FALSE)  + 
  scale_fill_gradientn(colors = 'black',na.value = "white")+ 
  ggnewscale::new_scale_fill()+ 
  geom_stars(data = hs_stars, alpha = 0.4) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  ggnewscale::new_scale_fill() +
  geom_sf(data = data_sf_crop, aes(shape = stream, color = stream),
          alpha = 0,
          show.legend = F,
          size = 0) + 
  scale_shape_manual(values = c(4,1),
                     labels = c("No Stream", "Stream"))+
  scale_color_manual(values =c('red', 'blue'),
                     labels = c("No Stream", "Stream")) +
  #labs(title = 'Topo-only Model') +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

p8 


p9 <- ggplot() + 
  geom_stars(data = tc_stars, show.legend = FALSE)  + 
  scale_fill_gradientn(colors = '#D81B60',na.value = "white")+ 
  ggnewscale::new_scale_fill()+ 
  geom_stars(data = hs_stars, alpha = 0.4) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  ggnewscale::new_scale_fill() +
  geom_sf(data = data_sf_crop, aes(shape = stream, color = stream),
          alpha = 0,
          show.legend = F,
          size = 1) + 
  scale_shape_manual(values = c(4,1),
                     labels = c("No Stream", "Stream"))+
  scale_color_manual(values =c('red', 'blue'),
                     labels = c("No Stream", "Stream")) +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

p9  


p10 <- ggplot() + 
  geom_stars(data = nhd_stars, show.legend = FALSE)  + 
  scale_fill_gradientn(colors = '#1E88E5',na.value = "white")+ 
  ggnewscale::new_scale_fill()+ 
  geom_stars(data = hs_stars, alpha = 0.4) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  ggnewscale::new_scale_fill() +
  geom_sf(data = data_sf_crop, aes(shape = stream, color = stream),
          alpha = 0,
          show.legend = F,
          size = 1) + 
  scale_shape_manual(values = c(4,1),
                     labels = c("No Stream", "Stream"))+
  scale_color_manual(values =c('red', 'blue'),
                     labels = c("No Stream", "Stream")) +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

p10

library(patchwork)



toc_gbm_stars <- st_as_stars(crop(topo_only,vect(aoi))) 
tc_gbm_stars <- st_as_stars(crop(topoclim, vect(aoi)))
p11 <- ggplot() + 
  geom_stars(data = 1-toc_gbm_stars, show.legend = T)  + 
  scale_fill_gradientn(colors = hcl.colors(11,'RdBu'),
                       name = 'Probability',
                       na.value = "white")+ 
  ggnewscale::new_scale_fill()+ 
  geom_stars(data = hs_stars, alpha = 0.4) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  ggnewscale::new_scale_fill() +
  geom_sf(data = data_sf_crop, 
          show.legend = F,
          size = 0, alpha = 0) +
  cowplot::theme_nothing() +
  theme(legend.position = 'bottom',
        legend.justification = "left",
        legend.key.width = unit(1.25, 'cm'))

p11


p12 <- ggplot() + 
  geom_stars(data = 1-tc_gbm_stars, show.legend = F)  + 
  scale_fill_gradientn(colors = hcl.colors(11,'RdBu'),na.value = "white")+ 
  ggnewscale::new_scale_fill()+ 
  geom_stars(data = hs_stars, alpha = 0.4, show.legend = F) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  ggnewscale::new_scale_fill() +
  geom_sf(data = data_sf_crop,
          show.legend = F,
          size = 0, alpha = 0) +
  cowplot::theme_nothing() +
  theme(legend.position = 'bottom',
        legend.justification = "left",
        legend.key.width = unit(1.5, 'cm'))

p12

(p2|p3|p4)/(p5|p6|p1)/(p8|p9|p10)/(p11|p12|p7)


aoi2 <- mapedit::drawFeatures()

aoi2 <- aoi2 %>% st_transform(crs = st_crs(topo_only))