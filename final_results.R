library(ggtext)
library(tidyverse)
library(ggridges)
library(wesanderson)
source('utils.R')

# for cairo_view()

trace(grDevices::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)


# Figure 3 ----------------------------------------------------------------

all_results <- read_csv("data/all_results_final.csv")

ffs_best <- all_results %>% group_by(model,cluster) %>% slice_max(order_by = Accuracy)%>% ungroup()

all_models <- ffs_best %>% mutate(cluster_no = str_replace_all(cluster, '[.]', '_'),
                                  cluster_no = parse_number(cluster_no))

all_models_gbm <- all_models %>% filter(model == 'gbm')
all_model_vars <- all_models_gbm %>% separate(char_vars, into = letters[1:15], sep = ',') %>% 
  pivot_longer(cols = 44:54) %>% mutate(value = str_remove_all(value, ' ') %>% factor(),
                                        value = fct_recode(value, NPP = 'nppmmid30agg',
                                                           UAA = 'accum30', TPI = 'tpi30agg',
                                                           CWD = 'deficitRS', Deciduous = 'decid30RS',
                                                           CAD = 'cad30RS', TWI = 'twi30agg',
                                                           NDVI = 'NDVI_med', NPOL = 'npol30agg',
                                                           NIR = 'NIR_med', Green = 'Green_med',
                                                           `VV sd` = 'vvsd30agg', VV = 'vv30agg',
                                                           `CPG Precip` = 'cpg30precip'))


#figure in paper
all_model_vars %>% count(value,model, sort = T) %>% na.omit() %>% mutate(value = fct_reorder(value, n)) %>% 
  ggplot(aes(value, n)) + 
  geom_col(alpha = 0.5) + 
  geom_point()  +
  labs(x = "Features", y = "Count") +
  coord_flip() + 
  geom_text(aes(label = n), nudge_y = 3) + 
  custom_theme() + 
  theme(axis.title = element_text(size = 14), axis.text =element_text(size = 12))


#### generating table 2
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


# Figure 2  ---------------------------------------------------------------

## figure 2 in paper  


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
                        labels = c('Topography-only','Topoclimatic','NHDPlus HR')) +
  scale_color_manual(values = c('black', '#D81B60', '#1E88E5'),
                     labels = c('Topography-only','Topoclimatic','NHDPlus HR')) + 
guides(
  linetype = guide_legend(keywidth = 2.5),
  color = guide_legend(keywidth = 2.5,
                      override.aes = list(
                        linetype = c(1, 2, 3),
                        size = c(.5, .5, .5),
                      reverse = T))) + 
  labs(color = 'Model Type', linetype = 'Model Type', x = '') +
  facet_wrap(~name, scales = 'free')  +
  custom_theme() +
  theme(strip.background = element_rect(fill = 'white',
                                        color = 'black',
                                        size = 1))

ggsave('images/test.jpg',device='jpeg', dpi = 300, width = 30, height = 70, units = 'mm')
# Figures 4 and 5 ---------------------------------------------------------

### figure 4 and 5
library(pdp)
library(vip)
library(gbm)
# this is a tuned GBM model being used with 4 variables and 10 cluster hold out

GBM_k10_4_model <- readRDS("data/GBM_k10_4_model.rds")

inter <- pdp::partial(GBM_k10_4_model$finalModel, 
                      pred.var = c("accum30","deficitRS"),
                      type = 'classification',
                      n.trees = 1500,
                      train = as.data.frame(GBM_k10_4_model$trainingData),
                      trim.outliers = TRUE, prob = TRUE)
def_npp <- pdp::partial(GBM_k10_4_model$finalModel, pred.var = c("deficitRS","nppmmid30agg"),
                        type = 'classification',
                        n.trees = 1500,
                        train = as.data.frame(GBM_k10_4_model$trainingData),
                        trim.outliers = TRUE, prob = TRUE) 

inter <- inter %>% dplyr::mutate(yhat = 1-yhat)
pdp_uaa_def <- plotPartial(inter,
                    levelplot = F,
                    plot.pdp = F,
                    rug = TRUE,
                    zlab = "yhat",
                    drape=T,
                    smooth = TRUE,
                    colorkey = T,
                    contour = T, 
                    #main = "Topoclimatic GBM Interactions (CWD and UAA)",
                    screen = list(z = 120, x = -50, y = 0),
                    scales = list(arrows = F),
                    xlab = 'UAA',
                    col.regions = rev(hcl.colors(100, 'Grays')),
                    ylab = 'CWD',
                    pretty = T, 
                    chull = TRUE)
pdp_uaa_def

def_npp <- def_npp %>% dplyr::mutate(yhat = 1-yhat)

pdp_def_npp <- plotPartial(def_npp,
                    levelplot = FALSE,
                    zlab = "yhat",
                    drape=T,
                    rug = TRUE,
                    smooth = TRUE,
                    colorkey = T,
                    contour = TRUE, 
                    #main = "Topoclimatic GBM Interactions (CWD and NPP)",
                    screen = list(z = -10, x = -75, y = -30),
                    col.regions = rev(hcl.colors(100, 'Grays')),
                    scales = list(arrows = F),
                    arrows = F,
                    pretty = T, 
                    xlab = 'CWD',
                    ylab = 'NPP')


pdp_def_npp

library(cowplot)
plot_grid(pdp_uaa_def, pdp_def_npp, 
          label_size = 12,
          labels = c('A)', 'B)'),
          hjust = -.25, vjust = 10.5)
#now make the single pdp plots 

pdp_uaa <- pdp::partial(GBM_k10_4_model$finalModel,
                        pred.var = "accum30",
                        trim.outliers = TRUE,
                        train = as.data.frame(GBM_k10_4_model$trainingData),
                        n.trees = 1500,
                        prob = TRUE,
                        rug = TRUE)
pdp_uaa <- pdp_uaa %>% dplyr::mutate(Variable = 'UAA') %>% 
           dplyr::rename(value = 'accum30')

pdp_npp <- pdp::partial(GBM_k10_4_model$finalModel,
                        pred.var = "nppmmid30agg",
                        trim.outliers = TRUE,
                        train = as.data.frame(GBM_k10_4_model$trainingData),
                        n.trees = 1500,
                        prob = TRUE,
                        rug = TRUE)
pdp_npp <- pdp_npp %>% dplyr::mutate(Variable = 'NPP') %>% 
           dplyr::rename(value = 'nppmmid30agg')

pdp_cwd <- pdp::partial(GBM_k10_4_model$finalModel,
                        pred.var = "deficitRS",
                        trim.outliers = TRUE,
                        train = as.data.frame(GBM_k10_4_model$trainingData),
                        n.trees = 1500,
                        prob = TRUE,
                        rug = TRUE)
pdp_cwd <- pdp_cwd %>% dplyr::mutate(Variable = 'CWD') %>% 
           dplyr::rename(value = 'deficitRS')

pdp_tpi <- pdp::partial(GBM_k10_4_model$finalModel,
                        pred.var = "tpi30agg",
                        trim.outliers = TRUE,
                        train = as.data.frame(GBM_k10_4_model$trainingData),
                        n.trees = 1500,
                        prob = TRUE,
                        rug = TRUE)
pdp_tpi <- pdp_tpi %>% dplyr::mutate(Variable = 'TPI') %>% 
           dplyr::rename(value = 'tpi30agg')

pdp_all_together <- dplyr::bind_rows(pdp_uaa, pdp_cwd, pdp_npp, pdp_tpi) %>% 
                    dplyr::mutate(yhat = 1-yhat)

library(ggplot2)
pdp_all_together %>% 
  ggplot(aes(value, yhat)) +
  geom_point(alpha = 0.15) + 
  geom_line(alpha = 0.15, lwd = 1) +
  # geom_line(aes(linetype = 'Partial Dependence'),
  #           key_glyph = 'timeseries') + 
  geom_rug(aes(x = NULL), size = 0.5) + 
  geom_smooth(aes(linetype = 'LOESS'), key_glyph = 'timeseries') +
  scale_linetype_manual(values = c(1), name = 'Model Fit') +
  guides(linetype = guide_legend(keywidth = 2.5,
                                 override.aes = list(linetype = c(1),
                                                     color = c('blue')))) + 
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
aoi2 <- mapedit::editFeatures(aoi1)
st_write(aoi1, 'hw_study.gpkg', layer = 'aoi1')
st_write(aoi2, 'hw_study.gpkg', layer = 'aoi2', delete_layer = T)
aoi1 <- read_sf('hw_study.gpkg', layer = 'aoi1')
aoi2 <- read_sf('hw_study.gpkg', layer = 'aoi2')

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
data_sf <- data_sf %>% dplyr::rename(stream = 'copy_TWI_1') %>% 
           dplyr::mutate(stream = factor(stream)) %>% st_transform(st_crs(aoi))

data_sf_crop <- st_crop(data_sf, toc_thin)

## color pallete 'black' = topo-only, '#D81B60' = topoclimatic, '#1E88E5' = nhdplus

p1 <- ggplot() +     
  geom_stars(data = hs_stars, alpha = 0.4) + 
  scale_fill_distiller(palette = "Greys", guide = 'none', na.value = NA) + 
  geom_sf(data = data_sf_crop, aes(color = stream), size = 1, show.legend = F) + 
  scale_shape_manual(
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
  labs(title = 'NHDPlus HR') +
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
  geom_sf(data = data_sf_crop, aes(color = stream), size = 1) + 
  scale_shape_manual(
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


plot_prediction <- (p2|p3|p4)/(p5|p6|p1)/(p8|p9|p10)/(p11|p12|p7)
plot_prediction

top_plot <- plot_grid(p2, p3, p4,
                      p5, p6, p1,
                      p8, p9, p10,
                      p11, p12, p7, ncol = 3, rel_widths = c(5,5,5,1,1,1,1,1,1,2,2,2),
                      labels = 'AUTO')
top_plot

bottom_plot <- plot_grid(p8, p9, p10,
                         p11, p12, p7,
                         ncol = 3,
                         labels = 'B)')

plot_grid(top_plot, bottom_plot, ncol = 1, rel_widths = c(5,5,5,1,1,1,1,1,1,2,2,2))

aoi2 <- mapedit::drawFeatures()

aoi2 <- aoi2 %>% st_transform(crs = st_crs(topo_only))


#### figure 1 and parts of figure 6
# get the 
data34 <- read_csv('data34.csv') %>% st_as_sf(coords = c('utm1', 'utm2'), crs = proj.study)
def <- rast('deficitRS.tif')
def[def> 500] <- NA
def <- def %>% st_as_stars()


district <- read_sf('data/district_boundary.shp') %>% st_transform(4326)



main.plot <- ggplot() + 
  geom_stars(data = def) + 
  scale_fill_gradientn(colors = rev(hcl.colors(11,'RdYlGn')),na.value = NA)+ 
  #ggnewscale::new_scale_fill() +
  geom_sf(data = data34, aes(color = stream),
          shape = 21, size = 0.5) + 
  scale_color_manual(values = c('red', 'blue'), labels = c('No HWS', 'HWS')) +
  guides(color = guide_legend(override.aes = list(size = c(2, 2)))) + 
  labs(color = 'HWS Occurrence', 
       fill = 'Climatic Water Deficit (mm)') +
  labs(y = 'Latitude', x = 'Longitude', shape = '')+ 
  geom_sf(data = district , fill = NA,color = 'black',aes(shape = 'Modeling Extent'), lwd = 1) +
  theme_bw()+
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey60", "white"),
    text_family = "ArcherPro Book"
  ) +
  ggspatial::annotation_north_arrow(
    location = "br", which_north = "true",
    pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20",
      text_family = "ArcherPro Book"
    )
  )
main.plot

bbox <- stats::setNames(sf::st_bbox(sf::st_buffer(district, 0.05)), 
                        c("left", "bottom", "right", "top"))
basemap_satellite <- ggmap::get_map(maptype = 'satellite', location = bbox, 
                                    zoom = 7, source = 'google')

# Define a function to fix the bbox to be in EPSG:3857
ggmap_bbox <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")), 
                       c("ymin", "xmin", "ymax", "xmax"))
  
  # Coonvert the bbox to an sf polygon, transform it to 3857, 
  # and convert back to a bbox (convoluted, but it works)
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
  map
}
class(mt_area)
# Use the function:
map <- ggmap_bbox(basemap_satellite)
maoi <- mapedit::drawFeatures()
us <- AOI::aoi_get(state = 'conus') %>% st_transform(3857) %>% filter(name %in% c('Montana', 'Idaho', 'Washington')) %>% 
  st_crop(st_transform(maoi, st_crs(.))) 

def_proj <- stars::st_transform_proj(def, 3857)

mt_area <-  
  # ggmap::ggmap(map) +
  # coord_sf(crs = st_crs(3857)) +
  ggplot() + 
  geom_stars(data = def_proj, inherit.aes = F, show.legend = F, downsample = 4) + 
  scale_fill_gradientn(colors = rev(hcl.colors(11,'RdYlGn')),na.value = NA)+ 
  geom_sf(data = us, fill = NA, color = 'black', 
          lwd = .5, 
          inherit.aes = F) + 
           geom_sf(data = district %>% st_transform(3857), fill = NA, color = 'black', 
                   lwd = 0.5, 
                   inherit.aes = F) + theme_void()
mt_area

ggdraw(main.plot) +
  draw_plot(mt_area, .09, .07, .2, .2) +
  draw_label("Montana", .225, .16, .2, .2, color = "black", size = 10, angle = 0)+
  draw_label("Idaho", .15, .1, .2, .2, color = "black", size = 10, angle = 0)

fig_6 <- ggplot() + 
  geom_stars(data = def_proj, downsample = 3) + 
  scale_fill_gradientn(colors = rev(hcl.colors(11,'RdYlGn')),na.value = NA)+ 
  #ggnewscale::new_scale_fill() +
  geom_sf(data = aoi1 %>% st_transform(3857), fill = NA, color = 'black', lwd = 0.75)  +
  geom_text(
    data = aoi1,
    aes(label = 'A', geometry = geom),
    stat = "sf_coordinates"
  )+ 
  geom_sf(data = aoi2%>% st_transform(3857), fill = NA, color = 'black', lwd = 0.75) +
  geom_text(
    data = aoi2,
    aes(label = 'B', geometry = geom),
    stat = "sf_coordinates"
  ) + 
  labs(
    fill = 'Climatic Water Deficit (mm)') +
  labs(y = 'Latitude', x = 'Longitude', shape = '')+ 
  geom_sf(data = district , fill = NA,color = 'black', lwd = 1) +
  theme_bw()+
  theme( legend.position = 'bottom') + 
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("grey60", "white"),
    text_family = "ArcherPro Book"
  ) +
  ggspatial::annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20",
      text_family = "ArcherPro Book"
    )
  ) + theme_void() + theme(legend.position = 'bottom')
plot_grid(plot_prediction, fig_6, ncol = 2, rel_widths = c(2, 1))
((plot_prediction)|fig_6) + 
  plot_layout(widths = c(2, 1), heights = c(2,1))
