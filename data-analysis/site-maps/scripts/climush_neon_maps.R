library(fs)
library(grid)

## FILE PATHS ##########################################################################################################

## WORKING DIRECTORY ###################################################

working_dir <- path('/Users/carolyndelevich/main/github_repos/climush/data-analysis/site-maps/')
setwd(working_dir)


## INPUT ###############################################################

## DIRECTORIES ########

# mycopull project repository, where the mapping settings / functions are located
mycopull_mapping_scripts <- path('/Users/carolyndelevich/main/projects/hja-macrofungi/hja-macrofungi/data-analysis/mapping/scripts')

# main directory for site map data
sitemaps_data_main <- path(working_dir, 'data')

# main directory for site map figures
sitemaps_figs_main <- path(working_dir, 'figures')

# main reference data directory for NEON ecoregion mapping data
neon_sites_main <- path(sitemaps_data_main, 'neon-sites')

# reference data subdirectory with domain area shape files
neon_domains_dir <- path(neon_sites_main, 'domain-polygons')

# reference data subdirectory with NEON site boundary shape files
neon_sitebounds_dir <- path(neon_sites_main, 'site-boundaries')


# main reference data directory for LTER ecoregion mapping data
lter_sites_main <- path(sitemaps_data_main, 'lter-sites')

# directory with shapefiles of LTER site boundaries
lter_sitebounds_dir <- path(lter_sites_main, 'site-boundaries')


## FILES ##############

# mycopull mapping source file with mapping settings and functions
mycopull_srcfile_pathin <- path(mycopull_mapping_scripts, 'mapping-source.R')

# shapefile of the NEON ecoregion domain area polygons
neon_domains_pathin <- path(neon_domains_dir, 'NEON_Domains.shp')

# shapefile of the NEON site boundaries
neon_sitebounds_pathin <- path(neon_sitebounds_dir, 'terrestrialSamplingBoundaries.shp')

# shapefile of the LTER site boundaries
lter_sitebounds_pathin <- path(lter_sitebounds_dir, 'lter_site-boundaries.shp')


## OUTPUT ##############################################################

## DIRECTORIES ########

## FILES ##############

# filename for the output map of the climush NEON sites (legend on bottom)
map_pathout <- path(sitemaps_figs_main, 'climush_neon_maps.png')

# filename for  the output map of the climush NEON sites (legend on top)
map_toplegend_pathout <- path(sitemaps_figs_main, 'climush_neon_maps_legendtop.png')

########################################################################################################################


## SOURCE MAPPING SETTINGS #############################################################################################

# source the mycopull mapping file for settings and functions for mapping
source(mycopull_srcfile_pathin)

# reset working directory to override working directory set by src file
setwd(working_dir)

########################################################################################################################


## IMPORT SHAPEFILES ###################################################################################################

## NEON SITES ##########################################################

# import the shapefile of the NEON ecoregion domains
neon_domains_sf <- read_sf(neon_domains_pathin)

# import the shapefile of the NEON site boundaries
neon_sitebounds_sf <- read_sf(neon_sitebounds_pathin)


## LTER SITES ##########################################################

# import the shapefile of the LTER site boundaries
lter_sitebounds_sf <- read_sf(lter_sitebounds_pathin)

########################################################################################################################


## SUBSET CLIMUSH / NON-CLIMUSH NEON DOMAINS ###########################################################################

# list of domain numbers that contain a climush site, used to filter the domain polygons
climush_ecoregions_domainNumb <- c('D19', 'D14', 'D13', 'D03', 'D06', 'D01', 'D05', 'D16')

# ecoregions with a climush site (NEON or LTER)
domains_climush_sf <- neon_domains_sf %>%
  filter(domainID %in% climush_ecoregions_domainNumb) %>%
  mutate(domainName = case_when(
    domainName == 'Southern Rockies / Colorado Plateau' ~ 'Southern Rockies + Colorado Plateau',
    .default = domainName
  ))

# ecoregions without a climush site
domains_nonclimush_sf <- neon_domains_sf %>%
  filter(!domainID %in% climush_ecoregions_domainNumb)

########################################################################################################################


## FILTER + REFORMAT SITE BOUNDARIES ###################################################################################

## COLUMN RENAME MAPPING ###############################################

# the new column names to use for the combine site boundaries shapefile object
sitebounds_colnames <- c(
  'domain_code',
  'domain_name',
  'site_code',
  'site_name'
)

# create a renaming map of the old NEON site boundary column names with the new site boundary column names
neon_sitebounds_rename <- c('domainNumb', 'domainName', 'siteID', 'siteName')
names(neon_sitebounds_rename) <- sitebounds_colnames

# create a renaming map of the old LTER site boundary column names with the new site boundary column names
lter_sitebounds_rename <- c('domain_code', 'domain_name', 'SITE', 'NAME')
names(lter_sitebounds_rename) <- sitebounds_colnames


## FILTER + REFORMAT SITE BOUNDARIES ###################################

## NEON ###############

# list of site names for climush sites that are a NEON site, used to filter the site centroids
climush_neon_sites <- c(
  'Harvard Forest',                         # also an LTER site
  'Ordway Swisher Biological Station',
  'Konza Prairie Biological Station',       # also an LTER site
  'Niwot Ridge Mountain Research Station',  # also an LTER site
  'Santa Rita Experimental Range'
)

# filter the NEON site boundary shapefile to include only climush sites; drop extra columns + rename to match LTER
climush_neon_sitebounds_sf <- neon_sitebounds_sf %>%
  filter(siteName %in% climush_neon_sites) %>%             # filter by site name
  filter(siteType == 'Core Terrestrial') %>%               # Konza has both Core + Gradient Terrestrial; only keep Core
  select(c(domainNumb, domainName, siteID, siteName)) %>%  # drop unused columns
  rename(neon_sitebounds_rename)                           # rename remaining columns to match LTER shapefile


## LTER ###############

# list of site names for climush sites that are a LTER site, used to filter the site centroids
climush_lter_sites <- c(
  'Andrews',         # technically also a NEON site, but only watersheds (McRae Creek)
  'Bonanza Creek',
  'Cedar Creek',
  'Harvard Forest',  # also a NEON site
  'Konza Prairie',   # also a NEON site
  'Niwot Ridge'      # also a NEON site
)

# assign NEON domain codes for each of the LTER sites by name
climush_lter_domaincodes <- c('D16', 'D19', 'D05', 'D01', 'D06', 'D13')

# match the LTER site names to their corresponding NEON domain codes
names(climush_lter_domaincodes) <- climush_lter_sites

# match the NEON domain codes with their domain name
climush_neon_domains <- domains_climush_sf %>%
  filter(domainID %in% climush_lter_domaincodes) %>%
  arrange(domainID) %>%
  .$domainID
names(climush_neon_domains) <- domains_climush_sf %>%
  filter(domainID %in% climush_lter_domaincodes) %>%
  arrange(domainID) %>%
  .$domainName

# filter the LTER site boundary shapefile to only include climush sites; rename columns to match NEON
climush_lter_sitebounds_sf <- lter_sitebounds_sf %>%
  filter(NAME %in% climush_lter_sites) %>%
  arrange(NAME) %>%
  mutate(domain_code = climush_lter_domaincodes) %>%
  arrange(domain_code) %>%
  mutate(domain_name = names(climush_neon_domains)) %>%  # already sorted by ascending domain code
  rename(lter_sitebounds_rename) %>%
  relocate(geometry, .after = last_col()) %>%
  relocate(c(site_code, site_name), .before = geometry) %>%
  mutate(site_name = case_when(
    site_name == 'Niwot Ridge' ~ 'Niwot Ridge Mountain Research Station',
    site_name == 'Konza Prairie' ~ 'Konza Prairie Biological Station',
    site_name == 'Andrews' ~ 'Andrews Forest',
    site_name == 'Cedar Creek' ~ 'Cedar Creek Ecosystem Science Reserve',
    .default = site_name)) %>%
  mutate(domain_name = case_when(
    domain_name == 'Southern Rockies / Colorado Plateau' ~ 'Southern Rockies & Colorado Plateau',
    .default = domain_name))

########################################################################################################################


## CONVERT SITE BOUNDARIES TO CENTROIDS ################################################################################

# must set to FALSE to avoid the error - Edge x has duplicate vertex with edge y
# sf_use_s2(FALSE)

## NEON ###############

# replace the multipolygon geometries with the centroids of each polygon
climush_neon_sitecentroids_sf <- climush_neon_sitebounds_sf %>%
  mutate(geometry = st_centroid(geometry))

## LTER ###############

# replace the multipolygon geometries with the centroids of each polygon
climush_lter_sitecentroids_sf <- climush_lter_sitebounds_sf %>%
  mutate(geometry = st_centroid(geometry))

# reset to original setting after calculating centroids
# sf_use_s2(TRUE)

########################################################################################################################


## ADD CENTROID FOR NON-LTER / NON-NEON SITE ###########################################################################

# x / y coordinates of the Mt. Pisgah grasslands site in OR (Pacific Northwest), which is neither NEON nor LTER
#   coords based on Google Maps pindrop approx on site
climush_misc_sitecentroids_sf <- st_sf(
  domain_code = 'D16',
  domain_name = 'Pacific Northwest',
  site_code = 'MTPS',
  site_name = 'Mount Pisgah Arboretum',
  geometry = st_sfc(st_point(c(-122.9396, 44.000)))  # longitude then latitude
)

########################################################################################################################


## COMBINE SITE COORDINATES ############################################################################################

# combine the site coordinates of the climush sites into one shapefile object
climush_sitecentroids_sf <- st_as_sf(
  rbind(
    as.data.frame(climush_neon_sitecentroids_sf),
    as.data.frame(climush_lter_sitecentroids_sf),
    as.data.frame(climush_misc_sitecentroids_sf)
  )
)

########################################################################################################################


## GET PLOT COORDINATE X / Y LIMS ######################################################################################

# create rough boundaries for continental US, for clipping main NEON shapefile
l48_lims <- c(-70, 25, -130, 50)
names(l48_lims) <- c('xmin', 'ymin', 'xmax', 'ymax')
l48_bbox <- st_bbox(l48_lims)

# crop the input domains shapefile to include only continental US; get x/y axes limits from this version of the shapefile
l48_axes_lims <- domains_climush_sf %>%
  st_crop(l48_bbox) %>%
  get_mapping_axes_limits(buffer_prop = 0.05, shrink_or_expand = 'expand')

# create rough boundaries for AK, for clipping main NEON shapefile
ak_lims <- c(st_bbox(neon_domains_sf)$xmin, 50, -120, st_bbox(neon_domains_sf)$ymax+10)
names(ak_lims) <- c('xmin', 'ymin', 'xmax', 'ymax')
ak_bbox <- st_bbox(ak_lims)

# crop the input domains shapefile to include only AK; get x/y axes limits for AK map from this version of the shapefile
AK_axes_lims <- neon_domains_sf %>%
  st_crop(ak_bbox) %>%
  get_mapping_axes_limits(buffer_prop = 0.05, shrink_or_expand = 'expand')

########################################################################################################################


## CREATE MAPPING LAYERS ###############################################################################################

# mapping layer of climush NEON ecoregions
domains_climush_layer <- geom_sf(
  data = domains_climush_sf,
  aes(fill = domainName),
  color = 'black',
  size = 0.25,
  alpha = 1,
  show.legend = TRUE)  # must set to TRUE in order to extract legend (below)

domains_other_layer <- geom_sf(
  data = domains_nonclimush_sf,
  fill = 'grey',
  color = 'dark grey',
  size = 0.25,
  alpha = 0.25,
  show.legend = FALSE)  # set to FALSE; legend will only show climush ecoregions

# mapping layer of centroids of sites in the NEON ecoregions of the lower 48
sitecentroids_layer <- geom_sf(
  data = climush_sitecentroids_sf,
  aes(fill = domain_name),
  color = 'black',
  pch = 21,
  size = 4,
  alpha = 1,
  show.legend = FALSE)

# extract the legend from the main climush NEON ecoregion figure
domains_src <- ggplot() +
  domains_climush_layer +
  scale_fill_manual(values = climush_ecoregions_colors) +
  guides(
    fill = guide_legend(
      ncol = 2,
      position = 'bottom',
      theme = theme(
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(0.25, 'in'),
        legend.background = element_blank(),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'in')
      )
    )
  )
domains_legend <- get_legend(domains_src)
if (is.null(domains_legend))
  stop('Legend was not successfully extracted from figure. Check that show.legend = TRUE when first creating ggplot layer.')


# create map of the climush ecoregions and climush NEON site points in lower 48 states
domains_continental_map <- ggplot() +
  domains_other_layer +
  domains_climush_layer +
  sitecentroids_layer +
  scale_fill_manual(values = climush_ecoregions_colors) +
  scale_y_continuous(limits = c(l48_axes_lims$ymin, l48_axes_lims$ymax)) +
  scale_x_continuous(limits = c(l48_axes_lims$xmin, l48_axes_lims$xmax)) +
  theme_void() +
  theme(
    legend.position = 'none',
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'in')
  )

# create map of the AK NEON climush sites + size points
domains_AK_map <- ggplot() +
  domains_other_layer +
  domains_climush_layer +
  sitecentroids_layer +
  scale_fill_manual(values = climush_ecoregions_colors) +
  scale_y_continuous(limits = c(AK_axes_lims$ymin, AK_axes_lims$ymax)) +
  scale_x_continuous(limits = c(AK_axes_lims$xmin, AK_axes_lims$xmax)) +
  theme_void() +
  theme(
    legend.position = 'none',
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'in')
  )

rect_fill <- rectGrob(
  width = unit(2.5, "in"),
  height = unit(2.5, "in"),
  gp = gpar(fill = "cornsilk", alpha = 0.25)
)
rect_outline <- rectGrob(
  width = unit(2.5, "in"),
  height = unit(2.5, "in"),
  gp = gpar(color = "black", alpha = 0.25, size = 0.25)
)

domains_l48plusAK_map <- ggdraw(domains_continental_map) +
  draw_grob(rect_outline,
            x = -0.325,
            y = -0.3) +
  draw_plot(domains_AK_map,
            x = -0.325,
            y = -0.3,
            # hjust = 0,
            # vjust = 0,
            scale = 0.25)


# add the legend to the bottom of the map of the continental US + AK
domains_map_wlegend <- ggdraw(domains_l48plusAK_map) +
  draw_plot(domains_legend,
            x = 0.15,
            y = -0.36,
            hjust = 0,
            vjust = 0
  )


ggsave(
  plot = domains_map_wlegend,
  filename = map_pathout,
  height = 10,
  width = 10,
  units = 'in',
  dpi = 600
)


rect_outline <- rectGrob(
  width = unit(2, "in"),
  height = unit(2, "in"),
  gp = gpar(color = "black", 
            alpha = 0.5, 
            lwd = 3)
)

# add AK to the top of the continental US map
domains_l48plusAK_maptop <- ggdraw(domains_continental_map) +
  draw_grob(rect_outline,
            x = -0.3,
            y = 0.26) +
  draw_plot(domains_AK_map,
            x = -0.3,
            y = 0.26,
            hjust = 0,
            vjust = 0,
            scale = 0.22)

# add legend to the top of the map of the continental US + AK
domains_map_wtoplegend <- ggdraw(domains_l48plusAK_maptop) +
  draw_plot(domains_legend,
            x = 0.15,
            y = 0.2275,
            hjust = 0,
            vjust = 0)


save_plot(filename = map_toplegend_pathout,
          plot = domains_map_wtoplegend,
          ncol = 1,
          nrow = 2,
          base_asp = 1.618,
          base_height = NULL,  # leave null; will calculate based on aspect ratio and base_width
          base_width = 10,
          unit = 'in',
          dpi = 600)

########################################################################################################################