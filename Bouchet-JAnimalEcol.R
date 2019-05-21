#' -----------------------------------------------------------
#' -----------------------------------------------------------
#
# Bouchet et al. Spatial dimensions of pelagic diversity in a geodiverse offshore marine park
# NESP Oceanic Shoals pelagic videography. Journal of Animal Ecology.
# 
#' -----------------------------------------------------------
#' -----------------------------------------------------------
#' 
#' Author: PJ Bouchet
#' Last update: May 16, 2019
#' 
#' -----------------------------------------------------------

#' ====================================
# LIBRARIES & SETTINGS ====
#' ====================================

# devtools::install_github("beckyfisher/FSSgam_package")

pack<-c("tidyverse", # Tidyverse for data science
        "reshape2", # Transform data between wide and long formats
        # "ecole", # For zero-adjusted Bray Curtis
        # "MASS", # For simulating from multivariate Normal
        "raster", # GIS
        "FSSgam", # A simple function for full subsets multiple regression in ecology with R
        # "rgdal", # GIS
        # "bioDist",
        # "mgcv", # GLM/GAM models
        # "pals",
        # "gstat",
        # "dismo",
        # "rgeos",
        # "statmod",
        # "AER",
        # "car",
        # "pscl",
        "iNEXT", # Species rarefaction
        # "robustbase",
        # "scales",
        "vegan", # Community ecology: ordination, disversity & dissimilarities
        # "fANCOVA", # LOESS smoothing
        "rfishbase")

#'---------------------------------------------
# Loops through the list and loads (and installs if required) all packages
#'---------------------------------------------

for (i in 1:length(pack)){
  p<-pack[i]
  if(as.character(p) %in% rownames(installed.packages())==TRUE){
    require(p, character.only = TRUE)
  }else{install.packages(p)
    require(p, character.only = TRUE)}
}

#'---------------------------------------------
# Set tibble options
#'---------------------------------------------

options(tibble.width = Inf) # All tibble columns shown
options(pillar.neg = FALSE) # No colouring negative numbers
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)

set.seed(32) # Reproducible results

#' ====================================
# FUNCTIONS ====
#' ====================================

source("R/Functions.R")

#' ====================================
# GIS ====
#' ====================================

utmOS <- sp::CRS("+proj=utm +zone=52 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
wgs84 <- sp::CRS("+proj=longlat +datum=WGS84")

#'---------------------------------------------
# AMP boundary
#'---------------------------------------------

amp_boundaries <- raster::shapefile("gis/OceanicShoals_AMP.shp")
amp_boundaries <- sp::spTransform(amp_boundaries, CRSobj = utmOS)

#'---------------------------------------------
# Line shapefile of submerged banks
#'---------------------------------------------

banks <- readRDS("gis/os_banks.rds")


#' ====================================
# ENVIRONMENTAL LAYERS ====
#' ====================================

#'---------------------------------------------
# Load the grid-scale data
#'---------------------------------------------

# See raster prep R file for details on how geomorphic layers
# were derived from the multibeam bathymetry grids

grids.stack <- readRDS("data/grid_stack.rds")
geom <- readr::read_csv("data/geom.csv")

#'---------------------------------------------
# Load the AMP-scale data
#'---------------------------------------------

amp.stack <- readRDS("data/AMP_stack.rds")

#' ====================================
# SPECIES DATA ====
#' ====================================

#'---------------------------------------------
# Load the data
#'---------------------------------------------

oshoals <- read.csv("data/oshoals_data.csv")
oshoals <- janitor::clean_names(oshoals)

#'---------------------------------------------
# Retrieve trophic levels (from FishBase)
#'---------------------------------------------

tl.1 <- oshoals %>% dplyr::group_by(full_name) %>% dplyr::mutate(TL = getTL(full_name))
tl.2 <- tl.1 %>% dplyr::mutate(TLgenus = ifelse(is.na(TL), getTLgenus(genusname = genus), TL))
tl.3 <- tl.2 %>% dplyr::mutate(TLfamily = ifelse(is.na(TLgenus), getTLfamily(familyname = family), TLgenus))
tl.3$TL <- tl.3$TLfamily
tl.3$TLgenus <- tl.3$TLfamily <- NULL

#'---------------------------------------------
# Manual input from published literature (see main)
#'---------------------------------------------

tl.3[tl.3$full_name=="Orcinus orca",]$TL <- 4.5 
tl.3[tl.3$full_name=="Lepidochelys olivacea",]$TL <- 3.1 
tl.3[tl.3$full_name=="Hydrophiidae sp",]$TL <- 4
tl.3[tl.3$waypoint=="WP71",]$TL <- 0 # No detections at wp71

#'---------------------------------------------
# Adjust values for species complexes
#'---------------------------------------------

tl.3[tl.3$full_name=="Carcharhinus limbatus/sorrah/amblyrhynchoides",]$TL <- 
  mean(c(getTL("Carcharhinus limbatus"),
         getTL("Carcharhinus sorrah"),
         getTL("Carcharhinus amblyrhynchoides")))

tl.3[tl.3$full_name=="Scomberomorus commerson/semifasciatus",]$TL <- 
  mean(c(getTL("Scomberomorus commerson"),
         getTL("Scomberomorus semifasciatus")))

oshoals <- dplyr::ungroup(tl.3)
rm(tl.1, tl.2, tl.3)

#'---------------------------------------------
# Calculate species richness for each site
#'---------------------------------------------

TL.threshold <- 3.5

srcalc(dframe = oshoals, name = waypoint, removeTL = TRUE, TLvalue = TL.threshold)

#'---------------------------------------------
# Save removed species to text file
#'---------------------------------------------

utils::capture.output(x = data.frame(species = removed.species[!removed.species==""]), 
                      file = "output/removed_species.txt")

#'---------------------------------------------
# Ensure MaxN = 0 when SR = 0
#'---------------------------------------------

dfres[dfres$SR==0,]$max_n <- 0

#'---------------------------------------------
# Retrieve species richness and sum(MaxN) for each
#'---------------------------------------------

srwpt <- dfres %>% 
  dplyr::mutate(wpt = waypoint) %>% 
  dplyr::group_by(wpt) %>% 
  dplyr::summarise(sr = mean(SR), max_n = sum(max_n), duration = mean(duration))

srwpt$max_n <- ifelse(is.na(srwpt$max_n), 0, srwpt$max_n)

#'---------------------------------------------
# Extract environmental data
#'---------------------------------------------

sampling.sites <- split(oshoals, as.factor(oshoals$grid))
names(sampling.sites) <- paste0("g", 1:3)

sampling.sites <- sampling.sites %>% 
  purrr::map(.x = ., .f = ~dplyr::select(., lon, lat, waypoint) %>% 
               dplyr::distinct() %>% 
               sp::SpatialPointsDataFrame(data = ., coords = cbind(.$lon, .$lat), 
                                          proj4string = wgs84) %>% 
               sp::spTransform(., CRSobj = utmOS))

raster.coordinates <- sampling.sites %>% 
  purrr::map(.x = ., .f = ~raster::as.data.frame(.x, xy=TRUE) %>% 
               dplyr::select(coords.x1, coords.x2) %>% 
               dplyr::rename(x = coords.x1) %>% 
               dplyr::rename(y = coords.x2))
  
env.data <- purrr::map2(.x = grids.stack, 
              .y = sampling.sites,
              .f = ~raster::extract(x = .x, y = .y)) %>% 
            purrr::map2(.x = ., 
              .y = raster.coordinates,
              .f = ~cbind(.x, .y))



#'---------------------------------------------
# Append to df
#'---------------------------------------------

osdata <- purrr::map(.x = sampling.sites,
                       .f = ~.x@data) %>% 
  purrr::map2(.x = .,
              .y = env.data,
              .f = ~cbind(.x, .y)) %>% 
  do.call(rbind, .) %>% 
  as_tibble(.)

#'---------------------------------------------
# Retrieve SR, MaxN, and grid info
#'---------------------------------------------

osdata <- dplyr::inner_join(x = osdata, y = srwpt, by = c('waypoint' = 'wpt'))
osdata <- dplyr::inner_join(x = osdata, y = oshoals[, c("waypoint", "grid")], by = 'waypoint')

par(mfrow=c(1,2))
plot(table(osdata$sr), main = "SR")
plot(table(osdata$max_n), main = "MaxN")

#' ====================================
# RAREFACTION ====
#' ====================================

# Need to re-run srcalc function to get species richness
# for whole assemblage 

srcalc(dframe = oshoals,
       name = all, 
       removeTL = F, 
       TLvalue = 0)

#'---------------------------------------------
# Uses package {iNEXT}
# Tutorial available online 
# @ https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html
# Input data can be either df or list of abundances
#'---------------------------------------------

species.counts <- dfres.all %>% 
  dplyr::filter(full_name%in%c(listres.all$splist$full_name)) %>% 
  dplyr::filter(., TL >= TL.threshold)

species.counts <- species.counts %>% 
  dplyr::group_by(full_name) %>% 
  dplyr::summarise(count = sum(max_n))

species.counts <- species.counts[complete.cases(species.counts),]

#'---------------------------------------------
# Rarefaction curve
#'---------------------------------------------

species.raref <- iNEXT::iNEXT(x = species.counts$count,
                              q = 0,
                              datatype = "abundance",
                              endpoint = 1000)

#'---------------------------------------------
# Estimate for 50 individuals
#'---------------------------------------------

ind.50 <- estimateD(species.counts$count, datatype = "abundance", level = 50)

#'---------------------------------------------
# Plotting the resulting curve
#'---------------------------------------------

myggiNEXT(species.raref, grey=TRUE)+
  theme(legend.position="none")+
  geom_point(data=ind.50[ind.50$order==0,], aes(x=m,y=qD),fill="white",
             colour="black", shape=23, size=4)+
  ggsidekick::theme_sleek()+
  theme(legend.position="none")+
  ylab("Species richness")+
  xlab("Number of individuals")+
  theme(axis.text.y   = element_text(size=13),
        axis.text.x   = element_text(size=13),
        axis.title.y  = element_text(size=14, 
                                     margin = unit(c(0, 5, 0, 0), "mm")),
        axis.title.x  = element_text(size=14,
                                     margin = unit(c(5, 0, 0, 0), "mm")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_x_continuous(labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE))

ggsave(filename = "output/rarefaction_curve.pdf", width = 5, height = 5)

#' ====================================
# PERMANOVA ====
#' ====================================

#'---------------------------------------------
# Prepare data for use with {vegan}
#'---------------------------------------------

oshoals.vegan <- dfres.all %>% 
  dplyr::filter(full_name%in%c(listres.all$splist$full_name)) %>% 
  dplyr::filter(., TL >= TL.threshold)

oshoals.vegan <- dplyr::select(oshoals.vegan, c("waypoint", "full_name", "max_n"))
oshoals.vegan <- reshape2::dcast(oshoals.vegan, waypoint ~ full_name,
                                 value.var = "max_n", fill = 0)

#'---------------------------------------------
# Recover deployment duration and geom factor
#'---------------------------------------------

oshoals.vegan <- dplyr::inner_join(x = oshoals.vegan,
                                   y = geom,
                                   by = c("waypoint" = "wpt"))


wpt.duration <- as_tibble(oshoals[,c("waypoint","duration")])
wpt.duration <- wpt.duration[!duplicated(wpt.duration),]

oshoals.vegan <- dplyr::inner_join(x = oshoals.vegan,
                                   y = wpt.duration,
                                   by = 'waypoint')

# oshoals.vegan <- oshoals.vegan %>% dplyr::filter(duration>180)

#'---------------------------------------------
# Response variables in a sample x species matrix
#'---------------------------------------------

oshoals.vegan.matrix <- as.matrix(oshoals.vegan[,2:(which(names(oshoals.vegan)=="lat")-1)])

#'---------------------------------------------
# Use fourth root to minimize influence of most abundant groups.
#'---------------------------------------------

oshoals.vegan.mat <- sqrt(sqrt(oshoals.vegan.matrix)) # Fourth root transform

#'---------------------------------------------
# Dissimilarity matrices
#'---------------------------------------------

oshoals.bray0 <- ecole::bray0(x = oshoals.vegan.mat) # Zero-adjusted Bray-Curtis
# oshoals.mGower5 <- vegdist(decostand(oshoals.vegan.mat, method = "log", logbase = 5, range.global = "range"), "altGower") # Anderson et al. (2006) version of Gower | Range standardization with "altGower" (that excludes double-zeros)


#'---------------------------------------------
# Run the permanova
#'---------------------------------------------

oshoals.permanova.BC0 <- vegan::adonis2(oshoals.bray0 ~ duration + geom, # Model formula
                                        data = oshoals.vegan, # df for independent vars
                                        permutations = 9999, #
                                        strata = "geom", # strata to constrain permutations
                                        by = NULL) # Type III sums of squares - most conservative
# See page 71 Anderson's manual

oshoals.permanova.BC0


#'---------------------------------------------
# PERMDISP - test of homogeneity of dispersions
#'---------------------------------------------

anv.perm <- anova(vegan::betadisper(oshoals.bray0, group = oshoals.vegan$geom))

#'---------------------------------------------
# SIMPER - similarity percentage analysis - which species contribute to differences?
#'---------------------------------------------

# https://www.tandfonline.com/doi/pdf/10.1080/00288330809509934
# When ANOSIM pairwise tests showed fish
# assemblages from two habitats differed significantly,
# analysis of similarities (SIMPER) was used to
# determine the best discriminating species. SIMPER identified a core
# group of discriminating species that contributed
# highly to both the within-habitat similarity and
# between-habitat dissimilarity, exceeding a "ratio"
# of 1.5 for each (Clarke & Warwick 2001). 

oshoals.sim <- vegan::simper(comm = oshoals.vegan.mat, 
                             group = oshoals.vegan$geom,
                             permutations = 9999)
summary(oshoals.sim)

#'---------------------------------------------
# Non-metric multidimensional scaling
#'---------------------------------------------

oshoalsMDS <- vegan::metaMDS(oshoals.bray0, # Dissimilarity matrix
                             k = 3, # Number of dimensions
                             try = 1500, # Random starts
                             autotransform = TRUE,
                             previous.best = T)

# GOF and Shepard Plot for Nonmetric Multidimensional Scaling

# Clarke (1993) suggests the following guidelines
# for acceptable stress values: <0.05 = excellent, <0.10 = good, <0.20 = usable,
# >0.20 = not acceptable

stressplot(oshoalsMDS)
nmds.stress <- oshoalsMDS$stress

mynmds <- plotnmds(nmds.object = oshoalsMDS, 
                   ndim = 2, 
                   size.pts = 8,
                   input.data = oshoals.vegan,
                   plot.factor = geom,
                   plot.ellipse = F)

utils::capture.output(x = oshoals.permanova.BC0, file = "output/permanova.txt")
utils::capture.output(x = anv.perm, file = "output/permdisp.txt")
utils::capture.output(x = summary(oshoals.sim), file = "output/simper.txt")
utils::capture.output(x = nmds.stress, file = "output/stress.txt")

#' ====================================
# EFFECTIVE DIVERSITY ====
#' ====================================

#'---------------------------------------------
# Retrieve species names
#'---------------------------------------------

lnames <- oshoals %>% 
  dplyr::filter(TL >= TL.threshold) %>% 
  dplyr::select(full_name) %>% 
  unique(.) %>% 
  dplyr::pull(.)

#'---------------------------------------------
# Filter and reshape the data
#'---------------------------------------------

esr <- dfres %>% 
  dplyr::filter(full_name%in%lnames) %>% 
  dplyr::select(waypoint, full_name, max_n)

# CHECK: sum(esr$adjMaxN) = 231

esrcast <- reshape2::dcast(data = esr, 
                           formula = waypoint ~ full_name, 
                           value.var = "max_n", 
                           fun.aggregate = sum)
rownames(esrcast) <- esrcast$waypoint
esrcast$waypoint <- NULL
esrcast <- as.matrix(esrcast)

#'---------------------------------------------
# Calculate the Shannon-Wiener index at each site
#'---------------------------------------------

shannon.index <- vegan::diversity(esrcast, index = "shannon", MARGIN = 1)

#'---------------------------------------------
# Compute ESR as exponential of SW
#'---------------------------------------------

esr.index <- as.data.frame(shannon.index) %>% 
  dplyr::mutate(waypoint = row.names(.)) %>% 
  dplyr::mutate(esr = exp(shannon.index)) %>% 
  as_tibble(.)

#'---------------------------------------------
# Add to master tibble
#'---------------------------------------------

osdata <- dplyr::inner_join(x = osdata,
                            y = esr.index,
                            by = 'waypoint')

