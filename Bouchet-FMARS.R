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

pacman::p_load(tidyverse,                  # Tidyverse for data science
               reshape2,                   # Transform data between wide and long formats
               ecole,                      # For zero-adjusted Bray Curtis
               MASS,                       # For simulating from multivariate Normal
               raster,                     # GIS
               FSSgam,                     # Full subsets multiple regression in ecology with R
               mgcv,                       # Generalised additive models
               iNEXT,                      # Species rarefaction
               vegan,                      # Community ecology: ordination & dissimilarities
               rfishbase,                  # Access to Fishbase data
               pals,                       # Colour ramps
               mgcViz,                     # Advanced visualisation for GAMs
               gratia,                     # Extra functionalities for GAMs
               ape,                        # Moran's I
               fields,                     # Euclidean distance matrices
               FSSgam,                     # Full subsets GAMs
               prettymapr,                 # Scale Bar, North Arrow, and Pretty Margins in R
               gridGraphics                # Redraw Base Graphics Using 'grid' Graphics
)

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
# Line shapefile of submerged banks
#'---------------------------------------------

banks <- readRDS("gis/os_banks.rds")

#' ====================================
# ENVIRONMENTAL LAYERS ====
#' ====================================

covariates <- c('depth', 'slope', 'curv', 'east', 'north', 'stde', 'tpiL', 'tpiS', 'bk')

#'---------------------------------------------
# Load the grid-scale data
#'---------------------------------------------

# See raster prep R file for details on how geomorphic layers
# were derived from the multibeam bathymetry grids

# raster::writeRaster(x = grids.stack$g1,
#                     file = 'data/grids.stack.g1.tif',
#                     options = 'INTERLEAVE=BAND')

grids.stack <- list()
grids.stack$g1 <- raster::stack("data/grids.stack.g1.tif")
grids.stack$g2 <- raster::stack("data/grids.stack.g2.tif")
grids.stack$g3 <- raster::stack("data/grids.stack.g3.tif")

names(grids.stack$g1) <- names(grids.stack$g2) <- names(grids.stack$g3) <- covariates

geom <- readr::read_csv("data/geom.csv")

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

osdata <- dplyr::left_join(x = osdata, y = srwpt, by = c('waypoint' = 'wpt'))
osdata <- dplyr::left_join(x = osdata, y = dplyr::distinct(oshoals[, c("waypoint", "grid")]), 
                           by = 'waypoint')

#'---------------------------------------------
# Factorise grid labels
#'---------------------------------------------

osdata$grid <- factor(osdata$grid)

#'---------------------------------------------
# Add x/y coordinates
#'---------------------------------------------

osdata <- osdata %>% 
  dplyr::select(lon, lat) %>% 
  sp::SpatialPoints(coords = ., proj4string = wgs84) %>% 
  sp::spTransform(., utmOS) %>% 
  sp::coordinates(.) %>% 
  tibble::as_tibble(.) %>% 
  dplyr::rename(x = lon, y = lat) %>% 
  dplyr::bind_cols(osdata, .) 


par(mfrow = c(1,2))
plot(table(osdata$sr), main = "SR")
plot(table(osdata$max_n), main = "MaxN")

#' ====================================
# RAREFACTION ====
#' ====================================

# Need to re-run srcalc function to get species richness
# for whole assemblage 

srcalc(dframe = oshoals,
       name = all, 
       removeTL = FALSE, 
       TLvalue = 0)

#'---------------------------------------------
# Uses package {iNEXT}
# Tutorial available online 
# @ https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html
# Input data can be either df or list of abundances
#'---------------------------------------------

species.counts <- dfres.all %>% 
  dplyr::filter(full_name%in%c(listres.all$splist$full_name)) 
# %>% 
#   dplyr::filter(., TL >= TL.threshold)

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
                              endpoint = 5000)

#'---------------------------------------------
# Estimate for 50 individuals
#'---------------------------------------------

ind.40 <- estimateD(species.counts$count, datatype = "abundance", level = 40)
ind.50 <- estimateD(species.counts$count, datatype = "abundance", level = 50)

#'---------------------------------------------
# Plotting the resulting curve
#'---------------------------------------------

myggiNEXT(species.raref, grey=TRUE)+
  theme(legend.position="none")+
  # geom_segment(data = df40, aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"))+
  # geom_segment(data = df40, aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"))+
  geom_point(data = ind.50[ind.50$order==0,], aes(x = m, y = qD), fill = "white",
             colour = "black", shape = 23, size = 4)+
  geom_point(data = ind.40[ind.40$order==0,], aes(x = m, y = qD), fill = "gray",
             colour = "black", shape = 23, size = 4)+
  ggsidekick::theme_sleek()+
  theme(legend.position = "none")+
  ylab("Species richness")+
  xlab("Number of individuals")+
  theme(axis.text.y   = element_text(size = 13),
        axis.text.x   = element_text(size = 13),
        axis.title.y  = element_text(size = 14, 
                                     margin = unit(c(0, 5, 0, 0), "mm")),
        axis.title.x  = element_text(size = 14,
                                     margin = unit(c(5, 0, 0, 0), "mm")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  scale_x_continuous(labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE))

ggsave(filename = "output/rarefaction_curve.pdf", width = 5, height = 5)

#' ====================================
# PERMANOVA ====
#' ====================================

#'---------------------------------------------
# Prepare data for use with {vegan}
#'---------------------------------------------

oshoals.vegan <- dfres.all %>% 
  dplyr::filter(full_name%in%c(listres.all$splist$full_name)) 
# %>% 
#   dplyr::filter(., TL >= TL.threshold)

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
                                        permutations = 9999,
                                        strata = "geom", # strata to constrain permutations
                                        by = "margin") # Type III sums of squares - most conservative
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
# When ANOSIM pairwise tests showed fish assemblages from two habitats differed significantly,
# analysis of similarities (SIMPER) was used to determine the best discriminating species.
# SIMPER identified a core group of discriminating species that contributed
# highly to both the within-habitat similarity and between-habitat dissimilarity, exceeding a "ratio"
# of 1.5 for each (Clarke & Warwick 2001). 
# 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3605449/
# When significant P values were obtained from pairwise tests using PERMANOVA for baited video data, similarity percentages (SIMPER) were used to identify significant distinguishing fish species. This criteria was based on dissimilarity to standard deviation ratios (Diss/SD) >2 and percentage contributions >10%. 

# https://www.frontiersin.org/articles/10.3389/fmars.2016.00241/full
# SIMPER analysis indicated that dissimilarity between communities in the two study regions was accounted for by small contributions (<2%) by a large number of taxa (Supplementary Table 4).
# 
# We considered species with %–δi >  3% and with –δi/SD(–δi) > 1 as being important contributors to dissimilarity between parks (Terlizzi et al. 2005).

oshoals.sim <- vegan::simper(comm = oshoals.vegan.mat, 
                             group = oshoals.vegan$geom,
                             permutations = 9999)
summary(oshoals.sim)

oshoals.sim$`bank_non-bank`$overall # overall between group dissimilarity


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

pdf(file = "output/stressplot.pdf")
stressplot(oshoalsMDS)
dev.off()

mynmds <- plotnmds(nmds.object = oshoalsMDS, 
                   ndim = 2, 
                   size.pts = 8,
                   input.data = oshoals.vegan,
                   plot.factor = geom,
                   plot.ellipse = F)

utils::capture.output(x = oshoals.permanova.BC0, file = "output/permanova.txt")
utils::capture.output(x = anv.perm, file = "output/permdisp.txt")
utils::capture.output(x = summary(oshoals.sim), file = "output/simper.txt")
utils::capture.output(x = oshoals.sim$`bank_non-bank`$overall, file = "output/overall_diss.txt")
utils::capture.output(x = oshoalsMDS$stress, file = "output/stress.txt")



#' ====================================
# COVARIATES ==== 
#' ====================================

#'---------------------------------------------
# Visual inspection of covariates
#'---------------------------------------------

for(p in which(names(osdata)%in%covariates)){
  par(mfrow = c(2,1))
  hist(pull(osdata[,p]), main = names(osdata)[p])
  plot(jitter(pull(osdata[,p])))
}

#'---------------------------------------------
# Univariate range of each covariate (excluding depth, bk in tensor product)
#'---------------------------------------------

envrg <- data.frame(apply(osdata[, covariates[2:8]], MARGIN = 2, FUN = function(x) range(x)))

#'---------------------------------------------
# Filter raster stacks
#'---------------------------------------------

grids.stack.filtered <- list()

grids.stack.filtered$g1 <- extrap.filter(input.raster.stack = grids.stack$g1, 
                                         covariate.rge = envrg, 
                                         dat = osdata, 
                                         convex = TRUE,
                                         concave = FALSE)

grids.stack.filtered$g2 <- extrap.filter(input.raster.stack = grids.stack$g2, 
                                         covariate.rge = envrg, 
                                         dat = osdata, 
                                         convex = TRUE,
                                         concave = FALSE)


grids.stack.filtered$g3 <- extrap.filter(input.raster.stack = grids.stack$g3, 
                                         covariate.rge = envrg, 
                                         dat = osdata,
                                         convex = TRUE,
                                         concave = FALSE)



#' ====================================
# ESR ==== 
#' ====================================

#'---------------------------------------------
# Calculates effective diversity
#'---------------------------------------------

# Extract species names (after filtering on TL)

# lnames <- purrr::map_lgl(listres, is.data.frame) %>% 
#   listres[.] %>% 
#   purrr::map(., ~dplyr::pull(., full_name)) %>% 
#   unlist()
# lnames <- sort(unique(lnames))
# 
# # Extract relevant data
# 
# os.esr <- dfres %>% dplyr::filter(full_name%in%lnames)
# os.esr <- os.esr %>% dplyr::select(., waypoint, full_name, max_n)
# 
# # Generate species x waypoint matrix
# 
# oscast <- reshape2::dcast(os.esr, waypoint~full_name, value.var = "max_n", fun.aggregate = sum) %>% 
#   tibble::as_tibble(.)
# 
# os.cast <- as.matrix(oscast[,2:ncol(oscast)])
# 
# # Calculating Shannon-wiener index for each site
# # Shannon or Shannon–Weaver (or Shannon–Wiener) index is defined as H = -sum p_i log(b) p_i, where p_i is the proportional abundance of species i and b is the base of the logarithm.
# 
# os.shannon <- vegan::diversity(os.cast, index = "shannon", MARGIN = 1)
# 
# # Effective SR is the exponential of shannon
# # http://www.loujost.com/Statistics%20and%20Physics/Diversity%20and%20Similarity/EffectiveNumberOfSpecies.htm
# # It is the first order hill number
# 
# os.esr <- exp(os.shannon) %>% 
#   tibble::tibble(esr = .) %>% 
#   dplyr::mutate(waypoint = oscast$waypoint)
# 
# # Combine with full dataset
# 
# osdata <- dplyr::left_join(x = osdata, y = os.esr, by = "waypoint")
# 
# # Replace NAs with zeros
# 
# osdata$esr <- ifelse(is.na(osdata$esr), 0, osdata$esr)
# osdata$esr <- round(osdata$esr)



#'---------------------------------------------
# Standarde the covariates
#'---------------------------------------------

# Not needed

# osdatast <- standardise.variables(df = osdata, 
#                       col.indices = which(names(osdata)%in%covariates),
#                       sd.factor = 2,
#                       append = FALSE)               
# 
# osdatast <- osdata %>% 
#   dplyr::select(lon, lat, waypoint, sr, max_n, duration, grid) %>% 
#   cbind(., osdatast) %>% 
#   as_tibble(.)

#' ====================================
# SPATIAL MODELS ==== 
#' ====================================

#'---------------------------------------------
# Fit test models
#'---------------------------------------------

gam.sr <- mgcv::gam(formula = sr ~ 
                      te(depth, bk) +
                      s(grid, bs = 're') +
                      s(curv) +
                      s(east) +
                      s(north) +
                      s(slope) +
                      s(tpiL) +
                      s(tpiS) +
                      s(stde) +
                      offset(log(duration)),
                    REML = TRUE,
                    data = osdata,
                    family = poisson(link = "log"))

gam.max_n <- mgcv::gam(formula = max_n ~ 
                         te(depth, bk) +
                         s(grid, bs = 're') +
                         s(curv) +
                         s(east) +
                         s(north) +
                         s(slope) +
                         s(tpiL) +
                         s(tpiS) +
                         s(stde) +
                         offset(log(duration)),
                       REML = TRUE,
                       data = osdata,
                       family = poisson(link = "log"))


#'---------------------------------------------
# Model validation (residual diagnostics) + spatial autocorrelation
#'---------------------------------------------

quick.validate(gam.sr, random.quantiles = TRUE)
quick.validate(gam.max_n, random.quantiles = TRUE)

gratia::appraise(gam.sr) # Does not work for Tweedie models
gratia::appraise(gam.max_n)

gratia::draw(gam.sr) # Partial plots
gratia::draw(gam.max_n) 

plot(gam.sr, scale = 0)
plot(gam.max_n, scale = 0)

coord.s <- cbind(osdata$x, osdata$y)      
w <-  fields::rdist(coord.s)  # point at which R runs into memory problems 

ape::Moran.I(x = residuals(gam.sr), w = w) 
ape::Moran.I(x = residuals(gam.max_n), w = w) 

# countreg::rootogram(object = gam.sr)

# Average distance between sites within a grid

osdata %>% 
  dplyr::filter(grid==1) %>% 
  dplyr::select(x, y) %>% 
  fields::rdist(.) %>% mean(.)

osdata %>% 
  dplyr::filter(grid==2) %>% 
  dplyr::select(x, y) %>% 
  fields::rdist(.) %>% mean(.)

osdata %>% 
  dplyr::filter(grid==3) %>% 
  dplyr::select(x, y) %>% 
  fields::rdist(.) %>% mean(.)

plot.variogram(gam.sr, osdata)
plot.variogram(gam.max_n, osdata)

# plot.tensor(gam.sr)
# plot.tensor(gam.max_n)

# DHARMa also an option

# sim_n <- simulate(test.gam.max_n, nsim = 10000)
# sim_res <- DHARMa::createDHARMa(simulatedResponse = sim_n, 
#                                observedResponse = osdata$max_n,
#                                fittedPredictedResponse = predict(test.gam.max_n, type = "response"),
#                                integerResponse = TRUE)
# 
# plot(sim_res)
# DHARMa::testSpatialAutocorrelation(simulationOutput = sim_res, x = osdata$x, y= osdata$y, alternative = "greater")


#'---------------------------------------------
# Define distributions and link functions
#'---------------------------------------------

# resp.vars.fams <- list("sr" = poisson(link = "log"), "max_n" = poisson(link = "log"))

#'---------------------------------------------
# Generate model sets
#'---------------------------------------------

model.set.sr <- list.model.sets(response.var = 'sr', 
                                dat = osdata, 
                                corr.cutoff = 0.7, 
                                max.predictors = 4)

model.set.max_n <- list.model.sets(response.var = 'max_n', 
                                   dat = osdata, 
                                   corr.cutoff = 0.7, 
                                   max.predictors = 4)

#'---------------------------------------------
# Fit all models
#'---------------------------------------------

gam.sr.list <- fit.model.set(model.set.list = model.set.sr, max.models = 500)
gam.max_n.list <- fit.model.set(model.set.list = model.set.max_n, max.models = 500)

#'---------------------------------------------
# Extract results
#'---------------------------------------------

gam.sr.fss <- FSS.results(gam.sr.list, plot.models = FALSE)
gam.max_n.fss <- FSS.results(gam.max_n.list, plot.models = FALSE)

#'---------------------------------------------
# Summarise results
#'---------------------------------------------

summary(gam.sr.fss$best.model$obj)
summary(gam.max_n.fss$best.model$obj)

quick.validate(gam.sr.fss$best.model$obj)
quick.validate(gam.max_n.fss$best.model$obj)

gratia::appraise(gam.sr.fss$best.model$obj) # Does not work for Tweedie models
gratia::appraise(gam.max_n.fss$best.model$obj)

gratia::draw(gam.sr.fss$best.model$obj) # Partial plots
gratia::draw(gam.max_n.fss$best.model$obj, select = 1)
ggsave(filename = "output/curv_sm.pdf")
gratia::draw(gam.max_n.fss$best.model$obj, select = 2) 
ggsave(filename = "output/tpiL_sm.pdf")

coord.s <- cbind(osdata$x, osdata$y)      
w <-  fields::rdist(coord.s)  # point at which R runs into memory problems 

ape::Moran.I(x = residuals(gam.sr.fss$best.model$obj), w = w) 
ape::Moran.I(x = residuals(gam.max_n.fss$best.model$obj), w = w) 

# countreg::rootogram(object = gam.sr)

plot.variogram(gam.sr.fss$best.model$obj, osdata)
plot.variogram(gam.max_n.fss$best.model$obj, osdata)

plot.tensor(gam.sr.fss$best.model$obj)
plot.tensor(gam.sr.fss$success.models$depth.te.bk)
# ggsave(filename = "output/tensor_sr.pdf")
# 
plot.tensor(gam.max_n.fss$best.model$obj)
plot.tensor(gam.max_n.fss$success.models$depth.te.bk)
# ggsave(filename = "output/tensor_max_n.pdf")

#' ====================================
# PREDICTIONS ==== 
#' ====================================

#'---------------------------------------------
# Generate predictions from the best models
#'---------------------------------------------

preds.sr <- purrr::map(.x = grids.stack.filtered,
                       .f = ~get.predictions(model = gam.sr.fss$best.model$obj, 
                                             proj.syst = utmOS, 
                                             pred.pts = .x,
                                             verbose = TRUE))

preds.max_n <- purrr::map(.x = grids.stack.filtered,
                          .f = ~get.predictions(model = gam.max_n.fss$best.model$obj, 
                                                proj.syst = utmOS, 
                                                pred.pts = .x,
                                                verbose = TRUE))

plot(preds.sr$g1)
plot(preds.sr$g2)
plot(preds.sr$g3)

plot(preds.max_n$g1)
plot(preds.max_n$g2)
plot(preds.max_n$g3)

#' ====================================
# BOOTSTRAP ==== 
#' ====================================

n.iter <- 100

#'-------------------------------------------------
# Nest the data
#'-------------------------------------------------

os.nest <- osdata %>% 
  dplyr::group_by(grid) %>% 
  tidyr::nest()

#'-------------------------------------------------
# Add names to list elements
#'-------------------------------------------------

os.nest$data <- os.nest$data %>% 
  purrr::set_names(., levels(osdata$grid))

#'-------------------------------------------------
# Generate bootstrap resamples, by group
#'-------------------------------------------------

os.boots <- purrr::map(.x = 1:n.iter, 
                       .f = ~purrr::map2(.x = os.nest$data, 
                                         .y = as.factor(c("1", "2", "3")),
                                        .f = ~dplyr::sample_n(tbl = ., 
                                                              size = nrow(.), 
                                                              replace = T # With replacement
                                        ) %>% dplyr::mutate(grid = .y))) %>% 
  purrr::map(.x = ., .f = ~do.call(rbind, .x))


#' ====================================
# MODEL FREQUENCIES ==== 
#' ====================================

#'-------------------------------------------------
# Run full subset on bootstrap datasets and record the best models
#'-------------------------------------------------

wenger.probs <- compute.mod.probs(bootstrap.data = os.boots)

# saveRDS(object = wenger.probs, file = "data/wenger_probs.rds")
wenger.probs <- readRDS(file = "data/wenger_probs.rds")

#'-------------------------------------------------
# Calculate model frequencies
#'-------------------------------------------------

wenger.probs <- purrr::map(.x = wenger.probs, 
                           .f = ~.x[!.x=="null"] %>% 
                             table(.) %>% 
                             tibble::as_tibble(.) %>% 
                             dplyr::arrange(-n) %>% 
                             dplyr::mutate(prob = n/sum(n)) %>% 
                             dplyr::rename(modname = '.'))

#'-------------------------------------------------
# Sample models according to their frequencies
#'-------------------------------------------------

wenger.mods <- purrr::map(.x = wenger.probs, 
                          .f = ~sample(x = .x$modname, 
                                       size = n.iter, 
                                       replace = TRUE, 
                                       prob = .x$prob))


#' ====================================
# UNCERTAINTY ==== 
#' ====================================

set.seed(32) # Reproducible results

#'-------------------------------------------------
# Fit the selected models to bootstrap data  
#'-------------------------------------------------

pb <- dplyr::progress_estimated(n.iter)

boot.mods.sr <- purrr::map2(.x = os.boots,
                            .y = wenger.mods$sr,
                            .f = ~{
                              pb$tick()$print()  
                              gam.form <- model.set.sr$mod.formula[[.y]] %>% 
                                as.character(.)
                              gam.form <- as.formula(paste0("sr~", gam.form[2]))
                              
                              mgcv::gam(formula = gam.form, 
                                        REML = TRUE,
                                        data = .x,
                                        family = poisson(link = 'log'))
                            })



pb <- dplyr::progress_estimated(n.iter)

boot.mods.max_n <- purrr::map2(.x = os.boots,
                               .y = wenger.mods$max_n,
                               .f = ~{
                                 pb$tick()$print()  
                                 gam.form <- model.set.max_n$mod.formula[[.y]] %>% 
                                   as.character(.)
                                 gam.form <- as.formula(paste0("max_n~", gam.form[2]))
                                 
                                 mgcv::gam(formula = gam.form, 
                                           REML = TRUE,
                                           data = .x,
                                           family = poisson(link = 'log'))
                               })


#'-------------------------------------------------
# Generate predictions including parameter uncertainty
#'-------------------------------------------------

# SR

pb <- dplyr::progress_estimated(n.iter * 3)

boot.preds.sr <- purrr::cross2(.x = boot.mods.sr, 
                               .y = grids.stack.filtered) %>% 
  purrr::map(.x = ., 
             .f = ~ {
               pb$tick()$print()
               get.predictions(model = .x[[1]], 
                               proj.syst = utmOS,
                               pred.pts = .x[[2]],
                               sample.coef = TRUE,
                               verbose = FALSE)})
                            
# MaxN

pb <- dplyr::progress_estimated(n.iter * 3)

boot.preds.max_n <- purrr::cross2(.x = boot.mods.max_n, 
                                  .y = grids.stack.filtered) %>% 
  purrr::map(.x = ., 
             .f = ~ {
               pb$tick()$print()
               get.predictions(model = .x[[1]], 
                               proj.syst = utmOS,
                               pred.pts = .x[[2]],
                               sample.coef = TRUE,
                               verbose = FALSE)})

#'-------------------------------------------------
# Stack predictions for each grid
#'-------------------------------------------------

boot.preds.sr$g1 <- boot.preds.sr[1:100] %>% raster::stack(.)
boot.preds.sr$g2 <- boot.preds.sr[101:200] %>% raster::stack(.)
boot.preds.sr$g3 <- boot.preds.sr[201:300] %>% raster::stack(.)
boot.preds.sr[1:300] <- NULL

boot.preds.max_n$g1 <- boot.preds.max_n[1:100] %>% raster::stack(.)
boot.preds.max_n$g2 <- boot.preds.max_n[101:200] %>% raster::stack(.)
boot.preds.max_n$g3 <- boot.preds.max_n[201:300] %>% raster::stack(.)
boot.preds.max_n[1:300] <- NULL

#'-------------------------------------------------
# Calculating robust statistics (median and robust CV)
#'-------------------------------------------------

# rCV = nIQR/median, with Normalised IQR (nIQR) is IQR * 0.7413, which converts to an estimate of SD
# https://arxiv.org/pdf/1907.01110.pdf

boot.sr <- boot.max_n <- list()

boot.sr$median <- purrr::map(.x = boot.preds.sr,
                             .f = ~raster::calc(x = .x, fun = median))

boot.max_n$median  <- purrr::map(.x = boot.preds.max_n,
                             .f = ~raster::calc(x = .x, fun = median))                     

boot.sr$rcv <- purrr::map(.x = boot.preds.sr,
                          .f = ~raster::calc(x = .x, 
                                             fun = function(x){0.7413*IQR(x, na.rm = TRUE)/median(x, na.rm = TRUE)}))

boot.max_n$rcv <- purrr::map(.x = boot.preds.max_n,
                             .f = ~raster::calc(x = .x, 
                                                fun = function(x){0.7413*IQR(x, na.rm = TRUE)/median(x, na.rm = TRUE)}))                        
                            

#' ====================================
# BIVARIATE ==== 
#' ====================================

biv.rasters <- purrr::map2(.x = preds.sr, 
                      .y = preds.max_n,
                      .f = ~bivariate.index(percentile = 10,
                                            raster.one = .x,
                                            raster.two = .y,
                                            variable.one = 'sr',
                                            variable.two = 'max_n'))

biv.rasters.boot <- purrr::map2(.x = boot.sr$median, 
                           .y = boot.max_n$median,
                           .f = ~bivariate.index(percentile = 10,
                                                 raster.one = .x,
                                                 raster.two = .y,
                                                 variable.one = 'sr',
                                                 variable.two = 'max_n'))

b <- purrr::map(.x = 1:3,
           .f = ~{
             biv.rasters.boot[[paste0('g', .x)]]$bivariate %>% 
             raster::as.data.frame(., xy = TRUE) %>% 
               ggplot2::ggplot(data = ., aes(x, y, fill = z)) +
               geom_tile(colour = 'black') +
               scale_fill_gradientn(colours = pals::viridis(10)) +
               theme_minimal() + labs(x = "", y = "") +
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())})

cowplot::plot_grid(b[[1]], b[[2]], b[[3]], ncol = 3)
# ggsave(filename = "output/biv_distr.pdf", height = 2.5, width = 10)


#' ====================================
# MAPPING ==== 
#' ====================================

#'---------------------------------------------
# Get raster summaries (min, max) to help parameterise plot.grid function
#'---------------------------------------------

sr.fivenum <- purrr::map(.x = preds.sr, .f = ~sumrast(.x)) %>% 
  purrr::map(.x = ., .f = ~as.list(.x))

sr.fivenum.boot <- purrr::map(.x = boot.sr$median, .f = ~sumrast(.x)) %>% 
  purrr::map(.x = ., .f = ~as.list(.x))

sr.fivenum.cv <- purrr::map(.x = boot.sr$rcv, .f = ~sumrast(.x)) %>% 
  purrr::map(.x = ., .f = ~as.list(.x))

max_n.fivenum <- purrr::map(.x = preds.max_n, .f = ~sumrast(.x)) %>% 
  purrr::map(.x = ., .f = ~as.list(.x))

max_n.fivenum.boot <- purrr::map(.x = boot.max_n$median, .f = ~sumrast(.x)) %>% 
  purrr::map(.x = ., .f = ~as.list(.x))

max_n.fivenum.cv <- purrr::map(.x = boot.max_n$rcv, .f = ~sumrast(.x)) %>% 
  purrr::map(.x = ., .f = ~as.list(.x))

#'---------------------------------------------
# Map predictions from the best models
#'---------------------------------------------

mapcol <- sapply(X = pals::parula(100), FUN = function(x) hexa2hex(input.colour = x, transparency = 0.8))

plot.grid(input.raster = biv.rasters.boot$g1$raster, 
          response.variable = "sr", 
          log.response = FALSE,
          grid.number = 1, 
          add.banks = TRUE, 
          add.bruvs = FALSE,
          zmin = 0,
          zmax = 1,
          z.int = 0.1, 
          skew = 0,
          max.pt.size = 3,
          col.palette = mapcol,
          save.to.pdf = TRUE,
          file.name = "bivr")

plot.grid(input.raster = biv.rasters.boot$g2$raster, 
          response.variable = "sr", 
          log.response = FALSE,
          grid.number = 2, 
          add.banks = TRUE, 
          add.bruvs = FALSE,
          zmin = 0,
          zmax = 1,
          z.int = 0.1, 
          skew = 0,
          max.pt.size = 3,
          col.palette = mapcol,
          save.to.pdf = TRUE,
          file.name = "bivr")

plot.grid(input.raster = biv.rasters.boot$g3$raster, 
          response.variable = "sr", 
          log.response = FALSE,
          grid.number = 3, 
          add.banks = TRUE, 
          add.bruvs = FALSE,
          zmin = 0,
          zmax = 1,
          z.int = 0.1, 
          skew = 0,
          max.pt.size = 3,
          col.palette = mapcol,
          save.to.pdf = TRUE,
          file.name = "bivr")


# Grid 1 ====

plot.grid(input.raster = preds.sr$g1, 
          response.variable = "sr", 
          log.response = FALSE,
          grid.number = 1, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = sr.fivenum$g1$Max.,
          z.int = 0.2, 
          skew = 0,
          max.pt.size = 3,
          col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
          save.to.pdf = TRUE,
          file.name = "preds")

plot.grid(input.raster = preds.max_n$g1, 
          response.variable = "max_n", 
          log.response = TRUE,
          grid.number = 1, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = max_n.fivenum$g1$Max.,
          z.int = 1,
          skew = 0,
          col.palette = rev(brewer.spectral(100)[1:100]),
          save.to.pdf = TRUE,
          file.name = "preds")

# ................................

plot.grid(input.raster = boot.sr$median$g1, 
          response.variable = "sr", 
          log.response = FALSE,
          grid.number = 1, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = sr.fivenum.boot$g1$Max.,
          z.int = 0.2, 
          skew = 0,
          max.pt.size = 3,
          col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
          save.to.pdf = TRUE,
          file.name = "preds_boot")

plot.grid(input.raster = boot.max_n$median$g1, 
          response.variable = "max_n", 
          log.response = TRUE,
          grid.number = 1, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = max_n.fivenum.boot$g1$Max.,
          z.int = 5,
          skew = 0,
          col.palette = rev(brewer.spectral(100)[10:100]),
          save.to.pdf = TRUE,
          file.name = "preds_boot")

# ................................

plot.grid(input.raster = boot.sr$rcv$g1, 
          response.variable = "sr", 
          log.response = FALSE,
          cv.over.one = F,
          grid.number = 1, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = 1,
          z.int = 0.2, 
          skew = 0,
          max.pt.size = 3,
          col.palette = pals::viridis(100),
          save.to.pdf = TRUE,
          file.name = "preds_cv")

plot.grid(input.raster = boot.max_n$rcv$g1, 
          response.variable = "max_n", 
          log.response = FALSE,
          cv.over.one = TRUE,
          grid.number = 1, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = 1,
          z.int = 0.2,
          skew = 0,
          over.one.col = "darkorange",
          col.palette = pals::viridis(100),
          save.to.pdf = TRUE,
          file.name = "preds_cv")


# Grid 2 ====

plot.grid(input.raster = preds.sr$g2, 
          response.variable = "sr", 
          log.response = FALSE,
          grid.number = 2, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = sr.fivenum$g2$Max.,
          z.int = 0.2, 
          skew = 0,
          max.pt.size = 3,
          col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
          save.to.pdf = TRUE,
          file.name = "preds")


plot.grid(input.raster = preds.max_n$g2, 
          response.variable = "max_n", 
          log.response = TRUE,
          grid.number = 2, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = max_n.fivenum$g2$Max.,
          z.int = 5,
          skew = 0,
          col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
          save.to.pdf = TRUE,
          file.name = "preds")

# ................................

plot.grid(input.raster = boot.sr$median$g2, 
          response.variable = "sr", 
          log.response = FALSE,
          grid.number = 2, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = sr.fivenum.boot$g2$Max.,
          z.int = 0.2, 
          skew = 0,
          max.pt.size = 3,
          col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
          save.to.pdf = TRUE,
          file.name = "preds_boot")

plot.grid(input.raster = boot.max_n$median$g2, 
          response.variable = "max_n", 
          log.response = TRUE,
          grid.number = 2, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = max_n.fivenum.boot$g2$Max.,
          z.int = 5,
          skew = 0,
          col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
          save.to.pdf = TRUE,
          file.name = "preds_boot")

# ................................

plot.grid(input.raster = boot.sr$rcv$g2, 
          response.variable = "sr", 
          log.response = FALSE,
          grid.number = 2, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = sr.fivenum.cv$g2$Max.,
          z.int = 0.2, 
          skew = 0,
          max.pt.size = 3,
          col.palette = pals::viridis(100),
          save.to.pdf = TRUE,
          file.name = "preds_cv")

plot.grid(input.raster = boot.max_n$rcv$g2, 
          response.variable = "max_n", 
          log.response = FALSE,
          cv.over.one = TRUE, 
          over.one.col = "darkorange",
          grid.number = 2, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = 1,
          z.int = 0.2,
          skew = 0,
          col.palette = pals::viridis(100),
          save.to.pdf = TRUE,
          file.name = "preds_cv")

# Grid 3 ====

plot.grid(input.raster = preds.sr$g3, 
          response.variable = "sr", 
          log.response = FALSE,
          grid.number = 3, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = sr.fivenum$g3$Max.,
          z.int = 0.2, 
          skew = 0,
          max.pt.size = 3,
          col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
          save.to.pdf = TRUE,
          file.name = "preds")


plot.grid(input.raster = preds.max_n$g3, 
          response.variable = "max_n", 
          log.response = TRUE,
          grid.number = 3, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = max_n.fivenum$g3$Max.,
          z.int = 5,
          skew = 0,
          col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
          save.to.pdf = TRUE,
          file.name = "preds")

# ................................

plot.grid(input.raster = boot.sr$median$g3, 
          response.variable = "sr", 
          log.response = FALSE,
          grid.number = 3, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = sr.fivenum.boot$g3$Max.,
          z.int = 0.2, 
          skew = 0,
          max.pt.size = 3,
          col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
          save.to.pdf = TRUE,
          file.name = "preds_boot")

plot.grid(input.raster = boot.max_n$median$g3, 
          response.variable = "max_n", 
          log.response = TRUE,
          grid.number = 3, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = max_n.fivenum.boot$g3$Max.,
          z.int = 5,
          skew = 0,
          col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
          save.to.pdf = TRUE,
          file.name = "preds_boot")

# ................................

plot.grid(input.raster = boot.sr$rcv$g3, 
          response.variable = "sr", 
          log.response = FALSE,
          grid.number = 3, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = 1,
          z.int = 0.2, 
          skew = 0,
          max.pt.size = 3,
          col.palette = pals::viridis(100),
          save.to.pdf = TRUE,
          file.name = "preds_cv")

plot.grid(input.raster = boot.max_n$rcv$g3, 
          response.variable = "max_n", 
          cv.over.one = TRUE,
          over.one.col = "darkorange",
          log.response = FALSE,
          grid.number = 3, 
          add.banks = TRUE, 
          add.bruvs = TRUE,
          zmin = 0,
          zmax = 1,
          z.int = 0.2, 
          skew = 0,
          col.palette = pals::viridis(100),
          save.to.pdf = TRUE,
          file.name = "preds_cv")


# S and N per grid

srcalc(dframe = oshoals, name = grid, removeTL = F, TLvalue = 0)

cat(paste0("S1: ", listres$spcount.1))
cat(paste0("S2: ", listres$spcount.2))
cat(paste0("S3: ", listres$spcount.3))

dfres %>% filter(grid==1) %>% filter(full_name%in%listres$splist.1$full_name) %>% summarise(N = sum(max_n))
dfres %>% filter(grid==2) %>% filter(full_name%in%listres$splist.2$full_name) %>% summarise(N = sum(max_n))
dfres %>% filter(grid==3) %>% filter(full_name%in%listres$splist.3$full_name) %>% summarise(N = sum(max_n))

# Variation in S and N across sites

srcalc(dframe = oshoals, name = waypoint, removeTL = F, TLvalue = 0)

fournums <- function(x){list(minimum = min(x, na.rm = TRUE), maximum = max(x, na.rm = TRUE), avg = mean(x, na.rm = TRUE), std = sd(x, na.rm = TRUE))}

listres[grepl(pattern = "spcount.", x = names(listres))] %>% unlist() %>% fournums()

dfres %>% group_by(waypoint) %>% summarise(max_n = sum(max_n)) %>% pull(max_n) %>% fournums()


oshoals %>% group_by(full_name) %>% summarise(tot = sum(max_n)) %>% arrange(-tot) %>% mutate(ctot = cumsum(tot)) %>% mutate(cperc = ctot/sum(tot, na.rm=T))

oshoals %>% group_by(full_name) %>% summarise(prev = length(unique(waypoint))) %>% arrange(-prev) %>% mutate(perc = prev/116)

# Number of species per grid

n1 <- dfres %>% filter(grid==1) %>% pull(full_name) %>% unique() %>% sort()
n1 <- n1[!n1%in%c("Carangidae sp", "Carcharhinus sp")]

n2 <- dfres %>% filter(grid==2) %>% pull(full_name) %>% unique() %>% sort()
n2 <- n2[!n2%in%c("Carangidae sp", "Carcharhinus sp", "Sphyrna sp")]

n3 <- dfres %>% filter(grid==3) %>% pull(full_name) %>% unique() %>% sort()
n3 <- n3[!n3%in%c("Carangidae sp", "Carcharhinus sp", "", "Carangoides sp", "Monacanthidae sp", "Psenes sp")]

oshoals %>% filter(grid==1) %>% filter(full_name%in%n1) %>% summarise(max_n = sum(max_n))
oshoals %>% filter(grid==2) %>% filter(full_name%in%n2) %>% summarise(max_n = sum(max_n))
oshoals %>% filter(grid==3) %>% filter(full_name%in%n3) %>% summarise(max_n = sum(max_n))



dfres %>% 
  dplyr::mutate(wpt = waypoint) %>% 
  dplyr::group_by(wpt) %>% 
  dplyr::summarise(sr = mean(SR), max_n = sum(max_n), duration = mean(duration)) %>% pull(sr) %>% max(., na.rm=T)


