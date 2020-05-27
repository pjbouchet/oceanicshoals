#' -------------------------------------------------------------------------------
#' -------------------------------------------------------------------------------
#'
#' R code accompanying the article: Submerged carbonate banks aggregate pelagic 
#' megafauna in offshore tropical Australia
#' 
#' -------------------------------------------------------------------------------
#' -------------------------------------------------------------------------------

# Authors: Bouchet PJ, Letessier TB, Caley JM, Nichol SL, Hemmi JM, Meeuwig JJ
# Journal: Frontiers in Marine Science (2020)

#'---------------------------------------------
# Functions to calculate overdispersion
#'---------------------------------------------

# https://cran.r-project.org/web/packages/bbmle/vignettes/quasi.pdf

overd <- function(object) {
  with(object,sum((weights * residuals^2)[weights > 0])/df.residual)
}

ovd <- function(model){
  p.res <- resid(model, type = "pearson")
  ovd <- sum(p.res^2)/(model$df.res)
  return(ovd)
}

#'---------------------------------------------
# Function to retrieve feeding guild from FishBase
#'---------------------------------------------

get.guild <- function(speciesname){
  
  guild <- rfishbase::ecology(speciesname)
  return(guild$FeedingType)
  
}

#'---------------------------------------------
# Wrapper around RFishbase functions to retrieve the trophic
# level of species of interest (takes mean of genus or family if necessary)
#'---------------------------------------------

getTL<-function(speciesname){
  
  troph <- rfishbase::ecology(speciesname, 
                              fields = c("DietTroph","FoodTroph"))
  
  troph <- troph %>% 
    dplyr::mutate(troph_val= ifelse(is.na(troph$FoodTroph),
                                    troph$DietTroph,
                                    troph$FoodTroph))
  result <- troph$troph_val
  return(result)}

getTLgenus<-function(genusname){
  
  # List of FAO areas where all species of the genus are known to occur
  fao_sp <- rfishbase::faoareas(species_list(Genus=genusname))
  
  # Extracts sp from eastern IO and Wester Central Pacific
  fao_sp <- unique(fao_sp[fao_sp$Status%in%c("native", "endemic") & fao_sp$AreaCode%in%c(57,71),]$Species)
  
  troph <- rfishbase::ecology(fao_sp, fields = c("DietTroph","FoodTroph"))
  
  troph <- troph %>% 
    dplyr::mutate(troph_val= ifelse(is.na(troph$FoodTroph),
                                    troph$DietTroph,
                                    troph$FoodTroph))
  
  result <- mean(troph$troph_val, na.rm=T)
  return(result)}

getTLfamily<-function(familyname){
  
  # List of FAO areas where all species of the genus are known to occur
  fao_sp <- rfishbase::faoareas(species_list(Family = familyname))
  
  # Extracts sp from eastern IO and Wester Central Pacific
  fao_sp <- unique(fao_sp[fao_sp$Status%in%c("native", "endemic") & fao_sp$AreaCode%in%c(57,71),]$Species)
  
  troph <- rfishbase::ecology(fao_sp, fields = c("DietTroph","FoodTroph"))
  troph <- troph %>% 
    dplyr::mutate(troph_val= ifelse(is.na(troph$FoodTroph),
                                    troph$DietTroph,
                                    troph$FoodTroph))
  
  result <- mean(troph$troph_val, na.rm=T)
  return(result)}

#'---------------------------------------------
# Function to calculate species richness at different hierarchy levels
#'---------------------------------------------

srcalc <- function(dframe, name, removeTL, TLvalue){
  
  # Defactorises ---
  
  dframe$full_name<-as.character(dframe$full_name)
  dframe$common_name<-as.character(dframe$common_name)
  
  if(removeTL==TRUE){
    
    assign(x = "removed.species", 
           value = sort(unique(as.character(dframe[dframe$TL<TLvalue,]$full_name))),
           env = .GlobalEnv)
    
    cat("Removing species:")
    cat("\n")
    cat("------------------------------")
    cat("\n")
    cat(paste(removed.species[!removed.species==""], collapse="\n"))
    cat("\n")
    cat(paste("------------------------------", collapse="\n"))
    cat("\n")
    
    dframe[dframe$TL<TLvalue,]$full_name <- ""
    dframe[dframe$TL<TLvalue,]$family <- ""
    dframe[dframe$TL<TLvalue,]$genus <- ""
    dframe[dframe$TL<TLvalue,]$species <- ""
    dframe[dframe$TL<TLvalue,]$common_name <- ""
    
    # Also need to set the corresponding MaxN to 0
    dframe[dframe$TL<TLvalue,]$max_n <- 0
  }
  
  # Level at which species richness should be calculated
  
  name <- deparse(substitute(name))
  
  # Total number species for whole dataset (all sites)
  
  if(name=="all"){coln <- 0}else{coln <- which(names(dframe)==name)}
  
  print("Calculating richness values ...")
  
  # INPUTS:
  # coln = index of column containing the variable to be used to split the data
  # e.g. grid or region or waypoint
  
  listres <- listres.all <- list() # Defines list container
  dfres <- dfres.all <- tibble() # Defines a data.table container
  
  # Using conditional statements enables the user to get a species list for the entire dataset by setting coln to 0 when calling the function
  if(coln==0){lvls<-data.frame(1)}else{# Determines the levels of the variable used to split the data
    lvls <- data.frame(unique(dframe[,coln]))}
  
  # For each level (grouping) of the variable of interest e.g. waypoint ID
  for (x in 1:nrow(lvls)){
    
    # If species list for entire dataset is required, assigns data to temp
    # Otherwise, determines which records correspond to the chosen level of the variable 
    if(coln==0){
      temp <- dframe
    }else{
      pos <- as.vector(as.data.frame(dframe[,coln])==as.character(lvls[x,]))
      
      # Filters data by level
      temp <- dframe[pos,]}
    
    # Calculates family/genus count totals
    # = Number of times any family/genus appears in the data
    
    fam.c <-  temp %>% dplyr::group_by(family) %>% count()
    gen.c <-  temp %>% dplyr::group_by(genus) %>% count()
    
    # fam.c<-temp[,length(Project), by=c("Family")]
    # gen.c<-temp[,length(Project), by=c("Genus")]
    
    # Fills in the original table with the above values
    temp$fam.c<-as.numeric(as.data.frame(apply(temp,1,function(e)
    {fam.c[fam.c$family==e[which(names(temp)=="family")],2]})))
    
    temp$gen.c<-as.numeric(as.data.frame(apply(temp,1,function(e)
    {gen.c[gen.c$genus==e[which(names(temp)=="genus")],2]})))
    
    
    # INCLUSION RULES:
    # (1) If an animal has a CAAB code, then it is a unique species: "SP"
    # (2) If no CAAB code, but family only appears once, then classified as species 
    # (e.g. sea snake in Timor was the only representative of its family)
    # (3) If no CAAB code, family appears only once, but genus is nondescript: "FAM"
    # (4) If no CAAB code, family appears more than once but genus only once: "SP"
    # (5) If no CAAB code, family and genus appear more than once - then count the number of
    # unique full names within that range. If this number is one, then "SP", if greater than 1
    # then "GEN"
    # Example of this: Decapterus sp. sighted at multiple waypoints within one grid hence
    # multiple rows of Decapterus sp, but no other Decapterus species identified within that grid.
    # Decapterus sp should therefore be counted as a different species.
    
    temp <- data.frame(temp)
    
    for(i in 1:nrow(temp)){
      if(is.na(temp[i,"code"])==FALSE) {temp$sp.all[i] <- "SP" # Has CAAB code
      }else{
        if(temp$fam.c[i]==1){temp$sp.all[i] <- "SP"
        }else{
          if(temp$genus[i]==""){temp$sp.all[i] <- "FAM"
          }else{
            if(temp$gen.c[i]==1){temp$sp.all[i] <- "SP"
            }else{
              
              name_temp <- temp[temp$family==temp$family[i] & temp$genus==temp$genus[i],]
              full_name_count <- length(unique(name_temp$full_name))
              
              if(full_name_count==1){temp$sp.all[i] <- "SP"
              }else{
                if(grepl("sp", temp$species[i])==TRUE){
                  temp$sp.all[i]<-"GEN"
                }else{temp$sp.all[i]<-"SP"}
              }
            }}}}}
    
    # Adjust MaxN values
    
    temp$max_n <- ifelse(!temp$sp.all=="SP", 0, temp$max_n)

    # Use the species inclusion rules to filter out all the non-"SP" records
    # Retrieves the list of species and counts how many there are.
    # Outputs everything in an alphabetically sorted data.table
    
    tt <- data.frame(unique(temp[temp$sp.all=="SP" & temp$full_name!="" & temp$full_name!="Unknown sp",c("full_name","common_name")]))
    
    tt <- tt[order(rank(tt$full_name)),]
    
    l <- length(sort(unique(temp[temp$sp.all=="SP" & temp$full_name!="" & temp$full_name!="Unknown sp",]$full_name)))
    
    
    # Combines all results
    results <- list(tt,l)
    
    # Assigns names to list elements
    if(coln==0){
      
      names(results) <- c("splist", "spcount")
      
      listres.all <- c(listres.all, results) # Combines results
      
      # Now brings these results back into the original data table
      temp$SR <- l
      temp <- data.frame(temp)
      dfres.all <- rbind(dfres.all, temp) # Combines results
      
    }else{
      
      names(results) <- c(paste0("splist.",as.character(lvls[x,])),
                          paste0("spcount.",as.character(lvls[x,])))
      
      
      listres <- c(listres, results) # Combines results
      
      # Now brings these results back into the original data table
      temp$SR <- l
      temp <- data.frame(temp)
      dfres <- rbind(dfres, temp) # Combines results
      
    }
    
  }
  
  
  if(coln==0){
    print("Saving results to global environment [listres.all, dfres.all] ...")
    assign("listres.all", listres.all, envir = .GlobalEnv)
    assign("dfres.all", dfres.all,envir = .GlobalEnv)
    
  }else{
    print("Saving results to global environment [listres, dfres] ...")
    assign("listres", listres, envir = .GlobalEnv) # Brings the variable into the global env
    assign("dfres", dfres,envir = .GlobalEnv)}
  
}

#'---------------------------------------------
# Function to produce rarefaction curves
#'---------------------------------------------

# Package iNEXT has a ggplot wrapper function for plotting rarefaction/extrapolation curves
# However, currently returns an error about factor levels being duplicated.
# The below is a slight tweak using the unique() command
# See https://github.com/JohnsonHsieh/iNEXT/issues/25

myggiNEXT<-function (x, type = 1, se = TRUE, 
                     facet.var = "none", 
                     color.var = "site",
                     grey = FALSE)
{
  TYPE <- c(1, 2, 3)
  SPLIT <- c("none", "order", "site", "both")
  if (is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  if (is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var, SPLIT) ==
      -1)
    stop("invalid facet variable")
  if (is.na(pmatch(color.var, SPLIT)) | pmatch(color.var, SPLIT) ==
      -1)
    stop("invalid color variable")
  type <- pmatch(type, 1:3)
  facet.var <- match.arg(facet.var, SPLIT)
  color.var <- match.arg(color.var, SPLIT)
  if (facet.var == "order")
    color.var <- "site"
  if (facet.var == "site")
    color.var <- "order"
  options(warn = -1)
  z <- fortify(x, type = type)
  options(warn = 0)
  if (ncol(z) == 7) {
    se <- FALSE
  }
  datatype <- unique(z$datatype)
  if (color.var == "none") {
    if (levels(factor(z$order)) > 1 & "site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites and orders, change setting as both")
      color.var <- "both"
      z$col <- z$shape <- paste(z$site, z$order, sep = "-")
    }
    else if ("site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object consists multiple orders, change setting as order")
      color.var <- "site"
      z$col <- z$shape <- z$site
    }
    else if (levels(factor(z$order)) > 1) {
      warning("invalid color.var setting, the iNEXT object consists multiple sites, change setting as site")
      color.var <- "order"
      z$col <- z$shape <- factor(z$order)
    }
    else {
      z$col <- z$shape <- rep(1, nrow(z))
    }
  }
  else if (color.var == "order") {
    z$col <- z$shape <- factor(z$order)
  }
  else if (color.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- z$site
  }
  else if (color.var == "both") {
    if (!"site" %in% names(z)) {
      warning("invalid color.var setting, the iNEXT object do not consist multiple sites, change setting as order")
      z$col <- z$shape <- factor(z$order)
    }
    z$col <- z$shape <- paste(z$site, z$order, sep = "-")
  }
  z$lty <- factor(z$method, levels=unique(c("interpolated", "observed", "extrapolated"),
                                          c("interpolation", "interpolation", "extrapolation")))
  z$col <- factor(z$col)
  data.sub <- z[which(z$method == "observed"), ]
  g <- ggplot(z, aes_string(x = "x", y = "y", colour = "col")) +
    geom_point(aes_string(shape = "shape"), size = 5, data = data.sub)
  g <- g + geom_line(aes_string(linetype = "lty"), lwd = 0.8) +
    guides(linetype = guide_legend(title = "Method"), colour = guide_legend(title = "Guides"),
           fill = guide_legend(title = "Guides"), shape = guide_legend(title = "Guides")) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          text = element_text(size = 18))
  if (type == 2L) {
    g <- g + labs(x = "Number of sampling units", y = "Sample coverage")
    if (datatype == "abundance")
      g <- g + labs(x = "Number of individuals", y = "Sample coverage")
  }
  else if (type == 3L) {
    g <- g + labs(x = "Sample coverage", y = "Species diversity")
  }
  else {
    g <- g + labs(x = "Number of sampling units", y = "Species diversity")
    if (datatype == "abundance")
      g <- g + labs(x = "Number of individuals", y = "Species diversity")
  }
  if (se)
    g <- g + geom_ribbon(aes_string(ymin = "y.lwr", ymax = "y.upr",
                                    fill = "factor(col)", colour = "NULL"), alpha = 0.1)
  if (facet.var == "order") {
    if (length(levels(factor(z$order))) == 1 & type != 2) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple orders.")
    }
    else {
      g <- g + facet_wrap(~order, nrow = 1)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides",
                                              ncol = length(levels(factor(z$order))), byrow = TRUE),
                        fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (facet.var == "site") {
    if (!"site" %in% names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites.")
    }
    else {
      g <- g + facet_wrap(~site, nrow = 1)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides",
                                              nrow = length(levels(factor(z$order)))), fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (facet.var == "both") {
    if (length(levels(factor(z$order))) == 1 | !"site" %in%
        names(z)) {
      warning("invalid facet.var setting, the iNEXT object do not consist multiple sites or orders.")
    }
    else {
      g <- g + facet_wrap(site ~ order)
      if (color.var == "both") {
        g <- g + guides(colour = guide_legend(title = "Guides",
                                              nrow = length(levels(factor(z$site))), byrow = TRUE),
                        fill = guide_legend(title = "Guides"))
      }
    }
  }
  if (grey) {
    g <- g + theme_bw(base_size = 18) + 
      scale_fill_grey(start = 0, end = 0.2) +
      scale_colour_grey(start = 0, end = 0.2) +
      guides(linetype = guide_legend(title = "Method"),
             colour = guide_legend(title = "Guides"), 
             fill = guide_legend(title = "Guides"),
             shape = guide_legend(title = "Guides")) + theme(legend.position = "bottom",
                                                             legend.title = element_blank())
  }
  g <- g + theme(legend.box = "vertical")
  return(g)
}

#'---------------------------------------------
# Function to create 2D and 3D NMDS plots
#'---------------------------------------------

plotnmds <- function(nmds.object, ndim=2, size.pts, input.data, plot.factor, plot.ellipse = TRUE){
  
  NMDS1 <- nmds.object$points[,1]
  NMDS2 <- nmds.object$points[,2]
  if(ndim==3) NMDS3 <- nmds.object$points[,3]
  
  nmds.plot<-cbind(input.data, NMDS1, NMDS2)
  if(ndim==3) nmds.plot<-cbind(input.data, NMDS1, NMDS2, NMDS3)
  
  plot.factor <- deparse(substitute(plot.factor))
  
  plot.lvls <- levels(nmds.plot[,names(nmds.plot)==plot.factor])
  nmds.plot$colplot <- ifelse(nmds.plot[,names(nmds.plot)==plot.factor]=="bank", 
                              "#ffbf00", "#009acd")
  
  if(ndim==2){
    
    if(plot.ellipse==T){
      ggplot(nmds.plot, 
             aes(NMDS1, NMDS2, size = size.pts,
                 color=nmds.plot[,names(nmds.plot)==plot.factor]))+
        geom_point(position=position_jitter(.1), shape=16)+##separates overlapping points
        stat_ellipse(type='t',size =1)+ ##draws 95% confidence interval ellipses
        theme_minimal()+
        scale_colour_manual(values = c("#ffbf00", "#009acd"))+
        theme(legend.title=element_blank())
      
    }else{
      
      ggplot(nmds.plot, 
             aes(NMDS1, NMDS2, size = size.pts,
                 color=nmds.plot[,names(nmds.plot)==plot.factor]))+
        geom_point(position=position_jitter(.1), shape=16)+##separates overlapping points
        theme_minimal()+
        scale_colour_manual(values = c("#ffbf00", "#009acd"))+
        theme(legend.title=element_blank())} # End If plot.ellipse
    
  }else{
    
    plot_ly(x=NMDS1, y=NMDS2, z=NMDS3, 
            data = nmds.plot,
            type="scatter3d",
            group_by = factor(nmds.plot[,names(nmds.plot)==plot.factor]),
            mode="markers", 
            marker = list(
              color = nmds.plot$colplot 
            ))
  }
} # End function

#'---------------------------------------------
# Function to perform quick model diagnostics
#'---------------------------------------------

quick.validate <- function(gam.model, random.quantiles = TRUE){
  
  message(paste0("Overdispersion: ", overd(gam.model), " / " , ovd(gam.model)))
  print(summary(gam.model)) # Model summary
  
  par(mfrow = c(2,2))
  mgcv::gam.check(gam.model)
  abline(0, 1, col = "orange")
  
  if(random.quantiles){
    dsm::rqgam.check(gam.model)
    abline(0, 1, col = "orange")}
  
  # QQplot with confidence interval
  
  test.fit <- mgcViz::getViz(gam.model)
  mgcViz::qq(test.fit, rep = 50, showReps = FALSE, CI = "normal", a.qqpoi = list("shape" = 19),
             a.replin = list("alpha" = 0.2))
  
}

#'---------------------------------------------
# Function to plot semivariogram (check for spatial autocorrelation)
#'---------------------------------------------

plot.variogram <- function(model, dat){
  
  E <- resid(model, type = "pearson")
  mydata <- data.frame(E, dat$x, dat$y)
  names(mydata) <- c("res", "x", "y")
  sp::coordinates(mydata) <- c("x", "y")
  
  par(mfrow = c(2,2))
  purrr::map(c(1000, 2500, 5000, 10000),
             ~gstat::variogram(res ~ 1, mydata, cutoff = .)) %>% 
    purrr::map(., ~plot(x = .$dist, 
                        y = .$gamma,
                        xlab = "Distance (km)", 
                        ylab = "Semi-variogram", 
                        xlim = c(.$dist[2], max(.$dist)),
                        pch = 16,
                        cex = 2 *.$np / max(.$np)))}

#'---------------------------------------------
# Function to create plot for tensor product bk x depth
#'---------------------------------------------

plot.tensor <- function(gam.model){
  
  if(sum(grepl(pattern = "te\\(depth,bk\\)", x = names(gam.model$coefficients)))==0) stop("No tensor product found")
  
  # Convert to gamViz object
  
  b <- mgcViz::getViz(gam.model)
  
  # Determine which term is te(bk, depth)
  
  te.pos <- as.character(gam.model$formula)[3]
  te.pos <- as.list(strsplit(te.pos, '\\+')[[1]]) 
  te.pos <- which(purrr::map_lgl(.x = te.pos, .f = ~grepl(pattern = "depth, bk", x = .x)))
  
  # Create plot
  
  plot(sm(b, te.pos)) + l_fitRaster() + l_fitContour() 
  
}

#'---------------------------------------------
# Function to compile results from full subsets
#'--------------------------------------------- 

FSS.results <- function(FSS.list, plot.models, AIC = FALSE, verbose = TRUE){
  
  fss.obj <- list()
  
  # Failed and successful models
  
  fss.obj$failed.models <- FSS.list$failed.models
  fss.obj$success.models <- FSS.list$success.models
  
  if(verbose) message(paste0('Number of failed models: ', length(fss.obj$failed.models), ' (', 100 * round(length(fss.obj$failed.models)/(length(fss.obj$failed.models)+length(fss.obj$success.models))), '%)'))
  
  if(verbose) message(paste0('Number of successful models: ', length(fss.obj$success.models), ' (', 100 * round(length(fss.obj$success.models)/(length(fss.obj$failed.models)+length(fss.obj$success.models))), '%)'))
  
  # Comparative table
  
  fss.obj$model.table <- FSS.list$mod.data.out
  
  if(AIC)  fss.obj$model.table <-   fss.obj$model.table[order(fss.obj$model.table$AIC),]
  if(!AIC)  fss.obj$model.table <-   fss.obj$model.table[order(fss.obj$model.table$AICc),]
  
  # Variable importance
  
  fss.obj$var.importance <- FSS.list$variable.importance$aic$variable.weights.raw
  names(fss.obj$var.importance) <- gsub(pattern = "depth", replacement = "te(depth, bk)", x = names(fss.obj$var.importance))
  fss.obj$var.importance <- sort(fss.obj$var.importance[!names(fss.obj$var.importance)=="bk"], decreasing = TRUE)

  # All models within 2 AICc units
  
  if(AIC) fss.obj$top.models <- fss.obj$model.table[which(fss.obj$model.table$delta.AIC<2),]
  if(!AIC) fss.obj$top.models <- fss.obj$model.table[which(fss.obj$model.table$delta.AICc<2),]
  
  # Model with minimum AICc
  
  if(AIC)   fss.obj$best.model <- list(summary = fss.obj$model.table[which(fss.obj$model.table$AIC==min(fss.obj$model.table$AIC, na.rm = TRUE)),])
  if(!AIC)   fss.obj$best.model <- list(summary = fss.obj$model.table[which(fss.obj$model.table$AICc==min(fss.obj$model.table$AICc, na.rm = TRUE)),])

  fss.obj$best.model$obj <- fss.obj$success.models[[as.character(fss.obj$best.model$summary$modname)]]
  
  if(verbose) message(paste0('Top-ranking model: ', fss.obj$best.model$summary$modname))
  if(verbose) message(paste0('Most important variable: ', names(fss.obj$var.importance)[1]))
  
  if(plot.models){
    
    par(oma = c(1,1,4,1))
    
    for(r in 1:nrow(fss.obj$top.models)){
      
      best.model.name <- as.character(fss.obj$top.models$modname[r])
      best.model <- FSS.list$success.models[[best.model.name]]
      
      plot(best.model, all.terms = TRUE, pages = 1, residuals = FALSE, scale = 0, shade = TRUE)
      mtext(side = 3, text = resp.vars[i], outer = TRUE)
    }
  }
  
  return(fss.obj)
  
}
  
#'---------------------------------------------
# Function to update list of model sets
#'--------------------------------------------- 

update.formulae <- function(model.set, k = NULL){
  
  # Identify all models with tensor products
  
  mod.depth.bk <- names(model.set$mod.formula)[stringr::str_detect(string = names(model.set$mod.formula), pattern = "depth.te.bk", negate = FALSE)] 
  
  # Identify all other models
  
  mod.others <- names(model.set$mod.formula)[stringr::str_detect(string = names(model.set$mod.formula), pattern = "depth.te.bk", negate = TRUE)]
  
  # Only retain those without single terms for depth and bk
  
  mod.others <- mod.others[stringr::str_detect(string = mod.others, 
                                               pattern = "depth", 
                                               negate = TRUE)] # All models with depth terms
  
  mod.others <- mod.others[stringr::str_detect(string = mod.others, 
                                               pattern = "bk", 
                                               negate = TRUE)] # All models with bk terms
  
  # Put it all back together
  
  allmod <- c(mod.others, mod.depth.bk)
  
  model.set$mod.formula <- model.set$mod.formula[which(names(model.set$mod.formula)%in%allmod)]
  model.set$n.mods <- length(allmod)
  
  # Change basis dimension of tensor product
  
  # model.set$mod.formula <- purrr::map(.x = model.set$mod.formula, 
  #                                     .f = ~Reduce(paste, deparse(.)) %>%
  #                                       gsub(pattern = "bk, k = 5,", 
  #                                            replacement = paste0("bk, k = c(", k.tensor, ",", k.tensor, "), "), x = .) %>% 
  #                                       as.formula(.))
  
  if(is.null(k)){
    
    model.set$mod.formula <- purrr::map(.x = model.set$mod.formula, 
                                        .f = ~Reduce(paste, deparse(.)) %>%
                                          gsub(pattern = "k = 5,", 
                                               replacement = paste0(" "), x = .) %>% 
                                          as.formula(.)) 
    
  }
  
  return(model.set)
  
}
  
#'---------------------------------------------
# Wrapper around generate.model.set
#'--------------------------------------------- 

list.model.sets <- function(response.var, 
                            dat, 
                            corr.cutoff = 0.7, 
                            max.predictors = 3, 
                            verbose = TRUE,
                            k = NULL){
  
  # Prepare the data
  
  if(verbose) message('Preparing data ...')
  
  use.dat <- dat %>% 
    dplyr::select(response.var, covariates, duration, grid) %>% 
    dplyr::rename(response = response.var)
  
  # Initial model
  
  if(verbose) message('Building null model ...')
  
  # Build model formula

  if(is.null(k)){
    
    null.model <- mgcv::gam(formula = as.formula(paste0("response ~ te(depth, bk) + s(grid, bs = 're') + offset(log(duration)) + ", paste0('s(', covariates[!covariates%in%c('depth', 'bk')], ') ', collapse = '+ '), collapse = "+")), data = use.dat, REML = TRUE, family = poisson(link = "log"))
    
  }else{
    
    null.model <- mgcv::gam(formula = as.formula(paste0("response ~ te(depth, bk) + s(grid, bs = 're') + offset(log(duration)) + ", paste0('s(', covariates[!covariates%in%c('depth', 'bk')], ', k = ', k,' ) ', collapse = '+ '), collapse = "+")), data = use.dat, REML = TRUE, family = poisson(link = "log"))
    
  }
  
#     null.model <- mgcv::gam(formula = as.formula(paste0("response ~ te(depth, bk, k = c(", k.tensor, ",", k.tensor, ")) + s(grid, bs = 're') + offset(log(duration)) + ", paste0('s(', covariates[!covariates%in%c('depth', 'bk')], ', k = 5) ', collapse = '+ '), collapse = "+")), data = use.dat, REML = TRUE, family = poisson(link = "log"))


  # Tensor product of depth and bk needs to be specified as smooth.smooth.interactions
  # in order to have it included in variable importance calculations
  # This in turn requires pred.vars.cont to contain all covariates
  # But has the downside that depth and bk are included as single terms too
  # Use custom function update.formulae to correct model formulae accordingly
  
  if(verbose) message('Generating model set ...')
  
  if(is.null(k)){
    
    model.set <- generate.model.set(use.dat = use.dat,
                                    test.fit = null.model,
                                    cov.cutoff = corr.cutoff,
                                    non.linear.correlations = FALSE,
                                    pred.vars.cont = covariates,
                                    # pred.vars.cont = covariates[!covariates%in%c('depth', 'bk')],
                                    pred.vars.fact = NA, 
                                    smooth.smooth.interactions = c("depth", "bk"),
                                    max.predictors = max.predictors,
                                    bs.arg = "'tp'",
                                    null.terms = "s(grid, bs = 're') + offset(log(duration))") %>% 
      update.formulae(model.set = ., k = k)
    
  }else{
  
  model.set <- generate.model.set(use.dat = use.dat,
                                  test.fit = null.model,
                                  cov.cutoff = corr.cutoff,
                                  non.linear.correlations = FALSE,
                                  pred.vars.cont = covariates,
                                  # pred.vars.cont = covariates[!covariates%in%c('depth', 'bk')],
                                  pred.vars.fact = NA, 
                                  smooth.smooth.interactions = c("depth", "bk"),
                                  max.predictors = max.predictors,
                                  k = k,
                                  bs.arg = "'tp'",
                                  null.terms = "s(grid, bs = 're') + offset(log(duration))") %>% 
    update.formulae(model.set = ., k = k)
  
  }
  
  if(verbose) message('Done!')
  return(model.set)
  
}

#'---------------------------------------------
# Function to remove raster values outside covariate ranges
#'--------------------------------------------- 

extrap.filter <- function(input.raster.stack, 
                          covariate.rge, 
                          dat, 
                          method = "convex", 
                          verbose = TRUE){
  
  if(verbose) message('Converting stack to data.frame ...')
  
  # Store projection system
  
  proj.syst <- sp::proj4string(input.raster.stack)
  
  # Convert raster stack to df
  
  temp.df <- raster::as.data.frame(x = input.raster.stack, xy = TRUE, na.rm = TRUE)
  
  # Discard data outside to univariate ranges
  
  if(verbose) message('Discarding data outside univariate range ...')
  
  for (g in 1:ncol(covariate.rge)){
    extracted.var <- temp.df[,which(names(temp.df)==names(covariate.rge)[g])]
    temp.df <- temp.df[which(extracted.var > covariate.rge[,g][1] & extracted.var < covariate.rge[,g][2]),]
  }
  
  # Discard data outside convex/concave hull (depth x bk)
  
  df.pts <- temp.df %>% 
    dplyr::select(depth, bk) %>% 
    sp::SpatialPoints(coords = .)
  
  locs.depth.bk <- dat %>% 
    dplyr::select(depth, bk)
  
  if(method == "convex"){
    
    if(verbose) message('Computing convex hull ...')
    
    hpts <- grDevices::chull(locs.depth.bk)
    convex.hull <- locs.depth.bk[c(hpts, hpts[1]), ]
    convex.hull <- sp::SpatialPolygons(list(Polygons(list(Polygon(convex.hull)), ID = 1))) 
    
  }else if (method == "concave"){
    
    if(verbose) message('Computing concave hull ...')
    
    locs.depth.bk <- locs.depth.bk %>% sp::SpatialPoints(coords = .)
    concave.hull <- concaveman::concaveman(locs.depth.bk)
    concave.hull <- broom::tidy(concave.hull) %>% 
      dplyr::select(long, lat) %>% 
      sp::SpatialPoints(coords = .)
    concave.hull <- sp::SpatialPolygons(list(Polygons(list(Polygon(concave.hull)), ID = 1)))

  }else if (method == "alpha"){
    
    if(verbose) message("Computing alpha hull ...")
    
    alpha.hull <- ConR::EOO.computing(XY = locs.depth.bk[, c(2,1)],
                                 exclude.area = FALSE,
                                 export_shp = TRUE,
                                 method.range = "alpha.hull",
                                 alpha = 5,
                                 buff.alpha = 0.001,
                                 write_results = FALSE,
                                 show_progress = FALSE)
    alpha.hull <- alpha.hull[[2]]
  }
    

  if(method == "convex"){
    
    temp.df$inside <- sp::over(x = df.pts, y = convex.hull)
  
  }else if (method =="concave"){
      
    temp.df$inside <- sp::over(x = df.pts, y = concave.hull)
    
  }else if(method=="alpha"){
      
    temp.df$inside <- sp::over(x = df.pts, y = alpha.hull)
    
    }
  
    # Filter the data
    
    temp.df <- temp.df %>% 
      dplyr::filter(inside == 1) %>% 
      dplyr::select(., -inside) %>% 
      tibble::as_tibble(.)
  
    if(verbose) message('Preparing rasters ...')
  
  # Create rasters from filtered data
  
  col.indices <- which(!names(temp.df)%in%c("x", "y"))
  
  r.filt <- purrr::map(.x = col.indices, 
                       .f = ~temp.df %>% 
                         dplyr::select(., x, y, .x) %>% 
                         raster::rasterFromXYZ(xyz = ., crs = proj.syst)) %>% 
    purrr::set_names(., names(temp.df[col.indices]))
  
  # Stack all rasters
  
  r.filt <- raster::stack(r.filt)

  if(verbose) message('Done!')
  return(r.filt)
}

#'---------------------------------------------
# Modified version of FSSgam::extract.mod.dat that includes AIC
#'--------------------------------------------- 

extract.model.dat <- function (mod.fit, r2.type. = r2.type) {
  mod.dat = list(AICc = NA, AIC = NA, BIC = NA, r2.vals = NA, r2.vals.unique = NA, 
                 edf = NA, edf.less.1 = NA)
  if (class(mod.fit)[[1]] != "try-error") {
    # mod.dat$AIC = MuMIn::QAIC(mod.fit, chat = mod.fit$deviance/mod.fit$df.residual)
    mod.dat$AIC = AIC(mod.fit)
    mod.dat$AICc = AICc(mod.fit)
    mod.dat$BIC = BIC(mod.fit)
    tempOut = NA
    if (class(mod.fit)[1] == "gam" & r2.type. == "dev") {
      tempOut = summary(mod.fit)$dev.expl
    }
    if (class(mod.fit)[1] == "gam" & r2.type. == "r2") {
      tempOut = summary(mod.fit)$r.sq
    }
    if (class(mod.fit)[1] == "gam" & r2.type. == "r2.lm.est") {
      tempOut = summary(lm(mod.fit$y ~ predict(mod.fit)))$r.sq
    }
    if (class(mod.fit)[[1]] == "gamm4" & r2.type. == "dev") {
      tempOut = summary(mod.fit$gam)$dev.expl
      if (length(tempOut) == 0) {
        tempOut = NA
      }
    }
    if (class(mod.fit)[[1]] == "gamm4" & r2.type. == "r2") {
      tempOut = summary(mod.fit$gam)$r.sq
    }
    if (class(mod.fit)[[1]] == "gamm" & r2.type. == "r2") {
      tempOut = summary(mod.fit$gam)$r.sq
    }
    if (class(mod.fit)[[1]] == "gamm4" & r2.type. == "r2.lm.est") {
      tempOut = summary(lm(attributes(mod.fit$mer)$frame$y ~ 
                             predict(mod.fit[[1]], re.form = NA, type = "response")))$r.sq
    }
    if (is.null(tempOut)) {
      tempOut = NA
    }
    mod.dat$r2.vals = round(tempOut, 5)
    if (class(mod.fit)[1] == "gam") {
      edf.m = summary(mod.fit)$edf
      p.coeff.m = summary(mod.fit)$p.coeff
    }
    else {
      edf.m = mod.fit$gam$edf
      p.coeff.m = mod.fit$gam$p.coeff
    }
    edf.m[which(edf.m < 1)] = 1
    mod.dat$edf = round(sum(c(edf.m, length(p.coeff.m))), 
                        2)
    if (class(mod.fit)[1] == "gam") {
      edf.m = summary(mod.fit)$edf
    }
    else {
      edf.m = mod.fit$gam$edf
    }
    mod.dat$edf.less.1 = length(which(edf.m < 0.25))
  }
  return(mod.dat)
}

#'---------------------------------------------
# Modified version of FSSgam::fit.mod.set that includes AIC
#'--------------------------------------------- 

fit.model.set <- function (model.set.list, max.models = 200, save.model.fits = T, parallel = F, n.cores = 4, r2.type = "r2.lm.est", report.unique.r2 = F) 
  {
    use.datModSet <- model.set.list$used.data
    n.mods = length(model.set.list$mod.formula)
    mod.formula = model.set.list$mod.formula
    test.fit = model.set.list$test.fit
    included.vars = model.set.list$included.vars
    if (n.mods > max.models) {
      warning(paste("You have ", n.mods, " models. Individual models fits will not be saved.\n        If you want to save all the model fits all of these you need to\n        increase 'max.models' from ", 
                    max.models, ".", sep = ""))
      save.model.fits = F
    }
    require(MuMIn)
    if (save.model.fits == T) {
      pb <- txtProgressBar(max = length(mod.formula), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      if (parallel == T) {
        require(doSNOW)
        cl = makeCluster(n.cores)
        registerDoSNOW(cl)
        opts <- list(progress = progress)
        out.dat <- foreach(l = 1:length(mod.formula), .packages = c("mgcv", 
                                                                    "gamm4", "MuMIn", "FSSgam"), .errorhandling = "pass", 
                           .options.snow = opts) %dopar% {
                             fit.mod.l(mod.formula[[l]], test.fit. = test.fit, 
                                       use.dat = use.datModSet)
                           }
        close(pb)
        stopCluster(cl)
        registerDoSEQ()
      }
      else {
        out.dat <- vector("list", length(mod.formula))
        for (l in 1:length(mod.formula)) {
          mod.l = fit.mod.l(mod.formula[[l]], test.fit. = test.fit, 
                            use.dat = use.datModSet)
          out.dat[[l]] = mod.l
          setTxtProgressBar(pb, l)
        }
      }
      close(pb)
      names(out.dat) = names(mod.formula)[1:n.mods]
      model.success = lapply(lapply(out.dat, FUN = class), 
                             FUN = function(x) {
                               length(grep("gam", x)) > 0
                             })
      failed.models = mod.formula[which(model.success == F)]
      success.models = out.dat[which(model.success == T)]
      if (length(success.models) == 0) {
        stop("None of your models fitted successfully. Please check your input objects.")
      }
      var.inclusions = build.inclusion.mat(included.vars = included.vars, 
                                           formula.list = success.models)
      mod.data.out = data.frame(modname = names(success.models))
      mod.data.out$formula = unlist(lapply(success.models, 
                                           FUN = function(x) {
                                             as.character(formula(x)[3])
                                           }))
      mod.data.out = cbind(mod.data.out, do.call("rbind", lapply(success.models, 
                                                                 FUN = function(x) {
                                                                   unlist(extract.model.dat(x, r2.type. = r2.type))
                                                                 })))
    }
    else {
      var.inclusions = build.inclusion.mat(included.vars = included.vars, 
                                           formula.list = mod.formula)
      mod.data.out = data.frame(modname = names(mod.formula))
      mod.data.out$formula = unlist(lapply(mod.formula, FUN = function(x) {
        as.character(formula(x))[2]
      }))
      pb <- txtProgressBar(max = length(mod.formula), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      if (parallel == T) {
        require(doSNOW)
        cl = makeCluster(n.cores)
        registerDoSNOW(cl)
        opts <- list(progress = progress)
        mod.dat <<- foreach(l = 1:length(mod.formula), .packages = c("mgcv", 
                                                                     "gamm4", "MuMIn", "FSSgam"), .options.snow = opts) %dopar% 
          {
            unlist(extract.model.dat(fit.mod.l(mod.formula[[l]], 
                                               test.fit. = test.fit, use.dat = use.datModSet), 
                                     r2.type. = r2.type))
          }
        close(pb)
        stopCluster(cl)
        registerDoSEQ()
      }
      else {
        mod.dat = vector("list", length(mod.formula))
        for (l in 1:length(mod.formula)) {
          mod.l = fit.mod.l(mod.formula[[l]], test.fit. = test.fit, 
                            use.dat = use.datModSet)
          out = unlist(extract.model.dat(mod.l, r2.type. = r2.type))
          mod.dat[[l]] = out
          setTxtProgressBar(pb, l)
        }
      }
      close(pb)
      names(mod.dat) = names(mod.formula[1:n.mods])
      mod.data.out = cbind(mod.data.out, do.call("rbind", mod.dat))
      failed.models = mod.formula[which(is.na(mod.data.out$AICc) == 
                                          T)]
      success.models = mod.formula[which(is.na(mod.data.out$AICc) == 
                                           F)]
    }
    mod.data.out$delta.AIC = round(mod.data.out$AIC - min(mod.data.out$AIC, 
                                                          na.rm = T), 3)
    
    mod.data.out$delta.AICc = round(mod.data.out$AICc - min(mod.data.out$AICc, 
                                                            na.rm = T), 3)
    mod.data.out$delta.BIC = round(mod.data.out$BIC - min(mod.data.out$BIC, 
                                                          na.rm = T), 3)
    mod.data.out$wi.AIC = round(wi(mod.data.out$AIC), 3)
    mod.data.out$wi.AICc = round(wi(mod.data.out$AICc), 3)
    mod.data.out$wi.BIC = round(wi(mod.data.out$BIC), 3)
    if (report.unique.r2 == T) {
      null.r2 = mod.data.out$r2.vals[which(mod.data.out$modname == 
                                             "null")]
      mod.data.out$r2.vals.unique = mod.data.out$r2.vals - 
        null.r2
    }
    mod.data.out = cbind(mod.data.out, var.inclusions)
    min.mods = min(colSums(mod.data.out[, included.vars]))
    var.weights = unlist(lapply(included.vars, FUN = function(x) {
      sum(sort(mod.data.out$wi.AICc[which(mod.data.out[, x] == 
                                            1)], decreasing = T)[1:min.mods])
    }))
    var.weights = unlist(lapply(included.vars, FUN = function(x) {
      sum(sort(mod.data.out$wi.AIC[which(mod.data.out[, x] == 
                                           1)], decreasing = T)[1:min.mods])
    }))
    names(var.weights) = included.vars
    variable.weights.raw = var.weights
    aic.var.weights = list(variable.weights.raw = variable.weights.raw)
    var.weights = unlist(lapply(included.vars, FUN = function(x) {
      sum(sort(mod.data.out$wi.BIC[which(mod.data.out[, x] == 
                                           1)], decreasing = T)[1:min.mods])
    }))
    names(var.weights) = included.vars
    variable.weights.raw = var.weights
    bic.var.weights = list(variable.weights.raw = variable.weights.raw)
    return(list(mod.data.out = mod.data.out, failed.models = failed.models, 
                success.models = success.models, variable.importance = list(aic = aic.var.weights, 
                                                                            bic = bic.var.weights)))
  }

#'---------------------------------------------
# Function to generate model predictions, taking parameter uncertainty into account
#'--------------------------------------------- 

get.predictions <- function(model, 
                             proj.syst,
                             pred.pts = NULL, 
                             boot.dat,
                             method = "convex",
                             sample.coef = FALSE,
                             verbose = FALSE){
  
  # Define univariate range of covariates
  
  env.rg <- data.frame(apply(boot.dat[, covariates[2:8]], MARGIN = 2, FUN = function(x) range(x)))
  
  # Filter out areas outside range
  
  stack.df <- suppressMessages(extrap.filter(input.raster.stack = pred.pts, 
                            covariate.rge = env.rg, 
                            dat = boot.dat,
                            method = method,
                            verbose = FALSE))
  
  r.boundaries <- c(apply(coordinates(pred.pts), 2, min)[2], apply(coordinates(pred.pts), 2, max)[2])
  stack.df <- raster::as.data.frame(x = stack.df, xy = TRUE, na.rm = TRUE)
    

  whichgrid <- as.numeric(which(sapply(grid.centroids$y, 
                                       FUN = function(x) {x > r.boundaries[1] & x < r.boundaries[2]})))

  # Random effect term
  
  stack.df$grid <- factor(as.character(whichgrid), levels = levels(osdata$grid))
  
  # Set offset to mean duration
  
  stack.df$duration <- mean(boot.dat$duration)
  
  # Get model predictions
  
  # Because GAMs have an underlying parametric representation,
  # it is possible to obtain a 'prediction matrix', Xp, say, which maps the model
  # parameters, to the predictions of the linear predictor,
  # One approach is to parameterize the spline in terms of its values at the knots.
  # When type="lpmatrix" then a matrix is returned which yields the values of the linear predictor (minus any offset) when post-multiplied by the parameter vector (in this case se.fit is ignored).
  
  if(verbose) message('Generating model predictions ...')
  
  preds <- mgcv::predict.gam(object = model,
                             newdata = stack.df,
                             exclude = "s(grid)",
                             type = "lpmatrix")
  
  
  # Other way of doing it
  # https://stackoverflow.com/questions/51634953/gam-with-gp-smoother-predict-at-new-locations/51635487#51635487
  
  # preds <- mgcv::predict.gam(object = model,
  #                            newdata = stack.df,
  #                            exclude = "s(grid)",
  #                            type = "terms")
  # 
  # ## linear predictor without random effect - simplify to vector
  # lp <- rowSums(preds) + attr(preds, "constant")
  # lp <- lp + log(mean(osdata$duration))
  # lp <- family(model)$linkinv(lp)
  
  # Model coefficients and variance-covariance Matrix
  
  beta <- coef(object = model)
  Vb <- vcov(object = model)
  
  # Draw from posterior distribution of parameters, if required 
  
  if(sample.coef){
    if(verbose) message('Sampling from posterior ...')
    beta <- MASS::mvrnorm(n = 1, mu = beta, Sigma = Vb)}
  
  # Predictions on scale of linear predictor
  
  preds.lp <- preds %*% beta
  
  # Add offset
  
  preds.lp <- preds.lp + log(mean(boot.dat$duration))
  
  
  # Use the inverse link function to obtain predictions on the response scale 
  
  stack.df$preds <- family(model)$linkinv(preds.lp)
  
  # Create raster
  
  stack.df <- stack.df %>% dplyr::select(x, y, preds)
  
  if(class(pred.pts) == "RasterStack"){
    
    if(verbose) message('Mapping predictions to raster ...')
    
    final.r <- stack.df %>%
      raster::rasterFromXYZ(., crs = proj.syst)  
    
    message('')
    return(final.r)
    
  }else{
    
    return(stack.df)
    
  }
  
}

#'---------------------------------------------
# Summary statistics for raster objects
#'---------------------------------------------

sumrast <- function(input.raster){
  summary(raster::getValues(input.raster), na.rm = TRUE)
}

#'---------------------------------------------
# Function to plot predictions over sampling grids
#'---------------------------------------------

plot.all <- function(log.response.sr = FALSE,
                     log.response.max_n = FALSE,
                     log.cv.sr = FALSE,
                     log.cv.max_n = FALSE,
                     cv.over.one.sr = FALSE,
                     cv.over.one.max_n = FALSE,
                     z.int.sr = 1,
                     z.int.max_n = 5,
                     z.int.sr.cv = 0.2,
                     z.int.max_n.cv = 0.2)
  {
  
  
  # Grid 1 ----
  
  plot.grid(input.raster = boot.sr$m$g1, 
            response.variable = "sr", 
            log.response = log.response.sr,
            grid.number = 1, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = sr.fivenum.boot$g1$Max.,
            z.int = z.int.sr, 
            skew = 0,
            max.pt.size = 3,
            col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
            save.to.pdf = TRUE, 
            plot.figure = TRUE,
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_boot")
  
  plot.grid(input.raster = boot.max_n$m$g1, 
            response.variable = "max_n", 
            log.response = log.response.max_n,
            grid.number = 1, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = max_n.fivenum.boot$g1$Max.,
            z.int = z.int.max_n,
            skew = 0,
            col.palette = rev(brewer.spectral(100)[10:100]),
            save.to.pdf = TRUE, 
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_boot")
  
  # ................................
  
  plot.grid(input.raster = boot.sr$rcv$g1, 
            response.variable = "sr", 
            log.response = log.cv.sr,
            cv.over.one = cv.over.one.sr,
            grid.number = 1, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = 1,
            z.int = z.int.sr.cv, 
            skew = 0,
            max.pt.size = 3,
            col.palette = pals::viridis(100)[10:100],
            save.to.pdf = TRUE, 
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_cv")
  
  plot.grid(input.raster = boot.max_n$rcv$g1, 
            response.variable = "max_n", 
            log.response = log.cv.max_n,
            cv.over.one = cv.over.one.max_n,
            grid.number = 1, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = 1,
            z.int = z.int.max_n.cv,
            skew = 0,
            over.one.col = "darkorange",
            col.palette = pals::viridis(100)[10:100],
            save.to.pdf = TRUE, 
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_cv")
  
  # plot.grid(input.raster = boot.sr$median$g1, 
  #           response.variable = "sr", 
  #           log.response = log.response.sr,
  #           grid.number = 1, 
  #           add.banks = FALSE, 
  #           add.bruvs = FALSE,
  #           zmin = 0,
  #           zmax = sr.fivenum.boot$g1$Max.,
  #           z.int = z.int.sr, 
  #           skew = 0,
  #           max.pt.size = 3,
  #           col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
  #           save.to.pdf = TRUE, 
  #           plot.figure = FALSE,
  #           add.arrow = FALSE, 
  #           add.scale = FALSE,
  #           plot.legend = TRUE,
  #           file.name = "preds_boot_colourscale")
  # 
  # plot.grid(input.raster = boot.max_n$median$g1, 
  #           response.variable = "max_n", 
  #           log.response = log.response.max_n,
  #           grid.number = 1, 
  #           add.banks = FALSE, 
  #           add.bruvs = FALSE,
  #           zmin = 0,
  #           zmax = max_n.fivenum.boot$g1$Max.,
  #           z.int = z.int.max_n, 
  #           skew = 0,
  #           max.pt.size = 3,
  #           col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
  #           save.to.pdf = TRUE, 
  #           plot.figure = FALSE,
  #           add.arrow = FALSE, 
  #           add.scale = FALSE,
  #           plot.legend = TRUE,
  #           file.name = "preds_boot_colourscale")
  # 
  # plot.grid(input.raster = boot.sr$rcv$g1, 
  #           response.variable = "sr", 
  #           log.response = log.cv.sr,
  #           cv.over.one = cv.over.one.sr,
  #           grid.number = 1, 
  #           add.banks = TRUE, 
  #           add.bruvs = FALSE,
  #           zmin = 0,
  #           zmax = 1,
  #           z.int = z.int.sr.cv, 
  #           skew = 0,
  #           max.pt.size = 3,
  #           col.palette = pals::viridis(100)[10:100],
  #           save.to.pdf = TRUE, 
  #           plot.figure = FALSE,
  #           add.arrow = FALSE, 
  #           add.scale = FALSE,
  #           plot.legend = TRUE,
  #           file.name = "preds_cv_colourscale")
  # 
  # plot.grid(input.raster = boot.max_n$rcv$g1, 
  #           response.variable = "max_n", 
  #           log.response = log.cv.max_n,
  #           cv.over.one = cv.over.one.max_n,
  #           grid.number = 1, 
  #           add.banks = TRUE, 
  #           add.bruvs = FALSE,
  #           zmin = 0,
  #           zmax = 1,
  #           z.int = z.int.max_n.cv,
  #           skew = 0,
  #           over.one.col = "darkorange",
  #           col.palette = pals::viridis(100)[10:100],
  #           save.to.pdf = TRUE, 
  #           plot.figure = FALSE,
  #           add.arrow = FALSE, 
  #           add.scale = FALSE,
  #           plot.legend = TRUE,
  #           file.name = "preds_cv_colourscale")
  
  
  
  # Grid 2 ----
  
  plot.grid(input.raster = boot.sr$m$g2, 
            response.variable = "sr", 
            log.response = log.response.sr,
            grid.number = 2, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = sr.fivenum.boot$g2$Max.,
            z.int = z.int.sr, 
            skew = 0,
            max.pt.size = 3,
            col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
            save.to.pdf = TRUE, 
            plot.figure = TRUE,
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_boot")
  
  plot.grid(input.raster = boot.max_n$m$g2, 
            response.variable = "max_n", 
            log.response = log.response.max_n,
            grid.number = 2, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = max_n.fivenum.boot$g2$Max.,
            z.int = z.int.max_n,
            skew = 0,
            col.palette = rev(brewer.spectral(100)[10:100]),
            save.to.pdf = TRUE, 
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_boot")
  
  # ................................
  
  plot.grid(input.raster = boot.sr$rcv$g2, 
            response.variable = "sr", 
            log.response = log.cv.sr,
            cv.over.one = cv.over.one.sr,
            grid.number = 2, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = 1,
            z.int = z.int.sr.cv, 
            skew = 0,
            max.pt.size = 3,
            col.palette = pals::viridis(100)[10:100],
            save.to.pdf = TRUE, 
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_cv")
  
  plot.grid(input.raster = boot.max_n$rcv$g2, 
            response.variable = "max_n", 
            log.response = log.cv.max_n,
            cv.over.one = cv.over.one.max_n,
            grid.number = 2, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = 1,
            z.int = z.int.max_n.cv,
            skew = 0,
            over.one.col = "darkorange",
            col.palette = pals::viridis(100)[10:100],
            save.to.pdf = TRUE, 
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_cv")
  
  
  
  
  
  
  # Grid 3 ----
  
  plot.grid(input.raster = boot.sr$m$g3, 
            response.variable = "sr", 
            log.response = log.response.sr,
            grid.number = 3, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = sr.fivenum.boot$g3$Max.,
            z.int = z.int.sr, 
            skew = 0,
            max.pt.size = 3,
            col.palette = rev(pals::brewer.spectral(100)[c(15:100, rep(100, 10))]),
            save.to.pdf = TRUE, 
            plot.figure = TRUE,
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_boot")
  
  plot.grid(input.raster = boot.max_n$m$g3, 
            response.variable = "max_n", 
            log.response = log.response.max_n,
            grid.number = 3, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = max_n.fivenum.boot$g3$Max.,
            z.int = z.int.max_n,
            skew = 0,
            col.palette = rev(brewer.spectral(100)[10:100]),
            save.to.pdf = TRUE, 
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_boot")
  
  # ................................
  
  plot.grid(input.raster = boot.sr$rcv$g3, 
            response.variable = "sr", 
            log.response = log.cv.sr,
            cv.over.one = cv.over.one.sr,
            grid.number = 3, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = 1,
            z.int = z.int.sr.cv, 
            skew = 0,
            max.pt.size = 3,
            col.palette = pals::viridis(100)[10:100],
            save.to.pdf = TRUE, 
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_cv")
  
  plot.grid(input.raster = boot.max_n$rcv$g3, 
            response.variable = "max_n", 
            log.response = log.cv.max_n,
            cv.over.one = cv.over.one.max_n,
            grid.number = 3, 
            add.banks = TRUE, 
            add.bruvs = FALSE,
            zmin = 0,
            zmax = 1,
            z.int = z.int.max_n.cv,
            skew = 0,
            over.one.col = "darkorange",
            col.palette = pals::viridis(100)[10:100],
            save.to.pdf = TRUE, 
            add.arrow = FALSE, 
            add.scale = FALSE,
            plot.legend = TRUE,
            file.name = "preds_cv")
  
}

plot.grid <- function(input.raster,
                      grid.number = NULL,
                      response.variable = "sr",
                      log.response = FALSE,
                      cv.over.one = FALSE,
                      over.one.col = 'grey',
                      add.banks = TRUE,
                      add.bruvs = TRUE,
                      zmin = NULL,
                      zmax = NULL,
                      z.int = 1,
                      skew = 0,
                      pts.outline = TRUE,
                      col.palette, 
                      no.pts.legend = 4, 
                      multiply = 1, 
                      min.pt.size = 1,
                      max.pt.size = 4,
                      add.arrow = TRUE,
                      add.scale = TRUE,
                      plot.axes = FALSE,
                      plot.figure = TRUE,
                      plot.legend = TRUE,
                      col.bruvs = "black",
                      save.to.pdf = FALSE,
                      file.name = NULL){

  
  if(save.to.pdf) {
    if(is.null(file.name)) stop("Must specify output file name")
    if(grid.number==1) pdf(file = paste0('output/', file.name, '_', response.variable, '_grid1.pdf'), width = 5, height = 6)
    if(grid.number==2) pdf(file = paste0('output/', file.name, '_', response.variable, '_grid2.pdf'), width = 6, height = 5)
    if(grid.number==3) pdf(file = paste0('output/', file.name, '_', response.variable, '_grid3.pdf'), width = 6, height = 4)
  }
  
  #'---------------------------------------------
  # Multiply is for CVs etc to change z scale
  #'--------------------------------------------- 
  
  input.raster <- raster::projectRaster(input.raster*multiply, crs = wgs84)
  input.raster[input.raster<0] <- NA

  #'---------------------------------------------
  # Modify Z scale if needed
  #'--------------------------------------------- 
  
  min.z <- zmin
  max.z <- zmax
  
  if(is.null(zmin)) min.z <- min(raster::getValues(input.raster), na.rm = TRUE)
  if(is.null(zmax)) max.z <- max(raster::getValues(input.raster), na.rm = TRUE)
  
  #'---------------------------------------------
  # Axis labels
  #'--------------------------------------------- 
  
  xat <- pretty(range(coordinates(input.raster)[,1]))
  xlab <- paste0(xat, "E")
  
  yat <- pretty(range(coordinates(input.raster)[,2]))
  ylab <- paste0(abs(yat), "S")
  
  # if(quantile.breaks){
  #   
  #   brks <- raster::getValues(input.raster)
  #   brks <- brks[!is.na(brks)]
  #   brks <- quantile_breaks(xs = brks, n = no.quantiles)}
  
  #'---------------------------------------------
  # Plot raster
  #'--------------------------------------------- 
  
  grab_grob <- function(){
    grid.echo()
    grid.grab()
  }
  
  if(cv.over.one){
    
    if(log.response){
      over1.g <- log(input.raster)
      over1.g[over1.g<=0] <- NA 
    }else{
    over1.g <- input.raster
    over1.g[over1.g<=1] <- NA
    }
  }
  
  if(plot.figure){
    
  if(log.response){
    
    if(plot.axes){
      plot(log(input.raster), 
           col = col.palette, 
           legend = FALSE, 
           axes = FALSE)
    }else{
    plot(log(input.raster), 
         col = col.palette, 
         legend = FALSE, 
         axes = FALSE, 
         box = FALSE)}
    
  }else{
    if(plot.axes){
      plot(input.raster, 
           col = col.palette, 
           legend = FALSE, 
           zlim = c(min.z, max.z), 
           axes = FALSE)
    }else{
    plot(input.raster, 
         col = col.palette, 
         legend = FALSE, 
         zlim = c(min.z, max.z), 
         axes = FALSE, box = FALSE)}
  }
  }else{
    
    empty.raster <- input.raster
    values(empty.raster) <- NA
    plot(empty.raster, legend = FALSE, axes = FALSE, box = FALSE)
    
}
  
  #'---------------------------------------------
  # Adds axes
  #'--------------------------------------------- 
  
  if(plot.axes){
  axis(1, at = xat, labels = xlab)
  axis(2, at = yat, labels = ylab)}
  
  if(log.response){
    
    if(skew>0){

      ticks.seq <- seq(min.z, max.z, by = z.int)
      
      skew.threshold <- as.numeric(quantile(ticks.seq, 1-skew))
      n.skew <- (ticks.seq-1<=skew.threshold)
      if(n.skew==1) n.skew = n.skew + 1
      
      tickslabels <- pretty(c(ticks.seq[1], ticks.seq[n.skew]), n = max(round(z.int*skew), n = round(z.int*(1-skew))))
      
      tickslabels <- c(tickslabels, pretty(c(skew.threshold, ticks.seq[length(ticks.seq)]), n = max(length(ticks.seq)-length(tickslabels), length(tickslabels)-length(ticks.seq))))
      
      tickslabels <- tickslabels[!duplicated(tickslabels)] 
      
    }else{
      
      tickslabels <- seq(min.z, max.z, by = z.int)
      tickslabels <- c(seq(0, 1, by = 0.1), seq(1, 10, by = 1), seq(10, 100, 10))
      tickslabels <- tickslabels[!duplicated(tickslabels)] 
      
    }

  
    whereticksat <- log(tickslabels)
    whereticksat[1] <- log(0.001)
    
    tickslabels[!tickslabels%in%c(0, 0.1, 1, 10, 100)] <- NA
    tickslabels <- as.character(tickslabels)
    
  }else{
    
    whereticksat <- seq(min.z, max.z, by = z.int)
    tickslabels <- as.character(whereticksat) 
  }

  #'---------------------------------------------
  # Add legend
  #'--------------------------------------------- 
  if(plot.legend){
  if(log.response){
   
    plot(input.raster, 
         col = col.palette,
         legend.only = TRUE, 
         legend.shrink = 1, 
         legend.width = 1.5,
         bty="n",
         zlim = c(log(min.z+0.001), log(max.z)),
         axis.args = list(at = whereticksat,
                          labels = tickslabels))
     
  }else{
    
  plot(input.raster, 
         col = col.palette,
         legend.only = TRUE, 
         legend.shrink = 1, 
         legend.width = 1.5,
         bty="n",
         zlim = c(min.z, max.z),
         axis.args = list(at = whereticksat,
                          labels = tickslabels))
  }
}
  #'---------------------------------------------
  # Add banks
  #'--------------------------------------------- 
  
  if(add.banks){plot(sp::spTransform(banks[[grid.number]], CRSobj = wgs84), add = TRUE)}
  
  #'---------------------------------------------
  # Add bruvs
  #'--------------------------------------------- 
  
  if(add.bruvs){
    
    #'---------------------------------------------
    # Filter relevant data
    #'--------------------------------------------- 
    
    bruvs <- osdata %>% 
      dplyr::filter(grid == as.character(grid.number)) %>% 
    dplyr::select(., lon, lat, response.variable)
    
    bruvs.values <- bruvs %>% dplyr::pull(response.variable)
    
    #'---------------------------------------------
    # Parameters for plotting
    #'--------------------------------------------- 
    
    bruvs.min <- min(bruvs.values, na.rm = TRUE)
    bruvs.max <- max(bruvs.values,  na.rm = TRUE)
    
    bruvs.labels <- as.character(pretty(bruvs.min:bruvs.max, n = no.pts.legend))
  
    
    # if(log.response){
    #   
    #   bruvs.max <- log(max(plot.cex+1), base = 2)
    #   plot.cex <- log(plot.cex+1, base = 2)
    # 
    # }else{
    #   
    #   plot.cex <- plot.cex/cex.coef
    #   }
    # 
    # bruvs.legend <- pretty(bruvs.min:bruvs.max, n = no.pts.legend)
    
    #'---------------------------------------------
    # Convert to spatialpts
    #'--------------------------------------------- 
    
    lm.df <- data.frame(x = range(as.numeric(bruvs.labels)),
                        y = seq(min.pt.size, max.pt.size, length.out = 2))
    
    ndat <- data.frame(x = bruvs[, response.variable])
    names(ndat) <- "x"
    
    bruvs$pts.cex <- predict(object = lm(y~x, 
                                         data = lm.df),
                              newdata = ndat)
    
    bruvs.sp <- bruvs %>% 
      sp::SpatialPointsDataFrame(coords = cbind(.$lon, .$lat), data = ., proj4string = wgs84)
    
    #'---------------------------------------------
    # Plot points and legend
    #'--------------------------------------------- 

    points(bruvs.sp, pch = 16, cex = bruvs.sp@data$pts.cex, 
           col = alpha(col.bruvs, 0.25), lwd = 1)
    
    points(bruvs.sp, pch = 1, cex = bruvs.sp@data$pts.cex, col = "black", lwd = 1)
    
    legend("topright", 
           legend = bruvs.labels, 
           pch = 21, 
           pt.bg = alpha(col.bruvs, 0.25),
           pt.cex = seq(min.pt.size, max.pt.size, length.out = length(bruvs.labels)))

  }

  
  if(cv.over.one) plot(raster::projectRaster(over1.g, crs = wgs84), 
                       col = over.one.col, add = TRUE, bty="n", legend = FALSE, axes = FALSE)
 
  #'---------------------------------------------
  # Function to add north arrow
  #'---------------------------------------------
  
  # http://www.flutterbys.com.au/stats/tut/tut5.4.html
  
  north.arrow <-  function(x, y, h, lab = "N", lab.pos = "above") {
    polygon(c(x, x, x + h/2), c(y - (1.5*h), y, y - (1 + sqrt(3)/2) * h), col = "black", border = NA)
    polygon(c(x, x + h/2, x, x - h/2), c(y - (1.5*h), y - (1 + sqrt(3)/2) * h, y, y - (1 + sqrt(3)/2) * h), col = "black")
    if(lab.pos=="below") text(x, y-(2.5*h), lab, adj = c(0.5, 0), cex = 1)
    else text(x, y+(0.5*h), lab, adj = c(0.5, 0), cex = 1)
  }
  
  if(add.arrow){
    
    # prettymapr::addnortharrow(pos = "bottomright", scale = 0.5)
    if(grid.number==1) north.arrow(127.37, -12.075, 0.0085)
    if(grid.number==2) north.arrow(126.9, -11.785, 0.0085)
    if(grid.number==3) north.arrow(127.075, -11.38, 0.0085)
   
  }

  #'---------------------------------------------
  # Scale bar
  #'---------------------------------------------
  
 if(add.scale) {
   if(grid.number==1) custom.scalebar(pos = "bottomright")
   if(grid.number==2) custom.scalebar(pos = "bottomleft")
   if(grid.number==3) custom.scalebar(pos = "bottomright")
   
 }
  
  # raster::scalebar(1, type = "bar", divs = 1, below = "km", 
  #                  lonlat = TRUE, label = c("0", "", "1"), 
  #                  adj = c(1, -0.75), 
  #                  lwd = 10)
  
  
  if(save.to.pdf) {
    dev.off()
  }
  
  
} # End plot.grid

#'---------------------------------------------
# Modified version of prettymap::add.scalebar
#'---------------------------------------------

custom.scalebar <- function (plotunit = NULL, plotepsg = NULL, widthhint = 0.25, 
            unitcategory = "metric", htin = 0.1, padin = c(0.15, 0.15), 
            style = "bar", bar.cols = c("black", "white"), lwd = 1, linecol = "black", 
            tick.cex = 0.7, labelpadin = 0.08, label.cex = 0.8, label.col = "black", 
            pos = "bottomleft") 
  {
    params <- scalebarparams(plotunit = plotunit, plotepsg = plotepsg, 
                             widthhint = widthhint, unitcategory = unitcategory)
    extents <- params$extents
    bottomin <- graphics::grconvertY(extents[3], from = "user", 
                                     to = "inches")
    leftin <- graphics::grconvertX(extents[1], from = "user", 
                                   to = "inches")
    topin <- graphics::grconvertY(extents[4], from = "user", 
                                  to = "inches")
    rightin <- graphics::grconvertX(extents[2], from = "user", 
                                    to = "inches")
    ht <- graphics::grconvertY(bottomin + htin, from = "inches", 
                               to = "user") - extents[3]
    paduser <- graphics::grconvertX(leftin + labelpadin, from = "inches", 
                                    to = "user") - extents[1]
    if (pos == "bottomleft") {
      x <- graphics::grconvertX(leftin + padin[1], from = "inches", 
                                to = "user")
      y <- graphics::grconvertY(bottomin + padin[2], from = "inches", 
                                to = "user")
      adj <- c(0, 0)
      textadj <- c(0, 0.5)
      textx <- x + params$widthplotunit + paduser
      texty <- y + 0.5 * ht
    }
    else if (pos == "topleft") {
      x <- graphics::grconvertX(leftin + padin[1], from = "inches", 
                                to = "user")
      y <- graphics::grconvertY(topin - padin[2], from = "inches", 
                                to = "user")
      adj <- c(0, 1)
      textadj <- c(0, 0.5)
      textx <- x + params$widthplotunit + paduser
      texty <- y - 0.5 * ht
    }
    else if (pos == "topright") {
      x <- graphics::grconvertX(rightin - padin[1], from = "inches", 
                                to = "user")
      y <- graphics::grconvertY(topin - padin[2], from = "inches", 
                                to = "user")
      adj <- c(1, 1)
      textadj <- c(1, 0.5)
      textx <- x - params$widthplotunit - paduser
      texty <- y - 0.5 * ht
    }
    else if (pos == "bottomright") {
      x <- graphics::grconvertX(rightin - padin[1], from = "inches", 
                                to = "user")
      y <- graphics::grconvertY(bottomin + padin[2], from = "inches", 
                                to = "user")
      adj <- c(1, 0)
      textadj <- c(1, 0.5)
      textx <- x 
      texty <- y + 2.5 * ht
    }
    plotscalebar(x, y, ht, params, adj = adj, style = style, 
                 lwd = lwd, linecol = linecol, bar.cols = bar.cols, tick.cex = tick.cex)
    graphics::text(textx, texty, params$labeltext, adj = textadj, 
                   cex = label.cex, col = label.col)
  }


#'---------------------------------------------
# Function to fit models to bootstrap data
#'---------------------------------------------

compute.mod.probs <- function(bootstrap.data, k = NULL){
  
  message('Generating model sets ...')
  
  #'---------------------------------------------
  # Generate model sets
  #'---------------------------------------------
  
  pb <- dplyr::progress_estimated(length(bootstrap.data))
  
  blist.sr <- purrr::map(.x = bootstrap.data, 
                         .f = ~{
                           pb$tick()$print()
                           list.model.sets(response.var = 'sr', 
                                               dat = .x, 
                                               corr.cutoff = 0.7,
                                               max.predictors = 4,
                                               k = k,
                                               verbose = FALSE)})
  
  pb <- dplyr::progress_estimated(length(bootstrap.data))
  
  blist.max_n <- purrr::map(.x = bootstrap.data, 
                            .f = ~{
                              pb$tick()$print()
                              list.model.sets(response.var = 'max_n', 
                                                  dat = .x, 
                                                  corr.cutoff = 0.7,
                                                  max.predictors = 4,
                                                  k = k,
                                                  verbose = FALSE)})
  
  message('Fitting models ...')
  
  #'---------------------------------------------
  # Fit all models
  #'---------------------------------------------

  bmod.sr <- purrr::map(.x = blist.sr, 
                        .f = ~fit.model.set(model.set.list = .x, 
                                                max.models = 500, 
                                                save.model.fits = FALSE))
  
  bmod.max_n <- purrr::map(.x = blist.max_n, 
                           .f = ~fit.model.set(model.set.list = .x, 
                                                   max.models = 500, 
                                                   save.model.fits = FALSE))
  
  #'---------------------------------------------
  # Identify best models
  #'---------------------------------------------
  
  message('Identifying best models ...')
  
  bmod.sr.AIC <- purrr::map(.x = bmod.sr, 
                            .f = ~FSS.results(.x, AIC = FALSE,
                                              plot.models = FALSE, 
                                              verbose = FALSE)) %>% 
    purrr::map_chr(.x = ., .f = ~as.character(.x$best.model$summary$modname))
  
  bmod.max_n.AIC <- purrr::map(.x = bmod.max_n, 
                               .f = ~FSS.results(.x, AIC = FALSE,
                                                 plot.models = FALSE, 
                                                 verbose = FALSE)) %>% 
    purrr::map_chr(.x = ., .f = ~as.character(.x$best.model$summary$modname))
  
  message('Done!')
  
  return(list(sr = bmod.sr.AIC, max_n = bmod.max_n.AIC))
  
}

#'---------------------------------------------
# Function to divide a range into specified quantiles/percentiles
#'---------------------------------------------

cutrange <- function(input.raster, variable.name, percent){
  
  r <- raster::as.data.frame(input.raster, xy = TRUE, na.rm = TRUE)
  
  names(r) <- c("x", "y", "preds")
  
  percent.probs <- seq(0, 1, by = percent/100)
  percent.labels <- paste0(variable.name, "Q", 1:as.character(100/percent))
  
  # Creates df for raster plotting
  
  rlab <- cut(x = r$preds, 
              breaks = quantile(x = r$preds, 
                              probs = percent.probs, 
                              right = TRUE,
                              na.rm = TRUE),
              include.lowest = TRUE,
              labels = percent.labels)
  
  rnum <- cut(x = r$preds,
              breaks = quantile(r$preds, 
                              probs = percent.probs, 
                              right = TRUE,
                              na.rm = TRUE),
              include.lowest = TRUE)
  # 
  # 
  # quantval <<- data.frame(perc = sort(unique(rlab)),
  #                 yrast = seq(1/(100/percent),1, by = 1/(100/percent)))
  # assign(x = paste0('quantval.', variable.name), value = quantval, envir = .GlobalEnv)
  
  res <- cbind(r, rlab, rnum)
  names(res) <- c(names(r), paste0(variable.name, ".perc"),
                    paste0(variable.name, ".percnum"))
  
  return(res)
  
}

#'---------------------------------------------
# Function to calculate combined index of richness and abundance
#'---------------------------------------------

bivariate.index <- function(percentile = 10,
                            raster.one,
                            raster.two,
                            variable.one,
                            variable.two
                            ){

  proj.syst <- sp::proj4string(raster.one)
  
  #'-------------------------------------------------
  # Partition predictions into percentiles 
  #'-------------------------------------------------
  
  var.one <- cutrange(input.raster = raster.one, 
                      variable.name = variable.one, 
                      percent = percentile)
  
  var.two <- cutrange(input.raster = raster.two, 
                      variable.name = variable.two, 
                      percent = percentile)
  
  names(var.one)[3] <- "preds.one"
  names(var.two)[3] <- "preds.two"
  
  #'-------------------------------------------------
  # Combine percentile ranks and calculate normalised scores
  #'-------------------------------------------------

  # biv <- cbind(var.one, var.two[,!names(var.two)%in%c("x","y")]) %>% 
    
    biv <- dplyr::full_join(x = var.one, y = var.two, by = c("x", "y")) %>% 
    
    dplyr::mutate(comb = paste0(get(paste0(variable.one, '.perc')), 
                                '-', get(paste0(variable.two, '.perc')))) %>% 
    
    dplyr::mutate(score = as.numeric(sub(pattern = paste0(variable.one, "Q"), 
                                         replacement = "", 
                                         x = get(paste0(variable.one, '.perc')))) + 
                    as.numeric(sub(pattern = paste0(variable.two, "Q"),
                                   replacement = "", 
                                   x = get(paste0(variable.two, '.perc'))))) %>% 
    dplyr::mutate(nscore = (score - min(score, na.rm = TRUE))/ (max(score, na.rm = TRUE)-min(score, na.rm = TRUE)))

  r <- biv %>% dplyr::select(x, y, nscore) %>% 
    raster::rasterFromXYZ(xyz = ., crs = proj.syst)
  
  #'-------------------------------------------------
  # Compute number of combinations
  #'-------------------------------------------------
  
  biv.freqs <- data.frame(table(biv$comb)/nrow(biv))
  names(biv.freqs) <- c("comb", "freq")
  
  #'-------------------------------------------------
  # List all possible combinations
  #'-------------------------------------------------
  
  freq.grid <- expand.grid(sort(unique(var.one[,paste0(variable.one, '.perc')])), 
                           sort(unique(var.two[,paste0(variable.two, '.perc')])))
  
  names(freq.grid) <- c(paste0(variable.one, '.perc'), paste0(variable.two, '.perc'))
  
  #'-------------------------------------------------
  # Retrieve those combs that are in the data
  #'-------------------------------------------------
  
  freq.grid <- freq.grid %>% 
    dplyr::mutate(comb = paste0(get(paste0(variable.one, '.perc')), 
                                '-', get(paste0(variable.two, '.perc')))) %>% 
    dplyr::left_join(x = ., y = biv.freqs, by = 'comb') 
  
  freq.grid <- freq.grid %>% 
    dplyr::mutate(xcoord = (1/percentile) * as.numeric(sub(pattern = paste0(variable.one, "Q"), replacement = "", x = get(paste0(variable.one, '.perc'))))) %>% 
    dplyr::mutate(ycoord = (1/percentile) * as.numeric(sub(pattern = paste0(variable.two, "Q"), replacement = "", x = get(paste0(variable.two, '.perc')))))
  
  
  #'-------------------------------------------------
  # Convert to raster
  #'-------------------------------------------------
  
  rbiv <- data.frame(x = freq.grid$xcoord, 
                          y = freq.grid$ycoord, 
                          z = freq.grid$freq) %>% 
    raster::rasterFromXYZ(xyz = .) 
  
  return(list(raster = r, bivariate = rbiv))
  
}

#'---------------------------------------------
# Function to find the HEX colour code corresponding to an input colour 
# with a set opacity level (i.e. emulate transparency)
#'---------------------------------------------

hexa2hex <- function(input.colour, 
                     transparency, 
                     bg.colour = "white"){
  
  #'---------------------------------------------
  # PARAMETERS
  #'---------------------------------------------
  #' @param input.colour Initial colour.
  #' @param transparency Desired level of transparency (number between 0 and 1).
  #' @param bg.colour Colour of the background. Defaults to 'white'.
  #'---------------------------------------------
  
  # White background
  
  bg <- grDevices::col2rgb(bg.colour, alpha = FALSE)
  
  # Convert input colour to RGB
  
  rgbcol <- grDevices::col2rgb(input.colour, alpha = FALSE)
  
  # Calculate red, green, blue values corresponding to input colour at chosen transparency level
  
  rc <- (1 - transparency) * bg[1,] + transparency * rgbcol[1,]
  gc <- (1 - transparency) * bg[2,] + transparency * rgbcol[2,]
  bc <- (1 - transparency) * bg[3,] + transparency * rgbcol[3,]
  
  # Convert back to hex
  
  rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
  return(rgb2hex(r = rc, g = gc, b = bc))
}


#'---------------------------------------------
# Functions to stack rasters with different extents
#'---------------------------------------------

extend_all <- function(rasters){
  # https://gis.stackexchange.com/questions/252862/sum-rasters-with-different-extent-in-r
  extent(Reduce(extend,rasters))
}

stack_boots <- function(input.list){
  
  input.list$g1 <- purrr::map(.x = input.list[1:n.iter], 
                              .f = ~raster::extend(.x, y = extend_all(input.list[1:n.iter]))) %>% 
    raster::stack(.)
  
  input.list$g2 <- purrr::map(.x = input.list[(1+n.iter):(n.iter*2)], 
                              .f = ~raster::extend(.x, y = extend_all(input.list[(1+n.iter):(n.iter*2)]))) %>% 
    raster::stack(.)
  
  input.list$g3 <- purrr::map(.x = input.list[(1+(n.iter*2)):(n.iter*3)], 
                              .f = ~raster::extend(.x, y = extend_all(input.list[(1+(n.iter*2)):(n.iter*3)]))) %>% 
    raster::stack(.)
  
  input.list[1:(3*n.iter)] <- NULL
  
  return(input.list)
  
} # End stack_boots()

#'---------------------------------------------
# Raster statistics
#'---------------------------------------------

calc_mCV <- function(input.list, method = "median"){
  
  res <- list()
  if(method == "median"){
  res$m <- purrr::map(.x = input.list, 
                           .f = ~raster::calc(x = .x, fun = function(x) median(x, na.rm = TRUE)))}
  
  if(method == "mean"){
    res$m <- purrr::map(.x = input.list, 
                             .f = ~raster::calc(x = .x, fun = function(x) mean(x, na.rm = TRUE)))}
  
  res$rcv <- purrr::map(.x = input.list,
                        .f = ~raster::calc(x = .x, 
                                           fun = function(x){0.7413*IQR(x, na.rm = TRUE)/median(x, na.rm = TRUE)}))
  
  return(res)
  
} # End calc_mCV