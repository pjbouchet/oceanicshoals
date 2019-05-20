#' ====================================
# FUNCTIONS ====
#' ====================================
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
# Calculates species richness at different hierarchy levels
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
  
  # For each level (grouping) of the variable of interest
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
    
    # Adjusts MaxN values
    temp$adjMaxN <- temp$max_n
    temp$adjMaxN <- ifelse(!temp$sp.all=="SP", 0, temp$adjMaxN)
    
    
    # Uses the species inclusion rules to filter out all the non-"SP" records
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
# Rarefaction curves
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
# 2D and 3D NMDS plots
#'---------------------------------------------

plotnmds <- function(nmds.object, ndim=2, size.pts, input.data, plot.factor, plot.ellipse=T){
  
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