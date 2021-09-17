#' Mapping univariate and multivariate NSPOT results
#'
#' @description This function creates maps of univariate and multivariate
#' differences between q100 positive and negative values of the indices,
#' based on results from a call to the nspot() function.
#' @param sts A matrix or data frame containing (all) the stations information.
#' Must have, at least: ID, NAME, LON, LAT.
#' @param ind_names names of the indices.
#' @param res results object (list) from a call to the nspot() function.
#' @param path Directory where the results in pdf format will be saved.
#' @param brks Vector (numeric). Breaks in which indices values (not time) will
#' be categorized for representation purposes.
#' @param brks_time Vector (numeric). Breaks in which the temporal index values
#' (temporal trends) will be categorized for representation purposes.
#' @param country uses administrative boundaries of provinces (or equivalent) of
#' a specific country. If NULL, only country boundaries are shown.
#' @return collection of PDF files representing the univariate and multivariate
#' trends, by indices. 
#'
#' @import ggplot2
#' @import rnaturalearth
#' @import sf
#' @import raster
#' @import dplyr
#' @import reshape
#' @export

mapq  <- function(sts, ind_names, res, path, 
                  brks = c(1.25, 2, 3, 4), 
                  brks_time = c(1.1,1.2,1.3,1.5),
                  country = NULL) {
  
  # Preliminary issues ----
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  inds_fig <- c(ind_names, 'tme')
  sts <- sts[, c('ID', 'LON', 'LAT')]
  inds_nams <- c(paste0(inds_fig, 'sign'), 
                 paste0(inds_fig,'q100pos'), 
                 paste0(inds_fig,'q100neg'))
  sts[, inds_nams] <-  NA
  sts_uni <- sts_multi <- sts
  
   # Extracting data from res ----
  for (i in 1:length(res$id)) {
    w <- which(as.character(sts$ID) == res$id[i])
    if(length(w) == 0) next else{
      sts_uni[w, grep('sign', names(sts_uni))] <- res$m1sign[, i] <= 0.05
      sts_uni[w, grep('q100pos', names(sts_uni))] <- res$m1q100pos[, i]
      sts_uni[w, grep('q100neg', names(sts_uni))] <- res$m1q100neg[, i]
      sts_multi[w, grep('sign', names(sts_multi))] <- res$m2sign[, i] <= 0.05
      sts_multi[w, grep('q100pos', names(sts_multi))] <- res$m2q100pos[, i]
      sts_multi[w, grep('q100neg', names(sts_multi))] <- res$m2q100neg[, i]
    }
  }
  
  # Map constants ----
  if(is.null(country)){
    cnt <- ne_countries(scale = 10, type = "countries", returnclass = "sf")
  } else{
    if(class(country)!='character'){
      stop('You must provide a valid country name')
    } else{
      cnt <- ne_states(country = tolower(country), returnclass = "sf")
      }
    }
  sts_uni <- sts_uni[complete.cases(sts_uni), ]
  sts_multi <- sts_multi[complete.cases(sts_multi), ]
  sts_sf <- st_as_sf(sts_uni, coords = c("LON", "LAT"), crs = 4326)
  esa <- extent(sts_sf) * 1.20
  
  
  # Make the plots (Santiago) ----
  sts_uni2 <- merge(
    sts_uni %>%
      melt(id.vars = c('ID', 'LON', 'LAT'),
         measure.vars = paste0(inds_fig,'q100pos')) %>%
      dplyr::rename(q100pos = value) %>%
      dplyr::mutate(across(contains('variable'), ~gsub('q100pos', '', .))),
    sts_uni %>%
      melt(id.vars = c('ID', 'LON', 'LAT'),
         measure.vars = paste0(inds_fig,'q100neg')) %>%
      dplyr::rename(q100neg = value) %>%
      dplyr::mutate(across(contains('variable'), ~gsub('q100neg', '', .)))
  ) %>%
  dplyr::mutate(q100ratio = ifelse(q100pos>q100neg,
                                   q100pos/q100neg, -q100neg/q100pos))
  sts_uni2 <- merge(
    sts_uni %>%
      melt(id.vars = c('ID', 'LON', 'LAT'),
      measure.vars = paste0(inds_fig, 'sign')) %>%
      dplyr::rename(sign = value) %>%
      dplyr::mutate(across(contains('variable'), ~gsub('sign', '', .))),
    sts_uni2) %>%
    dplyr::mutate(
      class = factor(ifelse(q100pos > q100neg, 'up', 'down'))
    )
  #head(sts_uni2)
  #summary(sts_uni2)
  
  sts_multi2 <- merge(
    sts_multi %>%
      melt(id.vars = c('ID', 'LON', 'LAT'),
           measure.vars = paste0(inds_fig,'q100pos')) %>%
      dplyr::rename(q100pos = value) %>%
      dplyr::mutate(across(contains('variable'), ~gsub('q100pos', '', .))),
    sts_multi %>%
      melt(id.vars = c('ID', 'LON', 'LAT'),
           measure.vars = paste0(inds_fig,'q100neg')) %>%
      dplyr::rename(q100neg = value) %>%
      dplyr::mutate(across(contains('variable'), ~gsub('q100neg', '', .)))
  ) %>%
    dplyr::mutate(q100ratio = ifelse(q100pos>q100neg,
                                     q100pos/q100neg, -q100neg/q100pos))
  sts_multi2 <- merge(
    sts_multi %>%
      melt(id.vars = c('ID', 'LON', 'LAT'),
           measure.vars = paste0(inds_fig, 'sign')) %>%
      dplyr::rename(sign = value) %>%
      dplyr::mutate(across(contains('variable'), ~gsub('sign', '', .))),
    sts_multi2) %>%
    dplyr::mutate(
      class = factor(ifelse(q100pos > q100neg, 'up', 'down'))
    )
  
  # a function to make the plots
  myplot <- function(mode, variable, sts_uni2, sts_multi2) {
    if (mode == 'univariate'){
      data <- sts_uni2
      hidden <- sts_multi2
    } else if(mode == 'multivariate'){
      data <- sts_multi2
      hidden <- sts_uni2
    }
    hidden$LON <- hidden$LON*10000
    hidden$LAT <- hidden$LAT*10000
    w <- data$variable == variable
    ggplot(data[w,]) +
      # reference map
      geom_sf(data = cnt, fill= "antiquewhite") +
      coord_sf(crs = "+proj=longlat +datum=WGS84 +no_defs", 
               xlim=c(esa[1], esa[2]), ylim=c(esa[3],esa[4])) +
      # triangles
      geom_point(aes(x = LON, y = LAT,
                     size = abs(q100ratio),
                     fill = class,
                     shape = class,
                     color = sign,
                     alpha = sign)) +
      geom_point(data = hidden[w,],
                 aes(x = LON, y = LAT,
                     size = abs(q100ratio))) +
      scale_size_continuous(name = 'N times', #labels = brkk, breaks = brkk,
                            range = c(1,10), limits = c(NA,NA)) +
      scale_fill_manual(values=c('blue', 'red')) +
      scale_shape_manual(name='Sign', values=c(25, 24), labels = c('neg.','pos.')) +
      scale_color_manual(values=c('NA', 'white')) +
      scale_alpha_manual(values=c(0.5, 1)) +
      theme_bw() +
      theme(legend.position = 'bottom', 
            panel.grid.major = element_line(color = gray(.5), 
                                            linetype = 'dashed', 
                                            size = 0.1),
            panel.background = element_rect(fill = 'aliceblue'),
            axis.title = element_blank()) +
      guides(size = guide_legend(override.aes=list(shape = 24)),
             shape = guide_legend(override.aes=list(shape = c(25,24),
                                                    fill = c('blue','red'))),
             fill = 'none',
             color = 'none',
             alpha = 'none')
  }
  
  # loop across variables
  for(ind in inds_fig) {
    ggsave(filename = paste(path, paste0('univariate_',ind,'.pdf'),sep='/'),
           plot = myplot(mode = 'univariate', variable = ind, sts_uni2, sts_multi2),
           device = cairo_pdf)
    ggsave(filename = paste(path, paste0('multivariate_',ind,'.pdf'),sep='/'),
           plot = myplot(mode = 'multivariate', variable = ind, sts_uni2, sts_multi2),
           device = cairo_pdf)
  }
}
  
