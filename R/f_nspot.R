#' Non-Stationary Peak-over-threshold analysis (NSPOT)
#'
#' @description This function analyzes the individual and combined influence of climatic indices on extreme precipitation events.
#' @param pcp Vector containing the ID of the station (first element) and the daily precipitation data series (rest of the elements).
#' @param ind Matrix or data.frame with daily indices, by columns.
#' @param sts Matrix or data.frame containing (all) the stations information. Must have, at least: ID, NAME, LON, LAT.
#' @param tstep Time step considered for the analysis. It can be a character ("annual", default) or a numeric vector containing the months.
#' @param thres numeric vector indicating the percentile (from 0 to 1) over which events will be considered for the analysis.
#' @param var Variable to analyze computed from the extreme precipitation events. Possible values are: "magnitude" (total precipitation of the event); "intensity" (mean daily precipitation of the event); "maximum_intensity" (maximum daily value of the event) and "duration" (number of days of the event).
#' @param path Directory where results in pdf will be saved. If NULL, no files will be created.
#' @param isEvent Logical. If TRUE, precipitation events (consecutive rainy days) will be considered instead of single rainy days (FALSE).
#' @param method Method used for fitting variables. At this moment, only "gpd" (Maximum-likelihood Fitting for the GPD Model) is available.
#' @param ind_names Character vector with the names of the indices to consider in the analysis. The names must match with one or more column names in "ind".
#' @return A list with resulting parameters is returned. If path is defined, a pdf file is created. 
#'
#'
#'
#' @examples
#'
#
#'
#' @import ggplot2
#' @import stats
#' @import ismev
#' @importFrom evd qgpd
#' @importFrom quantreg rq
#' @importFrom rnaturalearth ne_countries
#' @importFrom sf st_as_sf
#' @importFrom raster extent
#' @export



nspot <- function(pcp, ind, sts, dates, tstep = 'annual', thres = 0.9, 
                  var = 'intensity', path=NULL, isEvent = TRUE, 
                  method = 'gpd', ind_names) {
  
  # Preliminary checks ----
  
  if(!(is.matrix(ind) | is.data.frame(ind))){
    stop("'ind' must be a matrix or a data.frame.")
  } else {
    if(!all(ind_names%in%colnames(ind))){
      stop("'ind_names' must coincide with one or more colnames of 'ind'.")
    }
  }
  if(length(which(sts$ID == pcp[1])) == 0){
    stop("Station ID is not in stations list")
  }
  if(length(dates) != length(pcp[-1])){
    stop("Length of dates does not coincide with length of the data.")
  }
  if(is.character(tstep)){
    if(tstep != "annual"){
      stop("'tstep' only can be 'annual' or a numeric vector with months.")
    }
  } else{
    if(!is.numeric(tstep)){
      stop("'tstep' only can be 'annual' or a numeric vector with months.")
    } else{
      if(min(tstep) < 1 | max(tstep) > 12){
        stop("'tstep' must be between 1 and 12")
      }
    }
  }
  if(!is.numeric(thres)){
    stop("'thres' must be a numeric element between 0 and 1.")
  } else{
    if(min(thres) < 0.000000001 | max(thres) > 1){
      stop("'thres' must be a numeric element between 0 and 1.")
    }
  }
  if(!any(var == c('intensity', 'magnitude', "maximum_intensity", 
                   "duration"))){
    stop("'var' is not correct.")
  }
  if(!is.null(path)){
    if(!is.character(path)){
      stop("'path' must be a path to a directory.")
    }
  }
  
  # Setting station name and data
  est <- pcp[1]
  pcp <- as.numeric(as.character(pcp[-1]))
  stsname <- as.character(sts$NAME[which(sts$ID == est)])
  stsname <- gsub('/','_',stsname)
  stsname <- gsub('[^ -~]','_',stsname)
  years <- as.numeric(unique(substr(dates, 1, 4)))

  
  # Create event series ----
  if (isEvent) {
    # find events
    evs <- findevents(pcp)
    if(nrow(evs) < length(unique(substr(dates, 1, 4)))){
      return(list(
        id = est,
        nam = stsname,
        m0=NA, m0aic=NA, m0q100=NA,
        m1=NA, m1aic=NA, dtest1=NA, m1q100pos=NA, m1q100neg=NA, m1sign=NA,
        m2=NA, m2aic=NA, dtest2=NA, m2q100pos=NA, m2q100neg=NA, m2sign=NA))
    }
    dat <- data.frame(day = rep(NA, nrow(evs)), intmax = NA, mag = NA, 
                      dur = NA, int = NA, mon=NA, yea=NA)
    dat[, c(ind_names, paste0('u.', ind_names), 'u.tme')] <- NA
    mind <- match(ind_names, names(dat))
    
    for(i in 1:nrow(evs)) {
      a <- evs[i,1]
      b <- evs[i,2]
      dat$day[i] <- as.character(dates[a])
      dat$intmax[i] <- max(pcp[a:b])
      dat$mag[i] <- sum(pcp[a:b])
      dat$dur[i] <- b - a + 1
      dat$int[i] <- dat$mag[i]/dat$dur[i]
      for(h in 1:length(ind_names)){
        wtel <- which(ind_names[h] == colnames(ind))
        dat[i, mind[h]] <- mean(ind[a:b, wtel])
        dat[i, (mind[h] + length(mind))] <- 
          ifelse(!is.na(max(ind[a:b, wtel])),
                 ind[a:b,
                     wtel][which(abs(ind[a:b, 
                                         wtel]) == max(abs(ind[a:b, 
                                                               wtel])))], NA)
      }
    }
    # dat <- dat[!is.na(dat$mag),]
  }
  
  # or else create daily series
  if (!isEvent) {
    var <- 'intensity'
    dat <- data.frame(day = dates, int = pcp, mon = NA, yea = NA)
    dat[, c(ind_names, paste0('u.', ind_names), 'u.tme')] <- NA
    dat <- dat[dat[, 2] > 0, ]
  }
  
  dat$mon <- as.numeric(substr(dat$day, 6,7))
  dat$yea <- as.numeric(substr(dat$day, 1,4))
  dat$yea_p <- scale(dat$yea)[,1]
  
  # Select events ----
  if(!is.character(tstep)){
    # keeps only selected months
    fw <- function(x, y) which(x == y)
    dat <- dat[sort(unlist(sapply(tstep, fw, y = dat$mon))), ]
  }
  
  # Constants ----
  
  # create output
  dir.create(paste0(path,'/stations/'), recursive = TRUE, showWarnings = FALSE)
  # names of the dependent variables
  deps <- c('magnitude', 'intensity', 'maximum_intensity', 'duration')
  # their column numbers
  depsn <- c(3, 5, 2, 4)
  # their units
  depu <- c('mm', 'mm day-1', 'mm day-1', 'days')
  # column of the dependent variable being analysed
  yn <- depsn[which(deps == var)]
  # names of covariates
  covs <- c(ind_names, 'Tme')
  # their column numbers
  covsn <- c(match(covs[1:(length(covs) - 1)], names(dat)), which(names(dat) == 'yea'))
  # linear versions of teleconnection indices and their column numbers
  covs_p <- paste0(covs[1:(length(covs) - 1)], '_p')
  dat[, covs_p] <- NA
  covsn_p <- c(match(covs_p[1:length(covs_p)], names(dat)), which(names(dat) == 'yea_p'))
  for(jj in 1:(length(covsn_p)-1)){
    # dat[, covsn_p[jj]] <- pnorm(dat[, c(mind, max(mind) + (1:length(mind)))[jj]]) - 0.5
    dat[, covsn_p[jj]] <- pnorm(dat[, c(mind, max(mind) + (1:length(mind)))[jj]])
  }
  # names of thresholds
  thrs <- names(dat)[grep('u\\.', names(dat))]
  # their column numbers
  thrsn <- grep('u\\.', names(dat))
  
  # outputs
  m0aic <- NA
  m1aic <- m2aic <- vector('numeric', length(covs))
  m1q100pos <- m1q100neg <- vector('numeric', length(thrs))
  m2q100pos <- m2q100neg <- vector('numeric', (length(ind_names) + 1))
  dtest1 <- dtest2 <- m1scale <- m1sign <- vector('numeric', length(covs))
  
  # location
  loc <- quantile(dat[,yn], thres)[[1]]
  # rate (Poisson's lambda)
  rate  <- sum(dat[,yn]>loc) / (sum(!is.na(pcp))/365.25)
  
  ####
  # Start pdf ----
  ####
  pdf(paste(path,'/stations/',stsname,'_',est,'.pdf',sep=''), width=10)
  
    # Location map
    cnt <- ne_countries(scale = 10, type = "countries", returnclass = "sf")
    sts_sf <- st_as_sf(sts, coords = c("LON", "LAT"), crs = 4326)
    sts_sf_ref <- sts_sf[-which(sts$ID == est),]
    sts_sf_cand <- sts_sf[which(sts$ID == est),]
    esa <- extent(sts_sf) * 1.10
    p <- ggplot() + 
      geom_sf(data = cnt,fill= "antiquewhite") +
      geom_sf(data = sts_sf_ref, fill = 'grey', color = 'grey', alpha = 0.5, 
              size=2) +
      geom_sf(data = sts_sf_cand, fill = 'red', color = 'red', size = 3) +
      ggtitle(paste(stsname,', ',est,sep='')) +
      coord_sf(crs = "+proj=longlat +datum=WGS84 +no_defs", 
               xlim = c(esa[1], esa[2]), ylim = c(esa[3], esa[4])) +
      theme_bw() +
      theme(legend.position = 'none', 
            panel.grid.major = element_line(color = gray(.5), 
                                            linetype = 'dashed', size = 0.1),
            panel.background = element_rect(fill = 'aliceblue'))
    print(p)
  
  # m0: stationary model ----
  
  par(mfrow=c(1,2),mar=c(4,4,2,2))
  # MRL plot
  mrl.plot(dat[, yn], quantile(dat[, yn], 0.80), quantile(dat[, yn], 0.99), 
           nint = 250)
  title(paste('Threshold (u = q',thres,' = ',round(loc,2),', n = ', 
              sum(dat[, yn] > loc),')', sep = ''))
  abline(v = loc, lty = 'dotted')
  # fit it
  m0 <- gpd.fit(dat[, yn], loc, rate, show = FALSE)
  if (method == 'gpd') {
    sca <- m0$mle[1]
    sha <- m0$mle[2]
  } else {
    m0 <- pp.fit(dat[, yn], loc, rate, muinit = loc[[1]], 
                 siginit = m0$mle[1], show=FALSE)
    sha <- m0$mle[3]
    sca <- m0$mle[2] + m0$mle[3] * (loc[[1]] - m0$mle[1])
  }
  m0aic <- round(2 * length(m0$mle) + 2 * m0$nllh, 0)
  m0q100 <-  qgpd(1 - 1 / (rate * 100), loc, sca, sha)
  
  # plot: IDF curve
  # .. model
  y <- log10(1.1:100.1)
  pexc <- 1 - (1 / (rate * (10^y)))
  y <- y[-pexc < 0]
  pexc <- pexc[-pexc < 0]
  z <- qgpd(pexc, loc, sca, sha)
  plot(y~z, type = 'l', yaxt = 'n', xlab = paste(var,' (',depu,')',sep=''),
       ylab='Return period (years)')
  axis(2, log10(c(1, 2, 5, 10, 25, 50, 100)), c(1, 2, 5, 10, 25, 50, 100))
  grid()
  abline(h = log10(c(1, 2, 5, 10, 25, 50, 100)), col = 'lightgrey', 
         lty = 'dotted')
  title(paste('Stationary model (AIC=', 
              round(2 * length(m0$mle) + 2 * m0$nllh, 0), ')', sep = ''))
  # .. empirical return periods
  emp <- ecdf(dat[dat[, yn] > loc, yn])
  yy <- log10(1 / rate / (1 - emp(dat[dat[, yn] > loc, yn])))
  zz <- dat[dat[, yn] > loc, yn]
  points(yy ~ zz)
  lines(y ~ z)
  
  # m1: NSEV, univariate ----
  
  # threshold values
  # for each index
  tar <- paste0('u.', covs)
  thrModels <- list()
  for(idx in 1:length(covs_p)){
    u.mod <- rq(as.formula(paste0(colnames(dat)[yn], "~", covs_p[idx])),
                tau = thres, data = dat)
    dat[, which(names(dat) == tar[idx])] <- u.mod$fitted.values
    thrModels[[idx]] <- u.mod
  }
  # for time
  u.mod <- rq(as.formula(paste0(colnames(dat)[yn], "~yea_p")), 
              tau = thres, data = dat)
  dat$u.tme <- u.mod$fitted.values
  thrModels[[(length(thrModels) + 1)]] <- u.mod
  names(thrModels) <- tar
  
  # plot it
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
  for (cov in covs) {
    k <- which(covs == cov)
    # covariate column
    xn <- covsn[k]
    if (cov != 'Tme') {
      covariate <- as.matrix(pnorm(dat[,xn]))
    } else {
      covariate <- as.matrix(dat$yea_p)
    }
    # threshold column
    un <- thrsn[k]
    # threshold
    threshold <- dat[,un]
    # threshold model
    thrModel <- thrModels[k][[1]]
    
    # fit it
    if (method == 'gpd') {
      m1 <- gpd.fit(dat[, yn], threshold, rate, covariate, 
                    sigl = 1, siglink = exp, show = FALSE)
      sca <- m1$mle[1:2]
      sha <- m1$mle[3]
    } else {
      m1 <- pp.fit(dat[, yn], threshold, rate, covariate, 
                   mul = 1, maxit = 1000)
      sha <- m1$mle[4]
      sca <- exp(log(m1$mle[3]) - m1$mle[4] * log(rate))
    }
    m1aic[k] <- round(2 * length(m1$mle) + 2 * m1$nllh, 0)
    xs <- c(-2.5,-2,-1,-0.5,0,0.5,1,2,2.5)
    if (cov != 'Tme') {
      newd <- data.frame(newd = pnorm(xs))
    } else {
      newd <- data.frame(newd = unique(dat$yea_p))
    }
    names(newd) <- names(thrModel$coefficients)[2]
    loc <- predict(thrModel, newd)[c(1, nrow(newd))]
    m1q100pos[k] <-  qgpd(1-1/(rate*100), loc[2], exp(sca[1] + pnorm(2.5) * sca[2]), sha)
    m1q100neg[k] <-  qgpd(1-1/(rate*100), loc[1], exp(sca[1] + pnorm(-2.5) * sca[2]), sha)
    kk <- summary(thrModel)$coeff
    m1sign[k] <- ifelse(dim(kk)[2] == 4, kk[2, 4], 1)
    m1scale[k] <- m1$mle[2]
    dtest1[k] <- D.test(m1, m0, sign = 0.05)
    
    # plot it
    # threshold vs. index
    w <- dat[, yn] > threshold
    plot(dat[, yn] ~ covariate, xaxt='n',
         ylab = paste0(var, ' (', depu[which(var == deps)], ')'), 
         xlab = paste0(cov, ' (AIC = ',round(AIC(thrModel)[1]),
                       ifelse(m1sign[k]<=0.05,'*',''),')'))
    grid()
    points(dat[w, yn] ~ covariate[w], col = 'red')
    abline(thrModel, lty = 'dashed', col = 'red')
    if (cov!='Tme') {
      axis(1)
    } else {
      ss <- round(seq(0, length(years), length(years)/5))
      axis(1, unique(covariate)[c(1,ss[-1])], ss)
    }
    text(max(covariate), max(dat[w,yn]), paste0('n=',sum(w)),
         pos=2, col='red')
    
    # scale vs. index
    if (cov != 'Tme') {
      x <- seq(-2.5, 2.5, 0.25)
      y <- exp(sca[1] + pnorm(x) * sca[2])
    } else {
      x <- newd[,1]
      y <- exp(sca[1] + x * sca[2])
    }
    plot(y~x, type = 'l', ylab = 'Scale parameter', xaxt = 'n',
         xlab = paste0(cov,' (AIC=', m1aic[k],
                       ifelse(m1sign[k + length(covs)], '*', ''), ')'))
    grid()
    if (cov!='Tme') {
      axis(1)
    } else {
      ss <- round(seq(0, length(years), length(years)/5))
      axis(1, newd[, 1][c(1, ss[-1])], ss)
    }
    
    # IDF curves for values of covariate
    y <- log10(1.1:100.1)
    pexc <- 1-(1/(rate*(10^y)))
    y <- y[-pexc < 0]
    pexc <- pexc[-pexc < 0]
    z <- qgpd(pexc, predict(thrModel, newd)[1],
              exp(sca[1] + newd[,1]*sca[2]), sha)
    if (cov!='Tme') {
      xl <- c(
        qgpd(pexc[c(1, length(pexc))], predict(thrModel, newd)[3],
             exp(sca[1] + pnorm(-2) * sca[2]), sha),
        qgpd(pexc[c(1, length(pexc))], predict(thrModel, newd)[4],
             exp(sca[1] + pnorm(2) * sca[2]), sha)
      )
    } else {
      xl <- c(
        qgpd(pexc[c(1, length(pexc))], predict(thrModel, newd)[length(loc)],
             exp(sca[1] + newd[length(loc), 1] * sca[2]), sha),
        qgpd(pexc[c(1, length(pexc))], predict(thrModel, newd)[1],
             exp(sca[1] + newd[1,1] * sca[2]), sha)
      )
    }
    xl <- c(min(xl), max(xl))
    plot(y~z, type = 'n', yaxt = 'n', 
         xlab = paste(var, ' (', depu[which(var == deps)],')', sep = ''),
         ylab = 'Return period (years)',
         xlim = xl)
    axis(2, log10(c(1, 2, 5, 10, 25, 50, 100)), c(1, 2, 5, 10, 25, 50, 100))
    grid(ny = 0)
    abline(h = log10(c(1, 2, 5, 10, 25, 50, 100)), 
           col = 'lightgrey', lty = 'dotted')
    if (cov!='Tme') {
      for (i in c(2:8)) {
        z <- qgpd(pexc, predict(thrModel, newd)[i],
                  exp(sca[1] + pnorm(xs[i]) * sca[2]), sha)
        lines(y~z, type = 'l')
        text(z[length(pexc)], 2.04, xs[i], cex = 0.75)
      }
    } else {
      for (i in c(1,length(loc))) {
        z <- qgpd(pexc, predict(thrModel, newd)[i],
                  exp(sca[1] + newd[i,1] * sca[2]), sha)
        lines(y~z, type = 'l')
        text(z[length(pexc)], 2.04, paste0('t=', names(loc[i])), cex = 0.75)
      }
    }
    
    # bivariate plot
    if (xn < 10) {
      x <- seq(-2.5, 2.5, 0.25)
    } else {
      x <- seq(-5, 5, 0.25)
    }
    if (cov == 'Tme') {
      x <- newd[, 1]
    }
    y <- log10(1.1:100.1)
    z <- NULL
    if (cov != 'Tme') {
      newd <- data.frame(newd = pnorm(x))
      names(newd) <- names(thrModel$coefficients)[2]
    }
    for (i in 1:length(x)) {
      z <- c(z, qgpd(1-(1/(10^y)), predict(thrModel, newd)[i],
                     exp(sca[1] + newd[i, 1] * sca[2]), sha))
    }
    z <- matrix(z, ncol = length(x))
    z <- t(z)
    image(x, y, z, xlab=cov, yaxt = 'n', xaxt = 'n',
          ylab = 'Return period (years)', col = cm.colors(50))
    if (cov != 'Tme') {
      axis(1)
    } else {
      ss <- round(seq(0, length(years), length(years)/5))
      axis(1, unique(covariate)[c(1, ss[-1])], ss)
    }
    axis(2, log10(c(1, 2, 5, 10, 25, 50, 100)), c(1, 2, 5, 10, 25, 50, 100))
    contour(x, y, z, add = TRUE, nlevels = 20, box = 'plot')
    grid(ny = NA, col = 'grey')
    abline(h = log10(c(1, 2, 5, 10, 25, 50, 100)), 
           col = 'grey', lty = 'dotted')
    par(new = TRUE)
    plot(y~1, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
    
  } # next cov
    
  # m2: NSEV, multivariate ----
  
  # threshold values
  covs_p <- paste0(ind_names, '_p')
  mrq <- rq(as.formula(paste0(colnames(dat)[yn], "~", 
                              paste(c(covs_p, "yea_p"), collapse = '+'))), 
            tau = thres, data = dat)
  dat$u.multi <- mrq$fitted.values
  
  # fit it
  inds_fig <- c(ind_names, 'Tme')
  if (method == 'gpd') {
    threshold <- dat$u.multi
    covariate <- as.matrix(dat[, c(covs_p, 'yea_p')])
    m2 <- gpd.fit(dat[, yn], threshold, rate, covariate, 
                  sigl = c(1:ncol(covariate)), # number of covs + yea_p
                  siglink = exp,
                  siginit = c(m1$mle[1], m1scale[c(1:length(covs_p), length(m1scale))]),
                  shinit = m1$mle[3], show = FALSE)
    sca <- m2$mle[1:(length(covs_p) + 2)]
    sha <- m2$mle[length(m2$mle)]
  } else {
    #m2 <- pp.fit(dat[,yn], threshold, rate, covariate, mul=c(1:3), maxit=1000)
    #sha <- m1$mle[1:5]
    #sca <- exp(log(m1$mle[3])-m1$mle[4]*log(rate))
  }
  m2aic <- round(2 * length(m2$mle) + 2 * m2$nllh, 0)
  kk <- summary(mrq)$coeff
  if (dim(kk)[2] == 5) {
    m2sign <- kk[2:(length(inds_fig) + 1), 4]
  } else {
    m2sign <- as.numeric(rowSums(sign(kk[2:(length(inds_fig) + 1), 2:3])) == 0)
  }
  # m2sign[(length(inds_fig) + 1):(length(inds_fig) * 2)] <- pnorm(0, abs(m2$mle[2:(length(inds_fig)+1)]), m2$se[2:(length(inds_fig)+1)])
  
  # insertado por Roberto ----
  summ <- summary(mrq, se = "boot")
  m2sign <- summ$coefficients[2:(length(inds_fig) + 1),4]
  
  dtest2 <- D.test(m2, m0, sign=0.05)
  
  # plot it
  x <- y <- seq(-2.5, 2.5, 0.25)
  tt <- c(1,seq(0, length(years), 5)[-1])
  newd <- eval(parse(text = paste0("expand.grid(", 
                                   paste(rep("x", 
                                             length(ind_names)), 
                                         collapse = ','), 
                                   ",tt)")))
  names(newd) <- c(ind_names,'yea')
  covs_p <- c(paste0(ind_names[1:length(ind_names)], '_p'), 'yea_p')
  newd[, covs_p] <- NA
  for(jj in 1:(length(covs_p)-1)){
    newd[, which(covs_p[jj] == names(newd))] <- pnorm(newd[, jj])
  }
  newd$yea_p <- scale(newd$yea)[,1]
  newd$thres <- predict(mrq, newd)
  nn <- paste0("newd$", covs_p)
  newd$scale <- eval(parse(text = paste0("exp(", sca[1],'+', 
                                         paste(paste(nn,
                                                     sca[2:length(sca)], 
                                                     sep = '*'), 
                                               collapse = '+'), ")")
               )
         )
  pexc <- 1 - (1 / (rate * (100)))
  newd$p100 <- qgpd(pexc, newd$thres, newd$scale, sha)
  midyear <- tt[floor(length(tt) / 2)] # mid year
  mind <- match(ind_names, names(newd)) #columns with vals
  m2q100pos <- m2q100neg <- numeric()
  for(i in mind){
    # pos
    w2 <- newd[which(newd[, mind[i]] == 2), ]
    w0 <- w2[which(rowSums(abs(w2[,mind[-i]])) == 0),]
    m2q100pos[i] <- w0$p100[which(w0$yea == midyear)]
    # neg
    wm2 <- newd[which(newd[, mind[i]] == -2), ]
    wm0 <- wm2[which(rowSums(abs(wm2[,mind[-i]])) == 0),]
    m2q100neg[i] <- wm0$p100[which(wm0$yea == midyear)]
  }
  w <- which(rowSums(abs(newd[,mind])) == 0)
  newd_yea <- newd[w,]
  m2q100pos <- c(m2q100pos, newd_yea$p100[which(newd_yea$yea == tt[length(tt)])])
  m2q100neg <- c(m2q100neg, newd_yea$p100[which(newd_yea$yea == tt[1])])
  
  # close pdf ----
  dev.off()
  
  # return a list with important values
  return(list(
    id = est,
    nam = stsname,
    m0=m0, m0aic=m0aic, m0q100=m0q100,
    m1=m1, m1aic=m1aic, dtest1=dtest1, m1q100pos=m1q100pos, m1q100neg=m1q100neg, m1sign=m1sign,
    m2=m2, m2aic=m2aic, dtest2=dtest2, m2q100pos=m2q100pos, m2q100neg=m2q100neg, m2sign=m2sign))
  
}
