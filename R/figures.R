
require(akima)
require(smatr)

defaultColours <- function() {
  c("#577A33", "#C951CA", "#678BBD", "#CC4830", "#5F3422", "#54345B", "#71D49A",
    "#C58B37", "#C6BF93", "#475C54", "#8CC9CC", "#CDD24E", "#C27B72", "#7C6AC9",
    "#CE99C1")
}

default_colour_pallete <- function(i = 1) {
  if (i == 1) {
    p <- list(x = c(rep(-1.1, 4), rep(-0.75, 4), rep(-0.6, 4), rep(-0.3, 4)),
      y = rep(c(1.4, 1.2, 0.8, 0.3), 4), cols = c("#3bd02b", "#197902", "#82603e",
        "#82603e", "#1fc6d0", "#185299", "#8456b3", "#8456b3", "#1fc6d0",
        "#185299", "#8456b3", "#8456b3", "#f8fb20", "#ff1212", "#ff29c0",
        "#ff29c0"))
  }
  p
}

color_2D_pallete <- function(xo, yo, space = default_colour_pallete()) {

  # translate cols into RGB
  z <- col2rgb(space$cols)
  # get new rgb colors by interpolating between existing colours

  xo <- pmax(pmin(xo, max(space$x)), min(space$x))
  yo <- pmax(pmin(yo, max(space$y)), min(space$y))

  R <- interp(space$x, space$y, z[1, ], xo, yo, linear = TRUE, extrap = FALSE)$z
  G <- interp(space$x, space$y, z[2, ], xo, yo, linear = TRUE, extrap = FALSE)$z
  B <- interp(space$x, space$y, z[3, ], xo, yo, linear = TRUE, extrap = FALSE)$z
  rgb(R, G, B, maxColorValue = 255)
}

getColours <- function(obj) {
  sapply(seq_len(nrow(obj$strat)), function(i) color_2D_pallete(log10(obj$strat$Tr0[i]),
    log10(obj$strat$Tr1[i])))
}


generatePaths <- function(base.dir, T, P, stringFormat = "[%s,%s]", times = FALSE) {
  paths <- dirs <- titles <- NULL
  for (i in 1:length(T)) for (j in 1:length(P)) {
    k <- matrixIndexAsVectorIndex(i, j, length(P))
    dirs[k] <- file.path(base.dir, sprintf(stringFormat, T[i], P[j]))
    x <- getFitnessFilesInDir(dirs[k], 200, 400, 2, "")
    if (length(x$paths) > 0){
      paths[k] <- last(x$paths)
      titles[k] <- paste(P[j],T[i])
    }
  }
  titles[is.na(titles)] <- ""
  list(dirs = dirs, paths = paths, titles = titles)
}


# load resident data
loadResident <- function(DIR) {
  x <- list()
  x$strat <- read.table(file.path(DIR, "Strategy.txt"), head = TRUE)
  x$pAge <- read.table(file.path(DIR, "patch_age.txt"), head = TRUE)[, 1:3]
  x$SpAge <- read.table(file.path(DIR, "patch_age.txt"), head = TRUE)[, 1:3]

  tmp <- readLines(file.path(DIR, "params.m"))[-c(1)]
  tmp2 <- do.call(rbind, strsplit(gsub(";", "", tmp), "=", fixed = TRUE))
  x$params <- data.frame(value = as.numeric(tmp2[, 2]), row.names = tmp2[, 1],
    stringsAsFactors = FALSE)

  x$nSpp <- length(x$strat[, 1])
  x$Age <- x$Size <- x$Density <- x$Popn <- x$Lambda <- NULL
  for (i in 1:x$nSpp) {
    x$Age[[i]] <- read.table(file.path(DIR, paste0(i - 1, "_age.txt")), head = FALSE)
    x$Size[[i]] <- read.table(file.path(DIR, paste0(i - 1, "_bound_m.txt")),
      head = FALSE)
    x$Density[[i]] <- read.table(file.path(DIR, paste0(i - 1, "_bound_n.txt")),
      head = FALSE)
    x$Pop[[i]] <- read.table(file.path(DIR, paste0(i - 1, "_popn.txt")), head = TRUE)
    x$Lambda[[i]] <- read.table(file.path(DIR, paste0(i - 1, "_coh_lam.txt")),
      head = FALSE)
  }
  x
}

# reorders Residents in community default is to order by trait LMA: trait=
# 'Tr0'. Change to trait='Tr1' to order by Height alternatively you can impose
# an entirely new ordering by specifying new.order

sortResidentsByTrait <- function(Res, trait = "Tr0", decreasing = FALSE, new.order = order(Res$strat[[trait]],
  decreasing = decreasing)) {
  x <- Res
  x$strat <- x$strat[new.order, ]

  for (i in 1:x$nSpp) for (v in c("Age", "Size", "Pop", "Density", "Lambda")) x[[v]][[i]] <- Res[[v]][[new.order[i]]]  # fix here based on order
  x
}

axisInfoFn <- function(parToGet, dim = 2) {
  # returns x-,y-,z- axis ticks, ticklabels, titles args: dim: indicates which
  # trait is evolving, 0 ->lcc, 1 -->hsp, 2--> both

  # lcc
  lccLab <- c(0.01, 0.1, 1, 10, 100)
  lccTick <- log10(lccLab)
  lccTitle <- expression(paste("Leaf mass per unit leaf area (kg ", m^-2, " )"))


  # hsp
  hspLab <- c(0.4, 2, 10, 50)
  hspTick <- log10(hspLab)
  hspTitle <- expression(paste("Height at maturation (m)"))

  # fitness
  fTick <- seq(-6, 6, by = 2)
  fLab <- fTick
  fTitle <- "Fitness"

  out.lcc <- list(ytick = fTick, ylab = fLab, ytitle = fTitle, xtick = lccTick,
    xlab = lccLab, xtitle = lccTitle)

  out.hsp <- list(ytick = fTick, ylab = fLab, ytitle = fTitle, xtick = hspTick,
    xlab = hspLab, xtitle = hspTitle)

  out <- list(ytick = hspTick, ylab = hspLab, ytitle = hspTitle, xtick = lccTick,
    xlab = lccLab, xtitle = lccTitle)
  # lcc evolving plot
  if (dim == 0)
    out <- out.lcc
  # hsp evolving plot
  if (dim == 1)
    out <- out.hsp

  out[[parToGet]]
}

fitnessFileNameFormat <- function(Dir, time, dim, Tsuffix = "") {

  D <- 2
  if (dim %in% c(0, 1))
    D <- 1

  # string format used for fitness files
  paste0(Dir, "/T", Tsuffix, "-", time, "-", D, "D.txt")
}

getFitnessFilesInDir <- function(Dir, from = 200, step = 200, dim = 2, Tsuffix = "") {
  # searches directory Dir for all fitness landscape files


  time <- 0
  if (file.exists(fitnessFileNameFormat(Dir, time, dim, Tsuffix)))
    times <- 0 else times <- NULL

  time <- from
  while (file.exists(fitnessFileNameFormat(Dir, time, dim, Tsuffix))) {
    times <- c(times, time)
    time <- time + step
  }

  # generate pathnames for time sequence
  if (is.null(times))
    paths <- NULL else paths <- fitnessFileNameFormat(Dir, times, dim, Tsuffix)

  list(times = times, paths = paths)
}

loadFitnessDataFromFile <- function(filename = filename, zOffset = 0, zlim = c(-Inf,
  Inf)) {

  Data <- read.table(filename, skip = 1, header = FALSE)

  names(Data) <- switch(as.character(dim(Data)[2]), `2` = c("x", "z"), `3` = c("x",
    "y", "z"))

  # apply offset
  Data$z <- Data$z - zOffset
  # trim within specified limits
  Data$z[Data$z > zlim[2]] <- zlim[2]
  Data$z[Data$z < zlim[1]] <- zlim[1]


  Data
}


plotFitnessLandscapeFromFile <- function(filename, dim = 2, ...) {

  if (!file.exists(filename)) {
    plot.new()  # file not available, plot blank
  } else {

    zlim <- c(-6, 6)

    Data <- try(loadFitnessDataFromFile(filename = filename, zlim = zlim),
      silent = TRUE)
    # use try, because sometimes load fails, if errors in file
    if (inherits(Data, "try-error"))
      plot.new()  # file not available, plot blank
 else {
      Residents <- loadResidnetTraitFromFitnessFile(filename = filename)

      plotFitnessLandscape(Data, Residents, dim = dim, ...)
    }
  }
}


plotFitnessLandscape <- function(Data, Residents, title, dim = 2, Xtick = TRUE,
  Xlab = TRUE, Xtitle = TRUE, Ytick = TRUE, Ylab = TRUE, Ytitle = TRUE, cex = 1,
  outer = FALSE, line = 3, ParseTitle = FALSE, ...) {

  if (dim == 2)
    plotFitnessSurf2D(Data, Residents, ...) else plotFitnessSurf1D(Data, Residents, dim, ...)

  # add title
  if(!is.na(title) && !is.null(title) && title != "" && length(title) > 0){
    if(ParseTitle)
      title <- parse(text=title)
    mtext(title, side = 3, line = 0.5, cex = cex)
  }
  # x axis - add ticks, tick labels, title

  if (Xtick & Xlab)
    axis(1, at = axisInfoFn("xtick", dim), labels = axisInfoFn("xlab", dim),
      las = 1)
  if (Xtick & !Xlab)
    axis(1, at = axisInfoFn("xtick", dim), labels = FALSE, las = 1)
  if (Xtitle)
    mtext(axisInfoFn("xtitle", dim), side = 1, line = line, outer = outer,
      cex = cex)

  # yaxis - add ticks, tick labels, title
  if (Ytick & Ylab)
    axis(2, at = axisInfoFn("ytick", dim = dim), labels = axisInfoFn("ylab",
      dim = dim), las = 1)
  if (Ytick & !Ylab)
    axis(2, at = axisInfoFn("ytick", dim = dim), labels = FALSE, las = 1)
  if (Ytitle)
    mtext(axisInfoFn("ytitle", dim = dim), side = 2, line = line, outer = outer,
      cex = cex)

}

plotFitnessSurf1D <- function(Data, Residents, dim, zlim = c(-6, 6), xlab = "",
  ylab = "", main = "", xlim = par("usr")[1:2], col_pallete_fn = color_2D_pallete,
  ...) {
  # plots trait vs fitness landscape for data stored in filename

  # make plot
  plot(Data$x, Data$z, type = "l", ylim = zlim, xlim = xlim, xlab = xlab, ylab = ylab,
    main = main, xaxs = "i", yaxs = "i", las = 1, ann = FALSE, xaxt = "n",
    yaxt = "n", ...)

  # upper and lower limits of x-axis
  points(xlim, xlim * 0, type = "l", lty = "dashed")

  # add resident location interpolate to get fitness at resident point
  if (nrow(Residents) > 0) {
    cols <- sapply(seq_len(nrow(Residents)), function(i) col_pallete_fn(Residents[i,
      1], Residents[i, 2]))
    points(Residents[, dim + 1], approx(Data$x, Data$z, Residents[, dim + 1])$y,
      type = "p", pch = 21, bg = cols, cex = 1.5)
  }

}


plotFitnessSurf2D <- function(Data, Residents, zlim = c(-6, 6), add = FALSE, xlim = range(Data$x,
  finite = TRUE), ylim = range(Data$y, finite = TRUE), extinctThreshold = 0.02,
  col_pallete_fn = color_2D_pallete) {

  Data <- collapse.grid(Data)

  # load colorMap
  myColorMap <- loadColorMap("R/colormap.csv", zlim = zlim)

  # make plot
  if (!add) {
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, "", xaxs = "i", yaxs = "i", las = 1)
    box()
  }

  .filled.contour(Data$x, Data$y, Data$z, levels = myColorMap$levels, col = myColorMap$col)

  # add resident location
  if (nrow(Residents) > 0) {
    cols <- sapply(seq_len(nrow(Residents)), function(i) col_pallete_fn(Residents[i,
      1], Residents[i, 2]))

    i <- (Residents[, 3] > extinctThreshold)
    points(Residents[i, 1], Residents[i, 2], type = "p", pch = 21, bg = cols[i],
      cex = 1.5)
  }

  return(myColorMap)
  # add residents
}


loadColorMap <- function(filename = "R/colormap.csv", zlim = c(-1, 1)) {
  # load custom colormap read colormap
  tmp <- t(read.csv(filename, h = TRUE))

  # structure to store output
  out <- list()
  out$col <- rgb(tmp[1, ], tmp[2, ], tmp[3, ], alpha = 1)

  # set number of levels to match colormap - should be one more
  nlevels <- length(out$col) + 1
  out$levels <- seq(from = zlim[1], to = zlim[2], length.out = nlevels)
  out
}


panelPlot <- function(datafiles, nrow, ncol, titles = list(), dim = 2, cex = 0.5, omi = c(0.5,
  0.5, 0.1, 0.1), mai = c(0.3, 0.2, 0.2, 0), main = NULL, ...) {

  par(mfrow = c(nrow, ncol), omi = omi, mai = mai)
  for (i in 1:nrow) for (j in 1:ncol) {
    plotFitnessLandscapeFromFile(datafiles[matrixIndexAsVectorIndex(i, j, ncol)],
      titles[[matrixIndexAsVectorIndex(i, j, ncol)]], dim = dim, Xtick = TRUE,
      Xlab = (i == nrow), Xtitle = ((i == nrow) & (j == ceiling(ncol/2))),
      Ytick = TRUE, Ylab = (j == 1), Ytitle = ((j == 1) & (i == ceiling(nrow/2))),
      cex = cex, outer = TRUE, line = 1, ...)
  }
  if (is.null(main))
    mtext(main, side = 3, line = 1, outer = TRUE, cex = cex * 2)
}

plotAssemblyFinalPanel <- function(Dirs, titles = Dirs, from = 200, step = 200,
  ncol = 5, dim = 1, Tsuffix = "", cex = 1, omi = c(0.5, 0.5, 0.1, 0.1), mai = c(0.3,
    0.2, 0.2, 0), main = NULL) {
  if (length(Dirs) > 0) {
    x <- do.call("rbind", lapply(Dirs, function(x) {
      last(getFitnessFilesInDir(x, from = from, step = step, dim = dim, Tsuffix = Tsuffix))
    }))
    panelPlot(x$paths, ceiling(length(x$paths)/ncol), ncol, titles = titles,
      dim = dim, cex = cex, omi = omi, mai = mai, main = main)
  }
}


plotAssembly <- function(Dir, from = 200, step = 200, ncol = 5, dim = 1, Tsuffix = "",
  cex = 1, omi = c(0.5, 0.5, 0.1, 0.1), mai = c(0.3, 0.2, 0.2, 0), main = NULL,
  x = NULL, ...) {
  if (is.null(x))
    x <- getFitnessFilesInDir(Dir, from = from, step = step, dim = dim, Tsuffix = Tsuffix)
  if (length(x) > 0)
    panelPlot(x$paths, ceiling(length(x$paths)/ncol), ncol, titles = paste("Step",
      x$times), dim = dim, cex = cex, omi = omi, mai = mai, main = main,
      ...)
}


loadResidnetTraitFromFitnessFile <- function(filename, ncol = 3, removeExtinct = TRUE,
  extinctThreshold = 0.01) {
  con <- file(filename, "r", blocking = FALSE)
  Res.str <- strsplit(readLines(con, 1), "%")[[1]][2]
  close(con)
  Res <- as.numeric(unlist(strsplit(Res.str, "\t")))

  out <- matrix(Res, ncol = ncol, nrow = length(Res)/3, byrow = TRUE)

  if (ncol == 3 & removeExtinct)
    out <- subset(out, out[, 3] > extinctThreshold)
  out
}


loadStochFile <- function(path, filename = "Stoch.txt") {
  con <- file(file.path(path, filename), "r", blocking = FALSE)
  lines <- readLines(con)
  close(con)

  split <- lapply(lines, function(line) {
    as.numeric(strsplit(line, "\t")[[1]])
  })

  # start at back, so when item set to NULL only affects list elements behind you
  for (i in rev(1:length(split))) if (length(split[[i]]) == 0)
    split[[i]] <- NULL
  # removes dud lines, occur sometimes due to printing errors in file

  split
}


plotSeedRain <- function(res, xlim = c(0, 175), ylim = c(0, 0.05), cols = 1:100,
  ytck = seq(0, 0.05, 0.01), xtck = seq(0, 200, by = 50), xlab = expression(paste("Time since disturbance (y)")),
  ylab = expression(paste("Seed production (", y^-1, ")"))) {

  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", ann = FALSE, xlim = xlim, ylim = ylim,
    xaxs = "i", yaxs = "i")
  axis(1, at = xtck, las = 1)
  axis(2, at = ytck, las = 1)
  mtext(xlab, 1, 3)
  mtext(ylab, 2, 3)

  x <- res$pAge$X.Age
  for (i in seq_len(res$nSpp)) {
    y <- res$Pop[[i]]$Seed_Rain_out
    y_norm <- y/trapz(x, y)
    points(x, y_norm, type = "l", col = cols[i])
  }

}

# rescales log(Density) to fall in range from 0-1, based on limits defined in
# zlim.
scaleDensityForPlot <- function(x, zlim = c(-6, 6), log = TRUE) {
  Density <- log10(x)
  Density[Density < zlim[1]] <- zlim[1]
  Density[Density > zlim[2]] <- zlim[2]
  (Density - zlim[1])/diff(zlim)
}

plotSpeciesDensity2D <- function(Age, Size, Density, col = "red", zlim = c(-6,
  6)) {

  density_scaled <- scaleDensityForPlot(Density, zlim = zlim)

  for (i in 2:(dim(Size)[1] - 1)) {
    j.max <- max(which(!is.na(Size[i, ]))) - 1
    for (j in seq_len(j.max)) {
      x <- Age[c(i, i + 1, i + 1, i, i)]
      y <- as.matrix(Size[i + 0:1, j + 0:1])[c(1, 2, 4, 3, 1)]
      z <- sum(density_scaled[i:(i + 1), j:(j + 1)]) * 0.25  # average density of 4 corners of polygon
      polygon(x, y, col = make.transparent(col, z), border = FALSE)
    }
  }
}

# gives height as a function of leaf mass and LMA
Height <- function(ml, LMA, a1, B1) {
  A <- ml/LMA
  H <- a1 * A^B1
  H
}

dHdM <- function(ml, LMA, a1, B1) {
  a1 * B1/LMA * (ml/LMA)^(B1 - 1)
}


plotSizeDensityDistribution2D <- function(res, cols = 1:10, xlim = c(0, 180), ylim = c(0,
  35), xtck = seq(0, 200, by = 50), ytck = seq(0, 40, by = 10), xlab = expression(paste("Time since disturbance (y)")),
  ylab = expression(paste("Height (", m, ")"))) {

  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", ann = FALSE, xlim = xlim, ylim = ylim,
    xaxs = "i", yaxs = "i")
  axis(1, at = xtck, las = 1)
  axis(2, at = ytck, las = 1)
  mtext(xlab, 1, 3)
  mtext(ylab, 2, 3)

  for (n in 1:res$nSpp) {
    Age <- res$Age[[n]][, 1]
    Size <- Height(res$Size[[n]], LMA = res$strat$Tr0[n], a1 = res$params["p.a1",
      1], B1 = res$params["p.B1", 1])
    Density <- res$Density[[n]]/dHdM(res$Size[[n]], LMA = res$strat$Tr0[1],
      a1 = res$params["p.a1", 1], B1 = res$params["p.B1", 1])
    plotSpeciesDensity2D(Age, Size, Density, col = cols[n])
  }

}

# plots size distributions using same technique as used in plant
plotSizeDensityDistribution2D_plant <- function(res, cols = 1:10, xlim = c(0, 180), ylim = c(0,
  30)) {

  ## Arbitrarily clamp this at xmin and convert to relative 0/1
  rel <- function(x, xmin) {
    x[x < xmin] <- xmin
    xmax <- max(x, na.rm=TRUE)
    (x - xmin) / (xmax - xmin)
  }

  plot(NA, xlim=xlim, ylim=ylim,
       las=1, xlab="Time since disturbance (years)", ylab="Height (m)")

  for (i in 1:res$nSpp) {
    Age <- res$Age[[i]][, 1]
    Size <- Height(res$Size[[i]], LMA = res$strat$Tr0[i], a1 = res$params["p.a1",
      1], B1 = res$params["p.B1", 1])
    Density <- as.matrix(res$Density[[i]]/dHdM(res$Size[[i]], LMA = res$strat$Tr0[1],
      a1 = res$params["p.a1", 1], B1 = res$params["p.B1", 1]))

    density_scaled <- rel(scaleDensityForPlot(Density, zlim = c(-6,6)), 0)
    density_scaled[is.na(density_scaled)] <- 0
    col1 <- matrix(make.transparent(cols[i], density_scaled), nrow(density_scaled))

    n <- length(Age)

    h2 <- Size
    x <- matrix(rep(Age, ncol(h2)), nrow(h2))

    segments(x[-1, ], h2[-1, ], x[-n, ], h2[-n, ], col=col1[-n, ], lend="butt")
  }
}



get.Moment <- function(dat, trait, n = 1) {

  trait.mean <- function(z, trait) {
    sum(z$Traits[, trait] * z$PopSize)/sum(z$PopSize)
  }


  moment.nth <- function(z, trait, n = 2) {
    moment.1 <- trait.mean(z, trait)
    moment.2 <- sum((z$Traits[, trait] - moment.1)^n * z$PopSize)/sum(z$PopSize)
    moment.2
  }

  x <- sapply(dat, function(x) x$time)
  if (n == 1)
    y <- sapply(dat, trait.mean, trait = trait) else if (n > 1)
    y <- sapply(dat, moment.nth, trait = trait, n = n) else stop(paste("bad value given to moment function", n))

  list(x = x, y = y)
}


getEvolvedCommunity <- function(Dir, filename = "Stoch.txt", step = NA) {

  dat <- lapply(loadStochFile(Dir, filename), parseLine)
  if (is.na(step))
    step <- length(dat)
  dat[[step]]
}

parseLine <- function(x, nTraits = 4) {
  nSpecies <- x[2]

  Traits <- PopSize <- NULL
  if (nSpecies > 0) {
    Traits <- matrix(x[2 + 1:(nSpecies * nTraits)], ncol = nTraits, nrow = nSpecies,
      byrow = TRUE)
    PopSize <- x[2 + nSpecies * nTraits + 1:nSpecies]
  }

  list(time = x[1], nSpecies = nSpecies, nTraits = nTraits, Traits = Traits,
    PopSize = PopSize)
}

make.grid <- function(xvec, yvec, value = 0) {
  ncol <- length(xvec) - 1
  nrow <- length(yvec) - 1
  Z <- matrix(value, ncol = ncol, nrow = nrow)
  X <- matrix(xvec[col(Z)], ncol = ncol, nrow = nrow)
  Y <- matrix(yvec[row(Z)], ncol = ncol, nrow = nrow)
  list(X = X, Y = Y, Z = Z)
}

getTraitDensity <- function(community, grid = make.grid(seq(-2, 1, by = 0.05),
  seq(-1, 2, by = 0.05))) {
  Z <- grid$Z * 0
  for (i in seq_len(community$nSpecies)) {
    loc <- locatePointInGrid(community$Traits[i, 1], community$Traits[i, 2],
      grid)
    Z[loc$y, loc$x] <- Z[loc$y, loc$x] + community$PopSize[i]
  }
  Z
}

locatePointInGrid <- function(x, y, grid) {
  list(x = locatePointInVec(x, grid$X[1, ]), y = locatePointInVec(y, grid$Y[,
    1]))
}

locatePointInVec <- function(x, xvec) {
  i <- (x >= xvec[-length(xvec)]) & (x < xvec[-1])
  if (!any(i))
    warning("Trait ", x, " outside grid with range ", paste(range(xvec), collapse = "-"))
  which(i)
}

collectDensity <- function(dirs, grid) {
  Z <- grid$Z * 0
  for (d in dirs) {
    community <- getEvolvedCommunity(d)
    Z <- Z + getTraitDensity(community, grid)
  }
  Z
}

removeExtinct <- function(community, threshold = 0.01) {
  i <- (community$PopSize > threshold)
  community$Traits <- community$Traits[i, , drop = FALSE]
  community$PopSize <- community$PopSize[i]
  community$nSpecies <- sum(i)
  community
}

figure1a_getPatchCoordinates <- function(i, patchSize, patchSpace, y0 = 0) {
  x0 <- i * patchSize[1] + (i - 1) * patchSpace
  list(x = c(x0, x0 + patchSize[1], x0 + patchSize[1], x0), y = c(y0, y0, patchSize[2],
    patchSize[2]))
}


figure1c_getPatchCoordinates <- function(i, patchSize, patchSpace, y0 = 0) {
  ncol <- 4
  ii <- floor((i - 1)/ncol) + 1
  jj <- i - (ii - 1) * ncol
  x0 <- (ii - 1) * (patchSize[1] + patchSpace)
  y0 <- (jj - 1) * (patchSize[2] + patchSpace)
  # offset even numbers
  if (jj%%2 == 0)
    x0 <- x0 + 0.5 * (patchSize[1] + patchSpace)

  data.frame(x = x0 + c(0, patchSize[1], patchSize[1], 0), y = y0 + c(0, 0, patchSize[2],
    patchSize[2]))
}


figure2 <- function() {

  cex <- 1.5

  dir1 <- "output/data/[60,1]"
  dir2 <- "output/finalCommunities"
  res.lcc <- sortResidentsByTrait(loadResident(file.path(dir2, "[60,1]", "lcc",
    "[1.3725]", "base", "T-0")), trait = "Tr0")
  landscape.lcc <- file.path(dir1, "lcc/[1.3725]/base/T0-1000-1D.txt")
  res.hsp <- sortResidentsByTrait(loadResident(file.path(dir2, "[60,1]", "hsp",
    "[0.625]", "base", "T-0")), trait = "Tr1", decreasing = TRUE)
  landscape.hsp <- file.path(dir1, "hsp/[0.625]/base/T1-3500-1D.txt")
  res.both1 <- sortResidentsByTrait(loadResident(file.path(dir2, "[60,1]", "2trait",
    "base", "T-0")), trait = "Tr0")
  res.both <- sortResidentsByTrait(res.both1, new.order = c(1:5, 7, 8, 6, 9:12))
  landscape.both <- file.path(dir1, "2trait/base/T-5000-2D.txt")

  mylabel <- function(label, x, y = 1.2) {
    label(-0.2, y, label, cex = 2, xpd = NA, font = 2)
    label(-0.1, y, x, cex = 2, xpd = NA, font = 1)
  }

  HLIM <- log10(c(0.4, 60))

  m <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), ncol = 3, byrow = TRUE)
  layout(m, widths = c(0.33, 0.33, 0.33), heights = c(0.33, 0.33, 0.33))

  par(oma = rep(2, 4), mar = c(5.1, 4.1, 4, 3.1))
  cols <- getColours(res.lcc)
  xlim <- c(0, 1)
  plotSizeDensityDistribution2D(res.lcc, cols = cols, ylim = c(0, 40), xlab = "")
  mylabel("A", "Only leaf mass per unit leaf area evolving")
  plotSeedRain(res.lcc, cols = cols, xlab = "", ylim = c(0, 0.04), ytck = seq(0,
    0.04, 0.01), )
  plotFitnessLandscapeFromFile(landscape.lcc, dim = 0, title = "", zlim = c(-2,
    2), xlim = c(-2, 1))

  # ------------
  cols <- getColours(res.hsp)
  plotSizeDensityDistribution2D(res.hsp, cols = cols, ylim = c(0, 40), xlab = "")
  mylabel("B", "Only height at maturation evolving")
  plotSeedRain(res.hsp, cols = cols, ylim = c(0, 0.4), ytck = seq(0, 0.4, 0.1),
    xlab = "")
  plotFitnessLandscapeFromFile(landscape.hsp, dim = 1, title = "", zlim = c(-2,
    2), xlim = HLIM)

  # ------------
  cols <- getColours(res.both)
  plotSizeDensityDistribution2D(res.both, cols = cols, ylim = c(0, 40))
  mylabel("C", "Both traits evolving together")
  plotSeedRain(res.both, cols = cols, ylim = c(0, 0.16), ytck = seq(0, 0.16,
    0.04))
  plotFitnessLandscapeFromFile(landscape.both, dim = 2, title = "", ylim = HLIM)

}


figure3 <- function() {

  dir1 <- "output/data/[60,1]"
  dir2 <- "output/finalCommunities"

  path <- c(file.path(dir2, "[60,1]", "lcc", "[1.3725]"), file.path(dir2, "[60,1]",
    "hsp", "[0.625]"), file.path(dir2, "[60,1]", "2trait"))

  viability <- file.path(path[3], "base", "Viability.txt")
  res.both1 <- sortResidentsByTrait(loadResident(file.path(dir2, "[60,1]", "2trait",
    "base", "T-0")), trait = "Tr0")
  res.both <- sortResidentsByTrait(res.both1, new.order = c(1:5, 7, 8, 6, 9:12))
  landscape.both <- file.path(dir1, "2trait/base/T-5000-2D.txt")

  plotMixturesFromPaths <- function(Dirs, grid = make.grid(xvec = seq(-2, 2,
    length.out = 200), yvec = seq(-1, 2, length.out = 200))) {

    grid$Z <- collectDensity(Dirs, grid)
    points(grid$X, grid$Y, cex = (grid$Z > 0.01), pch = 16, col = "black")
    box()
  }

  emptyTraitsPlot <- function(xlim = c(-2.3, 1.3), ylim = log10(c(1, 70)), xlab = axisInfoFn("xtitle",
    2), ylab = axisInfoFn("ytitle", 2)) {

    plot(0, 0, type = "n", xaxt = "n", yaxt = "n", ann = FALSE, xaxs = "i",
      yaxs = "i", xlim = xlim, ylim = ylim)

    axis(1, at = axisInfoFn("xtick", 2), labels = axisInfoFn("xlab", 2), las = 1)
    mtext(xlab, side = 1, line = 3, cex = 0.8)

    axis(2, at = axisInfoFn("ytick", dim = 2), labels = axisInfoFn("ylab",
      2), las = 1)
    mtext(ylab, side = 2, line = 3, cex = 0.8)
  }

  plotViability <- function(viability) {
    dat <- read.table(viability, header = FALSE, sep = "\t")
    usr <- par("usr")
    polygon(usr[c(1, 1, 2, 2, 1)], usr[c(3, 4, 4, 3, 3)], col = "grey", border = NA)
    polygon(dat[, 2], dat[, 3], col = "white", border = NA)
  }

  mylabel <- function(label, x, y = 1.1) {
    label(-0.1, y, label, cex = 1.05, xpd = NA, font = 2)
    label(-0, y, x, cex = 1.05, xpd = NA, font = 1)
  }

  HLIM <- log10(c(0.5, 80))

  par(mfrow = c(3, 1), oma = c(3,3,0,1), mar = c(2, 2, 3, 0))

  emptyTraitsPlot(ylim = HLIM, xlab = "", ylab="")
  plotViability(viability)
  plotMixturesFromPaths(Dirs = file.path(dir1, "lcc", dir(file.path(dir1, "lcc")),
    "base"))
  mylabel("A", "Only leaf mass per unit leaf area evolving")
  box()

  emptyTraitsPlot(ylim = HLIM, xlab = "")
  plotViability(viability)
  plotMixturesFromPaths(Dirs = file.path(dir1, "hsp", dir(file.path(dir1, "hsp")),
    "base"))
  mylabel("B", "Only height at maturation evolving")
  box()

  emptyTraitsPlot(ylim = HLIM, ylab = "")
  plotViability(viability)
  points(log10(res.both$strat$Tr0), log10(res.both$strat$Tr1), type = "p", col = "black",
    pch = 16, cex = 1.5)
  mylabel("C", "Both traits evolving together")
  box()

}


clusterStrategies <- function(community, h=0.05){

  out <-community

  if(community$nSpecies > 1){

    # run clustering algorithm
    hc <- hclust(dist(community$Traits))
    # cut tree at given height
    membership <- cutree(hc,h=h)  # can cut tree at different heights, or specificy

    #average trait values
    out$nSpecies <- length(unique(membership))
    out$Traits <- matrix(0, nrow = out$nSpecies, ncol = ncol(community$Traits))
    out$PopSize <- rep(0, out$nSpecies)

    for(i in seq_len(community$nSpecies)){
      j <- membership[i]
      out$Traits[j,] = out$Traits[j,] + community$Traits[i,]*community$PopSize[i]
      out$PopSize[j] = out$PopSize[j] + community$PopSize[i]
    }

    # Divide by total population size for each of new strategies
    for(j in seq_len(out$nSpecies))
      out$Traits[j,] = out$Traits[j,] / out$PopSize[j]
  }

  out
}


figure_seed_rain <- function() {

  dir1 <- "output/data/[60,1]"
  Dirs <- file.path(dir1, "hsp", dir(file.path(dir1, "hsp")), "base")

  LMA <- HMAT <- PopSize <- rep(NA_real_, length(Dirs))
  for (i in seq_along(Dirs)) {
    # load community, running through clustering routine
    community <- clusterStrategies(removeExtinct(getEvolvedCommunity(Dirs[i])), h=0.1)
    # find tallest species
    h <- community$Traits[,2]
    ii <- which(h == max(h))
    HMAT[i] <- h[ii]
    PopSize[i] <- community$PopSize[ii]
    LMA[i] <-  community$Traits[1,1]
  }
  par(oma = rep(2, 4), mar = c(5.1, 4.1, 4, 3.1))

  plot(LMA, PopSize, ylim = c(0,400), xlim = c(-2, 1), xlab = "", ylab = "",
     xaxs = "i", yaxs = "i", las = 1, ann = FALSE, xaxt = "n", yaxt = "n")

  axis(1, at = axisInfoFn("xtick", 2), labels = axisInfoFn("xlab", 2),las = 1)
  mtext(expression(paste("Leaf mass per unit leaf area (kg ", m^-2, ")")), line = 3, side=1)
  axis(2, at = c(0, 100, 200, 300, 400), las = 1)
  mtext(expression(paste("Seed rain (", m^-2, ")")), line = 3, side=2)
}

add_fitness_colorbar <- function() {

  myColorMap <- loadColorMap("R/colormap.csv", zlim = -c(-6, 6))
  nticks <- 7

  min=min(myColorMap$levels)
  max=max(myColorMap$levels)
  ticks <- seq(min, max, len=nticks)
  scale <- (length(myColorMap$col)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n',
      xaxs='i', yaxs='i',
      xlab='', yaxt='n', ylab='')
  axis(4, ticks, las=1)
  mtext('Fitness', side=3, line=1)

  for (i in 1:(length(myColorMap$col)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=myColorMap$col[i], border=NA)
  }
  box()

  px <- 2.5
  py <- 0.25

  usr <- par("usr")
  x <- usr[1] + px * (usr[2] - usr[1])
  y <- usr[3] + c(0.5 - py, 0.5, 0.5 + py) * (usr[4] - usr[3])
  dy <- py * (usr[4] - usr[3])
  text(x, y, labels = c(" fail ", "neutral", "invade"), xpd = NA, cex = 1.5,
    font = 1)
  for (i in c(-1, 1)) arrows(x, y0 = y[2] + i * 0.2 * dy, y1 = y[2] + i * 0.8 *
    dy, length = 0.08, angle = 30, code = 2, lwd = 2, xpd = NA)
}


figure4 <- function() {

  steps <- c(0, 1, 5, 10, 200, 2000, 3000, 5000)

  V <- read.table("output/data/[60,1]/2trait/base/Viability.txt", sep = "\t")[, 2:4]


  par(oma = c(5, 5, 1, 1), mar = c(1, 1, 3, 1))

  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), 3, 3, byrow = TRUE), c(1, 1,
    1), c(1, 1,1,1))
  # layout.show(7) ## show the regions that have been allocated to each plot

  for (i in 1:8) {

    plotFitnessLandscapeFromFile(sprintf("output/data/[60,1]/2trait/base/T-%d-2D.txt",
      steps[i]), dim = 2, title = "", zlim = c(-6, 6), Ylab = (i %in% c(4)),
    Ytitle = 0, Xtitle = 0, Xlab = (i > 6), line = 3, ylim = log10(c(0.5,
      80)), xlim = c(-2.5, 1.45), col_pallete_fn = function(x, y) {
      "white"
    })
    label(-0, 1.15, toupper(letters[i]), cex = 1.5, xpd = NA, font = 2)
    label(0.1, 1.15, paste("Step", steps[i]), xpd = NA, cex = 1.5)
    # add viability curve
    if (i > 1)
      points(V[, 1], V[, 2], lty = "dashed", col = "white", type = "l")
  }

  # add color bar
  mar <- par("mar")
  mar[c(2, 4)] <- c(2, 12)
  par(mar = mar)
  add_fitness_colorbar()

  mtext(axisInfoFn("xtitle"), 1, line = 3, outer = TRUE)
  mtext(axisInfoFn("ytitle"), 2, line = 3, outer = TRUE)
}


figS_Assembly_hsp <- function(path,
            times = c(0:20, seq(40, 200, by = 20),
                      seq(200, 1000, by = 200)),
            xlim= c(-1, 2),
            zlim=c(-4,4), ...){
  x <- list(times = times, paths = fitnessFileNameFormat(path,
    times, dim = 1, Tsuffix = "1"))
  plotAssembly("", ncol = 5, dim = 1, xlim = xlim, x = x, zlim=zlim, col_pallete_fn = function(x,
    y) {
    "white"
  }, ...)
}

figS_Assembly_hsp_C <- function(...) {
  figS_Assembly_hsp(..., omi = c(0.5, 0.75, 0.75, 0.1))
  mtext("C", 3 , at=-0.05, cex = 1.35, adj=0, line = 3, xpd = NA, font = 2, outer =TRUE)
  mtext("Fitness landscapes", 3 , at = 0, cex = 1.35, adj=0, line = 3,xpd = NA, font = 1, outer =TRUE)
}

figS_Assembly_lcc <- function(path,
            times = c(0:20, seq(40, 200, by = 20),
                    seq(200, 1000, by = 200)), ...){
  x <- list(times = times, paths = fitnessFileNameFormat(path,
    times, dim = 0, Tsuffix = "0"))
  plotAssembly("", ncol = 5, dim = 0, xlim = c(-2, 1), x = x, col_pallete_fn = function(x,
    y) {
    "white"
  }, ...)
}

figS_Assembly_lcc_C <- function(...) {
  figS_Assembly_lcc(..., omi = c(0.5, 0.75, 0.75, 0.1))
  mtext("C", 3 , at=-0.05, cex = 1.35, adj=0, line = 3, xpd = NA, font = 2, outer =TRUE)
  mtext("Fitness landscapes", 3 , at = 0, cex = 1.35, adj=0, line = 3,xpd = NA, font = 1, outer =TRUE)
}

# vary parameters by 10%
figure_sensitivity_1 <- function() {

  T <- c("a1", "a3", "a4", "c_acc", "c_d0", "c_d2", "c_d3", "c_ext", "c_Rl",
    "c_Rr", "c_Rs", "pi_0", "seed", "theta", "wood_dens")
  labels <- c("a[1]", "a[3]", "a[4]", "c[acc]", "d[I]", "c[d2]", "c[d3]", "c[ext]", "r[Rl0]",
    "r[r]", "r[s]", "S[D]", "s", "theta", "rho")

  P <- c(0.9, 1.1)

  x <- generatePaths("output/data/[60,1]/2trait", T, P, "[%s,%s]")

 titles <- x$titles
 for(i in seq_along(T))
    titles <- gsub(T[i], labels[i], titles, fixed=TRUE)
 titles <- gsub(" ", "~", titles, fixed=TRUE)

  panelPlot(x$paths, nrow = 5, ncol = 6, gsub(" (5000)", "", titles, fixed = TRUE),
    ParseTitle = TRUE,
    dim = 2, cex = 1.2, omi = c(1, 1, 0.5, 0.5), col_pallete_fn = function(x,
      y) {
      "white"
    })
}

# check sensitivity to reproductive allocation
figure_sensitivity_2 <- function() {

  T <- c("c_r1", "c_r2")
  labels <- c("c[r1]", "c[r2]")

  P <- c(0.2, 0.4, 0.6, 0.8, 1)

  x <- generatePaths("output/data/[60,1]/2trait", T, P, "[%s,%s]")
  x$paths[5] <- x$paths[10] <- file.path("output/data/[60,1]/2trait", "base",
    "T-5000-2D.txt")
  x$titles[5] <- "1.0 c_r1"
  x$titles[10] <- "1.0 c_r2"

 titles <- x$titles
 for(i in seq_along(T))
    titles <- gsub(T[i], labels[i], titles, fixed=TRUE)
 titles <- gsub(" ", "~", titles, fixed=TRUE)

  panelPlot(x$paths, nrow = length(T), ncol = length(P), titles, dim = 2, cex = 1.2,
    ParseTitle = TRUE,
    col_pallete_fn = function(x, y) {
      "white"
    })
}


figure_sensitivity_3 <- function() {

  dirs <- c("[c_r1,0.2]","[c_r1,0.4]","[c_r1,0.6]","[c_r1,0.8]","base")

  r1 <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  last_step <- c(860, 700, 720, 1000, 1000)

  paths <- file.path("output/data/[60,1]/hsp/[0.54]", dirs, paste0("T1-", last_step, "-1D.txt"))

  HLIM <- log10(c(0.4, 60))
  par(mfrow=c(3,2), oma = c(4,4,1,1), mar = c(1, 2, 4, 1))
  for(i in seq_along(paths)) {
    plotFitnessLandscapeFromFile(paths[i], dim = 1, title = "", 
      zlim = c(-2,2), xlim = log10(c(0.4, 60)),
      Ylab = i %in% c(1,3,5), Ytitle = i==3, Xtitle = i>=4, 
      Xlab = i >=4, line = 3)

    label(-0.10, 1.15, toupper(letters[i]), cex = 1.5, xpd=NA, font=2)
    label(0.5, 1.15, parse(text=paste0(r1[i],"~c[r1]")), cex = 1.5, xpd=NA)
  }
}


getCellData <- function(data, fun, trait, steps = 1000) {

  processMomentData <- function(dat, fun, trait, steps) {
    x <- 10^tail(dat[[trait]]$y, steps)

    out <- NA
    if (length(x) > 1)
      out <- fun(x)
    out
  }
  vars <- data$vars
  vars$z <- sapply(data$out, processMomentData, fun = fun, trait = trait, steps = steps)

  names(vars) <- c("y", "x", "z")
  vars
}

color.range <- function(start = 0, end = 1, gamma = 1, alpha = 1, rgb = c(0, 1,
  0), n = 100) {
  x <- seq.int(from = start^gamma, to = end^gamma, length.out = n)^(1/gamma)
  rgb(x * rgb[1], x * rgb[2], x * rgb[3], alpha = alpha, names = NULL, maxColorValue = 1)
}

plotCells <- function(Data, zlim = NULL, xlabels = FALSE, ylabels = FALSE, col = rev(color.range(start = 0.3))) {

  # fit bivariate polynomial regression to data Used to interpolate to higher
  # resolution along smoothed surface
  i <- !is.na(Data$z)
  polyfit <- lm(z ~ polym(log(y), x, degree = 5, raw = TRUE), data = Data[i,
    ])

  # interpolate data onto finer grid using fitted polynomial
  data_hires <- expand.grid(y = seq(min(Data$y), max(Data$y), length.out = 100),
    x = seq(0.6, 1.4, length.out = 100))
  suppressWarnings({
    data_hires$z <- predict(polyfit, newdata = data_hires)
  })
  data_image <- collapse.grid(data_hires)
  if (is.null(zlim))
    zlim <- range(data_hires$z)
  with(data_image, image(x, y, z, zlim = zlim, xaxs = "i", yaxs = "i", las = 1,
    log = "y", main = "", xaxt = "n", yaxt = "n", col = col, ylim = c(240,
      7.5), xlim = range(x, na.rm = TRUE), asp = 1))
  axis(1, at = seq(0.6, 1.4, by = 0.2), labels = xlabels)
  axis(2, at = c(10, 20, 50, 100, 200), labels = ylabels, las = 1)
  box()

  # fit contours through interpolated data
  ccs <- with(data_image, contourLines(x, y, z, nlevels = 5))
  for (dat in ccs) {
    points(dat$x, dat$y, type = "l", col = "black")
  }

}

figure5_data <- function(f = lapply) {

  # averages across whole data set
  d <- c(7.5, 10.6, 15, 21.2, 30, 42.4, 60, 84.8, 120, 169.7, 240)
  p <- c(0.6, 0.8, 1, 1.2, 1.4)
  vars <- expand.grid(d, p)
  dirs <- paste0("output/data/[", vars[, 1], ",", vars[, 2], "]/2trait/base")

  # get first and second moments from data for both traits
  out <- f(seq_len(length(dirs)), function(i) {
    cat(paste(i, " "))
    dat <- tail(lapply(loadStochFile(dirs[i]), parseLine), 500)
    list(lcc.1 = get.Moment(dat, trait = 1, n = 1), hsp.1 = get.Moment(dat,
      trait = 2, n = 1))
  })

  list(dirs = dirs, vars = vars, out = out)
}

figure5 <- function(data) {

  mylabel <- function(label, x, y = 1.15) {
    label(-0.15, y, label, cex = 1.5, xpd = NA, font = 2)
    label(-0.05, y, x, cex = 1.5, xpd = NA, font = 1)
  }

  # specific example communities
  x <- c(1, 0.6, 0.8, 0.8, 1, 1.4)
  y <- c(10.6, 60, 120, 21.2, 60, 120)
  labels <- c("C", "D", "E", "F", "G", "H")

  dirs <- paste0("output/finalCommunities/[", y, ",", x, "]/2trait/base/T-0")
  nspecies <- sapply(seq_len(length(dirs)), function(i) {
    nrow(read.table(file.path(dirs[i], "Strategy.txt"), head = TRUE))
  })

  titles <- paste(c("Riverbank,", "Arid shrubland,", "Temperate forest,", "Woodland,",
    "Tropical rain forest,", "Tall tropical rain forest,"), "N =", nspecies)

  col1 <- rev(color.range(start = 0, end = 0.7, gamma = 1, rgb = c(0, 1, 0)))
  col2 <- rev(color.range(start = 0, end = 1, gamma = 1, rgb = c(217, 95, 14)/256))

  layout(matrix(c(1, 2, 3, 4, 5, 5, 6, 6), ncol = 4, byrow = FALSE), heights = c(0.5,
    0.5), widths = c(0.25, 0.05, 0.35, 0.35))
  par(oma = c(6, 6, 6, 6), mar = c(1, 1, 5, 1))

  dat <- getCellData(data, fun = mean, trait = "lcc.1")
  zlim.1 <- range(dat$z, na.rm = TRUE)
  plotCells(dat, ylabels = TRUE, col = col1)
  points(x, y, pch = 16, cex = 3, xpd = NA, col = "grey20")
  text(x, y, labels, col = "white", xpd = NA)
  polygon(c(0.6, 0.7, 0.7, 0.6, 0.6), c(7.5, 7.5, 12, 12, 7.5), col = "white",
    border = NA)
  text(0.65, 9, labels = "Not\nviable", cex = 0.8)
  mylabel("A", "Average leaf mass per area")

  dat <- getCellData(data, fun = mean, trait = "hsp.1")
  zlim.2 <- range(dat$z, na.rm = TRUE)
  plotCells(dat, ylabels = TRUE, xlabels = TRUE, col = col2)
  points(x, y, pch = 16, cex = 3, xpd = NA, col = "grey20")
  text(x, y, labels, col = "white", xpd = NA)

  polygon(c(0.6, 0.7, 0.7, 0.6, 0.6), c(7.5, 7.5, 12, 12, 7.5), col = "white",
    border = NA)
  text(0.65, 9, labels = "Not\nviable", cex = 0.8)
  mylabel("B", "Average height at maturation")

  mtext("Site-productivity index", 1, line = 4, cex = 1.25)
  mtext("                                                     Time since disturbance (y)",
    side = 1, line = 3, outer = TRUE, cex = 1.25)
  mtext("Average disturbance interval (y)", side = 2, line = 3, outer = TRUE,
    cex = 1.25)

  par(mar = c(1, 1, 5, 2))
  add.colorbar(col1, zlim.1)
  label(1.7, 1.05, bquote("kg" ~ m^-2 ~ ""), cex = 0, xpd = NA, font = 1)
  add.colorbar(col2, zlim.2)
  label(1.7, 1.05, bquote("m"), cex = 0, xpd = NA, font = 1)

  par(mar = c(1, 8, 5, 1))
  plot(0, 1, type = "n", frame.plot = FALSE, axes = FALSE, xlim = c(0, 1), ann = FALSE)
  mylabel("C", titles[1], y = 1.05)
  mylabel("D", titles[2], y = 0.7)
  mylabel("E", titles[3], y = 0.32)

  plot(0, 1, type = "n", frame.plot = FALSE, axes = FALSE, xlim = c(0, 1), ann = FALSE)
  mylabel("F", titles[4], y = 1.05)
  mylabel("G", titles[5], y = 0.7)
  mylabel("H", titles[6], y = 0.32)

}


add.colorbar <- function(col, limits) {

  plot(c(0, 10), limits, type = "n", bty = "n", xaxt = "n", xaxs = "i", yaxs = "i",
    xlab = "", yaxt = "n", ylab = "")
  y <- pretty(limits)
  axis(4, at = y, labels = y, las = 1)

  scale <- (length(col) - 1)/diff(limits)
  for (i in 1:(length(col) - 1)) {
    y <- (i - 1)/scale + limits[1]
    rect(0, y, 10, y + 1/scale, col = col[i], border = NA)
  }
  box()
}


figure5_single <- function(data) {

  layout(matrix(c(1, 2), ncol = 2, byrow = FALSE), heights = c(0.5,
    0.5), widths = c(0.25, 0.05))

  col1 <- rev(color.range(start = 0, end = 0.7, gamma = 1, rgb = c(0, 1, 0)))

  par(oma = c(4, 4, 0, 4), mar = c(1, 1, 1, 1))
  dat <- getCellData(data, fun = mean, trait = "lcc.1")
  zlim.1 <- range(dat$z, na.rm = TRUE)
  plotCells(dat,  ylabels = TRUE, xlabels = TRUE, col = col1)
  polygon(c(0.6, 0.7, 0.7, 0.6, 0.6), c(7.5, 7.5, 12, 12, 7.5), col = "white",
    border = NA)
  # text(0.65, 9, labels = "Not\nviable", cex = 0.8)
  mtext("Site-productivity index", 1, line = 3, cex = 2)
  mtext("Average disturbance interval (y)", side = 2, line = 3, cex = 2)

  par(mar = c(1, 1, 5, 2))
  add.colorbar(col1, zlim.1)
  label(1.7, 1.05, bquote("kg" ~ m^-2 ~ ""), cex = 0, xpd = NA, font = 1)
}


figure_get_example_patch <- function(d, prod, scale) {
  path <- paste0("output/finalCommunities/[", d, ",", prod, "]/2trait/base/T-0")

  Res <- sortResidentsByTrait(loadResident(path), trait = "Tr0")

  # scale patches for nice viewing. The approach taken here was to sample a
  # larger area from communities with shorter disturbance intervals.

  nPatches <- 50
  # size of patch with d 120, the longest interval considered
  max_size <- c(20/nPatches, 5)
  max_area <- prod(max_size) * nPatches

  # area of this patch, in proportion to disturbance interval
  true_area <- max_area * d/120
  true_size <- c(true_area/5/nPatches, 5)

  # determine how much want to up scale area by. note we keep relative sizing of
  # patches by up-scaling to area slightly less than max area (scale)
  area_scaling <- max_area * scale/true_area
  scaled_size <- max_size * sqrt(area_scaling)
  scaled_area <- prod(scaled_size) * nPatches
  vadj <- sqrt(scaled_area/max_area)
  # used to adjust lines lengths in plot

  # generate data
  ages <- 10^seq(log10(0.5), log10(3 * d), length.out = nPatches)
  spatial <- getAgeSequence(Res, ages, patchSize = scaled_size, patchSpace = 0,
    scaling = c(1, 1), col = adjustHSV(getColours(Res), c(0, 0, 0.05)))

  list(spatial = spatial, ages = ages, scaled_size = scaled_size, d = d, prod = prod,
    scale = scale, vadj = sqrt(scaled_area/max_area)  # used to adjust lines lengths in plot
)
}

figure_plot_patch <- function(data) {

  plot_stand(data$spatial, size = c(1, 1, 2000, 2000))

  nPatches <- length(data$spatial$coords)

  # add box around plot
  coords <- data$spatial$coords[[nPatches]]
  coords$x[c(1, 4)] <- data$spatial$coords[[1]]$x[c(1, 4)]
  plotPatchOutline(coords)

  # add scale bar in top left corner, 5m tall
  y0 <- data$scaled_size[2]
  add_height_bar(x=coords$x[1], y=data$scaled_size[2],vadj = data$vadj )

  # add labels along base - rescale mean disturbance onto new axis
  i_di <- which(abs(data$ages - data$d) == min(abs(data$ages - data$d)))
  lab <- c(data$ages[1], as.numeric(format(c(data$d, data$ages[length(data$ages)]),
    digits = 0)))
  patch_index <- c(1, i_di, length(data$ages))
  for (ii in 1:3) {
    t <- mean(data$spatial$coords[[patch_index[ii]]]$x)
    lines3d(rep(t, 2), c(0, -1) * data$vadj, c(0, 0), col = "black")
    text3d(t, -3 * data$vadj, 0, lab[ii], adj = c(0.5, 0), col = "black")
  }
}

add_height_bar <- function(x,y,vadj=1, H = 5, lab="5m") {
  lines3d(rep(x, 2), c(y, y), c(0, H), col = "red")
  lines3d(rep(x, 2) - c(0, 0.5) * vadj, c(y, y), c(0, 0), col = "red")
  lines3d(rep(x, 2) - c(0, 0.5) * vadj, c(y, y), c(H, H), col = "red")
  text3d(x - 1 * vadj, y, H, lab, adj = c(0.75, 0), col = "black",
    cex = 0.75)
}

fig_2D_landscapes_all <- function() {

  T <- c(7.5, 10.6, 15, 21.2, 30, 42.4, 60, 84.8, 120, 169.7, 240)
  P <- c(0.6, 0.8, 1, 1.2, 1.4)
  x <- generatePaths("output/data", T, P, "[%s,%s]/2trait/base")

  dim <- 2
  cex <- 1

  nrow <- length(T)
  ncol <- length(P)

  par(mfrow = c(nrow, ncol), oma = c(6, 6, 6, 7), mar = c(0.5, 0.5, 0, 0))
  for (i in 1:nrow) for (j in 1:ncol) {
    plotFitnessLandscapeFromFile(x$paths[matrixIndexAsVectorIndex(i, j, ncol)],
      "", dim = dim, Xtick = TRUE, Xlab = (i == nrow), Xtitle = ((i == nrow) &
        (j == ceiling(ncol/2))), Ytick = TRUE, Ylab = (j == 1), Ytitle = ((j ==
        1) & (i == ceiling(nrow/2))), cex = cex, outer = TRUE, line = 4,
      ylim = log10(c(1, 85)), col_pallete_fn = function(x, y) {
        "white"
      })
    if (i == 1)
      mtext(P[j], side = 3, line = 1, col = "darkgreen")
    if (j == ncol)
      mtext(T[i], side = 4, line = 1, las = 1, col = "darkgreen")
  }
  text(4, 12, srt = -90, pos = 1, labels = "Disturbance interval (y)", xpd = NA,
    cex = 1.6, col = "darkgreen")
  mtext("Site-productivity index", 3, line = 4, outer = TRUE, col = "darkgreen")
}


figure1_data <- function(path = "output/finalCommunities/[60,1]/2trait/base/T-0") {
  Res <- sortResidentsByTrait(loadResident(path),
    trait = "Tr0")
  Res <- sortResidentsByTrait(Res, new.order = c(1:5, 7, 8, 6, 9:12))

  list(Res = Res, col = adjustHSV(getColours(Res), c(0, 0, -0.05)), col.hi = adjustHSV(getColours(Res),
    c(0, 0, 0.05)), ages = 10^seq(log10(1), log10(200), length.out = 50))
}

figure_1b_data <- function(all) {
  nPatches <- 16
  ages <- drawRandomAges(all$Res$pAge, nPatches)
  getAgeSequence(all$Res, ages, patchSize = c(3, 3), patchSpace = 4, patchCoordinatesFn = figure1c_getPatchCoordinates,
    rotateBy = pi/4, col = all$col.hi)
}

figure1_data_single <- function(path = "output/finalCommunities/[60,1]/2trait/single/T-0") {
  Res <- loadResident(path)

  list(Res = Res, col = adjustHSV(getColours(Res), c(0, 0, -0.05)), col.hi = adjustHSV(getColours(Res),
    c(0, 0, 0.05)), ages = 10^seq(log10(1), log10(200), length.out = 50))
}

figure1_data_lcc <- function(path = "output/finalCommunities/[60,1]/lcc/[1.3725]/base/T-0") {
  Res <- loadResident(path)

  list(Res = Res, col = adjustHSV(getColours(Res), c(0, 0, -0.05)), col.hi = adjustHSV(getColours(Res),
    c(0, 0, 0.05)), ages = 10^seq(log10(1), log10(200), length.out = 50))
}

figure1_data_hsp <- function(path = "output/finalCommunities/[60,1]/hsp/[0.625]/base/T-0") {
  Res <- loadResident(path)

  # subsample --> order by abundance and take take top 5 species
  i <- order(Res$strat[, "Tr1"])
  i <- order(Res$strat[, "X_end"], decreasing=TRUE)[1:5]

  x <- Res
  x$strat <- Res$strat[i, ]
  x$nSpp <- length(i)
  for (v in c("Age", "Size", "Pop", "Density", "Lambda")) {
    x[[v]] <- list()
    for(j in seq_along(i)) {
      x[[v]][[j]] <- Res[[v]][[i[j]]]
    }
  }

  list(Res = x, col = adjustHSV(getColours(x), c(0, 0, -0.05)), col.hi = adjustHSV(getColours(x),
    c(0, 0, 0.05)), ages = 10^seq(log10(1), log10(200), length.out = 50))
}

# gather LAI data - added cumulatively
extract_LAI_Data <- function(Res) {
  nSpp <- Res$nSpp
  LAI <- Res$Pop[[1]][, c("X.Age", "LAI")]
  if (nSpp > 1)
    for (i in 2:nSpp) LAI <- cbind(LAI, LAI[, i] + Res$Pop[[i]][, "LAI"])
  LAI[, 1 + 1:nSpp] <- LAI[, 1 + 1:nSpp]
  names(LAI) <- c("Age", as.character(1:nSpp))
  LAI
}


figure1_LAI <- function(Res, cols = defaultColours()) {

  LAI <- extract_LAI_Data(Res)

  mylabel <- function(label, x, y = 1.2) {
    label(-0.2, y, label, cex = 2, xpd = NA, font = 2)
    label(-0.1, y, x, cex = 2, xpd = NA, font = 1)
  }

  m <- matrix(c(1, 1, 2, 3), ncol = 2, byrow = FALSE)
  layout(m, widths = c(0.5, 0.5), heights = c(0.5, 0.5))

  par(oma = c(4, 2, 4, 1))
  par(mar = c(2, 8, 18, 8))
  ymax <- 5
  xlim <- c(1, 200)

  plot(NA, type = "n", xaxt = "n", yaxt = "n", ann = FALSE, xaxs = "i", yaxs = "i",
    xlim = xlim, ylim = c(0, ymax), log = "x")
  mylabel("A", "Temporal dynamics within a patch", 2.225)


  xtick <- c(1, 5, 10, 20, 40, 60, 120, 180)
  axis(1, at = xtick, labels = xtick, las = 1)
  axis(3, at = xtick, labels = NA, line = 1.5, outer = FALSE, tcl = -par("tcl"))

  mtext("Time since disturbance (y)", side = 1, cex = 1, line = 3)
  axis(2, at = 0:5, labels = 0:5, las = 1)
  mtext("Leaf area per ground area", side = 2, cex = 1, line = 3)
  axis(4, at = c(0, 0.25, 0.5, 0.75, 1) * ymax, labels = c(0, 0.25, 0.5, 0.75, 1),
    las = 1)
  text(500, 2.5, srt = -90, pos = 1, labels = "Probability patch\nremains undisturbed", xpd = TRUE,
    cex = 1.3)

  for (i in (Res$nSpp + 1):2) polygon(c(LAI[, 1], LAI[nrow(LAI), 1]), c(LAI[,
    i], 0), col = cols[i - 1], border = NA)
  points(Res$pAge[, 1], (Res$pAge[, 2]/Res$pAge[1, 2] * ymax), type = "l", lwd = 2)


  par(mar = c(2, 8, 2, 8))
  plot(1, 1, type = "n", frame.plot = FALSE, axes = FALSE, ann = FALSE)
  mylabel("B", "Spatial distribution of patches", 1.1)

  plot_color_pallette(default_colour_pallete(1))
  points(log10(Res$strat$Tr0), log10(Res$strat$Tr1), pch = 21, bg = "white",
    cex = 1.5)
  mylabel("C", "Evolved trait mixture")

  mtext(axisInfoFn("xtitle"), side = 1, cex = 1, line = 3)
  mtext(axisInfoFn("ytitle"), side = 2, cex = 1, line = 3)
}

plot_color_pallette <- function(space = default_colour_pallete(), xlim = c(-1.5,
  0), ylim = c(0.4, 1.6)) {
  xo <- seq(xlim[1], xlim[2], length.out = 100)
  yo <- seq(ylim[1], ylim[2], , length.out = 100)
  xy <- expand.grid(xo, yo)
  plot(xy[, 1], xy[, 2], col = color_2D_pallete(xo, yo, space), pch = 15, cex = 3,
    xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", ann = FALSE, xaxt = "n",
    yaxt = "n")

  x.tick <- 0.01 * 2^(1:10)
  y.tick <- 1 * 2^(1:10)
  axis(1, labels = x.tick, at = log10(x.tick), las = 1)
  axis(2, labels = y.tick, at = log10(y.tick), las = 1)
}

BCIplot <- function() {
  dat <- read.table("~/Dropbox/_research/active/Falster-CompnForLight/data/Wright2010.txt",
    header = TRUE, sep = "\t")
  par(oma = c(2, 2, 1, 1))
  plot(dat$MRT25SAP, dat$RGR95SAP, pch = 16, col = make.transparent("black",
    0.5), ann = FALSE, log = "", lwd = 0, cex = 1, las = 1)
  sm <- sma(dat$RGR95SAP ~ dat$MRT25SAP, log = "")
  plot(sm, type = "l", col = "darkgreen", lwd = 2, add = TRUE)
  mtext("Maximum relative growth rate of saplings", 2, outer = FALSE, line = 3)
  mtext("Mortality rate of slow-growing individuals", 1, outer = FALSE, line = 3)
  legend("topleft", legend = paste0("Site: BCI Panama\nn = ", sm$n[[1]], " species\nr2 = ",
    format(sm$r2[[1]], digits = 2)), bty = "n", cex = 0.8)
}


load_R_growth_model <- function() {
  e <- new.env()
  source("R/growthModel-params.R", local = e)
  source("R/growthModel.R", local = e)
  e
}

blank_plot <- function(xlim = c(0.005, 2), ylim = c(0.04, 15), log = "xy") {
  plot(1, 1, type = "n", xlim = xlim, ylim = ylim, log = log, ann = FALSE, las = 1,
    xaxs = "i", yaxs = "i")
}

LMA_plot <- function() {

  # data downloaded from http://dx.doi.org/10.1038/nature02403
  dat.Wright <- process_wright_2004()

  G <- dat.Wright[["dataset"]]
  kl <-  dat.Wright[["k_l"]]
  lma <- dat.Wright[["lma"]]

  # only use sites with n > 5
  site.n <- data.frame(table(G))
  i <- !is.na(lma) & !is.na(kl) & site.n[match(G, as.character(site.n[, 1])), 2] > 5

  sm1 <- sma(kl[i] ~ lma[i] * G[i], log = "xy")
  blank_plot(ylim = c(0.01, 60))
  points(lma, kl, col = make.transparent("grey", 0.6), pch = 16, cex = 0.9,
    type = "p")
  plot(sm1, add = T, col = "darkgreen", type = "l", lwd = 1, p.lines.transparent = 0.15)

  x <- 10^seq(log10(0.001), log10(3), by = 0.01)
  points(x, 0.0286 * x^-1.71, type = "l", col = "black", lwd = 2)

  mtext(expression(paste("Leaf mass per unit leaf area (kg ", m^-2, ")")), line = 3,
    side = 1, cex = 1)
  mtext(expression(paste("Leaf turnover rate (", y^-1, ")")), line = 3, side = 2,
    cex = 1)
}


growth_plot <- function() {

  e <- load_R_growth_model()

  # baseline traits
  traits <- list()
  traits$lma <- 10^-0.3
  traits$rho <- 608
  traits$hmat <- 20

  # returns vector from lo to hi with multiplication steps of incr.
  seqLog <- function(from, to, n, base = 10) {
    base^(seq(log(from, base), log(to, base), length.out = n))
  }

  x <- seqLog(0.005, 2, 100)
  y <- e$getLMAvHeightGrowth(x, traits, h = 0.25, env = 1)
  y2 <- -log(e$getLMAvMaxShading(x, traits, h = 0.25))/0.5

  ymax <- max(y) * 1.25
  blank_plot(log = "x", xlim = range(x), ylim = c(0, ymax))
  points(x, y, type = "l", col = "red")
  points(x, y2/5 * ymax, type = "l", col = "darkblue")

  axis(4, at = seq(0, ymax, by = ymax/5), labels = 0:5, las = 1)

  mtext(expression(paste("Leaf mass per unit leaf area (kg ", m^-2, ")")), line = 3,
    side = 1, cex = 1)
  mtext(expression(paste("Height growth rate (m ", y^-1, ")")), line = 3, side = 2,
    cex = 1, col = "red")
  text(8, 0.4, srt = -90, pos = 1, labels = expression(paste("Shade tolerance (",
    m^2 ~ ~m^-2, ")")), xpd = NA, cex = 1.1, col = "darkblue")

}

plotleaf <- function() {

  mylabel <- function(x) label(-0.2, 1.15, x, adj = c(0, 1), xpd = NA, cex = 1,
    font = 2)

  par(mfrow = c(1, 2), oma = c(2, 2, 2, 6))
  LMA_plot()
  mylabel("A")
  growth_plot()
  mylabel("B")
}


figure_ReproductiveAllocation <- function() {

  mylabel <- function(x) label(-0.2, 1.2, x, adj = c(0, 1), xpd = NA, cex = 1,
    font = 2)

  par(mfrow = c(1, 2), oma = c(2, 2, 2, 3))

  e <- load_R_growth_model()

  # p.c_r1 <- 1 p.c_r2 <- 50

  x <- seq(0, 2, by = 0.01)
  y <- e$ReproductiveAllocation(1, x)

  plot(x, y, ylim = c(0, 1.1), ann = FALSE, las = 1, xaxs = "i", yaxs = "i",
    xaxt = "n", type = "l", col = "red", lwd = 2)

  abline(v = 1, lty = "dashed")
  points(1, 0.5, pch = 16)

  axis(1, at = c(0, 1), labels = c(0, expression(paste(H[m]))))
  mtext("Plant height (m)", line = 3, side = 1, cex = 1)
  mtext("Fraction of net dry mass\nallocated to reproduction", line = 3, side = 2,
    cex = 1)
  mylabel("A")
  plot(0, 0, ann = FALSE, xaxt = "n", yaxt = "n", type = "n", axes = FALSE)
  mylabel("B")
}

figure_competition <- function() {

  par(mfrow=c(1,3))
  path <- "output/finalCommunities/[60,1]/2trait/base/T-0"
  res <- sortResidentsByTrait(loadResident(path), trait = "Tr0")
  cols <- getColours(res)

  eta <- res$params["p.eta",1]

  leaf.pdf <- function(z, h) {
    tmp <- (z / h)^eta
    ret <- 2 * eta * (1 - tmp) * tmp / z
    ret[z>h] <- 0
    ret
  }

  leaf.cdf <- function(x, h) {
    ret <- ((1 - sqrt(x))^(1/eta)) * h
    ret[z>h] <- 0
    ret
  }

  mylabel <- function(x) label(-0.2, 1.05, x, adj = c(0, 1), xpd = NA, cex = 1,
    font = 2)

  z <- seq(0, 1.1, length.out=500)
  plot(leaf.pdf(z, 1), z, type="n", las=1, axes=FALSE, xlab=expression(paste("Probability density of leaf area, q (",
    m^-1, ")")), ylab="Height, z (m)")
  polygon(leaf.pdf(z, 1), z, col="black")
  box()
  axis(2, at=c(0,1), labels = c(0, "H"), las=1)
  axis(1, at=0)
  mylabel("A")

  plot(leaf.cdf(z, 1), z, type="l", las=1, axes=FALSE, xlab=expression(paste("Fraction of leaf area, Q")), ylab="", xlim=c(0, 1.1))
  box()
  axis(2, at=c(0,1), labels = c(0, "H"), las=1)
  axis(1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), las=1)
  mylabel("B")


  Env <- read.table(file.path(path, "patch_age.txt"), head = TRUE)
  H_I <- grepl("H", names(Env))
  E_I <- grepl("E", names(Env))

  Ages <- c(2, 5, 10, 20, 30, 60, 100, 180)
  i <- seq_len(length(Env[,1])-1)
  rows <- sapply(Ages, function(a) which(Env[i,1] < a & Env[i+1,1] > a))
  plot(NA, type="l", xlim=c(0,1), ylim=c(0, 30), xlab = "Canopy openness, E", ylab = "", las=1)
  for(row in rows) {
    E <- as.numeric(Env[row, E_I])
    H <- as.numeric(Env[row, H_I])
    lines(E, H, col="black")
    lines(c(1,1), c(max(H), 35), col="black")
  }
  row <- rows[7]
  E <- as.numeric(Env[row, E_I])
  H <- as.numeric(Env[row, H_I])
  lines(E,H)

  lines(c(1,1), c(max(H), 35))
  mylabel("C")

}

figure_evolution_2_traits <- function(path){

  dat <- lapply(loadStochFile(path), parseLine)

  par(mfrow=c(3,2), oma = c(4,3,2,3), mar= c(3,5,2,1))

  landscape <- getLastLandscape(path)

  mylabel <- function(x) label(-0.15, 1.05, x, adj = c(0, 1), xpd = NA, cex = 1,
    font = 2)

  panel_evolution_raw(dat, step =10, trait=1, xlab = "", cex=0.5)
  mtext(axisInfoFn("xtitle",0), 3, line=2)
  mylabel("A")

  panel_evolution_raw(dat, step =10, trait=2, xlab = "", ylab = "", cex=0.75)
  mtext(axisInfoFn("xtitle",1), 3, line=2)
  mylabel("D")

  panel_evolution_mean(dat, trait =1, xlab = "")
  mylabel("B")

  panel_evolution_mean(dat, trait =2, xlab = "", ylab = "")
  mylabel("E")

  panel_evolution_CV(dat, trait =1, xlab = "")
  mylabel("C")
  mtext("Steps", 1, line=3, xpd=NA)

  panel_evolution_CV(dat, trait =2, xlab = "", ylab = "")
  mylabel("F")
  mtext("Steps", 1, line=3, xpd=NA)
}


figure_evolution_lcc <- function(path){
  figure_evolution_1_trait(path, 1)
}

figure_evolution_hsp <- function(path){
  figure_evolution_1_trait(path, 2)
}

figure_evolution_1_trait <- function(path, trait){

  dat <- lapply(loadStochFile(path), parseLine)

  par(mfrow=c(3,1), oma = c(4,3,2,3), mar= c(3,5,2,1))

  landscape <- getLastLandscape(path)

  mylabel <- function(x) label(-0.15, 1.05, x, adj = c(0, 1), xpd = NA, cex = 1,
    font = 2)

  panel_evolution_raw(dat, step =10, trait=trait, xlab = "", cex=0.5)
  mtext(axisInfoFn("xtitle",trait-1), 3, line=2)
  mylabel("A")

  panel_evolution_mean(dat, trait=trait, xlab = "")
  mylabel("B")

  panel_evolution_CV(dat,trait=trait, xlab = "")
  mylabel("C")
  mtext("Steps", 1, line=3, xpd=NA)
}


figure_evolution_1_trait_PIP <- function(path, val, trait, title=""){

  path <- gsub("XX", val, path)

  mylabel <- function(label, x, y = 1.1, x0 = -0.125, x02 = 0.05) {
    label(x0, y, label, cex = 1.5, xpd = NA, font = 2)
    label(x0+x02, y, x, cex = 1.5, xpd = NA, font = 1)
  }

  dat <- lapply(loadStochFile(path), parseLine)

  m <- matrix(c(1, 2, 3, 4), ncol = 1, byrow = FALSE)
  layout(m, widths = 1, heights = c(0.4, 0.2, 0.2, 0.2))

  par(oma = c(3,3,6,1), mar= c(8,11.75,1,11.75))

  if(trait == 1)
    figure_PIP_lcc(file.path(path, "PIP"))
  else
    figure_PIP_hsp(file.path(path, "PIP"))
  mylabel("A", "Pairwise invasibility plot", 1.2, x0=-0.75, x02 = 0.105)

  par(mar= c(3,3,1,1))

  panel_evolution_raw(dat, step =10, trait=trait, xlab = "", cex=0.5)
  mylabel("B", "Stochastic assembly", y=1.2)
  panel_evolution_mean(dat, trait=trait, xlab = "")
  panel_evolution_CV(dat,trait=trait, xlab = "")
  mtext("Steps", 1, line=3, xpd=NA)

}

panel_evolution_mean <- function(dat,trait =1, xlab = "Time", ylab = "Mean"){

  out <- get_moment(dat, trait, n=1)

  if(trait ==1){
    ylim= c(-1.5, 0.25)
    lab <- c(0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2)
   }
  else{
    ylim= c(0.1, log10(60))
    lab <- c(1.25, 2.5, 5, 10, 20, 40, 80)
  }

  plot(out$x, out$y, type='l',  xlab = "", ylab = "", las=1, ylim=ylim, col="red",
    yaxt = "n")
  mtext(xlab, 1, line=3)
  mtext(ylab, 2, line=3)

  axis(2, at = log10(lab), labels = lab, las = 1)
}


get_moment <- function(dat, trait, n=1){

  trait.mean <- function(z, trait){
    sum(z$Traits[,trait]*z$PopSize)/sum(z$PopSize)
  }


 moment.nth <- function(z, trait, n=2){
    moment.1 <- trait.mean(z, trait)
    moment.2 <- sum((z$Traits[,trait]-moment.1)^n*z$PopSize)/sum(z$PopSize)
    moment.2
  }

  x <- sapply(dat, function(x) x$time)
  if(n==1)
    y <- sapply(dat, trait.mean, trait=trait)
  else if (n > 1)
    y <- sapply(dat, moment.nth, trait=trait, n=n)
  else
    stop(paste("bad value given to moment function", n))

  list(x=x, y=y)
}

panel_evolution_CV <- function(dat,trait =1, xlab = "Time", ylab = "CV (%)"){

  out.var <- get_moment (dat, trait, n=2)
  out.mn <- get_moment (dat, trait, n=1)
  y <- abs(sqrt(out.var$y)/out.mn$y*100)

  plot(out.var$x, y, type='l', xlab = "", ylab = "", las=1, ylim = c(0, max(100, tail(y, 500), na.rm=TRUE)), col="red")
  mtext(xlab, 1, line=3)
  mtext(ylab, 2, line=3)
}

getLastLandscape <- function(path){
  files <- list.files(path, "D.txt")
  x <- as.numeric(sapply(strsplit(files,"-"), function(x) x[2]))
  i <- which(x==max(x))
  files[i]
}

panel_evolution_raw <- function(dat, cex=0.5, pch=15, trait =1, step=10, xlab = "Time", ylab="Raw values"){

  if(trait ==1){
    ylim= c(-1.5, 0.25)
    lab <- c(0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2)
   }
  else{
    ylim= c(0.1, log10(60))
    lab <- c(1.25, 2.5, 5, 10, 20, 40, 80)
  }

  plot(0,0, type='n', xlab = "", ylab = "", las=1, ylim=ylim, xlim = c(0, length(dat)),
    yaxt = "n")
  mtext(xlab, 1, line=3)
  mtext(ylab, 2, line=3)

  axis(2, at = log10(lab), labels = lab, las = 1)

  for(i in seq(1, length(dat), by=step)){
    x <- dat[[i]]$time
    y <- dat[[i]]$Traits[,trait]
    z <- dat[[i]]$PopSize
    points(x+0*y, y, pch=pch, cex=cex, col=make.transparent("red", getColorSclaedDensity(z)))
  }
}

getColorSclaedDensity <- function(x, Min=0.02, Max=10){
  x <- x -Min
  x[x<0] <-  0
  x[x > Max] <- Max
  x <- x/Max*0.9
  x+0.1
}

figure_PIP_lcc <- function(..., xlim=c(-1.8, 0.8)) {

  figure_PIP(..., xlim=xlim)

  lab <- c(0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2)
  for(i in 1:2)
    axis(i, at = log10(lab), labels = lab, las = 1)
}

figure_PIP_hsp <- function(..., xlim=c(0, log10(50))) {

  figure_PIP(..., xlim=xlim)

  dim = 1

 lab <- c(1.25, 2.5, 5, 10, 20, 40, 80)
 for(i in 1:2)
  axis(i, at = log10(lab), labels = lab, las = 1)
}

figure_PIP <- function(path, xlim = NULL, zlim = c(-6, 6)) {

  df <-  read.table(path)
  names(df) <- c("x", "y", "z")
  df <- collapse.grid(df)

  df$z[df$z<zlim[1]] <- zlim[1]
  df$z[df$z>zlim[2]] <- zlim[2]

  myColorMap <- loadColorMap("R/colormap.csv", zlim = zlim)

  if(is.null(xlim))
    xlim = range(x)
  plot(NA, xlim = xlim, ylim = xlim, axes=FALSE, ann = FALSE, xaxs = "i", yaxs = "i", las = 1, asp=1)
  box()

  image(df$x, df$y, df$z, zlim=range(myColorMap$levels), col = myColorMap$col, add=TRUE)
  lines(df$x,df$x)

  mtext("Resident trait", 1, line=3)
  mtext("Variant trait", 2, line=3)

}


