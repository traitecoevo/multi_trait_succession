
library(rgl)
library(standViz)


addXYCoordinates <- function(treeData, xlim = c(1, 10), ylim = c(0, 10)) {
  noTrees <- nrow(treeData)

  if (!is.null(noTrees)) {
    treeData$x <- runif(noTrees, min = xlim[1], max = xlim[2])
    treeData$y <- runif(noTrees, min = ylim[1], max = ylim[2])
  }

  treeData
}

trial_patch <- function(patchSize = c(1, 1), age = 15, x0 = 0, y0 = 0) {

  treeData <- getStandData(age, Res, siteArea = prod(patchSize))

  treeData <- addXYCoordinates(treeData, xlim = c(x0, x0 + patchSize[1]), ylim = c(y0,
    y0 + patchSize[2]))

  makeStand(treeData)
}


plotPatchOutline <- function(coords) {
  M <- cbind(coords$x, coords$y, 0)
  M <- rbind(M, M[1, ])
  lines3d(M, col = "darkgrey")
}

getAgeSequence <- function(Res, ages = seq(0, 150, by = 5), patchSize = c(0.1,
  5), patchSpace = 0.1, patchCoordinatesFn = figure1a_getPatchCoordinates, rotateBy = 0,
  scaling = c(1, 1), col = defaultColours()) {

  treeData <- lapply(ages, getStandData, Res = Res, siteArea = prod(patchSize),
    col = col)

  coords <- lapply(1:length(treeData), patchCoordinatesFn, patchSize = patchSize,
    patchSpace = patchSpace)

  for (i in 1:length(treeData)) if (!is.null(treeData[[i]]))
    treeData[[i]] <- addXYCoordinates(treeData[[i]], xlim = range(coords[[i]]$x),
      ylim = range(coords[[i]]$y))

  if (abs(rotateBy) > 0) {
    for (i in 1:length(treeData)) {
      coords[[i]] <- rotateData(coords[[i]], theta = rotateBy)
      treeData[[i]] <- rotateData(treeData[[i]], theta = rotateBy)
    }
  }

  allData <- do.call(rbind, treeData)
  list(stand = makeStand(allData, zvals = c(0, seq(0.5, 0.7, by = 0.1), seq(0.75,
    0.9, by = 0.05), seq(0.91, 1, by = 0.02), 1), scaling = scaling), coords = coords)
}

# rotate plot by 45deg (http://en.wikipedia.org/wiki/Rotation_matrix)
rotateData <- function(data, theta) {
  data.new <- data
  data.new$x <- data$x * cos(theta) - data$y * sin(theta)
  data.new$y <- data$x * sin(theta) + data$y * cos(theta)

  data.new
}

sampleStandDensity <- function(Res, age, ...) {
  # find closest row matching desired age
  i <- which.min(abs(Res$pAge[, 1] - age))
  OUT <- NULL
  for (n in 1:Res$nSpp) {
    ml <- sampleSpeciesFromStand(Res, i, n, ...)
    if (!is.null(ml))
      OUT <- rbind(OUT, data.frame(Species = n, ml = ml))
  }
  OUT
}

sampleSpeciesFromStand <- function(Res, row, n, siteArea) {
  i <- row
  j <- !is.na(Res$Lam[[n]][i, ])
  Size <- rev(t(Res$Size[[n]][i, j]))
  Lam <- rev(t(Res$Lambda[[n]][i, j]))

  # note, Reverse order of x and y so x is in increasing size

  # calculate total number individuals
  nTot <- sum(Lam)
  # calculate number of indivs to sample
  Nsample <- floor(nTot * siteArea)

  sizes <- NULL
  if (Nsample > 0)
    sizes <- sample(Size, Nsample, replace = TRUE, prob = Lam/nTot)
  sizes
}

# return stand density at age Returns list of trees sampled from continuous
# density distribution
getStandData <- function(age, Res, col = defaultColours(), ...) {
  out <- sampleStandDensity(Res, age, ...)
  treeData <- NULL
  if (!is.null(out)) {
    # calculate heights, dbh, LA
    lma <- Res$strat$Tr0[out$Species]
    LfAr <- out$ml/lma
    topHeight <- Height(out$ml, LMA = lma, a1 = Res$params["p.a1", 1], B1 = Res$params["p.B1",
      1])
    HDB.A <- 0.03062551
    HDB.B <- 0.5068199  # parameters relating DBH to leaf area in Coweeta dataset
    treeData <- data.frame(crownColor = col[out$Species], topHeight = topHeight,
      heightCrownBase = 0, dbh = 2 * HDB.A * LfAr^HDB.B, crownWidth = 0.3 *
        topHeight, crownShape = "yokozawa", eta = 13, stringsAsFactors = FALSE)
  }
  treeData
}

drawRandomAges <- function(ageDistribution, nPatches) {

  # draw patch ages from distribution
  x <- ageDistribution[, 1]
  y <- ageDistribution[, 2]
  sampleFromCDF(x, CDF = c(0, cum_trapz(x, y))/trapz(x, y), nPatches)
}

# sample using Inversion
# sampling(http://en.wikipedia.org/wiki/Inversion_sampling) uses computed CDF
# of a function, invert it, then sample a random number and interpolate to get
# sample

sampleFromCDF <- function(x, CDF, n) {
  sample.fn <- approxfun(CDF, x, rule = 1)
  sample.fn(runif(n))
}


plot_stand <- function(data, size = c(1, 1, 1500, 500)) {
  plotStand(makeStand(NULL), size = size)
  plotStand(data$stand, add = TRUE)
}
