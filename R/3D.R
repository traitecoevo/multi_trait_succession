
library(rgl)

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


#' Takes a dataframe of tree values, checks for essential variables and fills missing
#' values with defaults
#' @title Takes a dataframe of tree values, checks for essential variables and fills missing
#' values with defaults
#' @param treeData A dataframe of values for the stand.
#' @param crownShape Shape of crown. See \code{make3dShape} for allowable values.
#' @param eta Canopy shape factor, required when shape is set as \code{yokozawa}.
#' See \code{yokozawaDensity} for details.
#' @param crownColor Colour of the crown
#' @param stemColor Colour of the stem
#' @return A dataframe of values for the stand.
makeStandDataframe <- function(treeData, crownShape = "yokozawa", eta = 5, crownColor = "forestgreen",
  stemColor = "brown") {

  noTrees <- nrow(treeData)

  if (!is.null(noTrees)) {
    # Fill missing values as needed
    for (v in c("crownShape", "eta", "crownColor", "stemColor")) if (is.null(treeData[[v]]))
      treeData[[v]] <- rep(get(v), noTrees)

    treeData[["species"]] <- as.character(as.factor(paste(treeData$crownShape,
      treeData$eta, treeData$crownColor)))

    essentialVariables <- c("x", "y", "topHeight", "heightCrownBase", "crownWidth",
      "dbh")
    for (v in essentialVariables) if (is.null(treeData[[v]]))
      stop(paste(v, "values needed to make stand"))
  }
  treeData
}

#' Creates a 3D shapeMatrix suitable for passing to \code{rgl::plot.triangles}, and used
#' by \code{plot3dShape}) This object has height 1, width 1, and is located at (0,0,0). Use
#' the function \code{resize3dShape} to scale and move to a different set of coordinates
#' @title Creates a 3D shapeMatrix suitable
#' @param shape Desired shape of object
#' @param nz Number of vertical points used to generate the shape. Only used if zvals=NA.
#' @param nalpha Number of angles at which nodes in the canopy are created
#' @param zvals Heights to calculate points. Overrides specification of \code{nz}.
#' @param eta Canopy shape factor, required when shape is set as \code{yokozawa}.
#' See \code{yokozawaDensity} for details.
#' @return A matrix of coordinates with columns [x,y,z] defining the shape
#' @export
make3dShape <- function(shape = c("cone", "ellipsoid", "halfellipsoid", "paraboloid",
  "cylinder", "yokozawa"), nz = 25, nalpha = 25, zvals = NA, eta = 13) {
  # Canopy shape
  shape <- match.arg(shape)
  # Define canopy shape function. Gives relative distance from centre at
  # different relative heights within the crown
  if (shape == "cone")
    distfun <- function(z) (1 - z)
  if (shape == "ellipsoid")
    distfun <- function(z) sqrt(1 - ((z - 1/2)^2)/((1/2)^2))
  if (shape == "halfellipsoid")
    distfun <- function(z) sqrt(1 - z^2)
  if (shape == "paraboloid")
    distfun <- function(z) sqrt(1 - z)
  if (shape == "cylinder")
    distfun <- function(z) 1
  if (shape == "yokozawa") {
    distfun <- function(z) yokozawaDensity(z, 1, eta)
  }
  # Heights to calculate points, if not already defined
  if (is.na(zvals[1]))
    zvals <- seq(1, 0, length = nz) else nz <- length(zvals)
  # Make matrix of x,y,z values, by choosing points with given radius and angle
  # at a given height
  radius <- rep(0.5 * distfun(zvals)/max(distfun(zvals)), times = nalpha)  #scaled to specified object width
  angles <- rep(seq(0, 2 * pi, length = nalpha), each = nz)
  x <- radius * cos(angles)
  y <- radius * sin(angles)
  z <- rep(zvals, times = nalpha)
  m <- matrix(cbind(x, y, z), ncol = 3)
  # determine ordering of points to make triangles. Each row of tm gives indices
  # for the vertices of a single triangle traverse each vertical spines and
  # define two triangles between each set of four points (adjacent splines and
  # heights)
  tm <- matrix(nrow = 3, ncol = (nalpha - 1) * 2 * (nz - 1))
  r <- 1
  for (ai in 1:(nalpha - 1)) {
    # do for all angles do for each height
    for (zi in 1:(nz - 1)) {
      # first (upper) triangle
      tm[1, r] <- (ai - 1) * nz + zi
      tm[2, r] <- (ai) * nz + zi
      tm[3, r] <- (ai) * nz + zi + 1
      # second (lower) triangle
      tm[1, r + 1] <- (ai - 1) * nz + zi
      tm[2, r + 1] <- (ai - 1) * nz + zi + 1
      tm[3, r + 1] <- (ai) * nz + zi + 1
      # increment index
      r <- r + 2
    }
  }
  # return matrix of triangle coordinates (note, includes duplicate points
  # because each used in multiple triangles)
  m[tm, 1:3]
}
#' Creates a stand of trees to plot in 3D
#' @title Creates a stand of trees to plot in 3D
#' @param treeData A dataframe of values for the stand. Required values are defined in code{makeStandDataframe}
#' @param scaling a vector of length two, giving scaling of heights and widths in the plot
#' @param ... Additional arguments defining crown, to pass through to \code{make3dShape}
#' @export
makeStand <- function(treeData, scaling = c(1, 1), ...) {

  myTreeData <- makeStandDataframe(treeData)

  stand <- list()
  t <- 0

  for (s in unique(myTreeData$species)) {

    treeData.sp <- subset(myTreeData, myTreeData$species == s)

    # Make basic shapes for crown and stem for reuse
    crownShapeMatrix <- make3dShape(shape = as.character(treeData.sp$crownShape[1]),
      treeData.sp$eta[1], ...)

    stemShapeMatrix <- make3dShape(shape = "cone", nz = 2, nalpha = 5)

    noTrees <- nrow(treeData.sp)

    if (!is.null(noTrees)) {
      for (i in 1:noTrees) stand[[t + i]] <- makeTree(x = treeData.sp$x[i],
        y = treeData.sp$y[i], topHeight = treeData.sp$topHeight[i], heightCrownBase = treeData.sp$heightCrownBase[i],
        crownWidth = treeData.sp$crownWidth[i], dbh = treeData.sp$dbh[i],
        crownShapeMatrix = crownShapeMatrix, stemShapeMatrix = stemShapeMatrix,
        crownColor = treeData.sp$crownColor[i], stemColor = treeData.sp$stemColor[i],
        scaling = scaling)
      t <- t + i
    }

  }
  class(stand) <- "standVizStand"
  stand
}
#' Creates a single tree with specified characteristics. Uses existing values for \code{stemShapeMatrix}
#' and \code{crownShapeMatrix} if provided, leading
#' to faster computing times. Otherwise recalculates these.
#' @title Creates a single tree with specified characteristics
#' @param x Location in x dimension
#' @param y Location in x dimension
#' @param topHeight Top height of tree
#' @param heightCrownBase Height of crown base
#' @param crownWidth Width of crown
#' @param dbh Diameter at breast height
#' @param crownShapeMatrix A matrix with columns [X,Y,Z] describing shape of the crown
#' @param stemShapeMatrix A matrix with columns [X,Y,Z] describing shape of the stem
#' @param crownShape One of 'cone', 'elipsoid', 'ellipsoid', 'round', 'halfellipsoid', 'paraboloid', 'cylinder', 'yokozawa'
#' @param eta Canopy scaling factor. See \code{yokozawaDensity} for details.
#' @param crownColor Colour of the crown
#' @param stemColor Colour of the stem
#' @param nz Number of vertical points used to generate the shape. Only used if zvals=NA.
#' @param nalpha Number of angles at which nodes in the canopy are created
#' @param zvals Heights to calculate points. Overrides specification of \code{nz}.
#' @param scaling a vector of length two, giving scaling of heights and widths in the plot
#' @param ... Additional arguments defining crown, to pass through to \code{make3dShape}
#' @export
#' @seealso \code{\link{plotStand}}, \code{\link{plotTree}}
makeTree <- function(x = 0, y = 0, topHeight = 1, heightCrownBase = 0, crownWidth = 1,
  dbh = 0.01, crownShapeMatrix = NA, crownShape = c("cone", "elipsoid", "ellipsoid",
    "round", "halfellipsoid", "paraboloid", "cylinder", "yokozawa"), eta = 13,
  crownColor = "forestgreen", stemShapeMatrix = NA, stemColor = "brown", nz = 25,
  nalpha = 25, zvals = NA, scaling = c(1, 1), ...) {

  # Makes a basic crown shape, with topHeight 1.0, width 1.0, located at (0,0,0)
  if (is.na(crownShapeMatrix[1]))
    crownShapeMatrix <- make3dShape(shape = match.arg(crownShape), eta = eta,
      nz = nz, nalpha = nalpha, zvals = zvals, ...)

  # Makes a basic stemShape, with topHeight 1.0, width 1.0, located at (0,0,0)
  if (is.na(stemShapeMatrix[1]))
    stemShapeMatrix <- make3dShape(shape = "cone", nz = 2, nalpha = nalpha,
      zvals = zvals)

  myTree <- list(crown = resize3dShape(crownShapeMatrix, height = (topHeight -
    heightCrownBase) * scaling[1], width = crownWidth * scaling[2], x = x,
    y = y, z = heightCrownBase * scaling[1]), crownColor = crownColor, stem = resize3dShape(stemShapeMatrix,
    height = topHeight * scaling[1], width = dbh * scaling[2], x = x, y = y,
    z = 0), stemColor = stemColor)

  class(myTree) <- "standVizTree"

  myTree
}
#' Creates a new rgl scene with desired dimensions
#' @title Creates a new rgl scene with desired dimensions
#' @param size A vector of dimensions for the plot, as specified in \code{newRgl}
#' @export
newRgl <- function(size = c(1, 1, 500, 500)) {
  rgl::rgl.open()
  rgl::par3d(windowRect = size)
  rgl::par3d(mouseMode = rep("none", 3))
  rgl::bg3d("white")
  rgl::rgl.viewpoint(0, -65, fov = 0, scale = c(1, 1, 1), zoom = 1)
}
#' Adds a shape node to the current scene
#' @title Adds a shape node to the current scene
#' @param shapeMatrix Matrix with columns [X,Y,Z]
#' @param ... Additional arguments passed through to \code{rgl::rgl.triangles}
#' @export
plot3dShape <- function(shapeMatrix, ...) {
  rgl::rgl.triangles(shapeMatrix[, 1], shapeMatrix[, 2], shapeMatrix[, 3], ...)
}
#' Plot a stand of trees. To optimsie plotting to rgl scene, the function groups together everything
#' with similar colour and plots these in a single call.
#' @title Plot a stand of trees
#' @param stand A list containing trees, as generated \code{makeStand}
#' @param add Add to existing plot?
#' @param size A vector of dimensions for the plot, as specified in \code{newRgl}
#' @param ... Additional arguments passed through to \code{rgl::rgl.triangles}
#' @export
plotStand <- function(stand, add = FALSE, size = c(1, 1, 500, 500), ...) {

  if (!add)
    newRgl(size = size)

  crownColors <- sapply(stand, function(tree) tree[["crownColor"]])
  stemColors <- sapply(stand, function(tree) tree[["stemColor"]])

  getItem <- function(n, stand, v = "crown") stand[[n]][[v]]

  for (i in unique(crownColors)) plot3dShape(do.call(rbind, lapply(which(i ==
    crownColors), getItem, stand = stand)), col = i, ...)

  for (i in unique(stemColors)) plot3dShape(do.call(rbind, lapply(which(i ==
    stemColors), getItem, stand = stand, v = "stem")), col = i, ...)
}
#' Plot a single tree in 3D
#' @title Plot a single tree in 3D
#' @param tree Top height of plant
#' @param ... Additional arguments passed through to \code{rgl::rgl.triangles}
#' @export
plotTree <- function(tree, ...) {
  plot3dShape(tree$crown, col = tree$crownColor, ...)
  plot3dShape(tree$stem, col = tree$stemColor, ...)
}
#' Scale and reposition a shape.
#' @title Scale and reposition a shape
#' @param shapeMatrix Initial shape, defined by a shapeMatrix with columns [X,Y,Z]
#' @param height Desired height
#' @param width Desired width
#' @param x Desired offset in X dimension
#' @param y Desired offset in Y dimension
#' @param z Desired offset in Z dimension
#' @return Resized and repositioned shape, defined by a shapeMatrix with columns [X,Y,Z]
#' @export
resize3dShape <- function(shapeMatrix, height = 1, width = 1, x = 0, y = 0, z = 0) {
  shapeMatrix[, 1] <- x + 0.5 * width * (shapeMatrix[, 1]/max(shapeMatrix[, 1]))
  shapeMatrix[, 2] <- y + 0.5 * width * (shapeMatrix[, 2]/max(shapeMatrix[, 2]))
  shapeMatrix[, 3] <- z + height * (shapeMatrix[, 3]/max(shapeMatrix[, 3]))
  shapeMatrix
}
#' Function describing vertical distribution of leaf area within the crown of a plant.
#' First described by Yokozawa 1995 doi: 10.1006/anbo.1995.1096
#' @title Function describing vertical distribution of leaf area within the crown
#' @param z Height within canopy of the plant
#' @param h Top height of plant
#' @param eta Canopy shape factor
#' @return A vector of values, with same length as z
#' @export
#' @examples
#' # Plot the crown profile function for Yokozawa's function:
#' z <- seq(0, 1, 0.01)
#' q <- yokozawaDensity(z, 1, 12)
#' plot(q, z, type = 'l')
yokozawaDensity <- function(z, h, eta) {
  q <- 2 * eta * (1 - (z/h)^eta) * (z/h)^eta/z
  q[is.na(q)] <- 0
  q
}
