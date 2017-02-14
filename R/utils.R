
to.dev <- function(expr, dev, filename, ..., verbose = TRUE) {
  if (!file.exists(dirname(filename)))
    dir.create(dirname(filename), recursive = TRUE)
  if (verbose)
    cat(sprintf("Creating %s\n", filename))
  dev(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

to.pdf <- function(expr, filename, ..., verbose = TRUE) {
  to.dev(expr, pdf, filename, ..., verbose = verbose)
}

to.png <- function(expr, filename, ..., verbose = TRUE) {
  to.dev(expr, png, filename, ..., verbose = verbose)
}

# convert to png
convert_pdf_to_png <- function(filename) {
  system(paste("convert -density 300 -background white", filename, gsub(".pdf", ".png", filename, fixed = TRUE)))
}

# make png version without white background
convert_png_to_transparent_png <- function(filename) {
  system(paste("convert", filename, "-transparent white ", filename))
}

# return last element in dataframe or vector
last <- function(x) {
  tail(x, n = 1)
}

# make colours semitransparent:
make.transparent <- function(col, opacity = 0.5) {
  if (length(opacity) > 1 && any(is.na(opacity))) {
    n <- max(length(col), length(opacity))
    opacity <- rep(opacity, length.out = n)
    col <- rep(col, length.out = n)
    ok <- !is.na(opacity)
    ret <- rep(NA, length(col))
    ret[ok] <- Recall(col[ok], opacity[ok])
    ret
  } else {
    tmp <- col2rgb(col)/255
    rgb(tmp[1, ], tmp[2, ], tmp[3, ], alpha = opacity)
  }
}

# returns index of 2D matrix m[i,j] when read as 1D vector
matrixIndexAsVectorIndex <- function(i, j, ncol) {
  (i - 1) * ncol + j
}

collapse.grid <- function(data) {
  # does reverse of expand.grid function takes a data frame with up to 2 or 3
  # cols and converts into x,y,z vectors such that data[,c('x','y)] =
  # expand.grid(x,y)

  # find dimensions of X and Y vectors
  ncol <- match(FALSE, data$y[2:length(data$y)] > data$y[1:(length(data$y) - 1)])
  nrow <- as.integer(length(data$x)/ncol)
  # extract x and y vectors
  out <- NULL
  out$y <- data$y[1:ncol]
  out$x <- data$x[ncol * seq(0, (nrow - 1)) + 1]

  # reshape Z data into matrix
  out$z <- matrix(data$z[seq_len(ncol*nrow)], ncol = ncol, nrow = nrow, byrow = TRUE)
  out
}

# position label at a fractional x/y position on a plot
label <- function(px, py, lab, ..., adj = c(0, 1)) {
  usr <- par("usr")
  x <- usr[1] + px * (usr[2] - usr[1])
  y <- usr[3] + py * (usr[4] - usr[3])

  if (par("ylog"))
    y <- 10^y
  if (par("xlog"))
    x <- 10^x

  text(x, y, lab, adj = adj, ...)
}

# mixes two colours with fraction p
mix2Colours <- function(cols, col2, p) {
  m <- col2rgb(cols)
  m2 <- col2rgb(rep(col2, length = length(cols)))
  m3 <- (m * p + m2 * (1 - p))/255
  rgb(m3[1, ], m3[2, ], m3[3, ])
}

# convert colour into HSV space
col2hsv <- function(cols) rgb2hsv(col2rgb(cols))

# adjust HSV values of cols, p is matrix of Hue,Saturation,Brightness (or
# Value) each in range 0-1 returns hexadecimal colours
adjustHSV <- function(cols, p = c(0, 0, 0)) {
  m <- col2hsv(cols) + p

  # truncate new values at allowable bounds
  for (i in 1:3) {
    m[i, m[i, ] < 0] <- 0
    m[i, m[i, ] > 1] <- 1
  }

  hsv(m[1, ], m[2, ], m[3, ])
}


linear.rescale <- function(x, range, scale = range(x)) {
  p <- (x - scale[[1]])/(scale[[2]] - scale[[1]])
  range[[1]] + p * (range[[2]] - range[[1]])
}


# simple trapezoidal numerical integration using (x, y) values in vectors x and
# y
trapz <- function(x, y) sum(diff(x) * (y[-length(y)] + y[-1])/2)

# cumulative cumulative integral using Simple trapezoidal numerical integration
cum_trapz <- function(x, y) cumsum(diff(x) * (y[-length(y)] + y[-1])/2)
