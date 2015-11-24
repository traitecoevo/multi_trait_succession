library(xlsx)
library(downloader)

process_wright_2004 <- function() {

  filename <- "output/Wright2004.xls"

  if(!file.exists(filename)) {
    url <- "http://www.nature.com/nature/journal/v428/n6985/extref/nature02403-s2.xls"
     download(url, filename, mode="wb")
  }

  d <- read.xlsx2(filename, sheetIndex = 1, startRow = 11, stringsAsFactors = FALSE, check.names = FALSE)

  ## Do some name translations:
  tr <- c("Code"="Code",
          "Dataset"="Dataset",
          "BIOME"="Biome",
          "Species"="Species",
          "GF"="GrowthForm",
          "Decid/E'green"="Deciduous",
          "Needle/Broad lf"="Needle",
          "C3C4"="C3",
          "N2-fixer"="N2fixer",
          "log LL"="LogLeafLifespan",
          "log LMA"="LogLMA",
          "log Nmass"="Log.N.mass",
          "log Narea"="Log.N.area",
          "log Pmass"="Log.P.mass",
          "log Parea"="Log.P.area",
          "log Amass"="Log.A.mass",
          "log Aarea"="Log.A.area",
          "log Gs"="Log.Gs",
          "log Rdmass"="Log.Rd.mass",
          "log Rdarea"="Log.Rd.area",
          "Ca - Ci"="CaCi")
  names(d)[match(names(tr), names(d))] <- tr

  ## Drop blank columns
  d <- d[names(d) != " "]

  ## Data tweaking:
  d[["Code"]] <- as.integer(d[["Code"]])
  d[["CaCi"]] <- as.numeric(d[["CaCi"]])

  names(d) <- gsub("Log\\.", "Log", names(d))
  re <- "Log"
  i_log <- grep(re, names(d))
  d[i_log] <- lapply(d[i_log], as.numeric)
  d_unlogged <- as.data.frame(10^d[i_log])
  names(d_unlogged) <- sub(re, "", names(d_unlogged))

  data<- cbind(d[-c(i_log)], d_unlogged)

  # lowercase names
  names(data) <- tolower(names(data))

  ## Convert to Kg from g
  data[["lma"]] <- data[["lma"]]/1000

  ## Convert to years from month
  data[["leaflifespan"]] <- data[["leaflifespan"]]/12
  ## Convert to 1/year from year
  data[["k_l"]] <- 1/data[["leaflifespan"]]

  data
}
