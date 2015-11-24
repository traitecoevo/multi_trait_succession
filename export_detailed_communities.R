#!/usr/bin/Rscript

source("R/figures.R")

printCommunityAsInputFile<-function(community, path, filename){

  traits <- community$Traits
  seed.rain <- community$PopSize
  n <- nrow(traits)
  tmp <- c(`%T`=0, NumRes=n, c(t(traits)), seed.rain)
  names(tmp)[seq(3, length=n*4)] <-
    sprintf("y%d", seq(0, length=n*4))
  names(tmp)[seq(3+n*4, length.out=n)] <-
    sprintf("p%d", seq(0, length=n))
  write.table(data.frame(t(tmp), check.names=FALSE),
              file.path(path, filename), sep="\t", row.names=FALSE,
              quote=FALSE)
  file.path(path, filename)
}

clusterStrategies<-function(community, plotIt =TRUE, h=0.01){

  out <-community

  if(community$nSpecies > 1){

    # run clusterin algorithm
    hc <- hclust(dist(community$Traits))
    # cut tree at given height
    membership <- cutree(hc,h=h)  # can cut tree at different heights, or specificy

    if(plotIt){
      par(mfrow=c(2,1))
      plot(hc)
      plot(community$Traits[,1],community$Traits[,2], col = membership)
      par(mfrow=c(1,1))
    }

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

printDetailedCommunity<-function(community, input.Dir, output.dir, name = "Stoch.txt", trait = 2,  landscape = NA, reprocess = FALSE, resolution =0, viability=1){


  if(!file.exists(output.dir) | reprocess){
    if(!file.exists(output.dir)) #create directory if it doesn't exist
      dir.create(output.dir, recursive=TRUE)

    # copy params file
    file.copy(from = file.path(input.Dir,"params.m"), to = file.path(output.dir,"params.m"))

    #export average community
    filename <- printCommunityAsInputFile(community, output.dir, name)

    #Absolute path to program
    program <- normalizePath("src/evolve")

    # Change directory Run evolve, change back on exit
    owd <- setwd(output.dir)
    on.exit(setwd(owd))

    if(is.na(landscape))
      cmd <- sprintf("%s -L 0 -V %d -P 1 -r %d -f . -s 1", program,  viability, resolution)
    else
      cmd <- sprintf("%s -L 1 -t %d -V %d -P 1  -r %d -f . -s 1", program,  viability, resolution, trait)

    print(cmd)
    system(cmd)
    }

  #return directory
  file.path(output.dir, "res")
}



# for each final community, makes a copy to newDir and prints two versions, 1 based on original file and other using clusterring algorithm to reduce number of stratgies
printFinalCommunity <- function(inputDir, outputDir, step=NA, h=0.01, threshold = 1E-2, ...){
  community <- getEvolvedCommunity(inputDir, step=step)
  community <- removeExtinct(clusterStrategies(community,  plotIt = FALSE, h=h), threshold = threshold)
  printDetailedCommunity(community, inputDir,
          outputDir, ...)
}

base.dir <- "output/data"

for(p in file.path(base.dir,
      c("[60,1]","[10.6,1]","[21.2,0.8]","[60,0.6]","[120,0.8]","[120,1.4]")
      ,"2trait", "base")) {
  printFinalCommunity(inputDir = p,
    outputDir = gsub("output/data", "output/finalCommunities",p),
    landscape=NA, resolution=1, reprocess=FALSE, viability=0, step=5001)
}

for(p in c(file.path(base.dir,"[60,1]","lcc", "[1.3725]/base"),
		 file.path(base.dir,"[60,1]", "hsp", "[0.625]/base")) ) {
	printFinalCommunity(inputDir = p,
		outputDir = gsub("output/data", "output/finalCommunities",p),
		landscape=NA, resolution=1, reprocess=FALSE, viability=0, step=1001,  h=0.02, threshold=0.5)
}
