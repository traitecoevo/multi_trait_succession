# Trait-based tradeoffs generate competitive coexistence and evolved neutrality in forests

Contact: [Daniel Falster](http://danielfalster.com/)

This repository contains all the code used to produce figures in the manuscript:

Falster, D.S., Brännström, Å., Westoby, M. & Dieckmann, U. (2015) Multi-trait eco-evolutionary dynamics explain niche diversity and evolved neutrality in forests. *bioRxiv*, 014605. [10.1101/014605](http://doi.org/10.1101/014605)

currently in review at *Science Advances*. Similar code was used to generate results of 

Falster, D.S., Brännström, Å., Dieckmann, U. & Westoby, M. (2011) Influence of four major plant traits on average height, leaf-area cover, net primary productivity, and biomass density in single-species forests: a theoretical investigation. *Journal of Ecology*, 99, 148–164. doi: [10.1111/j.1365-2745.2010.01735.x](http://doi.org/10.1111/j.1365-2745.2010.01735.x)

The main simulation tool is written in C++ while figures are generated using `R`. Code is released under the [GNU GPL](LICENSE).

Please note, this code is supplied so that readers can reproduce the results of the published analysis, if that is desired. But it is not intended for use by others. **If you are interested in running this type of analysis, do not use this code**. Instead checkout the [`plant` package](https://github.com/traitecoevo/plant) for an improved and far-superior implementation. 

The following instructions show you how to reproduce figures from the 2015 paper. First you should download this code repository. 

## Compiling program

First you need to build the C++ program in the folder `src` by running `make` within that directory. You'll need the [GNU Scientific Library](http://www.gnu.org/software/gsl/) installed.

## Run simulations 

The file `launch_all.bash` contains code for launching the simulations. Warning, this will take a **long** time. Time estimates for completion are indicated in the launch file. 

After the initial simulations have completed, you'll need to generate some additional outputs using those results:

1. Some additional snapshots of the fitness landscape for the assemblies, using `landscape_print_sequence.bash`.

2. Export hi-resolution versions of the final communities, using `export_detailed_communities.R`.

For those wishing to use already generated data, skip steps 1-3 by downloading simulation output from the [releases page](https://github.com/traitecoevo/evolved_neutrality/releases).

## Figures

To generate figures you'll need to install the following packages:

```
install.packages(c("rgl", "xlsx", "downloader", "devtools", "akima", "smatr"))
devtools::install_github("dfalster/standViz")
```

You can then run material in `analysis.R` to generate figures.
