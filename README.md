# Multitrait successional forest dynamics enable diverse competitive coexistence

Contact: [Daniel Falster](http://danielfalster.com/)

This repository contains all the code used to produce figures in the paper:

Falster, D.S., Brännström, Å., Westoby, M. & Dieckmann, U. (2017) Multitrait successional forest dynamics enable diverse competitive coexistence. *Proceedings of the National Academy of Sciences USA*: doi: [10.1073/pnas.1610206114](http://doi.org/10.1073/pnas.1610206114).

also available as a preprint at

Falster, D.S., Brännström, Å., Westoby, M. & Dieckmann, U. (2015) Multi-trait eco-evolutionary dynamics explain niche diversity and evolved neutrality in forests. *bioRxiv*, 014605. [10.1101/014605](http://doi.org/10.1101/014605)

Similar code was used to generate results of

Falster, D.S., Brännström, Å., Dieckmann, U. & Westoby, M. (2011) Influence of four major plant traits on average height, leaf-area cover, net primary productivity, and biomass density in single-species forests: a theoretical investigation. *Journal of Ecology*, 99, 148–164. doi: [10.1111/j.1365-2745.2010.01735.x](http://doi.org/10.1111/j.1365-2745.2010.01735.x)

The main simulation tool is written in C++ while figures are generated using `R`.

All code is released under the [GNU GENERAL PUBLIC LICENSE](LICENSE) (a requirement of using the GNU Scientific Library).

Please note, this code is supplied so that readers can reproduce the results of the published analysis, if that is desired. But it is not intended for use by others. **If you are interested in running this type of analysis, do not use this code**. Instead checkout the [`plant` package](https://github.com/traitecoevo/plant) for an improved implementation of trait evolution in size- and patch-structured metapopulations.

The following instructions show you how to reproduce figures from the 2017 paper.

First you should download this code repository.

## Compiling program

The main evolutionary simulation are run in C++. Build the C++ program `evolve` within the folder `src` by running `make` within that directory. To compile the program you'll need the [GNU Scientific Library](http://www.gnu.org/software/gsl/) installed.

## Run simulations 

The file `launch_all.bash` contains code for launching the simulations. Warning, this step will take a **long** time. Time estimates for completion are indicated in the launch file but in short: the entire collection of results takes > 1000 days of computer time to complete.

After the initial simulations have completed, you'll need to generate some additional outputs using those results:

1. Some additional snapshots of the fitness landscape for the assemblies, using `landscape_print_sequence.bash`.

2. Export hi-resolution versions of the final communities, using `export_detailed_communities.R`.

## Download my simulation output

For those wishing to use already generated data, skip steps 1-3 by downloading simulation output from the [releases page](https://github.com/traitecoevo/multi_trait_succession/releases). This is the data used to produce the figures in the paper. Once downloaded and unzipped the folder `output` should be placed in the root of this project.

## Figures

You can then run material in `analysis.R` to generate figures in the folder `output/figures`.

To generate figures you'll need to install the following `R` packages:

```
install.packages(c("rgl", "xlsx", "downloader", "devtools", "akima", "smatr"))
```
