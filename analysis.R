source("R/figures.R")
source("R/utils.R")
source("R/wright_2004.R")
source("R/3D.R")

# Figures -----------------------------
dir.create("output/figures", FALSE)

# Figure 1
f1_data <- figure1_data()

# a) 3D plot of stand development
scaling <- 0.75
data <- getAgeSequence(f1_data$Res, f1_data$ages, patchSize = c(30/50, 6), patchSpace = 0,
  scaling = rep(scaling, 2), col = f1_data$col.hi)
plot_stand(data, size = c(1, 1, 2000, 2000))
add_height_bar(0.5,6,vadj=0.75, H=5*scaling)
rgl.viewpoint(0, -65, fov = 0, scale = c(1, 1, 1), zoom = 0.7)
snapshot3d("output/figures/fig1a.png")
rgl.close()

# b) 3D plot of metapopulation
data <- figure_1b_data(f1_data)
plot_stand(data, size = c(1, 1, 2000, 2000))
tmp <- lapply(data$coords, plotPatchOutline)
rgl.viewpoint(0, -65, fov = 0, scale = c(1, 1, 1), zoom = 0.7)
snapshot3d("output/figures/fig1b.png")
rgl.close()

# c) plot traits & LAI with respect to time
to.pdf(figure1_LAI(f1_data$Res, cols = f1_data$col), "output/figures/fig1c.pdf",
  height = 7, width = 15)

# Figure 2
to.pdf(figure2(), "output/figures/fig2.pdf", height = 8, width = 15)

# Figure 3
to.pdf(figure3(), "output/figures/fig3.pdf", height = 7, width = 4)

# Figure 4
to.pdf(figure4(), "output/figures/fig4.pdf", height = 8, width = 8)

# Figure 5
# a) panels of community average LMA and height
out <- figure5_data(f = parallel::mclapply)
to.pdf(figure5(out), "output/figures/fig5a.pdf", height = 7, width = 12)

# b) 3D pictures for select plots
D <- c(10.6, 60, 120, 21.2, 60, 120)
Prod <- c(1, 0.6, 0.8, 0.8, 1, 1.4)
Scale <- c(60, 100, 120, 80, 100, 120)/120

for (i in seq_along(D)) {
  data <- figure_get_example_patch(D[i], Prod[i], Scale[i])
  figure_plot_patch(data)
  rgl.viewpoint(0, -65, fov = 0, scale = c(1, 1, 0.8), zoom = 0.8)
  filename <- file.path("output/figures/fig5b-", paste0("[",D[i],",",Prod[i],"].png"))
  snapshot3d(filename)
  rgl.close()
}

# Suppmat figures --------------------------------------------------

# Leaf
to.pdf(plotleaf(), paste0("output/figures/figS-leaf.pdf"), height = 6, width = 10)

# Height
to.pdf(figure_ReproductiveAllocation(), paste0("output/figures/figS-height-a.pdf"),  height = 5, width = 12)

# All fitness landscapes
to.pdf(fig_2D_landscapes_all(), "output/figures/figS_landscapes_all.pdf", height = 8, width = 6)

# check sensitivity to parameter changes
to.pdf(figure_sensitivity_1(), "output/figures/figS-elasticity-pars_2D.pdf", height = 12, width = 12)

to.pdf(figure_sensitivity_2(), "output/figures/figS-elasticity-repro_2D.pdf", height = 4, width = 10)

to.pdf(figure_sensitivity_3(), "output/figures/figS-elasticity-repro_1D.pdf", height = 8, width = 8)

# Seed rain
to.pdf(figure_seed_rain(), paste0("output/figures/figS-seed_rain_height.pdf"), height = 6, width = 6)

# Competitive feedbacks
to.pdf(figure_competition(), paste0("output/figures/figS-competition.pdf"), height = 4, width =10)

# Plot trait assembly and evolution
sanitise <- function(x) { gsub(".", "_", xx, fixed=TRUE)}

for(xx in c(0.8175, 1.3725)) {
  f <- sprintf("output/figures/figS-evolution_lcc_%s_1.pdf", sanitise(xx))
  to.pdf(figure_evolution_1_trait_PIP("output/data/[60,1]/lcc/[XX]/base", xx, 1), f, height = 10, width = 6)

  f <- sprintf("output/figures/figS-evolution_lcc_%s_2.pdf", sanitise(xx))
  to.pdf(
    figS_Assembly_lcc_C(sprintf("output/data/[60,1]/lcc/[%s]/base", xx))
    , f, height = 10, width = 10)
}

xx <- -0.2
f <- sprintf("output/figures/figS-evolution_hsp_%s_1.pdf", sanitise(xx))
to.pdf(
  figure_evolution_1_trait_PIP("output/data/[60,1]/hsp/[XX]/base", xx, 2)
  , f,height = 10, width = 6)

f <- sprintf("output/figures/figS-evolution_hsp_%s_2.pdf", sanitise(xx))
to.pdf(
  figS_Assembly_hsp_C(sprintf("output/data/[60,1]/hsp/[%s]/base", xx))
  , f, height = 10, width = 10)

xx  <- 0.625
f <- sprintf("output/figures/figS-evolution_hsp_%s_1.pdf", sanitise(xx))
to.pdf(
  figure_evolution_1_trait_PIP("output/data/[60,1]/hsp/[XX]/base", xx, 2)
  , f,height = 10, width = 6)

f <- sprintf("output/figures/figS-evolution_hsp_%s_2.pdf", sanitise(xx))
to.pdf(
  figS_Assembly_hsp_C(sprintf("output/data/[60,1]/hsp/[%s]/base", xx),
    times = c(seq(0, 20, by = 5), seq(40, 180, by = 20),
      seq(200, 3600, by = 200)))
  , f, height = 10, width = 10)


# Two trait evolution
to.pdf(figure_evolution_2_traits("output/data/[60,1]/2trait/base/"), paste0("output/figures/figS-assembly_60.pdf"), height = 10, width = 10)

# Run with lower mutation rate
xx <- 0.625
f <- sprintf("output/figures/figS-evolution_hsp_%s_low_1.pdf", sanitise(xx))
to.pdf(
  figure_evolution_1_trait_PIP("output/data/[60,1]/hsp/[XX]/I[0.1]/base", xx, 2)
  , f, height = 10, width = 6)

f <- sprintf("output/figures/figS-evolution_hsp_%s_low_2.pdf", sanitise(xx))
to.pdf(
  figS_Assembly_hsp_C(sprintf("output/data/[60,1]/hsp/[%s]/I[0.1]/base", xx),
    times = seq(0, 6000-1, by = 200))
  , f, height = 10, width = 10)

# Resolve mixture with no extra mutation
to.pdf(
  figure_evolution_1_trait_PIP("output/data/[60,1]/hsp/[XX]/resolve", xx, 2)
  , sprintf("output/figures/figS-evolution_hsp_%s_resolve_1.pdf", sanitise(xx)),
  height = 10, width = 6)

to.pdf(
  figS_Assembly_hsp_C(sprintf("output/data/[60,1]/hsp/[%s]/resolve", xx),
    times = seq(0, 6000-1, by = 200))
  , sprintf("output/figures/figS-evolution_hsp_%s_resolve_2.pdf", sanitise(xx)),
  height = 10, width = 10)
