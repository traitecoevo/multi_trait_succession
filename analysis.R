
source("R/figures.R")
source("R/utils.R")
source("R/3D.R")
source("R/wright_2004.R")

# Figure 1
f1_data <- figure1_data()

# a) 3D plot of stand development
data <- getAgeSequence(f1_data$Res, f1_data$ages, patchSize = c(30/50, 6), patchSpace = 0,
  scaling = c(0.75, 0.75), col = f1_data$col.hi)
plot_stand(data, size = c(1, 1, 2000, 500))
rgl.viewpoint(0, -65, fov = 0, scale = c(1, 1, 1), zoom = 1)
snapshot3d("output/figures/fig1a.png")
rgl.close()

# b) 3D plot of metapopulation
data <- figure_1b_data(f1_data)
plot_stand(data)
tmp <- lapply(data$coords, plotPatchOutline)
rgl.viewpoint(0, -65, fov = 0, scale = c(1, 1, 1), zoom = 0.7)
snapshot3d("output/figures/fig1b.png")
rgl.close()

# c) plot traits & LAI with respect to time
to.pdf(figure1_LAI(f1_data$Res, cols = f1_data$col), "output/figures/fig1c.pdf",
  height = 7, width = 15)

# Figure 1 - 1 species
f1_data <- figure1_data_single()
data <- getAgeSequence(f1_data$Res, f1_data$ages, patchSize = c(30/50, 6), patchSpace = 0,
  scaling = c(0.75, 0.75), col = f1_data$col.hi)
plot_stand(data, size = c(1, 1, 2000, 500))
rgl.viewpoint(0, -65, fov = 0, scale = c(1, 1, 1), zoom = 1)
snapshot3d("output/figures/fig1a_1sp.png")
rgl.close()

data <- figure_1b_data(f1_data)
plot_stand(data)
tmp <- lapply(data$coords, plotPatchOutline)
rgl.viewpoint(0, -65, fov = 0, scale = c(1, 1, 1), zoom = 0.7)
snapshot3d("output/figures/fig1b_1sp.png")
rgl.close()


# Figure 2
to.pdf(figure2(), "output/figures/fig2.pdf", height = 11, width = 13)
convert_pdf_to_png("output/figures/fig2.pdf")

# Figure 3
to.pdf(figureSI(), "output/figures/fig3.pdf", height = 5, width = 8)

# Figure 4
to.pdf(figure3(), "output/figures/fig4.pdf", height = 6, width = 12)

# Figure 5
# a) panels of community average LMA and height
out <- figure5_data(f = parallel::mclapply)
to.pdf(figure5(out), "output/figures/fig5/fig5a.pdf", height = 7, width = 12)

to.pdf(figure5_single(out), "output/figures/fig5/fig5_single.pdf", height = 6, width = 7)

# b) 3D pictures for select plots
D <- c(10.6, 60, 120, 21.2, 60, 120)
Prod <- c(1, 0.6, 0.8, 0.8, 1, 1.4)
Scale <- c(60, 100, 120, 80, 100, 120)/120

for (i in seq_along(D)) {
  data <- figure_get_example_patch(D[i], Prod[i], Scale[i])
  figure_plot_patch(data)
  rgl.viewpoint(0, -65, fov = 0, scale = c(1, 1, 0.8), zoom = 1)
  filename <- file.path("output/figures/fig5", paste0("[",D[1],",",Prod[1],"]-2trait-base-T-0.png"))
  snapshot3d(filename)
  rgl.close()
  convert_png_to_transparent_png(filename)
}

# Suppmat figures --------------------------------------------------

# Leaf
to.pdf(plotleaf(), paste0("output/figures/figS-leaf.pdf"), height = 6, width = 10)

# Height
to.pdf(figure_ReproductiveAllocation(), paste0("output/figures/figS-height-a.pdf"),  height = 5, width = 12)

# 1D trait evolution
to.pdf(figS_Assembly_lma(), "output/figures/figS-1D_assembly_lma.pdf", height = 10, width = 10)
to.pdf(figS_Assembly_height("output/data/[60,1]/hsp/[0.625]/base"), "output/figures/figS-1D_assembly_hmat_high.pdf", height = 10, width = 10)
to.pdf(figS_Assembly_height("output/data/[60,1]/hsp/[0.54]/[c_r1,0.4]"), "output/figures/figS-1D_assembly_hmat_low.pdf", height = 10, width = 10)

# All fitness landscapes
to.pdf(fig_2D_landscapes_all(), "output/figures/figS_landscapes_all.pdf", height = 8, width = 6)
convert_pdf_to_png("output/figures/figS_landscapes_all.pdf")

# check sensitivity to parameter changes
to.pdf(figure_sensitivity_1(), "output/figures/figS-elasticity-pars_2D.pdf", height = 12, width = 12)
convert_pdf_to_png("output/figures/figS-elasticity-pars_2D.pdf")

to.pdf(figure_sensitivity_2(), "output/figures/figS-elasticity-repro_2D.pdf", height = 4, width = 10)
convert_pdf_to_png("output/figures/figS-elasticity-repro_2D.pdf")

to.pdf(figure_sensitivity_3(), "output/figures/figS-elasticity-repro_1D.pdf", height = 8, width = 8)

# Seed rain
to.pdf(figure_seed_rain(), paste0("output/figures/figS-seed_rain_height.pdf"), height = 6, width = 6)
