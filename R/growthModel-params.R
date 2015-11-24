# canopy shape parameters from Yokozawa et al 1995
p.eta <- 12

# ratio leaf area to sapwood area
p.theta <- 4669
p.theta <- 10000


# height - leaf mass scaling
p.a1 <- 5.44
p.B1 <- 0.306

# leaf area -stem volume scaling
p.a2 <- 6.67e-05
p.B2 <- 1.75

# root leaf scaling
p.a3 <- 0.07

# scaling of leaf turnover(/yr) to LMA tropical rate
p.a4 <- 0.0286 * 0.5
p.B4 <- 1.71

# diameter - total mass scaling
p.a5 <- 27.46
p.B5 <- 2.745

p.b <- 0.17

# nitrogen concentrations & photosynthesis leaf kg/m2
p.n_area <- 0.00187
p.c_p1 <- 150.36
p.c_p2 <- 0.19

# respiration rates mol / kg / yr
p.c_Rl <- 21000
# mol / m3 / yr
p.c_Rs <- 4012
# mol / kg / yr
p.c_Rr <- 217
p.c_Rb <- 2 * p.c_Rs

# carbon conversion parameter
p.Y <- 0.7
p.c_bio <- 0.012/0.49

# turnover
p.k_b <- 0.2
p.k_r <- 1

# REPRODUCTION
p.c_r1 <- 1
p.c_r2 <- 50
p.c_acc <- 4
