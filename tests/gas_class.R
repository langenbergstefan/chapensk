library(chapensk)

CH4 <- Gas("methane")
Ar <- Gas("argon")
# virial coefficient 
stopifnot(round(1E6*CH4$B(T=300),0) == -41)
# viscosity
stopifnot(round(1E6*CH4$viscosity(T=300),1) == 11.2)
# density 
stopifnot(round(CH4$density(T=300),3) == 0.652)
# diffusion coefficient
stopifnot(round(1E4*CH4$diffusion(T=300),3) == 0.199)
# thermal conductivity
stopifnot(round(Ar$thermal_conductivity(T=300),3)==0.018)


