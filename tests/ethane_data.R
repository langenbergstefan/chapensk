library(chapensk)

# graphics parameters
width <- 10
height <- 5
mar <- c(5.1, 4.6, 0.2, 0.6)
mai <- c(1.02, 0.82, 0.04, 0.04)
family <-"Helvetica"
pointsize <-12
lwd <- 1

c2h6 <- Gas("ethane")

if (interactive()) {
if (exists("outputDevice")) {
	if (outputDevice == "pdf") pdf(file="man/figures/ethane.pdf",width=width,height=height,family=family,pointsize=pointsize)
	if (outputDevice == "svg") svg(file="man/figures/ethane.svg",width=width,height=height,family=family,pointsize=pointsize)
}
 split.screen(figs=c(1,2))
}

ethane_viscosity <- subset(ethane_data, (property=="viscosity") & (T>=90) & (T<=675), select=c(T, value))
sigma0 <- c2h6$sigma
epsk0 <- c2h6$epsk
c2h6$zeta <- 0
print("Literature data")
print(c2h6$sigma)
print(c2h6$epsk)
ethane_viscosity$value <- 1E-6*ethane_viscosity$value
ethane_B <- subset(ethane_data, (property=="BQFH") & (T>= 220) & (T<=623))
ethane_B$value <- 1E-6*ethane_B$value
fit <- c2h6$fit_B_viscosity_data(ethane_B,ethane_viscosity,log=FALSE)
sigma1 <- c2h6$sigma
epsk1 <- c2h6$epsk
zeta1 <- c2h6$zeta
print("Fit of B and viscosity data")
print(c2h6$sigma)
print(c2h6$epsk)
print(c2h6$zeta)
stopifnot(round(c2h6$sigma,2) == 4.36)
stopifnot(round(c2h6$epsk,1) == 244.8)
c2h6$sigma <- sigma0
c2h6$epsk <- epsk0
c2h6$zeta <- 0
if (interactive()) {screen(1)}
par(mai=mai,mar=mar)
c2h6$fit_viscosity_data(ethane_viscosity)
lines(ethane_viscosity$T, 1E6*c2h6$viscosity(ethane_viscosity$T,sigma=sigma1, epsk=epsk1, zeta=zeta1), col="red", lty=2)
if (interactive()) {
	legend("top", bty="n", legend="(a)")	
	screen(2)
}
print("Fit of viscosity data")
print(c2h6$sigma)
print(c2h6$epsk)
stopifnot(round(c2h6$sigma,2) == 4.38)
stopifnot(round(c2h6$epsk,1) == 235.7)
par(mai=mai,mar=mar)
c2h6$zeta <- 0
c2h6$fit_B_data(ethane_B)
lines(ethane_B$T, 1E6*c2h6$B(ethane_B$T,sigma=sigma1, epsk=epsk1), col="red", lty=2)
if (interactive()) {
  legend("top", bty="n", legend="(b)")
	close.screen()
	dev.off()
}	
print("Fit of B data")
print(c2h6$sigma)
print(c2h6$epsk)
stopifnot(round(c2h6$sigma,2) == 4.95)
stopifnot(round(c2h6$epsk,1) == 202.0)





