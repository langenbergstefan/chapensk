library(chapensk)

# graphics parameters
width <- 10
height <- 5
mar <- c(5.1, 4.6, 0.2, 0.6)
mai <- c(1.02, 0.82, 0.04, 0.04)
family <-"Helvetica"
pointsize <-12
lwd <- 1

if (interactive()) {
	if (exists("outputDevice")) {
		if (outputDevice == "pdf") pdf(file="man/figures/critical_data.pdf",width=width,height=height,family=family,pointsize=pointsize)
		if (outputDevice == "svg") svg(file="man/figures/critical_data.svg",width=width,height=height,family=family,pointsize=pointsize)
}
 split.screen(figs=c(1,2))
}

if (interactive()) {screen(1)}
par(mai=mai,mar=mar)
np <- subset(gas, dipole_moment==0)
plot(np$epsk, np$Tc, xlim=c(0,600),ylim=c(0,700),
	xlab=expression(epsilon/italic(k) ~ "[K]"), 
	ylab="Critical temperature / [K]",
)
text(np$epsk,  np$Tc, np$name, pos=4, cex=0.7, col="blue")
abline(a=0,b=1.321)
if (interactive()) {
	legend("top", bty="n", legend="(a)")	
	screen(2)
}
x <- 0.001/pkg.env$Na*((1e-10*np$sigma)^(-3))
y <- np$rhoc 
par(mai=mai,mar=mar)
plot(x, np$rhoc, 
xlim=c(0,60), 
ylim=c(0,16),
  xlab=expression(1/(sigma^3)/Na ~ "[mol/l]"), 
  ylab="Critical density / [mol/l]",
)
text(x, np$rhoc, np$name, pos=4, cex=0.7, col="blue")
abline(a=0,0.316)

if (interactive()) {
  legend("top", bty="n", legend="(b)")
	close.screen()
	dev.off()
}	