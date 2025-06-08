library(chapensk)

# graphics parameters
width <- 10
height <- 10
mar <- c(5.1, 4.6, 0.2, 0.6)
mai <- c(1.02, 0.82, 0.04, 0.04)
family <-"Helvetica"
pointsize <-12
lwd <- 1

Diffusion <- function(T, D0, b) {
	return <- D0*(T/pkg.env$T0)^b
}

helium <- Gas("helium")
nitrogen <- Gas("nitrogen")
argon <- Gas("argon")
methane <- Gas("methane")
ethane <- Gas("ethane")
propane <- Gas("propane")
butane <- Gas("n-butane")
fluoromethane <- Gas("fluoromethane")
difluoromethane <- Gas("difluoromethane")
trifluoromethane <- Gas("trifluoromethane")

if (interactive()) {
	if (exists("outputDevice")) {
  	if (outputDevice == "pdf") pdf(file="man/figures/binary_diffusion.pdf",width=width,height=height,family=family,pointsize=pointsize)
		if (outputDevice == "svg") svg(file="man/figures/binary_diffusion.svg",width=width,height=height,family=family,pointsize=pointsize)
}
 split.screen(figs=c(2,2))
}

nitrogen_in_helium <- subset(binary_diffusion,(gas=="nitrogen" & bath_gas=="helium"))
argon_in_helium <- subset(binary_diffusion,(gas=="argon" & bath_gas=="helium"))

if (interactive()) {screen(1)}
par(mai=mai,mar=mar)

plot(nitrogen_in_helium$T, nitrogen_in_helium$D,
	log="xy",
	xlab=expression(paste(italic(T), " / [K]")), 
	ylab=expression(paste(italic(D), "/ [cm2/s]")),
	pch=0,col="blue")
Dcalc <- nitrogen$binary_diffusion(T=nitrogen_in_helium$T, bathGas=helium)
nitrogen_in_helium$Dcalc <- 1E4*Dcalc
lines(nitrogen_in_helium$T,1E4*Dcalc,col="blue")
nls_n2_he <- nls(D ~ Diffusion(T, D0, b), data=nitrogen_in_helium, weights=1/U_D, start=list(D0=0.62, b=1.7), trace=TRUE)
print(summary(nls_n2_he))
points(argon_in_helium$T, argon_in_helium$D,pch=1,col="red")
Dcalc <- argon$binary_diffusion(T=argon_in_helium$T, bathGas=helium)
argon_in_helium$Dcalc <- 1E4*Dcalc
lines(argon_in_helium$T,1E4*Dcalc,col="red")
nls_ar_he <- nls(D ~ Diffusion(T, D0, b), data=argon_in_helium, weights=1/U_D, start=list(D0=0.62, b=1.7), trace=TRUE)
print(summary(nls_ar_he))
legend("bottomright",
	legend=c("nitrogen","argon"),
	col=c("blue","red"),
	pch=c(0,1)
)
if (interactive()) {
  legend("top", bty="n", legend="(a)")
	close.screen()
}	

methane_in_helium <- subset(binary_diffusion,(gas=="methane" & bath_gas=="helium"))
ethane_in_helium <- subset(binary_diffusion,(gas=="ethane" & bath_gas=="helium"))
propane_in_helium <- subset(binary_diffusion,(gas=="propane" & bath_gas=="helium"))
butane_in_helium <- subset(binary_diffusion,(gas=="butane" & bath_gas=="helium"))

if (interactive()) {screen(2)}
par(mai=mai,mar=mar)

plot(methane_in_helium$T, methane_in_helium$D,
	log="xy",
	xlab=expression(paste(italic(T), " / [K]")), 
	ylab=expression(paste(italic(D), "/ [cm2/s]")),
	pch=0,col="blue",
	ylim=c(0.3,3))
Dcalc <- methane$binary_diffusion(T=nitrogen_in_helium$T, bathGas=helium)
lines(methane_in_helium$T,1E4*Dcalc,col="blue")
methane_in_helium$Dcalc <- 1E4*Dcalc
nls_ch4_he <- nls(D ~ Diffusion(T, D0, b), data=methane_in_helium, weights=1/U_D, start=list(D0=0.62, b=1.7), trace=TRUE)
print(summary(nls_ch4_he))
points(ethane_in_helium$T, ethane_in_helium$D,pch=1,col="red")
Dcalc <- ethane$binary_diffusion(T=nitrogen_in_helium$T, bathGas=helium)
lines(ethane_in_helium$T,1E4*Dcalc,col="red")
ethane_in_helium$Dcalc <- 1E4*Dcalc
nls_c2h6_he <- nls(D ~ Diffusion(T, D0, b), data=ethane_in_helium, start=list(D0=0.62, b=1.7), trace=TRUE)
print(summary(nls_c2h6_he))
points(propane_in_helium$T, propane_in_helium$D,pch=2,col="dark green")
Dcalc <- propane$binary_diffusion(T=propane_in_helium$T, bathGas=helium)
lines(propane_in_helium$T,1E4*Dcalc,col="dark green")
propane_in_helium$Dcalc <- 1E4*Dcalc
nls_c3h8_he <- nls(D ~ Diffusion(T, D0, b), data=propane_in_helium, start=list(D0=0.62, b=1.7), trace=TRUE)
print(summary(nls_c3h8_he))
points(butane_in_helium$T, butane_in_helium$D,pch=3,col="orange")
Dcalc <- butane$binary_diffusion(T=butane_in_helium$T, bathGas=helium)
lines(butane_in_helium$T,1E4*Dcalc,col="orange")
butane_in_helium$Dcalc <- 1E4*Dcalc
nls_c4h10_he <- nls(D ~ Diffusion(T, D0, b), data=butane_in_helium, start=list(D0=0.62, b=1.7), trace=TRUE)
print(summary(nls_c4h10_he))
legend("bottomright",
	legend=c("methane","ethane","propane","butane"),
	col=c("blue","red","dark green","orange"),
	pch=c(0,1,2,3)
)
if (interactive()) {
  legend("top", bty="n", legend="(b)")
	close.screen()
}

if (interactive()) {screen(3)}
par(mai=mai,mar=mar)
methane_in_nitrogen <- subset(binary_diffusion,(gas=="methane" & bath_gas=="nitrogen"))
ethane_in_nitrogen <- subset(binary_diffusion,(gas=="ethane" & bath_gas=="nitrogen"))
propane_in_nitrogen <- subset(binary_diffusion,(gas=="propane" & bath_gas=="nitrogen"))
butane_in_nitrogen <- subset(binary_diffusion,(gas=="butane" & bath_gas=="nitrogen"))
plot(methane_in_nitrogen$T, methane_in_nitrogen$D,
	log="xy",
	xlab=expression(paste(italic(T), " / [K]")), 
	ylab=expression(paste(italic(D), "/ [cm2/s]")),
	pch=0,col="blue",
	ylim=c(0.1,1))
Dcalc <- methane$binary_diffusion(T=methane_in_nitrogen$T, bathGas=nitrogen)
methane_in_nitrogen$Dcalc <- 1E4*Dcalc
lines(methane_in_nitrogen$T,1E4*Dcalc,col="blue")
nls_ch4_n2 <- nls(D ~ Diffusion(T, D0, b), data=methane_in_nitrogen, weights=1/U_D, start=list(D0=0.62, b=1.7), trace=TRUE)
print(summary(nls_ch4_n2))
points(ethane_in_nitrogen$T, ethane_in_nitrogen$D,pch=1,col="red")
Dcalc <- ethane$binary_diffusion(T=ethane_in_nitrogen$T, bathGas=nitrogen)
lines(ethane_in_nitrogen$T,1E4*Dcalc,col="red")
ethane_in_nitrogen$Dcalc <- 1E4*Dcalc
nls_c2h6_n2 <- nls(D ~ Diffusion(T, D0, b), data=ethane_in_nitrogen, weights=1/U_D, start=list(D0=0.62, b=1.7), trace=TRUE)
print(summary(nls_c2h6_n2))
points(propane_in_nitrogen$T, propane_in_nitrogen$D,pch=2,col="dark green")
Dcalc <- propane$binary_diffusion(T=propane_in_nitrogen$T, bathGas=nitrogen)
lines(propane_in_nitrogen$T,1E4*Dcalc,col="dark green")
propane_in_nitrogen$Dcalc <- 1E4*Dcalc
nls_c3h8_n2 <- nls(D ~ Diffusion(T, D0, b), data=propane_in_nitrogen, weights=1/U_D, start=list(D0=0.62, b=1.7), trace=TRUE)
print(summary(nls_c3h8_n2))
points(butane_in_nitrogen$T, butane_in_nitrogen$D,pch=3,col="orange")
Dcalc <- butane$binary_diffusion(T=butane_in_nitrogen$T, bathGas=nitrogen)
lines(butane_in_nitrogen$T,1E4*Dcalc,col="orange")
butane_in_nitrogen$Dcalc <- 1E4*Dcalc
nls_c4h10_n2 <- nls(D ~ Diffusion(T, D0, b), data=butane_in_nitrogen, weights=1/U_D, start=list(D0=0.62, b=1.7), trace=TRUE)
print(summary(nls_c4h10_n2))
legend("topleft",
	legend=c("methane","ethane","propane","butane"),
	col=c("blue","red","dark green","orange"),
	pch=c(0,1,2,3)
)
if (interactive()) {
  legend("top", bty="n", legend="(c)")
	close.screen()
}

if (interactive()) {screen(4)}
par(mai=mai,mar=mar)
fluoromethane_in_nitrogen <- subset(binary_diffusion,(gas=="fluoromethane" & bath_gas=="nitrogen"))
difluoromethane_in_nitrogen <- subset(binary_diffusion,(gas=="difluoromethane" & bath_gas=="nitrogen"))
trifluoromethane_in_nitrogen <- subset(binary_diffusion,(gas=="trifluoromethane" & bath_gas=="nitrogen"))
plot(fluoromethane_in_nitrogen$T, fluoromethane_in_nitrogen$D,
	log="xy",
	xlab=expression(paste(italic(T), " / [K]")), 
	ylab=expression(paste(italic(D), "/ [cm2/s]")),
	pch=0,col="blue",
	ylim=c(0.1,0.6))
fluoromethane_chisq = function(p) {
	fluoromethane$sigma <- p[1]
	fluoromethane$zeta <- 0
	fluoromethane$epsk <- p[2]
	Dcalc <- 1E4*fluoromethane$binary_diffusion(T=fluoromethane_in_nitrogen$T, bathGas=nitrogen)
	return <- sum(((fluoromethane_in_nitrogen$D-Dcalc)/fluoromethane_in_nitrogen$U_D)^2)
}
optim_fluoromethane <- optim(par=c(4,180),fluoromethane_chisq,method="L-BFGS-B",control=list(trace=1,parscale=c(1,100)))
print(optim_fluoromethane$par[1]) 
print(optim_fluoromethane$par[2]) 
fluoromethane$sigma <- optim_fluoromethane$par[1]
fluoromethane$zeta <- 0
fluoromethane$epsk <- optim_fluoromethane$par[2]
Dcalc <- fluoromethane$binary_diffusion(T=fluoromethane_in_nitrogen$T, bathGas=nitrogen)
lines(fluoromethane_in_nitrogen$T,1E4*Dcalc,col="blue")
fluoromethane_in_nitrogen$Dcalc <- 1E4*Dcalc

points(difluoromethane_in_nitrogen$T, difluoromethane_in_nitrogen$D,pch=1,col="red")
difluoromethane_chisq = function(p) {
	difluoromethane$sigma <- p[1]
	difluoromethane$zeta <- 0
	difluoromethane$epsk <- p[2]
	Dcalc <- 1E4*difluoromethane$binary_diffusion(T=difluoromethane_in_nitrogen$T, bathGas=nitrogen)
	return <- sum(((difluoromethane_in_nitrogen$D-Dcalc)/difluoromethane_in_nitrogen$U_D)^2)
}
optim_difluoromethane <- optim(par=c(4,180),difluoromethane_chisq,method="L-BFGS-B",control=list(trace=1,parscale=c(1,100)))
print(optim_difluoromethane$par[1]) 
print(optim_difluoromethane$par[2])
difluoromethane$sigma <- optim_difluoromethane$par[1]
difluoromethane$zeta <- 0
difluoromethane$epsk <- optim_difluoromethane$par[2]
Dcalc <- difluoromethane$binary_diffusion(T=difluoromethane_in_nitrogen$T, bathGas=nitrogen)
lines(difluoromethane_in_nitrogen$T,1E4*Dcalc,col="red")
difluoromethane_in_nitrogen$Dcalc <- 1E4*Dcalc

points(trifluoromethane_in_nitrogen$T, trifluoromethane_in_nitrogen$D,pch=2,col="dark green")
trifluoromethane_chisq = function(p) {
	trifluoromethane$sigma <- p[1]
	trifluoromethane$zeta <- 0
	trifluoromethane$epsk <- p[2]
	Dcalc <- 1E4*trifluoromethane$binary_diffusion(T=trifluoromethane_in_nitrogen$T, bathGas=nitrogen)
	return <- sum(((trifluoromethane_in_nitrogen$D-Dcalc)/trifluoromethane_in_nitrogen$U_D)^2)
}
optim_trifluoromethane <- optim(par=c(4,180),trifluoromethane_chisq,method="L-BFGS-B",control=list(trace=1,parscale=c(1,100)))
print(optim_trifluoromethane$par[1]) 
print(optim_trifluoromethane$par[2]) 
trifluoromethane$sigma <- optim_trifluoromethane$par[1]
trifluoromethane$zeta <- 0
trifluoromethane$epsk <- optim_trifluoromethane$par[2]
Dcalc <- trifluoromethane$binary_diffusion(T=trifluoromethane_in_nitrogen$T, bathGas=nitrogen)
lines(trifluoromethane_in_nitrogen$T,1E4*Dcalc,col="dark green")
trifluoromethane_in_nitrogen$Dcalc <- 1E4*Dcalc
legend("bottomright",
	legend=c("fluoromethane","difluoromethane","trifluoromethane"),
	col=c("blue","red","dark green"),
	pch=c(0,1,2)
)

if (interactive()) {
  legend("top", bty="n", legend="(d)")
	close.screen()
}

dev.off()
