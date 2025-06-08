# Test calculation of collision integrals
# reference data: Table 4 of Kim and Monroe (2014) 

library(chapensk)

# graphics parameters
width <- 5
height <- 5
mar <- c(5.1, 4.6, 0.2, 0.6)
mai <- c(1.02, 0.82, 0.04, 0.04)
family <-"Helvetica"
pointsize <-12
lwd <- 1
tempcol <- c("blueviolet",
	"darkblue",
	"deepskyblue4",
	"cornflowerblue",
	"chartreuse4",
	"chocolate4",
	"red",
	"blue",
	"green",
	"lavender")

O11 <- CollisionIntegral(l=1,s=1)
O12 <- CollisionIntegral(l=1,s=2)
O13 <- CollisionIntegral(l=1,s=3)
O22 <- CollisionIntegral(l=2,s=2)
O23 <- CollisionIntegral(l=2,s=3)
O24 <- CollisionIntegral(l=2,s=4)
O25 <- CollisionIntegral(l=2,s=5)
O26 <- CollisionIntegral(l=2,s=6)
O44 <- CollisionIntegral(l=4,s=4)

print(round(O11$Omega(0.3),digits=4))
print(round(O12$Omega(0.3),digits=4))
print(round(O13$Omega(0.3),digits=4))
print(round(O22$Omega(0.3),digits=4))
print(round(O24$Omega(0.3),digits=4))
print(round(O25$Omega(0.3),digits=4))
print(round(O26$Omega(0.3),digits=4))
print(round(O44$Omega(0.3),digits=4))

print(round(O11$Omega(400),digits=5))
print(round(O12$Omega(400),digits=5))
print(round(O13$Omega(400),digits=5))
print(round(O22$Omega(400),digits=5))
print(round(O24$Omega(400),digits=5))
print(round(O25$Omega(400),digits=5))
print(round(O26$Omega(400),digits=5))
print(round(O44$Omega(400),digits=5))

np <- 20
lT1 <- log10(0.4)
lT2 <- log10(399)
dt <- (lT2 - lT1) / np
Theta <- 10^(lT1 + c(0:np)*dt)
if (interactive()) {
	if (exists("outputDevice")) {
		if (outputDevice == "pdf") pdf(file="man/figures/collisionintegral.pdf",width=width,height=height,family=family,pointsize=pointsize)
		if (outputDevice == "svg") svg(file="man/figures/collisionintegral.svg",width=width,height=height,family=family,pointsize=pointsize)
	}
}
par(mai=mai,mar=mar)
plot(Theta,O11$Omega(Theta),log="xy",type="l",xlab=expression(Theta),ylab=expression(Omega(l,s)),col=tempcol)
lines(Theta,O12$Omega(Theta),col=tempcol[2])
lines(Theta,O13$Omega(Theta),col=tempcol[3])
lines(Theta,O22$Omega(Theta),col=tempcol[4])
lines(Theta,O24$Omega(Theta),col=tempcol[5])
lines(Theta,O26$Omega(Theta),col=tempcol[6])
lines(Theta,O44$Omega(Theta),col=tempcol[7])
legend("topright",legend=c("(1,2)","(1,3)","(2,2)","(2,4)","(2,6)","(4.4)"), text.col=tempcol, col=tempcol, box.lwd=0)
dev.off()