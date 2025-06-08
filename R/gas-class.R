Gas <- setRefClass("Gas",
	fields = list(
		name="character",
		M="numeric",
		m="numeric",
		sigma="numeric",
		zeta="numeric",
		epsk="numeric",
		dipole_moment="numeric",
		polarizability="numeric",
		O11="ANY",
		O12="ANY",
		O13="ANY",
		O22="ANY",
		O23="ANY",
		O24="ANY",
		O25="ANY",
		O26="ANY",
		O44="ANY",
		bathGas="ANY"
	),
	methods = list(

	   initialize = function (name="") {
	   	name1 <- name
	   	O11 <<- CollisionIntegral(l=1,s=1)
	   	O12 <<- CollisionIntegral(l=1,s=2)  
	   	O13 <<- CollisionIntegral(l=1,s=3)  	
	   	O22 <<- CollisionIntegral(l=2,s=2)
	   	O23 <<- CollisionIntegral(l=2,s=3)
	   	O24 <<- CollisionIntegral(l=2,s=4)
	   	O25 <<- CollisionIntegral(l=2,s=5)
	   	O26 <<- CollisionIntegral(l=2,s=6)
	   	O44 <<- CollisionIntegral(l=4,s=4)
	   	gas_data <- subset(gas,gas$name==name1)
	   	if (nrow(gas_data)==1)
	   	{
				name <<- name
				M <<- gas_data$M
				m <<- M / 1000 / pkg.env$Na
				sigma <<- gas_data$sigma
				epsk <<- gas_data$epsk
				zeta <<- 0
				polarizability <<- gas_data$polarizability
				dipole_moment <<- gas_data$dipole_moment  	
	   	}
	   },
    	
		B = function (T, sigma=.self$sigma, epsk=.self$epsk, zeta=.self$zeta) {
		  xsigma <- 1e-10*(sigma + T*zeta)
		  xepsk <- epsk*(1e-10*sigma/xsigma)^6			
			Theta <- T/xepsk
			u <- 1/(2*Theta)
			return <- pkg.env$Na * (2/3)*u*sqrt(2)*pi^2*(xsigma)^3*exp(u)*(
				BesselI(u,-3/4) + BesselI(u,3/4) - BesselI(u,1/4) - BesselI(u,-1/4)						
			)
		},   	
		
		fit_B_data = function(B_df) {
			plot(B_df$T, 
			  1E6*B_df$value,
			  xlab=expression(paste(italic(T), " / [K]")),
 				ylab=expression(paste(italic(B), " / [cm3/mol]")),
 				main="")
			if (zeta == 0) {
				nls1 <- nls(value ~ .self$B(T, sigma, epsk, zeta=0),  
				data=B_df,trace=TRUE,start=list(sigma=.self$sigma, epsk=.self$epsk))
				print(summary(nls1))	
				sigma <<- coef(nls1)["sigma"]
				names(sigma) <<- NULL
				epsk <<- coef(nls1)["epsk"]
				names(epsk) <<- NULL	
				lines(B_df$T, 1E6*.self$B(B_df$T, sigma, epsk, 0))	
			} else {
				nls1 <- nls(value ~ .self$B(T, sigma, epsk, zeta),  
				data=B_df,trace=TRUE,start=list(sigma=.self$sigma, epsk=.self$epsk, zeta=.self$zeta))
				print(summary(nls1))	
				sigma <<- coef(nls1)["sigma"]
				names(sigma) <<- NULL
				epsk <<- coef(nls1)["epsk"]
				names(epsk) <<- NULL
				zeta <<- coef(nls1)["zeta"]	
				names(zeta) <<- NULL
				lines(B_df$T, 1E6*.self$B(B_df$T, sigma, epsk, zeta))	
			}
			return <- nls1
		},
		 	   	
		viscosity = function(T, sigma=.self$sigma, epsk=.self$epsk, zeta=.self$zeta, third_order_correction=TRUE) {
		
			f3_viscosity = function(Theta) {
	      stopifnot(length(Theta) == 1)
	    	Omega22 <- O22$Omega(Theta)
				Omega23 <- O23$Omega(Theta)
				Omega24 <- O24$Omega(Theta)
				Omega25 <- O25$Omega(Theta)
				Omega26 <- O26$Omega(Theta)
				Omega44 <- O44$Omega(Theta)
				b <- matrix(ncol=3,nrow=3)
				b[1,1] <- 4*Omega22
				b[1,2] <- 7*Omega22 - 8*Omega23
				b[2,1] <- b[1,2]
				b[2,2] <- (301/12)*Omega22 - 28*Omega23 + 20*Omega24
				b[1,3] <- (63/8)*Omega22 - 18*Omega23 + 10*Omega24
				b[3,1] <- b[1,3]
				b[2,3] <- (1365/32)*Omega22 - (321/4)*Omega23 + (125/2)*Omega24 - 30*Omega25
				b[3,2] <- b[2,3]
				b[3,3] <- (25137/256)*Omega22 - (1755/8)*Omega23 + (1905/8)*Omega24 - 135*Omega25 + (105/2)*Omega26 + 12*Omega44
				y <- b[1,2]^2/(b[1,1]*b[2,2] - b[1,2]^2) + b[1,1]*(b[1,2]*b[2,3]-b[2,2]*b[1,3])^2/((b[1,1]*b[2,2]-b[1,2]^2)*det(b))
				return <- 1 + y
    	}
	
			xsigma <- 1e-10*(sigma + T*zeta)	
			xepsk <- epsk*(1e-10*sigma/xsigma)^6			
			Theta <- T/xepsk
		  if (third_order_correction) {
				f3v <- Vectorize(f3_viscosity,c("Theta"))
				f3 <- f3v(Theta)
 			} else {f3 <- 1}
			return <- f3*5*sqrt(pi*pkg.env$k*T*m) / (16*pi*(xsigma)^2*O22$Omega(Theta))
		},
		
		fit_viscosity_data = function(viscosity_df) {
			plot(viscosity_df$T, 
			  1E6*viscosity_df$value,
			  xlab=expression(paste(italic(T), " / [K]")),
 				ylab=expression(paste(italic(eta), " / [uPa.s]")),
 				main="")
			if (.self$zeta == 0) {
				nls1 <- nls(value ~ .self$viscosity(T, sigma, epsk, zeta=0, 
				third_order_correction=TRUE),  
				data=viscosity_df,trace=TRUE,start=list(sigma=.self$sigma, epsk=.self$epsk))
				print(summary(nls1))	
				sigma <<- coef(nls1)["sigma"]
				names(sigma) <<- NULL
				epsk <<- coef(nls1)["epsk"]
				names(epsk) <<- NULL	
				lines(viscosity_df$T, 1E6*.self$viscosity(viscosity_df$T, sigma, epsk, 0, third_order_correction=TRUE))	
			} else {
				nls1 <- nls(value ~ .self$viscosity(T, sigma, epsk, .self$zeta, 
				third_order_correction=TRUE),  
				data=viscosity_df,trace=TRUE,start=list(sigma=.self$sigma, epsk=.self$epsk, zeta=.self$zeta))
				print(summary(nls1))	
				sigma <<- coef(nls1)["sigma"]
				epsk <<- coef(nls1)["epsk"]	
				lines(viscosity_df$T, 1E6*.self$viscosity(viscosity_df$T, sigma, epsk, zeta, third_order_correction=TRUE))	
			}
			return <- nls1
		},
		
		fit_B_viscosity_data = function(B_df, viscosity_df,log=FALSE) {
			scaleP2 <- 100
			scaleP3 <- 0.0001
			Tmin <- min(min(B_df$T),min(viscosity_df$T))
			delta_b2v <- abs( max(B_df$value) - min(B_df$value) )
  		delta_viscosity <- abs( max(viscosity_df$value) - min(viscosity_df$value))
      if (.self$zeta != 0) {
	      if (log == TRUE) {
	  			chisq <- function(p) {	
	  			  chib <- sum( 
	  					(log(abs(B_df$value) 
	  						/ (abs(.self$B(B_df$T,p[1],p[2],p[3])))))^2 
	  				)
	  				chiv <- sum( 
	  				  (log(viscosity_df$value 
	  				  	/ .self$viscosity(viscosity_df$T,p[1],p[2],p[3])))^2
	  				)
	  				return <- chib + chiv
		  		}	
	 			}	else {
	 				chisq <- function(p) {
	 				  chib <- sum(
	 						((abs(.self$B(B_df$T,p[1],p[2],p[3])) 
	 							- abs(B_df$value))/delta_b2v)^2
	 					)
						chiv <- 
						sum(abs(B_df$value)
							((.self$viscosity(viscosity_df$T,p[1],p[2],p[3]) 
								- viscosity_df$value)/delta_viscosity)^2
						)					
						return <- chib + chiv		 	
	 				}
	 			}
	 			# nlm1 <- nlm(chisq,p=c(.self$sigma,.self$epsk/100, .self$zeta*1E3),steptol=1e-12,hessian=TRUE,stepmax=1e-6)
	 			nlm1 <- optim(par=c(.self$sigma,.self$epsk,.self$zeta),chisq,method="L-BFGS-B",lower=c(0,Tmin/0.3,0),control=list(trace=1,parscale=c(1,scaleP2,scaleP3)))
	 			sigma <<- nlm1$par[1]
	 			epsk <<- nlm1$par[2]
	 			zeta <<- nlm1$par[3]
	 			return <- nlm1
	 			} else {
	 			  if (log == TRUE) {
	  			chisq <- function(p) {	
	  			  chib <- sum( 
	  					(log(abs(B_df$value) 
	  						/ (abs(.self$B(B_df$T,p[1],p[2],zeta=0)))))^2 
	  				)
	  				chiv <- sum( 
	  				  (log(viscosity_df$value 
	  				  	/ .self$viscosity(viscosity_df$T,p[1],p[2],zeta=0)))^2
	  				)
	  				return <- chib + chiv
		  		}	
	 			}	else {
	 				chisq <- function(p) {	
	 				  chib <- sum(
	 						((abs(.self$B(B_df$T,p[1],p[2],zeta=0)) 
	 							- abs(B_df$value))/delta_b2v)^2
	 					)
						chiv <- 
						sum(
							((.self$viscosity(viscosity_df$T,p[1],p[2],zeta=0) 
								- viscosity_df$value)/delta_viscosity)^2
						)					
						return <- chib + chiv		 	
	 				}
	 			}
	 			# nlm1 <- nlm(chisq,p=c(.self$sigma ,.self$epsk, .self$sigmat),steptol=1e-12,hessian=TRUE,stepmax=1e-6)
	 			# nlm1 <- nlm(chisq,p=c(.self$sigma, .self$epsk/100, 0),steptol=1e-12,hessian=TRUE)
	 			nlm1 <- optim(par=c(.self$sigma,.self$epsk),chisq,method="L-BFGS-B",control=list(trace=0,parscale=c(1,scaleP2)))
	 			sigma <<- nlm1$par[1]
	 			epsk <<- nlm1$par[2]
	 			zeta <<- 0
	 			return <- nlm1
	 		}		 		
		},
		
		thermal_conductivity = function(T, third_order_correction=TRUE) {
		
			f3_thermal_conductivity = function(Theta) {
	      stopifnot(length(Theta) == 1)
	    	Omega22 <- O22$Omega(Theta)
				Omega23 <- O23$Omega(Theta)
				Omega24 <- O24$Omega(Theta)
				Omega25 <- O25$Omega(Theta)
				Omega26 <- O26$Omega(Theta)
				Omega44 <- O44$Omega(Theta)
				a <- matrix(ncol=3,nrow=3)
				a[1,1] <- 4*Omega22
				a[1,2] <- 7*Omega22 - 8*Omega23
				a[2,1] <- a[1,2]
				a[2,2] <- (77/4)*Omega22 - 28*Omega23 + 20*Omega24
				a[1,3] <- (63/8)*Omega22 - 18*Omega23 + 10*Omega24
				a[3,1] <- a[1,3]
				a[2,3] <- (945/32)*Omega22 - (261/4)*Omega23 + (125/2)*Omega24 - 30*Omega25
				a[3,2] <- a[2,3]
				a[3,3] <- (14533/256)*Omega22 - (1215/8)*Omega23 + (1565/8)*Omega24 - 135*Omega25 + (105/2)*Omega26 + 12*Omega44
				y <- a[1,2]^2/(a[1,1]*a[2,2] - a[1,2]^2) + a[1,1]*(a[1,2]*a[2,3]-a[2,2]*a[1,3])^2/((a[1,1]*a[2,2]-a[1,2]^2)*det(a))
				return <- 1 + y
    	}
	
			xsigma <- 1e-10*(sigma + T*zeta)	
			xepsk <- epsk*(1e-10*sigma/xsigma)^6			
			Theta <- T/xepsk
		  if (third_order_correction) {
				f3v <- Vectorize(f3_thermal_conductivity,c("Theta"))
				f3 <- f3v(Theta)
 			} else {f3 <- 1}
			return <- f3*75*pkg.env$k*sqrt(pkg.env$k*T/(pi*m)) / (64*(xsigma)^2*O22$Omega(Theta))
		},
		
    density = function(p=pkg.env$p0, T=pkg.env$T0) {
    	B <- .self$B(T)
    	return <- ifelse(B < 1E-10, p*.self$M/(1000*pkg.env$R*T), .self$M / (2*1000*B) * (sqrt(1 + 4*p*B/(pkg.env$R*T)) - 1))	
    },

		diffusion = function(p=pkg.env$p0,T=pkg.env$T0,second_order_correction=TRUE) {
		
			f2_diffusion <- function(Theta) {
				Omega11 <- O11$Omega(Theta)
				Omega22 <- O22$Omega(Theta)
				Omega12 <- O12$Omega(Theta)
				Omega13 <- O13$Omega(Theta)
				A <- Omega22 / Omega11
				B <- (5*Omega12 - 4*Omega13)	/ Omega11
				C <- Omega12 / Omega11 		
			}		
		
			stopifnot(!is.null(.self$sigma))
			stopifnot(!is.null(.self$epsk))	
			stopifnot(!is.null(.self$m)) 
			xsigma <- 1e-10*(sigma + T*zeta)	
			xepsk <- epsk*(1e-10*sigma/xsigma)^6			
			Theta <- T/xepsk
			Omega11 <- O11$Omega(Theta)
			if (second_order_correction) {
				f2 <- f2_diffusion(Theta)
			} else {
				f2 <- 1			
			}
			return <- 3*sqrt(pi*m*pkg.env$k*T) / (8*pi*xsigma^2*density(p,T)*Omega11)	* f2
		
		},
		
		binary_diffusion = function(p=pkg.env$p0, T=pkg.env$T0, bathGas) {
		  stopifnot(class(bathGas)[1]=="Gas") 
			stopifnot(!is.null(.self$sigma))
			stopifnot(!is.null(.self$epsk))	
			stopifnot(!is.null(.self$m)) 
			stopifnot(!is.null(bathGas$sigma))
			stopifnot(!is.null(bathGas$epsk))
			stopifnot(!is.null(bathGas$m)) 
			sigma1 <- 1e-10 * (.self$sigma + T*.self$zeta)
			sigma2 <- 1e-10 * (bathGas$sigma + T*bathGas$zeta)
			eps1 <- pkg.env$k*.self$epsk*(1e-10*.self$sigma/sigma1)^6
		  eps2 <- pkg.env$k*bathGas$epsk*(1e-10*bathGas$sigma/sigma2)^6
			dipole_moment1 <- pkg.env$Debye*.self$dipole_moment
			r_dipole_moment2 <- dipole_moment1^2 / (eps1 * sigma1^3)
			r_polarizability <- 1E-30 * bathGas$polarizability / sigma2^3
			xi <- 1 + r_polarizability*r_dipole_moment2 / (16*pi*pkg.env$eps0) * sqrt(eps1/eps2)
			eps <- xi^2 * sqrt(eps1 * eps2)
			epsk2 <- eps / pkg.env$k
			xsigma <- xi^(-1/6)*(sigma1 + sigma2) / 2
			mu <- .self$m * bathGas$m / (.self$m + bathGas$m)
		  Theta <- T/epsk2
			return <- (3/16)*sqrt(2*pi*pkg.env$k*T/mu) * pkg.env$k*T/(pi*xsigma^2*O11$Omega(Theta)*p)	
		}
				
	)
)	