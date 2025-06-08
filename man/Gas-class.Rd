\name{Gas-class}
\Rdversion{1.1}
\docType{class}
\alias{Gas-class}
\alias{Gas}

\title{Class \code{"Gas"}}
\description{
Reference class for gas with methods for the estimation of transport properties. The equations for the higher order corrections are taken from Kim and Monroe 2014. The second virial coefficient
is calculated using the equation given by Vargas \Remph{et al.} 2001. 
\if{text}{See the HTML and PDF versions of this document for details.}
}
\references{
   Bechtel S, Bayer B, Vidakovic-Koch T, Wiser A, Vogel H, Sundmacher K. Precise determination of LJ parameters and Eucken correction factors for a more accurate modeling 
   of transport properties in gases. \Remph{Heat and Mass Transfer} 2020;56:2515-27. \doi{10.1007/s00231-020-02871-4}.
  	
	Brown NJ, Bastien LAJ, Price PN. Transport properties for combustion modeling. \Remph{Progress in Energy and Combustion Science}
   2011;37:56582. \doi{10.1016/j.pecs.2010.12.001}.
   
   Kim SU, Monroe CW. High-accuracy calculations of sixteen collision integrals for Lennard-Jones 
   (12-6) gases and their interpolation to parameterize neon, argon, and krypton. 
   \Remph{Journal of Computational Physics 2014} 273:358-73, 
   \doi{10.1016/j.jcp.2014.05.018}.
  
  	Marrero TR, Mason EA. Gaseous Diffusion Coefficients. \Remph{J. Phys. Chem. Ref. Data} 1972;1:3-118.
  	\doi{10.1063/1.3253094}.
	
  	Vargas P, Munoz E, Rodriguez L. Second virial coefficient for the Lennard-Jones potential. 
  	\Remph{Physica A: Statistical Mechanics and Its Applications} 2001;290:92-100. \doi{10.1016/s0378-4371(00)00362-9}.
  	
  	Zarkova L, Hohm U. Effective (n-6) Lennard-Jones Potentials with Temperature-Dependent Parameters Introduced for Accurate Calculation of Equilibrium and Transport 
  	Properties of Ethene, Propene, Butene, and Cyclopropane. \Remph{Journal of Chemical & Engineering Data} 2009;54:164855. \doi{10.1021/je800733b}.
 	
}
\examples{\dontrun{
# Second virial coefficient of methane at 300 K and standard pressure
CH4 <- Gas("methane")
print(CH4$B(T=300))
# Self-diffusion coefficient at 300 K
print(CH4$diffusion(T=300))
# create an instance of Gas for a molecule not listed in data_frame gas
Hg <- Gas("mercury")
# relative molecular mass
Hg$M <- 200.59
# mass of 1 molecule in kg
Hg$m <- Hg$M / pkg.env$Na / 1000
Hg$sigma <- 2.969 
Hg$epsk <- 750
print(Hg$thermal_conductivity(T=700))
}}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{name}:}{Chemical name}
    \item{\code{M}:}{relative molecular mass \eqn{M}}
    \item{\code{m}:}{Mass \eqn{m} of one gas particle in \code{kg}}
    \item{\code{sigma}:}{Lennard-Jones Parameter \Rsigma at \eqn{T = 0}}
    \item{\code{zeta}:}{Change of \Rsigma with temperature. For the noble gases \eqn{\zeta = 0}. Due to the vibrational excitation of molecules, \Rsigma is not independent of temperature but increasing with \eqn{T}.
    If \eqn{\sigma_0} is \Rsigma at \eqn{T = 0}, \Rsigma at \eqn{T} is given by \eqn{\sigma = \sigma_0 + \zeta T}.}
    \item{\code{epsk}:}{Well depth of the Lennard-Jones potential \eqn{\varepsilon/k} at \eqn{T = 0} in \code{K}.}
    \item{\code{dipole_moment}:}{dipole moment \eqn{\mu} in \code{Debye}.}
    \item{\code{polarizability}:}{polarizability \eqn{\alpha} in \code{Ao3}.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{B(T)}:}{Second virial coefficient \eqn{B(T)} in units \code{m3/mol}. The second virial coefficient provides systematic corrections
    to the ideal gas law. 
    	\if{latex}{
				\deqn{\frac{pV}{nRT} = 1 + \frac{n}{V} B + \dots}
    	}
	  	The second virial coefficient \eqn{B} depends only on the pair interaction between the particles (Vargas et. al. 2001). 
	  	\if{latex,html}{It is calculated from the modified Bessel function \eqn{I} (\link[Bessel]{BesselI}):
				\deqn{B(\Theta) = \frac{\sqrt{2}\pi^2\sigma^3}{3\Theta} \left[ 
						I_{-3/4}\left(\frac{1}{2\Theta}\right) + 
						I_{3/4}\left(\frac{1}{2\Theta}\right) - 
						I_{1/4}\left(\frac{1}{2\Theta}\right) - 
						I_{-1/4}\left(\frac{1}{2\Theta}\right)
					\right]
				}
	  	}
    }
    \item{\code{density(p=p0,T=T0)}:}{Gas density in units \code{kg/m3} incorporating the second virial correction.
			\if{latex,html}{
				\deqn{\rho = \frac{M}{2B(\Theta)} \left(\sqrt{1 + \frac{4pB(\Theta)}{RT}} - 1\right)}			
			}
    }
    \item{\code{diffusion(p=p0, T=T0, second_order_correction=TRUE)}:}{
    		Calculates the self diffusion coefficient in \code{m2/s}.
			\if{latex,html}{
				\deqn{D = \frac{3 \sqrt{\pi m kT}}{8\pi\sigma^2 \rho(p,T)\Omega^{(1,1)}(\Theta)} f_{D}^{(n)}}	
				\eqn{f_{D}^{(n)}} is the second order correction term. If \code{third_order_correction=FALSE} this term is set to 1.		
			}
    }
	\item{\code{thermal_conductivity(T, third_order_correction=TRUE)}:}{
    	Calculates the thermal conductivity for a monoatomic gas in \code{W/(m.K)}.
    	\if{latex,html}{ 
    		\deqn{\kappa = \frac{75k}{64 \sigma^2 \Omega^{(2,2)}(\Theta)}	\sqrt{\frac{kT}{\pi m}} f_{\kappa}^{(n)}}
    		\eqn{f_{\kappa}^{(n)}} is the third order correction term. If \code{third_order_correction=FALSE} this term is set to 1.
    		For a polyatomic gas the thermal conductivity can be calculated from the expression for a monoatomic gas using the classical Eucken correction (Bechtel \Remph{et al.} 2020)
    		\deqn{\kappa_p = \kappa \left(\frac{4 c_v}{15 R} + \frac{3}{5} \right)}
    		where \eqn{c_v} is the molar heat capacity of the polyatomic gas at constant volume.
    	}	
    }
    \item{\code{viscosity(T,third_order_correction=TRUE)}:}{
		\if{latex,html}{	
 			The viscosity \eqn{\eta} is given gas as
 			\deqn{\eta = \frac{5\sqrt{\pi m T}}{16 \pi \sigma^2 \Omega^{(2,2)}(\Theta)} f_{\eta}^{(n)}.}
 			\eqn{f_{\eta}^{(n)}} is the third order correction term. If \code{third_order_correction=FALSE} this term is set to 1.}
 		}
 		\item{\code{binary_diffusion(p=p0, T=T0, bathGas)}:}{
			Binary Diffusion coefficient of gas 1 in a bath gas 2 \code{bathGas}. The gas 1 may be polar or nonpolar. The bath gas 2 must be nonpolar.			
			\if{latex,html}{
 				The binary diffusion coefficient \eqn{D_{12}} is given as
 				\deqn{D_{12} = \frac{3}{16} \sqrt{\frac{2 \pi k T (m_1 + m_2)}{m_1 m_2}} 
				\left(\frac{kT}{\pi \sigma_{12}^2 \Omega^{(1,1)}(\Theta)p} \right)}
	 	  		The influence of the dipole moment of the gas is treated by the scaling parameter (Brown \Remph{et al.} 2011)
	 	  		\deqn{\xi = 1 + \frac{\bar{\alpha}\bar{\mu}^2}{16 \pi \epsilon_0} \sqrt{\frac{\varepsilon_1}{\varepsilon_2}},}
	 	  		where the reduced dipole moment is given by
	 	  		\deqn{\bar{\mu}^2=\frac{\mu^2}{\varepsilon_1 \sigma_1^3}}
	 	  		and the reduced polarizability of the bath gas is given by
	 	  		\deqn{\bar{\alpha} = \frac{\alpha}{\sigma_2^3}.}
	 	  		The scaling parameter \eqn{\xi} is used in the following combination rules to calculate the well depth of the effective interaction potential
	 	  		\deqn{\varepsilon_{12} = \xi^2 \sqrt{\varepsilon_1 \varepsilon_2}}
	 	  		and the collision diameter
	 	  		\deqn{\sigma_{12} = \xi^{-1/6} \frac{\sigma_1 + \sigma_2}{2}.}
	 	  		The second order correction (Marrero and Mason 1972) is not considered yet. Therefore, the diffusion coefficient of a polar gas in a non-polar bath gas is equal
	 	  		to the diffusion coefficient of a non-polar gas in a polar bath gas. For diffusion of a polar gas in a polar bath gas, this formula cannot be used.
	  		}
 		}
 		\item{\code{fit_B_data(B_df)}:}{
 		  Determination of Lennard-Jones parameters \Rsigma and \Rvarepsilon by nonlinear regression from a data frame of second virial coefficient data. The data frame must contain the columns 
 		  \code{T} for the temperature and \code{value} for the viscosity in units \code{cm3/mol}.
 		}
 		\item{\code{fit_viscosity_data(viscosity_df)}:}{
 		  Determination of Lennard-Jones parameters \Rsigma and \Rvarepsilon by nonlinear regression from a data frame of viscosity data. The data frame must contain the columns \code{T} for the 
 		  temperature and \code{value} for the viscosity in units \code{uPa.s}.
 		}	
 		\item{\code{fit_B_viscosity_data(B_df, viscosity_df, log=FALSE)}:}{
 		  Determination of Lennard-Jones parameters \Rsigma and \Rvarepsilon by simultaneous nonlinear regression from a data frame of viscosity data and a data frame of second virial 
 		  coefficient data. 
 		  \if{latex,html}{
 		    If \code{log=FALSE} the function to be minimized is given by (Bechtel \Remph{et. al.} 2020)
         \deqn{\chi^2(\sigma, \varepsilon) = \sum{\left(\frac{\eta(\sigma, \varepsilon, T) - \eta_{\mathrm{exp}}}{\Delta\eta_{\mathrm{exp}}}\right)^2}
           + \sum{\left(\frac{B(\sigma, \varepsilon, T) - B_{\mathrm{exp}}}{\Delta B_{\mathrm{exp}}} \right)^2}}
          where
          \deqn{\Delta\eta_{\mathrm{exp}} = | \mathrm{max}(\eta_{\mathrm{exp}}) - \mathrm{min}(\eta_{\mathrm{exp}})  |}
          and \deqn{\Delta B_{\mathrm{exp}} = | \mathrm{max}( B_{\mathrm{exp}}) - \mathrm{min}(B_{\mathrm{exp}}) |}
          If \code{log=TRUE} the function to be minimized is given by (Zarkova and Hohm 2009)
          \deqn{\chi^2(\sigma, \varepsilon) = \sum\left[\ln\left(\frac{\eta_{\mathrm{exp}}}{\eta(\sigma, \varepsilon, T)} \right) \right]^2 
           + \sum\left[\ln\left(\frac{|B_{\mathrm{exp}}|}{|B(\sigma, \varepsilon, T)|} \right) \right]^2 }	  
 		   }
 		}	
 	}
}
