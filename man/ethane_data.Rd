\name{ethane_data}
\alias{ethane_data}
\docType{data}
\title{
Viscosity and second virial coefficient of ethane
}
\description{
Reference values for the second virial coefficient and viscosity of ethane from a intermolecular potential energy surface. The second virial coefficient was
calculated semiclassically by means of the Mayer-sampling Monte Carlo technique, while the transport properties were obtained using the classical kinetic theory of polyatomic gases.
The computed thermophysical property values are in excellent agreement with the best available experimental data and are recommended as reference values.  
\if{text}{See the HTML and PDF versions of this document for details.}
}
\usage{data("ethane_data")}
\format{
  A data frame with 107 observations on the following 4 variables.
  \describe{
    \item{\code{T}}{Temperature in \code{K}}
    \item{\code{property}}{Type of property. \code{B}: classically calculated second virial coefficient, \code{BQFH}: second virial coefficient calculated by a modification of the pair potential known as the quadratic Feynman-Hibbs (QFH)
    effective pair potential. \code{viscosity}: gas phase viscosity}
    \item{\code{unit}}{a factor with levels \code{cm3/mol}, \code{uPa.s}}
    \item{\code{value}}{numerical value of property. For the second virial coefficient the calculated data are supported by experimental data in the temperature range 220 - 623 K. For the 
    viscosity the calculated data are supported by experimental data in the temperature range 90 - 675 K.}
  }
}
\details{
	\if{html}{\figure{ethane.svg}{options: width=1024 alt=""}}
	\if{latex}{\figure{ethane.pdf}{options: width=14cm}}
	\if{latex,html}{\cr\Remph{(a) Viscosity of ethane. (b) second virial coefficient of ethane. The black solid curves are fits using viscosity data and second virial coefficient data
		respectively. The red dotted curves are simultaneous fits against the viscosity and second virial coefficient data.}
	}
	\tabular{lll}{
		Type of fit                                      \tab \Rsigma    \tab     \Rvarepsilon\eqn{/k} \cr
		vs. viscosity data                               \tab  4.38(2)   \tab     235(6) \cr
		vs. second virial coefficient data               \tab  4.95(6)   \tab     202(3) \cr    
		vs. viscosity and second virial coefficient data \tab  4.35      \tab     244 \cr	
	}
}
\source{
Hellmann R. Reference Values for the Second Virial Coefficient and Three Dilute Gas Transport Properties of Ethane from a State-of-the-Art Intermolecular Potential Energy Surface.
\Remph{Journal of Chemical & Engineering Data} 2018;63:470-81. \doi{10.1021/acs.jced.7b01069}.
}
\examples{
\dontrun{	
c2h6 <- Gas("ethane")
# estimate LJ-coefficients from viscosity data
ethane_viscosity <- subset(ethane_data, 
  (property=="viscosity") & (T>=90) & (T<=675), select=c(T, value)
)
c2h6$dsigma_dt <- 0
ethane_viscosity$value <- 1E-6*ethane_viscosity$value
c2h6$fit_viscosity_data(ethane_viscosity)
print(c2h6$sigma)
print(c2h6$epsk)

# estimate LJ-coefficients from second virial coefficient 
ethane_B <- subset(ethane_data, (property=="BQFH") & (T>= 220) & (T<=623))
ethane_B$value <- 1E-6*ethane_B$value
c2h6$fit_B_data(ethane_B)
print(c2h6$sigma)
print(c2h6$epsk)
}
}
\keyword{datasets}
