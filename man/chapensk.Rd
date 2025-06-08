\name{chapensk-package}
\alias{chapensk-package}
\alias{chapensk}
\alias{chapensk-init}
\alias{pkg.env}
\docType{package}
\title{
\packageTitle{chapensk}
}
\description{
\packageDescription{chapensk}
}
%\details{
%	\packageIndices{chapensk}
%}
\section{Introduction}{
Transport properties, such as viscosity, diffusion and thermal conductivity, play a crucial role in the modelling of combustion processes and chemical 
reactions. They depend on the intermolecular potential. In practice it is not necessary to have a detailed calculation of the intermolecular potential
for the calculation of transport properties.	
	
The interaction between spherical gas particles without a dipole moment can be described by the Lennard-Jones potential. 
\if{latex,html}{It is given by the following equation:
\deqn{U(r) = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]}{E(r) = 4/r)^12 - (6]}
where \eqn{r} is the distance between two interacting particles, \Rvarepsilon is the depth of the potential well and \Rsigma is the distance at which the 
particle-particle potential energy \eqn{U} is zero. The theory of van der Waals interaction gives the exponent 6 for the attractive term (London, 1937). 
The 12 exponent in the repulsive term is largely a matter of computational practicality, though it does represent the general nature of Pauli repulsion.} 

The Lennard-Jones potential only can be used for non-polar molecules, but sometimes is also used for polar molecules. However, for the
latter the Stockmayer (12-6-3) potential is more appropriate (Mourits and Rummens, 1977).

\if{latex,html}{
For monoatomic gases \Rsigma and \Rvarepsilon are independent of tmeperature. However, for non-monoatomic gases, averaging over different orientations 
and vibrational states gives temperature dependent parameters \Rsigma and \Rvarepsilon (Zarkova and Hohm, 2002).  The temperature dependency is simplified 
in this package as follows by a linear temperature coefficient \Rzeta:

\deqn{\sigma(T) = \sigma + \zeta T}
\deqn{\varepsilon(T) = \varepsilon\left(\frac{\sigma}{\sigma + \zeta T}\right)^6}
}

The Chapman-Enskog theory is a theoretical framework used to describe the transport properties of gases, such as viscosity, thermal conductivity, 
and diffusion coefficients. The theory is based on the idea that the properties of a gas can be related to the collisional interactions between 
individual gas molecules (Chapman and Larmor 1918). 

Collision integrals are mathematical expressions that arise in the Chapman-Enskog theory. They quantify the effects of molecular collisions on the transport properties of a gas.
}
\section{Methods}{
An object-oriented framework has been developed to calculate transport properties from potential parameters and vice versa. A class \link{Gas} has been defined for the 
calculation of the properties of gas. To facilitate calculations, a \link{gas} data set is provided for the properties of some common gases. The \link{CollisionIntegral}
class is used to calculate collision integrals using an interpolation function and fit parameters from datataset \link{coefficients_collisionintegral}.
}
\section{Results}{
Lennard-Jones parameters of non-polar molecules can be estimated using high quality of viscosity and second virial coefficients. This is demontrated for ethane, see data 
set \link{ethane_data}. For polar molecules Lennard-Jones parameters for the Van der Walls interaction part can be estimated from measurements of binary diffusion coefficients, 
see data set \link{binary_diffusion}.   
}
\section{Nomenclature}{
	\tabular{llll}{
	 Symbol					\tab Description                               												\tab Unit					\tab Global variable \cr				
	 \eqn{b}             \tab temperature coefficient of diffusion            										\tab \code{-}				\cr			
	 \eqn{B}             \tab second virial coefficient            													\tab \code{m3}\cr	
	 \eqn{D}					\tab \href{https://doi.org/10.1351/goldbook.D01719}{diffusion coefficient}			\tab \code{m2/2}\cr
	 \eqn{k}             \tab \href{https://doi.org/10.1351/goldbook.B00695}{Boltzmann constant}				\tab \code{J/K}         \tab \code{pkg.env$k}\cr
	 \eqn{m}					\tab molecular mass                       													\tab \code{kg}\cr
	 \eqn{M}             \tab \href{https://doi.org/10.1351/goldbook.R05271}{relative molecular mass}		\tab \code{-}\cr
	 \eqn{N_a}				\tab \href{https://doi.org/10.1351/goldbook.A00543}{Avogadro constant}				\tab \code{1/mol}       \tab \code{pkg.env$Na} \cr
	 \eqn{n}             \tab \href{https://doi.org/10.1351/goldbook.M03980}{mole}								\tab \code{mol}\cr	
	 \eqn{p}					\tab \href{https://doi.org/10.1351/goldbook.P04819}{pressure} 							\tab \code{Pa}\cr
	 \eqn{p_0}           \tab \href{https://doi.org/10.1351/goldbook.S05921}{standard pressure 101325 Pa} \tab \code{Pa}          \tab \code{pkg.env$p0} \cr
	 \eqn{p_c}				\tab \href{https://doi.org/10.1351/goldbook.C01397}{critical pressure} 				\tab \code{Pa}\cr
	 \eqn{R}             \tab \href{https://doi.org/10.1351/goldbook.G02579}{gas constant}  					\tab \code{J/(K.mol)}   \tab \code{pkh.env$R} \cr
    \eqn{T}      			\tab \href{https://doi.org/10.1351/goldbook.T06261}{temperature} 						\tab \code{K}\cr
    \eqn{T_0}      		\tab standard temperature 273.15 K            												\tab \code{K}           \tab \code{pkg.env$T0} \cr
    \eqn{T_c}				\tab \href{https://doi.org/10.1351/goldbook.C01402}{critical temperature}			\tab \code{K}\cr
    \eqn{V}					\tab gas volume																						\tab \code{m3}\cr
    \eqn{\alpha}        \tab \href{https://doi.org/10.1351/goldbook.P04711}{polarizability} 					\tab \code{Ao3}\cr
    \eqn{\bar{\alpha}}	\tab reduced polarizability																		\tab \code{-}\cr
    \eqn{\eta}          \tab \href{https://doi.org/10.1351/goldbook.D01877}{dynamic viscosity}				\tab \code{Pa.s}\cr
    \eqn{\epsilon}		\tab \href{https://doi.org/10.1351/goldbook.P04508}{permittivity of vacuum} 		\tab \code{F/m}         \tab \code{pkg.env$eps0} \cr
    \eqn{\varepsilon} 	\tab depth of potential well	  	         													\tab \code{J}\cr 
    \eqn{\kappa}        \tab \href{https://doi.org/10.1351/goldbook.T06299}{thermal conductivity}   		\tab \code{W/(m.K)}\cr
    \eqn{\mu}           \tab \href{https://doi.org/10.1351/goldbook.E01929}{dipole moment}               \tab \code{D}\cr
    \eqn{\bar{\mu}}		\tab reduced dipole moment  																		\tab \code{-}\cr
    \eqn{\Omega}			\tab reduced collision integral 	         													\tab \code{-}\cr
    \eqn{\rho}          \tab \href{https://doi.org/10.1351/goldbook.D01590}{gas density} 						\tab \code{kg/m3}\cr
    \eqn{\rho_c}        \tab \href{https://doi.org/10.1351/goldbook.C01396}{critical density}  				\tab \code{mol/l}\cr
    \eqn{\Theta} 			\tab reduced temperature        	         													\tab \code{-}\cr
    \eqn{\sigma}        \tab distance at which the potential energy is zero										\tab \code{Ao}\cr
    \eqn{\xi}				\tab scaling parameter      		         													\tab \code{-}\cr
    \eqn{\zeta}         \tab temperature coefficient of \eqn{\sigma}          									\tab \code{Ao/K}
  }
  Physical units are displayed in \href{https://ucum.org/}{UCUM} notation.
}

\author{
\packageAuthor{chapensk}

Maintainer: \packageMaintainer{chapensk}
}
\references{
Brown NJ, Bastien LAJ, Price PN. Transport properties for combustion modeling. 
\Remph{Progress in Energy and Combustion Science} 2011;37:565-82. \doi{10.1016/j.pecs.2010.12.001}.

Chapman S, Larmor J. V. On the kinetic theory of a gas. Part II. A composite monatomic gas: diffusion, viscosity, and thermal conduction.
\Remph{Philosophical Transactions of the Royal Society of London. Series A} 1918;217:11597. \doi{10.1098/rsta.1918.0005}.
   
London F. The general theory of molecular forces. \Remph{Transactions of the Faraday Society} 1937;33:8b. \doi{10.1039/tf937330008b}.
   
Mourits FM, Rummens FHA. A critical evaluation of Lennard-Jones and Stockmayer potential parameters and of some correlation methods. 
\Remph{Can. J. Chem.} 1977;55:300720. \doi{10.1139/v77-418}.

Zarkova L, Hohm U. pVT-Second Virial Coefficients \eqn{B(T)}, Viscosity \Reta and Self-Diffusion \eqn{\rho D(T)} of the Gases: BF3, CF4, 
SiF4, CCl4, SiCl4, SF6, MoF6, WF6, UF6, C(CH3)4, and Si(CH3)4 Determined by Means of an Isotropic Temperature-Dependent Potential. 
\Remph{Journal of Physical and Chemical Reference Data} 2002;31:183216. \doi{10.1063/1.1433462}.
}
\keyword{ package }

