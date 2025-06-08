\name{gas}
\alias{gas}
\docType{data}
\title{
Gas data
}
\description{
	Physical properties of some gas species. 
}
\usage{data("gas")}
\format{
  A data frame with properties of 43 gases on the following 12 variables.
  \describe{
    \item{\code{formula}}{chemical formula}
    \item{\code{name}}{chemical name}
    \item{\code{CAS}}{Chemical abstracts registry number}
    \item{\code{group}}{\href{https://doi.org/10.1351/goldbook.P04703}{point group}}
    \item{\code{M}}{\href{https://doi.org/10.1351/goldbook.R05271}{relative molecular mass}}
    \item{\code{dipole_moment}}{\href{https://doi.org/10.1351/goldbook.E01929}{electric dipole moment} in \code{Debye}}
    \item{\code{polarizability}}{\href{https://doi.org/10.1351/goldbook.E01933}{electric polarizability} in \code{Ao}}
    \item{\code{IE}}{\href{https://doi.org/10.1351/goldbook.I03199}{ionization energy} in \code{eV}}
    \item{\code{Tc}}{\href{https://doi.org/10.1351/goldbook.C01402}{critical temperature} in \code{K}}
    \item{\code{pc}}{\href{https://doi.org/10.1351/goldbook.C01397}{critical pressure} in \code{bar}}
    \item{\code{rhoc}}{critical density in \code{mol/l}}
    \item{\code{sigma}}{distance at which the intermolecular potential between the two particles is zero in \code{Ao}}
    \item{\code{epsk}}{Well depth \eqn{\varepsilon/k} of the Lennard-Jones potential in \code{K}}
    \item{\code{DOI}}{Data source of Lennard-Jones parameters. If not specified otherwise, they are taken from Poling \Remph{et al.} (2004)}
  }
}
\details{
	Properties of simple gases for the calculation of transport properties. Not all properties are given for all molecules.
	\if{text}{See the HTML and PDF versions of this document for details.}
	\if{latex,html}{ 
	However, correlations exists between the Lennard-Jones parameters and critical data as derived from numerous numerical simulations 
	of the Lennard-Jones fluid (Stephan \Remph{et al.} 2019): \Rvarepsilon and \Rsigma can be determined from critical temperature 
	\deqn{T_c = (1.321 \pm 0.007) (\varepsilon/k)} 
	critical density
	\deqn{\rho_c = (0.316 \pm 0.005) / \sigma^3} 
	and critical pressure 
	\deqn{p_c = (0.129 \pm 0.005) (\varepsilon/\sigma^3)}
	where \eqn{k} is the Boltzmann constant.}
	\if{html}{\figure{critical_data.svg}{options: width=1024 alt=""}}
	\if{latex}{\figure{critical_data.pdf}{options: width=14cm}}
	\if{latex,html}{\cr\Remph{Correlation of Lennard-Jones parameters of nonpolar gases with critical data: The solid lines are the 
	expectations from numerical simulations of the Lennard-Jones Fluid. 
	(a) critical temperature vs. well depth \eqn{\varepsilon}. (b) critical density vs. \eqn{1/\sigma^3.}}}
}
\references{
  NIST Chemistry WebBook. NIST Standard Reference Database, vol. 69, 2023. \doi{10.18434/T4D303}.
  
  NIST. Experimental Polarizabilites. in: III RDJ, editor. NIST Computational Chemistry Comparison and Benchmark Database, NIST; 2020. \doi{10.18434/T47C7Z}.
     
  Poling BE, Prausnitz JM, OConnell JP. The Properties of Gases and Liquids. 5 ed. New York: McGraw-Hill; 2004. 
  
  Stephan S, Thol M, Vrabec J, Hasse H. Thermophysical Properties of the Lennard-Jones Fluid: Database and Data Assessment.
  \Remph{Journal of Chemical Information and Modeling} 2019;59:424865. \doi{10.1021/acs.jcim.9b00620}.
}
\keyword{datasets}
