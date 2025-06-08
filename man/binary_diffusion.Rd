\name{binary_diffusion}
\alias{binary_diffusion}
\docType{data}
\title{
	Binary gas phase diffusion coefficients of methane, ethane, propane and butane in helium and nitrogen; fluoromethane, difluoromethane and trifluoromethane
	in nitrogen measured by reverse-flow gas chromatography 
}
\description{
	Reference data of binary diffusion coefficients for comparison with calculated diffusion coefficients. Diffusion coefficients were determined using a 
	reversed-flow gas chromatography system. 
}
\usage{data("binary_diffusion")}
\format{
  A data frame with 91 observations on the following 6 variables.
  \describe{
    \item{\code{doi}}{DOI of data source}
    \item{\code{bath_gas}}{Bath gas \code{helium} or \code{nitrogen}}
    \item{\code{gas}}{Diffusing species \code{C2H6}, \code{C3H8}, \code{C4H10}, \code{CH2F2}, \code{CH3F}, \code{CH4}, \code{CHF3}}
    \item{\code{T}}{Temperature in \code{K}}
    \item{\code{D}}{Diffusion coefficient in \code{cm2/s}}
    \item{\code{U_D}}{Uncertainty of diffusion coefficient in \code{cm2/s}}
  }
}
\details{
	\if{text}{See the HTML and PDF versions of this document for details.}
	\if{html}{\figure{binary_diffusion.svg}{options: width=1024 alt=""}}
	\if{latex}{\figure{binary_diffusion.pdf}{options: width=14cm}}
	\if{latex,html}{\cr\Remph{Plot of experimental diffusion coefficient vs. temperature. (a) diffusion of nitrogen and
	argon in helium. (b) diffusion of methane, ethane, propane and butane in helium. (c) diffusion of methane, ethane, propane and butane in nitrogen. (d) diffusion of fluoromethane,
	difluoromethane and trifluoromethane in nitrogen. The solid lines are calculated using the Lennard-Jones model. For figures (a)-(c) the Lennard-Jones parameters are taken from 
	\link{gas}, for figure (d) the Lennard-Jones parameters are estimated by nonlinear regression using \link[stats]{optim} from experimental data.}}
	
	The diffusion coefficient \eqn{D} as function of pressure in a narrow temperature range close to the reference temperature \eqn{T_0} is usually expressed as 
	(Langenberg \Remph{et al.} 2020)
	\deqn{D = D_0 \left(\frac{p_0}{p}\right)\left(\frac{T}{T_0}\right)^b}
	For the experimental data, the temperature coefficient \eqn{b} is obtained from the fit. For the calculated diffusion coefficients, the temperature coefficient is calculated by
	\deqn{b = \left(\frac{\partial D}{\partial T}\right)_{T_0} \left(\frac{T_0}{D_0}\right).}
	The diffusion coefficients \eqn{D_\mathrm{calc}} are calulated using \link{Gas-class}. The deviation is calculated by
	\deqn{\frac{D_\mathrm{exp} - D_\mathrm{calc}}{D_\mathrm{exp}}.}
	\tabular{llllllr}{Gas 		\tab 	Bath gas 	\tab  Experimental \tab       \tab Calculated  	\tab 		  	\tab	Deviation\cr
                              \tab              \tab 	\eqn{D_0} / \code{[cm2/s]}\tab \eqn{b} 	\tab 	\eqn{D_0} / \code{[cm2/s]}\tab 	\eqn{b} \cr
							nitrogen \tab 	helium 		\tab 	0.605(3) \tab 1.664(8)	\tab 	0.620		\tab 	1.68	\tab 	-3\%\cr
		               argon 	\tab 	helium 		\tab 	0.630(2) \tab 1.665(6)	\tab 	0.640		\tab 	1.68 	\tab 	-2\%\cr
		               methane	\tab 	helium 		\tab 	0.575(3) \tab 1.675(7) 	\tab	0.597 	\tab 	1.68	\tab 	-4\%\cr
		               ethane 	\tab 	helium 		\tab 	0.421(5) \tab 1.68(2) 	\tab 	0.446 	\tab 	1.70	\tab 	-6\%\cr
		               propane 	\tab 	helium 		\tab 	0.341(7)	\tab 1.67(2)	\tab 	0.361 	\tab 	1.70	\tab 	-7\%\cr
		               butane 	\tab 	helium 		\tab 	0.294(6)	\tab 1.65(2)	\tab 	0.368 	\tab 	1.74 	\tab -32\%\cr
		               methane 	\tab 	nitrogen 	\tab 	0.201(2)	\tab 1.74(2)	\tab 	0.186 	\tab 	1.83	\tab 	7\%\cr
		               ethane 	\tab 	nitrogen 	\tab 	0.136(2)	\tab 1.70(2)	\tab 	0.123 	\tab 	1.87 	\tab 	7\%\cr
		               propane 	\tab 	nitrogen 	\tab 	0.106(2)	\tab 1.72(3)	\tab	0.094 	\tab 	1.88 	\tab 	7\%\cr
		               butane 	\tab 	nitrogen 	\tab 	0.090(1)	\tab 1.72(2)	\tab 	0.084 	\tab 	1.97 	\tab -8\%\cr
	}
	With the exception of the diffusion of butane in helium, the calculated diffusion coefficients well resemble the measured diffusion coefficients within an error limit of < 10\%.
	For butane in helium more advanced combining rules need to be appplied (Li \Remph{et al.} 2023).
 
   The experimental data for the diffusion coefficients of difluoromethane and trifluoromethane can in turn be used to estimate the Lennard-Jones parameters for the Van der Waals
   interaction. The values for \eqn{\sigma} obtained are smaller than \eqn{\sigma} obtained from data of viscosity measurements (Shibasaki-Kitakawa \emph{et. al.} 1995,
   Clifford \emph{et al.} 1979).
   \if{latex}{\out{\pagebreak}} 	
	\tabular{lllllll}{Gas 				\tab  \eqn{D_0}/ \code{[cm2/s]} \tab \eqn{b} \tab	Viscosity \tab \tab Diffusion \cr
												\tab 			    \tab 			  \tab  \eqn{\sigma} / \code{[Ao]}\tab \eqn{\varepsilon/k} / \code{[K]} \tab \eqn{\sigma} / \code{[Ao]}\tab \eqn{\varepsilon/k} / \code{[K]} \cr
						fluoromethane 		\tab 0.1576(7)  \tab 1.784(8)	  \tab  --  \tab --   \tab 3.5 \tab 174 \cr
						difluoromethane 	\tab 0.133(2) 	 \tab 1.76(2)    \tab 4.9 \tab 204 \tab 3.9 \tab 153 \cr
						trifluoromethane 	\tab 0.123(2)   \tab 1.73(2)    \tab 4.4 \tab 182 \tab 4.5 \tab 63 \cr
	}
	This is due to the fact that both molecules have a dipole moment. This is why the intermolecular interaction of polar molecules cannot be described in terms of the 
	Lennard-Jones potential.
}
\references{
Clifford AA, Gray P, Scott AC. Viscosities of CFCl3, CF3Cl, CHFCl2, CHF2Cl and CHF3 from 373 to 570 K. \Remph{J. Chem. Soc., Faraday Trans. 1},
1979;75:1752. \doi{10.1039/F19797501752}
   
Langenberg S, Carstens T, Hupperich D, Schweighoefer S, Schurath U. Technical note: Determination of binary gas-phase diffusion coefficients of unstable
and adsorbing atmospheric trace gases at low temperature  arrested flow and twin tube method.
\Remph{Atmospheric Chemistry and Physics} 2020;20:366982. \doi{10.5194/acp-20-3669-2020}.
  
Li Y, Gui Y, You X. On the binary diffusion coefficients of n-alkanes in He/N2. \Remph{Combustion and Flame} 2023;257:112795. \doi{10.1016/j.combustflame.2023.112795}.

Shibasaki-Kitakawa N, Takahashi M, Yokoyama C, Takahashi S. Gas Viscosity of Difluoromethane from 298.15 to 423.15 K and up to 10 MPa 
\Remph{J. Chem. Eng. Data} 1995; 40:900-902 \doi{10.1021/je00020a036}
}
\source{
McGivern WS, Manion JA. Extending reversed-flow chromatographic methods for the measurement of diffusion coefficients to higher temperatures. 
\Remph{J. Chromatogr. A} 2011;1218:8432-42. \doi{10.1016/j.chroma.2011.09.035}.
	
McGivern WS, Manion JA. Hydrocarbon binary diffusion coefficient measurements for use in combustion modeling. \Remph{Combustion and Flame} 
2012;159:3021-6. \doi{10.1016/j.combustflame.2012.04.015}.

McGivern WS, Manion J. Binary Diffusion Coefficients for Methane and Fluoromethanes in Nitrogen. 
\Remph{Journal of Chemical & Engineering Data} 2021. \doi{10.1021/acs.jced.1c00161}.
}
\examples{
# binary diffusion data of nitrogen in bath gas helium
nitrogen_in_helium <- subset(binary_diffusion,(gas=="nitrogen" & bath_gas=="helium"))
print(nitrogen_in_helium)
}
\keyword{datasets}
