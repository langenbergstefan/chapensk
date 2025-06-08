\name{CollisionIntegral-class}
\Rdversion{1.1}
\docType{class}
\alias{CollisionIntegral-class}
\alias{CollisionIntegral}

\title{Class \code{"CollisionIntegral"}}
\description{
	Class for the empirical interpolation functions for 16 reduced Lennard-Jones (12-6) collision integrals
	in the range 0.3 <= \RTheta <= 400.
	\if{html}{\figure{collisionintegral.svg}{options: width=640 alt=""}}
	\if{latex}{\figure{collisionintegral.pdf}{options: width=7cm}}
	\if{latex,html}{\cr\Remph{Collision integrals as function of reduced temperature \RTheta}
	\if{text}{See the HTML and PDF versions of this document for details.}
	}
}
\references{
Kim SU, Monroe CW. High-accuracy calculations of sixteen collision integrals for Lennard-Jones (12-6) gases and their interpolation to parameterize neon, argon, and krypton.
\Remph{Journal of Computational Physics} 2014;273:35873. \doi{10.1016/j.jcp.2014.05.018}.
}
\examples{
O11 <- CollisionIntegral(l=1,s=1)
O11$Omega(0.3)
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{l}:}{order of collision integral \eqn{l}}
    \item{\code{s}:}{order of collision integral \eqn{s}}
    \item{\code{coeff}:}{coefficients for collision integrals for order \eqn{l} and \eqn{s}, see \code{\link{coefficients_collisionintegral}}}
    \item{\code{A}:}{model parameter \Remph{A}}
    \item{\code{B}:}{model parameter \Remph{B\sub{1} ... B\sub{6}}}
    \item{\code{C}:}{model parameter \Remph{C\sub{1} ... C\sub{6}}}
  }
}
\section{Methods}{
  \describe{
    \item{\code{initialize(s, l)}:}{Initialize reduced collision integral of order \eqn{l} and \eqn{s}}
    \item{\code{Omega(Theta)}:}{Calculates reduced collision integral \ROmega(l,s) at reduced temperature \eqn{\Theta} using the empirical interpolation function
        \if{latex}{ 
    		\deqn{\Omega^{(l,s)} = A^{(l,s)} + \sum_{k=1}^6 \left[ \frac{B_k^{(l,s)}}{\Theta^k} 
    		+ C_k^{(l,s)}(\ln \Theta)^k\right]}}}
  }
}
