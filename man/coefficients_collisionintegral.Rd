\name{coefficients_collisionintegral}
\alias{coefficients_collisionintegral}
\docType{data}
\title{
  Coefficients to calculate collision integrals for the Lennard Jones (12-6) potential
}
\description{
	Coefficients for \code{"\linkS4class{CollisionIntegral}"}, Table 1 from Kim and Monroe (2014). 
}
\usage{data("coefficients_collisionintegral")}
\format{
    \describe{
    \item{\code{l}}{order of collision integral \eqn{s}}
    \item{\code{s}}{order of collision integral \eqn{l}}
    \item{\code{A}}{model parameter \eqn{A}}
    \item{\code{B1}}{model parameter \eqn{B_1}}
    \item{\code{B2}}{model parameter \eqn{B_2}}
    \item{\code{B3}}{model parameter \eqn{B_3 \times 10}}
    \item{\code{B4}}{model parameter \eqn{B_4 \times 10}}
    \item{\code{B5}}{model parameter \eqn{B_5 \times 100}}
    \item{\code{B6}}{model parameter \eqn{B_6 \times 1000}}
    \item{\code{C1}}{model parameter \eqn{C_1}}
    \item{\code{C2}}{model parameter \eqn{C_2 \times 10}}
    \item{\code{C3}}{model parameter \eqn{C_3 \times 10}}
    \item{\code{C4}}{model parameter \eqn{C_4 \times 100}}
    \item{\code{C5}}{model parameter \eqn{C_5 \times 1000}}
    \item{\code{C6}}{model parameter \eqn{C_6 \times 10000}}
  }
}
\source{
Kim SU, Monroe CW. High-accuracy calculations of sixteen collision integrals for Lennard-Jones (12-6) gases and their interpolation to parameterize neon, argon, and krypton.
\Remph{Journal of Computational Physics} 2014;273:35873. \doi{10.1016/j.jcp.2014.05.018}.
}
\keyword{datasets}
