pkg.env <- NULL

.onLoad <- function(libname, pkgname) {

	# defining global constants in package
	
	pkg.env <<- new.env(parent=emptyenv())	
	
	# Boltzmann constant
	assign('k', 1.380649E-23, envir=pkg.env)
	# Avogadro constant
	assign('Na', 6.02214076E23, envir=pkg.env)
	# gas constant
	assign("R", pkg.env$k* pkg.env$Na, envir=pkg.env)
	# Vacuum permittivity
	assign('eps0', 8.8541878128E-12, envir=pkg.env)
	# Debye unit
	assign('Debye',  3.33564E-30, envir=pkg.env)
	# Standard Pressure
	assign('p0', 101325, envir=pkg.env)
	# Standard Temperature
	assign('T0', 273.15, envir=pkg.env)
	
}

