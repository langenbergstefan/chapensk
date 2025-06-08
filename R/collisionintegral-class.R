CollisionIntegral <- setRefClass("CollisionIntegral",
	fields = list(
		l="integer",
		s="integer",
		coeff="data.frame",
		A="numeric",
		B="numeric",
		C="numeric"
	),
	methods = list(
		initialize = function(l,s) {
			l1 <- l
		   s1 <- s
			coeff <<- subset(coefficients_collisionintegral,(l==l1) & (s==s1))
			stopifnot(nrow(coeff) == 1)
			A <<- coeff$A
			B <<- c(coeff$B1, coeff$B2, 0.1*coeff$B3, 0.1*coeff$B4, 0.01*coeff$B5, 0.001*coeff$B6)
			C <<- c(coeff$C1, 0.1*coeff$C2, 0.1*coeff$C3, 0.01*coeff$C4, 0.001*coeff$C5, 0.0001*coeff$C6)
		},
		Omega = function(Theta) {
		  stopifnot(Theta >= 0.3, Theta <= 400)		
			sum <- A
			for(k in 1:6) {
				sum <- sum + B[k] / Theta^k + C[k] * log(Theta)^k		
			}
			return <- sum
		}
	)
)	