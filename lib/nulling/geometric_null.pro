FUNCTION GEOMETRIC_NULL, th_star, bases, lambda, ULD=uld

; Fixed parameters
m2r = 4.848136811D-9     ; Converts milli arc seconds --> radians
geom_null = (!Dpi*bases*th_star*m2r/(4*lambda))^2*(1-7*uld/15)*(1-uld/3)

RETURN, geom_null
END