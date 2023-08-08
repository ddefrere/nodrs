; +
; NAME: SINC
;
; PURPOSE:
;   Returns SINC(X) defined as SIN(X*!PI)/(X*!PI)
;   Also handles the case for X=0
;
; INPUTS:
;   X          : The value to compute
;
; MODIFICATION HISTORY:
;   Version 1.0, 18-SEP-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu

FUNCTION SINC, X

ON_ERROR,2
eps=1D-10
y=!DPI*x
n=N_ELEMENTS(y)
IF n EQ 1 THEN BEGIN
  IF ABS(y) LE eps THEN f=1D0 ELSE BEGIN
    f=SIN(y)/y
    IF ABS(f) LT eps THEN f=0D0
  ENDELSE
ENDIF ELSE BEGIN
  f=DBLARR(n)+1D0
  w=WHERE(ABS(y) GT eps, n)
  IF n GT 0 THEN f(w)=SIN(y(w))/y(w)
  w=WHERE(ABS(f) LT eps, n)
  IF n GT 0 THEN f(w)=0D0
ENDELSE

RETURN, f
END
