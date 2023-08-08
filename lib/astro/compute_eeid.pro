;+
; NAME: COMPUTE_EEID
;
; PURPOSE:
;   Simple function to compute the Earth Equivalent Irradiation Distance (EEID).
;
; INPUTS
;   rad  : the stellar radius in mas
;   dist : the stellar distance in pc
;   teff : the effective temperature of the star in K
; 
; KEYWORDS:
;   LUM  : on output, the stellar luminosity
;   
; OUTPUT
;   The EEID in arcseconds
;
; MODIFICATION HISTORY:
;   Version 1.0, 02-DEC-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (adapted from former routine 'nomic_null.pro')

FUNCTION COMPUTE_EEID, rad, dist, teff, LUM=lum

; Fixed paramaters
m2r  = 4.848136811D-9
pc   = 3.0856776e+16        ; [m]  !CONST.PARSEC 
SB   = 5.6703D-8            ; [W/m2/K4]
Lsun = 3.9D+26              ; [W]

; Compute radius in m
rad_m = rad*m2r*dist*pc

; Compute stellar luminosity (relative to Lsun)
lum = 4*!Dpi*rad_m^2*SB*double(teff)^4/Lsun

; Compute the EEID
eeid_as = SQRT(lum)/dist

RETURN, eeid_as
END