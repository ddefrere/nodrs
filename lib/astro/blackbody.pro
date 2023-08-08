;+
; NAME: BLACKBODY
; 
; DESCRIPTION:
;   Returns the Planck blackbody surface brightness in Jy/sr,
;   or in W/m2/m/sr if the keyword STANDARD is set.
;   The Lambertian assumption is implicit in the omission of the factor pi!
;
; INPUTS:
;   T_eff:  effective temperature in K
;   lambda: wavelength in m
;
; KEYWORD:
;   standard: 0 (default) - result in Jy / sr
;             1           -           W / m2 / m / sr
;             2           -           ph / s / m2 / m / sr
;             3           -           ph / s / m2 / Hz / sr
;
; MODIFICATION HISTORY:
;   Version 1.0, 10-JAN-2002, by   Olivier Absil / Genie Team, Olivier.Absil@obsprm.fr
;   Version 1.1, 10-APR-2002, by   Roland den Hartog, ESA / SCI-ST / Genie Team, rdhartog@rssd.esa.int
;   Version 1.2, 11-MAR-2003, OA:  Directly compute the brightness in Jy/sr.
;   Version 1.3, 30-OCT-2003, RdH: Checked with Wien's displacement law
;   Version 1.4, 26-NOV-2003, RdH: Implemented standard application: output in W/m2/m/sr
;   Version 1.5, 17-SEP-2004, RdH: Implemented second standard: output in ph/s/m2/m/sr
;   Version 1.6, 19-FEV-2008, OA:  Implemented third standard: output in ph/s/m2/Hz/sr
;   Version 1.7, 02-DEC-2015, DD:  Cleanup comments

FUNCTION BLACKBODY, T_eff, lambda, STANDARD=standard
ON_ERROR, 2

h=6.6262D-34  ; J s
c=2.9979D+8   ; m/s
k=1.3807D-23  ; J/K
W2Jy=1D26     ; Convertion factor W --> Jy

; Surface brightness in Jy/sr because Jansky is 1E-26 W / m2 / Hz

IF NOT KEYWORD_SET(standard) THEN standard=0
CASE standard OF
  1:    RETURN, 2D0*h*c^2 / lambda^5 / (EXP(h*c/lambda/k/T_eff)-1D0)      ; W / m2 / m / sr
  2:    RETURN, 2D0*c / lambda^4 / (EXP(h*c/lambda/k/T_eff)-1D0)          ; ph / s / m2 / m / sr
  3:    RETURN, 2D0 / lambda^2 / (EXP(h*c/lambda/k/T_eff)-1D0)            ; ph / s / m2 / Hz / sr
  ELSE: RETURN, W2Jy * 2D0*h*c / lambda^3 / (EXP(h*c/lambda/k/T_eff)-1D0) ; Jy /sr
ENDCASE
END

;----------------------------------------------------------------
; Test harness

PRO TEST_BB

; Wien's displacement law is different for effectively using frequencies l_peak * T = 5.098E-3
dl=10D-6
l=1D-6 + dl*DINDGEN(101)/100.
s=BLACKBODY(1000., l, /STANDARD)
PLOTXY, l, s*!DPI*1D-26*1D-6, /HOLD, /YLOG;, YRANGE=[2., 20000.]
s=BLACKBODY(800., l, /STANDARD)
PLOTXY, l, s*!DPI*1D-26*1D-6, /ADD, /YLOG, SYMBOL=-5
s=BLACKBODY(289.78, l, /STANDARD)
PLOTXY, l, s*!DPI*1D-26*1D-6, /ADD, /YLOG, SYMBOL=-4

s=BLACKBODY(1000., l)
PLOTXY, l, s*!DPI*1D-26*1D-6, /ADD, /YLOG, SYMBOL=-3
s=BLACKBODY(800., l)
PLOTXY, l, s*!DPI*1D-26*1D-6, /ADD, /YLOG, SYMBOL=-2
s=BLACKBODY(289.78, l)
PLOTXY, l, s*!DPI*1D-26*1D-6, /ADD, /YLOG, SYMBOL=-1

; Stefan-Boltzmann's law
sigma=5.6703D-8 ; W m^-2 T^-4
T=100.
bw=100.D-6
l=1D-6 + bw*DINDGEN(1001)/1000.
s=BLACKBODY(T, l, /STANDARD)
PLOTXY, l, s
PRINT, INT_TABULATED(l, s), sigma*T^4

END