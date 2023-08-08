;+
; NAME: GET_PRM
; 
; PURPOSE:
;   Initializes global astronomical and physical constants
;
; OUTPUTS:
;   prm: structure containing all relevant parameters.
;
; MODIFICATION HISTORY:
;   Version 1.0,  30-AUG-2005, by Roland den Hartog, ESA / ESTEC / Darwin team, rdhartog@rssd.esa.int
;   Version 1.1,  03-MAR-2015, DD: updated to use IDL CONST structure
;   Version 1.1,  21-DEC-2015, DD: removed CONST which doesn't exist before IDL 8.2

PRO GET_PRM, prm

; Parameters that are used everywhere
h=     6.6262D-34           ; [J s]
c=     299792458D0          ; [m/s]
k=     1.3807D-23           ; [J/K]
G=     6.6725985d-11        ; [m�/kg/s�] Gravitational const.
SB=    5.6703D-8            ; [W/m�/K^4]
w2jy=  1D26                 ;            1 Jy is 1E-26 W / m� / Hz
pi=    !Dpi
twopi= 2*pi
m2r=   4.848136811D-9       ; Convert milli arc seconds to radians
r2m=   206264806.246D0      ;         radians --> mas
m2d=   2.7777777777778D-7   ;         mas --> degrees
r2d=   57.2957795128D0      ;         radians --> degrees
d2r=   1.74532925199D-2     ;         degrees --> radians
sr2as2=4.25451703D+10       ;         sr --> arcsec�
s2day= 1D0/8.64D+4          ;         sec --> day
d2a=   0.88622692545D0      ; =0.5*SQRT(pi): converts telescope diameter into area, after squaring
AU=    1.4959787e+11        ; [m]
pc=    3.0856776e+16        ; [m]
; Solar parameters
msun=  1.9884159e+30        ; [kg]
tsun=  5777.                ; [K]
rsun=  6.96D8               ; [m]     solar radius
lsun=  3.9D+26              ; [W]
; Earth parameters
day =  86400D0              ; [s]
year=  31558118.4D0         ; [s]
mear=  5.9721864e+24        ; [kg]
tear=  265.                 ; [K]
rear=  6378136.6            ; [m]
lear=  1.                   ; [W]
eear=  0.017                ;         eccentricity

prm={H:h, C:c, K:k, G:g, SB:sb, W2JY:w2jy, PI:pi, TWOPI:twopi, M2R:m2r, R2M:r2m, M2D:m2d, R2D:r2d, D2R:d2r, SR2AS2:sr2as2, S2DAY:s2day, D2A:d2a, AU:au, PC:pc, $
     MSUN:msun, TSUN:tsun, RSUN:rsun, LSUN:lsun, DAY:day, YEAR:year, MEAR:mear, TEAR:tear, REAR:rear, LEAR:lear, EEAR:eear}

RETURN
END
