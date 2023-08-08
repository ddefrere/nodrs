;+
; NAME: FIZEAU_PSF
; 
; PURPOSE:
;   Create a synthetic Fizeau PSF image
;
; INPUTS:
;   diam      :  Time stamps corresponding to data (in sec)
;   lam       :  Data series
;   base      :  Corresponding error
;
; KEYWORDS
;   BANDWIDTH :  Set this keyword to the bandwidth of the observations (0 by default)
;   CEN_OBS   :  Set this keyword to the size of the central obscuration (0 by default)
;   PIX_SIZE  :  Set this keyword to the pixel size in mas. If not set, the code returns the automatic value.
;   VERBOSE   :
;
; MODIFICATION HISTORY:
;   Version 1.0, 16-APR-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 21-OCT-2014, DD: normalize output and minor bug corrected

FUNCTION FIZEAU_PSF, diam, lam, base, BANDWIDTH=bandwidth, CEN_OBS=cen_obs, PIX_SIZE=pix_size, VERBOSE=verbose

; Fixed paramaters
 m2r = 4.848136811D-9
 
; Running paramaters  
img_size = 2*1.22*2*lam/diam    ; image in rad (2 times the size of the airy pattern)
ratio    = cen_obs/diam         ; secondary to primary size ratio

; Keyword sanity check 
IF NOT KEYWORD_SET(CEN_OBS)   THEN cen_obs   = 0 ELSE IF cen_obs GT diam THEN MESSAGE, 'Central obscuration must be smaller than the aperture diameter'
IF NOT KEYWORD_SET(PIX_SIZE)  THEN pix_size  = img_size/128 ELSE pix_size *= m2r
IF NOT KEYWORD_SET(BANDWIDTH) THEN bandwidth = 0

; Derive the nimber of pixels (force an odd number of pixels)
n_pix = FIX(img_size/pix_size)
IF n_pix MOD 2 EQ 0 THEN n_pix -= 1

; Compute single-aperture PSF
r    = pix_size*SHIFT(DIST(n_pix,n_pix),FLOOR(0.5*n_pix),FLOOR(0.5*n_pix)) > 10E-10   ; 10E-10 to avoid dividing by 0
z    = !DPI*diam*r/lam
airy = 1./(1.-ratio^2)*(2.*BESELJ(z,1)/z-2.*ratio*BESELJ(ratio*z,1)/z)

; Compute fringe pattern
r    = (-FIX(0.5*n_pix) + DINDGEN(n_pix))*pix_size
tmp  = 1 + COS(2*!DPI*base*SIN(r)/lam)*SINC(!Dpi*base*SIN(r)*bandwidth/lam^2)
coh  = DBLARR(n_pix, n_pix, /NOZERO)
FOR i=0, n_pix-1 DO coh[0,i] = tmp

; Compute final image
fizeau_psf = (airy*coh)^2

RETURN, fizeau_psf/MAX(fizeau_psf)
END