;+
; NAME:
;   ALPHA_R_V61
;   
; PURPOSE:
;   Evaluates the self-attenuation in (ADI) image processing  
;   
; CALLING SEQUENCE:
;   ALPHA_R_V61, fwhm, imgfin_nofcp, n_br, rad_c, dim, DISPLAY=display
;                        
; DESCRIPTION:
;   Compares the results of data processing with an without fake companions added to the image cube,
;   to derive the amount of self-subtraction encountered in the process.
;  
; INPUT VARIABLES:
;   fwhm    : FWHM of the telescope PSF in pixels
;   imgfin_nofcp : final median image after processing for the cube without fake companions  
;   n_br    : number of radial branches for the injected companions
;   rad_c   : radial positions of the fake companions [arcsec]
;   dim     : size of the input image in number of pixels
;   
; KEYWORD PARAMETERS:
;   display : set this keyword to display intermediate and final results
;  
; COMMON BLOCK:
;   common_name : contains target name => implies that the data is available in datadir, 
;                 in the form 'img_'+tg_name+'_dc.fits'
;                 the parallactic angle and photometry vectors (optional) should also be available in datadir 
;                 in the form 'vec_'+tg_name+'_paral.fits' and 'vec_'+tg_name+'_photometry.fits'
; 
; OUTPUTS:
;   Three fits files containing
;      - the 2D map of the attenuation, under procdir+'img_'+tg_name_bas+'_alpha_r_2D.fits'
;      - the radial attenuation at the location of the fake companions, under procdir+'vec_'+tg_name_bas+'_alpha_r.fits'
;      - an interpolation of the radial attenuation on half the image size, under procdir+'vec_'+tg_name_bas+'_inter_alpha.fits'
;   
; DEPENDENCIES:
;   READFITS, WRITEFITS (AstroLib)
;   SHIFTI : shift a 1-d or 2-d array by fractional pixels using interpolation
;   The routine expects to have input PSF files named "fake_companion_INST_FILT_(PLSC).fits" in the routine directory
;     (with INST the instrument name, FILT the filter name and PLSC the plate scale in case various settings are possible for the considered instrument)
;
; LIMITATIONS:
;   - Assumes that the fake companions are placed on radial branches, with the first branch being horizontal and oriented towards the x positive direction.
;   - The 2D map created by the routine is filled with zeros beyond a radius equal to half the image size
; 
; MODIFICATION HISTORY:
;   Version 1.0, XX-YYY-ZZZZ, by Dimitri Mawet, ESO - JPL, originally developped for Keck pipeline
;   (...)
;   Version 6.0, 21-MAY-2013, DM: adapted for NACO/AGPM data 
;   Version 6.1, 24-JUL-2013, OA: cleaned up code, added header
;-

PRO ALPHA_R_V61, fwhm, imgfin_nofcp, n_br, rad_c, dim, DISPLAY=display

COMMON common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir


; Read the input images and arrays
; --------------------------------
imgfin_fcp = READFITS(procdir+'img_'+tg_name_dyn+'_median.fits')
fcp_array = READFITS(procdir+'img_'+tg_name_bas+'_fcparray.fits') ;/1000d ; level of injected fake companions
IF (SIZE(fcp_array))[1] NE dim THEN MESSAGE, 'Conflicting array size!!!' 

fcp_coords = READFITS(procdir+'vec_'+tg_name_bas+'_coords.fits')
n_rad = (SIZE(fcp_coords))[3]


; Derive the local attenuation
; ----------------------------
alpha_fcp = FLTARR(n_br,n_rad) ; Mean value of alpha on the whole fwhm of each fake companion
alpha = FLTARR(n_rad) ; Mean azimutal value of alpha for the different radii
alpha_2d = FLTARR(dim,dim)

ratio = MEDIAN((imgfin_fcp - imgfin_nofcp) / fcp_array, fwhm)
;ratio = (imgfin_fcp - imgfin_nofcp) / fcp_array
IF KEYWORD_SET(display) THEN BEGIN
  WSET, 3
  TVSCL, congrid(ratio,500,500)
ENDIF
WRITEFITS, procdir+'ratio_test.fits', ratio


; Perform "aperture photometry"
; -----------------------------
anti_aliasing_f = 3
dim_aa = anti_aliasing_f*dim
w_t = SHIFT(DIST(dim_aa), dim_aa/2, dim_aa/2)
w_t = TEMPORARY(w_t) LE CEIL(anti_aliasing_f*1.0*FLOAT(fwhm)/2) ; define an aperture of diameter = FWHM
w_t = SHIFT(REBIN(DOUBLE(TEMPORARY(w_t)), dim, dim), -dim/2, -dim/2) ; aperture centred on pixel (0,0)
WRITEFITS, 'aperture.fits', w_t

FOR i=0, n_br-1 DO BEGIN
	FOR j=0, n_rad-1 DO BEGIN
	  xc = (fcp_coords[*,i,j])[0];-dim/2
    yc = (fcp_coords[*,i,j])[1];-dim/2
    aper_tmp = SHIFTI(w_t, xc, yc) ; shift the aperture to the companion position
    alpha_fcp[i,j] = TOTAL((ratio*aper_tmp)[WHERE(aper_tmp GT 0)]) / TOTAL(aper_tmp)	; aperture "photometry" for the ratio
	ENDFOR
ENDFOR


; Azimutal mean of alpha --> alpha(r)
; ------------------------------------
alpha = TOTAL(alpha_fcp,1)/n_br
IF KEYWORD_SET(display) THEN BEGIN
  WSET, 0
  PLOT, alpha
ENDIF


; Interpolation of alpha(r)
; -------------------------
x_interp_out = FINDGEN(dim/2)
x_interp = REFORM(fcp_coords[0,0,*]) - dim/2 ; x coordinates of the fake companions along (first) the horizontal branch in the pattern, relative to the image center
interp_alpha = INTERPOL(alpha, x_interp, x_interp_out) ;, /spline)
IF KEYWORD_SET(display) THEN BEGIN
  WSET, 1
  PLOT, interp_alpha
ENDIF


; Creation of the 2D representation
; ---------------------------------
FOR x = 0d, dim-1 DO BEGIN
  FOR y = 0d, dim-1 DO BEGIN
    rect_coord = [x-dim/2, y-dim/2]
		polar_coord = CV_COORD(FROM_RECT=rect_coord, /TO_POLAR, /DEGREES) 
		IF ABS(polar_coord[1]) GT dim/2-1 THEN alpha_2d[x,y] = 0 ELSE alpha_2d[x,y] = interp_alpha[ABS(polar_coord[1])]
	ENDFOR
ENDFOR


; Final outputs in fits files
; ---------------------------
WRITEFITS, procdir+'img_'+tg_name_bas+'_alpha_r_2D.fits', alpha_2d
WRITEFITS, procdir+'vec_'+tg_name_bas+'_alpha_r.fits', alpha
WRITEFITS, procdir+'vec_'+tg_name_bas+'_inter_alpha.fits', interp_alpha

END