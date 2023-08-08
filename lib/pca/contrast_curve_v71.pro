;+
; NAME:
;   CONTRAST_CURVE_V71
;   
; PURPOSE:
;   Compute contrast curve from a single reduced image 
;   
; CALLING SEQUENCE:
;   CONTRAST_CURVE_V71, img, cx, cy, fwhm, plsc, dim, normal, sigma, BOX=box, GAUSSFILT=gaussfilt, $
;                       LEVEL_FILE=level_file, ALPHA=alpha, AGPM=agpm, STUDENT=student, DISPLAY=display
;                        
; DESCRIPTION:
;   Compute the noise level in the image as a function of the radial distance, using either the standard deviation on concentric annuli,
;   or a median of the local standard deviations computed in square boxes.
;   The noise level is then multiplied by the requested (gaussian) confidence level in terms of sigma to obtain the contrast curve.
;  
; INPUT VARIABLES:
;   img     : input image on which the contrast curve is to be computed 
;   cx      : x position of the star in the input image cube [pixels]
;   cy      : y position of the star in the input image cube [pixels]
;   fwhm    : FWHM of the telescope PSF in pixels
;   plsc    : plate scale [as/pixel]
;   dim     : size of the input image in number of pixels
;   sigma   : level of the injected companions relative to the noise  
;   
; KEYWORD PARAMETERS:
;   box        : if set to a number, compute noise using the azimuthal median of the standard deviations computed on small boxes with size = box*FWHM
;                if not set, a simple azimuthal RMS will be used to compute the noise in concentric annuli
;   gaussfilt  : if set, a gaussian filter will be applied to the image before estimating the noise 
;   level_file : if set, will look for already computed radial noise levels in predefined FITS file (procdir+'vec_'+tg_name_ori+'_level_r.fits')
;   alpha      : if set, the routine will look for a file (procdir+'vec_'+tg_name_bas+'_inter_alpha.fits') containing the self-attenuation as a function of pixel number
;   agpm       : if set, the fwhm of the PSF will be adapted to account for the use of the undersized Lyot stop
;   student    : if set, Student's t law will be used to compute the noise level at short angular separations 
;   display    : set this keyword to display intermediate and final results
;  
; COMMON BLOCK:
;   common_name : contains target name => implies that the data is available in datadir, 
;                 in the form 'img_'+tg_name+'_dc.fits'
;                 the parallactic angle and photometry vectors (optional) should also be available in datadir 
;                 in the form 'vec_'+tg_name+'_paral.fits' and 'vec_'+tg_name+'_photometry.fits'
; 
; OUTPUTS:
;   An output, a fits file is saved unedr the name procdir+'vec_'+tg_name_dyn+'_contrast.fits', with a dim_cc*3 table containing:
;     - the angular separations where the contrast has been computed
;     - the contrast curve (including self-attenuation, if the keyword alpha is set)
;     - the contrast curve before self-attenuation
;   Another fits file (procdir+'vec_'+tg_name_bas+'_level_r.fits') containing the radial noise level is saved if the keyword LEVEL_FILE is not set.
;   
; DEPENDENCIES:
;   READFITS, WRITEFITS, FILTER_IMAGE, ROBUST_SIGMA (AstroLib)
;   IMAGE_STDDEV : computes the local standard deviation in an image (using square boxes)

; MODIFICATION HISTORY:
;   Version 1.0, XX-YYY-ZZZZ, by Dimitri Mawet, ESO - JPL, originally developped for Keck pipeline
;   (...)
;   Version 7.0, 21-MAY-2013, DM: adapted for AGPM and implemented Student correction 
;   Version 7.1, 30-JUL-2013, OA: added header, cleaned up code, added keyword GAUSSFILT, now accepts image as input and saves radial noise level in a fits file
;-

PRO CONTRAST_CURVE_V71, img, cx, cy, fwhm, plsc, dim, sigma, BOX=box, GAUSSFILT=gaussfilt, $
                        LEVEL_FILE=level_file, ALPHA=alpha, AGPM=agpm, STUDENT=student, DISPLAY=display

COMMON common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

img_local = img
dim_cc = MIN([cx,cy,dim-cx,dim-cy]) ; define the minimum FOV radius in the image

; Load self-attenuation file if needed
IF KEYWORD_SET(alpha) THEN alpha_r = READFITS(procdir+'vec_'+tg_name_bas+'_inter_alpha.fits',header) ELSE alpha_r = 1d


; Compute the contrast curve
; --------------------------
IF KEYWORD_SET(level_file) THEN BEGIN ; If the contrast curve has already been computed, use it
  level_r = READFITS(procdir+'vec_'+tg_name_ori+'_level_r.fits',header)  ; load the radial noise level file
  contrast_img = sigma*TEMPORARY(level_r)/alpha_r  
ENDIF ELSE BEGIN ; If not, recompute the contrast curve
  level_r_th = FLTARR(360)
  level_r = FLTARR(dim_cc)
	; Remove NaNs from image, if any
	nans = WHERE(FINITE(img_local,/NAN) EQ 1)
	IF nans[0] NE -1 THEN img_local[nans]=0
	IF (SIZE(img_local))[1] NE dim THEN MESSAGE, 'Conflicting array size!!!'
	; Apply gaussian filter with width = fwhm (so that RMS is computed on resels hereafter)
	IF KEYWORD_SET(gaussfilt) THEN img_local = FILTER_IMAGE(img_local,fwhm_gaussian=fwhm)
	; Compute the RMS in boxes if needed
	IF KEYWORD_SET(box) THEN img_local = IMAGE_STDDEV(img_local,box/2D0*fwhm)
	; Loop on radii 
	FOR rad = 0d, dim_cc-1 DO BEGIN
	  ; Transform image in polar coordinates
		FOR th = 0d, 360-1 DO BEGIN
  		polar_coord = [th,rad]
  		rect_coord = CV_COORD(FROM_POLAR=polar_coord, /TO_RECT,/DEGREES) ; convert to polar coordinates to allow azimuthal treatment
  		level_r_th[th] = img_local[round(rect_coord[0])+cx,round(rect_coord[1])+cy]
    ENDFOR
  	r_fwhm = DOUBLE(rad)/fwhm
  	; Take into account the effect of the AGPM if needed
  	IF KEYWORD_SET(agpm) THEN BEGIN 
  		r_fwhm = r_fwhm*0.85 ; correct the effective fwhm
  		coro_thr = 1-EXP(-0.831*(r_fwhm)^2) ; transmission profile
  	ENDIF ELSE coro_thr = 1d
  	; Apply Student's correction if needed
  	IF KEYWORD_SET(student) THEN sigma_r = (1.0+(1.3291/(r_fwhm-0.64663))^1.1595) ELSE sigma_r=1d
  	; Compute the noise at the specified radius
  	IF KEYWORD_SET(box) THEN level_r[rad] = 1d/coro_thr*sigma_r*MEDIAN(level_r_th) $
  	ELSE level_r[rad] = 1d/coro_thr*sigma_r*ROBUST_SIGMA(level_r_th)
	ENDFOR
	WRITEFITS, procdir+'vec_'+tg_name_bas+'_level_r.fits', level_r ; save the radial noise level into a fits file if no noise level file was specified on input
	contrast_img = sigma*TEMPORARY(level_r)/alpha_r
ENDELSE


; Plot contrast curve
; -------------------
IF KEYWORD_SET(display) THEN BEGIN
  wset, 3
  plot, findgen(dim_cc)*plsc, contrast_img, /ylog, /xlog, yrange=[1e-6,1e-1],xrange=[1e-1,7],xtitle='Angular separation [arcsec]', ytitle='Contrast',$
        linestyle=0, XTicklen=1.0, YTicklen=1.0,XGridStyle=1, YGridStyle=1
  oplot, findgen(dim_cc)*plsc, contrast_img*alpha_r, linestyle=1
ENDIF

; Fixed parameter
m2r=   4.848136811D-9       ; Convert milli arc seconds to radians
r2m=   206264806.246D0      ;         radians --> mas
m2d=   2.7777777777778D-7   ;         mas --> degrees
r2d=   57.2957795128D0      ;         radians --> degrees
d2r=   1.74532925199D-2     ;         degrees --> radians
sr2as2=4.25451703D+10       ;         sr --> arcsec�
s2day= 1D0/8.64D+4          ;         sec --> day
d2a=   0.88622692545D0      ; =0.5*SQRT(pi): converts telescope diameter into area, after squaring
AU=    1.495978D+11         ; [m]
pc=    3.085678D+16         ; [m]

; Target parameters
dist  = 39.                 ; [pc]
sep_b = 68.*AU/(dist*pc)*r2m/1000.                 ; AU, Orbital separation (Marois 2010)
sep_c = 38.*AU/(dist*pc)*r2m/1000. 
sep_d = 24.*AU/(dist*pc)*r2m/1000. 
sep_e = 14.5*AU/(dist*pc)*r2m/1000. 
mag_b = 12.66-2.24          ; Marois et al. 2008
mag_c = 11.74-2.24
mag_d = 11.56-2.24
mag_e = 9.37                ; Marois et al. 2010

; Convert magnitude to contrast
flx_b = 10.^(-mag_b/2.5)
flx_c = 10.^(-mag_c/2.5)
flx_d = 10.^(-mag_d/2.5)
flx_e = 10.^(-mag_e/2.5)

; Telescope parameter
lambda  = 3.3*1D-6
dif_lim = ATAN(lambda/8.4)*dist*pc/au   ; diffraction limit

; Compute values of secondary axes
smin = 1e-1
smax = 7
cmax = 1e-6
cmin = 1e-2
mmin = -2.5*ALOG10(cmin)
mmax = -2.5*ALOG10(cmax)
amin = smin/AU*(dist*pc)/r2m*1000.
amax = smax/AU*(dist*pc)/r2m*1000.

mydevice = !D.NAME
SET_PLOT, 'PS'
IF keyword_set(alpha) THEN DEVICE, FILENAME=procdir+'plot_'+tg_name_dyn+'_contrast_alpha.eps', /ENCAPSULATE $
ELSE DEVICE, FILENAME=procdir+'plot_'+tg_name_dyn+'_contrast.eps', /ENCAPSULATE
plot, findgen(dim_cc)*plsc, contrast_img, /ylog, /xlog, yrange=[cmax,cmin], xrange=[smin,smax], xtitle='Angular separation [arcsec]', ytitle='5-!7r !1 Lp-band contrast',$
      linestyle=0, XTicklen=1.0, YTicklen=1.0,XGridStyle=1, YGridStyle=1, THICK=4, XTHICK=4, YTHICK=4, XSTYLE=9, YSTYLE=9, CHARTHICK=4
oplot, findgen(dim_cc)*plsc, contrast_img*alpha_r, linestyle=1, THICK=4

; Add planets
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, 1.4*COS(A), 1.4*SIN(A), /FILL, COLOR=220
t1 = 1 - 0.15 & t2 = 1 - 0.20
LOADCT, 13, /SILENT
OPLOT, [sep_b], [flx_b], PSYM=8 & XYOUTS, t1*[sep_b], t2*[flx_b], 'b', COLOR=235, CHARTHICK=3.0           ; b
OPLOT, [sep_c], [flx_c], PSYM=8 & XYOUTS, t1*[sep_c], t2*[flx_c], 'c', COLOR=235, CHARTHICK=3.0           ; c
OPLOT, [sep_d], [flx_d], PSYM=8 & XYOUTS, t1*[sep_d], t2*[flx_d], 'd', COLOR=235, CHARTHICK=3.0           ; d
OPLOT, [sep_e], [flx_e], PSYM=8 & XYOUTS, t1*[sep_e], t2*[flx_e], 'e', COLOR=235, CHARTHICK=3.0           ; e
LOADCT, 0, /SILENT

; Add axis
AXIS, YAXIS=1, YSTYLE=1, YLOG=0, YTHICK=4, CHARTHICK=4, YRANGE=[mmax,mmin], YTITLE='Magnitude contrast',  /SAVE
AXIS, XAXIS=1, XSTYLE=1, /XLOG, XTHICK=4, CHARTHICK=4, XRANGE=[amin,amax], XTITLE='Linear separation [AU]', /SAVE
DEVICE, /CLOSE
SET_PLOT, mydevice


; Write fits output
; -----------------
res_fin=fltarr(dim_cc,3)
res_fin[*,0]=findgen(dim_cc)*plsc
res_fin[*,1]=contrast_img
res_fin[*,2]=contrast_img*alpha_r
if keyword_set(alpha) then begin
  in = fix(0.15/plsc)
  out = fix(0.5/plsc)
  if mean(alpha_r[in:out]) ge 0.1 then begin
    print, '***************************************************************'
    print, '***************************************************************'
    print, 'Target name + reduction params: ', tg_name_dyn 
    print, 'Throughput-corrected X-sigma contrast between 0".15 and 0".5', median(res_fin[in:out,1])
    print, '***************************************************************'
    print, '***************************************************************'
  endif else begin
    print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print, '!!! Signal Attenuation > 10 !!!'
    print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  endelse
endif
WRITEFITS, procdir+'vec_'+tg_name_dyn+'_contrast.fits', res_fin

END