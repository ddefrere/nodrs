;+
; NAME:
;   FCP_V83
;   
; PURPOSE:
;   Fake companion injection into an image cube
;   
; CALLING SEQUENCE:
;   FCP_V83, fwhm, sigma, rad_c, n_br, plsc, dim, camera_filter, right_handed, ADI=adi, KECK=keck, LBT=lbt, NACO=naco, NICI=nici, $
;            EVAL_NOISE=eval_noise, NORMAL=normal, GAUSSFILT=gaussfilt, BOX=box, AGPM=agpm, STUDENT=student, DISPLAY=display
;                        
; DESCRIPTION:
;   Injection of fake companion into an image cube, to assess the attenuation of an image processing method.
;   Companions are injected at various radial distances and on several branches (spread on 360�).
;  
; INPUT VARIABLES:
;   fwhm    : FWHM of the telescope PSF in pixels
;   sigma   : level of the injected companions relative to the noise  
;   rad_c   : radial positions of the fake companions [arcsec]
;   n_br    : number of radial branches for the injected companions
;   plsc    : plate scale [as/pixel]
;   dim     : size of the input image in number of pixels
;   camera_filter : name of the filter used in the camera
;   right_handed  : for right_handed WCS (rare, e.g. NICI) 
;   
; KEYWORD PARAMETERS:
;   adi        : set this keyword if the input cube is taken in ADI mode 
;                --> will look for a file '..._paral.fits' where the parallactic angle is stored for each image of the cube
;   keck       : set this keyword to display intermediate and final results
;   naco       : display various outputs/status in the console
;   nici       : set this keyword for WCS data
;   eval_noise : if set, will (re)evaluate the radial noise level and save it into a FITS file (procdir+'vec_'+tg_name_ori+'_level_r.fits')
;   normal     : photometric normalisation factor. If not set, looks for 'vec_'+tg_name+'_photometry.fits'.
;   gaussfilt  : if set, a gaussian filter will be applied to the image before estimating the noise 
;   box        : if set to a number, compute noise using azimuthal median of RMS computed on small boxes with size = box*FWHM
;                if not set, a simple azimuthal RMS will be used to compute the noise in annuli
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
;   A new fits image cube containing the fake companions is saved under procdir+'img_'+tg_name_dyn+'_dc.fits'.
;   In addition, three other fits outputs: 
;      - the noise level as a function of radius under procdir+'vec_'+tg_name_bas+'_level_r.fits'
;      - the image of the fake companion "pattern" under procdir+'img_'+tg_name_bas+'_fcparray.fits'
;      - the coordinates of the fake companions under procdir+'vec_'+tg_name_bas+'_coords.fits'
;   
; DEPENDENCIES:
;   READFITS, WRITEFITS, SXPAR, FILTER_IMAGE (AstroLib)
;   IMAGE_STDDEV : computes the local standard deviation in an image (using square boxes)
;   SHIFTI : shift a 1-d or 2-d array by fractional pixels using interpolation
;   CONTRAST_CURVE_V71 : computation of the radial noise level (in case EVAL_NOISE is set)
;   The routine expects to have input PSF files named "fake_companion_INST_FILT_(PLSC).fits" in the routine directory
;      (with INST the instrument name, FILT the filter name and PLSC the plate scale in case various settings are  
;      possible for the considered instrument)
;      
; LIMITATIONS:
;   This routine assumes the star to be located at the center of the image.
;   The first input branch is along the x axis, and all other branches are evenly spread in angles. 
;
; MODIFICATION HISTORY:
;   Version 1.0, XX-YYY-ZZZZ, by Dimitri Mawet, ESO - JPL, originally developped for Keck pipeline
;   (...)
;   Version 8.0, 21-MAY-2013, DM: adapted for AGPM and implemented Student correction 
;   Version 8.1, 24-MAY-2013, OA: added header, cleaned up code 
;   Version 8.2, 23-JUL-2013, OA: added keyword gaussfilt for image filtering, removed rotation of companion PSF in ADI cube
;   Version 8.3, 30-JUL-2013, OA: transformed into a procedure, removed radial noise computation (deferred to CONTRAST_CURVE), removed input image (useless)
;   Version 8.4, 27-NOV-2013, DD: updated for LBT
;-


PRO FCP_V83, fwhm, sigma, rad_c, n_br, plsc, dim, camera_filter, right_handed, ADI=adi, KECK=keck, LBT=lbt, NACO=naco, NICI=nici, $
             EVAL_NOISE=eval_noise, NORMAL=normal, GAUSSFILT=gaussfilt, BOX=box, AGPM=agpm, STUDENT=student, DISPLAY=display

COMMON common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir


;Read input data
;---------------
obj_tmp=readfits(procdir+'img_'+tg_name_dyn+'_dc.fits',header)
nref=sxpar(header,'NREF') 

nobj=(size(obj_tmp))[3]; Number of frames in the input object datacube
if (size(obj_tmp))[1] ne dim then MESSAGE, 'Conflicting array size!!!'

IF KEYWORD_SET(normal) THEN normal_vec = fltarr(nobj) + normal $ 
ELSE normal_vec = readfits(procdir+'vec_'+tg_name_bas+'_photometry.fits')


;Read the appropriate input PSF
;------------------------------
r = ROUTINE_INFO('fcp_v83', /SOURCE)
path = STRMID(r.path,0,STRLEN(r.path)-STRLEN('fcp_v81.pro'))

IF KEYWORD_SET(keck)+KEYWORD_SET(nici)+KEYWORD_SET(naco) GT 1 THEN MESSAGE, 'Please choose only one instrument among Kec, NICI, NACO.'

IF KEYWORD_SET(keck) THEN BEGIN
  ; Keck plate scale is unique and = 10 mas/pix
  if camera_filter eq 'Jj' then fcp=readfits(path+'fake_companion_keck_J.fits')
  if camera_filter eq 'Hh' then fcp=readfits(path+'fake_companion_keck_H.fits')
  if camera_filter eq 'Kk' then fcp=readfits(path+'fake_companion_keck_K.fits')
  if camera_filter eq 'Ll' then fcp=readfits(path+'fake_companion_keck_L.fits')
ENDIF ELSE IF KEYWORD_SET(nici) THEN BEGIN
  if camera_filter eq 'Hh' then fcp=readfits(path+'fake_companion_nici_H_s18.fits')
ENDIF ELSE IF KEYWORD_SET(naco) THEN BEGIN
  ; NACO plate scale is either 13 or 27 mas/pix
  if plsc eq 0.013 then begin
    if camera_filter eq 'Jj' then fcp=readfits(path+'fake_companion_naco_J_s13.fits')
    if camera_filter eq 'Hh' then fcp=readfits(path+'fake_companion_naco_H_s13.fits')
    if camera_filter eq 'Kk' then fcp=readfits(path+'fake_companion_naco_K_s13.fits')
  endif
  if plsc eq 0.0271 then begin
    if camera_filter eq 'Jj' then fcp=readfits(path+'fake_companion_naco_J_s27.fits')
    if camera_filter eq 'Hh' then fcp=readfits(path+'fake_companion_naco_H_s27.fits')
    if camera_filter eq 'Kk' then fcp=readfits(path+'fake_companion_naco_K_s27.fits')
    if camera_filter eq 'Ll' then fcp=readfits(path+'fake_companion_naco_L_l27.fits')
    if camera_filter eq 'apo' then fcp=readfits(path+'fake_companion_naco_Lp_apo_l27.fits')
  endif
ENDIF ELSE IF KEYWORD_SET(lbt) THEN BEGIN
  ; 2x2 binned LBT/LMIRcam plate scale
  if plsc eq 0.0213600 then begin
    if camera_filter eq 'Lp' then fcp=readfits(path+'fake_companion_lbt_Lp_s21.fits')
  endif
ENDIF


; Compute the radial noise levels if needed, and load it
; ------------------------------------------------------
IF KEYWORD_SET(eval_noise) THEN BEGIN
  img = READFITS(procdir+'img_'+tg_name_dyn+'_pca_median.fits')
  CONTRAST_CURVE_V71, img, dim/2, dim/2, fwhm, plsc, dim, sigma, /GAUSSFILT, /AGPM, BOX=box, DISPLAY=display
ENDIF
level_r = READFITS(procdir+'vec_'+tg_name_ori+'_level_r.fits')  ; load the radial noise level file
  
IF KEYWORD_SET(display) THEN BEGIN
	wset, 2
	plot, findgen(dim/2)*plsc, sigma*level_r, /ylog, /xlog, yrange=[1e-8,1e-2],xrange=[1e-1,7],xtitle='Angular separation [arcsec]', ytitle='FCP level',linestyle=0, XTicklen=1.0, YTicklen=1.0,XGridStyle=1, YGridStyle=1
ENDIF


;Compute coordinates and flux level of fake companions
;-----------------------------------------------------
arr_final = fltarr(dim,dim)
fcp_arr = fltarr(dim,dim)
sz_fcp = (size(fcp))[1]/2
fcp_arr[dim/2-sz_fcp,dim/2-sz_fcp] = fcp

n_rad = (size(rad_c))[1] ; Number of different radii where we want to place our fake companions
coords = intarr(2,n_br,n_rad) ; Coordinates of the center of all the fake companions

FOR j = 0, n_br-1 DO BEGIN
  theta = j*2*!DPI/n_br
  FOR i = 0, n_rad-1 DO BEGIN
    r_pix = DOUBLE(rad_c[i])/plsc
    xpos = r_pix * COS(theta)
    ypos = r_pix * SIN(theta)
    coords[0,j,i] = dim/2 + xpos
    coords[1,j,i] = dim/2 + ypos
    arr_final = arr_final + SHIFTI(fcp_arr, xpos, ypos) * sigma * level_r[FIX(r_pix)]
  ENDFOR
ENDFOR


;Add the fake companions to the input image cube
;---------------------------------------------------
IF KEYWORD_SET(adi) then begin ; derotate the companion position in case of an ADI sequence
  paral = readfits(datadir+'vec_'+tg_name_bas+'_paral.fits')
	FOR i=0,nobj-1 DO BEGIN
    arr_tmp = fltarr(dim,dim)
  	IF right_handed EQ 1 then paral[i]=-paral[i]
    FOR j = 0, n_br-1 DO BEGIN
      theta = j*2*!DPI/n_br
      FOR k = 0, n_rad-1 DO BEGIN
        r_pix = DOUBLE(rad_c[k])/plsc
        xpos = r_pix * COS(theta-paral[i]*!DPI/180D0)
        ypos = r_pix * SIN(theta-paral[i]*!DPI/180D0)
        arr_tmp = arr_tmp + SHIFTI(fcp_arr, xpos, ypos) * sigma * level_r[FIX(r_pix)]
      ENDFOR
    ENDFOR
  	obj_tmp[*,*,i] = obj_tmp[*,*,i] + normal_vec[i]*arr_tmp
	ENDFOR
ENDIF ELSE BEGIN
	FOR i = 0, nobj-nref-1 DO obj_tmp[*,*,i] = obj_tmp[*,*,i] + normal_vec[i]*arr_final
ENDELSE


;Save the new image cube and ancillary data
;------------------------------------------
tg_name_dyn=tg_name_dyn + '_fcp'

writefits, procdir+'img_'+tg_name_dyn+'_dc.fits', obj_tmp, header
writefits, procdir+'img_'+tg_name_bas+'_fcparray.fits', Arr_final
writefits, procdir+'vec_'+tg_name_bas+'_coords.fits', coords

END
