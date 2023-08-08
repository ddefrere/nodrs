;+
; NAME:
;   CALL_PCA
;   
; PURPOSE:
;   Perform a (smart or classical) PCA analysis on an ADI image cube and derive detection limits (contrast curves).
;   
; CALLING SEQUENCE:
;   CALL_PCA, dim, klip, ZONE=zone, DELTA=delta, SMART=smart, FILTER=filter, NOCURVE=nocurve, DISPLAY=display
;                        
; DESCRIPTION:
;   This is a high-level routine calling a sequence of specialized routines to perform a full (s)PCA analysis of an ADI cube.
;   The main steps are the following:
;      - image filtering of the whole cube (if needed)
;      - (s)PCA analysis of the cube
;      - computation of the contrast curve
;      - injection of fake companions
;      - image filtering of the whole cube with fake companions (if needed)
;      - (s)PCA analysis of the cube with fake companions
;      - evaluation of the fake companion attenuation
;      - computation of the final contrast curve after correcting the attenuation
;  
; INPUT VARIABLES:
;   dim     : dimensions of the ouptut image, on which (s)PCA will be performed [pixels]
;   klip    : number of principal components to be used in the (s)PCA
;             
;   Additional inputs are hard-coded in the routine:
;     - datadir : the path to the directory where input data are stored  
;     - tg_name_ori : the routine will look in datadir for files named 
;                     * 'img_tg_name_ori_dc.fits': the input image cube
;                     * 'vec_tg_name_ori_paral.fits': parallactic angles associated to the images of the cube
;     - parameters of the input image cube (plsc, diam, lambda, dim_init, cx, cy, right_handed, normal)
;     - additional PCA parameters (rin_init, step_init) 
;   
; INPUT KEYWORDS:
;   zone    : [only in classical PCA] square zone on which the PCA will be perfomed [in pixels]. Default = dim. 
;   delta   : [only in smart PCA] parallactic angle exclusion region for building the PSF library in sPCA mode, defined as the minimum 
;             angular separation in units of $\lambda/D$ between current companion position other images of the cube. Default = 1*lambda/D.
;   smart   : if set, annulus-wise smart PCA will be used instead of classical PCA
;   filter  : if set, filtering will be applied on the image cube before PCA 
;   nocurve : if set, skip all the computations related to the contrast curves (including the 2nd PCA with fake companions) 
;   display : set this keyword to display intermediate and final results
;  
; COMMON BLOCK:
;   common_name : contains target name => implies that the data is available in datadir, in the form 'img_'+tg_name+'_dc.fits';
;                 the parallactic angle and photometry vectors (optional) should also be available in datadir, 
;                 in the form 'vec_'+tg_name+'_paral.fits' and 'vec_'+tg_name+'_photometry.fits'
; 
; OUTPUTS:
;   All outputs are in the form of FITS and postscript files dumped into the working directory. They include:
;     - tg_name_ori_rsz_fcp_pca_dc.fits : image cube after subtraction of the principal components
;     - img_tg_name_ori_rsz_fcp_pca_median.fits : final image obtained as the median of the derotated processed images of the cube
;     - vec_tg_name_ori_rsz_fcp_pca_contrast.fits : 
;     - plot_tg_name_ori_rsz_fcp_pca_contrast.ps : contrast curve obtained after the 1st sPCA pass 
;     - plot_tg_name_ori_rsz_fcp_pca_contrast_alpha.ps : final contrast curve taking into account fake companion attenuation 
;   
;  DEPENDENCIES:
;    merging_window_v5, spca_adi_v13, pca_adi_v38, contrast_curve_v71, fcp_v83, alpha_r_v61
;
; MODIFICATION HISTORY:
;   Version 1.0, 06-JUN-2013, by Olivier Absil (ULg), adapted from a top-level routine developped by Dimitri Mawet for Keck pipeline
;   Version 1.1, 26-JUL-2013, OA: upgraded various routines, implemented the 2D contrast and SNR maps
;-

PRO CALL_PCA, dim, klip, ZONE=zone, DELTA=delta, SMART=smart, FILTER=filter, DISPLAY=display, NOCURVE=nocurve

COMMON common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir


; Specify the directory where the data are stored (datadir) and the generic name of the image cube (tg_name_ori)
; --------------------------------------------------------------------------------------------------------------
;datadir='/Users/dimitrimawet/data/AGPM/2013-01-31_betaPic_AGPM_collapsed/'
;datadir='D:\Doc_Oli\IDL\Corono\NACO\'
datadir='/Volumes/LaCie/nodrs/results/lmircam/l1_fits/2016-10-16/'
;datadir='/home/absil/IDL/Corono/NACO/'
tg_name_ori = 'hd179218'
;tg_name_ori = 'bpic_batch1_rec'

procdir=datadir
tg_name_dyn = tg_name_ori
tg_name_bas = tg_name_ori
;spawn, 'rm '+datadir+'*.'+tg_name_bas+'_level_r.*'
;spawn, 'rm '+datadir+'*.'+tg_name_bas+'_coords.*'


; Parameters of the input image cube
; ----------------------------------
plsc = 0.010707   ; 0.01768    ; plate scale [arcsec]
diam = 8.4        ; telescope diameter [m]
lambda = 3.8e-6   ; 11.1e-6   ; filter central wavelength [m]
dim_init = 256    ; Input image size (not used, unless WCS data are considered)
cx = 128          ; xposition of the star in input image cube [pixels] 
cy = 128          ; y position of the star in input image cube [pixels]
right_handed = 0  ; 0 by default
normal = 77092.0  ; normalization by photometric variations (TBC)
fwhm = lambda/diam*206265/plsc ; fhwm of the psf, assuming no Lyot stop
camera_filter = 'Lp'

; Additional parameters
; ---------------------
sigma = 19        ; required significance level in contrast curves
cutoff_l = 4*fwhm ; parameter of high pass filter


; Initialise the output windows, if necessary
; -------------------------------------------
IF KEYWORD_SET(display) THEN BEGIN
  window, 0, xs=512, ys=512
  window, 2, xs=512, ys=512
  window, 3, xs=512, ys=512
  window, 1, xs=512, ys=512
ENDIF


; Specific PCA input parameters
; -----------------------------
IF NOT KEYWORD_SET(zone)  THEN zone=dim
IF NOT KEYWORD_SET(delta) THEN delta=0.5
rin_init  = 1.0 ; inner radius where PCA starts, in resolution elements (fwhm) 
step_init = 1.5 ; width of the annuli for sPCA, in resolution elements (fwhm)
;onezone = 4.7 ; if non-zero, use only one region for the sPCA, optimised at this location [in fwhm]
;box = 3 ; if non-zero, noise will be computed as the median of RMS in boxes of this width [in fwhm]


; Fake companion throughput calibration parameters
; ------------------------------------------------
sep_fcp = 4. ; radial angular separation between fake companions in the radial direction [in fwhm] 
n_rad = ROUND(dim/2/(sep_fcp*fwhm)) > 4  ;  number of fake companions along the radius (minimum 4)
n_br = 3 ; Number of branches of the pattern used to place our fake companions
pseudo_sigma_fcp = 20 ; SNR level at which the fake companions will be introduced
rad_c_in = (rin_init+0.5 > 2.0)*fwhm*plsc ; put the first FCP at 0.5 lambda/D from the PCA inner rim, but in any case beyond 2 lambda/D
rad_c_out_min = dim/2*plsc - 4.0*fwhm*plsc ; outer FCP at 4 lambda/D from the edge to prevent side effects
rad_c = rad_c_in + findgen(n_rad)*(rad_c_out_min-rad_c_in)/(n_rad-1)


;Run smart PCA on input cube (initial pass, no fake companion)
;-------------------------------------------------------------
;Crop array around center
MERGING_WINDOW_V5, cx, cy, cx, cy, dim, /resize_only

;Spatial Filter
IF KEYWORD_SET(filter) THEN MFILTER_V4, HIGHPASS=cutoff_l, LOWPASS=fwhm/2, /MEDIAN_F, DISPLAY=display

; Smart PCA without the fake companions
IF KEYWORD_SET(smart) THEN pca_img_median = SPCA_ADI_V13(cx, cy, dim, dim_init, fwhm, plsc, klip, rin_init, step_init, delta, right_handed, $
                                            NORMAL=normal, ONEZONE=onezone, /VERBOSE, DISPLAY=display) $
ELSE pca_img_median = PCA_ADI_V38(cx, cy, dim, dim_init, fwhm, klip, zone, rin_init, right_handed, NORMAL=normal, DISPLAY=display)

IF KEYWORD_SET(nocurve) THEN GOTO, THE_END

;Contrast curve (no normalization to alpha_r)
CONTRAST_CURVE_V71, pca_img_median, dim/2, dim/2, fwhm, plsc, dim, sigma, /GAUSSFILT, /AGPM, BOX=box, DISPLAY=display


;Second sPCA, inject fake companion at sigma * residual noise
;------------------------------------------------------------
;mask=readfits('mask_agpm2.fits')
;Restore original target name
tg_name_dyn = tg_name_ori+'_rsz'

;Injection of fake companions
FCP_V83, fwhm, pseudo_sigma_fcp, rad_c, n_br, plsc, dim, camera_filter, right_handed, $ 
         NORMAL=normal, /GAUSSFILT, /LBT, /ADI, BOX=box, /AGPM, DISPLAY=display

;Spatial Filter
IF KEYWORD_SET(filter) THEN MFILTER_V4, HIGHPASS=cutoff_l, LOWPASS=fwhm/2, /MEDIAN_F, DISPLAY=display

;PCA with fake companions
IF KEYWORD_SET(smart) THEN pca_fcp_img_median = SPCA_ADI_V13(cx, cy, dim, dim_init, fwhm, plsc, klip, rin_init, step_init, delta, right_handed, $
                                                NORMAL=normal, ONEZONE=onezone, /VERBOSE, DISPLAY=display) $
ELSE pca_fcp_img_median = PCA_ADI_V38(cx, cy, dim, dim_init, fwhm, klip, zone, rin_init, right_handed, NORMAL=normal, DISPLAY=display)

;Analyze fake companion throughput => alphar(r)
ALPHA_R_V61, fwhm, pca_img_median, n_br, rad_c, dim, DISPLAY=display

;Normalisation of the contrast curve : 5N(r)/alpha(r)
CONTRAST_CURVE_V71, 0, dim/2, dim/2, fwhm, plsc, dim, sigma, LEVEL_FILE=procdir+'vec_'+tg_name_ori+'_level_r.fits', $
                    /GAUSSFILT, /ALPHA, /AGPM, BOX=box

; Plot the 2D contrast map, based on a 3-pixel box
alpha_2d = READFITS(procdir+'img_'+tg_name_bas+'_alpha_r_2D.fits',header)
alpha_2d[WHERE(alpha_2d EQ 0)] = 1D0
noise_map = IMAGE_STDDEV(FILTER_IMAGE(pca_img_median, fwhm_gaussian=fwhm), 1.5*fwhm)/(alpha_2d>0.01)  ; use half the box size
PREP_PS & DEVICE, FILENAME='contrast_map.eps', /ENCAPS, /COLOR, BITS_PER_PIXEL=8, /PORTRAIT
LOADCT, 2, /SILENT 
TVIM, ALOG10(sigma*noise_map), /SCALE, XRANGE=[-dim/2,dim/2]*plsc, YRANGE=[-dim/2,dim/2]*plsc, $
      XTITLE='!4D!3x [arcsec]', YTITLE='!4D!3y [arcsec]', STITLE=STRING(sigma, FORMAT='(I0)')+' sigma detection limit (log)'
OPLOT, [-dim/2,dim/2]*plsc, [0,0]*plsc, LINESTYLE=1 & OPLOT, [0,0]*plsc, [-dim/2,dim/2]*plsc, LINESTYLE=1
LOADCT, 0, /SILENT 
DEVICE, /CLOSE & END_PS

; Plot the SNR in the image
PREP_PS & DEVICE, FILENAME='SNR_map.eps', /ENCAPS, /COLOR, BITS_PER_PIXEL=8
LOADCT, 2, /SILENT 
TVIM, FILTER_IMAGE(pca_img_median, fwhm_gaussian=fwhm)/noise_map, /SCALE, XRANGE=[-dim/2,dim/2]*plsc, YRANGE=[-dim/2,dim/2]*plsc, $
      XTITLE='!4D!3x [arcsec]', YTITLE='!4D!3y [arcsec]', STITLE='SNR', RANGE=[-5,5]
OPLOT, [-dim/2,dim/2]*plsc, [0,0]*plsc, LINESTYLE=1 & OPLOT, [0,0]*plsc, [-dim/2,dim/2]*plsc, LINESTYLE=1
DEVICE, /CLOSE
DEVICE, FILENAME='SNR_fcp_map.eps', /ENCAPS, /COLOR, BITS_PER_PIXEL=8
TVIM, FILTER_IMAGE(pca_fcp_img_median, fwhm_gaussian=fwhm)/noise_map, /SCALE, XRANGE=[-dim/2,dim/2]*plsc, YRANGE=[-dim/2,dim/2]*plsc, $
      XTITLE='!4D!3x [arcsec]', YTITLE='!4D!3y [arcsec]', STITLE='SNR', RANGE=[-5,5]
OPLOT, [-dim/2,dim/2]*plsc, [0,0]*plsc, LINESTYLE=1 & OPLOT, [0,0]*plsc, [-dim/2,dim/2]*plsc, LINESTYLE=1
LOADCT, 0, /SILENT 
DEVICE, /CLOSE & END_PS

THE_END:

;STOP
END
