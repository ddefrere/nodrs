;+
; NAME:
;   FPNEG_PCA_ADI
;   
; PURPOSE:
;   Uses the negative fake companion technique to determine the photometry and astrometry of a faint off-axis companion in an ADI data cube.
;   
; CALLING SEQUENCE:
;   Result = FPNEG_PCA_ADI(dim, klip, ZONE=zone, DELTA=delta, SMART=smart, FILTER=filter)
;                        
; DESCRIPTION:
;   This is a high-level routine calling a sequence of specialized routines to perform the (s)PCA analysis of an ADI cube
;   and derive the photometry/astrometry of a faint off-axis companion.
;   The main steps are the following:
;      - add a negative fake companion at a predefined position and with the predefined flux 
;      - image filtering of the whole cube (if needed)
;      - (s)PCA analysis of the cube
;      - computation of a figure of merit in a predefined region around the companion position
;      - use the simplex-ameoba algorithm to optimise this figure of merit by varying the position and flux of the negative fake companion
;   This routine has been developped for NACO/AGPM_L data sets. It could easily be adapted to other instruments.
;  
; INPUT VARIABLES:
;   dim     : dimensions of the image on which (s)PCA will be performed [pixels]
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
;  
; COMMON BLOCKS:
;   common_name : contains target name => implies that the data is available in datadir, in the form 'img_'+tg_name+'_dc.fits';
;                 the parallactic angle and photometry vectors (optional) should also be available in datadir, 
;                 in the form 'vec_'+tg_name+'_paral.fits' and 'vec_'+tg_name+'_photometry.fits'
;   pca_amoeba  : contains variables that need to be passed to the cost function in the amoeba optimisation process
; 
; OUTPUT:
;   A three-element vector containing the following variables:
;     - [0]: the best-fit astrometric position of the companion along the x dimension of the input images
;     - [1]: the best-fit astrometric position of the companion along the y dimension of the input images
;     - [2]: the best-fit photometry of the companion
;   
; DEPENDENCIES:
;    spca_adi_v13, pca_adi_v38, neg_fcp_v11, mfilter_v4
;
; MODIFICATION HISTORY:
;   Version 1.0, 30-JUL-2013, by Olivier Absil (ULg), adapted from a top-level routine developped by Dimitri Mawet for Keck pipeline
;-

; *********************************************************************************************************
; ****************************************** COST FUNCTION ************************************************
; *********************************************************************************************************

FUNCTION COST_FPNEG_PCA_ADI, fpneg

; fpneg : a 3-element vector containing (i) the x position [pixels], (ii) the y position [pixels] and (iii) the flux of the negative fake companion  

COMMON common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir
COMMON PCA_AMOEBA, plsc, dim_init, cx, cy, fwhm, right_handed, normal, camera_filter, cutoff_l, $
       dim_tmp, klip_tmp, smart_tmp, zone_tmp, delta_tmp, rin_init, step_init, mask_merit

tg_name_dyn = tg_name_ori

;Inject negative fake companion
NEG_FCP_V11, fpneg[0], fpneg[1], fpneg[2], fwhm, plsc, dim_tmp, camera_filter, right_handed, /ADI, /NACO, NORMAL=normal

;Spatial Filter
IF KEYWORD_SET(filter) THEN MFILTER_V4, highpass=cutoff_l, lowpass=fwhm/2, /median_f

; Classical or smart PCA on the cube containing the negative fake companion
IF KEYWORD_SET(smart_tmp) THEN BEGIN
  onezone = SQRT(fpneg[0]^2 + fpneg[1]^2) / fwhm ; angular separation of the negative fake companion in fwhm  
  pca_img_median = SPCA_ADI_V13(cx, cy, dim_tmp, dim_init, fwhm, plsc, klip_tmp, rin_init, step_init, delta_tmp, right_handed, ONEZONE=onezone, NORMAL=normal, /verbose)
ENDIF ELSE pca_img_median = PCA_ADI_V38(cx, cy, dim_tmp, dim_init, fwhm, klip_tmp, zone_tmp, rin_init, right_handed, NORMAL=normal)

figureofmerit = DOUBLE(TOTAL(ABS(pca_img_median*mask_merit*normal)^2)/TOTAL(mask_merit))

print, '***********************'
print, 'x, y, flux', fpneg[0], fpneg[1], fpneg[2]
print, 'sep, PA', SQRT(fpneg[0]^2+fpneg[1]^2)*plsc, ATAN(-fpneg[0],fpneg[1])*180/!DPI-104.8+360
print, 'figure of merit: ', figureofmerit
print, '***********************'

RETURN, figureofmerit
END


; *********************************************************************************************************
; ************************************* HIGH LEVEL FUNCTION ***********************************************
; *********************************************************************************************************

FUNCTION FPNEG_PCA_ADI, dim, klip, ZONE=zone, DELTA=delta, SMART=smart, FILTER=filter

COMMON common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir
COMMON PCA_AMOEBA, plsc, dim_init, cx, cy, fwhm, right_handed, normal, camera_filter, cutoff_l, $
       dim_tmp, klip_tmp, smart_tmp, zone_tmp, delta_tmp, rin_init, step_init, mask_merit


; Directory where the data are stored (datadir) and the generic name of the image cube (tg_name_ori)
; --------------------------------------------------------------------------------------------------
;datadir='/Users/dimitrimawet/data/AGPM/2013-01-31_betaPic_AGPM_collapsed/'
datadir='D:\Doc_Oli\IDL\Corono\NACO\'
;datadir='/home/absil/IDL/Corono/NACO/'
tg_name_ori = 'bpic_batch1_rec'
procdir=datadir
tg_name_dyn = tg_name_ori
tg_name_bas = tg_name_ori

; Parameters of the input image cube
; ----------------------------------
plsc = 0.0271 ; plate scale [arcsec]
diam = 8.2 ; telescope diameter [m]
lambda = 3.8e-6 ; filter central wavelength [m]
dim_init = 300 ; Input image size (not used, unless WCS data are considered)
cx = 150 ; xposition of the star in input image cube [pixels] 
cy = 150 ; y position of the star in input image cube [pixels]
right_handed = 0 ; 0 by default
normal = 77092.0 ; normalization by photometric variations (TBC)
fwhm = lambda/diam*206265/plsc ; fwhm of the psf, assuming no undersized Lyot stop [pixels]
camera_filter = 'apo'
cutoff_l = 4*fwhm ; parameter of high pass filter
dim_tmp = dim
klip_tmp = klip

;Crop array around center
MERGING_WINDOW_V5, cx, cy, cx, cy, dim_tmp, /resize_only
tg_name_ori = tg_name_dyn

; Specific PCA input parameters
; -----------------------------
IF KEYWORD_SET(smart) THEN smart_tmp = 1 ELSE smart_tmp = 0
IF NOT KEYWORD_SET(zone) THEN zone_tmp = dim ELSE zone_tmp = zone
IF NOT KEYWORD_SET(delta) THEN delta_tmp = 0.5 ELSE delta_tmp = delta
rin_init = 0.9 ; inner radius where PCA starts, in resolution elements (fwhm) 
step_init = 1.5 ; width of the annuli for sPCA, in resolution elements (fwhm)

; Simplex-amoeba parameters
; -------------------------
Ftol = 1e-3 ; fractional tolerance to be achieved in the function value
P0 = [11, 12, 1e-3] ; initial starting point for the optimisation of the three parameters
scale = [0.5, 0.5, 3e-4] ; caracteristic scale for each dimension

; Define the mask to compute the figure of merit
; ----------------------------------------------
osamp = 3 ; oversampling
mask1 = SHIFT(DIST(dim_tmp*osamp),dim_tmp/2*osamp,dim_tmp/2*osamp) ; create a mask centered on the star
inring = SQRT(P0[0]^2+P0[1]^2)*osamp - 1.22*fwhm*osamp/0.85 ; in pixels, taking into account the effet of the APO165 stop
outring = SQRT(P0[0]^2+P0[1]^2)*osamp + 1.22*fwhm*osamp/0.85 ; in pixels, taking into account the effet of the APO165 stop
mask1 = mask1 LE outring AND mask1 GT inring ; define an annulus
mask2 = SHIFT(DIST(dim_tmp*osamp),(dim_tmp/2+P0[0])*osamp,(dim_tmp/2+P0[1])*osamp) ; create a mask centered on the negative companion position
mask2 = mask2 LE 6*1.22*fwhm*osamp/0.85
mask_merit = REBIN(mask1 AND mask2, dim_tmp, dim_tmp)

result = AMOEBA(Ftol, FUNCTION_NAME='COST_FPNEG_PCA_ADI', P0=P0, SCALE=scale, FUNCTION_VALUE=fval, NMAX=1000)

RETURN, result
END