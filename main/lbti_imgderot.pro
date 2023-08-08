;+
; NAME: LBTI_IMGDEROT
; 
; PURPOSE:
;   This procedure performs the ADI or PCA processing (depending on drs.adi_mode).
;
; INPUTS:
;   file        : Input FITS file with the image cube to process
;   data_file   : Data file with relevant information about the frames in file (in the first extension of file by default)
;
; KEYWORDS:
;   MEDIAN      : If set, the median of the selected image cube is removed from each image
;   PLOT        : Set this keyword to plot the data
;   VERBOSE     : Define the level of information printed to screen (see main routine for more info)
;
; MODIFICATION HISTORY:
;   Version 1.0, 02-FEB-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 14-JUN-2018, DDe: now save both ADI and MEDIAN images

PRO LBTI_IMGDEROT, file, data_file, MEDIAN=median, PSF_FILE=psf_file, PLOT=plot, VERBOSE=verbose

COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs, log
  
; Keyword check
IF NOT KEYWORD_SET(RIN_INIT) THEN rin_init = 0

; Read the data file  (or look in the first extension for the data)
IF KEYWORD_SET(DATA_FILE) THEN data_in = MRDFITS(data_file, 1, RANGE=range, /SILENT) ELSE data_in = MRDFITS(file, 1, RANGE=range, /SILENT)
img_in = READFITS(file, header, /SILENT)

; Compute psf's FWHM
psf_rad  = 1.028*FXPAR(HEADER, 'WAVELENG')/cnf.tel_diam     ; [rad], psf fwhm 
psf_fwhm = psf_rad*prm.r2m/(FXPAR(HEADER, 'PIXSCALE')*1D+3) ; [pixels], psf fwhm 

; Size of the input image
n_img = (size(img_in))[3]
dim   = (size(img_in))[2]

; Keep only a given X position on the detector
IF KEYWORD_SET(drs.X_RANGE) THEN BEGIN
  xcen = data_in.xcen
  idx_range = WHERE(xcen GT drs.x_range[0] AND xcen LE drs.x_range[1], n_range)
  IF n_range NE n_img THEN BEGIN
    img_in  = img_in[*,*,idx_range]
    data_in = data_in[idx_range]
  ENDIF
  PRINT, 'Number of rejected frames based on X position : ', n_img-n_range
  n_img = n_range
ENDIF

; Size of the input image
n_img = (size(img_in))[3]
dim   = (size(img_in))[2]

; Extract time and parallactic angles
time  = data_in.mjd_obs
time  = (time-time[0])*(24.*60.*60.)
paral = data_in.lbt_para
IF drs.right_handed THEN paral=-paral

; Perform ADI
IF drs.adi_mode EQ 0 THEN BEGIN  
  ; Compute median image 
  IF KEYWORD_SET(PSF_FILE) THEN BEGIN
    IF KEYWORD_SET(verbose)  THEN PRINT, 'Now performing PSF fitting'
    img_cen = PSF_FIT(img_in, 0.5*dim, 0.5*dim, 0, FIT_METHOD=6, PSF_FILE=psf_file, ZOOM_RATIO=2, BCK_FLX=bck_flx, BCK_ERR=bck_err, IMG_FIT=img_fit, $
                      STR_FLX=str_flx, STR_ERR=str_err, XCEN=xcen, YCEN=ycen, /NO_SHIFT)
    img_in  = img_cen - img_fit
    img_psf = 0   ; just because it's used later with other options of the code
  ENDIF ELSE BEGIN
    IF KEYWORD_SET(verbose)  THEN PRINT, 'Now performing median subtraction'
    IF KEYWORD_SET(MEDIAN) THEN RESISTANT_MEAN, img_in, 10, img_psf, DIMENSION=3, /SILENT  ;img_psf = MEDIAN(img_in, dim=3) 
  ENDELSE
   
  ; Image filter
  ;img_in=BANDPASS_FILTER(img_in, 1/(4.*psf_fwhm), 1., butterworth=1) 
  
  ; Derotate image and combine
  FOR i=0, n_img-1 DO BEGIN
    IF KEYWORD_SET(verbose)  THEN PRINT, 'Derotate by: ', paral[i]
    img_in[0,0,i]=ROT(img_in[*,*,i],-paral[i],1.0,dim/2,dim/2,CUBIC=-0.5,/PIVOT)
    ;img_in[0,0,i]=BANDPASS_FILTER(img_in[*,*,i], 1/(4.*psf_fwhm), 1., butterworth=1)
  ENDFOR

  ; Compute median of derotated image
  ;img_med = MEDIAN(img_in, dim=3)
  RESISTANT_MEAN, img_in, 10, img_med, DIMENSION=3, /SILENT
  ;img_con = CONVOL_FFT(img_med, psf)

  ; Mask center (for coronagraphy or saturated images)
  mask_t = SHIFT(DIST(dim),dim/2,dim/2)
  IF drs.rin_init THEN BEGIN
    mask   = mask_t ge (psf_fwhm*drs.rin_init)
    mask(where(mask eq 0))=0
    img_med = img_med * mask
  ENDIF

  ; Outer mask
  mask    = mask_t le dim/2
  mask(where(mask eq 0))=0
  img_med = img_med * mask

  ; Unsharp mask the data
  ;img_med = SIGMA_FILTER(img_med, 5, N_SIGMA=3, /ITERATE)
  ;img_med = UNSHARP_MASK(img_med, RADIUS=psf_fwhm)

  ; Write image to disk
  IF NOT KEYWORD_SET(no_save) THEN BEGIN
    ; First update header
    SXADDPAR, header, 'NAXIS3', 1
    ; Save files
    i_pos     = STRPOS(file, '_SEL.fits')
    WRITEFITS, STRMID(file, 0, i_pos) + '-DEROT-MED.fits', img_med, header
  ENDIF
   
  ; Now subtract PSF (or median)
  FOR i=0, n_img-1 DO BEGIN
    img_in[0,0,i]-=ROT(img_psf,-paral[i],1.0,dim/2,dim/2,CUBIC=-0.5,/PIVOT)
    ;img_in[0,0,i]=BANDPASS_FILTER(img_in[*,*,i], 1/(4.*psf_fwhm), 1., butterworth=1)   
  ENDFOR
  
  ; Compute median of derotated image
  ;img_med = MEDIAN(img_in, dim=3)
  RESISTANT_MEAN, img_in, 10, img_med, DIMENSION=3, /SILENT
  ;img_con = CONVOL_FFT(img_med, psf)
  
  ; Mask center (for coronagraphy or saturated images)
  mask_t = SHIFT(DIST(dim),dim/2,dim/2)
  IF drs.rin_init THEN BEGIN
    mask   = mask_t ge (psf_fwhm*drs.rin_init)
    mask(where(mask eq 0))=0
    img_med = img_med * mask
  ENDIF
  
  ; Outer mask
  mask    = mask_t le dim/2
  mask(where(mask eq 0))=0
  img_med = img_med * mask
   
  ; Unsharp mask the data
  ;img_med = SIGMA_FILTER(img_med, 5, N_SIGMA=3, /ITERATE)
  ;img_med = UNSHARP_MASK(img_med, RADIUS=psf_fwhm)
  
  ; Write image to disk
  IF NOT KEYWORD_SET(no_save) THEN BEGIN
    ; First update header
    SXADDPAR, header, 'NAXIS3', 1
    ; Save files
    i_pos     = STRPOS(file, '_SEL.fits')
    WRITEFITS, STRMID(file, 0, i_pos) + '-DEROT-ADI.fits', img_med, header
  ENDIF  
ENDIF

; Run PCA
IF drs.adi_mode EQ 1 THEN BEGIN
  ; Find data_path
  i_pos     = STRPOS(file, '/', /REVERSE_SEARCH)
  data_path = STRMID(file, 0, i_pos + 1)
  
  ; First write data in the format of Dimi's pipeline
  tgt_name = STRCOMPRESS(STRTRIM(STRING(FXPAR(header, 'OBJNAME'))), /REMOVE_ALL)
  writefits, data_path + 'vec_' + tgt_name + '_paral.fits', paral
  writefits, data_path + 'img_' + tgt_name + '_dc.fits', img_in
    
  ; Set common block for PCA routines
  COMMON common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir
  datadir      = data_path
  tg_name_ori  = tgt_name
  
  ; Run PCA
  CALL_PCA, dim, drs.pca_klip, ZONE=pca_zone, DELTA=delta, SMART=drs.pca_smart, FILTER=filter, NOCURVE=nocurve, /DISPLAY
ENDIF

END