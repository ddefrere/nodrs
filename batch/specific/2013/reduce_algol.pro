PRO REDUCE_ALGOL, REDUCE=reduce

; Infos
; I measured 3.5 mas of RMS tip/tilt on the PSF data

; Running parameter
plot         = 0
info         = 3
aper_rad     = 0.
obstype      = 2.     ; 2 is for photometric data (1 for interferometric data)
no_gpu       = 1.     ; Damned, gpuFFT does not work with the free version!!!
no_overlap   = 1
n_bin        = 2
n_clip       = 256
ob_mode      = 2
tgt_name     = 'Algol'
psf_name     = tgt_name
flag         = 'SCI'
date         = '131024'
date_lng     = '2013-10-24'
right_handed = 0
rin_init     = 0.2
no_center    = 0
log_file     = 1
fit_method   = 5            ; use Moffat model fitting
verbose      = 1

; Declare path
DECLARE_PATH, pth, LMIRCAM=lmircam
data_path = pth.l1fits_path + date_lng + pth.sep

IF KEYWORD_SET(reduce) THEN BEGIN
  ; reduction parameters
  bck_irad = 15
  bck_orad = 25
 
  ; Pre-crop the frames
  IF date EQ '131024' THEN BEGIN
    ; LBTI_FITSUTIL, date, [0,1040], CROP=[256,256,512,768], LMIRCAM=lmircam
      
    ; PSF data
    data_idx = [830,1819]
    ;nod_idx  = [[0,17],[18,40],[41,60],[61,260],[261,287],[288,524],[525,680],[681,700],[701,800],[801,820],[821,920],[921,940],[941,1040]]
    LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx, BAD_IDX=bad_idx, $   ; Date to be reduced
              DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                            ; Path to data directories
              XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, NIL_DIC=nil_dic, $                              ; Optional input keywords superseeding keywords in the header
              APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, OFFSET=offset, OB_MODE=ob_mode, NO_GPU=no_gpu, $
              FIT_METHOD=fit_method, NO_CENTER=no_center, SKIP_CAL=skip_cal,  $                                                                        ; Optional input keywords critical for data reduction
              LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                    ; Optional input keywords for output information    
   ENDIF  
ENDIF

; Files and bumber of files
files   = FILE_SEARCH(data_path + date_lng + '*' + flag + '*' + tgt_name + '*FULL_IMG.fits')
n_files = N_ELEMENTS(files)

; Data files
data_files = FILE_SEARCH(data_path + date_lng + '*' + flag + '*' + tgt_name + '*FULL_DATA.fits')

; PSF file
psf_files  = FILE_SEARCH(data_path + date_lng +  '*PSF*' + psf_name + '*FULL_MED.fits')
  
; Loop over the number of files
FOR i_f = 0, n_files-1 DO BEGIN  
  ; PSF correlation
  ;img     = READFITS(files[i_f], hdr_img)
  ;img_psf = READFITS(psf_files[i_f])
  ;nobj    = (size(img))[3]
  ;dim     = (size(img))[2]
  ;img_con = DBLARR(dim,dim,nobj)
  ;FOR i=0, nobj-1 DO BEGIN
  ;  img_con[0,0,i] = CONVOL_FFT(img[*,*,i], img_psf)
  ;ENDFOR
  ;p=STRPOS(files[i_f],'.fits')
  ;IF p GT 0 THEN files[i_f]=STRMID(files[i_f], 0, p) + '_CON.fits' ELSE MESSAGE, 'Invalid file format.'
  ;MWRFITS, img_con, files[i_f], hdr_img, /CREATE, /SILENT
  
  ; Frame selection
  n_coadd = 5
  sig_sl  = 10
  sig_pos = 10
  ;range   = [600,1000]
  ;y_range = [256,768] 
  LBTI_MERGEL1IMG, files[i_f], NCOADD=n_coadd, DATA_FILE=data_files[i_f], RANGE=range, Y_RANGE=y_range, SIGMA_SLOPE=sig_sl, SIGMA_POSITION=sig_pos, PLOT=plot;, /MEDIAN
  
  ; Read image file
  dimi_file = FILE_SEARCH(data_path + date_lng + '*' + flag + '*' + tgt_name + '*IMG_SEL.fits')
  objt      = READFITS(dimi_file)

  ; Find corresponding PARAL file
  i_pos = STRPOS(files[i_f], 'IMG')
  file  = STRMID(files[i_f], 0, i_pos)
  med   = READFITS(file + 'MED.fits')
  
  ; Extract parralactic angle
  data  = MRDFITS(dimi_file, 1, /SILENT)
  paral = data.lbt_para
  
  ; Correct parralactic angle
  idx_par = WHERE(paral LT 0)
  paral[idx_par] = paral[idx_par] + 360
  
 ; Create file for Dimi's pipeline
 data         = MRDFITS(dimi_file, 1, /SILENT)
 writefits, data_path + 'vec_' + tgt_name + '_paral.fits', data.lbt_para
 ;FILE_MOVE, dimi_file, data_path + 'img_' + tgt_name + '_dc.fits', /OVERWRITE, /ALLOW_SAME
  
 ; Run PCA
 ;dim          = 128
 ;klip         = 10
 ;zone         = 128
 ;smart        = 0
 ;CALL_PCA, dim, klip, ZONE=zone, DELTA=delta, SMART=smart, FILTER=filter, NOCURVE=nocurve, /DISPLAY
    
  ; Subtract median and derotate image
  nobj = (size(objt))[3]
  dim  = (size(objt))[2]
  FOR i=0, nobj-1 DO BEGIN
    if keyword_set(verbose) then print, 'Derotate by: ', paral[i]
    if right_handed eq 1 then paral[i]=-paral[i]
    objt[*,*,i]=rot(objt[*,*,i]-med, -paral[i],1.0,dim/2,dim/2,cubic=-0.5,/pivot)
  ENDFOR
  
  ; Compute median of derotated image
  img_med = MEDIAN(objt, dim=3)
  img_con = CONVOL_FFT(img_med, med)
  
  ;Mask center (for coronagraphy or saturated images)
  IF rin_init NE 0 THEN BEGIN
    mask_t=shift(dist(dim),dim/2,dim/2)
    mask= mask_t ge (rin_init)
    mask(where(mask eq 0))=0
    img_med = img_med * mask
    img_con = img_con * mask
  ENDIF
  
  ; Plot if requested
  IF plot GT 1 THEN PLOTXY, img_med
  ;IF plot GT 1 THEN PLOTXY, img_con
  
  ; Write image to disk
  IF NOT KEYWORD_SET(no_save) THEN BEGIN
    WRITEFITS, file + 'MED_DEROT.fits', img_med
    WRITEFITS, file + 'CON_DEROT.fits', img_con
  ENDIF
ENDFOR

END