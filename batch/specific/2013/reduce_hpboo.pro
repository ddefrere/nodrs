PRO REDUCE_HPBOO, REDUCE=reduce

; Running parameter
plot     = 1
info     = 3
log_file = 1
no_gpu   = 1.     ; Damned, gpuFFT does not work with the free version!!!
verbose  = 1

; Observation info
tgt_name = 'hp_boo'
flag     = 'SCI'
date     = '130524'
date_lng = '2013-05-24'

; Reduction parameters
aper_rad     = 0.
obstype      = 2.     ; 2 is for photometric data (1 for interferometric data)
n_clip       = 400
right_handed = 0
pca          = 0
fit_method   = 5

; Read and derotate image
DECLARE_PATH, pth, LMIRCAM=lmircam
data_path = pth.l1fits_path + date_lng + pth.sep

IF KEYWORD_SET(reduce) THEN BEGIN
  ; June data (HP Boo -- full frame).
  data_idx = [0,848]            ; [190,199]  ; nod1
  nod_idx  = [[0,29],[30,94],[95,159],[160,224],[225,289],[290,354],[355,419],[420,484],[485,549],[550,614],[680,744],[745,829],[830,848]]   ; [300,309]
  ;bckg_idx = [[0,1] ,[2,848]]
  
  ;LBTI_FITSUTIL, date, [0,848], CROP=[320,256,384,768], LMIRCAM=lmircam
  
  ;xc = [0] & yc = [0]
  bck_irad    = 35.
  bck_orad    = 55.
  no_overlap  = 1
  n_bin       = 2
  nil_dic     = 'imaging'
  LBTI_DRS, date, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$                            ; Date to be reduced
            DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                  ; Path to data directories
            XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, NIL_DIC=nil_dic, $                     ; Optional input keywords superseeding keywords in the header
            APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, FIT_METHOD=fit_method, NO_GPU=no_gpu, $   ; Optional input keywords critical for data reduction
            LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                          ; Optional input keywords for output information
ENDIF

; Number of files
files   = FILE_SEARCH(data_path + date_lng + '_' + flag + '_' + tgt_name + '*FULL_IMG.fits')
n_files = N_ELEMENTS(files)

; Frame selection
n_coadd = 10
sig_sl  = 3
sig_pos = 3
;y_range = [256,512]
LBTI_MERGEL1IMG, files, NCOADD=n_coadd, DATA_FILE=data_files[i_f], RANGE=range, Y_RANGE=y_range, SIGMA_SLOPE=sig_sl, SIGMA_POSITION=sig_pos, PLOT=plot;, /MEDIAN

; Create file for Dimi's pipeline
dimi_file = FILE_SEARCH(data_path + date_lng + '*' + flag + '*' + tgt_name + '*IMG_SEL.fits')
data      = MRDFITS(dimi_file, 1, /SILENT)
objt      = READFITS(dimi_file)

; Compute median
img_med = MEDIAN(objt, dim=3)

; Derotate image
nobj = (size(objt))[3]
dim  = (size(objt))[2]
FOR i=0, nobj-1 DO BEGIN
  if keyword_set(verbose) then print, 'Derotate by: ', paral[i]
  if right_handed eq 1 then paral[i]=-paral[i]
  objt[*,*,i]=rot(objt[*,*,i]-img_med,-paral[i],1.0,dim/2,dim/2,cubic=-0.5,/pivot)
ENDFOR

; Compute median of derotated image
img_med = MEDIAN(objt, dim=3)
;img_con = CONVOL_FFT(img_med, psf)

;Mask center (for coronagraphy or saturated images)
rin_init     = 1.
fwhm         = 31.0
IF rin_init NE 0 THEN BEGIN
  mask_t=shift(dist(dim),dim/2,dim/2)
  mask= mask_t ge (fwhm*rin_init)
  mask(where(mask eq 0))=0
  img_med = img_med * mask
  img_con = img_con * mask
ENDIF

; Plot if requested
IF plot GT 1 THEN PLOTXY, img_med
;IF plot GT 1 THEN PLOTXY, img_con

; Write image to disk
IF NOT KEYWORD_SET(no_save) THEN BEGIN
  WRITEFITS, data_path + date_lng + '_' + flag + '_' + tgt_name + '_MED_DEROT.fits', img_med
;  WRITEFITS, data_path + date_lng + '_' + flag + '_' + tgt_name + '_CON_DEROT.fits', img_con
ENDIF
END