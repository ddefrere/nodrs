PRO REDUCE_EPSERI, REDUCE=reduce

; Running parameter
plot         = 1
info         = 3
aper_rad     = 0.
obstype      = 2.     ; 2 is for photometric data (1 for interferometric data)
no_gpu       = 1.     ; Damned, gpuFFT does not work with the free version!!!
n_clip       = 300
tgt_name     = 'eps_eri'
flag         = 'SCI'
date         = '131020'
date_lng     = '2013-10-20'
lmircam      = 0
right_handed = 0
pca          = 0
rin_init     = 0
verbose      = 1.
ob_mode      = 2
fit_method   = 5

; Read and derotate image
DECLARE_PATH, pth, LMIRCAM=lmircam
data_path = pth.l1fits_path + date_lng + pth.sep
psf_path  = pth.l1fits_path + '2012-11-05' + pth.sep

IF KEYWORD_SET(reduce) THEN BEGIN
  ; reduction parameters
  ; xc = [279] & yc = [297] ; 255
  bck_irad    = 15.
  bck_orad    = 25.
  no_overlap  = 1
  n_bin       = 2
  nil_dic     = 'imaging' 
  
  ; Pre-crop the frames
  ;LBTI_FITSUTIL, date, [10,19730], CROP=[256,256,512,512], LMIRCAM=lmircam
      
  ; Dark and flat frames
  data_idx = [2035,12034]  ;18555
  flag1    = 'SCI_P1'
;  LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$                                ; Date to be reduced
;            DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                                       ; Path to data directories
;            XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag1, OBSTYPE=obstype, NIL_DIC=nil_dic, $                                          ; Optional input keywords superseeding keywords in the header
;            APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, FIT_METHOD=fit_method, NO_GPU=no_gpu, OB_MODE=ob_mode,   $     ; Optional input keywords critical for data reduction
;            LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                               ; Optional input keywords for output information
;            
;  ; Dark and flat frames
;  data_idx = [12035,22034]  ;18555
  flag2    = 'SCI_P2'
;  LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$                                ; Date to be reduced
;            DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                                       ; Path to data directories
;            XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag2, OBSTYPE=obstype, NIL_DIC=nil_dic, $                                          ; Optional input keywords superseeding keywords in the header
;            APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, FIT_METHOD=fit_method, NO_GPU=no_gpu, OB_MODE=ob_mode,   $     ; Optional input keywords critical for data reduction
;            LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                               ; Optional input keywords for output information
            
    ; Dark and flat frames
;  data_idx = [22035,32034]  ;18555
  flag3    = 'SCI_P3'
;  LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$                                ; Date to be reduced
;            DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                                       ; Path to data directories
;            XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag3, OBSTYPE=obstype, NIL_DIC=nil_dic, $                                          ; Optional input keywords superseeding keywords in the header
;            APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, FIT_METHOD=fit_method, NO_GPU=no_gpu, OB_MODE=ob_mode,   $     ; Optional input keywords critical for data reduction
;            LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                               ; Optional input keywords for output information
;            
;  ; Dark and flat frames
;  data_idx = [32035,39034]  ;18555
  flag4    = 'SCI_P4'
;  LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$                                ; Date to be reduced
;            DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                                       ; Path to data directories
;            XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag4, OBSTYPE=obstype, NIL_DIC=nil_dic, $                                          ; Optional input keywords superseeding keywords in the header
;            APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, FIT_METHOD=fit_method, NO_GPU=no_gpu, OB_MODE=ob_mode,   $     ; Optional input keywords critical for data reduction
;            LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                               ; Optional input keywords for output information
;            
;  ; Dark and flat frames
;  data_idx = [39035,46556]  ;18555
  flag5    = 'SCI_P5'
;  LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$                                ; Date to be reduced
;            DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                                       ; Path to data directories
;            XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag5, OBSTYPE=obstype, NIL_DIC=nil_dic, $                                          ; Optional input keywords superseeding keywords in the header
;            APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, FIT_METHOD=fit_method, NO_GPU=no_gpu, OB_MODE=ob_mode,   $     ; Optional input keywords critical for data reduction
;            LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                               ; Optional input keywords for output information
            
  ; Merge the two parts
  file1 = FILE_SEARCH(data_path + date_lng + '*' + flag1 + '*' + tgt_name + '*FULL_IMG.fits')
  file2 = FILE_SEARCH(data_path + date_lng + '*' + flag2 + '*' + tgt_name + '*FULL_IMG.fits')
  file3 = FILE_SEARCH(data_path + date_lng + '*' + flag3 + '*' + tgt_name + '*FULL_IMG.fits')
  file4 = FILE_SEARCH(data_path + date_lng + '*' + flag4 + '*' + tgt_name + '*FULL_IMG.fits')
  file5 = FILE_SEARCH(data_path + date_lng + '*' + flag5 + '*' + tgt_name + '*FULL_IMG.fits')
  LBTI_MERGEL1FILES, [file1,file2,file3,file4,file5]
ENDIF

; Number of files
files      = FILE_SEARCH(data_path + date_lng + '_' + flag + '_' + tgt_name + '*FULL_IMG.fits')
data_files = FILE_SEARCH(data_path + date_lng + '_' + flag + '_' + tgt_name + '*FULL_DATA.fits')
n_files = N_ELEMENTS(files)

; Do PSF
psf   = READFITS(FILE_SEARCH(data_path + date_lng +  '_PSF_' + tgt_name + '*_MED.fits'))

; Loop over the number of files
FOR i_f = 0, n_files-1 DO BEGIN
  ; --- Perform ADI/PCA if requested
  IF KEYWORD_SET(pca) THEN BEGIN
    ; COMMON block for the PCA routine
    common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir
    dim          = 256
    dim_init     = 256
    fwhm         = 8.5/2.
    truncate_pca = 8
    Na           = 32
    normal       = 1
    datadir      = data_path
    procdir      = data_path
    tg_name_dyn  = tgt_name
    tg_name_bas  = tgt_name
    wcs          = 0
   ; status       = PCA_ADI_V35(cx_init, cy_init, dim, dim_init, fwhm, truncate_pca, Na, rin_init, right_handed, normal=normal, wcs=wcs, /Display, /verbose)
  ENDIF ELSE BEGIN
    ; Frame selection
    n_coadd = 20
    sig_sl  = 3
    sig_pos = 3
   ; LBTI_MERGEL1IMG, files[i_f], NCOADD=n_coadd, DATA_FILE=data_files[i_f], RANGE=range, SIGMA_SLOPE=sig_sl, SIGMA_POSITION=sig_pos, PLOT=plot;, /MEDIAN

    ; Create file for Dimi's pipeline
    dimi_file = FILE_SEARCH(data_path + date_lng + '*' + flag + '*' + tgt_name + '*IMG_SEL.fits')
    data      = MRDFITS(dimi_file, 1, /SILENT)
    objt      = READFITS(dimi_file)
  
    ; Find corresponding PARAL file
    paral = data.lbt_para
    
    ; Subtract median and derotate image
    nobj = (size(objt))[3]
    dim  = (size(objt))[2]
    FOR i=0, nobj-1 DO BEGIN
      if keyword_set(verbose) then print, 'Derotate by: ', paral[i]
      if right_handed eq 1 then paral[i]=-paral[i]
      objt[*,*,i]=rot(objt[*,*,i],-paral[i],1.0,dim/2,dim/2,cubic=-0.5,/pivot)
    ENDFOR
    
    ; Compute median of derotated image
    img_med = MEDIAN(objt, dim=3)
    ;img_con = CONVOL_FFT(img_med, psf)
    
    ;Mask center (for coronagraphy or saturated images)
    fwhm         = 12.0
    IF rin_init NE 0 THEN BEGIN
      mask_t=shift(dist(dim),dim/2,dim/2)
      mask= mask_t ge (fwhm*rin_init)
      mask(where(mask eq 0))=0
      img_med = img_med * mask
      ;img_con = img_con * mask
    ENDIF
    
    ; Plot if requested
    IF plot GT 1 THEN PLOTXY, img_med
    IF plot GT 1 THEN PLOTXY, img_con
    
    ; Write image to disk
    IF NOT KEYWORD_SET(no_save) THEN BEGIN
      i_pos = STRPOS(files[i_f], 'IMG')
      file  = STRMID(files[i_f], 0, i_pos)
      WRITEFITS, file + 'MED_DEROT.fits', img_med
      ;WRITEFITS, file + 'CON_DEROT.fits', img_con
    ENDIF
  ENDELSE
ENDFOR

END