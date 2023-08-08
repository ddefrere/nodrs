PRO REDUCE_FUTAU, REDUCE=reduce, ADI=adi, PCA=pca

; General running parameters
plot         = 1
info         = 1

; Reduction parameters
aper_rad     = 0.
bck_irad     = 45.
bck_orad     = 70.
fit_method   = 5
obstype      = 2.     ; 2 is for photometric data (1 for interferometric data)
no_gpu       = 1.     ; Damned, gpuFFT does not work with the free version!!!
n_clip       = 1024
n_bin        = 2
no_overlap   = 1
nil_dic      = 'imaging'
tgt_name     = 'fu_tau'
flag         = 'SCI'
date         = '121105'
date_lng     = '2012-11-05'
lmircam      = 1

; Frame selection and derotation parameters
right_handed = 0
rin_init     = 1
verbose      = 1
n_coadd      = 1
sig_sl       = 3
sig_pos      = 0
x_range      = 0
y_range      = 0

IF KEYWORD_SET(reduce) THEN BEGIN
  ; Reduction parameters
  ;xc = [0] & yc = [0]
  
  ; Some files have no fits header, write the information
  ;lmir_fw1 = ' LargeDualAperturte'
  ;lmir_fw2 = ' SX-Half-moon'
  ;lmir_fw3 = ' Lp'
  ;lmir_fw4 = ' Open'
  ;LBTI_FITSUTIL, date, [1992,2011], LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $
  ;               NOM_APW=nom_apw, PHA_FW1=pha_fw1, PHA_FW2=pha_fw2, PHA_IMG=pha_img, NIL_DIC=nil_dic, LMIRCAM=lmircam
  ;lmir_fw1 = ' LargeDualAperturte'
  ;lmir_fw2 = ' L-cont4'
  ;lmir_fw3 = ' Open'
  ;lmir_fw4 = ' Open'
  ;LBTI_FITSUTIL, date, [0,199], LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $
  ;               NOM_APW=nom_apw, PHA_FW1=pha_fw1, PHA_FW2=pha_fw2, PHA_IMG=pha_img, NIL_DIC=nil_dic, LMIRCAM=lmircam
  
  ; Read linearity data
  ;LBTI_DRS, '121105', DATA_IDX=[0,199], /LMIRCAM, /NO_CENTER, XCEN=548, YCEN=512, TGT_NAME='sky', APER_RAD=2, BCK_IRAD=3, BCK_ORAD=15, PLOT=3, FLAG='BCK', OBSTYPE=2
  
  ; Crop all frames
  ;LBTI_FITSUTIL, date, [1872,2521], CROP=[0,148,1024,728], LMIRCAM=lmircam
  ;LBTI_FITSUTIL, date, [180,199], CROP=[0,148,1024,728], LMIRCAM=lmircam
  ;
  ; Dark and flat frames
  flat_idx = [190,199]
  ;dark_idx = [3000,3099]  ; CDS data (dark automatically removed)
  
  ; Nov 2012 data (DH Tau -- full frame).
  ; Reference object
  ; tgt_name = 'hd_18881'
  ; flag     = 'CAL'
  ; Lp files
  ; data_idx = [2202,2521]
  ; nod_idx  = [[2202,2221],[2342,2361],[2362,2381],[2502,2521]]
  ; M files
  ; data_idx = [2202,2521];[1872,2201]
  ; nod_idx  = [[2222,2241],[2322,2341],[2382,2401],[2481,2501]]
  ;LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$                ; Date to be reduced
  ;  DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                               ; Path to data directories
  ;  XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, NIL_DIC=nil_dic, $                                  ; Optional input keywords superseeding keywords in the header
  ;  APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, NO_GPU=no_gpu, FIT_METHOD=fit_method,   $                                    ; Optional input keywords critical for data reduction
  ;  LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                       ; Optional input keywords for output information
  
  ; Lp files
  data_idx = [2545,2649]
  nod_idx  = [[2545,2549],[2590,2594],[2595,2599],[2645,2649]]
  ob_mode  = 2 
  LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$                        ; Date to be reduced
    DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                               ; Path to data directories
    XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, NIL_DIC=nil_dic, $                                  ; Optional input keywords superseeding keywords in the header
    APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, OB_MODE=ob_mode, NO_GPU=no_gpu, FIT_METHOD=fit_method,   $                                    ; Optional input keywords critical for data reduction
    LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                       ; Optional input keywords for output information
  ;
  ; M files
  data_idx = [2545,2649]
  nod_idx  = [[2550,2559],[2580,2589],[2600,2609],[2635,2644]]
  ob_mode  = 2
  LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$                         ; Date to be reduced
    DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                               ; Path to data directories
    XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, NIL_DIC=nil_dic, $                                  ; Optional input keywords superseeding keywords in the header
    APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, OB_MODE=ob_mode, NO_GPU=no_gpu, FIT_METHOD=fit_method,   $                                    ; Optional input keywords critical for data reduction
    LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                       ; Optional input keywords for output information
    
  ; Ls files
  data_idx = [2545,2649]
  nod_idx  = [[2560,2564],[2575,2579],[2610,2614],[2630,2634]] 
  ob_mode  = 0
  LBTI_DRS, date,  LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$                       ; Date to be reduced
    DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                               ; Path to data directories
    XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, NIL_DIC=nil_dic, $                                  ; Optional input keywords superseeding keywords in the header
    APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, OB_MODE=ob_mode, NO_GPU=no_gpu, FIT_METHOD=fit_method,   $                                    ; Optional input keywords critical for data reduction
    LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                       ; Optional input keywords for output information
    
  ; Ice
  data_idx = [2545,2649]
  nod_idx  = [[2565,2570],[2571,2574],[2615,2624],[2625,2629]] 
  ob_mode  = 2
  LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$     ; Date to be reduced
    DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                               ; Path to data directories
    XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, NIL_DIC=nil_dic, $                                  ; Optional input keywords superseeding keywords in the header
    APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, OB_MODE=ob_mode, NO_GPU=no_gpu, FIT_METHOD=fit_method,   $                                    ; Optional input keywords critical for data reduction
    LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                       ; Optional input keywords for output information
ENDIF

IF KEYWORD_SET(ADI) OR KEYWORD_SET(PCA) THEN LBTI_IMGPROCESS, date_lng, ADI=adi, PCA=pca, N_COADD=n_coadd, SIG_SL=sig_sl, SIG_POS=sig_pos, X_RANGE=x_range, Y_RANGE=y_range, RIN_INIT=rin_init, $
                                             TGT_NAME=tgt_name, LMIRCAM=lmircam, PLOT=plot, INFO=info
END