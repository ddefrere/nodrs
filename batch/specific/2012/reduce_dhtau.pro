PRO REDUCE_DHTAU, REDUCE=reduce

; General running parameters
plot         = 1
info         = 1

; Reduction parameters
aper_rad     = 0.
bck_irad     = 45.
bck_orad     = 70.
ob_mode      = 2
fit_method   = 5
obstype      = 0.   
n_clip       = 512.
n_bin        = 2
nil_dic      = 'imaging'
tgt_name     = 'dh_tau'
flag         = 'SCI'
date         = '121105'
overlap      = 1
lmircam      = 1
skip_null    = 1
skip_flx     = 1
skip_adi     = 0
skip_cal     = 1

; Frame selection and derotation parameters
right_handed = 0
rin_init     = 0
verbose      = 0
n_coadd      = 1
sig_sl       = 3
sig_pos      = 0
x_range      = 0
y_range      = 0

IF KEYWORD_SET(reduce) THEN BEGIN
  ; Reduction parameters
  ;xc = [0] & yc = [0]
  
  ; Some files have no fits header, write the information0
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
  ;data_idx = [1872,2161]
  ;nod_idx  = [[1872,1881],[1902,1911],[1912,1921],[1942,1961],[1982,1991],[1992,2011],[2032,2041],[2042,2051],[0,0],[2102,2111],[2152,2161]]  
  data_idx = [1942,2161]
  nod_idx  = [[1942,1961],[1982,1991],[1992,2011],[2032,2041],[2102,2111],[2152,2161]]
  
  ; M files
  ;data_idx = [1882,2151]
  ;nod_idx  = [[1882,1891],[1892,1901],[1922,1931],[0,0],[1962,1971],[1972,1981],[2012,2021],[2022,2031],[2052,2061],[2112,2121],[2142,2151]]  

  ; Ls files
  ;data_idx = [2122,2201]
  ;nod_idx  = [[2122,2126],[2137,2141],[2162,2166],[2177,2181],[2182,2186],[2197,2201]] 

  ; Ice
  ;data_idx = [2127,2196]
  ;nod_idx  = [[2127,2131],[2132,2136],[2167,2171],[2172,2176],[2187,2191],[2192,2196]]
  
  LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx, BAD_IDX=bad_idx, $   ; Date to be reduced and index of files
            DATA_PATH=data_path, BCKG_PATH=bckg_path, DARK_PATH=dark_path, FLAT_PATH=flat_path, $                                                    ; Path to data directories (superseed "*_idx" files)
            XCEN=xcen, YCEN=ycen, OVERLAP=overlap, PSF_FILE=psf_file, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, $
            LAMBDA_CEN=lambda_cen, BANDWIDTH=bandwidth, $
            LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $                          ; Optional input keywords superseeding keywords in the header
            MASTERLOG=masterlog, ADI_MODe=adi_mode, BCKG_MODE=bckg_mode, FLX_MODE=flx_mode, FIT_MODE=fit_mode, NULL_MODE=null_mode, OB_MODE=ob_mode, $
            APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, BCK_CEN=bck_cen, $
            KEEP_RATIO=keep_ratio,  N_BIN=n_bin, N_CLIP=n_clip, N_FRBCK=n_frbck, N_FROB=n_frob, N_TRANS=n_trans, OFFSET=offset,$
            NO_CENTER=no_center, NO_DARK=no_dark, NO_FLAT=no_flat, RIN_INIT=rin_init, SIG_POS=sig_pos, SIG_SL=sig_pl, SKIP_CAL=skip_cal, SKIP_ADI=skip_adi, SKIP_FLX=skip_flx, SKIP_NULL=skip_null, $   ; Optional input keywords critical for data reduction
            LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                    ; Optional input keywords for output information
ENDIF                                            
END