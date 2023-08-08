PRO BATCH_NOMIC, DATE=date, REDUCE=reduce, CALIBRATE=calibrate, MASTERLOG=masterlog, SKIP_CAL=skip_cal, FIT=fit, N_BIN=n_bin, REMOVE_OB=remove_ob, INFO=info, PLOT=plot, MOVIE=movie

; PURPOSE:
;   Batch routine to reduce LBTI/NOMIC data.
;
; MODIFICATION HISTORY:
;   Version 1.0,  24-NOV-2012, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  05-DEC-2012, DD: added keyword plot
;   Version 1.2,  25-JAN-2013, DD: updated call to 'nomic_null.pro' according to its new definition
;   Version 1.3,  25-APR-2013, DD: adapted to new reduction routine "lbti_drs.pro"
;   Version 1.4,  16-MAY-2013, DD: added keyword information
;   Version 1.4,  16-MAY-2013, DD: added DATE keyword
;   Version 1.5,  15-JAN-2014, DD: added MASTERLOG and SKIP_CAL keywords
;   Version 1.6,  15-FEB-2014, DD: added keyword REMOVE_OB

; Time the reduction and calibration parts
plot      = 1
movie     = 0
info      = 3 
log_file  = 1

IF KEYWORD_SET(reduce) THEN BEGIN
  ; General informatoin
  n_bin      = 0
  ;n_frob     = 1000
  aper_rad   = 15.
  bck_irad   = 20.
  bck_orad   = -1    ; -1 = go the channel edge (-1)
  method     = 0
  no_save    = 0
  n_clip     = 100
  ob_mode    = 0
  no_center  = 1
  fit_method = 1
  
  ; Define data path (26/02/2013)
  IF date EQ '130226' THEN BEGIN
    ; Filter information
    ubc_dxsp = 'Open'
    ubc_sxsp = 'Open'
    nic_fstp = 'Home'
    nic_beam = 'Trichroic'
    lmir_fw1 = 'N/A'
    lmir_fw2 = 'N/A'
    lmir_fw3 = 'N/A'
    lmir_fw4 = 'N/A'
    nom_fw1  = 'W-10145-9'
    nom_fw2  = 'Nprime'
    nom_apw  = 'Dual 3.05mm'
    pha_fw1  = 'N/A'
    pha_fw2  = 'N/A'
    pha_img  = 'N/A'
    nil_dic  = 'N/A'
    ;  lambda      = 10.9D-6
    ;  bandwidth   = 5.8D-6
    
    ; Reduce mu Uma
    data_idx = [102,4401]
    bckg_idx = [4403,4502]
    nod_idx  = data_idx
    ; dark_idx = [924,1453]
  
    ; Run computation
    xc = [76.5] & yc = [35] 
    LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced          
              DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories        
              XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME='mu_uma', FLAG='CAL', OBSTYPE=obstype, $
              LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $ ; Optional input keywords superseeding keywords in the header
              APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, METHOD=method, N_CLIP=n_clip, N_BIN=n_bin, N_FROB=n_frob,$          ; Optional input keywords critical for data reduction 
              OB_MODE=ob_mode, NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                        ; Optional input keywords for output information           
              
    
    ; Define data path (26/02/2013, alpa Boo)
    data_idx = [5503,7502]
    bckg_idx = [102,301]
    nod_idx  = data_idx
    
    ; Run computation
    xc = [33.] & yc = [21.] 
    LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced
              DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories        
              XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME='alp_boo', FLAG='SCI', OBSTYPE=obstype, $
              LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $              ; Optional input keywords superseeding keywords in the header
              APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, METHOD=method, N_CLIP=n_clip, N_BIN=n_bin, N_FROB=n_frob,$          ; Optional input keywords critical for data reduction 
              OB_MODE=ob_mode, NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                        ; Optional input keywords for output information      
        
    ; Define data path (26/02/2013, ups Boo)
    data_idx = [7503,8502]
    bckg_idx = [102,301]
    nod_idx  = data_idx
    
    ; Run computation
    xc = [32] & yc = [38.] 
    LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced
              DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories        
              XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME='ups_boo', FLAG='CAL', OBSTYPE=obstype, $
              LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2,$ ; Optional input keywords superseeding keywords in the header
              APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, METHOD=method, N_CLIP=n_clip, N_BIN=n_bin, N_FROB=n_frob,$          ; Optional input keywords critical for data reduction 
              OB_MODE=ob_mode, NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                        ; Optional input keywords for output information    
              
    ; Define data path (26/02/2013, Vega)
    data_idx = [25853,29651]
    bckg_idx = [29653,29852]
    nod_idx  = data_idx
  
    xc = [32.5] & yc = [32.] 
    LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced
              DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories        
              XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME='alf_Lyr', FLAG='SCI', OBSTYPE=obstype, $
              LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $; Optional input keywords superseeding keywords in the header
              APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, METHOD=method, N_CLIP=n_clip, N_BIN=n_bin, N_FROB=n_frob,$          ; Optional input keywords critical for data reduction 
              OB_MODE=ob_mode, NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                        ; Optional input keywords for output information    
   ENDIF
   IF date EQ '131019' THEN BEGIN
     ; Define data path
     ;data_idx = [28219,60273];,60273]
     ;bckg_folder = '/Volumes/nodrs/data/LBTI/nomic/131019/bckg_28ms/'
     data_idx = [83486,93485]
     bckg_idx = [93486,94485]
     nod_idx  = data_idx
     ;xc = [86] & yc = [60.5]
     n_frob   = 5000
     obstype  = 2
     
     ; Photometry
     ;data_idx = [0,9999]
     ;bckg_idx = [10000,10999]
     ;nod_idx  = data_idx
     ;xc = [32.0] & yc = [98.0] ; 28-34
     ;obstype  = 2
     ;n_frob   = 0
     ;n_clip   = 64
     
     LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced
               DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories
               XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME='alp_Per', FLAG='SCI', OBSTYPE=obstype, $
               LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $                               ; Optional input keywords superseeding keywords in the header
               APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, FIT_METHOD=FIT_method, N_CLIP=n_clip, N_BIN=n_bin, N_FROB=n_frob,$  ; Optional input keywords critical for data reduction
               NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                         ; Optional input keywords for output information
    ENDIF
    IF date EQ '131122' THEN BEGIN
      ; Define data path
      data_idx = [0,5999];,60273]
      nod_idx  = [[0,5060],[5061,5064],[5065,5999]]
      xc = [192] & yc = [119]
      fit_method = 5
      obstype    = 2
      
      LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced
                DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories
                XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME='NAC', FLAG='CAL', OBSTYPE=obstype, $
                LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $; Optional input keywords superseeding keywords in the header
                APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, FLX_METHOD=flx_method, N_CLIP=n_clip, FIT_METHOD=fit_method, N_BIN=n_bin, N_FROB=n_frob,$  ; Optional input keywords critical for data reduction
                NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                         ; Optional input keywords for output information
    ENDIF
    IF date EQ '131223' THEN BEGIN
      ; Define data path
      data_idx = [3600,9030];,60273]
      ;nod_idx  = data_idx
      bckg_idx = [9040,11039]
      xc = [197] & yc = [101]
      ;obstype    = 2
      
     ; LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced
     ;   DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories
     ;   XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME='iot_aur', FLAG='SCI', OBSTYPE=obstype, $
     ;   LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2,$ ; Optional input keywords superseeding keywords in the header
     ;   APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, FLX_METHOD=flx_method, N_CLIP=n_clip, FIT_METHOD=fit_method, N_BIN=n_bin, N_FROB=n_frob,$  ; Optional input keywords critical for data reduction
     ;   NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                         ; Optional input keywords for output information
        
      ; Define data path
      data_idx = [11041,16039];,60273]
      ;nod_idx  = data_idx
      bckg_idx = [16040,18000]
      ;xc = [197] & yc = [101]
      obstype    = 2
      
      LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced
                DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories
                XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME='bet_leo', FLAG='SCI', OBSTYPE=obstype, MASTERLOG=masterlog, $
                LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $ ; Optional input keywords superseeding keywords in the header
                APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, FLX_METHOD=flx_method, N_CLIP=n_clip, FIT_METHOD=fit_method, N_BIN=n_bin, N_FROB=n_frob,$  ; Optional input keywords critical for data reduction
                NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                         ; Optional input keywords for output information
    ENDIF
    IF date EQ '131230' THEN BEGIN
      ; Define data path
      data_idx = [2381,5380];,60273]
      nod_idx  = data_idx
      bckg_idx = [5381,6380]
      xc = [162] & yc = [190]
      n_frob    = 3000
      
      LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced
        DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories
        XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME='bet_leo', FLAG='SCI', OBSTYPE=obstype, $
        LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $; Optional input keywords superseeding keywords in the header
        APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, FLX_METHOD=flx_method, N_CLIP=n_clip, FIT_METHOD=fit_method, $
        N_BIN=n_bin, N_FROB=n_frob, SKIP_CAL=skip_cal, NO_CENTER=no_center, $                                                        ; Optional input keywords critical for data reduction
        NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                         ; Optional input keywords for output information
    ENDIF
;    IF date EQ '131231' THEN BEGIN
      ; Fix data number
      ; LBTI_FIXDATA, date

      ; Assign obstype and datatype
;      LBTI_FITSUTIL, date, [0,10666], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [100667,12666], OBSTYPE=0, DATATYPE=3
;      LBTI_FITSUTIL, date, [12667,16666], OBSTYPE=3, DATATYPE=3
;       
;      LBTI_FITSUTIL, date, [16667,21666], OBSTYPE=2, DATATYPE=0
;      LBTI_FITSUTIL, date, [21667,22166], OBSTYPE=3, DATATYPE=0
;      LBTI_FITSUTIL, date, [22167,27166], OBSTYPE=2, DATATYPE=0
;      LBTI_FITSUTIL, date, [27167,27666], OBSTYPE=3, DATATYPE=0
;      LBTI_FITSUTIL, date, [27667,32666], OBSTYPE=2, DATATYPE=0
;      LBTI_FITSUTIL, date, [32667,33166], OBSTYPE=3, DATATYPE=0
;      LBTI_FITSUTIL, date, [33167,38166], OBSTYPE=0, DATATYPE=0
;      LBTI_FITSUTIL, date, [38167,38666], OBSTYPE=3, DATATYPE=0
;      LBTI_FITSUTIL, date, [38667,43666], OBSTYPE=2, DATATYPE=0
;      LBTI_FITSUTIL, date, [43667,44166], OBSTYPE=3, DATATYPE=0
;      LBTI_FITSUTIL, date, [44167,49166], OBSTYPE=0, DATATYPE=0
;      LBTI_FITSUTIL, date, [49167,49666], OBSTYPE=3, DATATYPE=0  
;      LBTI_FITSUTIL, date, [49667,54666], OBSTYPE=2, DATATYPE=0     
;      LBTI_FITSUTIL, date, [54667,55166], OBSTYPE=3, DATATYPE=0, OBJNAME='gam Per'
;      
;      LBTI_FITSUTIL, date, [55167,60166], OBSTYPE=2, DATATYPE=3, OBJNAME='HD 14872', NOM_FW1='Nprime'
;      LBTI_FITSUTIL, date, [60167,60666], OBSTYPE=3, DATATYPE=3, OBJNAME='HD 14872', NOM_FW1='Nprime'
;      LBTI_FITSUTIL, date, [60667,65666], OBSTYPE=0, DATATYPE=3, OBJNAME='HD 14872', NOM_FW1='Nprime'
;      
;      LBTI_FITSUTIL, date, [65667,73718], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [73719,74218], OBSTYPE=3, DATATYPE=3
;      LBTI_FITSUTIL, date, [74219,79218], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [79219,79718], OBSTYPE=3, DATATYPE=3
;      LBTI_FITSUTIL, date, [79719,81718], OBSTYPE=0, DATATYPE=3
;      LBTI_FITSUTIL, date, [81719,82218], OBSTYPE=3, DATATYPE=3
; 
;      LBTI_FITSUTIL, date, [82219,87218], OBSTYPE=2, DATATYPE=0
;      LBTI_FITSUTIL, date, [87219,87718], OBSTYPE=3, DATATYPE=0
;      LBTI_FITSUTIL, date, [87719,92718], OBSTYPE=2, DATATYPE=0
;      LBTI_FITSUTIL, date, [92719,93218], OBSTYPE=3, DATATYPE=0
;      LBTI_FITSUTIL, date, [93219,95218], OBSTYPE=0, DATATYPE=0
;      LBTI_FITSUTIL, date, [95219,95718], OBSTYPE=3, DATATYPE=0
;      
;      LBTI_FITSUTIL, date, [95719,100718], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [100719,101218], OBSTYPE=3, DATATYPE=3
;      LBTI_FITSUTIL, date, [101219,106218], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [106219,106718], OBSTYPE=3, DATATYPE=3
;      LBTI_FITSUTIL, date, [106719,108718], OBSTYPE=0, DATATYPE=3
;      LBTI_FITSUTIL, date, [108719,109218], OBSTYPE=3, DATATYPE=3
;      
;      LBTI_FITSUTIL, date, [109219,114218], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [114219,114718], OBSTYPE=3, DATATYPE=3
;      LBTI_FITSUTIL, date, [114719,119718], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [119719,120218], OBSTYPE=3, DATATYPE=3
;      LBTI_FITSUTIL, date, [120219,122218], OBSTYPE=0, DATATYPE=3
;      LBTI_FITSUTIL, date, [122219,122718], OBSTYPE=3, DATATYPE=3
;      
;      LBTI_FITSUTIL, date, [122719,127718], OBSTYPE=2, DATATYPE=0
;      LBTI_FITSUTIL, date, [127719,128218], OBSTYPE=3, DATATYPE=0
;      LBTI_FITSUTIL, date, [128219,133218], OBSTYPE=2, DATATYPE=0
;      LBTI_FITSUTIL, date, [133219,133718], OBSTYPE=3, DATATYPE=0
;      LBTI_FITSUTIL, date, [133719,135718], OBSTYPE=0, DATATYPE=0
;      LBTI_FITSUTIL, date, [135719,136218], OBSTYPE=3, DATATYPE=0
;      
;      LBTI_FITSUTIL, date, [136219,141218], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [141219,141718], OBSTYPE=3, DATATYPE=3
;      LBTI_FITSUTIL, date, [141719,146718], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [146719,147218], OBSTYPE=3, DATATYPE=3
;      LBTI_FITSUTIL, date, [147219,149218], OBSTYPE=0, DATATYPE=3
;      LBTI_FITSUTIL, date, [149219,149718], OBSTYPE=3, DATATYPE=3
;      
;      LBTI_FITSUTIL, date, [149719,154718], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [154719,155218], OBSTYPE=3, DATATYPE=3
;      LBTI_FITSUTIL, date, [155219,160218], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [160219,160718], OBSTYPE=3, DATATYPE=3
;      LBTI_FITSUTIL, date, [160719,162718], OBSTYPE=0, DATATYPE=3
;      LBTI_FITSUTIL, date, [162719,163218], OBSTYPE=3, DATATYPE=3
;      
;      LBTI_FITSUTIL, date, [163219,168218], OBSTYPE=2, DATATYPE=0, OBJNAME='beta Leo', NOM_FW1='Nprime'
;      LBTI_FITSUTIL, date, [168219,168718], OBSTYPE=3, DATATYPE=0, OBJNAME='beta Leo', NOM_FW1='Nprime'
;      LBTI_FITSUTIL, date, [168719,173718], OBSTYPE=2, DATATYPE=0, OBJNAME='beta Leo', NOM_FW1='Nprime'
;      LBTI_FITSUTIL, date, [173719,174218], OBSTYPE=3, DATATYPE=0, OBJNAME='beta Leo', NOM_FW1='Nprime'
;      LBTI_FITSUTIL, date, [174219,179218], OBSTYPE=2, DATATYPE=0, OBJNAME='beta Leo', NOM_FW1='Nprime'
;      LBTI_FITSUTIL, date, [179219,179718], OBSTYPE=3, DATATYPE=0, OBJNAME='beta Leo', NOM_FW1='Nprime'
;      LBTI_FITSUTIL, date, [179719,181718], OBSTYPE=0, DATATYPE=0, OBJNAME='beta Leo', NOM_FW1='Nprime'
;      LBTI_FITSUTIL, date, [181719,182218], OBSTYPE=3, DATATYPE=0, OBJNAME='beta Leo', NOM_FW1='Nprime'
;      
;      LBTI_FITSUTIL, date, [182219,187218], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [187219,187718], OBSTYPE=3, DATATYPE=3
;      LBTI_FITSUTIL, date, [187719,192718], OBSTYPE=2, DATATYPE=3
;      LBTI_FITSUTIL, date, [192719,193718], OBSTYPE=3, DATATYPE=3
;      
;
;      LBTI_FITSUTIL, date, [193719,194118], OBSTYPE=0, DATATYPE=1, OBJNAME='Dome' 
;      LBTI_FITSUTIL, date, [194119,194518], OBSTYPE=0, DATATYPE=1, OBJNAME='Dome' 
;      LBTI_FITSUTIL, date, [194519,194918], OBSTYPE=0, DATATYPE=1, OBJNAME='Dome' 
;      LBTI_FITSUTIL, date, [194919,195318], OBSTYPE=0, DATATYPE=1, OBJNAME='Dome' 
                                                          
      ; Define data path
      ;data_idx = [6548,17047];,60273]
      ;nod_idx  = data_idx
      ;bckg_idx = [21048,22047]
      ;xc = [207] & yc = [182]
      ;n_frob    = 10500
      ;tgt_name = 'bet_and'
      
      ; Photometry frames
;      tgt_name = 'gam_per'
;      data_idx = [39048,56047];
;      nod_idx  = [[39048,39547],[39548,44547],[50048,50547],[50548,55547],[55548,56047]]
;       ;bckg_idx = [21048,22047]
;      xc = [200] & yc = [178]
;      n_frob     = 5000
;      fit_method = 5
;      obstype    = 2
;      
;      LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced
;        DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories
;        XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG='SCI', OBSTYPE=obstype, $
;        LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, NOM_APW=nom_apw, $ ; Optional input keywords superseeding keywords in the header
;        APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, FLX_METHOD=flx_method, N_CLIP=n_clip, FIT_METHOD=fit_method, N_BIN=n_bin, N_FROB=n_frob,$  ; Optional input keywords critical for data reduction
;        NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                         ; Optional input keywords for output information
      
;      tgt_name = 'gam_per'
;      data_idx = [23048,39547];
;      nod_idx  = [[23048,28047],[28048,28547],[28548,33547],[33548,34047],[34048,39047],[39048,39547]]
;      ;bckg_idx = [21048,22047]
;      xc = [200] & yc = [178]
;      n_frob    = 5000
;      obstype   = 1
;      
;      LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                                              ; Date to be reduced
;        DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                                       ; Path to data directories
;        XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG='SCI', OBSTYPE=obstype, $
;        LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, NOM_APW=nom_apw, $                    ; Optional input keywords superseeding keywords in the header
;        APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, FLX_METHOD=flx_method, N_CLIP=n_clip, FIT_METHOD=fit_method, N_BIN=n_bin, N_FROB=n_frob,$  ; Optional input keywords critical for data reduction
;        OB_MODE=ob_mode, NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                               ; Optional input keywords for output information
;        
;      tgt_name = 'gam_per2'
;      data_idx = [45048,61547];
;      nod_idx  = [[45048,50047],[50048,50547],[56048,61047],[61048,61547]]
;      ;bckg_idx = [21048,22047]
;      xc = [200] & yc = [178]
;      n_frob    = 5000
;      
;      LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced
;        DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories
;        XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG='SCI', OBSTYPE=obstype, $
;        LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, NOM_APW=nom_apw, $ ; Optional input keywords superseeding keywords in the header
;        APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, FLX_METHOD=flx_method, N_CLIP=n_clip, FIT_METHOD=fit_method, N_BIN=n_bin, N_FROB=n_frob,$  ; Optional input keywords critical for data reduction
;        OB_MODE=ob_mode, NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                         ; Optional input keywords for output information
      
;      data_idx = [168718,174217]
      ;nod_idx  = [[0,10666],[14667,16665],[10667,14666],[14667,16665],[16666,21665],[21666,22165],[22166,27165],[27166,27665],[27666,32665],[32666,33165],[33166,38165],$
      ;           [38166,38665],[38666,43665],[43666,44165],[44166,49165],[49166,49665],[49666,54665],[54666,55165],[55166,60165],[60166,60665],[60666,65665],[60166,60665]]    
      ;nod_idx  = [[55166,60165],[60166,60665],[60666,65665],[60166,60665]]
      ;nod_idx  = [[60166,60665],[60666,65665]];,[55166,60165],[60166,60665],[60666,65665],[60166,60665]]
      ;data_idx = [179718,193717]
      ;nod_idx  = [[179718,181717],[181718,182217],[182218,187217],[187218,187717],[187718,192717],[192718,193217],[193218,193717]] 
;      LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                                   ; Date to be reduced
;                DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                    ; Path to data directories
;                XCEN=xc, YCEN=yc, NO_OVERLAP=no_overlap, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, $
;                LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $                  ; Optional input keywords superseeding keywords in the header
;                APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, FLX_METHOD=flx_method, N_CLIP=n_clip, FIT_METHOD=fit_method, $
;                N_BIN=n_bin, N_FROB=n_frob, SKIP_CAL=skip_cal, NO_CENTER=no_center, $                                                            ; Optional input keywords critical for data reduction
;                NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                             ; Optional input keywords for output information
;    ENDIF
    ; Generic data reduction section
    IF KEYWORD_SET(date) THEN BEGIN
      ;data_idx=[12667,16666]
      ;nod_idx=[[12667,14665],[14666,16665]]
      data_idx = [6625,61330]
      LBTI_DRS, date, DATA_IDX=data_idx, NOD_IDX=nod_idx, BCKG_IDX=bckg_idx, $                                                               ; Date to be reduced
                DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories
                XCEN=xc, YCEN=yc, OVERLAP=overlap, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, $
                LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $                        ; Optional input keywords superseeding keywords in the header
                APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, FLX_METHOD=flx_method, N_CLIP=n_clip, FIT_METHOD=fit_method, N_BIN=n_bin, N_FROB=n_frob,$  ; Optional input keywords critical for data reduction
                OB_MODE=ob_mode, MASTERLOG=masterlog, SKIP_CAL=skip_cal, NO_CENTER=no_center,$
                NO_GPU=no_gpu, LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                         ; Optional input keywords for output information
    ENDIF
ENDIF

IF KEYWORD_SET(calibrate) THEN NULL_CALIB, date, CAL_METHOD=1, POLYDEG=0, REMOVE_OB=remove_ob, PLOT=plot
IF KEYWORD_SET(fit) THEN NULL_FIT, date, PLOT=plot, /PS
END
