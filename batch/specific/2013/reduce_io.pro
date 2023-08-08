PRO REDUCE_IO, REDUCE=reduce, MASTERLOG=masterlog
  
  ; WARNING. DATA from DEc. 24th 2013 have a lot of empty header fields (e.g., LBT_PARA, LBT_LST, etc)
  ; 
  
  ; General running parameters
  plot      = 1         ; Plot data on screen if larger than 1
  info      = 3
  log_file  = 1
  no_save   = 0
  
  ; Specific to Io data
  tgt_name     = 'io'
  flag         = 'SCI'
  date         = '131224'
  obstype      = 1        ; Coherent imaging
  skip_cal     = 1
  skip_flx     = 1
  skip_null    = 1
  skip_adi     = 0
  
  ; reduction parameters
  bckg_mode    = 0
  ob_mode      = 2
  fit_mode     = 0      ; Fit method gives the best results (more robust) but it's 50x more time consuming!
  flx_mode     = 0      ; 0 for aperture photometry and 1 for PSF fitting photometry
  null_mode    = 1      ; 0: mode, 1: best x%, 2: NSC
  adi_mode     = 0      ; 0: regular adi
  no_gpu       = 1.     ; Damned, gpuFFT does not work with the free version!!!
  keep_ratio   = 0.03
  n_clip       = 128.
  n_bin        = 1
  n_trans      = 5
  no_center    = 1      ; Don't recenter the photometric frames (centering computed only for tip/tilt)
  no_dark      = 1
  no_flat      = 1
  ;n_frbck      = 0
  ;n_frob       = 2500

  ; Frame selection and derotation parameters
  right_handed = 0
  rin_init     = 0
  n_coadd      = 1
  sig_sl       = 0
  sig_pos      = 0
  x_range      = 0
  y_range      = 0

  ; Perform data reduction
  ; **********************
  IF KEYWORD_SET(reduce) THEN BEGIN
    ; reduction parameters
    ;xc       = [108,108]  ; Position of the lower nod but code finds it fine 
    ;yc       = [65,65]    ; Position of the lower nod but code finds it fine 
    
    ; Pare object name properly
    ;LBTI_FITSUTIL, '131224', [1616,9367], OBJNAME='IO', /RAW_DATA
    ;LBTI_FITSUTIL, '131224', [0,1484], OBJNAME='ep_eri', /RAW_DATA
    
    ; Create PSF file
    ; hdr_tmp = HEADFITS(files[i_f])
    ; wav     = FXPAR(hdr_tmp, 'WAVELENG')
    ; bdwth   = FXPAR(hdr_tmp, 'BANDWIDT')
    ; GET_CNF, cnf, LAMBDA=wav, LMIRCAM=lmircam

    ; Create PSF file
    ; img_psf = FIZEAU_PSF(cnf.tel_diam, wav, cnf.base, BANDWIDTH=bdwth, CEN_OBS=cnf.cen_obs, PIX_SIZE=cnf.pix_size, VERBOSE=verbose)
    ; p=STRPOS(files[i_f],'_IMG.fits')
    ; IF p GT 0 THEN psf_file=STRMID(files[i_f], 0, p) + '_PSF.fits' ELSE MESSAGE, 'Invalid file format.'
    ; MKHDR, hdr_psf, img_psf
    ; MWRFITS, img_psf, psf_file, hdr_psf, /CREATE, /SILENT
    
    ; IO data
    ;lambda_cen   = 8.70D-6  ; 0 to read filter wheels
    ;bandwidth    = 1.22D-6  ; 0 to read filter wheels
    ;data_idx = [1868,9367]
    ;nod_idx  = [[1868,2413],[2414,2500],[3242,4100],[4242,4441],[4442,4941],[4942,5141],[5142,5941],[5942,6141],[6142,7124],[7125,7324],[7325,8324],[8325,8524],[8525,9167],[9168,9367]]
    
    lambda_cen   = 11.1D-6  ; 0 to read filter wheels
    bandwidth    = 2.5D-6  ; 0 to read filter wheels
    data_idx = [2842,3200]
    nod_idx  = [[2842,3041],[3042,3200]]
    
    LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx, BAD_IDX=bad_idx, $   ; Date to be reduced and index of files
              DATA_PATH=data_path, BCKG_PATH=bckg_path, DARK_PATH=dark_path, FLAT_PATH=flat_path, $                                                    ; Path to data directories (superseed "*_idx" files)        
              XCEN=xcen, YCEN=ycen, OVERLAP=overlap, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, $
              LAMBDA_CEN=lambda_cen, BANDWIDTH=bandwidth, $
              LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $                          ; Optional input keywords superseeding keywords in the header
              MASTERLOG=masterlog, ADI_MODe=adi_mode, BCKG_MODE=bckg_mode, FLX_MODE=flx_mode, FIT_MODE=fit_mode, NULL_MODE=null_mode, OB_MODE=ob_mode, $
              APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, BCK_CEN=bck_cen, $
              KEEP_RATIO=keep_ratio,  N_BIN=n_bin, N_CLIP=n_clip, N_FRBCK=n_frbck, N_FROB=n_frob, N_TRANS=n_trans, OFFSET=offset,$
              NO_CENTER=no_center, NO_DARK=no_dark, NO_FLAT=no_flat, RIN_INIT=rin_init, SIG_POS=sig_pos, SIG_SL=sig_pl, SKIP_CAL=skip_cal, SKIP_ADI=skip_adi, SKIP_FLX=skip_flx, SKIP_NULL=skip_null, $   ; Optional input keywords critical for data reduction 
              LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                    ; Optional input keywords for output information  
  ENDIF
  
END