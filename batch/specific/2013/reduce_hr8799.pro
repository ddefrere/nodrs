PRO REDUCE_HR8799, REDUCE=reduce, MASTERLOG=masterlog, SKIP_RED=skip_red, SKIP_FLX=skip_flx

; Infos
; I measured 3.5 mas of RMS tip/tilt on the PSF data

; General running parameters
plot        = 1         ; Plot data on screen if larger than 1
verbose     = 3
log_file    = 1
no_save     = 0
no_multi    = 0          ; Turn on to debug NSC reduction (multi-thread not easy to debug)
cfg_file    = 'leech.cfg'
date        = '131017'

  ; reduction parameters
  ;xc = [284] & yc = [257+44]  ; Ghost (Nod 1 PSF data)
  ;xc = [285] & yc = [257] 
  ;xc = [541] & yc = [513] ; if full frame
  ;offset = [2.52, -43.39]    ; as measured on the PSF data
  ;bck_irad = 15
  ;bck_orad = 25
 
  ; Pre-crop the frames
  ;LBTI_FITSUTIL, date, [10,19730], CROP=[256,256,512,512], LMIRCAM=lmircam
  
  ; Compute a master dark for Carlos
  ;flat_idx = [19301,19630]
  ;LBTI_DRS, date, LMIRCAM=lmircam, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx, BAD_IDX=bad_idx, $   ; Date to be reduced
  ;          DATA_PATH=data_folder, BCKG_PATH=bckg_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                            ; Path to data directories
  ;          XCEN=xc, YCEN=yc, OVERLAP=overlap, TGT_NAME=tgt_name, FLAG='BPM', OBSTYPE=obstype, $                                                     ; Optional input keywords superseeding keywords in the header
  ;          APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, OFFSET=offset, OB_MODE=ob_mode, NO_GPU=no_gpu, $
  ;          FIT_METHOD=fit_method, NO_CENTER=no_center, SKIP_CAL=skip_cal,  $                                                                        ; Optional input keywords critical for data reduction
  ;          LOG_FILE=log_file, INFO=info, MOVIE=movie, PLOT=plot, NO_SAVE=no_save                                                                    ; Optional input keywords for output information

   ; Dark and flat frames
  flat_idx = [19350,19380]
  dark_idx = [19631,19730]
  
  ; PSF data
;   data_idx = [18560,18959];,18559];[8960,18559]; [10,18559]  ;18555
  ; Science data
  data_idx = [10,18559];,18559];[8960,18559]; [10,18559]  ;18555
  nod_idx  = [[10,359],[360,559],[560,1759],[1760,1959],[1960,3959],[3960,4159],[4160,5159],[5160,5359],[5360,6359],[6360,6559],[6560,7559],[7560,7759],[7760,8759],[8760,8959],[8960,9959],$
              [9960,10159],[10160,11159],[11160,11359],[11360,12359],[12360,12559],[12560,13559],[13560,13759],$
              [13760,14759],[14760,14959],[14960,15959],[15960,16159],[16160,17159],[17160,17359],[17360,18359],[18360,18559]]
  bad_idx  = [36,152,153,154,155,156,157,158,159,160,4162+INDGEN(64),4360+INDGEN(138),4677+INDGEN(56)]
  LBTI_DRS, date, cfg_file, $
      BAD_IDX=bad_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, DATA_IDX=data_idx,  FLAT_IDX=flat_idx, NOD_IDX=nod_idx, OB_IDX=ob_idx, $   ; Date to be reduced and index of files
      DATA_PATH=data_path, BCKG_PATH=bckg_path, DARK_PATH=dark_path, FLAT_PATH=flat_path, $                                             ; Path to data directories (superseed "*_idx" files)
      MASTERLOG=masterlog, SKIP_ADI=skip, SKIP_RED=skip_red, SKIP_FLX=skip_flx, VERBOSE=verbose
END
