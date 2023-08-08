PRO NODRS, date, DATA_IDX=data_idx, BAD_IDX=bad_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, OB_IDX=ob_idx, REDUCE=reduce, CALIBRATE=calibrate, MASTERLOG=masterlog, $
                 SKIP_CAL=skip_cal, SKIP_ADI=skip_adi, SKIP_FLX=skip_flx, SKIP_NULL=skip_null, VERSION=version, $
                  QUICK=quick, FIT=fit, REMOVE_OB=remove_ob, VERBOSE=verbose, PLOT=plot

  ; General running parameters
  plot        = 1         ; Plot data on screen if larger than 1
  verbose     = 3
  log_file    = 1
  no_save     = 0
  no_multi    = 0          
  cfg_file    = 'fizeau.cfg'
  
  ;LBTI_FITSUTIL, '150308', [8245,11044], OBJNAME='Io'
  bad_idx = [9491,9492,12213]
  
  ; Run reduction
  LBTI_DRS, date, cfg_file[i], $
            BAD_IDX=bad_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, DATA_IDX=data_idx,  FLAT_IDX=flat_idx, NOD_IDX=nod_idx, OB_IDX=ob_idx, $   ; Date to be reduced and index of files
            DATA_PATH=data_path, BCKG_PATH=bckg_path, DARK_PATH=dark_path, FLAT_PATH=flat_path, $                                             ; Path to data directories (superseed "*_idx" files)
            MASTERLOG=masterlog, SKIP_ADI=skip_adi, SKIP_CAL=skip_cal, SKIP_FLX=skip_flx, SKIP_NULL=skip_null, $                              ; Running keywords
            LOG_FILE=log_file, NO_MULTI=no_multi, NO_SAVE=no_save, PLOT=plot, VERBOSE=verbose                                                 ; Optional input keywords for output information
END