;+
; NAME: REDUCE_GJ504
PRO REDUCE_GJ504, date, DATA_IDX=data_idx
             
; General running parameters
log_file = 1         ; Enable log files
no_multi = 0         ; Turn on to debug NSC reduction (multi-thread not easy to debug)
no_save  = 0         ; Prevent L1/L2 data saving 
plot     = 1         ; Plot data on screen if larger than 0
verbose  = 1

; Config files
dqtq_idx = [75403,76777] 

; Loop over the config files
LBTI_DRS, '180329', 'nomic.cfg' , $
         BAD_IDX=bad_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, DATA_IDX=data_idx,  FLAT_IDX=flat_idx, NOD_IDX=nod_idx, OB_IDX=ob_idx, $   ; Date to be reduced and index of files
         DATA_PATH=data_path, BCKG_PATH=bckg_path, DARK_PATH=dark_path, FLAT_PATH=flat_path, $                                             ; Path to data directories (superseed "*_idx" files)
         MASTERLOG=masterlog, RENEW=renew, SKIP_ADI=skip_adi, SKIP_FLX=skip_flx, SKIP_NULL=skip_null, SKIP_RED=skip_red, $                 ; Running keywords
         LOG_FILE=log_file, NO_MULTI=no_multi, NO_SAVE=no_save, PLOT=plot, VERBOSE=verbose                                                 ; Optional input keywords for output information   
                                  
END
