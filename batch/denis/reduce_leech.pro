;+
; NAME: REDUCE_LEECH
;   
; PURPOSE
;   Main batch routine to reduce LBTI LEECH data (see LBTI_DRS.pro for more information)
;
; MANDATORY INPUT
;   date          :  String with the UT date to be reduced based on the format 'yymmdd' (e.g., '130524')
; 
; OPTIONAL INPUT KEYWORDS (FILE INDEX)
;   BAD_IDX       :  Vector with the file number of frames to be removed
;   BCKG_IDX      :  Two-element vector with the lower and the upper file numbers of the background files (e.g., [10,19], used when BCKG_MODE = 4).
;   DARK_IDX      :  Two-element vector with the lower and the upper file numbers of the dark files (e.g., [0,9])
;   DATA_IDX      :  Two-element vector with the lower and the upper file numbers of the data files (e.g., [20,999])
;   FLAT_IDX      :  Two-element vector with the lower and the upper file numbers of the flat files (e.g., [10,19])
;   NOD_IDX       :  2x(number of nods) array with the lower and the upper file numbers of each nod position (e.g., [[120,239],[240,299],[300,399]])
;   
; OPTIONAL INPUT KEYWORDS (WHAT TO DO)
;   AGPM          :  Set this keyword if AGPM data
;   MASTERLOG     :  Set this keyword to force the creation of new masterlog file
;   RENEW         :  Set this keyword to force the code to create a new file (e.g., darks, flats, l0_fits, l1_fits, ...).
;                 :  Otherwise, the code doesn't process the data if the appropriate intermediate file already exists (with the same reduction parameters)
;   SKIP_ADI      :  Set this keyword to skip ADI processing (superseed the value in the config file)
;   SKIP_FLX      :  Set this keyword to skip flux computation (files restored from disk, superseed the value in the config file)
;   SKIP_RED      :  Set this keyword to skip image reduction (files restored from disk, superseed the value in the config file)
;       
; LAST MODIFICATION
;   17-APR-2015, Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   27-MAR-2016, DD: added keyword RENEW

PRO REDUCE_LEECH, date, BAD_IDX=bad_idx, BCKG_IDX=bckg_idx,  DARK_IDX=dark_idx, DATA_IDX=data_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx, $     ; Date to reduce and file index 
                  AGPM=agpm, MASTERLOG=masterlog, RENEW=renew, SKIP_ADI=skip_adi, SKIP_FLX=skip_flx, SKIP_RED=skip_red                       ; Control what to do
             
; General running parameters
log_file = 1         ; Enable log files
no_multi = 0         ; Turn on to debug NSC reduction (multi-thread not easy to debug)
no_save  = 0         ; Prevent L1/L2 data saving 
plot     = 1         ; Plot data on screen if larger than 1
verbose  = 3         ; Print onfo to screen if >0

; Config files
IF KEYWORD_SET(agpm) THEN cfg_file = 'agpm.cfg' ELSE cfg_file = 'leech.cfg'
  
; Perform data reduction   
LBTI_DRS, date, cfg_file, $
          BAD_IDX=bad_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, DATA_IDX=data_idx,  FLAT_IDX=flat_idx, NOD_IDX=nod_idx, OB_IDX=ob_idx, $   ; Date to be reduced and index of files
          DATA_PATH=data_path, BCKG_PATH=bckg_path, DARK_PATH=dark_path, FLAT_PATH=flat_path, $                                             ; Path to data directories (superseed "*_idx" files)
          MASTERLOG=masterlog, RENEW=renew, SKIP_ADI=skip_adi, SKIP_FLX=skip_flx, SKIP_NULL=skip_null, SKIP_RED=skip_red, $                 ; Running keywords
          LOG_FILE=log_file, NO_MULTI=no_multi, NO_SAVE=no_save, PLOT=plot, VERBOSE=verbose                                                 ; Optional input keywords for output information   
                                  
END
