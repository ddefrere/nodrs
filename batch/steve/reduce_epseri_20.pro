;+
; NAME: REDUCE_EPSERI_20
;   
; PURPOSE
;   Main batch routine to reduce LBTI nulling data (see LBTI_DRS.pro for more information)
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
;   OB_IDX        :  Two-element vector with the lower and upper OB ID numbers to process during the null computation
;   
; OPTIONAL INPUT KEYWORDS (WHAT TO DO)
;   CALIB_IMG     :  Set this keyword to turn on the null calibration using the L1 images (experimental) 
;   CALIB_NULL    :  Set this keyword to turn on the null calibration
;   FIT           :  Set this keyword to turn on the null fit (that also gives the final source null)
;   L1VERSION     :  Set this keyword to the version number of the L1 summary file to calibrate (last one if not set)
;   L2BOOTSPTRAP  :  Set this keyword to boostrap the L2 data 
;   LIVE          :  Set this keyword for real-time data reduction
;   MASTERLOG     :  Set this keyword to force the creation of new masterlog file
;   QUICK         :  Set this keyword to perform a quick reduction (e.g., real-time data reduction)
;   REDUCE        :  Set this keyword to reduce the data (L0 to L2)
;   REMOVE_OB     :  Set this keyword to an array with the number of the OBs to ignore in the null calibration
;   REMOVE_ID     :  Set this keyword to an array with the ID number of the files to ignore in the null calibration
;   RENEW         :  Set this keyword to force the code to create a new file (e.g., darks, flats, l0_fits, l1_fits, ...).
;                 :  Otherwise, the code doesn't process the data if the appropriate intermediate file already exists (with the same reduction parameters)
;   SKIP_ADI      :  Set this keyword to skip ADI processing (superseed the value in the config file)
;   SKIP_FLX      :  Set this keyword to skip flux computation (files restored from disk, superseed the value in the config file)
;   SKIP_NULL     :  Set this keyword to skip null computation (superseed the value in the config file)
;   SKIP_RED      :  Set this keyword to skip image reduction (files restored from disk, superseed the value in the config file)
;   VERBOSE       :  Define the level of information printed to screen:
;                       - 0: completely silent execution
;                       - 1: minimum level of information
;                       - 2: nominal level of information
;                       - 3: debugging use (print for instance the computation time of each step in the data reduction)
;       
; LAST MODIFICATION
;   17-APR-2015, Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   27-MAR-2016, DD: added keyword RENEW
;   14-FEB-2017, DD: added keyword VERBOSE
;   05-APR-2017, DD: added keyword LIVE
;   05-SEP-2018, DD: added keyword REMOVE_ID

PRO REDUCE_EPSERI_20, date, BAD_IDX=bad_idx, BCKG_IDX=bckg_idx,  DARK_IDX=dark_idx, DATA_IDX=data_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx, OB_IDX=ob_idx, $          ; Date to reduce and file index 
                  CALIB_IMG=calib_img, CALIB_NULL=calib_null, FIT=fit, L1VERSION=l1version, L2BOOTSTRAP=l2bootstrap, LIVE=live, MASTERLOG=masterlog, QUICK=quick, $ 
                  REDUCE=reduce, REMOVE_OB=remove_ob, REMOVE_ID=remove_id, RENEW=renew, SKIP_ADI=skip_adi, SKIP_FLX=skip_flx, SKIP_NULL=skip_null, SKIP_RED=skip_red, VERBOSE=verbose ; Control what to do
             
; General running parameters
log_file = 1         ; Enable log files
no_multi = 0         ; Turn on to debug NSC reduction (multi-thread not easy to debug)
no_save  = 0         ; Prevent L1/L2 data saving 
plot     = 0         ; Plot data on screen if larger than 0
verbose  = 1

; Set up real-time data reduction
IF KEYWORD_SET(LIVE) THEN BEGIN
  quick      = 1
  reduce     = 1
  calib_null = 1
  fit        = 1
ENDIF

; Config files
IF NOT KEYWORD_SET(QUICK) THEN cfg_file = 'epseri20.cfg' ELSE cfg_file = 'quicknull.cfg'
;cfg_file = ['hosts_a8.cfg','hosts_a30.cfg'];,'hosts_a36.cfg']

; Run LIVE loop
run_live:

; Loop over the config files
FOR i=0, N_ELEMENTS(cfg_file)-1 DO BEGIN
  
  ; --- Perform data reduction
  IF KEYWORD_SET(reduce) THEN BEGIN      
      LBTI_DRS, date, cfg_file[i], $
                BAD_IDX=bad_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, DATA_IDX=data_idx,  FLAT_IDX=flat_idx, NOD_IDX=nod_idx, OB_IDX=ob_idx, $   ; Date to be reduced and index of files
                DATA_PATH=data_path, BCKG_PATH=bckg_path, DARK_PATH=dark_path, FLAT_PATH=flat_path, $                                             ; Path to data directories (superseed "*_idx" files)
                MASTERLOG=masterlog, RENEW=renew, SKIP_ADI=skip_adi, SKIP_FLX=skip_flx, SKIP_NULL=skip_null, SKIP_RED=skip_red, $                 ; Running keywords
                LOG_FILE=log_file, NO_MULTI=no_multi, NO_SAVE=no_save, PLOT=plot, VERBOSE=verbose                                                 ; Optional input keywords for output information   
  ENDIF
  
  ; --- Data calibration          
  IF KEYWORD_SET(calib_null) THEN BEGIN
    NULL_CALIB, date, cfg_file[i], REMOVE_ID=remove_id, REMOVE_OB=remove_ob, PLOT=plot, VERSION=l1version, /NO_INSET, /RUNBIAS             ; run calibration on backgorund bias 
    NULL_CALIB, date, cfg_file[i], REMOVE_ID=remove_id, REMOVE_OB=remove_ob, PLOT=plot, VERSION=l1version, /NO_INSET, /CALPOB              ; L2 fits will be per OB
    NULL_CALIB, date, cfg_file[i], REMOVE_ID=remove_id, REMOVE_OB=remove_ob, PLOT=plot, VERSION=l1version, /LOG_FILE, /NO_INSET            ; L2 fits will be per pointing
  ENDIF
  
  ; --- Calibration using the images
  IF KEYWORD_SET(calib_img) THEN NULL_CALIBIMG, date, cfg_file[i], PLOT=plot, /LOG_FILE
  
  ; --- Null bootstrap
  ; IF KEYWORD_SET(l2bootstrap) THEN NULL_BOOTSTRAP, date, cfg_file, LOG_FILE=log_file, REMOVE_OB=remove_ob, PLOT=plot, VERSION=version
  
  ; --- Null fit 
  IF KEYWORD_SET(fit) THEN BEGIN
    NULL_FIT, date, cfg_file[i], PLOT=plot, /PS, L1TYPE='BIAS'
    NULL_FIT, date, cfg_file[i], PLOT=plot, /PS, /LOG_FILE;, L1TYPE='PT' 
  ENDIF
ENDFOR        

; Run LIVE loop
IF KEYWORD_SET(LIVE) THEN GOTO, run_live
                                  
END
