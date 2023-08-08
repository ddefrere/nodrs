PRO REDUCE_KAPAND, SKIP_RED=skip_red, SKIP_FLX=skip_flx, MASTERLOG=masterlog

; Note: data reduced with no sign of the planet but the parralactic angle given by the telescope seems fishy
date_obs = '131024'

; General running parameters
plot      = 1         ; Plot data on screen if larger than 1
verbose   = 3
log_file  = 1
no_save   = 0
cfg_file  = 'nomic.cfg';'grism.cfg';'leech.cfg'
masterlog = 0
renew     = 1

; Perform data reduction
;flat_idx = [19904, 20103]
;dark_idx = [20104, 20203]

; Pre-crop the frames
; Add nodding (DX only)
;n_nod = (829+1)/5
;FOR in = 0, n_nod-1 DO LBTI_FITSUTIL, date_obs, [in*5,(in+1)*5-1], LBT_RYOS=4.5*(in MOD 2)
;LBTI_FITSUTIL, date_obs, [2820,3819], DATATYPE=1
;LBTI_FITSUTIL, date_obs, [3820,4819], DATATYPE=2
;+ PRE_CROP  384,384,511,895  in config file

data_idx = [25,829]  ; Not good before that
LBTI_DRS, date_obs, cfg_file, $
          DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx, BAD_IDX=bad_idx, $   ; Date to be reduced and index of files
          DATA_PATH=data_path, BCKG_PATH=bckg_path, DARK_PATH=dark_path, FLAT_PATH=flat_path, $                             ; Path to data directories (superseed "*_idx" files)
          MASTERLOG=masterlog, RENEW=renew, SKIP_ADI=skip_adi, SKIP_RED=skip_red, SKIP_FLX=skip_flx, SKIP_NULL=skip_null, $ ; Running keywords
          LOG_FILE=log_file, NO_SAVE=no_save, PLOT=plot, VERBOSE=verbose                                                    ; Optional input keywords for output information


END