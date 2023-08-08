PRO REDUCE_HR8799_DIRECT, MASTERLOG=masterlog, RENEW=renew, SKIP_RED=skip_red, SKIP_FLX=skip_flx

; Infos
; I measured 3.5 mas of RMS tip/tilt on the PSF data

; General running parameters
plot        = 1         ; Plot data on screen if larger than 1
verbose     = 3
log_file    = 1
no_save     = 0
no_multi    = 0          ; Turn on to debug NSC reduction (multi-thread not easy to debug)
cfg_file    = 'leech.cfg'
date        = '131021'

; Dark and flat frames
flat_idx = [17029,17038]
dark_idx = [16689,16988]
  
; Science data
data_idx = [140,7362]
nod_idx  = [140,189] 
n_nod    = 122
FOR i=0, n_nod-1 DO nod_idx = [[nod_idx],[190+i*50,239+i*50]]
n_nod    = 20
FOR i=0, n_nod-1 DO nod_idx = [[nod_idx], [6363+i*50,6412+i*50]]  
;  bad_idx  = [36,152,153,154,155,156,157,158,159,160,4162+INDGEN(64),4360+INDGEN(138),4677+INDGEN(56)]

; PSF data
;data_idx = [7363,7402]
;nod_idx  = [[7363,7382],[7383,7402]]

; Next-star data
data_idx = [7403,7852]
nod_idx  = [7403,7452]
n_nod    = 7
FOR i=0, n_nod-1 DO nod_idx = [[nod_idx],[7453+i*50,7502+i*50]]


LBTI_DRS, date, cfg_file, $
      BAD_IDX=bad_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, DATA_IDX=data_idx,  FLAT_IDX=flat_idx, NOD_IDX=nod_idx, OB_IDX=ob_idx, $   ; Date to be reduced and index of files
      DATA_PATH=data_path, BCKG_PATH=bckg_path, DARK_PATH=dark_path, FLAT_PATH=flat_path, $                                             ; Path to data directories (superseed "*_idx" files)
      MASTERLOG=masterlog, RENEW=renew, SKIP_RED=skip_red, SKIP_FLX=skip_flx, VERBOSE=verbose     
END
