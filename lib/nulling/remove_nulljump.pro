;+
; NAME: REMOVE_NULLJUMP
; 
; PURPOSE:
;   Basic routine to remove fringe jumps.
;
; INPUTS:
;   data          :  1-D vector with the null mezasurements
;   threshold     :  Jump threshold between successive frames. When the difference between two frames is larger than that, a new sub-sequence is created
;   threshold2    :  This parameter is meant to keep all sub-sequences with a mean within 5% of the lowest one
;   min_fr        :  Minimum number of frames required to be a valid sub-sequences
;
; KEYWORDS:
;   IDX_OUT       :  Two-element vector with the lower and upper OB ID numbers to process during the null computation
;
; MODIFICATION HISTORY:
;   Version 1.0,  01-OCT-2016, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu 
;   Version 1.1,  05-AUG-2018, DD: added consistency check 

FUNCTION REMOVE_NULLJUMP, data, threshold, threshold2, min_fr, IDX_OUT=idx_out
  
  ; Constancy check
  n_data = N_ELEMENTS(data)
  IF n_data LT min_fr THEN BEGIN
    idx_out = LINDGEN(n_data)
    RETURN, data
  ENDIF
  
  ; Find jumps
  idx_jump = WHERE(ABS(data-SHIFT(data,1)) GT threshold, n_jump)
  
  ; If jump(s) found, extract sub sequences, compute the mean of each, and then keep only those below threshold2
  IF n_jump GT 0 THEN BEGIN
    ; 0 must be the starting index
    n   = N_ELEMENTS(data)-1
    IF idx_jump[0]        NE 0 THEN idx = [0, idx_jump] ELSE idx = idx_jump
    IF idx_jump[n_jump-1] NE n THEN idx = [idx, n]
    
    ; Compute the mean of each subsequence and mean of sequence with enough data
    n_seq = N_ELEMENTS(idx)-1
    seq   = DBLARR(n_seq)
    nfr   = INTARR(n_seq)
    FOR i=0, n_seq-1 DO BEGIN
      seq[i] = MEDIAN(data[idx[i]:idx[i+1]-1])
      nfr[i] = ABS(idx[i+1]-idx[i]) 
    ENDFOR
    
    ; Find lowest sequences of at least min_fr
    seq_min  = MIN(seq[WHERE(nfr GT min_fr)])
    idx_min  = WHERE(ABS(seq-seq_min) LT threshold2, n_ok)
    data_out = data[idx[idx_min[0]]:idx[idx_min[0]+1]-1] 
    idx_out  = idx[idx_min[0]] + LINDGEN(idx[idx_min[0]+1] - idx[idx_min[0]] + 1)
    FOR i = 1, n_ok-1 DO BEGIN
      data_out = [data_out, data[idx[idx_min[i]]:idx[idx_min[i]+1]]]
      idx_out  = [idx_out, idx[idx_min[i]] + LINDGEN(idx[idx_min[i]+1] - idx[idx_min[i]] + 1)]
    ENDFOR
  ENDIF ELSE BEGIN
    data_out = data
    idx_out  = LINDGEN(N_ELEMENTS(data))
  ENDELSE
  
  ; Apply second filter (because above code doesn't work well for stepping jumps)
  sig_out = 3
  AVGSDV, data_out, avg, rms, rmsm, KAPPA=sig_out
  idx_ok  = WHERE(ABS(data_out-MEDIAN(data_out)) LE sig_out*rms, n_ok)
  idx_out  = idx_out[idx_ok]
  data_out = data_out[idx_ok]

  RETURN, data_out
END

; TEST HARNESS
PRO TEST1
  ;data = MRDFITS('/Users/denis/IDLWorkspace84/Divers/LBTI/UT2016-10-16_ID019_CAL_HD209960_DIT-53ms_11um_NULL.fits', 1, hdr)  ; easy one
  ;data = MRDFITS('/Volumes/hosts_results/nomic/l1_fits/2016-11-14/UT2016-11-14_ID077_CAL_HD49968_DIT-68ms_11um_NULL.fits', 1, hdr) ; crazy one
  
  
  data = MRDFITS('/Volumes/hosts_results/nomic/l1_fits/2014-02-12_APR/UT2014-02-12_ID004_CAL_HD108522_DIT-85ms_11um_NULL.fits', 1, hdr) ; CG only
  pho1 = MRDFITS('/Volumes/hosts_results/nomic/l1_fits/2014-02-12_APR/UT2014-02-12_ID004_CAL_HD108522_DIT-85ms_11um_PHOT1.fits', 1, hdr) ; CG only
  pho2 = MRDFITS('/Volumes/hosts_results/nomic/l1_fits/2014-02-12_APR/UT2014-02-12_ID004_CAL_HD108522_DIT-85ms_11um_PHOT2.fits', 1, hdr) ; CG only
  
  phot_tot = MEAN(pho1.flx_tot[0]) + MEAN(pho2.flx_tot[0]) + 2*SQRT(MEAN(pho1.flx_tot[0])*MEAN(pho2.flx_tot[0]))
 
  null = data.flx_tot[0]/phot_tot
  
  null_ok = REMOVE_NULLJUMP(null, 0.18, 0.05)
  PLOT, null
  LOADCT, 3
  OPLOT, null_ok, COLOR=90
END