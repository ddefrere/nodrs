;+
; NAME: AVGERR
; 
; PURPOSE:
;   Basic routine to compute the average and error bar of a sample of measurements 
;
; INPUTS:
;   data          :  measurements
;   data_err      :  corresponding error bars
;
; KEYWORDS
;   KAPPA         :  On input, number of sigma for the sigma clipping
;   IDX_OUT       :  On output, data points that have been filter out 
; 
; MODIFICATION HISTORY:
;   Version 1.0,  03-AUG-2017, by Denis Defrère, University of Arizona, denis@lbti.org

PRO AVGERR, data0, data_err0, avg_unw, avg_wei, err_stat, err_sca_unw, err_sca_wei, exc_rms, KAPPA=kappa, IDX_OUT=idx_out, N_OUT=n_out, INFO=info

  ; Redefine data to preserve the input
  data = data0
  data_err = data_err0
  
  ; KEYWORD check
  IF NOT KEYWORD_SET(KAPPA) THEN kappa = 3
  
  ; 1. Remove outliers based on values
  AVGSDV, data, avg, rms, KAPPA=kappa
  idx_in1  = WHERE(ABS(data-avg) LE kappa*rms, COMPLEMENT=idx_out1, NCOMPLEMENT=n_out)
  data     = data[idx_in1]
  data_err = data_err[idx_in1] 
  IF KEYWORD_SET(INFO) THEN PRINT, 'Number of elements rejected based on value :', n_out
  
  ; 2. Remove outliers based on error bar
  AVGSDV, data_err, avg, rms, KAPPA=kappa
  idx_in2  = WHERE(ABS(data_err-avg) LE kappa*rms, COMPLEMENT=idx_out2, NCOMPLEMENT=n_out) & IF n_out GT 0 THEN idx_out1 = [idx_out1, idx_in1[idx_out2]]
  data     = data[idx_in2]
  data_err = data_err[idx_in2] 
  IF KEYWORD_SET(INFO) THEN PRINT, 'Number of elements rejected based on error bar :', n_out
  
  ; 3. Also remove point which are more than kappa sigma from the unweighted mean and median (using their own error!, if at least 5 points)
  ; Iterate once to compute the weighted mean
  n_data   = N_ELEMENTS(data)
  IF n_data GE 5 THEN BEGIN
    data_med = MEDIAN(data)
    AVGSDV, data, avg_unw, rms, err_sca_unw, KAPPA=kappa
    idx_in3  = WHERE(ABS(data-avg_unw) LE kappa*data_err OR ABS(data-data_med) LE kappa*data_err OR ABS(data-data_med) LE rms, n_data, COMPLEMENT=idx_out3, NCOMPLEMENT=n_out) 
    data_med = MEDIAN(data)
    AVGSDV, data[idx_in3], avg_unw, rms, err_sca_unw, KAPPA=kappa
    idx_in3  = WHERE(ABS(data-avg_unw) LE kappa*data_err OR ABS(data-data_med) LE kappa*data_err OR ABS(data-data_med) LE rms, n_data, COMPLEMENT=idx_out3, NCOMPLEMENT=n_out) & IF n_out GT 0 THEN idx_out1 = [idx_out1, idx_in1[idx_in2[idx_out3]]]
    IF KEYWORD_SET(INFO) THEN PRINT, 'Number of elements rejected based on distance to mean :', n_out
    data     = data[idx_in3]
    data_err = data_err[idx_in3]
  ENDIF  
  
  ; 4. Compute statistical error
  err_stat = SQRT(1./TOTAL(1./data_err^2))
  ; 5. Compute weighted mean and weighted error (and same for unweighted but using the same data point this time)
  AVGSDV, data, avg_unw, rms, err_sca_unw, KAPPA=kappa
  AVGSDV, data, avg_wei, rms_sca, err_sca_wei, KAPPA=kappa, WEIGHT=1./data_err^2
  ; 6. Compute excess variance for this pointing
  IF n_data GT 3 THEN exc_rms = SQRT(EXCESS_VARIANCE(data, data_err, MEANX=meanx, MSERR=mserr, VARX=varx, DIMENSION=dimension,  CHATTER=chatter, ERROR=error) > 0 ) $
                 ELSE exc_rms = 0
  ; If only one point, make sure it's 0
  IF n_data LE 1 THEN BEGIN
    err_sca_wei = 0.
    err_sca_unw = 0.
    avg_wei     = data
    avg_unw     = data
  ENDIF
  
  ; Sort output IDX
  idx_out = idx_out1[SORT(idx_out1)]
  idx_out = idx_out[WHERE(idx_out GE 0, n_out)]
END

PRO TEST_AVGERR  
  n_pt     = 100
  kappa    = 3
  avg_in   = 10
  err_in   = 1
  data     = avg_in + RANDOMN(seed, n_pt)
  data_err = err_in + 0.1*err_in*RANDOMN(seed2, n_pt)
  AVGERR, data, data_err, avg_unw, avg_wei, err_stat, err_sca_unw, err_sca_wei, exc_rms, KAPPA=kappa, IDX_OUT=idx_out, /INFO
  PRINT, 'Normal data'
  PRINT, avg_in, avg_unw, avg_wei
  PRINT, err_in/SQRT(n_pt), err_stat, err_sca_unw, err_sca_wei, exc_rms
  PRINT, N_ELEMENTS(idx_out)
  
  ; Produce irregular data for 5% of the frame
  idx_rndm = FIX(n_pt*RANDOMU(seed3, FIX(0.05*n_pt)))
  idx_rndm = idx_rndm[SORT(idx_rndm)]
  ;data[idx_rndm] += 200;FIX(10*err_in*RANDOMU(seed4, FIX(0.05*n_pt)))
  ;data_err[idx_rndm] += 200
  data_err[idx_rndm] = 0
  AVGERR, data, data_err, avg_unw, avg_wei, err_stat, err_sca_unw, err_sca_wei, exc_rms, KAPPA=kappa, IDX_OUT=idx_out, /INFO
  PRINT, ' Irregular data'
  PRINT, avg_in, avg_unw, avg_wei
  PRINT, err_in/SQRT(n_pt), err_stat, err_sca_unw, err_sca_wei, exc_rms
  PRINT, N_ELEMENTS(idx_out), N_ELEMENTS(idx_rndm)
  PRINT, idx_rndm
  PRINT, idx_out
  
END