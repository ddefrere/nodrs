;+
; NAME: MODE_WEIGHT
; 
; PURPOSE:
;   Estimate the (weighted) mode of the input data series by computing the underlying probability density using the kernel density estimator method.
;   The error bar is computed by bootstrapping and confidence intervals
;
; INPUTS:
;   data_in    :  The 1-D input data series
;
; KEYWORDS
;   MODE_RMS   :  Set this keyword to compute the standard deviation of the mode (by bootstrapping => time consuming)
;   WEIGHT     :  Set this keyword to the weight of each calue of data_in
;   INFO       :  Set this keyword to print info to screen
;   PLOT       :  Set this keyword to plot the results
;
; OUTPUT
;   A two element vector with the mode and corresponding error of the input data
;
; MODIFICATION HISTORY:
;   Version 1.0,  07-APR-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  09-APR-2014, DD: added keyword MODE_RMS
;   Version 1.2,  21-MAY-2014, DD: renamed MODE_WEIGHT to avoid conflicts with IDL routine MODE.pro
;   Version 1.3,  21-MAY-2014, DD: corrected minor bug

FUNCTION MODE_WEIGHT, data_in, MODE_RMS=mode_rms, WEIGHT=weight, INFO=info, PLOT=plot
    
  ; Derive the number of elements
  n_data   = N_ELEMENTS(data_in)
  IF n_data LE 1 THEN MESSAGE, 'Not enough data points'
  
  ; Critical parameters
  n_bin    = ROUND(SQRT(n_data))
  n_bstrp  = 100*n_bin
  n_rho    = 100*n_bin
   
  ; Keyword check
  IF NOT KEYWORD_SET(WEIGHT) THEN weight = REPLICATE(1D, n_data)
  
  ; Estimate the underlying probability density using the kernel density estimator method (Epanechnikov kernel by default).
  data_min   = MIN(data_in) 
  data_max   = MAX(data_in)
  data_range = data_max-data_min
  data_min   = data_min - 0.1*data_range
  data_max   = data_max + 0.1*data_range
  binsize    = (data_max-data_min)/n_rho
  rho_bin    = data_min + (0.5 + DINDGEN(n_rho))*binsize
  rho_data   = KDE(data_in,rho_bin,WEIGHT=weight,/GAUSSIAN)
  
  ; Find maximum of the probability distribution
  idx_max  = WHERE(rho_data EQ MAX(rho_data), n_max)
  IF n_max GT 1 THEN PRINT, 'WARNING: two maxima have been found. Using the smallest one'
  
  ; Make sure idx_max has one element
  idx_max  = idx_max[0]
  
  ; Compute the mode
  data_mode = rho_bin[MIN(idx_max)]
  
  ; Loop over the boostrap realizations to compute error by bootstrapping
  IF KEYWORD_SET(MODE_RMS) THEN BEGIN
    data_bstrp = DBLARR(n_bstrp)
    FOR i_b = 0, n_bstrp-1 DO BEGIN
      ; Draw random null points with repetition
      rndm_idx = FLOOR(RANDOMU(seed, n_data)*n_data)
      data_tmp = data_in[rndm_idx]
      wei_tmp  = weight[rndm_idx]
      ; Compute mode
      rho_tmp  = KDE(data_tmp,rho_bin,WEIGHT=wei_tmp)
      ; Find maximum of the histogram
      idx_max  = WHERE(rho_tmp EQ MAX(rho_tmp), n_max)
      ; Compute weighted average
      data_bstrp[i_b] = rho_bin[MIN(idx_max)]
    ENDFOR
    
    ; Compute the confidence interval (16% and 84% is the 1-sigma definition)
    data_sort = data_bstrp[SORT(data_bstrp)]
    data_low = data_sort[ROUND(0.16*n_bstrp)]
    data_up  = data_sort[ROUND(0.84*n_bstrp)]
    
    ; Be conservative, use the largest error as rms.
    mode_rms = 0.5*(ABS(data_up-data_mode) + ABS(data_low-data_mode))
    ;mode_rms = ABS(data_up-data_mode) > ABS(data_low-data_mode)
  ENDIF ELSE mode_rms = 0
  
  ; Plot if requested
  IF KEYWORD_SET(PLOT) THEN BEGIN
    ; Compute histogram
    data_hist = HISTOGRAM(data_in, NBINS=FIX(SQRT(n_data)), MIN=data_min, MAX=data_max, LOCATIONS=bin_hist)
    bin_hist  = bin_hist + 0.5*(bin_hist[1]-bin_hist[0])
    PLOT, bin_hist, data_hist, PSYM=10 
    LOADCT, 12, /SILENT
    OPLOT, rho_bin, rho_data/TOTAL(rho_data)/binsize*TOTAL(data_hist)*(bin_hist[1]-bin_hist[0]), COLOR=80
    LOADCT, 0, /SILENT
  ENDIF
  
RETURN, [data_mode, mode_rms]
END   
  