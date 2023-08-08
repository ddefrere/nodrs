;+
; NAME: MODSDV
;
; PURPOSE:
;   Calculate the mode value and corresponding error bar using high-density regions described by Hyndman, R. J. 1996, The American Statistician, 50, 120
;
; CALLING SEQUENCE:
;   MODSDV, data, mode, err_low, err_sup
;
; INPUTS:
;   indata: Input array. May be any type except string.
;
; OUTPUTS:
;   mode, error low, error sup
;
; RESTRICTIONS:
;   If N_ELEMENTS is one, then SDV is zero.
;
; MODIFICATION HISTORY:
;  Version 1.0, 28-MAY-2016, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;-

PRO MODSDV, indata, dist_mod, err_low, err_sup
  ON_ERROR,2
  
  ; Number of bins in the histogram
  n_data = N_ELEMENTS(indata)
  n_bin  = ROUND(SQRT(n_data))
  n_rho  = 100*n_bin  

;  ; Estimate the underlying probability density using the kernel density estimator method (Epanechnikov kernel by default).
;  data_min   = MIN(indata)
;  data_max   = MAX(indata)
;  data_range = data_max-data_min
;  data_min   = data_min - 0.1*data_range
;  data_max   = data_max + 0.1*data_range
;  binsize    = (data_max-data_min)/n_rho
;  rho_bin    = data_min + (0.5 + DINDGEN(n_rho))*binsize
;  rho_data   = KDE(indata, rho_bin, WEIGHT=weight, /GAUSSIAN)
;  
;  ; Find maximum of the probability distribution
;  idx_max  = WHERE(rho_data EQ MAX(rho_data), n_max)
;  IF n_max GT 1 THEN PRINT, 'WARNING: two maxima have been found. Using the smallest one'
;  dist_mod = rho_bin[MIN(idx_max[0])]
  
  ; Old way, with histograms (problematic if binsize too large, which is often the case)
  dist_data = HISTOGRAM(indata, NBINS=n_bin, LOCATIONS=rho_bin)
  bin_size  = rho_bin[1]-rho_bin[0]
  rho_bin   = 0.5*bin_size + rho_bin
  
  ; Compute mode and err bars
  PDFSDV, dist_data, rho_bin, dist_mod, err_low, err_sup
END
