;+
; NAME: BIN_NULLDATA
; 
; PURPOSE:
;   Perform null statistics and return a three element vector containing the null, the null error, and the photometric error on the null.
;
; INPUTS:
;   null          :  Input null measurements
;   null_err      :  Corresponding errors
;
; KEYWORDS
;   RATIO         :  Ratio of the input data to bin
;   MODE          :  If set, return the mode of the input data
;   INFO
;   PLOT
;
; OUTPUT
;   The weighted mean of the input null values or best input values if ratio is set
;
; MODIFICATION HISTORY:
;   Version 1.0,  01-OCT-2012, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  10-FEB-2014, DD: added keyword HISTOGRAM
;   Version 1.2,  07-APR-2014, DD: now compute the unbiased error on the null using AVGSDV
;   Version 1.3,  08-APR-2014, DD: added photometric error to output vector
;   Version 1.4,  09-APR-2014, DD: added keyword FLOOR to keep the ratio% around the minimum null value
;   Version 1.5,  10-APR-2014, DD: replaced HISTOGRAM keyword by the MODE keyword
;   Version 1.6,  10-APR-2014, DD: now computing the mode using Kernel density estimation (added call to MODE.pro)
;   Version 1.7,  15-APR-2014, DD: removed keyword FLOOR and implemented error bar computation by bootstrapping
;   Version 1.8,  16-APR-2014, DD: added keyword INFO

FUNCTION BIN_NULLDATA, null, null_err, RATIO=ratio, MODE=mode, INFO=info, PLOT=plot

; Keyword check
IF NOT KEYWORD_SET(INFO) THEN info = 0
IF NOT KEYWORD_SET(PLOT) THEN plot = 0
  
; Number of null
n_data0 = N_ELEMENTS(null)

IF NOT KEYWORD_SET(MODE) THEN BEGIN  
  ; If ratio is set, select the best ratio% frames either around the minimum value (if FLOOR is set) or around the mean value.
  IF KEYWORD_SET(ratio) AND n_data0 GT 1 THEN BEGIN
      ; Number of data point to keep data
      n_data   = FLOOR(ratio*n_data0)
      ; Sort the data by closest to minimum value
      out      = null[SORT(null)]
      out_err  = null_err[SORT(null)]
      ; Select the closest values to minimum or mean value
      null_out = out[INDGEN(n_data)]
      nerr_out = out_err[INDGEN(n_data)] 
  ENDIF ELSE BEGIN
      n_data   = n_data0
      null_out = null
      nerr_out = null_err
  ENDELSE
    
  ; Mean null computed as the weighted-average of 3-sigma clipped nulls
  AVGSDV, null_out, avg_null, rms_err, WEI=(1/nerr_out)^2, KAPPA=5.
  
  ; Compute error on the mean by bootstrapping
  n_bstrp    = 5000
  data_bstrp = DBLARR(n_bstrp)
  FOR i_b = 0, n_bstrp-1 DO BEGIN
    ; Draw random null points with repetition
    rndm_idx = FLOOR(RANDOMU(seed, n_data0)*n_data0)
    null_tmp = null[rndm_idx]
    err_tmp  = null_err[rndm_idx]
    ; Sort the data by closest to minimum value
    out      = null_tmp[SORT(null_tmp)]
    out_err  = err_tmp[SORT(null_tmp)]
    ; Select the closest values to minimum or mean value
    null_out = out[INDGEN(n_data)]
    nerr_out = out_err[INDGEN(n_data)] 
    ; Compute weighted average (and remove obvious outliers)
    AVGSDV, null_out, avg_tmp, rms_tmp, WEI=(1/nerr_out)^2, KAPPA=5.
    data_bstrp[i_b] = avg_tmp
  ENDFOR

  ; Compute the confidence interval (16% and 84% is the 1-sigma definition)
  data_sort = data_bstrp[SORT(data_bstrp)]
  data_low  = data_sort[ROUND(0.16*n_bstrp)]
  data_up   = data_sort[ROUND(0.84*n_bstrp)]

  ; Be conservative, use the largest error as rms.
  rms_null  = ABS(data_up-avg_null) > ABS(data_low-avg_null)
  ;rms_null  = rms_err/SQRT(N_ELEMENTS(null_out))
  
  ; Compute weighted photometric error on the mean
  rms_phot = 1./SQRT(TOTAL(1./nerr_out^2))
  
  ; Print info to screen
  IF info GT 1 THEN PRINT, 'Mean null, mean bootstrapped null, and corresponding bias :', avg_null, MEAN(data_bstrp), MEAN(data_bstrp)-avg_null
ENDIF ELSE BEGIN
  ; Compute mode
  weight   = 1/null_err^2
  tmp      = MODE_WEIGHT(null, WEIGHT=weight, PLOT=plot, /MODE_RMS)
  avg_null = tmp[0]
  rms_null = tmp[1]
  rms_phot = 0; 1./SQRT(TOTAL(1./null_err[idx_hist]^2))
ENDELSE

RETURN, [avg_null, rms_null, rms_phot]
END

; ---------
PRO TEST_BINNULLDATA
  ; Read data
  ;data = MRDFITS('/Volumes/nodrs/nodrs/results/nomic/l1_fits/2014-03-17/2014-03-17_ID9_SCI_alf_Lyr_DIT-27ms_11um_NULL.fits', 1, header)
  ;tmp  = MRDFITS('/Volumes/nodrs/nodrs/results/nomic/l1_fits/2014-03-17/2014-03-17_ID9_SCI_alf_Lyr_DIT-27ms_11um_NULL.fits', 0, header)
  data = MRDFITS('/Volumes/nodrs/nodrs/results/nomic/l1_fits/2014-02-12/2014-02-12_ID31_SCI_eta_crv_DIT-85ms_11um_NULL.fits', 1, header)
  tmp  = MRDFITS('/Volumes/nodrs/nodrs/results/nomic/l1_fits/2014-02-12/2014-02-12_ID31_SCI_eta_crv_DIT-85ms_11um_NULL.fits', 0, header)
  flx_tot = FXPAR(header, 'PHOT_AVG')
  
  ; Extract data
  flx_null = data.flx_tot/flx_tot
  flx_err  = data.flx_err/flx_tot
  
  ; Ratios
  ratio = (4*DINDGEN(25)+1)*0.01;[0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.75, 0.90, 1.00]
  n_r   = N_ELEMENTS(ratio)
  
  ; Compute null for various ratio
  null     = DBLARR(n_r)
  rms_null = null
  rms_phot = null
  FOR i_r = 0, n_r-1 DO BEGIN
    tmp           = BIN_NULLDATA(flx_null, flx_err, RATIO=ratio[i_r], INFO=2)
    null[i_r]     = tmp[0]
    rms_null[i_r] = tmp[1]
    rms_phot[i_r] = tmp[2]
    PRINT, ratio[i_r], null[i_r], rms_null[i_r], rms_phot[i_r]
  ENDFOR
  
  null_mode = BIN_NULLDATA(flx_null, flx_err, /MODE, /PLOT)
  PRINT, 'Mode', null_mode
  
  ; Plot results
  result_path = '/Volumes/nodrs/nodrs/results/nomic/'
  plotname    = result_path + 'null_analysis.eps'
  PREP_PS & DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=17.8, YSIZE=12.7
  LOADCT, 12, /SILENT
  xrange = [-0.1, 1.1]
  yrange = [0.01,0.2]
  PLOT, [0], [0], /YLOG, XTITLE='Ratio of preserved frames', YTITLE='Measured null', XSTYLE=1, YSTYLE=1, XRANGE=xrange, YRANGE=yrange, XTHICK=4, YTHICK=4, CHARTHICK=2.5, CHARSIZE=1.0
  USERSYM, [1,1,-1,-1,1]*0.2, [1,-1,-1,1,1]*0.2, THICK=1.5
  OPLOT, ratio, null, PSYM=-8, COLOR=100
  ERRPLOT, ratio, null-rms_null, null+rms_null, COLOR=100
  DEVICE, /CLOSE & END_PS
  
  plotname    = result_path + 'null_analysis2.eps'
  PREP_PS & DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=17.8, YSIZE=12.7
  LOADCT, 12, /SILENT
  xrange = [-0.1, 1.1]
  yrange = [0.9*(MIN(rms_null)<MIN(rms_phot)),1.1*(MAX(rms_null)>MAX(rms_phot))]
  PLOT, [0], [0], /YLOG, XTITLE='Ratio of preserved frames', YTITLE='Measured error on the mean', XSTYLE=1, YSTYLE=1, XRANGE=xrange, YRANGE=yrange, XTHICK=4, YTHICK=4, CHARTHICK=2.5, CHARSIZE=1.0
  OPLOT, ratio, rms_null, COLOR=100, THICK=4
  OPLOT, ratio, rms_phot, COLOR=0, THICK=4
  DEVICE, /CLOSE & END_PS
END