;+
; NAME: PDFSDV
;
; PURPOSE:
;   Calculate the mode value of an input PDF and corresponding error bar using high-density regions described by Hyndman, R. J. 1996, The American Statistician, 50, 120
;  
; CALLING SEQUENCE:
;   MODSDV, pdf, bins, mode, err_low, err_sup
;
; INPUTS:
;   indata: Input PDF
;
; OUTPUTS:
;   mode, error low, error sup
;
; MODIFICATION HISTORY:
;  Version 1.0, 28-MAY-2016, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;  Version 1.1, 12-FEB-2017, DD: added interpolation
;-

PRO PDFSDV, pdf_in, bins, dist_mod, err_low, err_sup, BAR_H=bar_h
  ON_ERROR,2
   
  ; Interpol over a finer grid
  nbins = N_ELEMENTS(bins)
  bins_fine = MIN(bins) + (MAX(bins)-MIN(bins))*DINDGEN(10*nbins)/(10*nbins-1)
  pdf_fine  = INTERPOL(pdf_in, bins, bins_fine, /QUADRATIC)
  
  ; Renormalize
  IF TOTAL(pdf_fine) NE 1 THEN scale = FLOAT(TOTAL(pdf_fine)) ELSE scale = 1.
  pdf_fine /= scale
  
  ; 1. Find mode
  max_pdf  = MAX(pdf_fine, idx)
  dist_mod = bins_fine[idx[0]]
  
  ; 2. Find unique values in PDF
  pdf_un = pdf_fine[UNIQ(pdf_fine)]
  pdf_un = REVERSE(pdf_un[SORT(pdf_un)])
  
  ; 2. Compute error bars by computing the cumulative probability density from max_pdf.
  ; The idea here is to use an horizontal bar and lower it from the maximum of the PDF. 
  ; Now use a bar and lower it
  n_stp = N_ELEMENTS(pdf_un)
  FOR i=0, n_stp-1 DO BEGIN
    bar_h = pdf_un[i]
    idx = WHERE(pdf_fine GT bar_h, n_ok)         ; where the PDF is higher than BAR (don't work if double peaked PDF) 
    IF TOTAL(pdf_fine[idx]) GE 0.6827 THEN i=n_stp
  ENDFOR
  
  ; Scale BAR 
  bar_h *= (scale/TOTAL(pdf_in))
  
  ; Normalize cumulative probability and find 68.2% value (use closest value superior to 68% to handle double-peaked distributions)
  binsize = 0.5*(bins[1]-bins[0]) 
  err_low = (dist_mod-MIN(bins_fine[idx])) > 0.5*binsize
  err_sup = (MAX(bins_fine[idx])-dist_mod) > 0.5*binsize
END

; TEST harnerss
PRO TEST_PDFSDV
  ; Create skewed distribution
  npt = 1D6
  d1 = 10*RANDOMN(seed,0.5*npt)
  d2 = RANDOMN(seed2, npt) - 5
  d  = d1;[d1, d2[WHERE(d2 LT 0)]]
  h = HISTOGRAM(d, NBINS=SQRT(npt), LOCATIONS=bins)
  h /= TOTAL(h)
  bins += 0.5*(bins[1]-bins[0])
  PDFSDV, h, bins, mode, err_low, err_sup, BAR_H=bar_o
  PRINT, mode, err_low, err_sup, bar_o
  PLOT, bins, h, PSYM=10
  OPLOT, -50 + INDGEN(100), bar_o + INTARR(100), LINESTYLE=1
END

; OLD APPROACH 
;; To avoid problems with double-peaked PDF, we remove first the part of the PDF that can crate more than two intersections with the bar
;; Extract bins below and above the mode
;idx_low = REVERSE(WHERE(bins LT dist_mod, n_low))
;idx_sup = WHERE(bins GT dist_mod, n_sup)
;
;; Keep only part where the PDF decreases (starting from the mode)
;tmp    = max_pdf
;idx_ok = idx[0]
;FOR i = 0, n_low-1 DO BEGIN
;  IF pdf_in[idx_low[i]] LT tmp THEN BEGIN
;    tmp = pdf_in[idx_low[i]]
;    idx_ok = [idx_ok, idx_low[i]]
;  ENDIF
;ENDFOR
;
;tmp    = max_pdf
;FOR i = 0, n_sup-1 DO BEGIN
;  IF pdf_in[idx_sup[i]] LT tmp THEN BEGIN
;    tmp = pdf_in[idx_sup[i]]
;    idx_ok = [idx_ok, idx_sup[i]]
;  ENDIF
;ENDFOR
;idx_srt = idx_ok[SORT(idx_ok)]
;pdf_tmp = pdf_in[idx_srt]
;bins_t  = bins[idx_srt]
;
;; Now use a bar and lower it
;n_stp   = 1D3
;cum_pr  = FLTARR(n_stp)
;lim_low = cum_pr
;lim_sup = cum_pr
;bar_h   = (1-DINDGEN(n_stp)/(n_stp-1))*max_pdf   ; bar height (starting from maximum)
;FOR i=0, n_stp-1 DO BEGIN
;  idx = WHERE(pdf_tmp GT bar_h[i], n_ok)         ; where the PDF is higher than BAR (don't work if double peaked PDF)
;  IF n_ok GT 0 THEN BEGIN
;    lim_low[i] = MIN(bins_t[idx])
;    lim_sup[i] = MAX(bins_t[idx])
;    idx = WHERE(bins GE lim_low[i] AND bins LE lim_sup[i])
;    cum_pr[i] = TOTAL(pdf_in[idx])
;  ENDIF
;ENDFOR
;
;; Normalize cumulative probability and find 68.2% value (use closest value superior to 68% to handle double-peaked distributions)
;idx68   = (WHERE(cum_pr GT 0.6827))[0]
;bar_o   = bar_h[idx68]
;err_low = (dist_mod-lim_low[idx68]) > 0.5*binsize
;err_sup = (lim_sup[idx68]-dist_mod) > 0.5*binsize
