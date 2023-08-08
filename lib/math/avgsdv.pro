;+
; NAME: AVGSDV
;
; PURPOSE:
;   Calculate the average value and the standard deviation of an array,
;   or calculate the average value and standard deviation over one dimension of ;   an array as a function of all the other dimensions.
;
; CALLING SEQUENCE:
;   AVGSDV, data, average, spread
;
; INPUTS:
;   data:
;     Input array. May be any type except string.
;
; OUTPUTS:
;   average, spread, spread on the mean:
;
; KEYWORDS:
;  CLIP:    If set to a value N, clip the N extremes from the data
;  DOUBLE:  Set this keyword for double precision calculations
;  KAPPA:   If set to a value K, remove all data points more than K standard deviations away
;           from the mean and recompute average and standard deviation
;  RUNBIN:  If set to a value M, with M < N_ELEMENTS(data), compute the average and standard deviation in bins of size M
;  WEI:     Can be used to provide weights to the data points.
;
; PROCEDURE:
;   For the computation of the standard deviation the corrected two-pass
;   algorithm is applied to minimize round-off errors: see section 14.1 of
;   Numerical Recipes in Fortran (second edition), publ. by Cambridge Univ. Pr.
;
; RESTRICTIONS:
;   If N_ELEMENTS is one, then SDV is zero.
;
; MODIFICATION HISTORY:
;   Based on the routine AVG by William Thompson, Applied Research Corporation
;   and the Numerical Recipes (2nd edition) routine MOMENTS.
;   Created 06-MAR-1996, by Roland den Hartog, ESA/ESTEC SCI-ST division, NL, rdhartog@rssd.esa.int
;   Version 11-JUN-1997, RdH: weighting implemented
;   Version 15-JUN-2004, RdH: running bins and kappa clipping implemented
;   Version 15-OCT-2014, DD:  added keyword DOUBLE
;   Version 15-MAR-2015, DD:  revised way to compute weighted standard deviation (now unbiased)
;   Version 07-AUG-2015, DD:  now output the standard error on the mean (weighted if specified)
;   Version 18-MAY-2017, DD:  corrected major bug in use of KAPPA!!!!
;-
PRO AVGSDV, indata, avg, sdv, sdv_mean, CLIP=clip, DOUBLE=double, KAPPA=kappa, RUNBIN=runbin, WEIGHT=wei
ON_ERROR,2

data=indata

; Clip extremes
IF KEYWORD_SET(clip) THEN BEGIN
  cnt=clip
  WHILE cnt GT 0 DO BEGIN
    n=N_ELEMENTS(data) & tmp=LONARR(n)
    m=MAX(data,i) & tmp(i)=1.
    m=MIN(data,i) & tmp(i)=1.
    w=WHERE(tmp LT 1.)
    data=data[w]
    cnt=cnt-1L
  ENDWHILE
ENDIF

; Compute data weighting
nd=N_ELEMENTS(data)
IF NOT KEYWORD_SET(wei) THEN wei=REPLICATE(1.,nd)
w=wei*DOUBLE(nd)/TOTAL(wei)

; Compute average and standard deviation
IF KEYWORD_SET(runbin) THEN BEGIN
  avg=data-data & sdv=avg & sdv_mean=avg
  num=(runbin < nd) > 1L
  nbin=nd-num+1L
  npad=LONG(0.5*num)
  FOR i=0L, nbin-1L DO BEGIN
    db=data[i:i+num-1L]
    wb=w[i:i+num-1L]

    av=TOTAL(db*wb) / num
    IF num GT 1 THEN BEGIN
      tmp=(db-av)*wb
      sd=SQRT((TOTAL(tmp^2)-TOTAL(tmp)^2/num)/(num-1D0))
    ENDIF ELSE sd=0

    ; Apply kappa clipping
    IF KEYWORD_SET(kappa) AND sd GT 0 THEN BEGIN
      k=WHERE(ABS(db-av) LE kappa*sd, n)
      IF n GT 1 THEN BEGIN
        n=DOUBLE(n)
        av=TOTAL(db[k]*wb[k]) / n
        tmp=(db[k]-av)*wb[k]
        sd=SQRT((TOTAL(tmp^2)-TOTAL(tmp)^2/n)/(n-1D0))
      ENDIF
    ENDIF

    j=i+npad
    IF i EQ 0 THEN BEGIN
      avg[0:j]=av & sdv[0:j]=sd
    ENDIF ELSE IF i GE nbin-1L THEN BEGIN
      avg[j:*]=av & sdv[j:*]=sd
    ENDIF ELSE BEGIN
      avg[j]=av & sdv[j]=sd
    ENDELSE
    
    ; Compute error on the mean
    sdv_mean[j] = sdv[j]*SQRT(TOTAL((wb/TOTAL(wb))^2))

  ENDFOR
ENDIF ELSE BEGIN
; When a running bin is not required do it quicker:
  num=nd
  avg=TOTAL(data*w, DOUBLE=double) / num
  IF num GT 1 THEN BEGIN
    ;tmp=(data-avg)*w
    ;sdv=SQRT((TOTAL(tmp^2, DOUBLE=double)-TOTAL(tmp, DOUBLE=double)^2/num)/(num-1D0)) ; old approach
    sdv=SQRT(TOTAL(wei*(data-avg)^2, DOUBLE=double)/((1-1D/num)*TOTAL(wei, DOUBLE=double)))
    sdv_mean=sdv*SQRT(TOTAL((wei/TOTAL(wei))^2))   ; Compute error on the mean
    ;xWbar = avg
    ;wbar = mean(wei)
    ;sdv_mean = SQRT(num/((num-1)*TOTAL(wei)^2)*(TOTAL((wei*data-wbar*xWbar)^2)-2*xWbar*TOTAL((wei-wbar)*(wei*data-wbar*xWbar))+xWbar^2*TOTAL((wei-wbar)^2)))
  ENDIF ELSE sdv=0

  ; Apply kappa clipping
  IF KEYWORD_SET(kappa) AND sdv GT 0 THEN BEGIN
    k=WHERE(ABS(data-avg) LE kappa*sdv, n)
    w=wei[k]*DOUBLE(n)/TOTAL(wei[k])
    IF n GT 1 THEN BEGIN
      n=DOUBLE(n)
      avg=TOTAL(data[k]*w, DOUBLE=double) / n
      ;tmp=(data[k]-avg)*w[k]
      ;sdv=SQRT((TOTAL(tmp^2, DOUBLE=double)-TOTAL(tmp, DOUBLE=double)^2/n)/(n-1D0)) ; old approach
      sdv=SQRT(TOTAL(wei[k]*(data[k]-avg)^2, DOUBLE=double)/((1-1D/n)*TOTAL(wei[k], DOUBLE=double)))
      ;xWbar = avg
      ;wbar = mean(wei[k])
      ;sdv_mean = SQRT(n/((n-1)*TOTAL(wei[k])^2)*(TOTAL((wei[k]*data[k]-wbar*xWbar)^2)-2*xWbar*TOTAL((wei[k]-wbar)*(wei[k]*data[k]-wbar*xWbar))+xWbar^2*TOTAL((wei[k]-wbar)^2)))
      sdv_mean = sdv*SQRT(TOTAL((wei[k]/TOTAL(wei[k]))^2))   ; Compute error on the mean
    ENDIF
  ENDIF  
ENDELSE

END
