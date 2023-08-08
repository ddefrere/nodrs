;+
; NAME: COMPUTE_PSD
; 
; PURPOSE:
;   This procedure compute the FFT, PSD, and cumulative RMS of the input data.
;
; MANDATORY INPUTS:
;   time       : Input data coordinates (must be the same length as data)
;   data       : Data to be Fourier transformed 
;
; OPTIONAL INPUT KEYWORDS
;   FRQ_MAX    : Set this keyword to the maximum frequency of the output FFT (the FFT will be extraqpolated if necessary)
;   REGRID     : Set this keyword to regrid the data on uniform data coordinates
;   WINDOW     : Set this keyword to enable the Hanning windowing
;      
; OUTPUTS:
;   FREQ      : The frequencies at which the FFT is computed
;   FRQ_EXT   : The extrapolated frequency range
;   FFT       : The computed FFT
;   PSD       : The corresponding PSD
;   CUM_RMS   : The corresponding cumulative RMS
;   SLOPE     : The slope of the PSD at high frequencies
;
; KEYWORD:
;   FRQ_LIM     : Maximum frequency
;   REGRID      : If set, regrid by interpolation the input data uniformly
;
; MODIFICATION HISTORY:
;   Version 1.0,  31-OCT-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  21-JUN-2014, DD: improved PSD normalization
;   Version 1.2,  02-SEP-2014, DD: added keyword REGRID
;   Version 1.3,  02-OCT-2014, DD: corrected bug in frequency computation!
;   Version 1.4,  02-OCT-2014, DD: added keyword FRQ_LIM
;   Version 1.5,  02-DEC-2015, DD: cleanup code
;   Version 1.6,  12-JAN-2016, DD: restored the WINDOW keyword which was lost at some point
;   Version 1.7,  16-FEB-2016, DD: minor bug corrected
;   Version 1.8,  16-FEB-2016, DD: now handle duplicate entries

PRO COMPUTE_PSD, time, data, FREQ=freq, FRQ_EXT=frq_ext, FRQ_LIM=frq_lim, FFT=fft1, PSD=psd1, CUM_RMS=cum_rms, SLOPE=coeff, $
                 FRQ_MAX=frq_max, REGRID=regrid, WINDOW=window, INFO=info

; Keyword consistency check
IF NOT KEYWORD_SET(INFO)    THEN info = 0

; Determine the number of data points
n_data = N_ELEMENTS(data)

; Remove duplicates
time0 = time
data0 = data

; Remove duplicates
idx_dup = WHERE(time EQ SHIFT(time, 1), n_dup, COMPLEMENT=idx_ok)
IF n_dup GT 0 THEN time = time[idx_ok]
IF n_dup GT 0 THEN data = data[idx_ok]

; The time series might not be evenly-spaced, determine the mean period and maximum frequency.
ptmp    = time-SHIFT(time, 1)
period  = MEAN(ptmp[1:*])
IF NOT KEYWORD_SET(FRQ_MAX) THEN frq_max = 1D/period

; Limiting frequency
IF NOT KEYWORD_SET(FRQ_LIM) THEN frq_lim = frq_max

; Resample the time on an evenly spaced grid and compute the frequency coordinates (see FFT manual) 
x = (FINDGEN((n_data - 1)/2) + 1)
IF n_data MOD 2 THEN freq = [0.0, x, -(n_data/2 + 1) + x]/(n_data*period) $
                ELSE freq = [0.0, x, n_data/2, -n_data/2 + x]/(n_data*period)

; Interpolate OPD over the new time grid 
data_in = data-MEAN(data)     ; First remove mean value
time_in = MIN(time) + (MAX(time)-MIN(time))*DINDGEN(n_data)/(n_data-1)
IF KEYWORD_SET(REGRID) THEN data_in = INTERPOL(data, time, time_in, /SPLINE)

; Apply windowing
IF KEYWORD_SET(WINDOW) THEN BEGIN
  filter = HANNING(n_data, /DOUBLE)
  data_in *= filter
ENDIF ELSE filter = 1

; Compute FFT
fft = FFT(data_in, -1)/TOTAL(filter)*n_data  ; compensate for the loss of energy due to windowing
psd = ABS(fft)^2/(frq_max/n_data)            ; see IDL definition of PSD in Ferree PhD 1999

; Compute one-sided FFT and PSD  (normalized properly)
frq_pos = LINDGEN(FLOOR(0.5*n_data)+1)
freq    = freq[frq_pos]
fft1    = 2.*fft[frq_pos]    ; Amplitude is 2*FFT
psd1    = 2.*psd[frq_pos] & psd1[0] = psd[0] & psd1[FLOOR(0.5*n_data)] = psd[FLOOR(0.5*n_data)]

; Compute slope of PSD
;psd_bin  = SMOOTH(psd, FLOOR(n_data/10)>1, /EDGE_TRUNCATE)
coeff   = ROBUST_LINEFIT(ALOG10(freq[1:*]), ALOG10(psd1[1:*])) ; avoid 0 frequency with log!

; Now compute the cumulative RMS from freq_lim.
; 1. define frequency domain for extrapolation
IF frq_lim GT frq_max THEN BEGIN
  n_frc   = FLOOR(frq_lim/(frq_max/n_data))
  frq_ext = frq_lim/n_frc*(DINDGEN(n_frc)+1)
  ; Extrapolate up to frq_lim
  psd_ext = 10^(coeff[1]*ALOG10(frq_ext[1:*])+coeff[0])        ; avoid 0 frequency with log!
  ; Parse the missing PSD up to freq limit (+ 0 frequency term)
  idx_ext  = WHERE(frq_ext GT frq_max, n_ext)
  IF n_ext GT 0 THEN psd_tmp = [psd1,psd_ext[idx_ext]] ELSE psd_tmp = psd1
ENDIF ELSE BEGIN
  frq_ext = freq
  n_frc   = N_ELEMENTS(freq)
  psd_tmp = psd1
ENDELSE

; Compute cumulative RMS
AVGSDV, data, avg_data, rms_data
cum_rms = REVERSE(SQRT((MAX(frq_ext)/n_frc)*TOTAL(REVERSE(psd_tmp), /CUMULATIVE)))
cum_rms *= rms_data/cum_rms[0]

; Consitency checks.
; 1. Energy consevration (Parseval theorem within 5%)
e1 = TOTAL(data_in^2)
e2 = TOTAL(psd1);*frq_max/n_data
IF info GT 0 THEN IF 2.*ABS(e1-e2)/(e1+e2) GT 0.05 THEN PRINT, 'Parseval identity  failed. There is a problem with the FFT computation.'; ELSE IF info GT 0 THEN PRINT, 'Parseval identity OK'

; 2. Cumulative RMS should be close to measured rms (within 5%)
;AVGSDV, data_in, avg_data, rms_data
;IF 2.*ABS(rms_data-cum_rms[0])/(rms_data+cum_rms[0]) GT 0.05 THEN MESSAGE, 'The cumulative RMS does not corresponds to the RMS of the time series' ELSE IF info GT 0 THEN PRINT, 'Consistency check on cumulative RMS OK'

; 
time = time0
data = data0
END

