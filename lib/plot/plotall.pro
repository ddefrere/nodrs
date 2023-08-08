;+
; NAME: PLOTALL
; 
; PURPOSE:
;   General plotting routine. Plot time series, histogram, FFT, and PSD
;
; INPUTS:
;   time     :  Time stamps corresponding to data (in sec)
;   data     :  Data series
;   error    :  Corresponding error (ignored if set to 0)
;
; KEYWORDS
;   NAME     :  Complete filename (without extension)
;   NO_FFT   :  If set, the FFT and corresponding plots are not computed
;   NO_HISTO :  If set, the histogram and corresponding plots are not computed
;   SCATTER  :  Set this keyword for scatter plot
;   TITLE    :  Title of the figures
;   XLOG     :  Set this keyword to plot the sequence on a log X scale
;   YLOG     :  Set this keyword to plot the sequence on a log Y scale
;   XRANGE   :  Set this keyword to force the XRANGE
;   YRANGE   :  Set this keyword to force the YRANGE
;   XTITLE   :  Title of the x axis
;   YTITLE   :  Title of the y axis
;   EPS      :  Send plots to EPS files rather than display them on screen
;
; MODIFICATION HISTORY:
;   Version 1.0, 16-APR-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 28-SEP-2014, DD: added keyword KERNEL
;   Version 1.2, 28-NOV-2014, DD: now properly ignore the cases where error = 0 and added keyword NO_FFT
;   Version 1.3, 12-FEB-2015, DD: nw works for negative 'time'
;   Version 1.4, 16-FEB-2015, DD: added keyword SCATTER
;   Version 1.5, 05-AUG-2015, DD: added keyword XLOG, YLOG, XRANGE, and YRANGE
;   Version 1.6, 25-AUG-2015, DD: added keyword BOLD
;   Version 1.7, 25-MAR-2016, DD: minor bug correction in file names
;   Version 1.8, 27-MAR-2016, DD: improved axis range computation
;   Version 1.9, 30-MAR-2016, DD: rearranged thye LOADCT commands to be after PLOTXY (which sdoes the SET_PLOT. LOADCT is very slow over SSH otherwise)
;   Version 2.0, 27-MAY-2016, DD: added keyword NO_HISTO
;   Version 2.1, 28-JUL-2016, DD: properly avoid negative values in plot if /XLOG or /YLOG
;   Version 2.2, 07-DEC-2016, DD: added CALL to MODSDV when /KERNEL is set

PRO PLOTALL, time, data, error, BOLD=bold, KERNEL=kernel, NAME=name, NO_FFT=no_fft, NO_HISTO=no_histo, TAG=tag, SCATTER=scatter, TITLE=title, $
             XLOG=xlog, YLOG=ylog, XRANGE=xrange, YRANGE=yrange, XTITLE=xtitle, YTITLE=ytitle, EPS=eps

 ; Turn off informational messages
 quiet0 = !quiet
 !quiet = 1
           
; Number of elements
n_time = N_ELEMENTS(time)
n_data = N_ELEMENTS(data)
n_err  = N_ELEMENTS(error)
             
; Input and keyword sanity checks
IF n_data NE n_time         THEN MESSAGE, 'Inputs should have the same size'
IF NOT KEYWORD_SET(TAG)     THEN tag_plot  = ''      ELSE tag_plot  = '_' + tag
IF NOT KEYWORD_SET(NAME)    THEN plot_file = 'image' ELSE plot_file = name + tag_plot
IF NOT KEYWORD_SET(XTITLE)  THEN xtitle    = 'x' 
IF NOT KEYWORD_SET(YTITLE)  THEN ytitle    = 'y' 
IF NOT KEYWORD_SET(SCATTER) THEN symbol    = 0 ELSE symbol = 2 ; 2 is a '+' sign

; Running parameters
time_exp = 0.1   ; factor by which the time range will be expanded
data_exp = 0.05  ; factor by which the data range will be expanded

; Prepare input data for later
;time0 = time-MIN(time)
data0 = data-MEAN(data)

; prepare time range
time_max = MAX(time)
time_min = MIN(time)
time_min -= 0.5*time_exp*(time_max-time_min)
time_max += 0.5*time_exp*(time_max-time_min)

; Prepare data range
data_max = MAX(data)
data_min = MIN(data)
data_min -= 0.5*data_exp*(data_max-data_min)
data_max += 0.5*data_exp*(data_max-data_min)

; Remove negative values if /XLOG or /YLOG are set
IF KEYWORD_SET(XLOG) THEN time_min = time_min > 0.5*MIN(time[WHERE(time GT 0)])
IF KEYWORD_SET(YLOG) THEN data_min = data_min > 0.5*MIN(data[WHERE(data GT 0)])


; 1. Compute histograms
; *********************

IF NOT KEYWORD_SET(NO_HISTO) THEN BEGIN
  ; Compute histogram
  n_bin      = ROUND(SQRT(N_ELEMENTS(data)))
  data_hist  = HISTOGRAM(data, NBINS=n_bin, LOCATIONS=bins)
  bins       = bins + 0.5*(bins[1]-bins[0]) ; Shift to bin center
  
  ; Fit histogram
  IF KEYWORD_SET(KERNEL) THEN BEGIN
    ; Sanity check
    IF n_data NE n_err THEN wei = 0 ELSE wei = 1./error^2 
    ; Error bar computation
    MODSDV, data, dat_mod, err_low, err_sup
    ; Compute kernel density
    n_rho      = 1D2*n_bin
    binsize    = (data_max-data_min)/n_rho
    rho_bin    = data_min + binsize*(DINDGEN(n_rho) + 0.5)
    rho_data   = KDE(data, rho_bin, WEIGHT=wei,/GAUSSIAN)
    rho_data   = rho_data/(TOTAL(rho_data)*binsize)*(TOTAL(data_hist)*(bins[1]-bins[0]))  ; normalize density
  ENDIF ELSE gauss_fit = GAUSSFIT(bins,data_hist,fit_param,NTERMS=3,SIGMA=sigma)  ; Gaussian fit
END

; 2. Compute FFT and PSD
; **********************

IF NOT KEYWORD_SET(NO_FFT) THEN COMPUTE_PSD, time, data, FREQ=freq, FRQ_EXT=frq_ext, FFT=fft, PSD=flx_psd, CUM_RMS=cum_rms, SLOPE=coeff, INFO=info, /REGRID


; 3. Plot figures
; ***************

; Initiate plot parameters
IF KEYWORD_SET(EPS) THEN fit = 20./1720. ELSE fit = 1
inv     = 0
!P.FONT = 0 
IF KEYWORD_SET(BOLD) THEN BEGIN
  thick  = 5.0
  xthick = 5.0
  ythick = xthick
  cthick = 3.5
  csize  = 1.3  
ENDIF ELSE BEGIN
  thick  = 4.0
  xthick = 3.5
  ythick = xthick
  cthick = 3
ENDELSE

; -- 3.1. Plot time series --
; ---------------------------

IF NOT KEYWORD_SET(XRANGE) THEN x_range = [time_min,time_max] ELSE x_range = xrange
IF NOT KEYWORD_SET(YRANGE) THEN y_range = [data_min,data_max] ELSE y_range = yrange
IF KEYWORD_SET(YLOG) THEN y_range[0] = y_range[0] > 0.001*(data_max-data_min) 
IF KEYWORD_SET(EPS) THEN PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plot_file + '_DATA.eps' $
                    ELSE PLOTXY, /INIT, INV=inv, /COLOR, WINDOW=[0, 0, 1200, 800]*fit                 
LOADCT, 0, /SILENT
PLOTXY, time, data, YLOG=ylog, /NEW, XRANGE=x_range, YRANGE=y_range, TITLE=title, XTITLE=xtitle, YTITLE=ytitle, GRID=0, $
        XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
LOADCT, 13, /SILENT
PLOTXY, time, data, /ADD, SYMBOL=symbol, COLOR=255, THICK=thick 
LOADCT, 0, /SILENT
PLOTXY, /FIN

; -- 3.2. Plot histograms --
; --------------------------

IF NOT KEYWORD_SET(NO_HISTO) THEN BEGIN
  xrange = [data_min,data_max];[-5,100]
  IF KEYWORD_SET(KERNEL) THEN yrange = [0.,MAX(rho_data)*1.2] ELSE yrange = [0.,MAX(gauss_fit)*1.2] 
  IF KEYWORD_SET(EPS) THEN PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plot_file + '_HIST.eps' $
                      ELSE PLOTXY, /INIT, INV=inv, /COLOR, WINDOW=[0, 0, 1200, 800]*fit
  PLOTXY, bins, data_hist, PSYM=10, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE=title, XTITLE=ytitle, YTITLE= 'Number of occurence', GRID=0, $
          XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
  LOADCT, 13, /SILENT
  PLOTXY, bins, data_hist, /ADD, PSYM=10, COLOR=90, THICK=thick
  IF NOT KEYWORD_SET(KERNEL) THEN BEGIN 
    PLOTXY, bins, gauss_fit, /ADD, COLOR=255, THICK=thick
    XYOUTS, xrange[0]+0.1*(xrange[1]-xrange[0]), yrange[1]-0.1*(yrange[1]-yrange[0]), 'AVG: ' + STRING(fit_param[1], FORMAT='(F8.2)'), CHARSIZE=1.0, CHARTHICK=2.5
    XYOUTS, xrange[0]+0.1*(xrange[1]-xrange[0]), yrange[1]-0.2*(yrange[1]-yrange[0]), 'RMS: ' + STRING(fit_param[2], FORMAT='(F8.2)'), CHARSIZE=1.0, CHARTHICK=2.5
  ENDIF ELSE BEGIN
    ;PLOTXY, rho_bin, rho_data, /ADD, COLOR=255, THICK=thick
    XYOUTS, xrange[0]+0.1*(xrange[1]-xrange[0]), yrange[1]-0.1*(yrange[1]-yrange[0]), 'MODE   : ' + STRING(dat_mod, FORMAT='(F8.2)'), CHARSIZE=1.0, CHARTHICK=2.5
    XYOUTS, xrange[0]+0.1*(xrange[1]-xrange[0]), yrange[1]-0.2*(yrange[1]-yrange[0]), 'ERR_LOW: ' + STRING(err_low, FORMAT='(F8.2)'), CHARSIZE=1.0, CHARTHICK=2.5
    XYOUTS, xrange[0]+0.1*(xrange[1]-xrange[0]), yrange[1]-0.3*(yrange[1]-yrange[0]), 'ERR_SUP: ' + STRING(err_sup, FORMAT='(F8.2)'), CHARSIZE=1.0, CHARTHICK=2.5
  ENDELSE
  LOADCT, 0, /SILENT
  PLOTXY, /FIN
ENDIF

; -- 3.2. Plot FFT amplitude, FFT phase, and PSD --
; -------------------------------------------------
IF NOT KEYWORD_SET(NO_FFT) THEN BEGIN
  ; Remove 0 frequency
  IF freq[0] EQ 0 THEN BEGIN
    freq    = freq[1:*]
    fft     = fft[1:*]
    flx_psd = flx_psd[1:*]
  ENDIF

  ; Compute FFT amplitude
  fft_amp = ABS(fft)
  yrange = [0.5*MIN(fft_amp), 1.5*MAX(fft_amp[1:*])]
  xrange = [MIN(freq),MAX(freq)]
  IF KEYWORD_SET(EPS) THEN PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plot_file + '_FFT_AMP.eps' $
                      ELSE PLOTXY, /INIT, INV=inv, /COLOR, WINDOW=[0, 0, 1200, 800]*fit
  PLOTXY, freq, fft_amp, /XLOG, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE=title, XTITLE='Frequency [Hz]', YTITLE='FFT amplitude', GRID=0, $
          XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
  LOADCT, 13, /SILENT
  PLOTXY, freq, fft_amp, /ADD, LINESTYLE=0, COLOR=255, THICK=thick
  LOADCT, 0, /SILENT
  PLOTXY, /FIN  
  
  ; Compute FFT phase
  fft_pha = ATAN(Imaginary(fft),Real_part(fft))
  yrange = [1.5*MIN(fft_pha), 1.5*MAX(fft_pha[1:*])]
  xrange = [MIN(freq[1:*]),MAX(freq)]
  IF KEYWORD_SET(EPS) THEN PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plot_file + '_FFT_PHA.eps' $
                      ELSE PLOTXY, /INIT, INV=inv, /COLOR, WINDOW=[0, 0, 1200, 800]*fit
  PLOTXY, freq, fft_pha, /XLOG, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE=title, XTITLE='Frequency [Hz]', YTITLE='FFT phase [rad]', GRID=0, $
    XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
  LOADCT, 13, /SILENT
  PLOTXY, freq, fft_pha, /ADD, LINESTYLE=0, COLOR=255, THICK=thick
  LOADCT, 0, /SILENT
  PLOTXY, /FIN
  
  ; Compute PSD
  yrange = [0.5*MIN(flx_psd) > 1D-6*MAX(flx_psd[1:*]), 1.5*MAX(flx_psd[1:*])]
  xrange = [MIN(freq[1:*]), MAX(freq)]
  IF KEYWORD_SET(EPS) THEN PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plot_file + '_PSD.eps' $
                      ELSE PLOTXY, /INIT, INV=inv, /COLOR, WINDOW=[0, 0, 1200, 800]*fit
  PLOTXY, freq, flx_psd, /XLOG, /YLOG, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE=title, XTITLE='Frequency [Hz]', YTITLE='PSD', GRID=0, $
    XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
  LOADCT, 13, /SILENT
  PLOTXY, freq, flx_psd, /ADD, LINESTYLE=0, COLOR=255, THICK=thick
  LOADCT, 0, /SILENT
  PLOTXY, /FIN
ENDIF

; -- 3.3. data versus error --
; ----------------------------

IF n_data EQ n_err AND MAX(data) NE 0 AND MAX(error) GT 0 THEN BEGIN
  xrange = [0.1, 2.*MAX(data)]
  yrange = [0.2*MIN(error[WHERE(error GT 0)]), 2.*MAX(error)]
  IF KEYWORD_SET(EPS) THEN PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plot_file + '_SEN.eps' $
                      ELSE PLOTXY, /INIT, INV=inv, /COLOR, WINDOW=[0, 0, 1200, 800]*fit
  PLOTXY, data, error, XLOG=xlog, YLOG=ylog, /NEW, XRANGE=xrange, YRANGE=yrange,TITLE=title, XTITLE=ytitle, YTITLE='Error', GRID=0, $
          XSTYLE=1, YSTYLE=1, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
  LOADCT, 13, /SILENT
  PLOTXY, data, error, /ADD, PSYM=1, COLOR=90, THICK=thick
  PLOTXY, [1D-4,1D+4], [1D-4,1D+4], /ADD, LINESTYLE=1, COLOR=0, THICK=thick
  LOADCT, 0, /SILENT
  PLOTXY, /FIN
ENDIF

; Restore informational messages
!quiet = quiet0
END




