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
;   Version 2.3, 20-NOC-2024, DD: prevented error with GAUSFFIT

pro PLOTALL, time, data, error, bold = bold, kernel = kernel, name = name, no_fft = no_fft, no_histo = no_histo, tag = tag, scatter = scatter, title = title, $
  xlog = xlog, ylog = ylog, xrange = xrange, yrange = yrange, xtitle = xtitle, ytitle = ytitle, eps = eps
  compile_opt idl2

  ; Turn off informational messages
  quiet0 = !quiet
  !quiet = 1

  ; Number of elements
  n_time = n_elements(time)
  n_data = n_elements(data)
  n_err = n_elements(error)

  ; Input and keyword sanity checks
  if n_data ne n_time then message, 'Inputs should have the same size'
  if not keyword_set(tag) then tag_plot = '' else tag_plot = '_' + tag
  if not keyword_set(name) then plot_file = 'image' else plot_file = name + tag_plot
  if not keyword_set(xtitle) then xtitle = 'x'
  if not keyword_set(ytitle) then ytitle = 'y'
  if not keyword_set(scatter) then symbol = 0 else symbol = 2 ; 2 is a '+' sign

  ; Running parameters
  time_exp = 0.1 ; factor by which the time range will be expanded
  data_exp = 0.05 ; factor by which the data range will be expanded

  ; Prepare input data for later
  ; time0 = time-MIN(time)
  data0 = data - mean(data)

  ; prepare time range
  time_max = max(time)
  time_min = min(time)
  time_min -= 0.5 * time_exp * (time_max - time_min)
  time_max += 0.5 * time_exp * (time_max - time_min)

  ; Prepare data range
  data_max = max(data)
  data_min = min(data)
  data_min -= 0.5 * data_exp * (data_max - data_min)
  data_max += 0.5 * data_exp * (data_max - data_min)

  ; Remove negative values if /XLOG or /YLOG are set
  if keyword_set(xlog) then time_min = time_min > 0.5 * min(time[where(time gt 0)])
  if keyword_set(ylog) then data_min = data_min > 0.5 * min(data[where(data gt 0)])

  ; 1. Compute histograms
  ; *********************

  if not keyword_set(no_histo) then begin
    ; Compute histogram
    n_bin = round(sqrt(n_elements(data)))
    data_hist = histogram(data, nbins = n_bin, locations = bins)
    bins = bins + 0.5 * (bins[1] - bins[0]) ; Shift to bin center
    n_terms = min(n_bin, 3)

    ; Fit histogram
    if keyword_set(kernel) then begin
      ; Sanity check
      if n_data ne n_err then wei = 0 else wei = 1. / error ^ 2
      ; Error bar computation
      MODSDV, data, dat_mod, err_low, err_sup
      ; Compute kernel density
      n_rho = 1d2 * n_bin
      binsize = (data_max - data_min) / n_rho
      rho_bin = data_min + binsize * (dindgen(n_rho) + 0.5)
      rho_data = kde(data, rho_bin, weight = wei, /gaussian)
      rho_data = rho_data / (total(rho_data) * binsize) * (total(data_hist) * (bins[1] - bins[0])) ; normalize density
    endif else gauss_fit = gaussfit(bins, data_hist, fit_param, nterms = n_terms, sigma = sigma) ; Gaussian fit
  end

  ; 2. Compute FFT and PSD
  ; **********************

  if not keyword_set(no_fft) then COMPUTE_PSD, time, data, freq = freq, frq_ext = frq_ext, fft = fft, psd = flx_psd, cum_rms = cum_rms, slope = coeff, info = info, /regrid

  ; 3. Plot figures
  ; ***************

  ; Initiate plot parameters
  if keyword_set(eps) then fit = 20. / 1720. else fit = 1
  inv = 0
  !p.font = 0
  if keyword_set(bold) then begin
    thick = 5.0
    xthick = 5.0
    ythick = xthick
    cthick = 3.5
    csize = 1.3
  endif else begin
    thick = 4.0
    xthick = 3.5
    ythick = xthick
    cthick = 3
  endelse

  ; -- 3.1. Plot time series --
  ; ---------------------------

  if not keyword_set(xrange) then x_range = [time_min, time_max] else x_range = xrange
  if not keyword_set(yrange) then y_range = [data_min, data_max] else y_range = yrange
  if keyword_set(ylog) then y_range[0] = y_range[0] > 0.001 * (data_max - data_min)
  if keyword_set(eps) then PlotXY, /init, inv = inv, /color, /eps, window = [0, 0, 1200, 800] * fit, filename = plot_file + '_DATA.eps' $
  else PlotXY, /init, inv = inv, /color, window = [0, 0, 1200, 800] * fit
  loadct, 0, /silent
  PlotXY, time, data, ylog = ylog, /new, xrange = x_range, yrange = y_range, title = title, xtitle = xtitle, ytitle = ytitle, grid = 0, $
    xstyle = 1, ystyle = 1, xthick = xthick, ythick = ythick, thick = thick, charthick = cthick, /nodata, /noerase, window = [150, 100, 1150, 700] * fit ; , INSET_UR='a.'
  loadct, 13, /silent
  PlotXY, time, data, /add, symbol = symbol, color = 255, thick = thick
  loadct, 0, /silent
  PlotXY, /fin

  ; -- 3.2. Plot histograms --
  ; --------------------------

  if not keyword_set(no_histo) then begin
    xrange = [data_min, data_max] ; [-5,100]
    if keyword_set(kernel) then yrange = [0., max(rho_data) * 1.2] else yrange = [0., max(gauss_fit) * 1.2]
    if keyword_set(eps) then PlotXY, /init, inv = inv, /color, /eps, window = [0, 0, 1200, 800] * fit, filename = plot_file + '_HIST.eps' $
    else PlotXY, /init, inv = inv, /color, window = [0, 0, 1200, 800] * fit
    PlotXY, bins, data_hist, psym = 10, /new, xrange = xrange, yrange = yrange, title = title, xtitle = ytitle, ytitle = 'Number of occurence', grid = 0, $
      xstyle = 1, ystyle = 1, xthick = xthick, ythick = ythick, thick = thick, charthick = cthick, /nodata, /noerase, window = [150, 100, 1150, 700] * fit ; , INSET_UR='a.'
    loadct, 13, /silent
    PlotXY, bins, data_hist, /add, psym = 10, color = 90, thick = thick
    if not keyword_set(kernel) then begin
      PlotXY, bins, gauss_fit, /add, color = 255, thick = thick
      xyouts, xrange[0] + 0.1 * (xrange[1] - xrange[0]), yrange[1] - 0.1 * (yrange[1] - yrange[0]), 'AVG: ' + string(fit_param[1], format = '(F8.2)'), charsize = 1.0, charthick = 2.5
      xyouts, xrange[0] + 0.1 * (xrange[1] - xrange[0]), yrange[1] - 0.2 * (yrange[1] - yrange[0]), 'RMS: ' + string(fit_param[2], format = '(F8.2)'), charsize = 1.0, charthick = 2.5
    endif else begin
      ; PLOTXY, rho_bin, rho_data, /ADD, COLOR=255, THICK=thick
      xyouts, xrange[0] + 0.1 * (xrange[1] - xrange[0]), yrange[1] - 0.1 * (yrange[1] - yrange[0]), 'MODE   : ' + string(dat_mod, format = '(F8.2)'), charsize = 1.0, charthick = 2.5
      xyouts, xrange[0] + 0.1 * (xrange[1] - xrange[0]), yrange[1] - 0.2 * (yrange[1] - yrange[0]), 'ERR_LOW: ' + string(err_low, format = '(F8.2)'), charsize = 1.0, charthick = 2.5
      xyouts, xrange[0] + 0.1 * (xrange[1] - xrange[0]), yrange[1] - 0.3 * (yrange[1] - yrange[0]), 'ERR_SUP: ' + string(err_sup, format = '(F8.2)'), charsize = 1.0, charthick = 2.5
    endelse
    loadct, 0, /silent
    PlotXY, /fin
  endif

  ; -- 3.2. Plot FFT amplitude, FFT phase, and PSD --
  ; -------------------------------------------------
  if not keyword_set(no_fft) then begin
    ; Remove 0 frequency
    if freq[0] eq 0 then begin
      freq = freq[1 : *]
      fft = fft[1 : *]
      flx_psd = flx_psd[1 : *]
    endif

    ; Compute FFT amplitude
    fft_amp = abs(fft)
    yrange = [0.5 * min(fft_amp), 1.5 * max(fft_amp[1 : *])]
    xrange = [min(freq), max(freq)]
    if keyword_set(eps) then PlotXY, /init, inv = inv, /color, /eps, window = [0, 0, 1200, 800] * fit, filename = plot_file + '_FFT_AMP.eps' $
    else PlotXY, /init, inv = inv, /color, window = [0, 0, 1200, 800] * fit
    PlotXY, freq, fft_amp, /xlog, /new, xrange = xrange, yrange = yrange, title = title, xtitle = 'Frequency [Hz]', ytitle = 'FFT amplitude', grid = 0, $
      xstyle = 1, ystyle = 1, xthick = xthick, ythick = ythick, thick = thick, charthick = cthick, /nodata, /noerase, window = [150, 100, 1150, 700] * fit ; , INSET_UR='a.'
    loadct, 13, /silent
    PlotXY, freq, fft_amp, /add, linestyle = 0, color = 255, thick = thick
    loadct, 0, /silent
    PlotXY, /fin

    ; Compute FFT phase
    fft_pha = atan(imaginary(fft), real_part(fft))
    yrange = [1.5 * min(fft_pha), 1.5 * max(fft_pha[1 : *])]
    xrange = [min(freq[1 : *]), max(freq)]
    if keyword_set(eps) then PlotXY, /init, inv = inv, /color, /eps, window = [0, 0, 1200, 800] * fit, filename = plot_file + '_FFT_PHA.eps' $
    else PlotXY, /init, inv = inv, /color, window = [0, 0, 1200, 800] * fit
    PlotXY, freq, fft_pha, /xlog, /new, xrange = xrange, yrange = yrange, title = title, xtitle = 'Frequency [Hz]', ytitle = 'FFT phase [rad]', grid = 0, $
      xstyle = 1, ystyle = 1, xthick = xthick, ythick = ythick, thick = thick, charthick = cthick, /nodata, /noerase, window = [150, 100, 1150, 700] * fit ; , INSET_UR='a.'
    loadct, 13, /silent
    PlotXY, freq, fft_pha, /add, linestyle = 0, color = 255, thick = thick
    loadct, 0, /silent
    PlotXY, /fin

    ; Compute PSD
    yrange = [0.5 * min(flx_psd) > 1d-6 * max(flx_psd[1 : *]), 1.5 * max(flx_psd[1 : *])]
    xrange = [min(freq[1 : *]), max(freq)]
    if keyword_set(eps) then PlotXY, /init, inv = inv, /color, /eps, window = [0, 0, 1200, 800] * fit, filename = plot_file + '_PSD.eps' $
    else PlotXY, /init, inv = inv, /color, window = [0, 0, 1200, 800] * fit
    PlotXY, freq, flx_psd, /xlog, /ylog, /new, xrange = xrange, yrange = yrange, title = title, xtitle = 'Frequency [Hz]', ytitle = 'PSD', grid = 0, $
      xstyle = 1, ystyle = 1, xthick = xthick, ythick = ythick, thick = thick, charthick = cthick, /nodata, /noerase, window = [150, 100, 1150, 700] * fit ; , INSET_UR='a.'
    loadct, 13, /silent
    PlotXY, freq, flx_psd, /add, linestyle = 0, color = 255, thick = thick
    loadct, 0, /silent
    PlotXY, /fin
  endif

  ; -- 3.3. data versus error --
  ; ----------------------------

  if n_data eq n_err and max(data) ne 0 and max(error) gt 0 then begin
    xrange = [0.1, 2. * max(data)]
    yrange = [0.2 * min(error[where(error gt 0)]), 2. * max(error)]
    if keyword_set(eps) then PlotXY, /init, inv = inv, /color, /eps, window = [0, 0, 1200, 800] * fit, filename = plot_file + '_SEN.eps' $
    else PlotXY, /init, inv = inv, /color, window = [0, 0, 1200, 800] * fit
    PlotXY, data, error, xlog = xlog, ylog = ylog, /new, xrange = xrange, yrange = yrange, title = title, xtitle = ytitle, ytitle = 'Error', grid = 0, $
      xstyle = 1, ystyle = 1, /nodata, /noerase, window = [150, 100, 1150, 700] * fit ; , INSET_UR='a.'
    loadct, 13, /silent
    PlotXY, data, error, /add, psym = 1, color = 90, thick = thick
    PlotXY, [1d-4, 1d+4], [1d-4, 1d+4], /add, linestyle = 1, color = 0, thick = thick
    loadct, 0, /silent
    PlotXY, /fin
  endif

  ; Restore informational messages
  !quiet = quiet0
end