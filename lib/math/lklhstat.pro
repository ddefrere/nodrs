;+
; NAME: LKLHSTAT
;
; PURPOSE:
;   Compute best-fit values (likelihood maximization and Bayesian approach) based on a input chi2 cube.
;   Compute the error bars using the statistics of the Pearson chi2 and some bayesian approach
;
; CALLING SEQUENCE:
;  LKLHSTAT, lklhcube, null_expl, mu_expl, sig_expl, nas, mu, sig, PLOT_PATH=plot_path, TAG=tag
;  
; INPUTS:
;
;
; OUTPUTS:
;   nas_chi2_coarse, mu_chi2_coarse, sig_chi2_coarse
;
; MODIFICATION HISTORY:
;  Version 1.0, 28-MAY-2016, by Denis Defrère, University of Arizona, denis@lbti.org
;-

PRO LKLHSTAT, lklhcube, p1, p2, p3, p4, p5, p1_out, p2_out, p3_out, p4_out, p5_out, PLOT_PATH=plot_path, TAG=tag

  ; 1. Find optimum position
  ; ************************

  tmp      = MAX(lklhcube, index)
  idx3D    = ARRAY_INDICES(lklhcube, index)
  p1_opt   = p1[idx3D[0]]
  p2_opt   = p2[idx3D[1]]
  p3_opt   = p3[idx3D[2]]
  p4_opt   = p4[idx3D[3]]
  p5_opt   = p5[idx3D[4]]

  ; 2. Bayesian approach
  ; ********************

  ; Now marginalize
  norm   = TOTAL(lklhcube)
  p1_pdf = TOTAL(TOTAL(TOTAL(TOTAL(lklhcube,5),4),3),2)/norm
  p2_pdf = TOTAL(TOTAL(TOTAL(TOTAL(lklhcube,5),4),3),1)/norm
  p3_pdf = TOTAL(TOTAL(TOTAL(TOTAL(lklhcube,5),4),2),1)/norm
  p4_pdf = TOTAL(TOTAL(TOTAL(TOTAL(lklhcube,5),3),2),1)/norm
  p5_pdf = TOTAL(TOTAL(TOTAL(TOTAL(lklhcube,4),3),2),1)/norm 

  ; Now compute mode and error bar using high-density regions described by Hyndman, R. J. 1996, The American Statistician, 50, 120
  PDFSDV, p1_pdf, p1, p1_bay, p1_err_low, p1_err_sup, BAR_H=p1_bar
  PDFSDV, p2_pdf, p1, p2_bay, p2_err_low, p2_err_sup, BAR_H=p2_bar
  PDFSDV, p3_pdf, p1, p3_bay, p3_err_low, p3_err_sup, BAR_H=p3_bar
  PDFSDV, p4_pdf, p1, p4_bay, p4_err_low, p4_err_sup, BAR_H=p4_bar
  PDFSDV, p5_pdf, p1, p5_bay, p5_err_low, p5_err_sup, BAR_H=p5_bar
  
  ; Output structure
  p1_out = {OPT: p1_opt, BAY: p1_bay, ERR_LOW: p1_err_low, ERR_SUP: p1_err_sup}
  p2_out = {OPT: p2_opt, BAY: p2_bay, ERR_LOW: p2_err_low, ERR_SUP: p2_err_sup}
  p3_out = {OPT: p3_opt, BAY: p3_bay, ERR_LOW: p3_err_low, ERR_SUP: p3_err_sup}
  p4_out = {OPT: p4_opt, BAY: p4_bay, ERR_LOW: p4_err_low, ERR_SUP: p4_err_sup}
  p5_out = {OPT: p5_opt, BAY: p5_bay, ERR_LOW: p5_err_low, ERR_SUP: p5_err_sup}

  ; Plot results
  IF KEYWORD_SET(PLOT_PATH) THEN BEGIN
; Plot name and parameters
    plotname = plot_path + 'PDF'
    IF NOT KEYWORD_SET(TAG) THEN tag = ' '  
    ; PDFs
    ; Initiate plot parameters
    fit = 20./1720.
    inv     = 0
    !P.FONT = 0
    thick  = 4.0
    xthick = 3.5
    ythick = xthick
    cthick = 3

    ; 1. NAS PDF
    x_range = [MIN(null_expl),MAX(null_expl)]*1D2 
    y_range = [MIN(nas_pdf),MAX(nas_pdf)] 
    PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plotname + '_NAS_' + tag + '.eps' 
    LOADCT, 0, /SILENT
    PLOTXY, 1D2*null_expl, nas_pdf, YLOG=ylog, /NEW, XRANGE=x_range, YRANGE=y_range, TITLE=title, XTITLE='Tested null [%]', YTITLE='Marginalized PDF', GRID=0, $
      XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
    LOADCT, 13, /SILENT
    PLOTXY, 1D2*null_expl, nas_pdf, /ADD, SYMBOL=symbol, COLOR=255, THICK=thick
    PLOTXY, 1D2*null_expl, bar_nas+INTARR(N_ELEMENTS(null_expl)), /ADD, SYMBOL=symbol, COLOR=0, THICK=thick, LINESTYLE=1
    LOADCT, 0, /SILENT
    PLOTXY, /FIN
    
    ; 1. mu PDF
    x_range = [MIN(mu_expl),MAX(mu_expl)]
    y_range = [MIN(mu_pdf),MAX(mu_pdf)]
    PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plotname + '_MU_' + tag + '.eps' 
    LOADCT, 0, /SILENT
    PLOTXY, mu_expl, mu_pdf, YLOG=ylog, /NEW, XRANGE=x_range, YRANGE=y_range, TITLE=title, XTITLE='Tested mean phase [rad]', YTITLE='Marginalized PDF', GRID=0, $
      XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
    LOADCT, 13, /SILENT
    PLOTXY, mu_expl, mu_pdf, /ADD, SYMBOL=symbol, COLOR=255, THICK=thick
    PLOTXY, mu_expl, bar_mu+INTARR(N_ELEMENTS(mu_expl)), /ADD, SYMBOL=symbol, COLOR=0, THICK=thick, LINESTYLE=1
    LOADCT, 0, /SILENT
    PLOTXY, /FIN
    
    ; 1. mu PDF
    x_range = [MIN(sig_expl),MAX(sig_expl)]
    y_range = [MIN(sig_pdf),MAX(sig_pdf)]
    PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plotname + '_SIG_' + tag + '.eps' 
    LOADCT, 0, /SILENT
    PLOTXY, sig_expl, sig_pdf, YLOG=ylog, /NEW, XRANGE=x_range, YRANGE=y_range, TITLE=title, XTITLE='Tested phase jitter [rad]', YTITLE='Marginalized PDF', GRID=0, $
      XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
    LOADCT, 13, /SILENT
    PLOTXY, sig_expl, sig_pdf, /ADD, SYMBOL=symbol, COLOR=255, THICK=thick
    PLOTXY, sig_expl, bar_sig+INTARR(N_ELEMENTS(sig_expl)), /ADD, SYMBOL=symbol, COLOR=0, THICK=thick, LINESTYLE=1
    LOADCT, 0, /SILENT
    PLOTXY, /FIN
  ENDIF

END