;+
; NAME: CHI2STAT
;
; PURPOSE:
;   Compute best-fit values (chi2 minimization and Bayesian approach) based on a input chi2 cube.
;   Compute the error bars using the statistics of the Pearson chi2 and some bayesian approach
;
; CALLING SEQUENCE:
;   CHI2STAT, chi2cube_coarse, null_expl_c, mu_expl_c, mu_expl_c, nas_chi2_coarse, mu_chi2_coarse, sig_chi2_coarse, DOF=dof, PLOT_PATH=plot_path, TAG=tag
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

PRO CHI2STAT, chi2cube, null_expl, mu_expl, sig_expl, nas, mu, sig, DOF=dof, PLOT_PATH=plot_path, TAG=tag

  ; Keyword sanity check (this assumes that the chi2 cube is unreduced)
  IF NOT KEYWORD_SET(DOF) THEN dof = 1
  
  ; 1. Find optimum position
  ; ************************
  
  tmp     = MIN(chi2cube, index)
  idx3D   = ARRAY_INDICES(chi2cube, index) 
  nas_opt = null_expl[idx3D[0]]
  mu_opt  = mu_expl[idx3D[1]]
  sig_opt = sig_expl[idx3D[2]]
  
  ; 2. Bayesian approach
  ; ********************
     
  ; Convert to non-reduced chi2
  chi2_cube = DOUBLE(chi2cube*dof)
  
  ; Computes proba density for chi-sq distribution
  q = !quiet
  !quiet = 1
  pro_den = CHI2_PDF(chi2_cube,DF=dof) ; chi2_pdf Computes proba density for chi-sq distribution
  !quiet = q
  
  ; Now marginalize
  nas_pdf = TOTAL(TOTAL(pro_den,3),2)/TOTAL(pro_den)
  mu_pdf  = TOTAL(TOTAL(pro_den,3),1)/TOTAL(pro_den)
  sig_pdf = TOTAL(TOTAL(pro_den,2),1)/TOTAL(pro_den)
  
  ; Now compute mode and error bar using high-density regions described by Hyndman, R. J. 1996, The American Statistician, 50, 120
  PDFSDV, nas_pdf, null_expl, nas_bay, nas_err_low, nas_err_sup, BAR_H=bar_nas
  PDFSDV, mu_pdf, mu_expl, mu_bay, mu_err_low, mu_err_sup, BAR_H=bar_mu
  PDFSDV, sig_pdf, sig_expl, sig_bay, sig_err_low, sig_err_sup, BAR_H=bar_sig
  
  ; Output structure
  nas = {NAS_OPT: nas_opt, NAS_BAY: nas_bay, NAS_ERR_LOW: nas_err_low, NAS_ERR_SUP: nas_err_sup}
  mu  = {mu_OPT: mu_opt, mu_BAY: mu_bay, mu_ERR_LOW: mu_err_low, mu_ERR_SUP: mu_err_sup}
  sig = {sig_OPT: sig_opt, sig_BAY: sig_bay, sig_ERR_LOW: sig_err_low, sig_ERR_SUP: sig_err_sup}
  
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