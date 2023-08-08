;+
; NAME: ZODI_SENSITIVITY
;   
; PURPOSE:
;   This routine computes the integration time needed to reach a given zodi level as a function of the stellar brightness.
;   This is assuming a solar-like star at 10pc.
;
; KEYWORDS
;   ZODI_GOAL : the zodi sensitivity goal (in zodi)
;
; MODIFICATION HISTORY:
;   Version 1.0,  16-JAN-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu

PRO ZODI_SENSITIVITY, ZODI_GOAL=zodi_goal

  trans = 0.15
  IF NOT KEYWORD_SET(zodi_goal) THEN zodi_goal = 3.

  ; Get fixed parameters
  GET_PRM, prm
  GET_CNF, cnf, INSTRUM='nomic'

  ; Convertion factor
  bdw_freq  = prm.c*cnf.bandwidth/cnf.lambda^2
  Jy2ph     = 1./(1D+26*prm.h*prm.c/cnf.lambda/(2*!Dpi*(0.5*cnf.tel_diam)^2)/bdw_freq)
  print, jy2ph
  
  ; Prepare scaling of stellar flux to the aperture radius used
  ratio = cnf.cen_obs/cnf.tel_diam                     ; secondary to primary ratio 
  n_pt  = 1D+3

  ; Scale flux to the aperture radius
  aper_rad    = 0.5*cnf.psf_fwhm
  theta       = aper_rad/cnf.psf_pix*cnf.psf_rad    ; radial distance from center (rad)
  r           = 2.*!Dpi/cnf.lambda*0.5*cnf.tel_diam*theta
  t           = (1+DINDGEN(n_pt))/n_pt*r
  energy_frac = 1./(1.-ratio^2)*(1.-BESELJ(r,0)^2-BESELJ(r,1)^2+ratio^2*(1-BESELJ(ratio*r,0)^2-BESELJ(ratio*r,1)^2)-$
                4.*ratio*INT_TABULATED(t,BESELJ(t,1)*BESELJ(ratio*t,1)/t))

  ; Background per pixel [Jy]
  bck_pix  = 1.
  zodi_ctr = 5.D-5 * zodi_goal
    
  ; Target flx
  n_mag   = 100.
  tgt_flx = DINDGEN(n_mag) + 1.   ; Jy
  
  ; Compute zodi contrast
  zodi_flx = zodi_ctr * tgt_flx 
  
  ; Convert zodi flux per pixel
  psf_pix   = !Dpi*aper_rad^2   ; number of pixels in PSF
  zodi_flx  = zodi_flx*energy_frac

  ; Convert to photon/s per pixel
  bck_cnt  = bck_pix * Jy2ph * trans
  zodi_cnt = zodi_flx * Jy2ph * trans
   
  ; Compute required integration time
  noise  = SQRT(psf_pix)*SQRT(bck_cnt)    ; over the PSF
  signal = zodi_cnt
  t_int  = noise^2/(signal)^2 
  
  ; Plot result
  LOADCT, 13, /SILENT
  fit     = 20./1720.
  yrange  = [0.001, MAX(t_int)]
  xrange  = [MIN(tgt_flx),MAX(tgt_flx)]
  plot_name = 'zodi_sensitivity.eps'
  PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME=plot_name
  PLOTXY, [0.,1.], [0.,1.], /YLOG, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='', XTITLE='Target magnitude [Jy]', YTITLE='Integration time for 3-zodi sensitivity [sec]', GRID=0, $
          CHARSIZE=1.0, CHARTHICK=2.5, THICK=4.5, XTHICK=4.5, YTHICK=4.5, XSTYLE=1, YSTYLE=1, /NODATA, /NOERASE, WINDOW=[250,100,1350,900]*fit;, INSET_UR='a.'

  ; Overplot the curve for this target (and subtract dark rms)
  PLOTXY, tgt_flx, t_int, /ADD, COLOR=0, THICK=4
  ;PLOTXY, tgt_flx, t_int*3., /ADD, COLOR=0, THICK=4, LINESTYLE=1

  ; Add legend
  ;LEGEND, [tgt_all], LINESTYLE=DBLARR(n_all), CHARTHICK=2.5, CHARSIZE=charsize, COLORS=[0,color], PTHICK=4.0, PSYM=DBLARR(n_all)+8, /TOP_LEGEND, /RIGHT_LEGEND
    
  ; Close plot
  LOADCT, 0, /SILENT
  PLOTXY, /FIN
END