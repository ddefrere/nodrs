;+
; NAME: PLOT_SENSITIVITY
; 
; PURPOSE:
;   Standard routine to plot sensitivity, transmission and backgrounds results
;
; INPUTS
;   date_in       :
;   tgt_name      :
;   int_time      :
;   det_mode      :
;   flx_tot       :  The measured flux
;   flx_err       :  The 1-sigma error on the measured flux
;   bck_flx       :  The background flux
;   bck_err       :  The 1-sigma error on the measured background flux
;
; KEYWORDS
;   LAMBDA        :  Set this keyword to the wavelength in m (10 um by default)
;
; MODIFICATION HISTORY:
;   Version 1.0,  01-MAR-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  24-JUN-2013, DD, added transmission computation and plots

PRO PLOT_SENSITIVITY, date_in, tgt_name, int_time, det_mode, flx_tot, flx_err, bck_flx, bck_err, aper_rad, lambda, bandwidth, nfrseq

; Recover the data path
IF !VERSION.OS_FAMILY EQ 'unix' THEN sep='/' ELSE sep='\'
path          = GET_PATH('lbti_drs.pro', N_DIR_UP=1)             ; Main path of LBTI software
result_path   = path + 'results' + sep + 'sensitivity'           ; Path to the result folder
IF NOT FILE_TEST(result_path) THEN FILE_MKDIR, result_path       ; Create directory if it does not exist
  
; Get fixed parameters
GET_PRM, prm
GET_CNF, cnf

; Efficiency and integration time parameter
dit_tot = 3.
eff     = 0.5
nfrseq  = nfrseq/nfrseq

; Prepare scaling of stellar flux to the aperture radius used
ratio = cnf.cen_obs/cnf.tel_diam                     ; secondary to primary ratio 
n_pt  = 1D+3
  
; Distinct dates
date_uniq = date_in(UNIQ(date_in, SORT(date_in)))      
n_date    = N_ELEMENTS(date_uniq)

; Distinct acquisition times
int_uniq = int_time(UNIQ(int_time, SORT(int_time)))      
n_int    = N_ELEMENTS(int_uniq)

; Compute sensitivity, background, and transmission
sen = DBLARR(n_date,n_int) & sen_err = sen & bck_jy = sen & bck_ejy = sen & trans = sen & adu2jy = trans
bck_date = DBLARR(n_date)  & trans_date = bck_date & sen_meas = bck_date & sen_exp = bck_date 
FOR i_d = 0, n_date-1 DO BEGIN

  ; Extract information for this date
  idx_date    = WHERE(date_in EQ date_uniq[i_d], nd)
  
  ; Get config parameters
  GET_CNF, cnf, LAMBDA=lambda[idx_date[0]], BANDWIDTH=bandwidth[idx_date[0]]
  
  ; Stellar flux (only working for one star per date right now)
  star      = GET_TGT(tgt_name[idx_date[0]])
  star_flux = BLACKBODY(star.temp, lambda[idx_date[0]], STANDARD=0) * !Dpi*(0.5*star.ldm*prm.m2r)^2        ; flux in Jy
  PRINT, 'Target name and flux [Jy]: ', tgt_name[idx_date[0]], star_flux
  
  ; Compute and plot sensitivity
  FOR i_i = 0, n_int-1 DO BEGIN
    ; Retrieve the file obtained at this acquisition time
    idx_int      = WHERE(int_time[idx_date] EQ int_uniq[i_i], nn)
    ; If acquisition time present, compute it. 0 otherwise.
    IF nn GT 0 THEN BEGIN
      idx_cur     = idx_date[idx_int]
      psf_pix     = (0.5*cnf.psf_fwhm)^2*!Dpi
      ; Scale flux to the aperture radius
      theta       = aper_rad[idx_cur[0]]/cnf.psf_pix*cnf.psf_rad    ; radial distance from center (rad)
      r           = 2.*!Dpi/lambda[idx_cur[0]]*0.5*cnf.tel_diam*theta
      t           = (1+DINDGEN(n_pt))/n_pt*r
      energy_frac = 1./(1.-ratio^2)*(1.-BESELJ(r,0)^2-BESELJ(r,1)^2+ratio^2*(1-BESELJ(ratio*r,0)^2-BESELJ(ratio*r,1)^2)-$
                    4.*ratio*INT_TABULATED(t,BESELJ(t,1)*BESELJ(ratio*t,1)/t))
      aper_flux   = star_flux*energy_frac
      ; Convert background to Jy
      bck_jy[i_d,i_i]  = MEAN(bck_flx[idx_cur])*aper_flux/MEAN(flx_tot[idx_cur])
      bck_ejy[i_d,i_i] = MEAN(bck_err[idx_cur])*aper_flux/MEAN(flx_tot[idx_cur])
      bck_ejy_min      = MEAN(bck_err[idx_cur])*aper_flux/MEAN(flx_tot[idx_cur]+flx_err[idx_cur])
      bck_ejy_max      = MEAN(bck_err[idx_cur])*aper_flux/MEAN(flx_tot[idx_cur]-flx_err[idx_cur])
      ; Compute sensitivity in Jy (per frame)
      sen[i_d,i_i]     = bck_ejy[i_d,i_i]*SQRT(psf_pix)
      sen_err[i_d,i_i] = (bck_ejy_max-bck_ejy_min)/2*SQRT(psf_pix)
      ; Compute transmission
      star_ph        = aper_flux*cnf.jy2ph*int_uniq[i_i]*2.*!Dpi*(0.5*cnf.tel_diam)^2
      IF det_mode[idx_cur[0]] EQ 0. THEN det_adu2ph = cnf.adu2ph[0] ELSE det_adu2ph = cnf.adu2ph[1]
      trans[i_d,i_i] = MEAN(flx_tot[idx_cur]*det_adu2ph/(star_ph))
      PRINT, 'flux [ADU]', MEAN(flx_tot[idx_cur])
      ; Compute conversion factor
      adu2jy[i_d,i_i] = aper_flux/(MEAN(flx_tot[idx_cur])/int_uniq[i_i])
    ENDIF 
  ENDFOR
  
  ; Compute mean background and date for this day
  idx_ok          = WHERE(bck_jy[i_d,*] NE 0., n_ok)
  bck_date[i_d]   = MEAN(bck_jy[i_d,idx_ok])
  trans_date[i_d] = MEAN(trans[i_d,idx_ok])
  
  ; Compute sensitivity in 10min for this date...now 3 hours with 50% efficiency (0.23 comes from noise correlation analysis)
  sen_meas[i_d]     = sen[i_d,idx_ok[n_ok-1]]*(int_uniq[idx_ok[n_ok-1]]/(eff*dit_tot*60.*60.))^0.23
  sen_exp[i_d]      = sen[i_d,idx_ok[n_ok-1]]*(int_uniq[idx_ok[n_ok-1]]/(eff*dit_tot*60.*60.))^0.50
  
  ; Print info to screen
  PRINT, ' - Fraction energy in aperture radius: ', energy_frac
  PRINT, ' - Estimated transmission [%]        : ', 100.*trans_date[i_d]
  PRINT, ' - Background [Jy]                   : ', bck_date[i_d]
  PRINT, ' - Sensitivity in 3 hours [mJy]      : ', 1D+3*sen_meas[i_d], 1D+3*sen_exp[i_d]
  PRINT, ' - ADU/s/pix to Jy conversion factor : ', MEAN(adu2jy[i_d,idx_ok])
ENDFOR


; 1. Plot the sensitivity curve 
; *****************************

; Initialize plot
LOADCT, 13, /SILENT
fit     = 20./1720.
color   = REVERSE(50 + 205.*INDGEN(n_date)/(n_date-1))
yrange  = [0.01, 10.]
xrange  = [0.5*MIN(int_time),1.5*MAX(int_time)]
IF n_date EQ 1 THEN plot_name = result_path + sep + date_in[0] + '_sensitivity_nomic.eps' ELSE plot_name = result_path + sep + 'NOMIC_sensitvity_vs_dit.eps'
PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME=plot_name
PLOTXY, [0.,1.], [0.,1.], /XLOG, /YLOG, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='LBTI/NOMIC sensitivity (9.81-12.41' + Greek('mu') + 'm)', XTITLE='Integration time [s]', YTITLE='1-sigma sensitivity [Jy/PSF]', GRID=0, $
        CHARSIZE=1.0, CHARTHICK=2.5, THICK=4.5, XTHICK=4.5, YTHICK=4.5, XSTYLE=1, YSTYLE=9, /NODATA, /NOERASE, WINDOW=[250,100,1350,900]*fit;, INSET_UR='a.'

; Overplot the curve for this date
FOR i_d=0, n_date-1 DO BEGIN
  idx_sen = WHERE(sen[i_d,*] GT 0, n_sen)
  IF n_sen GT 0. THEN BEGIN
    PLOTXY, [int_uniq[idx_sen]], [sen[i_d,idx_sen]], /ADD, LINESTYLE=0, COLOR=color[i_d], THICK=4
    PLOTXY, [int_uniq[idx_sen]], [sen[i_d,idx_sen]], /ADD, SYMBOL=29, COLOR=color[i_d], THICK=4
    ERRPLOT, [int_uniq[idx_sen]], [sen[i_d,idx_sen]-sen_err[i_d,idx_sen]], [sen[i_d,idx_sen]+sen_err[i_d,idx_sen]], COLOR=color[i_d], THICK=4
  ENDIF
ENDFOR

; Add secondary Y axis
psf_pix = (0.5*cnf.psf_fwhm)^2*!Dpi
AXIS, YAXIS=1, /YLOG, YRANGE=!Y.RANGE*SQRT(psf_pix), CHARTHICK=2.5, XTHICK=4.5, YTHICK=4.5, YSTYLE=1, YTITLE='1-sigma sensitivity [Jy/pix]'
      
; Add legend
LEGEND, date_uniq, PSYM=16, /FILL, CHARSIZE=1.0, CHARTHICK=2.5, THICK=4.5, PTHICK=4.5, COLORS = color, /TOP_LEGEND, /RIGHT_LEGEND

; Close plot
LOADCT, 0, /SILENT
PLOTXY, /FIN

; Sensitivity in 10 min
IF n_date GT 1 THEN BEGIN
  ; Plot background in Jy/pixel
  LOADCT, 13, /SILENT
  fit     = 20./1720.
  ;color   = 50 + 205.*INDGEN(n_date)/(n_date-1)
  color   = [250, 120, 130, 80]
  symsize = 0.9
  xrange  = [0., 1.]
  yrange  = [0.3*MIN(sen_exp),2.5*MAX(sen_meas)]*1D+3
  xval    = DINDGEN(n_date)/(n_date) + 1D/(n_date+5D)
  plot_name = result_path + sep + 'NOMIC_sen10_vs_date.eps'
  PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME=plot_name
  PLOTXY, xval, 1D+3*sen_exp, /YLOG, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='LBTI/NOMIC sensitivity (9.81-12.41' + Greek('mu') + 'm)', XTITLE='', YTITLE='Sensitivity in 3 hours [mJy/PSF]', GRID=0, LINESTYLE=0, $
    CHARSIZE=1.0, CHARTHICK=2.5, THICK=4.5, XTHICK=4., YTHICK=4., YSTYLE=9, XTICKS=n_date-1, XTICKV=xval, XTICKNAME=date_uniq, /NOERASE, WINDOW=[250,100,1350,900]*fit;, INSET_UR='a.'
  PLOTXY, xval, 1D+3*sen_exp, SYMBOL=29, COLOR=0, /ADD
  PLOTXY, xval, 1D+3*sen_meas, SYMBOL=29, COLOR=90, /ADD
  PLOTXY, xval, 1D+3*sen_meas, COLOR=90, /ADD
  
  ; Add Visir, Keck, and T-Rex sensitivity
  sen_visir = 5./10.*SQRT(1./dit_tot)*1.00/eff ; http://www.eso.org/sci/facilities/paranal/instruments/visir/inst/index.html
  PLOTXY, [0.,1.], [sen_visir,sen_visir], THICK=4.5, LINESTYLE=2, COLOR=255, /ADD
  sen_keck  = 4.*SQRT(3./dit_tot)  ; Bertrand private communication
  PLOTXY, [0.,1.], [sen_keck,sen_keck], THICK=4.5, LINESTYLE=2, COLOR=255, /ADD
  sen_trecs = 1.4/5.*SQRT(0.5/dit_tot)*0.50/eff ; GEMINI webpage
  PLOTXY, [0.,1.], [sen_trecs,sen_trecs], THICK=4.5, LINESTYLE=2, COLOR=255, /ADD
  
  ; Add LBTI requirements
  sen_lbti = 0.1
  PLOTXY, [0.,1.], [sen_lbti,sen_lbti], THICK=4.5, LINESTYLE=2, COLOR=90, /ADD
  
  ; Add secondary Y axis
  psf_pix = (0.5*cnf.psf_fwhm)^2*!Dpi
  AXIS, YAXIS=1, YRANGE=!Y.RANGE*SQRT(psf_pix), YSTYLE=1, YTITLE='Sensitivity in 3 hours [mJy/pixel]', CHARTHICK=2.5, YTHICK=4.
  
  ; Close plot
  LOADCT, 0, /SILENT
  PLOTXY, /FIN
ENDIF

; 2. Plot the background
; **********************

; Initialize plot
LOADCT, 13, /SILENT
fit     = 20./1720.
;color   = 50 + 205.*INDGEN(n_date)/(n_date-1)
color   = REVERSE(DINDGEN(n_date)/(n_date-1)*255.)
symsize = 0.9
yrange  = [0., 4.]
xrange  = [0.001,1.5*MAX(int_time)]
IF n_date EQ 1 THEN plot_name = result_path + sep + date_in[0] + '_background_nomic.eps' ELSE plot_name = result_path + sep + 'NOMIC_bckg_vs_dit.eps'
PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME=plot_name
PLOTXY, [0.,1.], [0.,1.], /XLOG, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='NOMIC background (9.81-12.41' + Greek('mu') + 'm)', XTITLE='Integration time [s]', YTITLE='Background flux [Jy/pixel]', GRID=0, $
        CHARSIZE=1.0, CHARTHICK=2.5, THICK=4.5, XTHICK=4., YTHICK=4., XSTYLE=1, YSTYLE=9, /NODATA, /NOERASE, WINDOW=[250,100,1350,900]*fit;, INSET_UR='a.'

FOR i_d=0, n_date-1 DO BEGIN
  idx_sen = WHERE(sen[i_d,*] GT 0, n_sen)
  IF n_sen GT 0. THEN BEGIN
    PLOTXY, [int_uniq[idx_sen]], [bck_jy[i_d,idx_sen]], /ADD, SYMBOL=7, COLOR=color[i_d], THICK=4, SYMSIZE=symsize
    ERRPLOT, [int_uniq[idx_sen]], [bck_jy[i_d,idx_sen]-bck_ejy[i_d,idx_sen]], [bck_jy[i_d,idx_sen]+bck_ejy[i_d,idx_sen]], COLOR=color[i_d], THICK=4
  ENDIF
ENDFOR

; Add secondary Y axis
psf_pix = !Dpi*cnf.psf_pix^2
AXIS, YAXIS=1, YRANGE=!Y.RANGE*psf_pix, YSTYLE=1, YTITLE='Background flux [Jy/PSF]', CHARTHICK=2.5, YTHICK=4.

; Add legend
LEGEND, date_uniq, LINESTYLE=DBLARR(n_date), CHARTHICK=1.8, PTHICK=4., CHARSIZE=charsize, COLORS=color, PSYM=DBLARR(n_date)+8, SYMSIZE=DBLARR(n_date)+symsize, /TOP_LEGEND, /RIGHT_LEGEND
    
; Close plot
LOADCT, 0, /SILENT
PLOTXY, /FIN

; 3. Plot the background vs date (Jy and ph/s)
; *********************************************
IF n_date GT 1 THEN BEGIN
  ; Plot background in Jy/pixel
  LOADCT, 13, /SILENT
  fit     = 20./1720.
  ;color   = 50 + 205.*INDGEN(n_date)/(n_date-1)
  color   = [250, 120, 130, 80]
  symsize = 0.9
  yrange  = [0., 4.]
  xrange  = [0.,1.]
  xval    = DINDGEN(n_date)/(n_date) + 1D/(n_date+5D)
  plot_name = result_path + sep + 'NOMIC_bckg-jy_vs_date.eps'
  PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME=plot_name
  PLOTXY, xval, bck_date, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='NOMIC background (9.81-12.41' + Greek('mu') + 'm)', XTITLE='', YTITLE='Background flux [Jy/pixel]', GRID=0, $
          CHARSIZE=1.0, CHARTHICK=2.5, THICK=4.5, XTHICK=4., YTHICK=4., YSTYLE=9, XTICKS=n_date-1, XTICKV=xval, XTICKNAME=date_uniq, SYMBOL=3, /NOERASE, WINDOW=[250,100,1350,900]*fit;, INSET_UR='a.'
                   
  ; Add secondary Y axis
  psf_pix = !Dpi*cnf.psf_pix^2
  AXIS, YAXIS=1, YRANGE=!Y.RANGE*psf_pix, YSTYLE=1, YTITLE='Background flux [Jy/PSF]', CHARTHICK=2.5, YTHICK=4.
  
  ; Close plot
  LOADCT, 0, /SILENT
  PLOTXY, /FIN   
  
;  ; Plot background in photon/pixels/s
;  LOADCT, 13, /SILENT
;  fit     = 20./1720.
;  ;color   = 50 + 205.*INDGEN(n_date)/(n_date-1)
;  color   = [250, 120, 130, 80]
;  symsize = 0.9
;  yrange  = [0., 4000.]
;  xrange  = [0.,1.]
;  xval    = DINDGEN(n_date)/(n_date) + 1D/(n_date+5D)
;  plot_name = result_path + sep + 'NOMIC_bckg-ph_vs_date.eps'
;  PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME=plot_name
;  PLOTXY, xval, bck_adu, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='NOMIC -- N-band background', XTITLE='', YTITLE='Background flux [ph/pixel/27ms]', GRID=0, $
;    CHARSIZE=1.0, CHARTHICK=2.5, THICK=4.5, XTHICK=4., YTHICK=4., YSTYLE=1, XTICKS=n_date-1, XTICKV=xval, XTICKNAME=date_uniq, SYMBOL=3, /NOERASE, WINDOW=[250,100,1350,900]*fit;, INSET_UR='a.'
;  
;  ; Close plot
;  LOADCT, 0, /SILENT
;  PLOTXY, /FIN           
ENDIF

; 4. Plot the transmission vs date
; ********************************
IF n_date GT 1 THEN BEGIN
  ; Initialize plot
  LOADCT, 13, /SILENT
  fit     = 20./1720.
  ;color   = 50 + 205.*INDGEN(n_date)/(n_date-1)
  color   = [250, 120, 130, 80]
  symsize = 0.9
  yrange  = [0., 0.3]
  xrange  = [0.,1.]
  xval    = DINDGEN(n_date)/(n_date) + 1D/(n_date+5D)
  plot_name = result_path + sep + 'NOMIC_trans_vs_date.eps'
  PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME=plot_name
  PLOTXY, xval, trans_date, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='NOMIC transmission (19.81-12.41' + Greek('mu') + 'm)', XTITLE='', YTITLE='Transmission', GRID=0, $
    CHARSIZE=1.0, CHARTHICK=2.5, THICK=4.5, XTHICK=4., YTHICK=4., YSTYLE=1, XTICKS=n_date-1, XTICKV=xval, XTICKNAME=date_uniq, SYMBOL=3, /NOERASE, WINDOW=[250,100,1350,900]*fit;, INSET_UR='a.'

  ; Close plot
  LOADCT, 0, /SILENT
  PLOTXY, /FIN
ENDIF
END