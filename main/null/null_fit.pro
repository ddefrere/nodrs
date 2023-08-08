;+
; NAME: NULL_FIT
; 
; PURPOSE:
;   This function fit null data given as input.
;
; INPUTS:
;   data_path     :  Path to the level 2 fits files (calibrated sequences of null).
;
; KEYWORDS
;    L1TYPE        : set this keyword to the type of L1 file to fit (OB, PT, or BIAS). Pt is default.
;    PLOT          : set this keyword to have plot the results
;    PS            : set this keyword to save the plot in an eps figure
;
; MODIFICATION HISTORY:
;   Version 1.0,  16-JAN-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  14-APR-2015, DD: corrected bug on statistical error computation
;   Version 1.2,  14-APR-2015, DD: modified for ORR version + added LOG_FILE
;   Version 1.3,  14-APR-2015, DD: now read HA from L! file rather than computing it
;   Version 1.4,  04-NOV-2016, DD: added keyword L1 type
;   Version 1.5,  04-APR-2017, DD: adapted for new formlaism of GET_TGT
;   Version 1.6,  03-AUG-2017, DD: improved output formalism
;   Version 1.7,  05-AUG-2017, DD: now use AVGERR instead of AVGSDV
;   Version 1.8,  29-SEP-2017, DD: updated touse the new keywords in the L2 files (i.e., NULL_EST and SIG_SCA)
;   Version 1.9,  29-JUL-2018, DD: added comments to error computation
;   Version 2.0,  31-JUL-2018, DD: improved consitency with NULL_CALIB.pro

PRO NULL_FIT, date, cfg_file, INFO=info, LOG_FILE=log_file, L1TYPE=l1type, PLOT=plot, PS=ps
   
; Start actual code
ON_ERROR, 0
IF NOT KEYWORD_SET(info) THEN info = 1
drs_version = 2.0
drs_date    = '31-JUL-2018'

; DEFINE GLOBAL VARIABLES
; ***********************

; Astronomical and physical constants
GET_PRM, prm

; Obtain the definition of the configuration
GET_CNF, cnf, INSTRUM='nomic'

; Read config file with reduction parameters
GET_DRS, drs, 'nodrs/cfg/' + cfg_file

; Recover the IDL running path
DECLARE_PATH, pth, INSTRUM='nomic'

; COMPUTE FILE PATH
; *****************

; Get long version of the date
date_lng    = '20' + STRMID(date, 0, 2) + '-' + STRMID(date, 2, 2) + '-' + STRMID(date, 4, 2)

; Define operational parameters
COMMON LBTI_PARAM, tel_diam, base, lambda, pix_size, psf_pix

; L2 directory
IF MAX(drs.null_rad) NE 0 AND NOT STRMATCH(drs.dir_label, '*_APR') THEN drs.dir_label = drs.dir_label + '_APR'
l2_dir = pth.l2fits_path + pth.sep + date_lng + drs.dir_label + pth.sep ;+ '_APR' + pth.sep
IF NOT FILE_TEST(l2_dir) THEN FILE_MKDIR, l2_dir

; Init log file
IF KEYWORD_SET(log_file) THEN BEGIN
  log_file =  l2_dir + date_lng + '.txt'
  OPENW, lun, log_file, /GET_LUN, WIDTH=800, /APPEND
  PRINTF,lun, ' '
  PRINTF,lun, '******************************************************************************************************* '
  PRINTF,lun, ' '
  PRINTF,lun, 'NULL_FIT.pro version ' + STRING(drs_version, FORMAT='(F3.1)') +  ' -- ' + drs_date + ' -- Denis Defrère - Steward Observatory (denis@lbti.org)'
  PRINTF,lun, 'Data fitted on ' + SYSTIME()
  PRINTF,lun, ' '
  PRINTF,lun, '******************************************************************************************************* '
  PRINTF,lun, ' '
ENDIF ELSE lun = -1


; START CALIBRATION
; *****************

; Search for L1 file and read it
IF drs.cal_mode EQ 1 THEN file_type = 'OB' ELSE file_type = 'PT'
IF KEYWORD_SET(L1TYPE) THEN file_type = file_type + '-' + l1type
l2file = l2_dir + 'UT' + date_lng + '_calib_' + file_type + '.fits'
IF NOT FILE_TEST(l2file) THEN RETURN ; NOT FOUND
l2data = MRDFITS(l2file, 1, header, /SILENT) & n_data = N_ELEMENTS(l2data.mjd_obs)
tmp    = MRDFITS(l2file, 0, header, /SILENT)

; Read header information
null_est = FXPAR(header, 'NULL_EST', /NOCONTINUE)
sig_sca  = FXPAR(header, 'SIG_SCA', /NOCONTINUE)
cal_mode = FXPAR(header, 'CAL_MODE', /NOCONTINUE)
;bck_mode = FXPAR(header, 'BCK_MODE', /NOCONTINUE)

; Print some info to screen
aper_rad = l2data.aper_rad
bck_irad = l2data.bck_irad
bck_orad = l2data.bck_orad
flag     = l2data.QUALITY_FLG
;IF info GT 0 THEN PRINT, 'Aperture radius : ' + STRING(aper_rad, FORMAT='(I0)')

; Read data
wav  = MEAN(l2data.wav_eff)
bdw  = MEAN(l2data.bandwidth)
tgt  = l2data.objname

; Wavelength range
n_lam       = 20
lambda_chnl = wav - 0.5*bdw + bdw * (1.+1./n_lam) * DINDGEN(n_lam)/(n_lam-1)

; Extract LBTI transmission
trans = 1.      ; Temporary implementation

; Find unique targets
tgt_uniq = tgt(UNIQ(tgt,  SORT(tgt)))
n_tgt    = N_ELEMENTS(tgt_uniq)

; Loop over the targets
FOR i_tgt = 0, n_tgt -1 DO BEGIN
  
  ; extract target
  idx_tgt      = WHERE(tgt EQ tgt_uniq[i_tgt], n_sci)
  sci_time     = l2data[idx_tgt].mjd_obs
  sci_para     = l2data[idx_tgt].lbt_para
  sci_ha       = l2data[idx_tgt].hr_angle
  idx_neg      = WHERE(sci_para GT -180 AND sci_para LT -90, nn)
  ;IF nn GT 0 THEN sci_para[idx_neg] += 360
  null_sci     = l2data[idx_tgt].null_cal ;- 0.005
  null_esyst   = l2data[idx_tgt].null_cal_sys
  null_estat   = l2data[idx_tgt].null_cal_sta
  null_etot    = l2data[idx_tgt].null_cal_err
  null_eraw    = l2data[idx_tgt].null_meas_err
  null_flo_err = l2data[idx_tgt].null_flo_err

  ; Extract target information
  GET_TGT, tgt_uniq[i_tgt], star, DATABASE= pth.input_path + 'hosts_db.dat'
  th_mean  = star.ldm
  th_min   = star.ldm - star.lde
  th_max   = star.ldm + star.lde
  uld      = star.uld
  
  ; Convert the time to UT hours
  time_mean = MEAN(sci_time)
  DAYCNV, time_mean+2400000.5D0, yr, mn, day, hr
  ;title_day = STRING(yr,FORMAT='(I0)') + '/' + STRING(mn,FORMAT='(I0)') + '/' + STRING(day,FORMAT='(I0)')
  JDCNV, yr, mn, day, 0., t0day
  t0day = t0day - 2400000.5D0
  scitime  = sci_time - t0day
    
  ; Generate stellar spectrum for broadband null computation
  star_flux = BLACKBODY(star.temp, lambda_chnl, STANDARD=3)             ; flux in photons/s/m²/Hz/sr  
  
  ; Compute the expected null for this star
  norm        = TOTAL((trans*star_flux)^2)
  null_lam    = GEOMETRIC_NULL(th_mean, cnf.base, lambda_chnl, ULD=uld)
  null_wb     = TOTAL((trans*star_flux)^2 * null_lam) / norm
  ; Compute minimum null (with minimum diameter)
  null_min    = GEOMETRIC_NULL(th_min, cnf.base, lambda_chnl, ULD=uld)
  null_wb_min = TOTAL((trans*star_flux)^2 * null_min) / norm
  ; Compute maximum null (with maximum diameter)
  null_max    = GEOMETRIC_NULL(th_max, cnf.base, lambda_chnl, ULD=uld)
  null_wb_max = TOTAL((trans*star_flux)^2 * null_max) / norm
  ; Derive the error bar on the estimated null
  null_err    = (null_wb_max - null_wb_min) / 2D0
   
  ; Compute statistical error bar on SCI OBs
  rms_stat   = 1./SQRT(TOTAL(1./null_eraw^2))                           ; statistical error on SCI measurements 
  
  ; Compute systematic and total error bar on TF
  rms_tf_tot = TOTAL(null_flo_err)/n_sci
  rms_syst   = SQRT(TOTAL(null_esyst/n_sci)^2 + null_err^2)               
  
  ; Compute weighted-average null and error bars
  IF n_sci LE 1 THEN BEGIN
    null_avg       = null_sci
    disp_m_null    = 0
    disp_m_null_we = 0
    rms_exc        = 0.
  ENDIF ELSE BEGIN
    AVGERR, null_sci, null_etot, null_avg_unw, null_avg_wei, rms_stat_null, disp_m_null, disp_m_null_we, rms_exc, KAPPA=sig_sca, IDX_OUT=idx_out, N_OUT=n_out
    IF null_est EQ 0 THEN null_avg = null_avg_unw ELSE null_avg = null_avg_wei
  ENDELSE
  
  ; Definition of systematic errors
  CASE (ABS(drs.err_mode) MOD 3) OF
    0: rms_sys = disp_m_null
    1: rms_sys = disp_m_null_we
    2: rms_sys = rms_exc
    ELSE: MESSAGE, 'Undefined error mode'
  ENDCASE
  
  ; Add errors (if drs.err_mode GT 3, use only the greater of the two)
  IF drs.err_mode GE 3 THEN rms_tot = SQRT(null_err^2 + rms_tf_tot^2 + rms_stat^2) > SQRT(null_err^2 + rms_tf_tot^2 + rms_sys^2) $
                       ELSE rms_tot = SQRT(null_err^2 + rms_tf_tot^2 + rms_stat^2 + rms_sys^2) ; null err is the error on the stellar diameter, rms_tf_tot is the total error on the TF at that position, rms_stat is the stat error on the SCI OBs, and rms_sys is the SYS error
  
  ; Compute stat error on TF
  rms_tf_stat = MEAN(SQRT(null_estat^2-null_eraw^2))
 
  ; Print results
  IF info GT 0 THEN BEGIN
    PRINT, tgt_uniq[i_tgt] + ' [%]:'
    PRINT, 'CALIBRATION PARAMETERS'
    PRINT, '   - Photometric aperture  :', aper_rad[i_tgt] 
    PRINT, '   - Calibration mode      :', cal_mode
    PRINT, '   - Null estimator        :', null_est 
    PRINT, '   - Sigma threshold limit :', sig_sca
    PRINT, 'CALIBRATED NULL'   
    PRINT, '   Estimated null for photopshere         : ', STRING((null_wb)*1D2, FORMAT='(F8.5)') + ' +/- ' + STRING(null_err*1D2, FORMAT='(F7.5)')
    PRINT, '   Calibrated null (photosphere corected) : ', STRING((null_avg-null_wb)*1D2, FORMAT='(F8.5)') + ' +/- ' + STRING(rms_tot*1D2, FORMAT='(F7.5)')
    PRINT, '   - Statistical error on SCI OBs      : ', STRING(rms_stat*1D2, FORMAT='(F7.5)')
    PRINT, '   - Scatter error on SCI OBs          : ', STRING(disp_m_null*1D2, FORMAT='(F7.5)')
    PRINT, '   - Error due to diameter uncertainty : ', STRING(null_err*1D2, FORMAT='(F7.5)')
    PRINT, '   - Error on TF                       : ', STRING(rms_tf_tot*1D2, FORMAT='(F7.5)')
    IF lun GT 0 THEN BEGIN
      PRINTF, lun, ' '
      PRINTF, lun, tgt_uniq[i_tgt] + ' [%]:'
      PRINTF, lun,  'CALIBRATION PARAMETERS'
      PRINTF, lun,  '   - Photometric aperture  :', aper_rad[i_tgt] 
      PRINTF, lun,  '   - Calibration mode      :', cal_mode
      PRINTF, lun,  '   - Null estimator        :', null_est 
      PRINTF, lun,  '   - Sigma threshold limit :', sig_sca
      PRINTF, lun,  'CALIBRATED NULL'   
      PRINTF, lun,  '   Estimated null for photopshere         : ', STRING((null_wb)*1D2, FORMAT='(F8.5)') + ' +/- ' + STRING(null_err*1D2, FORMAT='(F7.5)')
      PRINTF, lun,  '   Calibrated null (photosphere corected) : ', STRING((null_avg-null_wb)*1D2, FORMAT='(F8.5)') + ' +/- ' + STRING(rms_tot*1D2, FORMAT='(F7.5)')
      PRINTF, lun,  '   - Statistical error on SCI OBs      : ', STRING(rms_stat*1D2, FORMAT='(F7.5)')
      PRINTF, lun,  '   - Scatter error on SCI OBs          : ', STRING(disp_m_null*1D2, FORMAT='(F7.5)')
      PRINTF, lun,  '   - Error due to diameter uncertainty : ', STRING(null_err*1D2, FORMAT='(F7.5)')
      PRINTF, lun,  '   - Error on TF                       : ', STRING(rms_tf_tot*1D2, FORMAT='(F7.5)')
    ENDIF
    
    ; Write summary file
    IF KEYWORD_SET(log_file) THEN BEGIN
      ; Create log file
      ;log_file2 =  pth.l2fits_path + pth.sep + 'results_hosts.txt'
      log_file2 =  pth.l2fits_path + pth.sep + 'results_hosts2018.txt'
      IF NOT FILE_TEST(log_file2) THEN BEGIN
        OPENW, lun2, log_file2, /GET_LUN, WIDTH=800, /APPEND
        PRINTF, lun2, 'Column signification'
        PRINTF, lun2, '1: Date'
        PRINTF, lun2, '2: Target'
        PRINTF, lun2, '3: Aperture radius'
        PRINTF, lun2, '4: Background outer radius'
        PRINTF, lun2, '5: Background inner radius'
        PRINTF, lun2, '6: Null estimator mode'
        PRINTF, lun2, '7: Number of sigmas for scatter rejection of OBs'
        PRINTF, lun2, '8: Calibration mode (0 per pointing; 1 per OB)'
        PRINTF, lun2, '9: Calibrated null [%]'
        PRINTF, lun2, '10: Error on calibrated null [%]'
        PRINTF, lun2, '11: Quality flag'
        PRINTF, lun2, '12: L1 file version'
        PRINTF, lun2, '13: Execution date'
        PRINTF, lun2, ' 1 ', '  ', '   2  ','  ','                  3     ','  ','   4   ','  ','   5   ','  ','   6   ','  ','   7   ','  ','   8   ','  ','   9   ','  ','   10   ','   11   ','  ','   12   ', '   13   '
        PRINTF, lun2, '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
      ENDIF ELSE OPENW, lun2, log_file2, /GET_LUN, WIDTH=800, /APPEND
      PRINTF,lun2, date, '  ', tgt_uniq[i_tgt], aper_rad[idx_tgt[0]], '  ', bck_irad[idx_tgt[0]], '  ', bck_orad[idx_tgt[0]], '  ', null_est, '  ', sig_sca, '  ', cal_mode, '  ', $
                   STRING((null_avg-null_wb)*1D2, FORMAT='(F8.5)'), '  ', STRING(rms_tot*1D2, FORMAT='(F7.5)'),'  ',  flag[i_tgt], '  ', 0,  '  ', SYSTIME()
      CLOSE, lun2
      FREE_LUN, lun2
    ENDIF
  ENDIF
  
  ; Fit a binary model (very simple approach!!)
  n_tet = 1.   & tet_max = 500. * prm.m2r
  n_phi = 1.   & phi_max = !Dpi
  n_ctr = 100. & ctr_max = 0.20
  tet = 252.21 * prm.m2r;DINDGEN(n_tet)/(n_tet-1.) * tet_max  ; Angular separation
  phi = 244. * !DTOR ;DINDGEN(n_phi)/(n_phi-1.) * phi_max  ; Position angle
  ctr = DINDGEN(n_ctr)/(n_ctr-1.) * ctr_max
;  n_tet = 50. & tet_max = 500. * prm.m2r
;  n_phi = 50. & phi_max = !Dpi
;  n_ctr = 50. & ctr_max = 0.20
;  tet = DINDGEN(n_tet)/(n_tet-1.) * tet_max  ; Angular separation
;  phi = DINDGEN(n_phi)/(n_phi-1.) * phi_max  ; Position angle
;  ctr = DINDGEN(n_ctr)/(n_ctr-1.) * ctr_max
  
  ; Target visibilities
  v1 = (1-ABS(null_wb))/(1+ABS(null_wb))
  v2 = 1.00
  
  ; Measured baseline position angles
  phi_base = ATAN(l2data[idx_tgt].u_coord,l2data[idx_tgt].v_coord)
  
  ; Predict baseline position angle based on hour angle (consitancy check)
  ra  = star.ra
  dec = star.dec
  ra*=15
  precess, ra, dec, 2000.0, 2014.0
  ra/=15
  par_ang = parangle(sci_ha, dec, TEN(32,42,05.87))
  ;PRINT, par_ang-sci_para 
  ;PRINT, sci_para
  
  ; Loop over the parameters
  r_chi2 = DBLARR(n_tet,n_phi,n_ctr)
  FOR i_tet = 0, n_tet-1 DO BEGIN
      FOR i_phi = 0, n_phi-1 DO BEGIN
          FOR i_ctr = 0, n_ctr-1 DO BEGIN
             ;tm                        = SIN(!Dpi*cnf.base*tet[i_tet]/wav*COS(phi_meas-phi[i_phi]))^2*ctr[i_ctr] ; old approach
             vis                       = (v1 + ctr[i_ctr]*v2*EXP(-DCOMPLEX(0,2*!dpi*cnf.base*tet[i_tet]/wav*COS(phi[i_phi]-phi_base))))/(1+ctr[i_ctr])
             null                      = (1.-ABS(vis))/(1.+ABS(vis))
             r_chi2[i_tet,i_phi,i_ctr] = TOTAL((null_sci-null)^2/null_etot^2)/(N_ELEMENTS(null_sci)-2)
          ENDFOR
      ENDFOR
  ENDFOR
  idx_min = WHERE(r_chi2 EQ MIN(r_chi2), n_min)
  idx_chi = ARRAY_INDICES(r_chi2, idx_min)
  ;PRINT, 'Best chi2: ', MIN(r_chi2)
  ;PRINT, 'Best ang. sep.:', tet[idx_chi[0]] / prm.m2r
  ;PRINT, 'Best PA      :', phi[idx_chi[1]] / !DTOR
  ;PRINT, 'Best contrast:', ctr[idx_chi[2]]
  ctr_opt = ctr[idx_chi[2]]
  
  ; Compute error bar
  r_chi2 = DBLARR(n_tet,n_phi,n_ctr) 
  null_sci_up = null_sci + null_etot
  FOR i_tet = 0, n_tet-1 DO BEGIN
    FOR i_phi = 0, n_phi-1 DO BEGIN
      FOR i_ctr = 0, n_ctr-1 DO BEGIN
        ;tm                        = SIN(!Dpi*cnf.base*tet[i_tet]/wav*COS(phi_meas-phi[i_phi]))^2*ctr[i_ctr] ; old approach
        vis                       = (v1 + ctr[i_ctr]*v2*EXP(-DCOMPLEX(0,2*!dpi*cnf.base*tet[i_tet]/wav*COS(phi[i_phi]-phi_base))))/(1+ctr[i_ctr])
        null                      = (1.-ABS(vis))/(1.+ABS(vis))
        r_chi2[i_tet,i_phi,i_ctr] = TOTAL((null_sci_up-null)^2/null_etot^2)/(N_ELEMENTS(null_sci)-2)
      ENDFOR
    ENDFOR
  ENDFOR
  idx_min = WHERE(r_chi2 EQ MIN(r_chi2), n_min)
  idx_chi = ARRAY_INDICES(r_chi2, idx_min)
  ctr_up = ctr[idx_chi[2]]
  
  r_chi2 = DBLARR(n_tet,n_phi,n_ctr)
  null_sci_do = null_sci - null_etot
  FOR i_tet = 0, n_tet-1 DO BEGIN
    FOR i_phi = 0, n_phi-1 DO BEGIN
      FOR i_ctr = 0, n_ctr-1 DO BEGIN
        ;tm                        = SIN(!Dpi*cnf.base*tet[i_tet]/wav*COS(phi_meas-phi[i_phi]))^2*ctr[i_ctr] ; old approach
        vis                       = (v1 + ctr[i_ctr]*v2*EXP(-DCOMPLEX(0,2*!dpi*cnf.base*tet[i_tet]/wav*COS(phi[i_phi]-phi_base))))/(1+ctr[i_ctr])
        null                      = (1.-ABS(vis))/(1.+ABS(vis))
        r_chi2[i_tet,i_phi,i_ctr] = TOTAL((null_sci_do-null)^2/null_etot^2)/(N_ELEMENTS(null_sci)-2)
      ENDFOR
    ENDFOR
  ENDFOR
  idx_min = WHERE(r_chi2 EQ MIN(r_chi2), n_min)
  idx_chi = ARRAY_INDICES(r_chi2, idx_min)
  ctr_do = ctr[idx_chi[2]]
  
  ; Error bar
  PRINT, 'Error on contrast', 0.5*ABS(ctr_do-ctr_up)
  
  ; Use theoretical values
  ctr_opt = 0.0355
  ctr_do  = ctr_opt - 0.0022
  ctr_up  = ctr_opt + 0.0022
    
  ; Plot the calibrated nulls
  IF KEYWORD_SET(plot) THEN BEGIN
    result_path = pth.result_path + pth.sep + 'fit' + pth.sep
    IF NOT FILE_TEST(result_path) THEN FILE_MKDIR, result_path 
    plotname = result_path + 'UT' + date + '_' + STRTRIM(tgt_uniq[i_tgt],2) + '_calib.eps'
    IF KEYWORD_SET(ps) THEN BEGIN & PREP_PS, /BOLD & DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=17.8, YSIZE=13.2, /TIMES & ENDIF $
                       ELSE WINDOW, !D.WINDOW+1, TITLE='Calibrated null'
    LOADCT, 12, /SILENT
    pa_range = ABS(MAX(sci_para)-MIN(sci_para))
    xrange   = [MIN(sci_para)-0.1*pa_range,MAX(sci_para)+0.1*pa_range]
    yrange   = [-0.5,4.]                             
    PLOT, [0.], [0.], XTITLE='Parallactic angle [deg]', YTITLE='Calibrated null [%]', TITLE=' ', XSTYLE=1, YSTYLE=1, XRANGE=xrange, YRANGE=yrange, COLOR=0
    ; Overplot the estimated null values for calibrators measurements
    LOADCT, 13, /SILENT
    USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.7, THICK=3.0, /FILL
    OPLOT, [sci_para], 100.*[null_sci], PSYM=8, COLOR=255
    ERRPLOT, [sci_para], 100.*([null_sci]-[null_etot]), 100.*([null_sci]+[null_etot]), COLOR=255
    LOADCT, 12, /SILENT
    ; Overplot the stellar photosphere
    n_plot = 1000 & n_surf = 1000
    plot_gleak = DBLARR(n_plot,n_surf)                                                 ; Prepare surfaces for 3-sigma filled contour plot
    ysurf      = DINDGEN(n_surf)/(n_surf-1D0) * (yrange[1]-yrange[0]) + yrange[0]
    plot_pa    = DINDGEN(n_plot)/(n_plot-1D0) * (xrange[1]-xrange[0]) + xrange[0]
    FOR iplot = 0, n_plot-1D0 DO BEGIN
      plot_gleak[iplot,*] = (ABS(ysurf-100.*null_wb) LT 3.*100.*null_err)
    ENDFOR
    OPLOT, plot_pa, REPLICATE(100.*null_wb, n_plot), COLOR=0., LINESTYLE=2
    CONTOUR, plot_gleak, plot_pa, ysurf, /CELL_FILL, /OVERPLOT, C_COLORS=[100], LEVELS=[1]
    ; Overplot the best model  
    phi_plot = 0.5*!dpi+plot_pa*!DTOR
    idx_neg  = WHERE(plot_pa GT 0, n_neg)
    IF n_neg GT 0 THEN phi_plot[idx_neg] = -0.5*!dpi+plot_pa[idx_neg]*!DTOR
    vis      = (v1 + ctr_opt*v2*EXP(-DCOMPLEX(0,2*!dpi*cnf.base*tet[idx_chi[0]]/wav*COS(phi[idx_chi[1]]-phi_plot))))/(1+ctr_opt)
    bin_fit  = (1.-ABS(vis))/(1.+ABS(vis))    
    OPLOT, plot_pa, 100.*bin_fit, LINESTYLE=0, COLOR=100
    vis      = (v1 + ctr_up*v2*EXP(-DCOMPLEX(0,2*!dpi*cnf.base*tet[idx_chi[0]]/wav*COS(phi[idx_chi[1]]-phi_plot))))/(1+ctr_up)
    bin_fit  = (1.-ABS(vis))/(1.+ABS(vis))
    OPLOT, plot_pa, 100.*bin_fit, LINESTYLE=1, COLOR=100
    vis      = (v1 + ctr_do*v2*EXP(-DCOMPLEX(0,2*!dpi*cnf.base*tet[idx_chi[0]]/wav*COS(phi[idx_chi[1]]-phi_plot))))/(1+ctr_do)
    bin_fit  = (1.-ABS(vis))/(1.+ABS(vis))
    OPLOT, plot_pa, 100.*bin_fit, LINESTYLE=1, COLOR=100
    ; Print info to plot
    ;xpos0 = 0.04*(xrange[1]-xrange[0]) + xrange[0] 
    ;ypos0 = 0.95*(yrange[1]-yrange[0]) + yrange[0] 
    ;XYOUTS, xpos0, ypos0*0.98, 'Angular separation = ' + STRING(tet[idx_chi[1]]/m2r,FORMAT='(I0)') + ' mas'
    ;XYOUTS, xpos0, ypos0*0.92, 'Contrast = ' + STRING(ctr[idx_chi[2]],FORMAT='(F5.2)') + '%'
    ;XYOUTS, xpos0, ypos0*0.86, '!4v!3!dr!u2!n = ' + STRING(MIN(r_chi2),FORMAT='(F6.2)')
    LOADCT, 0, /SILENT
    IF KEYWORD_SET(ps) THEN BEGIN & DEVICE, /CLOSE & END_PS & ENDIF
  ENDIF
ENDFOR

IF KEYWORD_SET(log_file) THEN BEGIN
  PRINTF,lun, ' '
  PRINTF,lun, ' '
  CLOSE, lun
  FREE_LUN, lun
ENDIF 
END