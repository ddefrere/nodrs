; +
; NAME: NULL_BOOTSTRAP
;   
; PURPOSE:
;   This function bootstraps and calibrates null data given as input. It was used to extrapolate from 0.05% to 0.04% for the ORR but
;   has not been updated since then (so, probably still ok but not up-to-date regarding the calibration strategy).
;
; INPUTS:
;   date           :  The date to be calibrated
;   cfg_file       :  String with the name of the config file with the reduction parameters
;
; KEYWORDS
;   LOG_FILE       :
;   NO_INSET       : set this keyword to not display the bottom inset with the background null in the TF plot
;   REMOVE_OB      : set this keyword to the OB number to be removed from the data
;   PLOT           : set this keyword to have plot the results
;   VERSION        : set this keyword to fit an older version of the L1 summary file
;
; MODIFICATION HISTORY:
;   Version 1.0,  16-JAN-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  14-APR-2014, DD: removed high reduced chi2

PRO NULL_BOOTSTRAP, date, cfg_file, LOG_FILE=log_file, NO_INSET=no_inset, REMOVE_OB=remove_ob, INFO=info, PLOT=plot, VERSION=version

drs_version = 2.6
drs_date    = '13-MAR-2015'
  
; Start actual code
ON_ERROR, 0

; Keyword sanity check
IF NOT KEYWORD_SET(info)       THEN info       = 1
IF NOT KEYWORD_SET(no_inset)   THEN no_inset   = 0

; Define running and plotting paramaters
n_lam     = 20                     ;Number of wavelength bins within the bandwidth
n_time    = 1000                   ;Number of points in the interpolated TF curve
charsize  = 1.2
charthick = 3.0

; DEFINE GLOBAL VARIABLES
; ***********************

; Astronomical and physical constants
GET_PRM, prm

; Read config file with reduction parameters
GET_DRS, drs, 'nodrs/cfg/' + cfg_file
IF drs.null_weight GT 0 THEN wei_exp = 2 ELSE wei_exp = 0                      ; O to turn off weights (2 to turn on)

; Obtain the definition of the configuration
GET_CNF, cnf, INSTRUM=drs.instrum

; Recover the IDL running path
DECLARE_PATH, pth, INSTRUM=drs.instrum


; COMPUTE FILE PATH AND READ DATA
; *******************************

; Derive long version of the date
date_lng    = '20' + STRMID(date, 0, 2) + '-' + STRMID(date, 2, 2) + '-' + STRMID(date, 4, 2)
l1fits_path = pth.l1fits_path + date_lng + pth.sep
IF NOT FILE_TEST(l1fits_path) THEN MESSAGE, 'No L1 data for ' + date_lng

; Output directory
l2_dir = pth.l2fits_path + pth.sep + date_lng + pth.sep
IF NOT FILE_TEST(l2_dir) THEN FILE_MKDIR, l2_dir

; Init log file
IF KEYWORD_SET(log_file) THEN BEGIN
  log_file =  l2_dir + date_lng + '.txt'
  OPENW, lun, log_file, /GET_LUN, WIDTH=800, /APPEND
  PRINTF,lun, ' '
  PRINTF,lun, 'NULL_CALIB.pro version ' + STRING(drs_version, FORMAT='(F3.1)') +  ' -- ' + drs_date + ' -- Denis Defrère - Steward Observatory (ddefrere@email.arizona.edu)'
  PRINTF,lun, 'Proceeding at your own risk now, on '+ SYSTIME()
  PRINTF,lun, ' '
ENDIF ELSE lun = -1

; Search for L1 file and read it
IF KEYWORD_SET(VERSION) THEN vers = '_v' + STRING(version, FORMAT='(I0)') ELSE vers = ''
l1files = FILE_SEARCH(l1fits_path,'UT' + date_lng + vers + '.fits')
l1data  = MRDFITS(l1files[0], 1, hdr_col, /SILENT)   ; Read null data
detect  = MRDFITS(l1files[0], 2, /SILENT)            ; Read detector information data
phasec  = MRDFITS(l1files[0], 5, /SILENT)            ; Read PHASECam data
weather = MRDFITS(l1files[0], 6, /SILENT)            ; Read weather data
header  = HEADFITS(l1files[0], /SILENT)              ; Read main header
 
; If REMOVE_OB is set, remove the corresponding OBs from the data
n_data = N_ELEMENTS(l1data.mjd_obs)
IF KEYWORD_SET(REMOVE_OB) THEN BEGIN
  idx_bad = VALUE_LOCATE(remove_ob, l1data.ob_id) 
  idx_ok  = [0, WHERE(idx_bad EQ SHIFT(idx_bad, 1), n_data)]
  n_data  += 1
  IF n_data GT 0 THEN BEGIN
    l1data  = l1data[idx_ok]
    phasec  = phasec[idx_ok]
    detect  = detect[idx_ok]
    weather = weather[idx_ok]  
  ENDIF
ENDIF
;idx_pro = WHERE(STRMATCH(STRCOMPRESS(l1data.objname,/REMOVE_ALL), 'procyon') EQ 1)
;l1data[idx_pro].datatype = 'CAL'
;idx_sci = VALUE_LOCATE(l1data.ob_id,[1,3,5,7]) 
;l1data[idx_sci].datatype = 'SCI'

; Remove high chi2
n_data = N_ELEMENTS(l1data.mjd_obs)
IF drs.chi2_lim AND FXPAR(header, 'NUL_MODE') EQ 2 THEN BEGIN
  idx_ok = WHERE(l1data.NSC_CHI2 LT drs.chi2_lim, n_data)
  IF n_data GT 0 THEN BEGIN
    l1data  = l1data[idx_ok]
    phasec  = phasec[idx_ok]
    detect  = detect[idx_ok]
    weather = weather[idx_ok]
  ENDIF
ENDIF

; Make sure data is chronological
idx_srt  = SORT(l1data.mjd_obs)
l1data   = l1data[idx_srt]
phasec   = phasec[idx_srt]
detect   = detect[idx_srt]
weather  = weather[idx_srt]

; Assign pointing ID
pt_id = INTARR(n_data)
FOR i = 1, n_data-1 DO IF STRMATCH(l1data[i].objname, l1data[i-1].objname, /FOLD_CASE) EQ 0 THEN pt_id[i] = pt_id[i-1] + 1 ELSE pt_id[i] = pt_id[i-1]
n_pt  = MAX(pt_id) + 1

; SUBTRACT GEOMETRIC NULL FROM CALIBRATOR NULL
; ********************************************

; Find index of calibrator data
idx_cal = WHERE(STRMATCH(l1data.flag, 'CAL') EQ 1, n_caldat)
IF n_caldat LE 0 THEN BEGIN
  PRINT, 'No calibrator data. Use science data as calibrator!'
  IF lun GT 0 THEN PRINTF, lun, 'No calibrator data. Use science data as calibrator!'
  idx_cal = WHERE(STRMATCH(l1data.flag, 'SCI') EQ 1, n_caldat)
  IF n_caldat LE 0 THEN MESSAGE, 'No data available for null calibration'
ENDIF

; Read relevant header information
aper_rad = FXPAR(header, 'APER_RAD')
bck_irad = FXPAR(header, 'BCK_IRAD')
bck_orad = FXPAR(header, 'BCK_ORAD')
bck_mode = FXPAR(header, 'BCK_MODE')
null_cor = FXPAR(header, 'NULL_COR')
nsc_cube = FXPAR(header, 'NSC_CUBE') 
n_btstrp = FXPAR(header, 'N_BTSTRP')

; Print some info to screen
IF info GT 0 THEN BEGIN
  PRINT, '==== REDUCTION INFO ===='
  PRINT, 'Aperture radius : ' + STRING(aper_rad, FORMAT='(I0)')
  PRINT, 'Background inner radius : ' + STRING(bck_irad, FORMAT='(I0)')
  PRINT, 'Background outer radius : ' + STRING(bck_orad, FORMAT='(I0)')
  PRINT, 'Background mode         : ' + STRING(bck_mode, FORMAT='(I0)')
  PRINT, 'Null correction mode    : ' + STRING(null_cor, FORMAT='(I0)')
  PRINT, 'NSC cube size           : ' + nsc_cube
  PRINT, 'Number of bootstrap     : ' + STRING(n_btstrp, FORMAT='(I0)')
  PRINT, 'Reduced chi2 limit      : ' + STRING(drs.chi2_lim, FORMAT='(F3.1)')
  IF lun GT 0 THEN BEGIN
    PRINTF, lun,  '==== REDUCTION INFO ===='
    PRINTF, lun, 'Aperture radius : ' + STRING(aper_rad, FORMAT='(I0)')
    PRINTF, lun, 'Background inner radius : ' + STRING(bck_irad, FORMAT='(I0)')
    PRINTF, lun, 'Background outer radius : ' + STRING(bck_orad, FORMAT='(I0)')
    PRINTF, lun, 'Background mode         : ' + STRING(bck_mode, FORMAT='(I0)')
    PRINTF, lun, 'Null correction mode    : ' + STRING(null_cor, FORMAT='(I0)')
    PRINTF, lun, 'NSC cube size           : ' + nsc_cube
    PRINTF, lun, 'Number of bootstrap     : ' + STRING(n_btstrp, FORMAT='(I0)')
    PRINTF, lun, 'Reduced chi2 limit      : ' + STRING(drs.chi2_lim, FORMAT='(F3.1)')
  ENDIF
ENDIF

; Extract calibrator data
caldata    = l1data[idx_cal]
cal_name   = STRCOMPRESS(caldata.objname, /REMOVE_ALL)
cal_time   = caldata.mjd_obs
cal_wav    = caldata.wav_eff
cal_bdw    = caldata.bandwidth
cal_null   = caldata.null_meas
cal_nerr   = caldata.null_err
cal_bias   = caldata.bckg_bias
cal_ebias  = caldata.bckg_ebias
cal_fro    = caldata.nfr_ob
cal_seeing = caldata.seeing  

; Extract complementary data
pt_id_cal  = pt_id[idx_cal]
cal_smttau = weather[idx_cal].smttau
cal_wind   = weather[idx_cal].windspd
cal_fpcp   = phasec[idx_cal].fpc_pists*2.2/360.
cal_phavg  = phasec[idx_cal].pcphmean*2.2/360.
cal_phstd  = phasec[idx_cal].pcphstd*2.2/360.
cal_phavg_err = phasec[idx_cal].pcphmean_err*2.2/360.
cal_phstd_err = phasec[idx_cal].pcphstd_err*2.2/360.

; Define arrays
tf_wb       = DBLARR(n_caldat)
tf_cal      = tf_wb
tf_wb_estat = tf_wb
tf_wb_esyst = tf_wb
tf_wb_etot  = tf_wb
cal_flx     = tf_wb

; Derive number of unique wavelength
wav_uniq = cal_wav[UNIQ(cal_wav, SORT(cal_wav))]
n_wav    = N_ELEMENTS(wav_uniq)

; Loop over the wavelength
FOR i_wav = 0, n_wav-1 DO BEGIN  
  ; Extract data of this wavelength 
  idx_wav  = WHERE(cal_wav EQ wav_uniq[i_wav])
  calname  = cal_name[idx_wav]
  calbdw   = cal_bdw[idx_wav]
  
  ; Check whether there is only one bandwidth for this wavelength
  bdw_uniq = calbdw[UNIQ(calbdw, SORT(calbdw))]
  IF N_ELEMENTS(bdw_uniq) GT 1 THEN MESSAGE, 'Found two different bandwidths for the same wavelenth'
         
  ; Compute wavelength range
  lambda_chnl = wav_uniq[i_wav] - 0.5*bdw_uniq[0] + bdw_uniq[0] * (1.+1./n_lam) * DINDGEN(n_lam)/(n_lam-1)
  
  ; Extract LBTI transmission profile for this wavalength (N' only now)
  IF wav_uniq[i_wav] GE 1.1D-5 AND wav_uniq[i_wav] LE 1.2D-5 THEN BEGIN
    READ_TABLE, 'nodrs/input/n-band_thruput.txt', lam_tmp, thruput, FIRST=4, SEPARATOR='TAB'
    trans = INTERPOL(thruput, lam_tmp, lambda_chnl)
  ENDIF ELSE trans = 1
  
  ; Compute effective wavelength
  lam_eff  = TOTAL(trans*lambda_chnl)/TOTAL(trans)

  ; Compute unique calibrators for this wavelength
  nam_uniq = calname[UNIQ(calname, SORT(calname))]
  n_cal    = N_ELEMENTS(nam_uniq)
  
  ; Loop over calibrator data
  FOR i_cal = 0, n_cal - 1 DO BEGIN
    ; Extract data for this calibrator
    idx_cal = idx_wav[WHERE(calname EQ nam_uniq[i_cal])]
    
    ; Retrieve target information
    star    = GET_TGT(nam_uniq[i_cal])
    th_mean = star.ldm
    th_min  = star.ldm - star.lde
    th_max  = star.ldm + star.lde
    uld     = star.uld
    
    ; Generate stellar spectrum for broadband null computation
    star_flux        = BLACKBODY(star.temp, lambda_chnl, STANDARD=3)                       ; flux in photons/s/m²/Hz/sr
    flux_int         = TOTAL(trans*star_flux*!Dpi*(0.5*th_mean*prm.m2r)^2)/TOTAL(trans)    ; Integrated flux in ph/s/m2/Hz
    cal_flx[idx_cal] = flux_int*prm.h*prm.c*prm.w2jy/lam_eff                               ; Integrated flux in Jansky
    
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
    
    ; Derive the wide-band instrumental transfer function and its error bar
    tf_wb[idx_cal]       = cal_null[idx_cal] - null_wb                         ; correct for geometric leakage 
    tf_wb_estat[idx_cal] = cal_nerr[idx_cal]                                   ; statistical error bar
    tf_wb_esyst[idx_cal] = null_err                                            ; systematic error bar
    tf_wb_etot[idx_cal]  = SQRT(tf_wb_estat[idx_cal]^2+tf_wb_esyst[idx_cal]^2) ; total error bar
  ENDFOR
ENDFOR

; Compute non-weighted statistical error on TF and background (nearby empty region)
rms_stat_null = SQRT(1./TOTAL(1./tf_wb_estat^2))
rms_stat_bckg = SQRT(1./TOTAL(1./cal_ebias^2))
; Compute non-weighted calibrator null and background dispersion
AVGSDV, tf_wb, avg_tot_null, rms_tot_null, KAPPA=5
AVGSDV, cal_bias, avg_tot_bckg, rms_tot_bckg, KAPPA=5
; Compute weighted calibrator null dispersion
wei = 1./tf_wb_etot^wei_exp & AVGSDV, tf_wb, avg_tot_null_we, rms_tot_null_we, KAPPA=5, WEIGHT=wei
rms_tot_null_we_mean = rms_tot_null_we*SQRT(TOTAL((wei/TOTAL(wei))^2))
; Same for background region
wei = 1./cal_ebias^wei_exp & AVGSDV, cal_bias, avg_tot_bckg_we, rms_tot_bckg_we, KAPPA=5, WEIGHT=wei
rms_tot_bckg_we_mean = rms_tot_bckg_we*SQRT(TOTAL((wei/TOTAL(wei))^2))

; Search for science target fits files in the input folder
idx_sci = WHERE(STRMATCH(l1data.flag, 'SCI') EQ 1, n_scidat)
IF n_scidat GT 0 THEN scidata = l1data[idx_sci]

; Extract main science data in chronological order
sci_time   = scidata.mjd_obs
sxi_wav    = scidata.wav_eff
sci_null   = scidata.null_meas
sci_nerr   = scidata.null_err

; Define pointing ID
pt_id = INTARR(n_caldat)
FOR i = 1, n_caldat-1 DO IF STRMATCH(cal_name[i], cal_name[i-1], /FOLD_CASE) EQ 0 THEN pt_id[i] = pt_id[i-1] + 1 ELSE pt_id[i] = pt_id[i-1]
n_pt  = MAX(pt_id) + 1

; Add pointing
n_seq  = 8
n_time = 12
n_btsp = 100
rms_out = DBLARR(n_btsp, n_time)
rms_pt  = DBLARR(n_time)
rms_ept = rms_pt 
FOR i_t = 0, n_time-1 DO BEGIN
  n_pta = i_t + 1
  n_pt0 = 1 + n_pta
  FOR i=0, n_btsp-1 DO BEGIN
    
    ; Draw first random pointing (first cal pointing)
    idx_rndm     = FLOOR(RANDOMU(seed, n_seq)*n_caldat)
    pt_id0       = [INTARR(n_seq)]
    tf_wb0       = [tf_wb[idx_rndm]]
    tf_wb_etot0  = [tf_wb_etot[idx_rndm]]
    tf_wb_estat0 = [tf_wb_estat[idx_rndm]]
    idx_rndm     = FLOOR(RANDOMU(seed, n_seq)*n_scidat)
    sci_nerr0    = [sci_nerr[idx_rndm]]
    
    ; Loop over CAL pointings to add
    FOR j=0, n_pta-1 DO BEGIN
      idx_rndm    = FLOOR(RANDOMU(seed, n_seq)*n_caldat)
      pt_id0      = [pt_id0, 1+j+INTARR(n_seq)] 
      tf_wb0      = [tf_wb0, tf_wb[idx_rndm]]
      tf_wb_etot0 = [tf_wb_etot0, tf_wb_etot[idx_rndm]]
      tf_wb_estat0 = [tf_wb_estat0, tf_wb_estat[idx_rndm]]
    ENDFOR
       
    ; Loop over SCI pointings to add
    FOR j=0, n_pta-1 DO BEGIN
      idx_rndm  = FLOOR(RANDOMU(seed, n_seq)*n_scidat)
      sci_nerr0 = [sci_nerr0, sci_nerr[idx_rndm]]
    ENDFOR
    
    ; Compute STAT error
    rms_stat_null = SQRT(1./TOTAL(1./tf_wb_estat0^2))
    sci_stat_err  = SQRT(1./TOTAL(1./sci_nerr0^2))
    
    ; Compute weigthed instrumental null per pointing
    null_avg_pt = FLTARR(n_pt0) & null_err_pt = FLTARR(n_pt0)
    bias_avg_pt = FLTARR(n_pt0) & bias_err_pt = FLTARR(n_pt0)
    FOR i_pt = 0, n_pt0-1 DO BEGIN
      ; Derive OB of this pointing
      idx_pt = WHERE(pt_id0 EQ i_pt, n_ob)
      ; Compute weighted calibrator null dispersion
      wei = 1./tf_wb_etot0[idx_pt]^wei_exp & AVGSDV, tf_wb0[idx_pt], avg_null_pt_we, rms_null_pt_we, KAPPA=5, WEIGHT=wei
      null_avg_pt[i_pt] = avg_null_pt_we & null_err_pt[i_pt] = rms_null_pt_we*SQRT(TOTAL((wei/TOTAL(wei))^2))
    ENDFOR 
    
    ; Compute error between pointings (used RMS as systematic error)
    wei_null = 1D/null_err_pt^wei_exp & AVGSDV, null_avg_pt, avg_sys, rms_sys, WEIGHT=wei_null
    rms_mean = rms_sys*SQRT(TOTAL((wei_null/TOTAL(wei_null))^2))
    
    ; Compute error per OB
    wei = 1./tf_wb_etot0^wei_exp & AVGSDV, tf_wb0, avg_null_ob_we, rms_null_ob_we, KAPPA=5, WEIGHT=wei
    rms_mean2 = rms_null_ob_we*SQRT(TOTAL((wei/TOTAL(wei))^2))
    
    ; Take the worst of the 3
    rms_sys = rms_mean > rms_mean2 > rms_stat_null
    
    ; Error bar
    rms_out[i,i_t] = SQRT(rms_sys^2+sci_stat_err^2) 
  ENDFOR
  AVGSDV, rms_out[*,i_t], avg, rms
  rms_pt[i_t]  = avg*1D2
  rms_ept[i_t] = rms*1D2
ENDFOR

; Time vector
timppt = 25/60.
time   = [3*timppt + DINDGEN(n_time)*timppt*2]
time_fit = FINDGEN(100)/99*(MAX(time)+0.5)
pt     = time/timppt

; Fit polynomial
coeff = POLY_FIT(1./SQRT(time[INDGEN(3)]), rms_pt[INDGEN(3)], 1, MEASURE_ERRORS=rms_ept)
fit   = (coeff[1]/SQRT(time_fit)+coeff[0])

tf_path  = pth.result_path + 'TF' + pth.sep + date_lng + pth.sep
plotname =  tf_path + date_lng + '_BOOTSTRAP_APER' + STRING(aper_rad, FORMAT='(I0)') + '.eps'
PREP_PS, /BOLD & LOADCT, 12, /SILENT
DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=18.8, YSIZE=14.7
xrange = [1.0,10+0.5]
yrange = [0.,0.1]
PLOT, [0], [0], XTITLE='Sequence time [h]', YTITLE='Calibrated null uncertainty [%]', TITLE=title_day, XSTYLE=9, YSTYLE=1, $
      XRANGE=xrange, YRANGE=yrange, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize
AXIS, XAXIS=1, XTITLE='Number of pointings', XRANGE=[MIN(pt),MAX(pt)], XSTYLE=1
; Overplot the estimated null values for calibrators measurements
OPLOT, time_fit, fit, COLOR=0
OPLOT, xrange, [0.04,0.04], COLOR=0, LINESTYLE=1
; Overplot measured data
i_m = 2
OPLOT, [6*timppt], [0.053], PSYM=4, COLOR=25
;ERRPLOT, [time[i_m]]-0.5, [rms_pt[i_m]+rms_ept[i_m]], [rms_pt[i_m]-rms_ept[i_m]], COLOR=25
; Plot data subsets
idx = WHERE(time LT time[i_m])
OPLOT, time[idx], rms_pt[idx], PSYM=4, COLOR=200
ERRPLOT, time[idx], rms_pt[idx]+rms_ept[idx], rms_pt[idx]-rms_ept[idx], COLOR=200
; Plot data subsets
idx = WHERE(time GE time[i_m])
OPLOT, time[idx], rms_pt[idx], PSYM=4, COLOR=100
ERRPLOT, time[idx], rms_pt[idx]+rms_ept[idx], rms_pt[idx]-rms_ept[idx], COLOR=100
LEGEND, ['Data subsets','Measured value','Bootstrapped data'], COLOR=[200,25,100], PSYM=[4,4,4], THICK=1.5, PTHICK=5.5, CHARTHICK=charthick, CHARSIZE=charsize, /TOP_LEGEND, /RIGHT_LEGEND, /FILL
DEVICE, /CLOSE & END_PS

PRINT, 'Bootstrapped error and rms on error', 1D2*avg, 1D2*rms
END
