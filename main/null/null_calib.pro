;+
; NAME: NULL_CALIB
; 
; PURPOSE:
;   This function calibrates null data given as input.
;
; INPUTS:
;   date           :  The date to be calibrated
;   cfg_file       :  String with the name of the config file with the reduction parameters
;
; KEYWORDS
;   CALPOB         : set this keyword to calibrate per OB (per pointing by default)
;   LOG_FILE       :
;   NO_INSET       : set this keyword to not display the bottom inset with the background null in the TF plot
;   REMOVE_OB      : set this keyword to the OB number to be removed from the data
;   REMOVE_ID      : set this keyword to the file ID number to be removed from the data
;   RUNBIAS        : set this keyword to calibrate the background bias 
;   PLOT           : set this keyword to have plot the results
;   VERSION        : set this keyword to fit an older version of the L1 summary file
;   
; LIMITATION
;   Transmission profile of the filter only included for N' and 8.7um (used to compute the effective wavelength)
;
; MODIFICATION HISTORY:
;   Version 1.0,  16-JAN-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  15-MAY-2013, DD: adapted to new format (now binary fits) of input files
;   Version 1.2,  15-OCT-2013, DD: updated to the new format which consists in having all null measurements of a night in the same file
;   Version 1.3,  25-JAN-2014, DD: updated output file for a first example release to RMG and added keyword SPLIT_TF
;   Version 1.4,  25-FEB-2014, DD: added keyword REMOVE_OB
;   Version 1.5,  31-MAR-2014, DD: added keyword SPLIT_NOD
;   Version 1.6,  16-APR-2014, DD: added keyword version
;   Version 1.7,  25-MAY-2014, DD: added background bias plot
;   Version 1.8,  10-AUG-2014, DD: implemented filter transmission
;   Version 1.9,  13-AUG-2014, DD: added photometry plot
;   Version 2.0,  18-NOV-2014, DD: added aperture radius in plot filenames
;   Version 2.1,  22-NOV-2014, DD: improved speed and robustness
;   Version 2.2,  29-NOV-2014, DD: added keyword NO_INSET
;   Version 2.3,  02-DEC-2014, DD: now use AVGSDV for weighted average
;   Version 2.4,  07-JAN-2014, DD: added new diagnostic plots
;   Version 2.5,  23-JAN-2015, DD: added error on phasecam telemetry to output data (+corrected minor bug)
;   Version 2.6,  13-MAR-2015, DD: added log file and cfg_file input (removed associated keywords)
;   Version 2.7,  04-APR-2015, DD: improved output file header for major archiv release
;   Version 2.8,  10-APR-2015, DD: added unweighted results to output plot
;   Version 2.9,  12-APR-2015, DD: added chi2 limit
;   Version 3.0,  16-APR-2015, DD: added the possibility to split the TF at user-defined hours
;   Version 3.1,  13-MAY-2015, DD: added quality flag
;   Version 3.2,  26-MAY-2015, DD: now properly read the PID
;   Version 3.3,  03-AUG-2015, DD: updated for pipeline paper
;   Version 3.4,  06-AUG-2015, DD: major update (improved speed and cleaned code)
;   Version 3.5,  02-DEC-2015, DD: adapted for new column format of APER_RAD<M BCK_IRAD, and BCK_ORAD
;   Version 3.5,  06-DEC-2015, DD: corrected implementation of hour angle
;   Version 3.6,  12-JAN-2016, DD: updated axis labels in output plots (astrophysical null => null depth)
;   Version 3.7,  05-FEB-2016, DD: updated to write two L2 files (one with the results per pointing and one per OB)
;   Version 3.8,  15-APR-2016, DD: corrected minor bug with pointing ID
;   Version 3.9,  21-APR-2016, DD: added parameter NOD_COR
;   Version 4.0,  27-APR-2016, DD: now create a log file differnet for each L1 file version (and also print the null per OB)
;   Version 4.1,  12-MAY-2016, DD: added keyword CALPOB
;   Version 4.2,  06-JUL-2016, DD: improved display
;   Version 4.3,  12-OCT-2016, DD: corrected minor bug when only one OB per pointing
;   Version 4.4,  01-NOV-2016, DD: simplified code + removed bias computation
;   Version 4.5,  05-NOV-2016, DD: added keyword RUNBIAS
;   Version 4.6,  28-NOV-2016, DD: improved backwards compatibility
;   Version 4.7,  04-APR-2017, DD: adapted for new formlaism of GET_TGT
;   Version 4.8,  05-APR-2017, DD: now split TF when the aperture radius changes!
;   Version 4.9,  12-JUN-2017, DD: now split TF when one of the background radius change
;   Version 5.0,  03-AUG-2017, DD: updated to call AVGERR instead of repeating the same code over and over
;   Version 5.1,  29-SEP-2017, DD: implemented NULL_EST and SIG_SCA
;   Version 5.2,  11-OCT-2017, DD: updated to loop through the L1 summary files to find one with the right aperture radius
;   Version 5.3,  29-JUL-2018, DD: added comments to error computation
;   Version 5.4,  31-JUL-2018, DD: implemented unweighted POLY_FIT (based on drs.NULL_EST input parameter)
;   Version 5.5,  05-SEP-2018, DD: added keyword REMOVE_ID
;   Version 5.6,  01-NOV-2018, DD: now properly sort the file ID when more than ten L1 files are found + updated log
;   Version 5.7,  10-AUG-2020, DD: updated the name of the L2 file
;   Version 5.8,  15-SEP-2023, DD: added new filter options
;   Version 5.9,  15-OCT-2023, DD: Updated text for FRA_MODE=2 (i.e., PCA background subtraction)
;   Version 6.0,  25-OCT-2023, DD: Updated UV coordinate file name

PRO NULL_CALIB, date, cfg_file, CALPOB=calpob, LOG_FILE=log_file, NO_INSET=no_inset, REMOVE_ID=remove_id, REMOVE_OB=remove_ob, RUNBIAS=runbias, INFO=info, PLOT=plot, VERSION=version

; Calibration routine version, to be updated manually
drs_version = 5.7
drs_date    = '10-AUG-2020'
  
; Start actual code
ON_ERROR, 0

; Keyword sanity check
IF NOT KEYWORD_SET(info)       THEN info       = 1
IF NOT KEYWORD_SET(no_inset)   THEN no_inset   = 0

; Define running and plotting parameters
n_lam     = 20                     ; Number of wavelength bins within the bandwidth (used for broadband computation)
n_time    = 1000                   ; Number of points in the interpolated TF curve
charsize  = 1.2                    ; Character size in plots  
charthick = 3.0                    ; Character thickness in plots

; DEFINE GLOBAL VARIABLES
; ***********************

; Astronomical and physical constants
GET_PRM, prm

; Obtain the definition of the configuration
GET_CNF, cnf, INSTRUM='nomic'

; Read config file with reduction parameters
GET_DRS, drs, 'nodrs/cfg/' + cfg_file
IF KEYWORD_SET(RUNBIAS) THEN BEGIN
  IF KEYWORD_SET(CALPOB) THEN BEGIN
    id_name      = 'OB-BIAS'
    drs.cal_mode = 1
  ENDIF ELSE id_name = 'PT-BIAS'
ENDIF ELSE BEGIN
  IF KEYWORD_SET(CALPOB) THEN BEGIN
    id_name      = 'OB'
    drs.cal_mode = 1
  ENDIF ELSE id_name = 'PT'
ENDELSE

; Recover the IDL running path
DECLARE_PATH, pth, INSTRUM='nomic'


; COMPUTE FILE PATH AND READ DATA (+ basic heaser information)
; *******************************

; Derive long version of the date
date_lng = '20' + STRMID(date, 0, 2) + '-' + STRMID(date, 2, 2) + '-' + STRMID(date, 4, 2)
IF MAX(drs.null_rad) NE 0 AND NOT STRMATCH(drs.dir_label, '*_APR') THEN drs.dir_label = drs.dir_label + '_APR' ; If aperture radius is set, redefine dir_label
l1fits_path = pth.l1fits_path + date_lng + drs.dir_label + pth.sep
IF NOT FILE_TEST(l1fits_path) THEN BEGIN
  MESSAGE, 'No L1 data for ' + date_lng, /CONTINUE
  RETURN
ENDIF

; Output directory
l2_dir = pth.l2fits_path + pth.sep + date_lng + drs.dir_label + pth.sep
IF NOT FILE_TEST(l2_dir) THEN FILE_MKDIR, l2_dir

; File version
IF KEYWORD_SET(VERSION) THEN vers = '_v' + STRING(version, FORMAT='(I0)') ELSE vers = ''

; Init log file
IF KEYWORD_SET(log_file) THEN BEGIN
  log_file =  l2_dir + date_lng + vers + '.txt'
  OPENW, lun, log_file, /GET_LUN, WIDTH=800, /APPEND
  PRINTF,lun, ' '
  PRINTF,lun, '******************************************************************************************************* '
  PRINTF,lun, ' '
  PRINTF,lun, 'NULL_CALIB.pro version ' + STRING(drs_version, FORMAT='(F3.1)') +  ' -- ' + drs_date + ' -- Denis Defrère - Steward Observatory (denis@lbti.org)'
  PRINTF,lun, 'Data calibrated on ' + SYSTIME()
  PRINTF,lun, ' '
  PRINTF,lun, '******************************************************************************************************* '
  PRINTF,lun, ' '
ENDIF ELSE lun = -1

; Search for L1 summary file and read it
l1files = FILE_SEARCH(l1fits_path,'UT' + date_lng + vers + '.fits', COUNT=n_files)
IF NOT FILE_TEST(l1files[0]) THEN BEGIN
  MESSAGE, 'No L1 data for ' + date_lng, /CONTINUE
  RETURN
ENDIF ELSE BEGIN
  ; Search for the right aperture radius
  header   = HEADFITS(l1files[0], /SILENT)              ; Read main header
  aper_rad = FXPAR(header, 'NULL_RAD', /NOCONTINUE)
  IF aper_rad NE drs.null_rad THEN BEGIN
    ; Find the latest file with the right aperture
    l1files = FILE_SEARCH(l1fits_path,'UT' + date_lng + '_v*.fits', COUNT=n_files)
    m_time  = (FILE_INFO(l1files)).mtime
    l1files = l1files[REVERSE(SORT(m_time))]
    FOR i_f = 0, n_files-1 DO BEGIN
      file_id = i_f
      IF FXPAR(HEADFITS(l1files[file_id], /SILENT), 'NULL_RAD', /NOCONTINUE) EQ drs.null_rad THEN i_f = n_files
    ENDFOR
  ENDIF ELSE file_id = 0 
ENDELSE
l1data  = MRDFITS(l1files[file_id], 1, hdr_col, /SILENT)   ; Read null data
detect  = MRDFITS(l1files[file_id], 2, /SILENT)            ; Read detector information data
phasec  = MRDFITS(l1files[file_id], 5, /SILENT)            ; Read PHASECam data
weather = MRDFITS(l1files[file_id], 6, /SILENT)            ; Read weather data
header  = HEADFITS(l1files[file_id], /SILENT)              ; Read main header
 
; Remove duplicated OBs (temporary fix)
idx_ok = REM_DUP(l1data.ob_id)
IF N_ELEMENTS(idx_ok) GT 0 THEN BEGIN
  l1data  = l1data[idx_ok]
  phasec  = phasec[idx_ok]
  detect  = detect[idx_ok]
  weather = weather[idx_ok]
ENDIF ELSE MESSAGE, 'No OB to calibrate'

; If REMOVE_OB is set, remove the corresponding OBs from the data
n_data = N_ELEMENTS(l1data.mjd_obs)
IF KEYWORD_SET(REMOVE_OB) THEN BEGIN
  MATCH, l1data.ob_id, remove_ob, idx_out, idx_tmp, COUNT=n_out
  IF n_out GT 0 THEN BEGIN
    l1data[idx_out].ob_id = -1
    idx_ok = WHERE(l1data.ob_id NE -1, n_data)
    IF n_data GT 0 THEN BEGIN
      l1data  = l1data[idx_ok]
      phasec  = phasec[idx_ok]
      detect  = detect[idx_ok]
      weather = weather[idx_ok]
    ENDIF ELSE MESSAGE, 'No OB to calibrate'
  ENDIF
ENDIF

; If REMOVEIDB is set, remove the corresponding OBs from the data
n_data = N_ELEMENTS(l1data.mjd_obs)
;l1data.file_id = INDGEN(n_data)
IF KEYWORD_SET(REMOVE_ID) THEN BEGIN
  MATCH, l1data.file_id, remove_id, idx_out, idx_tmp, COUNT=n_out
  IF n_out GT 0 THEN BEGIN
    l1data[idx_out].file_id = -1
    idx_ok = WHERE(l1data.file_id NE -1, n_data)
    IF n_data GT 0 THEN BEGIN
      l1data  = l1data[idx_ok]
      phasec  = phasec[idx_ok]
      detect  = detect[idx_ok]
      weather = weather[idx_ok]
    ENDIF ELSE MESSAGE, 'No OB to calibrate'
  ENDIF
ENDIF

; Remove high chi2 if requested
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

; Make sure data is in chronological order
idx_srt  = SORT(l1data.mjd_obs)
l1data   = l1data[idx_srt]
phasec   = phasec[idx_srt]
detect   = detect[idx_srt]
weather  = weather[idx_srt]

; Read relevant header information
bck_mode = FXPAR(header, 'BCK_MODE', /NOCONTINUE)
bfl_mode = FXPAR(header, 'BFL_MODE', /NOCONTINUE)
file_ver = FXPAR(header, 'FILE_VER', /NOCONTINUE)
n_btstrp = FXPAR(header, 'N_BTSTRP', /NOCONTINUE)
nbin_fac = FXPAR(header, 'NBIN_FAC', /NOCONTINUE)
nsc_bins = FXPAR(header, 'NSC_BINS', /NOCONTINUE)
nsc_cube = FXPAR(header, 'NSC_CUBE', /NOCONTINUE)
nsc_mode = FXPAR(header, 'NSC_MODE', /NOCONTINUE)
nsc_omin = FXPAR(header, 'NSC_OMIN', /NOCONTINUE)
null_cor = FXPAR(header, 'NULL_COR', /NOCONTINUE)
null_lim = FXPAR(header, 'NULL_LIM', /NOCONTINUE)
null_mod = FXPAR(header, 'NUL_MODE', /NOCONTINUE)
aper_rad = FXPAR(header, 'NULL_RAD', /NOCONTINUE)
fra_mode = FXPAR(header, 'FRA_MODE', /NOCONTINUE)
img_mode = FXPAR(header, 'IMG_MODE', /NOCONTINUE)
fit_mode = FXPAR(header, 'FIT_MODE', /NOCONTINUE)
flx_mode = FXPAR(header, 'FLX_MODE', /NOCONTINUE)

; Read main HDU structure into variables
objname   = STRCOMPRESS(STRLOWCASE(l1data.objname), /REMOVE_ALL)
pid       = l1data.pid
ob_id     = l1data.ob_id
flag      = l1data.flag
time      = l1data.mjd_obs
wav       = l1data.wav_eff
bdw       = l1data.bandwidth
nfr_ob    = l1data.nfr_ob
nfr_rej   = l1data.nfr_rej
utc       = l1data.lbt_utc
lst       = l1data.lbt_lst
ra        = l1data.lbt_ra
dec       = l1data.lbt_dec
alt       = l1data.lbt_alt
az        = l1data.lbt_az
para      = l1data.lbt_para

; Null or bias
IF NOT KEYWORD_SET(RUNBIAS) THEN BEGIN
  null      = l1data.null_meas
  null_err  = l1data.null_meas_err 
ENDIF ELSE BEGIN
  null      = l1data.bckg_bias
  null_err  = l1data.bckg_bias_rms/SQRT(nfr_ob)    ; RMS is per frame and we want the error on the mean
  IF MAX(null) EQ 0 AND MAX(null_err) EQ 0 THEN BEGIN
    MESSAGE, 'No background bias information.', /CONTINUE
    RETURN
  ENDIF
ENDELSE

; Extract individual error terms
null_rms = l1data.null_meas_rms/l1data.BCKG_CNOD_RMS
err_phot = ABS(l1data.null_meas_avg)*(l1data.phot_err/l1data.phot_avg)       ; error on NULL due to photometric uncertainty (not included in NSC and correlated between all OBs of the same pointing!)
err_bckg = l1data.bckg_cnod_rms/SQRT(l1data.nfr_rej+l1data.nfr_ob)           ; error on background subtraction from complementary nod (this is only an approximation. Go to the L1 log to see the right results)
err_trm  = SQRT(err_phot^2+err_bckg^2)
err_nsc  = SQRT((l1data.null_meas_err^2-err_trm^2) > 0)                      ; error from NSC not saved but one can comptue it from the total error

; Compute expected error terms
IF bck_mode GE 0 THEN null_err_phot = SQRT(err_bckg^2 + l1data.bckg_ebias^2) ELSE null_err_phot = SQRT(2)*err_bckg
null_err_phase = SQRT(4*(l1data.nsc_phavg)^2*(l1data.nsc_phrms)^2 + 2*(l1data.nsc_phrms)^4)/4

; Backwards compatibility
IF TAG_EXIST(l1data, 'aper_rad') THEN aper  = l1data.aper_rad ELSE aper  = REPLICATE(FXPAR(header, 'APER_RAD', /NOCONTINUE), n_data)
IF TAG_EXIST(l1data, 'bck_irad') THEN birad = l1data.bck_irad ELSE birad = REPLICATE(FXPAR(header, 'BCK_IRAD', /NOCONTINUE), n_data)
IF TAG_EXIST(l1data, 'bck_orad') THEN borad = l1data.bck_orad ELSE borad = REPLICATE(FXPAR(header, 'BCK_ORAD', /NOCONTINUE), n_data)
IF NOT TAG_EXIST(l1data, 'null_offset') THEN STRUCT_ADD_FIELD, l1data, 'NULL_OFFSET', FLTARR(n_data)  

; Compute intensity mismatch
int_err  = 2.0*ABS(l1data.PHOTDX_AVG-l1data.PHOTSX_AVG)/(l1data.PHOTDX_AVG+l1data.PHOTSX_AVG)

; Read relevant detector information
xcen      = detect.xcen
ycen      = detect.ycen
dit       = detect.int_time

; Read relevant weather information
seeing    = weather.seeing
smttau    = weather.smttau
wind      = weather.windspd

; Read PHASECam information
IF TAG_EXIST(phasec, 'plc_wav')  THEN plc_wav = phasec.plc_wav*1D+6 ELSE plc_wav = REPLICATE(2.2, n_data)   ; Backward compatibility
IF TAG_EXIST(phasec, 'dith_per') THEN dith_per  = phasec.dith_per ELSE dith_per = REPLICATE(0, n_data)   ; Backward compatibility
fpcp      = phasec.fpc_pists*plc_wav/360.          ; Convert to um
phavg     = phasec.pcphmean*plc_wav/360.           ; Convert to um
phstd     = phasec.pcphstd*plc_wav/360.            ; Convert to um
phavg_err = phasec.pcphmean_err*plc_wav/360.       ; Convert to um
phstd_err = phasec.pcphstd_err*plc_wav/360.        ; Convert to um


; Compute hour angles for each OB
ha = DBLARR(n_data)
FOR i=0, n_data-1 DO BEGIN
  GET_TGT, objname[i], tgt, DATABASE= pth.input_path + drs.database
  ha[i] = HOUR_ANGLE(time[i], tgt.ra, tgt.dec, LAT=TEN(32,42,05.87), LON=-TEN(109,52,18.87), ALTITUDE=3267D0) * 12./180D0 ; in hours
ENDFOR

; Compute UV coordinates
u_coord = cnf.base/(206265.*wav)*COS(para*!DTOR)
v_coord = cnf.base/(206265.*wav)*SIN(-para*!DTOR)  ; u positive towards East

; Assign pointing ID (successive OB on the same target)
; Ensure backward compatibility (pointing ID was not present in the file before)
IF NOT TAG_EXIST(l1data, 'pt_id') THEN BEGIN
  pt_id = INTARR(n_data)
  FOR i = 1, n_data-1 DO IF STRMATCH(objname[i], objname[i-1], /FOLD_CASE) EQ 0 THEN pt_id[i] = pt_id[i-1] + 1 ELSE pt_id[i] = pt_id[i-1]
ENDIF ELSE pt_id = l1data.pt_id
pt_uniq = pt_id[UNIQ(pt_id, SORT(pt_id))]
n_pt    = N_ELEMENTS(pt_uniq)

; Derive the number of nod positions assuming that each nod is located in a different channel (in the Y direction)
; If nod are uncorrelated and split_nod is not set, force the number of nods to 1
IF drs.nod_cor EQ 1 OR drs.split_nod NE 0 THEN nod_pos = FIX(ycen/cnf.y_chan) ELSE nod_pos = INTARR(n_data)
nod_uniq = nod_pos[UNIQ(nod_pos, SORT(nod_pos))]
n_nod    = N_ELEMENTS(nod_uniq)

; Convert time to hour
DAYCNV, MEAN(time)+2400000.5D0, yr, mn, day, hr
title_day = STRING(yr,FORMAT='(I0)') + '/' + STRING(mn,FORMAT='(I0)') + '/' + STRING(day,FORMAT='(I0)')
JDCNV, yr, mn, day, 0., t0day
t0day       = t0day - 2400000.5D0
time        = time - t0day

; Define arrays
tf_wb       = null              ; Wideband transfer function
tf_wb_estat = null_err          ; Statistical error on wideband TF
tf_wb_etot  = null_err          ; Total error on wideband TF
tf_wb_esyst = FLTARR(n_data)    ; Systematic error on wideband TF
tgt_flx     = FLTARR(n_data)    ; Will contain the target flux in Jy
ob_out      = INTARR(n_data)    ; Will contain OBs to be removed
                                         
; Print some info to screen and log file
IF info GT 0 THEN BEGIN
  PRINT, '==== REDUCTION INFO ===='
  PRINT, ''
  PRINT, 'File version number     : ' + STRING(file_ver, FORMAT='(I0)')
  PRINT, 'Input aperture radius   : ' + STRING(aper_rad, FORMAT='(I0)') + '  (0 for EEID + FWHM)'
  ;PRINT, 'Background inner radius : ' + STRING(bck_irad, FORMAT='(I0)')
  ;PRINT, 'Background outer radius : ' + STRING(bck_orad, FORMAT='(I0)')
  PRINT, 'Image combination mode  : ' + STRING(img_mode, FORMAT='(I0)') + '  (0: for automatic)'
  PRINT, 'Background nod mode     : ' + STRING(bck_mode, FORMAT='(I0)') + '  (0: none, 1 nod pairs, 2 closest n frames, 3 dedicated, 3 chopping)'
  PRINT, 'Background floor mode   : ' + STRING(bfl_mode, FORMAT='(I0)') + '  (0: 5-sigma clipped, 1 median)'
  PRINT, 'Frame mode for flux     : ' + STRING(fra_mode, FORMAT='(I0)') + '  (0: mean-background subtracted, 1: processed raw, 2: pca-background subtracted)'
  PRINT, 'Flux mode               : ' + STRING(flx_mode, FORMAT='(I0)') + '  (0: aper. phot., 1: weigh. aper phot., 2: PSF-fitting)'
  PRINT, 'Centroid mode           : ' + STRING(fit_mode, FORMAT='(I0)') + '  (0: none, 1: CNTRD., 2 GCNTRD, 3: Gaussian, 4: Lorentzan, 5: Moffat)'
  PRINT, 'Null estimator          : ' + STRING(drs.null_est, FORMAT='(I0)')
  PRINT, 'Null error mode         : ' + STRING(drs.err_mode, FORMAT='(I0)')
  PRINT, 'Scatter-rejection sigma : ' + STRING(drs.sig_sca, FORMAT='(I0)')
  PRINT, 'Null correction mode    : ' + STRING(null_cor, FORMAT='(I0)') + "  (1: Denis', 2: Bertrand's method to subtract HF noise)"
  PRINT, 'NSC cube size           : ' + nsc_cube + '  (0 for automatic)'
  PRINT, 'NSC mode                : ' + STRING(nsc_mode, FORMAT='(I0)') + '  (0: mode, 1: %best, 2: NSC)' 
  PRINT, 'Acceptable null range   : ' + null_lim
  PRINT, 'Number of bin factor    : ' + STRING(nbin_fac , FORMAT='(I0)') + '  (Bin size for NSC reduction (0: constant, 1: variable)'
  PRINT, 'Number of bootstrap     : ' + STRING(n_btstrp, FORMAT='(I0)')
  PRINT, 'Calibration method      : ' + STRING(drs.cal_method, FORMAT='(I0)') + '  (0: interpolation nearest neighbors, 2: interpolation on all OBs)'
  PRINT, 'Degree of polynomial    : ' + STRING(drs.polydeg, FORMAT='(I0)') 
  PRINT, 'Calibration mode        : ' + STRING(drs.cal_mode, FORMAT='(I0)') + '  (0: per pointing, 1: per OB)'
  PRINT, 'Reduced chi2 limit      : ' + STRING(drs.chi2_lim, FORMAT='(i0)')
  PRINT, 'TF error mode           : ' + STRING(drs.err_mode, FORMAT='(I0)') + '  (0: unweighted null dispersion, 1: weighted null dispersion, 2: excess variance)'
  PRINT, 'Nod correlation         : ' + STRING(drs.nod_cor, FORMAT='(I0)') + '  (1 if nulls appear correlated per nod within a given pointing)'
  PRINT, 'Split TF                : ' + STRING(drs.split_tf, FORMAT='(I0)') + '  (1 to split the TF when there is a dead time longer than split_time between two null measurements)'
  PRINT, 'Split time              : ' + STRING(drs.split_time, FORMAT='(F3.1)') 
  PRINT, 'Split hour              : ' + STRING(drs.split_hour, FORMAT='(I0)') 
  PRINT, 'Split nod               : ' + STRING(drs.split_nod, FORMAT='(I0)') 
  PRINT, ''
  PRINT, 'Number of OBs           : ' + STRING(n_data, FORMAT='(I03)')
  PRINT, 'Number of pointings     : ' + STRING(n_pt, FORMAT='(I03)')
  ;PRINT, 'PLC wavelength          : ' + STRING(plc_wav, FORMAT='(F5.2)')
  PRINT, ''
  PRINT, ''
 ; Column signification
  PRINT, 'Column signification'
  PRINT, 'A: OB identification number'
  PRINT, 'B: Object name'
  PRINT, 'C: Flag (SCI/CAL)'
  PRINT, 'D: UT time (from LBT)'
  PRINT, 'E: Telescope elevation in degrees (from LBT)'
  PRINT, 'F: Parallactic angle in degrees (from LBT)'
  PRINT, 'G: DIMM seeing [arcsec]'
  PRINT, 'H: PWV [mm]'
  PRINT, 'I: DIT [ms]'
  PRINT, 'J: Dithering period [number of frames]'
  PRINT, 'K: Wavelength [microns]'
  PRINT, 'L: Bandwidths [microns]'
  PRINT, 'M: SNR on photometry (DX and SX)' 
  PRINT, 'N: Intensity mismatch (in %)'  
  PRINT, 'O: Tip/tilt error (DX and SX, in mas)'
  PRINT, 'P: RMS of photometry relative to constructive peak (DX and SX, in %)'  
  PRINT, 'Q: RMS of background relative to constructive peak (in %)'
  PRINT, 'R: Background bias relative to constructive peak (in %)'  
  PRINT, 'S: Null (in %)'
  PRINT, 'T: Null offset between optimum NSC value and bootstrapped or bayesian value (in %)'  
  PRINT, 'U: Null uncertainty due to error on constructive peak (in %)'
  PRINT, 'V: Null uncertainty due to external error on background floor (in %)'
  PRINT, 'W: Null uncertainty from NSC fit (in %)'  
  PRINT, 'X: Total null uncertainty (in %)'  
  PRINT, 'Y: Expected null uncertainty due to photometric errors (in %)'
  PRINT, 'Z: Expected null uncertainty due to phase variation (i.e., sqrt(4mu^2.sig^2 + 2*sig^4)/4, in %)'  
  PRINT, 'a: Best-fit mean PHASE from NSC (in microns)'
  PRINT, 'b: Best-fit PHASE jitter from NSC (in microns)'
  PRINT, 'c: Reduced chi2 from NSC'
  PRINT, 'd: Number of rejected null frames'
  PRINT, 'e: Preserved number of frames'  
  PRINT, ' '
  PRINT, ' A ', '  ', '    B   ','  ',' C ','  ','     D     ','  ','   E  ','  ','    F  ','  ','  G ','  ','  H ','  ','  I  ','  ','   J ','  ','  K  ','  ','      L     ','  ','   M   ','  ','     N     ','  ','    O     ','  ','  P  ','  ','  Q  ','  ','    R  ',' ','   S  ','  ','    T ','  ','   U  ',' ','   V  ','  ','   W  ','   ','   X  ',' ','  Y  ','  ','   Z  ','  ','   a   ','  ','  b  ','  ','  c  ','  ','  d  ','  ','  e  ','  '
  PRINT, '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
  FOR i = 0, n_data-1 DO PRINT, STRING(l1data[i].ob_id, FORMAT='(I03)'), '  ', STRING(objname[i], FORMAT='(A8)'), '  ', STRING(flag[i], FORMAT='(A3)'), '  ', STRING(utc[i], FORMAT='(A11)'), '  ', STRING(alt[i], FORMAT='(F6.2)'), '  ', STRING(para[i], FORMAT='(F7.2)'), '  ', STRING(seeing[i] ,FORMAT='(F4.2)'), '  ', STRING(smttau[i], FORMAT='(F4.2)'), $
                         '  ', STRING(1D+3*dit[i], FORMAT='(F5.1)'), '  ', STRING(dith_per[i], FORMAT='(I04)'), '  ', STRING(1D+6*wav[i], FORMAT='(F5.1)'), '  ', STRING(1D+6*bdw[i], FORMAT='(F4.1)'), ' | ', STRING(l1data[i].photdx_snr, FORMAT='(F5.1)'), ' ', STRING(l1data[i].photsx_snr, FORMAT='(F5.1)'),$
                         '  ', STRING(int_err[i], FORMAT='(F7.4)'), '  ', STRING(l1data[i].TTDX_RMS, FORMAT='(F5.1)'), ' ', STRING(l1data[i].TTSX_RMS, FORMAT='(F5.1)'), '  ', STRING(100.*l1data[i].photdx_rms/l1data[i].phot_avg, FORMAT='(F4.2)'), '  ', STRING(100.*l1data[i].photsx_rms/l1data[i].phot_avg, FORMAT='(F4.2)'), $
                         '  ', STRING(100.*l1data[i].bckg_cnod_rms, FORMAT='(F5.2)'), '  ', STRING(100.*l1data[i].bckg_bias, FORMAT='(F5.2)'), ' || ', STRING(100.*l1data[i].null_meas, FORMAT='(F5.2)'), '  ', STRING(100.*l1data[i].null_offset, FORMAT='(F5.2)'), ' || ', STRING(100.*err_phot[i], FORMAT='(F5.2)'), '  ', STRING(100.*err_bckg[i], FORMAT='(F5.2)'), '  ', STRING(100.*err_nsc[i], FORMAT='(F5.2)'), $
                         ' | ', STRING(100.*l1data[i].null_meas_err, FORMAT='(F5.2)'), ' || ',  STRING(100.*null_err_phot[i], FORMAT='(F5.2)'), ' ',  STRING(100.*null_err_phase[i], FORMAT='(F5.2)'), ' ||', STRING(1D6*l1data[i].nsc_phavg*wav[i]/(2*!dpi), FORMAT='(F6.2)'),'  ', STRING(1D6*l1data[i].nsc_phrms*wav[i]/(2*!dpi), FORMAT='(F6.2)'), '  ', STRING(l1data[i].nsc_chi2, FORMAT='(F6.2)'), $
                         '   ', STRING(l1data[i].nfr_rej, FORMAT='(I04)'), '  ', STRING(l1data[i].nfr_ob, FORMAT='(I04)')                               
  PRINT, ''
  IF lun GT 0 THEN BEGIN
    PRINTF, lun, '==== REDUCTION INFO ===='
    PRINTF, lun, ''
    PRINTF, lun, 'File version number     : ' + STRING(file_ver, FORMAT='(I0)')
    PRINTF, lun, 'Input aperture radius   : ' + STRING(aper_rad, FORMAT='(I0)') + '  (0 for EEID + FWHM)'
    ;PRINTF, lun, 'Background inner radius : ' + STRING(bck_irad, FORMAT='(I0)')
    ;PRINTF, lun, 'Background outer radius : ' + STRING(bck_orad, FORMAT='(I0)')
    PRINTF, lun, 'Image combination mode  : ' + STRING(img_mode, FORMAT='(I0)') + '  (0 for automatic)'
    PRINTF, lun, 'Background nod mode     : ' + STRING(bck_mode, FORMAT='(I0)') + '  (0 none, 1 nod pairs, 2 closest n frames, 3 dedicated, 3 chopping)'
    PRINTF, lun, 'Background floor mode   : ' + STRING(bfl_mode, FORMAT='(I0)') + '  (0 5-sigma clipped, 1 median)'
    PRINTF, lun, 'Frame mode for flux     : ' + STRING(fra_mode, FORMAT='(I0)') + '  (0 background subtracted, 1: processed raw)'
    PRINTF, lun, 'Flux mode               : ' + STRING(flx_mode, FORMAT='(I0)') + '  (0 aper. phot., 1: weigh. aper phot., 2: PSF-fitting)'
    PRINTF, lun, 'Centroid mode           : ' + STRING(fit_mode, FORMAT='(I0)') + '  (0 none, 1: CNTRD., 2 GCNTRD, 3: Gaussian, 4: Lorentzan, 5: Moffat)'  
    PRINTF, lun, 'Null estimator          : ' + STRING(drs.null_est, FORMAT='(I0)')
    PRINTF, lun, 'Null error mode         : ' + STRING(drs.err_mode, FORMAT='(I0)')
    PRINTF, lun, 'Scatter-rejection sigma : ' + STRING(drs.sig_sca, FORMAT='(I0)')
    PRINTF, lun, 'Null correction mode    : ' + STRING(null_cor, FORMAT='(I0)') + "  (1: Denis', 2: Bertrand's method to subtract HF noise)"
    PRINTF, lun, 'NSC cube size           : ' + nsc_cube + '  (0 for automatic)'
    PRINTF, lun, 'NSC mode                : ' + STRING(nsc_mode, FORMAT='(I0)') + '  (0: mode, 1: %best, 2: NSC)' 
    PRINTF, lun, 'Acceptable null range   : ' + null_lim
    PRINTF, lun, 'Number of bin factor    : ' + STRING(nbin_fac , FORMAT='(I0)') + '  (Bin size for NSC reduction (0: constant, 1: variable)'
    PRINTF, lun, 'Number of bootstrap     : ' + STRING(n_btstrp, FORMAT='(I0)')
    PRINTF, lun, 'Calibration method      : ' + STRING(drs.cal_method, FORMAT='(I0)') + '  (0: interpolation nearest neighbors, 2: interpolation on all OBs)'
    PRINTF, lun, 'Degree of polynomial    : ' + STRING(drs.polydeg, FORMAT='(I0)') 
    PRINTF, lun, 'Calibration mode        : ' + STRING(drs.cal_mode, FORMAT='(I0)') + '  (0: per pointing, 1: per OB)'
    PRINTF, lun, 'Reduced chi2 limit      : ' + STRING(drs.chi2_lim, FORMAT='(i0)')
    PRINTF, lun, 'TF error mode           : ' + STRING(drs.err_mode, FORMAT='(I0)') + '  ( 0: unweighted null dispersion, 1: weighted null dispersion, 2: excess variance)'
    PRINTF, lun, 'Nod correlation         : ' + STRING(drs.nod_cor, FORMAT='(I0)') + '  ( 1 if nulls appear correlated per nod within a given pointing)'
    PRINTF, lun, 'Split TF                : ' + STRING(drs.split_tf, FORMAT='(I0)') + '  (1 to split the TF when there is a dead time longer than split_time between two null measurements)'
    PRINTF, lun, 'Split time              : ' + STRING(drs.split_time, FORMAT='(F3.1)') 
    PRINTF, lun, 'Split hour              : ' + STRING(drs.split_hour, FORMAT='(I0)') 
    PRINTF, lun, 'Split nod               : ' + STRING(drs.split_nod, FORMAT='(I0)') 
    PRINTF, lun, ''
    PRINTF, lun, 'Number of OBs           : ' + STRING(n_data, FORMAT='(I03)')
    PRINTF, lun, 'Number of pointings     : ' + STRING(n_pt, FORMAT='(I03)')
    ;PRINTF, lun, 'PLC wavelength          : ' + STRING(plc_wav, FORMAT='(F5.2)')
    PRINTF, lun, ''
    PRINTF, lun,''
    PRINTF, lun,'==== OB-BASED INFORMATION ===='
    ; Column signification
    PRINTF, lun, 'Column signification'
    PRINTF, lun, 'A: OB identification number'
    PRINTF, lun, 'B: Object name'
    PRINTF, lun, 'C: Flag (SCI/CAL)'
    PRINTF, lun, 'D: UT time (from LBT)'
    PRINTF, lun, 'E: Telescope elevation in degrees (from LBT)'
    PRINTF, lun, 'F: Parallactic angle in degrees (from LBT)'
    PRINTF, lun, 'G: DIMM seeing [arcsec]'
    PRINTF, lun, 'H: PWV [mm]'
    PRINTF, lun, 'I: DIT [ms]'
    PRINTF, lun, 'J: Dithering period [number of frames]'
    PRINTF, lun, 'K: Wavelength [microns]'
    PRINTF, lun, 'L: Bandwidths [microns]'
    PRINTF, lun, 'M: SNR on photometry (DX and SX)' 
    PRINTF, lun, 'N: Intensity mismatch (in %)'  
    PRINTF, lun, 'O: Tip/tilt error (DX and SX, in mas)'
    PRINTF, lun, 'P: RMS of photometry relative to constructive peak (DX and SX, in %)'  
    PRINTF, lun, 'Q: RMS of background relative to constructive peak (in %)'
    PRINTF, lun, 'R: Background bias relative to constructive peak (in %)'  
    PRINTF, lun, 'S: Null (in %)'
    PRINTF, lun, 'T: Null offset between optimum NSC value and bootstrapped or bayesian value (in %)'  
    PRINTF, lun, 'U: Null uncertainty due to error on constructive peak (in %)'
    PRINTF, lun, 'V: Null uncertainty due to external error on background floor (in %)'
    PRINTF, lun, 'W: Null uncertainty from NSC fit (in %)'  
    PRINTF, lun, 'X: Total null uncertainty (in %)'  
    PRINTF, lun, 'Y: Expected null uncertainty due to photometric errors (in %)'
    PRINTF, lun, 'Z: Expected null uncertainty due to phase variation (i.e., sqrt(4mu^2.sig^2 + 2*sig^4)/4, in %)'  
    PRINTF, lun, 'a: Best-fit mean PHASE from NSC (in microns)'
    PRINTF, lun, 'b: Best-fit PHASE jitter from NSC (in microns)'
    PRINTF, lun, 'c: Reduced chi2 from NSC'
    PRINTF, lun, 'd: Number of rejected null frames'
    PRINTF, lun, 'e: Preserved number of frames'  
    PRINTF, lun, ' '
    PRINTF, lun, ' A ', '  ', '    B   ','  ',' C ','  ','     D     ','  ','   E  ','  ','    F  ','  ','  G ','  ','  H ','  ','  I  ','  ','   J ','  ','  K  ','  ','      L     ','  ','   M   ','  ','     N     ','  ','    O     ','  ','  P  ','  ','  Q  ','  ','    R  ',' ','   S  ','  ','    T ','  ','   U  ',' ','   V  ','  ','   W  ','   ','   X  ',' ','  Y  ','  ','   Z  ','  ','   a   ','  ','  b  ','  ','  c  ','  ','  d  ','  ','  e  ','  '
    PRINTF, lun, '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    FOR i = 0, n_data-1 DO PRINTF, lun, STRING(l1data[i].ob_id, FORMAT='(I03)'), '  ', STRING(objname[i], FORMAT='(A8)'), '  ', STRING(flag[i], FORMAT='(A3)'), '  ', STRING(utc[i], FORMAT='(A11)'), '  ', STRING(alt[i], FORMAT='(F6.2)'), '  ', STRING(para[i], FORMAT='(F7.2)'), '  ', STRING(seeing[i] ,FORMAT='(F4.2)'), '  ', STRING(smttau[i], FORMAT='(F4.2)'), $
                           '  ', STRING(1D+3*dit[i], FORMAT='(F5.1)'), '  ', STRING(dith_per[i], FORMAT='(I04)'), '  ', STRING(1D+6*wav[i], FORMAT='(F5.1)'), '  ', STRING(1D+6*bdw[i], FORMAT='(F4.1)'), ' | ', STRING(l1data[i].photdx_snr, FORMAT='(F5.1)'), ' ', STRING(l1data[i].photsx_snr, FORMAT='(F5.1)'),$
                           '  ', STRING(int_err[i], FORMAT='(F7.4)'), '  ', STRING(l1data[i].TTDX_RMS, FORMAT='(F5.1)'), ' ', STRING(l1data[i].TTSX_RMS, FORMAT='(F5.1)'), '  ', STRING(100.*l1data[i].photdx_rms/l1data[i].phot_avg, FORMAT='(F4.2)'), '  ', STRING(100.*l1data[i].photsx_rms/l1data[i].phot_avg, FORMAT='(F4.2)'), $
                           '  ', STRING(100.*l1data[i].bckg_cnod_rms, FORMAT='(F5.2)'), '  ', STRING(100.*l1data[i].bckg_bias, FORMAT='(F5.2)'), ' || ', STRING(100.*l1data[i].null_meas, FORMAT='(F5.2)'), '  ', STRING(100.*l1data[i].null_offset, FORMAT='(F5.2)'), ' || ', STRING(100.*err_phot[i], FORMAT='(F5.2)'), '  ', STRING(100.*err_bckg[i], FORMAT='(F5.2)'), '  ', STRING(100.*err_nsc[i], FORMAT='(F5.2)'), $
                           ' | ', STRING(100.*l1data[i].null_meas_err, FORMAT='(F5.2)'), ' || ',  STRING(100.*null_err_phot[i], FORMAT='(F5.2)'), ' ',  STRING(100.*null_err_phase[i], FORMAT='(F5.2)'), ' ||', STRING(1D6*l1data[i].nsc_phavg*wav[i]/(2*!dpi), FORMAT='(F6.2)'),'  ', STRING(1D6*l1data[i].nsc_phrms*wav[i]/(2*!dpi), FORMAT='(F6.2)'), '  ', STRING(l1data[i].nsc_chi2, FORMAT='(F6.2)'), $
                           '   ', STRING(l1data[i].nfr_rej, FORMAT='(I04)'), '  ', STRING(l1data[i].nfr_ob, FORMAT='(I04)')                               
    PRINTF, lun, ''
  ENDIF
ENDIF

; SUBTRACT GEOMETRIC NULL FROM CALIBRATOR NULL
; ********************************************

; Find index of calibrator data
idx_cal = WHERE(STRMATCH(flag, 'CAL') EQ 1, n_caldat)
IF n_caldat LE 0 THEN BEGIN
  PRINT, 'No calibrator data => no calibration'
  IF lun GT 0 THEN PRINTF, lun, 'No calibrator data => no calibration'
  GOTO, SKIP_CAL
ENDIF

; Derive number of unique wavelength
cal_wav  = wav[idx_cal]
wav_uniq = cal_wav[UNIQ(cal_wav, SORT(cal_wav))]
n_wav    = N_ELEMENTS(wav_uniq)

; Loop over the wavelengths
FOR i_wav = 0, n_wav-1 DO BEGIN  
  ; Extract data of this wavelength 
  idx_wav  = idx_cal[WHERE(cal_wav EQ wav_uniq[i_wav])]
  calname  = objname[idx_wav]
  calbdw   = bdw[idx_wav]
  
  ; Check whether there is only one bandwidth for this wavelength
  bdw_uniq = calbdw[UNIQ(calbdw, SORT(calbdw))]
  IF N_ELEMENTS(bdw_uniq) GT 1 THEN MESSAGE, 'Found two different bandwidths for the same wavelength'
         
  ; Compute wavelength range
  lambda_chnl = wav_uniq[i_wav] - 0.5*bdw_uniq[0] + bdw_uniq[0] * (1.+1./n_lam) * DINDGEN(n_lam)/(n_lam-1)
  
  ; Extract LBTI transmission profile for this wavalength (N' and 8.7um only now)
  IF wav_uniq[i_wav] GE 1.1D-5 AND wav_uniq[i_wav] LE 1.2D-5 THEN BEGIN
    READ_TABLE, 'nodrs/input/n-band_thruput.txt', lam_tmp, thruput, FIRST=4, SEPARATOR='TAB'
    trans = INTERPOL(thruput, lam_tmp, lambda_chnl*1D6)
  ENDIF ELSE BEGIN
    IF wav_uniq[i_wav] GE 8.6D-6 AND wav_uniq[i_wav] LT 8.8D-6 THEN BEGIN
      READ_TABLE, 'nodrs/input/f87_thruput.txt', lam_tmp, thruput, FIRST=4, SEPARATOR='TAB'
      trans = INTERPOL(thruput, lam_tmp, lambda_chnl*1D6)
    ENDIF ELSE BEGIN
      IF wav_uniq[i_wav] GE 8.8D-6 AND wav_uniq[i_wav] LE 9.0D-6 THEN BEGIN
        READ_TABLE, 'nodrs/input/n08909_thruput.txt', lam_tmp, thruput, FIRST=4, SEPARATOR='TAB'
        trans = INTERPOL(thruput, lam_tmp, lambda_chnl*1D6) 
      ENDIF ELSE BEGIN 
          IF wav_uniq[i_wav] GE 12D-6 AND wav_uniq[i_wav] LE 12.5D-6 THEN BEGIN
            READ_TABLE, 'nodrs/input/n12520_thruput.txt', lam_tmp, thruput, FIRST=4, SEPARATOR='TAB'
            trans = INTERPOL(thruput, lam_tmp, lambda_chnl*1D6) 
          ENDIF ELSE trans = 1
      ENDELSE
    ENDELSE 
  ENDELSE
  
  ; Compute effective wavelength
  lam_eff  = TOTAL(trans*lambda_chnl)/TOTAL(trans)

  ; Compute unique calibrators for this wavelength
  nam_uniq = calname[UNIQ(calname, SORT(calname))]
  n_cal    = N_ELEMENTS(nam_uniq)
  
  ; Loop over calibrator data
  FOR i_cal = 0, n_cal - 1 DO BEGIN
    ; Extract data for this calibrator
    idx_cur = idx_wav[WHERE(calname EQ nam_uniq[i_cal])]
    
    ; Retrieve target information
    GET_TGT, nam_uniq[i_cal], star, DATABASE= pth.input_path + drs.database
    th_mean = star.ldm
    th_min  = star.ldm - star.lde
    th_max  = star.ldm + star.lde
    uld     = star.uld
    
    ; Generate stellar spectrum for broadband null computation
    star_flux        = BLACKBODY(star.temp, lambda_chnl, STANDARD=3)                       ; flux in photons/s/m²/Hz/sr
    flux_int         = TOTAL(trans*star_flux*!Dpi*(0.5*th_mean*prm.m2r)^2)/TOTAL(trans)    ; Integrated flux in ph/s/m2/Hz
    tgt_flx[idx_cur] = flux_int*prm.h*prm.c*prm.w2jy/lam_eff                               ; Integrated flux in Jansky
    
    ; Compute the expected null for this star
    norm         = TOTAL((trans*star_flux)^2)
    null_lam     = GEOMETRIC_NULL(th_mean, cnf.base, lambda_chnl, ULD=uld)
    null_wb      = TOTAL((trans*star_flux)^2 * null_lam) / norm
    ; Compute minimum null (with minimum diameter)
    null_min     = GEOMETRIC_NULL(th_min, cnf.base, lambda_chnl, ULD=uld)
    null_wb_min  = TOTAL((trans*star_flux)^2 * null_min) / norm
    ; Compute maximum null (with maximum diameter)
    null_max     = GEOMETRIC_NULL(th_max, cnf.base, lambda_chnl, ULD=uld)
    null_wb_max  = TOTAL((trans*star_flux)^2 * null_max) / norm
    ; Derive the error bar on the estimated null
    null_err_sys = (null_wb_max - null_wb_min) / 2D0
    
    ; Derive the wide-band instrumental transfer function and its error bar
    ; Statistical error bar remains unchanged
    tf_wb[idx_cur]      -= null_wb                                             ; correct for geometric leakage 
    tf_wb_esyst[idx_cur] = null_err_sys                                        ; systematic error bar
    tf_wb_etot[idx_cur]  = SQRT(tf_wb_estat[idx_cur]^2+tf_wb_esyst[idx_cur]^2) ; total error bar
  ENDFOR
ENDFOR

; COMPUTE MEAN NULL AND ERROR PER POINTING/NOD
; ********************************************

; Number of unique pointings
IF info GT 0 THEN BEGIN
  PRINT, '==== POINTING-BASED NULLS (NO TIME GATING!) ===='
  PRINT, '   '
  PRINT, '  1. Instrumental null floor per pointing (pointing, name, flag, null[%], w_null[%], stat error[%], disp error[%], w_disp error[%], excess error [%])  '
  PRINT, '   '  
  IF lun GT 0 THEN BEGIN
    PRINTF, lun, '==== POINTING-BASED NULLS (NO TIME GATING!) ===='
    PRINTF, lun, '   '
    PRINTF, lun, '  1. Instrumental null floor per pointing (pointing, name, flag, null[%], w_null[%], stat error[%], disp error[%], w_disp error[%], excess error [%])'
    PRINTF, lun, '   '
  ENDIF
ENDIF

; Compute weigthed instrumental null per pointing and per nod
flag_pt          = STRARR(n_pt)       & rms_exc          = FLTARR(n_pt)       
time_avg_pt      = FLTARR(n_pt,n_nod) & null_avg_pt      = FLTARR(n_pt,n_nod) & null_err_pt       = FLTARR(n_pt,n_nod) 
null_err_stat_pt = FLTARR(n_pt,n_nod) & null_err_disp_pt = FLTARR(n_pt,n_nod) & null_err_wdisp_pt = FLTARR(n_pt,n_nod)
FOR i_pt = 0, n_pt-1 DO BEGIN
  ; Derive OB of this pointing and parse stats results
  idx_pt = WHERE(pt_id EQ pt_uniq[i_pt], n_ob)
  AVGERR, tf_wb[idx_pt], tf_wb_estat[idx_pt], avg_unw, avg_wei, err_stat, err_sca_unw, err_sca_wei, exc_rms, KAPPA=drs.sig_sca
  flag_pt[i_pt] = flag[idx_pt[0]]
  rms_exc[i_pt] = exc_rms
  ; Now loop over the nods
  FOR i_nod = 0, n_nod-1 DO BEGIN
    ; Derive nod of this pointing
    idx_npt = idx_pt[WHERE(nod_pos[idx_pt] EQ nod_uniq[i_nod], n_tmp)]
    ; Compute mean time
    time_avg_pt[i_pt,i_nod] = MEAN(time[idx_npt])
    ; Compute average value and error bars
    AVGERR, tf_wb[idx_npt], tf_wb_estat[idx_npt], avg_unw, avg_wei, err_stat, err_sca_unw, err_sca_wei, sys_ecc, KAPPA=drs.sig_sca
    IF drs.null_est EQ 0 THEN null_est = avg_unw ELSE null_est = avg_wei
    null_avg_pt[i_pt,i_nod]       = null_est
    null_err_stat_pt[i_pt,i_nod]  = err_stat
    null_err_disp_pt[i_pt,i_nod]  = err_sca_unw
    null_err_wdisp_pt[i_pt,i_nod] = err_sca_wei
    ; Systematic error
    CASE (ABS(drs.err_mode) MOD 3) OF 
      0: rms_sys = null_err_disp_pt[i_pt,i_nod]
      1: rms_sys = null_err_wdisp_pt[i_pt,i_nod]
      2: rms_sys = sys_ecc
      ELSE: MESSAGE, 'Undefined error mode'
    ENDCASE
    ; Compute total error (use maximum of the two if greater than 2)
    IF drs.err_mode GE 3 THEN err_tot = null_err_stat_pt[i_pt,i_nod] > rms_sys ELSE err_tot = SQRT(null_err_stat_pt[i_pt,i_nod]^2+rms_sys^2) 
    null_err_pt[i_pt,i_nod] = SQRT(err_tot^2+MEAN(tf_wb_esyst[idx_npt])^2)  ; in most cases, the systematic error due to the uncertainty of the stellar diameter is the same for all points
    ; Print info to screen
    IF info GT 0 THEN BEGIN
      PRINT, STRING(i_pt, FORMAT='(I02)'), '  ', STRING(objname[idx_npt[0]], FORMAT='(A8)'), '  ', flag[idx_npt[0]], '  ', STRING(1D2*avg_unw, FORMAT='(F7.4)'), '  ', STRING(1D2*avg_wei, FORMAT='(F7.4)'), '  ', STRING(1D2*null_err_stat_pt[i_pt,i_nod], FORMAT='(F6.4)'), '  ', STRING(1D2*err_sca_unw, FORMAT='(F6.4)'), '  ', STRING(1D2*err_sca_wei, FORMAT='(F6.4)'), '  ', STRING(1D2*sys_ecc, FORMAT='(F6.4)')
      IF lun GT 0 THEN PRINTF, lun, STRING(i_pt, FORMAT='(I02)'), '  ', STRING(objname[idx_npt[0]], FORMAT='(A8)'), '  ', flag[idx_npt[0]], '  ', STRING(1D2*avg_unw, FORMAT='(F7.4)'), '  ', STRING(1D2*avg_wei, FORMAT='(F7.4)'), '  ', STRING(1D2*null_err_stat_pt[i_pt,i_nod], FORMAT='(F6.4)'), '  ', STRING(1D2*err_sca_unw, FORMAT='(F6.4)'), '  ', STRING(1D2*err_sca_wei, FORMAT='(F6.4)'), '  ', STRING(1D2*sys_ecc, FORMAT='(F6.4)')
    ENDIF
  ENDFOR
ENDFOR

; Merge nods if requested
IF drs.split_nod EQ 0 AND n_nod GT 1 THEN BEGIN
  IF info GT 0 THEN BEGIN
    PRINT, '   '
    PRINT, '  2. Instrumental null floor per pointing (pointing, null[%], stat error[%], disp error[%], w_disp error[%], excess error [%], tot error [%])'
    IF lun GT 0 THEN PRINTF, lun,  '   '
    IF lun GT 0 THEN PRINTF, lun, '  2. Instrumental null floor per pointing (pointing, null[%], stat error[%], disp error[%], w_disp error[%], excess error [%], tot error[%])  '
  ENDIF
  time_avg_pt      = MEAN(time_avg_pt, DIMENSION=2)
  null_err_stat_pt = MEAN(null_err_stat_pt, DIMENSION=2)/SQRT(n_nod)
  null_avg_pt0 = null_avg_pt  & null_err_pt0 = null_err_pt  & null_err_disp_pt0 = null_err_disp_pt & null_err_wdisp_pt0 = null_err_wdisp_pt
  null_avg_pt  = FLTARR(n_pt) & null_err_pt  = FLTARR(n_pt) & null_err_disp_pt  = FLTARR(n_pt)     & null_err_wdisp_pt  = FLTARR(n_pt) 
  FOR i_pt = 0, n_pt-1 DO BEGIN  
    ; Compute dispersion errors 
    AVGSDV, null_avg_pt0[i_pt,*], avg_tmp1, rms_disp, rmsm_disp, KAPPA=drs.sig_sca
    AVGSDV, null_avg_pt0[i_pt,*], avg_tmp2, rms_wdisp, rmsm_wdisp, KAPPA=drs.sig_sca, WEIGHT=1./null_err_pt0[i_pt,*]^2
    IF drs.null_est EQ 0 THEN null_est = avg_tmp1 ELSE null_est = avg_tmp2
    null_avg_pt[i_pt]       = null_est
    null_err_disp_pt[i_pt]  = rmsm_disp
    null_err_wdisp_pt[i_pt] = rmsm_wdisp
    ; Systematic error. RMS_SYS is added quadratically to the statiscal error, unless err_mode is greater than 2 (in that case, the maximum of the 2 is used)
    CASE (ABS(drs.err_mode) MOD 3) OF
      0: rms_sys = null_err_disp_pt[i_pt]
      1: rms_sys = null_err_wdisp_pt[i_pt]
      2: rms_sys = rms_exc[i_pt]
      ELSE: MESSAGE, 'Undefined error mode'
    ENDCASE
    ; Compute total error (use maximum of the two if greater than 2)
    IF drs.err_mode GE 3 THEN null_err_pt[i_pt] = null_err_stat_pt[i_pt] > rms_sys ELSE null_err_pt[i_pt] = SQRT(null_err_stat_pt[i_pt]^2+rms_sys^2)
    IF info GT 0 THEN BEGIN
      PRINT, STRING(i_pt, FORMAT='(I02)'), '  ', STRING(1D2*null_avg_pt[i_pt], FORMAT='(F7.4)'), '  ', STRING(1D2*null_err_stat_pt[i_pt], FORMAT='(F6.4)'), '  ', STRING(1D2*null_err_disp_pt[i_pt], FORMAT='(F6.4)'), '  ', STRING(1D2*null_err_wdisp_pt[i_pt], FORMAT='(F6.4)'), '  ', STRING(1D2*rms_exc[i_pt], FORMAT='(F6.4)'), '  ', STRING(1D2*null_err_pt[i_pt], FORMAT='(F6.4)')
      IF lun GT 0 THEN PRINTF, lun, STRING(i_pt, FORMAT='(I02)'), '  ', STRING(1D2*null_avg_pt[i_pt], FORMAT='(F7.4)'), '  ', STRING(1D2*null_err_stat_pt[i_pt], FORMAT='(F6.4)'), '  ', STRING(1D2*null_err_disp_pt[i_pt], FORMAT='(F6.4)'), '  ', STRING(1D2*null_err_wdisp_pt[i_pt], FORMAT='(F6.4)'), '  ', STRING(1D2*rms_exc[i_pt], FORMAT='(F6.4)'), '  ', STRING(1D2*null_err_pt[i_pt], FORMAT='(F6.4)')
    ENDIF
  ENDFOR
ENDIF

; Extract CAL pointings (for later)
idx_cal_pt = WHERE(STRMATCH(flag_pt, 'CAL') EQ 1, n_caldat_pt)
caltime    = time[idx_cal]
caltime_pt = time_avg_pt[idx_cal_pt,*]

; READ SCIENCE DATA
; *****************

; Jumping point if no calibrator data
SKIP_CAL:

; Search for science target fits files in the input folder
idx_sci = WHERE(STRMATCH(flag, 'SCI') EQ 1, n_scidat)
IF n_scidat GT 0 THEN scidata = l1data[idx_sci] ELSE  GOTO, SKIP_INTERP 
idx_sci_pt = WHERE(STRMATCH(flag_pt, 'SCI') EQ 1, n_scidat_pt)

; Extract SCI times
scitime    = time[idx_sci] 
scitime_pt = time_avg_pt[idx_sci_pt,*]
sci_pt     = pt_uniq[idx_sci_pt]

; Compute number of nod positions to consider (otherwise, it uses nod_pos as defined above)
IF drs.split_nod EQ 0 THEN nod_pos = INTARR(n_data)
nod_uniq = nod_pos[UNIQ(nod_pos, SORT(nod_pos))]
n_nod    = N_ELEMENTS(nod_uniq)

; Define output arrays
IF drs.cal_mode EQ 1 THEN n_out = n_scidat ELSE n_out = n_scidat_pt
IF drs.cal_mode EQ 1 THEN n_bis = 1        ELSE n_bis = n_nod
tf_sci   = DBLARR(n_out,n_bis) & tf_estat       = tf_sci    & tf_esyst       = tf_sci    & tf_etot       = tf_sci    ; will contain the estimate of the TF at the science times (and error bars)
null_sci = DBLARR(n_out,n_bis) & null_sci_estat = null_sci  & null_sci_esyst = null_sci  & null_sci_etot = null_sci  ; will contain the final calibrated null and corresponding error bars


; INTERPOLATE THE TF OVER THE NIGHT
; *********************************

; Define the number of sub-TFs based on maximum dead time allowed
; This is badly coded...each element of idx_split gives the beginning of the new sub TF
IF drs.split_tf NE 0 AND drs.split_time NE 0 THEN BEGIN
  idx_split = WHERE(ABS(l1data[1:*].mjd_obs-l1data[0:(n_data-1)].mjd_obs)*24. GT drs.split_time, n_tf)
  IF n_tf GT 0 THEN idx_split = [0, idx_split+1, n_data] ELSE idx_split = [0, n_data]
ENDIF ELSE idx_split = [0, n_data]

; Add sub-TFs based on user-defined hours at which the TF has to be split
IF drs.split_hour[0] NE 0 THEN BEGIN
  idx_tmp = WHERE((l1data.mjd_obs- t0day)*24-drs.split_hour GT 0, n_tf)
  idx_tmp = idx_tmp[0]
  IF n_tf GT 0 THEN BEGIN
    idx_split = [idx_split, idx_tmp]
    idx_split = idx_split[SORT(idx_split)]
  ENDIF
ENDIF

; Now split if aper radius changes (when it's EEID dependend)
FOR i = 1, N_ELEMENTS(aper)-1 DO IF aper[i] NE aper[i-1] THEN idx_split = [idx_split, i]
idx_split = idx_split[UNIQ(idx_split,SORT(idx_split))]

; Also split if aper radius changes
; Inner radius
FOR i = 1, N_ELEMENTS(birad)-1 DO IF birad[i] NE birad[i-1] THEN idx_split = [idx_split, i]
idx_split = idx_split[UNIQ(idx_split,SORT(idx_split))]
; Outer radius (no, don't use it anymore because stars are likely not on the same spot of the detector anyway if 'borad' is different)
;FOR i = 1, N_ELEMENTS(borad)-1 DO IF borad[i] NE borad[i-1] THEN idx_split = [idx_split, i]
;idx_split = idx_split[UNIQ(idx_split,SORT(idx_split))]

; Derive number of sub-TFs
n_tf = N_ELEMENTS(idx_split)-1

; Prepare array for interpolated TF
tf_time    = FLTARR(n_tf, n_time)
tf_plot    = tf_time & tf_plot_err = FLTARR(n_tf)
tf_plot_pt = FLTARR(n_tf, n_time, n_nod) & tf_plot_pt_err = FLTARR(n_tf, n_nod)

; Loop over the TFs
FOR i_tf = 0, n_tf-1 DO BEGIN
  ; Extract data of this sub-TF
  IF i_tf NE n_tf-1 THEN idx_tf = WHERE(time GE time[idx_split[i_tf]] AND time LT time[idx_split[i_tf+1]], n_intf) $
                    ELSE idx_tf = WHERE(time GE time[idx_split[i_tf]], n_intf)
  flag_tf    = flag[idx_tf]
  idx_cal_tf = idx_tf[WHERE(STRMATCH(flag_tf, 'CAL') EQ 1, n_cal_tf)]
  idx_sci_tf = idx_tf[WHERE(STRMATCH(flag_tf, 'SCI') EQ 1, n_sci_tf)] 
  
  ; Skip this TF if no science and no cal frames
  IF n_sci_tf GT 0 AND n_cal_tf GT 0 THEN BEGIN
  
    ; Number of unique calibrators
    cal_name = objname[idx_cal_tf]
    cal_uniq = cal_name(UNIQ(cal_name,SORT(cal_name)))
    n_cal    = N_ELEMENTS(cal_uniq)
     
    ; Define plot range
    min_time        = 0.98*MIN(time[idx_tf])
    max_time        = 1.02*MAX(time[idx_tf])
    plottime        = DINDGEN(n_time)/(DOUBLE(n_time)-1D0) * (max_time-min_time) + min_time
    tf_time[i_tf,*] = plottime 
    
    ; Remove bad points (see rules within AVGERR routine)
    AVGERR, tf_wb[idx_cal_tf], tf_wb_estat[idx_cal_tf], avg_null, avg_null_we, rms_stat_null, disp_m_null, disp_m_null_we, rms_exc, KAPPA=drs.sig_sca, IDX_OUT=idx_out, N_OUT=n_out
    IF n_out GT 0 THEN ob_out[idx_cal_tf[idx_out]] = 1
    idx_in     = idx_cal_tf[WHERE(ob_out[idx_cal_tf] NE 1)]
    n_cal_tf_f = N_ELEMENTS(null_cal)
    
    ; Calibation by linear interpolation between neighbours (REMOVING bad points!)
    IF drs.cal_method EQ 0 THEN BEGIN
      ; To be implemented
      MESSAGE, 'Not yet implemented'
    ENDIF ELSE BEGIN  ; Calibation by polynomial interpolation on all calibrator data
      ; Interpolate the data at data and plot times
      polydeg = drs.polydeg
      IF polydeg GT n_cal THEN MESSAGE, 'Degree of polynomial fit should not exceed the number of calibrator data points.'
      polyfit = POLY_FIT(time[idx_in], tf_wb[idx_in], polydeg, /DOUBLE, MEASURE_ERRORS=tf_wb_estat[idx_in]^drs.null_est, CHISQ=chi2, COVAR=covar, SIGMA=polysig, YFIT=yfit, YBAND=yband, YERROR=yerror) ; include only statistical error bars 
      tf_data = 0D & tf_tmp = 0D
      FOR ipol=0, polydeg DO BEGIN
        tf_data += polyfit[ipol]*time[idx_tf]^ipol
        tf_tmp  += polyfit[ipol]*plottime^ipol
      ENDFOR
      tf_plot[i_tf,*] = tf_tmp
    ENDELSE  
    
    ; Calibrate null
    null_c = tf_wb[idx_tf] - tf_data
    
    ; Extract calibrator data
    null_cal     = null_c[WHERE(STRMATCH(flag_tf, 'CAL') EQ 1, ntmp)]
    null_cal_err = tf_wb_estat[idx_cal_tf]     
           
    ; Compute mean and errors on TF: OB BASED and also flag points too far from TF
    AVGERR, null_cal, null_cal_err, avg_null, avg_null_we, rms_stat_null, disp_m_null, disp_m_null_we, rms_exc, KAPPA=drs.sig_sca, IDX_OUT=idx_out, N_OUT=n_out
    IF n_out GT 0 THEN ob_out[idx_cal_tf[idx_out]] = 1
    n_cal_tf_f = N_ELEMENTS(null_cal)
                 
    ; Systematic error defined by the err_mode input parameter
    ; This term will be added to the statistcal error
    CASE (ABS(drs.err_mode) MOD 3) OF
      0: rms_sys = disp_m_null
      1: rms_sys = disp_m_null_we
      2: rms_sys = rms_exc
      ELSE: MESSAGE, 'Undefined error mode'
    ENDCASE
    
    ; Compute the systematic error on TF estimate (due to uncertainty on stellar diameter)
    ; Loop over the diffferent calibrators (systematic error is fully correlated between different OBs on the same calibrator)
    null_err_diam = 0D
    FOR ical = 0, n_cal-1 DO null_err_diam += (TOTAL(tf_wb_esyst[idx_cal_tf[WHERE(cal_name EQ cal_uniq[ical], n_c)]])/n_c)^2
    rms_sys_tot = SQRT(rms_sys^2 + null_err_diam^2)
    
    ; Compute total error (use maximum of the two if greater than 2)
    IF drs.err_mode GE 3 THEN null_err_tot = rms_stat_null > rms_sys_tot ELSE null_err_tot = SQRT(rms_stat_null^2+rms_sys_tot^2)
    tf_plot_err[i_tf] = null_err_tot
    
    ; Ok, now save the calibrated null and null floor in output arrays if CAL_MODE is 1
    IF drs.cal_mode EQ 1 THEN BEGIN
      idx_tmp = WHERE(STRMATCH(flag_tf, 'SCI') EQ 1, ntmp)
      FOR i_s = 0, n_sci_tf-1 DO BEGIN
        idx_s                 = WHERE(scitime EQ time[idx_sci_tf[i_s]], n_ok) & IF n_ok LE 0 THEN MESSAGE, 'Problem with science time computation!!'
        null_sci[idx_s]       = null_c[idx_tmp[i_s]]
        tf_sci[idx_s]         = tf_data[i_s]
        tf_estat[idx_s]       = rms_stat_null
        tf_esyst[idx_s]       = rms_sys_tot
        tf_etot[idx_s]        = null_err_tot
        null_sci_estat[idx_s] = SQRT(tf_wb_estat[idx_sci_tf[i_s]]^2+rms_stat_null^2)
        null_sci_esyst[idx_s] = rms_sys_tot
        null_sci_etot[idx_s]  = SQRT(null_sci_estat[idx_s]^2+null_sci_esyst[idx_s]^2)
      ENDFOR
    ENDIF
    
    ; Print diagnostic info to screen and log
    IF info GT 0 THEN BEGIN
      PRINT, ' '
      PRINT, '==== SUB-TF COMPUTATION (OB-based approach) ==== '
      PRINT, '================================================ '
      PRINT, ' '
      PRINT, '==== PART ' + STRING(i_tf+1, FORMAT='(I0)') + ' ==== '
      PRINT, '  Statistical error on TF [%] :'
      PRINT, '    - based on individual error bars    : ', STRING(rms_stat_null*1D+2, FORMAT='(F6.4)')
      PRINT, '    - standard deviation on the mean    : ', STRING(disp_m_null_we*1D+2, FORMAT='(F6.4)') + ' (' + STRING(disp_m_null*1D+2, FORMAT='(F6.4)') + ' if unweighted)'
      ;PRINT, '    - quadratic sum of the 2            : ', STRING(null_err_tot*1D+2, FORMAT='(F6.4)')
      PRINT, '  Systematic error on TF [%]            : '
      PRINT, '    - due to diameter uncertainty       : ', STRING(null_err_diam*1D+2, FORMAT='(F6.4)')
      PRINT, '    - sqrt(excess variance)             : ', STRING(rms_exc*1D+2, FORMAT='(F6.4)')
      PRINT, '  Total error on TF [%]                 : ', STRING(null_err_tot*1D+2, FORMAT='(F6.4)')  
      PRINT, '  Number of initial CAL OBs             : ', STRING(n_cal_tf, FORMAT='(I0)')
      PRINT, '  Number of preserved CAL OBs           : ', STRING(n_cal_tf_f, FORMAT='(I0)')
      PRINT, '  Number of initial CAL frames          : ', STRING(TOTAL(nfr_ob[idx_cal_tf]), FORMAT='(I0)') ;+ ' (' + STRING(TOTAL(sci_fro), FORMAT='(I0)') + ',' + STRING(TOTAL(cal_fro), FORMAT='(I0)') +')'
      PRINT, ' '
      IF lun GT 0 THEN BEGIN
        PRINTF, lun,  ' '
        PRINTF, lun,  '==== SUB-TF COMPUTATION (OB-based approach) ==== '
        PRINTF, lun,  '================================================ '
        PRINTF, lun,  ' '
        PRINTF, lun,  '==== PART ' + STRING(i_tf+1, FORMAT='(I0)') + ' ==== '
        PRINTF, lun,  '  Statistical error on TF [%] :'
        PRINTF, lun,  '    - based on individual error bars    : ', STRING(rms_stat_null*1D+2, FORMAT='(F6.4)')
        PRINTF, lun,  '    - standard deviation on the mean    : ', STRING(disp_m_null_we*1D+2, FORMAT='(F6.4)') + ' (' + STRING(disp_m_null*1D+2, FORMAT='(F6.4)') + ' if unweighted)'
        ;PRINTF, lun,  '    - quadratic sum of the 2            : ', STRING(null_err_tot*1D+2, FORMAT='(F6.4)')
        PRINTF, lun,  '  Systematic error on TF [%]            : '
        PRINTF, lun,  '    - due to diameter uncertainty       : ', STRING(null_err_diam*1D+2, FORMAT='(F6.4)')
        PRINTF, lun,  '    - sqrt(excess variance)             : ', STRING(rms_exc*1D+2, FORMAT='(F6.4)')
        PRINTF, lun,  '  Total error on TF [%]                 : ', STRING(null_err_tot*1D+2, FORMAT='(F6.4)')  
        PRINTF, lun,  '  Number of initial CAL OBs             : ', STRING(n_cal_tf, FORMAT='(I0)')
        PRINTF, lun,  '  Number of preserved CAL OBs           : ', STRING(n_cal_tf_f, FORMAT='(I0)')
        PRINTF, lun,  '  Number of initial CAL frames          : ', STRING(TOTAL(nfr_ob[idx_cal_tf]), FORMAT='(I0)') ;+ ' (' + STRING(TOTAL(sci_fro), FORMAT='(I0)') + ',' + STRING(TOTAL(cal_fro), FORMAT='(I0)') +')'
        PRINTF, lun,  ' '
      ENDIF
    ENDIF
    
    ; Now do the pointing approach (i.e., calibrate based on the pointing values
    ; **************************************************************************
    ; 
    ; Find pointings in this TF
    pt_id_tf   = pt_id[idx_tf]
    nod_pos_tf = nod_pos[idx_tf]
    pt_uniq_tf = pt_id_tf[UNIQ(pt_id_tf, SORT(pt_id_tf))] & n_pt_tf = N_ELEMENTS(pt_uniq_tf)
    idx_pt_tf  = VALUE_LOCATE(pt_uniq, pt_uniq_tf) 
    flag_pt_tf = flag_pt[idx_pt_tf]
    
    ; Prepare arrays
    time_avg_pt_tf = DBLARR(n_pt_tf, n_nod)
    null_avg_pt_tf = time_avg_pt_tf
    null_err_pt_tf = time_avg_pt_tf
    
    ; Recompute pointing-based value because they might be different from the ones computed above (because of specific TF splitting)
    FOR ip = 0, n_pt_tf-1 DO BEGIN
      ; Derive nod position for this TF
      idx_p       = WHERE(pt_id_tf EQ pt_uniq_tf[ip])  
      nod_pt_tf   = nod_pos_tf[idx_p]
      nod_uniq_tf = nod_pt_tf[UNIQ(nod_pt_tf, SORT(nod_pt_tf))]
      n_nod       = N_ELEMENTS(nod_uniq_tf) 
      FOR i_nod = 0, n_nod-1 DO BEGIN
        ; Extract nod data and parse to TF arrays (also ignore points to far from the TF)
        idx_pn =idx_tf[ WHERE(pt_id_tf EQ pt_uniq_tf[ip] AND nod_pos_tf EQ nod_uniq_tf[i_nod] AND ob_out[idx_tf] NE 1)]
        AVGERR, tf_wb[idx_pn], tf_wb_estat[idx_pn], avg_null, avg_null_we, rms_stat_null, disp_m_null, disp_m_null_we, rms_exc, KAPPA=drs.sig_sca, IDX_OUT=idx_out, N_OUT=n_out
        time_avg_pt_tf[ip,i_nod] = MEAN(time[idx_pn])
        IF drs.null_est EQ 0 THEN null_est = avg_null ELSE null_est = avg_null_we   
        null_avg_pt_tf[ip,i_nod] = null_est
        IF n_out GT 0 THEN ob_out[idx_pn[idx_out]]  = 1
        
        ; Systematic error (according to err_mode definition)
        CASE (ABS(drs.err_mode) MOD 3) OF
          0: rms_sys = disp_m_null
          1: rms_sys = disp_m_null_we
          2: rms_sys = rms_exc
          ELSE: MESSAGE, 'Undefined error mode'
        ENDCASE
        ; Compute total error as sum of the two or just the maximum one (keep error due to stellar diameter in both case)
        IF drs.err_mode GE 3 THEN null_err_pt_tf[ip,i_nod] = SQRT(rms_stat_null^2+MEAN(tf_wb_esyst[idx_pn])^2) > SQRT(MEAN(tf_wb_esyst[idx_pn])^2+rms_sys^2) ELSE null_err_pt_tf[ip,i_nod] = SQRT(rms_stat_null^2+MEAN(tf_wb_esyst[idx_pn])^2+rms_sys^2)
        ; Parse to global pointing arrays
        idx_pc = WHERE(pt_uniq EQ pt_uniq_tf[ip])
        time_avg_pt[idx_pc,i_nod] = time_avg_pt_tf[ip,i_nod]
        null_avg_pt[idx_pc,i_nod] = null_avg_pt_tf[ip,i_nod]
        null_err_pt[idx_pc,i_nod] = null_err_pt_tf[ip,i_nod]
      ENDFOR
    ENDFOR
    
    ; Find CAL and SCI pointings  
    idx_cal_pt_tf = WHERE(STRMATCH(flag_pt_tf, 'CAL') EQ 1, n_cal_pt_tf)
    idx_sci_pt_tf = WHERE(STRMATCH(flag_pt_tf, 'SCI') EQ 1, n_sci_pt_tf)
    
    ; Now loop over nods
    FOR i_nod = 0, n_nod-1 DO BEGIN
      ; Extract value for this nod
      time_avg_pt_tf_cal = time_avg_pt_tf[idx_cal_pt_tf,i_nod]    
      null_avg_pt_tf_cal = null_avg_pt_tf[idx_cal_pt_tf,i_nod] 
      null_err_pt_tf_cal = null_err_pt_tf[idx_cal_pt_tf,i_nod]
      
      ; Do interpolation on pointing nulls
      IF drs.cal_method EQ 0 THEN BEGIN
        ; To be implemented
        MESSAGE, 'Not yet implemented'
      ENDIF ELSE BEGIN  ; Calibation by polynomial interpolation on all calibrator data
        ; Interpolate the data at data and plot times
        polydeg = drs.polydeg
        IF polydeg GT n_cal THEN MESSAGE, 'Degree of polynomial fit should not exceed the number of calibrator data points.'
        polyfit = POLY_FIT(time_avg_pt_tf_cal, null_avg_pt_tf_cal, polydeg, /DOUBLE, MEASURE_ERRORS=null_err_pt_tf_cal^drs.null_est, CHISQ=chi2, COVAR=covar, SIGMA=polysig, YFIT=yfit, YBAND=yband, YERROR=yerror) ; include only statistical error bars 
        tf_data = 0D & tf_tmp = 0D
        FOR ipol=0, polydeg DO BEGIN
          tf_data += polyfit[ipol]*time_avg_pt_tf[idx_pt_tf,i_nod]^ipol
          tf_tmp  += polyfit[ipol]*plottime^ipol
        ENDFOR
        tf_plot_pt[i_tf,*,i_nod] = tf_tmp
      ENDELSE
      
      ; Calibrate null
      null_c_pt = null_avg_pt_tf[*,i_nod] - tf_data
  
      ; Extract calibrator data
      null_cal_pt     = null_c_pt[WHERE(STRMATCH(flag_pt_tf, 'CAL') EQ 1, ntmp)]
      null_cal_err_pt = null_err_pt_tf_cal
          
      ; Compute propagated statistical error on TF: POINTING BASED and renove large residues
      AVGERR, null_cal_pt, null_cal_err_pt, avg_null, avg_null_we, rms_stat_null_pt, disp_m_null_pt, disp_m_null_pt_we, rms_exc, KAPPA=drs.sig_sca, IDX_OUT=idx_out, N_OUT=n_out
                        
      ; Definition of systematic errors
      CASE (ABS(drs.err_mode) MOD 3) OF
        0: rms_sys_pt = disp_m_null_pt
        1: rms_sys_pt = disp_m_null_pt_we
        2: rms_sys_pt = rms_exc
        ELSE: MESSAGE, 'Undefined error mode'
      ENDCASE
     
      ; Systematic error is the same as the one per OB
       rms_sys_tot_pt = SQRT(rms_sys_pt^2 + null_err_diam^2) 
      
      ; Compute total error with systematic errors
      IF drs.err_mode GE 3 THEN null_err_tot_pt = SQRT(rms_stat_null_pt^2+null_err_diam^2)  > SQRT(rms_sys_pt^2+null_err_diam^2) ELSE null_err_tot_pt = SQRT(rms_stat_null_pt^2+rms_sys_tot_pt^2)
      tf_plot_pt_err[i_tf,i_nod] = null_err_tot_pt 
      
      ; Print diagnostic info to screen and log
      IF info GT 0 THEN BEGIN
        PRINT, ' '
        PRINT, '==== SUB-TF COMPUTATION (pointing-based approach) ==== '
        PRINT, '====================================================== '
        PRINT, ' '
        PRINT, '==== PART ' + STRING(i_tf+1, FORMAT='(I0)') + ' -- NOD ' + STRING(i_nod+1, FORMAT='(I0)') + ' ==== '
        PRINT, '  Statistical error on TF [%] :'
        PRINT, '    - based on individual error bars    : ', STRING(rms_stat_null_pt*1D+2, FORMAT='(F6.4)')
        PRINT, '    - standard deviation on the mean    : ', STRING(disp_m_null_pt_we*1D+2, FORMAT='(F6.4)') + ' (' + STRING(disp_m_null_pt*1D+2, FORMAT='(F6.4)') + ' if unweighted)'
        ;PRINT, '    - quadratic sum of the 2            : ', STRING(null_err_tot_pt*1D+2, FORMAT='(F6.4)')
        PRINT, '  Systematic error on TF [%]            : '
        PRINT, '    - due to diameter uncertainty       : ', STRING(null_err_diam*1D+2, FORMAT='(F6.4)')
        PRINT, '    - sqrt(excess variance)             : ', STRING(rms_sys_pt*1D+2, FORMAT='(F6.4)')
        PRINT, '  Total error on TF [%]                 : ', STRING(null_err_tot_pt*1D+2, FORMAT='(F6.4)')
        PRINT, '  Number of CAL pointings               : ', STRING(n_cal_pt_tf, FORMAT='(I0)')
        PRINT, ' '
        IF lun GT 0 THEN BEGIN
          PRINTF, lun, ' '
          PRINTF, lun, '==== SUB-TF COMPUTATION (pointing-based approach) ==== '
          PRINTF, lun, '====================================================== '
          PRINTF, lun, ' '
          PRINTF, lun, '==== PART ' + STRING(i_tf+1, FORMAT='(I0)') + ' -- NOD ' + STRING(i_nod+1, FORMAT='(I0)') + ' ==== '
          PRINTF, lun, '  Statistical error on TF [%] :'
          PRINTF, lun, '    - based on individual error bars    : ', STRING(rms_stat_null_pt*1D+2, FORMAT='(F6.4)')
          PRINTF, lun, '    - standard deviation on the mean    : ', STRING(disp_m_null_pt_we*1D+2, FORMAT='(F6.4)') + ' (' + STRING(disp_m_null_pt*1D+2, FORMAT='(F6.4)') + ' if unweighted)'
          ;PRINTF, lun, '    - quadratic sum of the 2            : ', STRING(null_err_tot_pt*1D+2, FORMAT='(F6.4)')
          PRINTF, lun, '  Systematic error on TF [%]            : '
          PRINTF, lun, '    - due to diameter uncertainty       : ', STRING(null_err_diam*1D+2, FORMAT='(F6.4)')
          PRINTF, lun, '    - sqrt(excess variance)             : ', STRING(rms_sys_pt*1D+2, FORMAT='(F6.4)')
          PRINTF, lun, '  Total error on TF [%]                 : ', STRING(null_err_tot_pt*1D+2, FORMAT='(F6.4)')
          PRINTF, lun, '  Number of CAL pointings               : ', STRING(n_cal_pt_tf, FORMAT='(I0)')
          PRINTF, lun,  ' '
        ENDIF
      ENDIF
      
      ; Ok, now save the calibrated null and null floor in output arrays if CAL_MODE is 1
      IF drs.cal_mode EQ 0 THEN BEGIN
        idx_tmp = WHERE(STRMATCH(flag_pt_tf, 'SCI') EQ 1)
        FOR i_s = 0, n_sci_pt_tf-1 DO BEGIN
          idx_s                       = WHERE(pt_uniq_tf[idx_tmp[i_s]] EQ sci_pt, n_ok) & IF n_ok LE 0 THEN MESSAGE, 'Problem with science pointing association!!'
          tf_sci[idx_s,i_nod]         = tf_data[idx_tmp[i_s]] 
          tf_estat[idx_s,i_nod]       = rms_stat_null_pt
          tf_esyst[idx_s,i_nod]       = rms_sys_tot_pt
          tf_etot[idx_s,i_nod]        = null_err_tot_pt  
          null_sci[idx_s,i_nod]       = null_c_pt[idx_tmp[i_s]]   
          null_sci_estat[idx_s,i_nod] = SQRT(null_err_pt_tf[idx_tmp[i_s],i_nod]^2+rms_stat_null_pt^2)
          null_sci_esyst[idx_s,i_nod] = rms_sys_tot_pt
          null_sci_etot[idx_s,i_nod]  = SQRT(null_sci_estat[idx_s]^2+null_sci_esyst[idx_s]^2)
        ENDFOR       
      ENDIF
    ENDFOR  
   ENDIF
ENDFOR

; PREPARE OUTPUT DATA (ALSO USED FOR PLOTS)
; *****************************************
IF drs.cal_mode EQ 1 THEN BEGIN
  idx_ok        = WHERE(null_sci NE 0, n_scidat)  ; 0 can happen when some SCI OB are not associated to any CAL OBs (e.g., if taken too far away from CAL OBs)
  sci_pid       = LONG(pid[idx_sci[idx_ok]])
  sci_name      = objname[idx_sci[idx_ok]]
  sci_time      = DOUBLE(scidata[idx_ok].mjd_obs)
  sci_utc       = utc[idx_sci[idx_ok]]
  sci_lst       = lst[idx_sci[idx_ok]]
  sci_ra        = ra[idx_sci[idx_ok]]
  sci_dec       = dec[idx_sci[idx_ok]]
  sci_alt       = FLOAT(alt[idx_sci[idx_ok]])
  sci_az        = FLOAT(az[idx_sci[idx_ok]])
  sci_para      = FLOAT(para[idx_sci[idx_ok]])
  sci_ha        = FLOAT(ha[idx_sci[idx_ok]])
  sci_ucoord    = FLOAT(u_coord[idx_sci[idx_ok]])
  sci_vcoord    = FLOAT(v_coord[idx_sci[idx_ok]])
  sci_wav       = FLOAT(wav[idx_sci[idx_ok]])
  sci_bdw       = FLOAT(bdw[idx_sci[idx_ok]])
  sci_aper      = FIX(aper[idx_sci[idx_ok]])
  sci_birad     = FIX(birad[idx_sci[idx_ok]])
  sci_borad     = FIX(borad[idx_sci[idx_ok]])
  sci_dit       = FLOAT(dit[idx_sci[idx_ok]])
  sci_fro       = LONG(nfr_ob[idx_sci[idx_ok]])
  sci_frr       = LONG(nfr_rej[idx_sci[idx_ok]])
  sci_seeing    = FLOAT(seeing[idx_sci[idx_ok]])
  sci_smttau    = FLOAT(smttau[idx_sci[idx_ok]])
  sci_fpcp      = FLOAT(fpcp[idx_sci[idx_ok]])
  qua_flag      = 2 + INTARR(n_scidat)  ; Quality flag (manual for now)
  null_meas     = DOUBLE(null[idx_sci[idx_ok]])
  null_meas_err = DOUBLE(null_err[idx_sci[idx_ok]])
  null_flo      = DOUBLE(tf_sci[idx_ok])
  null_flo_sta  = DOUBLE(tf_estat[idx_ok])
  null_flo_sys  = DOUBLE(tf_esyst[idx_ok])
  null_flo_err  = DOUBLE(tf_etot[idx_ok])
  null_cal      = DOUBLE(null_sci[idx_ok])
  null_cal_sta  = DOUBLE(null_sci_estat[idx_ok])
  null_cal_sys  = DOUBLE(null_sci_esyst[idx_ok])
  null_cal_err  = DOUBLE(null_sci_etot[idx_ok])
ENDIF ELSE BEGIN
  sci_pid  = LONARR(n_scidat_pt,n_nod) & sci_name = STRARR(n_scidat_pt,n_nod) & sci_time = DBLARR(n_scidat_pt,n_nod) & sci_utc = STRARR(n_scidat_pt,n_nod) & sci_lst  = STRARR(n_scidat_pt,n_nod)
  sci_ra   = STRARR(n_scidat_pt,n_nod) & sci_dec  = STRARR(n_scidat_pt,n_nod) & sci_alt  = FLTARR(n_scidat_pt,n_nod) & sci_az  = FLTARR(n_scidat_pt,n_nod) & sci_para = FLTARR(n_scidat_pt,n_nod)
  sci_ha   = FLTARR(n_scidat_pt,n_nod) & sci_wav  = FLTARR(n_scidat_pt,n_nod) & sci_bdw  = FLTARR(n_scidat_pt,n_nod) & sci_aper= INTARR(n_scidat_pt,n_nod) & sci_birad= INTARR(n_scidat_pt,n_nod)
  sci_borad= INTARR(n_scidat_pt,n_nod) & sci_dit  = FLTARR(n_scidat_pt,n_nod) & sci_fro  = LONARR(n_scidat_pt,n_nod) & sci_frr = LONARR(n_scidat_pt,n_nod) & sci_fpcp = FLTARR(n_scidat_pt,n_nod) 
  sci_seeing = FLTARR(n_scidat_pt,n_nod) & sci_smttau = FLTARR(n_scidat_pt,n_nod) 
  sci_ucoord = FLTARR(n_scidat_pt,n_nod) & sci_vcoord = FLTARR(n_scidat_pt,n_nod)
  qua_flag = 2 + INTARR(n_scidat_pt)
  FOR i_nod = 0, n_nod-1 DO BEGIN
    FOR i_pt = 0, n_scidat_pt-1 DO BEGIN
      idx_pt           = WHERE(pt_id EQ pt_uniq[idx_sci_pt[i_pt]] AND nod_pos EQ nod_uniq[i_nod], n_inpt)
      sci_pid[i_pt,i_nod]    = LONG(pid[idx_pt[0]])
      sci_name[i_pt,i_nod]   = objname[idx_pt[0]]
      sci_time[i_pt,i_nod]   = DOUBLE(MEAN(l1data[idx_pt[0]].mjd_obs))
      sci_utc[i_pt,i_nod]    = utc[idx_pt[n_inpt/2]]
      sci_lst[i_pt,i_nod]    = lst[idx_pt[n_inpt/2]]
      sci_ra[i_pt,i_nod]     = ra[idx_pt[n_inpt/2]]
      sci_dec[i_pt,i_nod]    = dec[idx_pt[n_inpt/2]]
      sci_alt[i_pt,i_nod]    = FLOAT(MEAN(alt[idx_pt]))
      sci_az[i_pt,i_nod]     = FLOAT(MEAN(az[idx_pt]))
      sci_para[i_pt,i_nod]   = FLOAT(MEAN(para[idx_pt]))
      sci_ha[i_pt,i_nod]     = FLOAT(MEAN(ha[idx_pt]))
      sci_ucoord[i_pt,i_nod] = FLOAT(MEAN(u_coord[idx_pt]))
      sci_vcoord[i_pt,i_nod] = FLOAT(MEAN(v_coord[idx_pt]))
      sci_wav[i_pt,i_nod]    = FLOAT(wav[idx_pt[0]])
      sci_bdw[i_pt,i_nod]    = FLOAT(bdw[idx_pt[0]])
      sci_aper[i_pt,i_nod]   = FIX(aper[idx_pt[0]])
      sci_birad[i_pt,i_nod]  = FIX(birad[idx_pt[0]]) 
      sci_borad[i_pt,i_nod]  = FIX(borad[idx_pt[0]]) 
      sci_dit[i_pt,i_nod]    = FLOAT(dit[idx_pt[0]])
      sci_fro[i_pt,i_nod]    = LONG(TOTAL(nfr_ob[idx_pt]))
      sci_frr[i_pt,i_nod]    = LONG(TOTAL(nfr_rej[idx_pt]))
      sci_seeing[i_pt,i_nod] = FLOAT(MEAN(seeing[idx_pt]))
      sci_smttau[i_pt,i_nod] = FLOAT(MEAN(smttau[idx_pt]))
      sci_fpcp[i_pt,i_nod]   = FLOAT(MEAN(fpcp[idx_pt]))
  ENDFOR
  ENDFOR  
  null_meas     = DOUBLE(REFORM(null_avg_pt[idx_sci_pt,*],1,n_scidat_pt*n_nod))
  null_meas_err = DOUBLE(REFORM(null_err_pt[idx_sci_pt,*],1,n_scidat_pt*n_nod))
  sci_pid       = REFORM(sci_pid,1,n_scidat_pt*n_nod)
  sci_name      = REFORM(sci_name,1,n_scidat_pt*n_nod)
  sci_time      = REFORM(sci_time,1,n_scidat_pt*n_nod)
  sci_utc       = REFORM(sci_utc,1,n_scidat_pt*n_nod)
  sci_lst       = REFORM(sci_lst,1,n_scidat_pt*n_nod)
  sci_ra        = REFORM(sci_ra,1,n_scidat_pt*n_nod)
  sci_dec       = REFORM(sci_dec,1,n_scidat_pt*n_nod)
  sci_alt       = REFORM(sci_alt,1,n_scidat_pt*n_nod)
  sci_az        = REFORM(sci_az,1,n_scidat_pt*n_nod)
  sci_para      = REFORM(sci_para,1,n_scidat_pt*n_nod)
  sci_ha        = REFORM(sci_ha,1,n_scidat_pt*n_nod)
  sci_wav       = REFORM(sci_wav,n_scidat_pt*n_nod)
  sci_bdw       = REFORM(sci_bdw,1,n_scidat_pt*n_nod)
  sci_dit       = REFORM(sci_dit,1,n_scidat_pt*n_nod)
  sci_fro       = REFORM(sci_fro,1,n_scidat_pt*n_nod)
  sci_frr       = REFORM(sci_frr,1,n_scidat_pt*n_nod)
  sci_seeing    = REFORM(sci_seeing,1,n_scidat_pt*n_nod)
  sci_smttau    = REFORM(sci_smttau,1,n_scidat_pt*n_nod)
  sci_fpcp      = REFORM(sci_fpcp,1,n_scidat_pt*n_nod)
  null_flo      = DOUBLE(REFORM(tf_sci,1,n_scidat_pt*n_nod))
  null_flo_sta  = DOUBLE(REFORM(tf_estat,1,n_scidat_pt*n_nod))
  null_flo_sys  = DOUBLE(REFORM(tf_esyst,1,n_scidat_pt*n_nod))
  null_flo_err  = DOUBLE(REFORM(tf_etot,1,n_scidat_pt*n_nod))
  null_cal      = DOUBLE(REFORM(null_sci,1,n_scidat_pt*n_nod))
  null_cal_sta  = DOUBLE(REFORM(null_sci_estat,1,n_scidat_pt*n_nod))
  null_cal_sys  = DOUBLE(REFORM(null_sci_esyst,1,n_scidat_pt*n_nod))
  null_cal_err  = DOUBLE(REFORM(null_sci_etot,1,n_scidat_pt*n_nod))
ENDELSE
          
; ESTIMATE TARGET FLUX BASED ON CALIBRATOR FLUX
; *********************************************

IF n_scidat GT 0 AND n_caldat GT 1 THEN BEGIN
  ; Individual photometry
  cal_adu_dx_avg = l1data[idx_cal].photdx_avg & cal_adu_dx_rms = l1data[idx_cal].photdx_rms
  cal_adu_sx_avg = l1data[idx_cal].photsx_avg & cal_adu_sx_rms = l1data[idx_cal].photsx_rms
  sci_adu_dx_avg = l1data[idx_sci].photdx_avg & sci_adu_dx_rms = l1data[idx_sci].photdx_rms
  sci_adu_sx_avg = l1data[idx_sci].photsx_avg & sci_adu_sx_rms = l1data[idx_sci].photsx_rms
  ; Combined photometry
  cal_adu_avg = 2*(cal_adu_dx_avg+cal_adu_sx_avg) & cal_adu_rms = SQRT(cal_adu_dx_rms^2+cal_adu_sx_rms^2)
  sci_adu_avg = 2*(sci_adu_dx_avg+sci_adu_sx_avg) & sci_adu_rms = SQRT(sci_adu_dx_rms^2+sci_adu_sx_rms^2)
  ; Compute estimated Jy to ADU conversion factor based on estimated calibrator flux
  jy2adu_avg   = tgt_flx[idx_cal]/cal_adu_avg & jy2adu_rms = tgt_flx[idx_cal]/cal_adu_avg^2*cal_adu_rms
  jy2adu_estat = 1./SQRT(TOTAL(1./jy2adu_rms^2))
  ; Compute mean over all measurements
  AVGSDV, jy2adu_avg, jy2adu, jy2adu_err, jy2adu_err_mean, KAPPA=5
  ; Number of different science objects
  sci_uniq = sci_name(UNIQ(sci_name,SORT(sci_name))) & n_sci = N_ELEMENTS(sci_uniq)
  IF info GT 0 THEN BEGIN
    PRINT, '==== PHOTOMETRY INFO ===='
    PRINT, 'Jy to ADU conversion               : ', STRING(jy2adu, FORMAT='(F8.6)')
    PRINT, 'Science target estimated flux [Jy] : '
    FOR i=0, n_sci-1 DO BEGIN
      idx_tmp = WHERE(sci_name EQ sci_uniq[i])
      PRINT, '    - ' + sci_uniq[i] + '         : ' + STRING(MEAN(jy2adu*sci_adu_avg[idx_tmp]), FORMAT='(F5.1)')
    ENDFOR
    PRINT, ' '
  ENDIF
ENDIF

IF info GT 0 AND FXPAR(header, 'NUL_MODE') EQ 2 THEN BEGIN
  PRINT, '==== NSC-RELATED INFO ===='
  PRINT, 'Mean phase offset [rad]               : ', STRING(MEAN(l1data.nsc_phavg), FORMAT='(F8.6)')
  PRINT, 'Mean error on mean phase offset [rad] : ', STRING(MEAN(l1data.nsc_phavg_err), FORMAT='(F8.6)')
  PRINT, 'Phase jitter [rad]                    : ', STRING(MEAN(l1data.nsc_phrms), FORMAT='(F8.6)')
  PRINT, 'Mean error on phase jitter [rad]      : ', STRING(MEAN(l1data.nsc_phrms_err), FORMAT='(F8.6)')
  PRINT, ' '
  IF lun GT 0 THEN BEGIN
    PRINTF, lun,  '==== NSC-RELATED INFO ===='
    PRINTF, lun, 'Mean phase offset for calibrator [rad]               : ', STRING(MEAN(l1data.nsc_phavg), FORMAT='(F8.6)')
    PRINTF, lun, 'Mean error on mean phase offset for calibrator [rad] : ', STRING(MEAN(l1data.nsc_phavg_err), FORMAT='(F8.6)')
    PRINTF, lun, 'Phase jitter for calibrator [rad]                    : ', STRING(MEAN(l1data.nsc_phrms), FORMAT='(F8.6)')
    PRINTF, lun, 'Mean error on phase jitter for calibrator [rad]      : ', STRING(MEAN(l1data.nsc_phrms_err), FORMAT='(F8.6)')
    PRINTF, lun, ' '
  ENDIF
ENDIF

; SAVE L2 FITS FILE
; *****************

IF NOT KEYWORD_SET(NO_SAVE) THEN BEGIN
  ; Save results to an oifits file
  outfile = l2_dir + 'UT' + date_lng + '_calib_' + id_name + '.fits'

  ; Create new FITS file with primary HDU
  FXHMAKE,  hdr, /INIT, /EXTEND, 0
  FXADDPAR, hdr, "COMMENT", "This is a FITS file of calibrated (level 2) data taken with LBTI/" + STRTRIM(FXPAR(header, 'INSTRUME')) + "."
  FXADDPAR, hdr, "COMMENT", "The raw data has been reduced using the nodrs software version " +  STRING(FXPAR(header, 'DRS_VERS'), FORMAT='(F3.1)') + "."
  FXADDPAR, hdr, "COMMENT", "Contact: lbtisupp@ipac.caltech.edu"
  FXADDPAR, hdr, "COMMENT", "Date and DRS information", BEFORE='DATE_OBS'
  FXADDPAR, hdr, 'DATE_OBS',  date_lng,                            'Date of observation'
  FXADDPAR, hdr, 'TELESCOP', 'LBT',                                'Telescope'
  FXADDPAR, hdr, 'INSTRUME',  FXPAR(header, 'INSTRUME'),           'Instrument'
  FXADDPAR, hdr, 'DRS_VERS',  FXPAR(header, 'DRS_VERS'),           'Version of the DRS'
  FXADDPAR, hdr, 'DRS_DATE',  FXPAR(header, 'DRS_DATE'),           'DRS version date'
  FXADDPAR, hdr, 'DATE_RED',  FXPAR(header, 'DATE_RED'),           'Date of data reduction'
  FXADDPAR, hdr, 'DATE_CAL',  systime(),                           'Date of data calibration'
  FXADDPAR, hdr, "COMMENT", "Image calibration parameters", AFTER='DATE_CAL'
  FXADDPAR, hdr, 'BCK_MODE',  FXPAR(header, 'BCK_MODE'),           '0=nod pairs, 1=0, 2=adjacent nods, 3=closest frames''
  FXADDPAR, hdr, 'N_BIN',     FXPAR(header, 'N_BIN'),              'Image re-binning size'
  FXADDPAR, hdr, 'NO_BPM',    FXPAR(header, 'NO_BPM'),             'Bad pixel correction off if 1 '
  FXADDPAR, hdr, 'NO_DARK',   FXPAR(header, 'NO_DARK'),            'Dark subtraction off if 1 '
  FXADDPAR, hdr, 'NO_FLAT',   FXPAR(header, 'NO_FLAT'),            'Flat correction off if 1'
  FXADDPAR, hdr, 'NOD_FRQ',   FXPAR(header, 'NOD_FRQ'),            'Mean nodding period [s]'
  FXADDPAR, hdr, 'N_FRBCK',   FXPAR(header, 'N_FRBCK'),            'Number of preserved frames for background subtraction'
  FXADDPAR, hdr, "COMMENT", "Flux computation parameters", AFTER='N_FRBCK'
  FXADDPAR, hdr, 'FLX_MODE',  FXPAR(header, 'FLX_MODE'),           '0=aperture phot, 1=PSF fitting'
  FXADDPAR, hdr, 'FIT_MODE',  FXPAR(header, 'FIT_MODE'),           '0=no centroid, 1-3=Gaussian fit, 4=Lorentzian fit, 5=Moffat fit'
  FXADDPAR, hdr, "COMMENT", "Null computation parameters", AFTER='FIT_MODE'
  FXADDPAR, hdr, 'OB_MODE',   FXPAR(header, 'OB_MODE'),            '0=all, 1=nod splitted, 2=user defined'
  FXADDPAR, hdr, 'NUL_MODE',  FXPAR(header, 'NUL_MODE'),           '0=mode, 1=%best, 2=statistical'
  FXADDPAR, hdr, 'R_FRNULL',  FXPAR(header, 'R_FRNULL'),           'Ratio of preserved frames for null computation'
  FXADDPAR, hdr, 'NULL_COR',  FXPAR(header, 'NULL_COR'),           'Subtract the null estimated from high-frequency phase noise'
  FXADDPAR, hdr, 'N_BTSTRP',  FXPAR(header, 'N_BTSTRP'),           'Number of bootstrap samples to compute the error bar'
  FXADDPAR, hdr, 'NSC_CUBE',  nsc_cube,                            'Cube size for NSC reduction (null, phase offset, phase rms)'
  FXADDPAR, hdr, 'NSC_BFAC',  nbin_fac,                            'Multiplier on number of histogram bins (nbins_fac*sqrt(Np))
  FXADDPAR, hdr, 'NSC_BINS' , nsc_bins,                            'Bin size for NSC reduction (0: constant, 1: variable)'
  FXADDPAR, hdr, 'NSC_OMIN',  nsc_omin,                            'Minimum number of occurences per bin for the fit'
  FXADDPAR, hdr, 'NULL_COR',  null_cor,                            'Correct for high frequency phase noise to each null measurement'
  FXADDPAR, hdr, 'NULL_LIM',  null_lim,                            'Acceptable raw null range (before NSC)'
  FXADDPAR, hdr, "COMMENT", "Null calibration parameters", AFTER='NULL_LIM'
  FXADDPAR, hdr, 'CAL_METH',  FIX(drs.cal_method),                 '0=linear interpolation between two nearest neighbors, 1=global polynomial interpolation'
  FXADDPAR, hdr, 'CAL_MODE',  FIX(drs.cal_mode),                   '0=calibrate per pointing, 1=calibrate per OB'
  FXADDPAR, hdr, 'CHI2_LIM',  FIX(drs.cal_mode),                   'Limit on acceptable chi2 (only with NSC)'
  FXADDPAR, hdr, 'NULL_EST',  FIX(drs.null_est),                   '0=unweighted, 1=weighted average of OBs'
  FXADDPAR, hdr, 'POLYDEG',   FIX(drs.polydeg),                    'Degree of the polynomial for method CAL_MODE=1'
  FXADDPAR, hdr, 'SIG_SCA',   FIX(drs.sig_sca),                    'Number of sigmas used for OB-based scatter OB removal'
  FXADDPAR, hdr, 'SPLT_NOD',  FIX(drs.split_nod),                  '0=calibrate both nod positions at once, 1=calibrate separately'
  FXADDPAR, hdr, 'SPLT_TF' ,  FIX(drs.split_tf),                   '1=split the TF where dead time > time_split'
  FXADDPAR, hdr, 'TIME_SPL',  drs.split_time,                      'Maximum time between 2 nulls for a single TF [hours]'
  FXWRITE, outfile, hdr
  
  ; Create header for the main table
  n_row = N_ELEMENTS(sci_wav)
  FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'RESULTS_SUMMARY', 'results of data calibration'

  ; Init column number
  n_col = 35
  col   = LINDGEN(n_col)+1L

  ; Fill extension header with column names
  FXBADDCOL, 1L, hdr, LONG(0)               , 'PID',          'Project ID number'
  FXBADDCOL, 2L, hdr, 'aaaaaaaaaaaaaaaaaaaa', 'OBJNAME',      'Object name'
  FXBADDCOL, 3L, hdr, DOUBLE(0)             , 'MJD_OBS',      'Modified Julian Date of observation'
  FXBADDCOL, 4L, hdr, '00:00:00.00'         , 'LBT_UTC',      'UTC from observatory'
  FXBADDCOL, 5L, hdr, '00:00:00.00'         , 'LBT_LST',      'LST from observatory'
  FXBADDCOL, 6L, hdr, '+00:00:00.000'       , 'LBT_RA',       'RA from observatory'
  FXBADDCOL, 7L, hdr, '+00:00:00.000'       , 'LBT_DEC',      'DEC from observatory'
  FXBADDCOL, 8L, hdr, FLOAT(0)              , 'LBT_ALT',      'ALT from observatory'
  FXBADDCOL, 9L, hdr, FLOAT(0)              , 'LBT_AZ',       'AZ from observatory'
  FXBADDCOL, 10L, hdr, FLOAT(0)             , 'LBT_PARA',     'Parralactic angle from obs.'
  FXBADDCOL, 11L, hdr, FLOAT(0)             , 'HR_ANGLE',     'Hour angle from obs.'
  FXBADDCOL, 12L, hdr, FLOAT(0)             , 'U_COORD',      'u coordinate', TUNIT='cycles/arcsec'
  FXBADDCOL, 13L, hdr, FLOAT(0)             , 'V_COORD',      'v coordinate', TUNIT='cycles/arcsec'
  FXBADDCOL, 14L, hdr, FLOAT(0)             , 'INT_TIME',     'Integration time', TUNIT='s'
  FXBADDCOL, 15L, hdr, FLOAT(0)             , 'WAV_EFF',      'Central wavelength', TUNIT='m'
  FXBADDCOL, 16L, hdr, FLOAT(0)             , 'BANDWIDTH',    'Bandwidth', TUNIT='m'
  FXBADDCOL, 17L, hdr, FIX(0)               , 'APER_RAD',     'Radius for aperture photometry', TUNIT='pix'
  FXBADDCOL, 18L, hdr, FIX(0)               , 'BCK_IRAD',     'Inner radius for background computation', TUNIT='pix'
  FXBADDCOL, 19L, hdr, FIX(0)               , 'BCK_ORAD',     'Outer radius for background computation', TUNIT='pix'
  FXBADDCOL, 20L, hdr, DOUBLE(0)            , 'NULL_MEAS',    'Measured raw null'
  FXBADDCOL, 21L, hdr, DOUBLE(0)            , 'NULL_MEAS_ERR','Error on the measured raw null'
  FXBADDCOL, 22L, hdr, DOUBLE(0)            , 'NULL_FLO',     'Null floor at SCI position'
  FXBADDCOL, 23L, hdr, DOUBLE(0)            , 'NULL_FLO_STA', 'Statistical error on null floor'  
  FXBADDCOL, 24L, hdr, DOUBLE(0)            , 'NULL_FLO_SYS', 'Systematic error on null floor'  
  FXBADDCOL, 25L, hdr, DOUBLE(0)            , 'NULL_FLO_ERR', 'Total error on null floor'  
  FXBADDCOL, 26L, hdr, DOUBLE(0)            , 'NULL_CAL',     'Calibrated null'
  FXBADDCOL, 27L, hdr, DOUBLE(0)            , 'NULL_CAL_STA', 'Statistical error on the calibrated null'
  FXBADDCOL, 28L, hdr, DOUBLE(0)            , 'NULL_CAL_SYS', 'Systematic error on the calibrated null'
  FXBADDCOL, 29L, hdr, DOUBLE(0)            , 'NULL_CAL_ERR', 'Total error on the calibrated null'
  FXBADDCOL, 30L, hdr, LONG(0)              , 'NFR_OB',       'Number of null frames'
  FXBADDCOL, 31L, hdr, LONG(0)              , 'NFR_REJ',      'Number of rejected null frames'
  FXBADDCOL, 32L, hdr, FLOAT(0)             , 'SEEING',       'Median seeing', TUNIT='arcsec'
  FXBADDCOL, 33L, hdr, FLOAT(0)             , 'SMTTAU',       'Precipitable water vapor estimate', TUNIT='mm'
  FXBADDCOL, 34L, hdr, FLOAT(0)             , 'FPC_PISTS',    'FPC RMS piston', TUNIT='um'
  FXBADDCOL, 35L, hdr, FIX(0)               , 'QUALITY_FLG',  'Quality flag'

  ; Write extension header to FITS file
  FXBCREATE, unit, outfile, hdr
  FXBWRITM,  unit, col, sci_pid, sci_name, sci_time, sci_utc, sci_lst, sci_ra, sci_dec, sci_alt, sci_az, sci_para, sci_ha, sci_ucoord, sci_vcoord, sci_dit, sci_wav, sci_bdw, $
                        sci_aper, sci_birad, sci_borad, null_meas, null_meas_err, null_flo, null_flo_sta, null_flo_sys, null_flo_err, null_cal, null_cal_sta, null_cal_sys, null_cal_err, sci_fro, $
                        sci_frr, sci_seeing, sci_smttau, sci_fpcp, qua_flag
  FXBFINISH, unit
  
  ; Copy file to new version
  FILE_COPY, outfile, l2_dir + 'UT' + date_lng + '_calib_' + id_name + '_aper' + STRING(aper_rad, FORMAT='(I0)') + '_L1v' + STRING(file_ver, FORMAT='(I0)') + '.fits', /OVERWRITE
ENDIF

; Jumping point if no science data
SKIP_INTERP:
IF n_scidat EQ 0 THEN BEGIN
  MESSAGE, 'WARNING, no science target in the input files', /CONTINUE
  RETURN
ENDIF

; Now start plotting....lots and lots of plots
; ********************************************

  
  IF KEYWORD_SET(RUNBIAS) THEN bias_tag = '_BIAS' ELSE bias_tag = ''
  ; 1. Plot the uv plane
  ; --------------------
  
  ;charthick = 4.0
  ;charsize  = 1.3
  uv_path  = pth.result_path + 'uv' + pth.sep
  IF NOT FILE_TEST(uv_path) THEN FILE_MKDIR, uv_path
  
  PREP_PS, /BOLD & DEVICE, FILENAME= uv_path + date_lng + '_uv_coord.eps', /ENCAPSULATE, /COLOR, XSIZE=15.4, YSIZE=14.7, /TIMES
  xrange = -[-10,10] & yrange = -xrange
  PLOT, [0,0], [0,0], XTITLE='u [cycles/arcsec]', YTITLE='v [cycles/arcsec]', TITLE='uv coordinates', $;',$ ;STRTRIM(sci_name[0],2), $
        XRANGE=xrange, YRANGE=yrange, XSTYLE=1, YSTYLE=1, /NODATA
  ;PLOTS, [MIN(xrange), MAX(xrange)], [0,0], LINESTYLE=2
  ;PLOTS, [0,0], [MIN(yrange), MAX(yrange)], LINESTYLE=2
  
  ; Then plot the good data
  LOADCT, 13, /SILENT
  FOR ik = 0, n_scidat-1 DO OPLOT, [-u_coord[idx_sci[ik]],u_coord[idx_sci[ik]]], [-v_coord[idx_sci[ik]],v_coord[idx_sci[ik]]], PSYM=1, COLOR=250
  ;FOR ik = 0, n_caldat-1 DO OPLOT, [-u_coord[idx_cal[ik]],u_coord[idx_cal[ik]]], [-v_coord[idx_cal[ik]],v_coord[idx_cal[ik]]], PSYM=1, COLOR=90
  
  ; Overplot the orientation of the main disk (+ 90 because pa is given for the photosphere)
  disk_ori = (114) * !Dpi / 180D0
  x = xrange[0] + (-xrange[0] + xrange[1]) * DINDGEN(1D+2)/(1D+2-1)
  y = TAN(disk_ori) * x
  OPLOT, y, x, LINESTYLE=2, THICK = 4          ; PA counted from North to (y,x) and not (x,y)
  
  ; Draw E-N small axis at the bootm right of the figure
  ARROW, 0.90*xrange[1], 0.90*yrange[0], 0.65*xrange[1], 0.90*yrange[0], THICK=5.5, HEAD_INDENT=0.1, /DATA
  ARROW, 0.90*xrange[1], 0.90*yrange[0], 0.90*xrange[1], 0.65*yrange[0], THICK=5.5, HEAD_INDENT=0.1, /DATA
  XYOUTS, 0.63*xrange[1], 0.87*yrange[0], "E", CHARTHICK=charthick, CHARSIZE=charsize
  XYOUTS, 0.85*xrange[1], 0.64*yrange[0], "N", CHARTHICK=charthick, CHARSIZE=charsize
  ;XYOUTS, 1.2, -2.4, "Moerchen midplane", ORIENTATION=70
  DEVICE, /CLOSE & END_PS
  
  ; 2. Plot data and transfer function (OB based)
  ; ---------------------------------------------
 
  ; Plot path (create a specific path if DIR_LABEL is set)
  ; The label must not constain _APR, which is treated differently
  IF STRMATCH(drs.dir_label, '*_APR') THEN label = STRMID(drs.dir_label, 0, STRPOS(drs.dir_label, '_APR')) $
                                      ELSE label = drs.dir_label
  ; If label is not empty, add label to the TF path
  IF STRLEN(STRCOMPRESS(label, /REMOVE_ALL)) EQ 0 THEN tf_path  = pth.result_path + 'TF' + pth.sep + date_lng + pth.sep $
                                                  ELSE tf_path  = pth.result_path + 'TF' + pth.sep + date_lng + pth.sep + label + pth.sep  
  
  IF NOT FILE_TEST(tf_path) THEN FILE_MKDIR, tf_path
  plotname =  tf_path + date_lng + '_TF-OB_APER' + STRING(aper_rad, FORMAT='(I0)') + bias_tag + '.eps'
  ; Init plot
  PREP_PS, /BOLD
  DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=17.8, YSIZE=13.2, /TIMES
  IF KEYWORD_SET(NO_INSET) THEN BEGIN
    !P.MULTI = [0,1,1]
    xtitle   = 'UT hour'
  ENDIF ELSE BEGIN
    !P.MULTI  = [0,1,2]
    xtickname = [' ',' ',' ',' ',' ',' ',' ',' ',' ']
    xtitle    = ' '
    postop    = [0.151,0.28,0.954,0.932]
    posbot1   = [0.151,0.10,0.954,0.28]
  ENDELSE
  LOADCT, 12, /SILENT
  fac        = 24D
  min_time   = 0.98*MIN(time)
  max_time   = 1.02*MAX(time)
  null_plot  = [-1,5]
  xrange     = [min_time,max_time]*fac
  null_range = MAX(null)-MIN(null)
  yrange     = [MIN(null)-0.10*null_range,MAX(null)+0.10*null_range]*1D2
  yrange[0]  = yrange[0] < null_plot[0]
  yrange[1]  = yrange[1] > null_plot[1]
  PLOT, [0], [0], XTITLE=xtitle, YTITLE='Measured null depth per OB [%]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
        POSITION=postop, XTICKNAME=xtickname, XRANGE=xrange, YRANGE=yrange
  ; Overplot the estimated null values for calibrators measurements
  IF n_caldat GT 0 THEN BEGIN
    USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.7, THICK=1.5
    OPLOT, caltime*fac, 100.*tf_wb[idx_cal], PSYM=8, COLOR=100
    ERRPLOT, caltime*fac, 100.*(tf_wb[idx_cal]-tf_wb_etot[idx_cal]), 100.*(tf_wb[idx_cal]+tf_wb_etot[idx_cal]), COLOR=100
    cal_null = null[idx_cal] & cal_nerr = null_err[idx_cal]
    USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.7, THICK=1.5, /FILL
    OPLOT, caltime*fac, 100.*cal_null, PSYM=8, COLOR=100
    ERRPLOT, caltime*fac, 100.*(cal_null-cal_nerr), 100.*(cal_null+cal_nerr), COLOR=100
  ENDIF
  ; Overplot the estimated null values for science measurements
  IF n_scidat GE 1 THEN BEGIN
    sci_null = null[idx_sci] & sci_nerr = null_err[idx_sci]
    USERSYM, [1,0,-1,0,1], [0,1,0,-1,0], THICK=1.5, /FILL
    OPLOT, scitime*fac, 100.*sci_null, PSYM=8, COLOR=200
    ERRPLOT, scitime*fac, 100.*(sci_null-sci_nerr), 100.*(sci_null+sci_nerr), COLOR=200
  ENDIF
  
  ; Overplot the estimated T2 values at science measurements (and confidence intervals)
  LOADCT, 0, /SILENT
  IF drs.cal_method EQ 1 AND n_scidat NE 0 THEN BEGIN
    FOR i_tf = 0, n_tf-1 DO BEGIN
      OPLOT, REFORM(tf_time[i_tf,*])*fac, 100.*REFORM(tf_plot[i_tf,*]), LINESTYLE=0, THICK=4
      OPLOT, REFORM(tf_time[i_tf,*])*fac, 100.*(REFORM(tf_plot[i_tf,*])-tf_plot_err[i_tf]), LINESTYLE=1, THICK=3.5     
      OPLOT, REFORM(tf_time[i_tf,*])*fac, 100.*(REFORM(tf_plot[i_tf,*])+tf_plot_err[i_tf]), LINESTYLE=1, THICK=3.5
    ENDFOR  
  ENDIF
  
  ; Plot background bias in the bottom inset if requested
  no_inset = 1
  IF NOT KEYWORD_SET(NO_INSET) THEN BEGIN
    yrange = [-0.004,0.004]*1D2 + avg_bias*1D2
    PLOT, [0], [0], XTITLE='UT hour', YTITLE='Bckg null [%]', TITLE=' ', XSTYLE=1, YSTYLE=1, $
          POSITION=posbot1, XRANGE=xrange, YRANGE=yrange, YMINOR=3
    LOADCT, 12, /SILENT
    ; Overplot the estimated null values for calibrators measurements
    IF n_caldat GT 0 THEN BEGIN
      USERSYM, [1,1,-1,-1,1]*0.5, [1,-1,-1,1,1]*0.5, THICK=1.5, /FILL
      OPLOT, caltime*fac, 1D2*bias[idx_cal], PSYM=8, COLOR=100
      ERRPLOT, caltime*fac, 1D2*(bias[idx_cal]-bias_err[idx_cal]), 1D2*(bias[idx_cal]+bias_err[idx_cal]), COLOR=100
    ENDIF
    IF n_scidat GT 0 THEN BEGIN
      USERSYM, [1,0,-1,0,1]*0.5, [0,1,0,-1,0]*0.5, THICK=1.5, /FILL
      OPLOT, scitime*fac, 1D2*bias[idx_sci], PSYM=8, COLOR=200
      ERRPLOT, scitime*fac, 1D2*(bias[idx_sci]-bias_err[idx_sci]), 1D2*(bias[idx_sci]+bias_err[idx_sci]), COLOR=200
    ENDIF
    AVGSDV, bias, avg_bias_we, rms_bias, rms_bias_mean, KAPPA=5, WEIGHT=1./bias_err^2
    OPLOT, [0,24], [avg_bias_we,avg_bias_we]*1D+2, LINESTYLE=0, THICK=4
    OPLOT, [0,24], [avg_bias_we+rms_bias_mean,avg_bias_we+rms_bias_mean]*1D2, LINESTYLE=1, THICK=4
    OPLOT, [0,24], [avg_bias_we-rms_bias_mean,avg_bias_we-rms_bias_mean]*1D2, LINESTYLE=1, THICK=4
  ENDIF
    
  ; Finish plot
  DEVICE, /CLOSE & END_PS
  !P.MULTI=[0,1,1]

  ; 2bis. Plot filtered data and transfer function (OB based)
  ; --------------------------------------------------------
  
  ; Filter data
  idx_ok       = WHERE(ob_out NE 1)
  time_f       = time[idx_ok]
  null_f       = null[idx_ok]
  null_err_f   = null_err[idx_ok]
  tf_wb_f      = tf_wb[idx_ok]
  tf_wb_etot_f = tf_wb_etot[idx_ok]
  flag_f       = flag[idx_ok]
  idx_cal_f    = WHERE(STRMATCH(flag_f, 'CAL') EQ 1, n_caldat_f)
  idx_sci_f    = WHERE(STRMATCH(flag_f, 'SCI') EQ 1, n_scidat_f)
  caltime_f    = time_f[idx_cal_f]
  scitime_f    = time_f[idx_sci_f]
  
  ; Plot path
  plotname =  tf_path + date_lng + '_TF-OB_APER' + STRING(aper_rad, FORMAT='(I0)') + bias_tag + '_FILT.eps'
  ; Init plot
  PREP_PS, /BOLD
  DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=17.8, YSIZE=13.2, /TIMES
  IF KEYWORD_SET(NO_INSET) THEN BEGIN
    !P.MULTI = [0,1,1]
    xtitle   = 'UT hour'
  ENDIF ELSE BEGIN
    !P.MULTI  = [0,1,2]
    xtickname = [' ',' ',' ',' ',' ',' ',' ',' ',' ']
    xtitle    = ' '
    postop    = [0.151,0.28,0.954,0.932]
    posbot1   = [0.151,0.10,0.954,0.28]
  ENDELSE
  LOADCT, 12, /SILENT
  PLOT, [0], [0], XTITLE=xtitle, YTITLE='Measured null depth per OB [%]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
    POSITION=postop, XTICKNAME=xtickname, XRANGE=xrange, YRANGE=yrange
  ; Overplot the estimated null values for calibrators measurements
  IF n_caldat_f GT 0 THEN BEGIN
    USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.7, THICK=1.5
    OPLOT, caltime_f*fac, 100.*tf_wb_f[idx_cal_f], PSYM=8, COLOR=100
    ERRPLOT, caltime_f*fac, 100.*(tf_wb_f[idx_cal_f]-tf_wb_etot_f[idx_cal_f]), 100.*(tf_wb_f[idx_cal_f]+tf_wb_etot_f[idx_cal_f]), COLOR=100
    cal_null = null_f[idx_cal_f] & cal_nerr = null_err_f[idx_cal_f]
    USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.7, THICK=1.5, /FILL
    OPLOT, caltime_f*fac, 100.*cal_null, PSYM=8, COLOR=100
    ERRPLOT, caltime_f*fac, 100.*(cal_null-cal_nerr), 100.*(cal_null+cal_nerr), COLOR=100
  ENDIF
  ; Overplot the estimated null values for science measurements
  IF n_scidat_f GE 1 THEN BEGIN
    sci_null = null_f[idx_sci_f] & sci_nerr = null_err_f[idx_sci_f]
    USERSYM, [1,0,-1,0,1], [0,1,0,-1,0], THICK=1.5, /FILL
    OPLOT, scitime_f*fac, 100.*sci_null, PSYM=8, COLOR=200
    ERRPLOT, scitime_f*fac, 100.*(sci_null-sci_nerr), 100.*(sci_null+sci_nerr), COLOR=200
  ENDIF

  ; Overplot the estimated T2 values at science measurements (and confidence intervals)
  LOADCT, 0, /SILENT
  IF drs.cal_method EQ 1 AND n_scidat NE 0 THEN BEGIN
    FOR i_tf = 0, n_tf-1 DO BEGIN
      OPLOT, REFORM(tf_time[i_tf,*])*fac, 100.*REFORM(tf_plot[i_tf,*]), LINESTYLE=0, THICK=4
      OPLOT, REFORM(tf_time[i_tf,*])*fac, 100.*(REFORM(tf_plot[i_tf,*])-tf_plot_err[i_tf]), LINESTYLE=1, THICK=3.5
      OPLOT, REFORM(tf_time[i_tf,*])*fac, 100.*(REFORM(tf_plot[i_tf,*])+tf_plot_err[i_tf]), LINESTYLE=1, THICK=3.5
    ENDFOR
  ENDIF

  ; Plot background bias in the bottom inset if requested
  no_inset = 1 ; Depreciated
  IF NOT KEYWORD_SET(NO_INSET) THEN BEGIN
    yrange = [-0.004,0.004]*1D2 + avg_bias*1D2
    PLOT, [0], [0], XTITLE='UT hour', YTITLE='Bckg null [%]', TITLE=' ', XSTYLE=1, YSTYLE=1, $
      POSITION=posbot1, XRANGE=xrange, YRANGE=yrange, YMINOR=3
    LOADCT, 12, /SILENT
    ; Overplot the estimated null values for calibrators measurements
    IF n_caldat GT 0 THEN BEGIN
      USERSYM, [1,1,-1,-1,1]*0.5, [1,-1,-1,1,1]*0.5, THICK=1.5, /FILL
      OPLOT, caltime_f*fac, 1D2*bias_f[idx_cal_f], PSYM=8, COLOR=100
      ERRPLOT, caltime_f*fac, 1D2*(bias_f[idx_cal_f]-bias_err_f[idx_cal_f]), 1D2*(bias_f[idx_cal_f]+bias_err_f[idx_cal_f]), COLOR=100
    ENDIF
    IF n_scidat GT 0 THEN BEGIN
      USERSYM, [1,0,-1,0,1]*0.5, [0,1,0,-1,0]*0.5, THICK=1.5, /FILL
      OPLOT, scitime_f*fac, 1D2*bias[idx_sci_f], PSYM=8, COLOR=200
      ERRPLOT, scitime_f*fac, 1D2*(bias[idx_sci_f]-bias_err_f[idx_sci_f]), 1D2*(bias_f[idx_sci_f]+bias_err_f[idx_sci_f]), COLOR=200
    ENDIF
    AVGSDV, bias, avg_bias_we, rms_bias, rms_bias_mean, KAPPA=5, WEIGHT=1./bias_err^2
    OPLOT, [0,24], [avg_bias_we,avg_bias_we]*1D+2, LINESTYLE=0, THICK=4
    OPLOT, [0,24], [avg_bias_we+rms_bias_mean,avg_bias_we+rms_bias_mean]*1D2, LINESTYLE=1, THICK=4
    OPLOT, [0,24], [avg_bias_we-rms_bias_mean,avg_bias_we-rms_bias_mean]*1D2, LINESTYLE=1, THICK=4
  ENDIF

  ; Finish plot
  DEVICE, /CLOSE & END_PS
  !P.MULTI=[0,1,1]
  
  ; 3. Plot data and transfer function (POINTING based)
  ; ---------------------------------------------------
  
  plotname =  tf_path + date_lng + '_TF-PTG_APER' + STRING(aper_rad, FORMAT='(I0)') + bias_tag + '.eps'
  PREP_PS, /BOLD & LOADCT, 12, /SILENT
  DEVICE, FILEN=plotname, /ENCAPS, /COLOR,  XSIZE=17.8, YSIZE=13.2, /TIMES
  xrange     = [min_time,max_time]*fac
  null_range = MAX(null_avg_pt)-MIN(null_avg_pt)
  yrange     = [MIN(null_avg_pt)-0.75*null_range,MAX(null_avg_pt)+0.75*null_range]*1D2
  yrange[0]  = yrange[0] < null_plot[0]
  yrange[1]  = yrange[1] > null_plot[1]
  PLOT, [0], [0], XTITLE='UT hour', YTITLE='Measured null depth per pointing [%]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
        XRANGE=xrange, YRANGE=yrange
  ; Overplot the estimated null values for calibrators measurements
  FOR i_nod = 0, n_nod-1 DO BEGIN
    IF n_caldat GE 1 THEN BEGIN
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, caltime_pt[*,i_nod]*fac, null_avg_pt[idx_cal_pt,i_nod]*1D2, PSYM=8, COLOR=100
      ERRPLOT, caltime_pt[*,i_nod]*fac, (null_avg_pt[idx_cal_pt,i_nod]+null_err_pt[idx_cal_pt,i_nod])*1D2, (null_avg_pt[idx_cal_pt,i_nod]-null_err_pt[idx_cal_pt,i_nod])*1D2, COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      USERSYM, [1,0,-1,0,1], [0,1,0,-1,0], THICK=1.5, /FILL
      OPLOT, scitime_pt[*,i_nod]*fac, null_avg_pt[idx_sci_pt,i_nod]*1D2, PSYM=8, COLOR=200
      ERRPLOT, scitime_pt[*,i_nod]*fac, (null_avg_pt[idx_sci_pt,i_nod]+null_err_pt[idx_sci_pt,i_nod])*1D2, (null_avg_pt[idx_sci_pt,i_nod]-null_err_pt[idx_sci_pt,i_nod])*1D2, COLOR=200
    ENDIF
  ENDFOR
  ; Overplot the estimated T2 values at science measurements (and confidence intervals)
  LOADCT, 0, /SILENT
  FOR i_nod = 0, n_nod-1 DO BEGIN
    IF drs.cal_method EQ 1 THEN BEGIN
      FOR i_tf = 0, n_tf-1 DO BEGIN
        OPLOT, REFORM(tf_time[i_tf,*])*fac, 100.*REFORM(tf_plot_pt[i_tf,*,i_nod]), LINESTYLE=0, THICK=4
        OPLOT, REFORM(tf_time[i_tf,*])*fac, 100.*(REFORM(tf_plot_pt[i_tf,*,i_nod])-tf_plot_pt_err[i_tf,i_nod]), LINESTYLE=1, THICK=3.5
        OPLOT, REFORM(tf_time[i_tf,*])*fac, 100.*(REFORM(tf_plot_pt[i_tf,*,i_nod])+tf_plot_pt_err[i_tf,i_nod]), LINESTYLE=1, THICK=3.5
      ENDFOR
    ENDIF
  ENDFOR
  DEVICE, /CLOSE & END_PS
  
  ; Now only if not RUNBIAS
  IF NOT KEYWORD_SET(RUNBIAS) THEN BEGIN
    ; 4. Plot the photometries and throughput (only computed if both SCI and CAL data)
    ; --------------------------------------------------------------------------------
    
    IF n_scidat GT 0 AND n_caldat GT 1 THEN BEGIN
      plotname =  tf_path + date_lng + '_PHOT_APER' + STRING(aper_rad, FORMAT='(I0)') + '.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=17.8, YSIZE=14.7, /TIMES
      xrange = [min_time,max_time]*fac
      yrange = [0.,1.5*MAX(cal_adu_sx_avg>cal_adu_dx_avg>sci_adu_sx_avg>sci_adu_dx_avg)]
      PLOT, [0], [0], XTITLE='UT time', YTITLE='Measured photometry [ADU]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
            XRANGE=xrange, YRANGE=yrange
      ; Overplot the estimated null values for calibrators measurements
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, caltime*fac, cal_adu_dx_avg, PSYM=8, COLOR=100
      ERRPLOT, caltime*fac, cal_adu_dx_avg+cal_adu_dx_rms, cal_adu_dx_avg-cal_adu_dx_rms, COLOR=100
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5
      OPLOT, caltime*fac, cal_adu_sx_avg, PSYM=8, COLOR=100
      ERRPLOT, caltime*fac, cal_adu_sx_avg+cal_adu_sx_rms, cal_adu_sx_avg-cal_adu_sx_rms, COLOR=100
      USERSYM, [1,0,-1,0,1], [0,1,0,-1,0], THICK=1.5, /FILL
      OPLOT, scitime*fac, [sci_adu_dx_avg], PSYM=8, COLOR=200
      ERRPLOT, scitime*fac, [sci_adu_dx_avg+sci_adu_dx_rms], [sci_adu_dx_avg-sci_adu_dx_rms], COLOR=200
      USERSYM, [1,0,-1,0,1], [0,1,0,-1,0], THICK=1.5
      OPLOT, scitime*fac, [sci_adu_sx_avg], PSYM=8, COLOR=200
      ERRPLOT, scitime*fac, [sci_adu_sx_avg+sci_adu_sx_rms], [sci_adu_sx_avg-sci_adu_sx_rms], COLOR=200
      DEVICE, /CLOSE & END_PS
         
      plotname =  tf_path + date_lng + '_TRANS_APER' + STRING(aper_rad, FORMAT='(I0)') + '.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=17.8, YSIZE=14.7, /TIMES
      xrange = [min_time,max_time]*fac
      yrange = [0.,1.5*MAX(jy2adu_avg)]
      PLOT, [0], [0], XTITLE='UT time', YTITLE='Jy to ADU', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
            XRANGE=xrange, YRANGE=yrange, XTHICK=xthick, YTHICK=ythick, CHARTHICK=charthick, CHARSIZE=charsize
      ; Overplot the estimated null values for calibrators measurements
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, caltime*fac, jy2adu_avg, PSYM=8, COLOR=100
      ERRPLOT, caltime*fac, jy2adu_avg+jy2adu_rms, jy2adu_avg-jy2adu_rms, COLOR=100
      DEVICE, /CLOSE & END_PS
    ENDIF
    
    ; 5. Plot NSC results
    ; -------------------
    
    IF FXPAR(header, 'NUL_MODE') EQ 2 THEN BEGIN
      ; a. Plot astrophysical nulls vs chi2
      chi2       = l1data.nsc_chi2
      chi2_range = MAX(chi2)-MIN(chi2)
      plotname =  tf_path + date_lng + '_NAS-VS-CHI2_APER' + STRING(aper_rad, FORMAT='(I0)') + '.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=17.8, YSIZE=14.7, /TIMES
      xrange     = [(MIN(chi2)-0.05*chi2_range)>0.1,MAX(chi2)+0.05*chi2_range]
      null_range = MAX(null)-MIN(null)
      yrange     = [MIN(null)-0.75*null_range,MAX(null)+0.75*null_range]*1D2
      PLOT, [0], [0], XTITLE='!9c!3!e2!ir!3', YTITLE='Best-fit null depth [%]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
            XRANGE=xrange, YRANGE=yrange,  /XLOG
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        cal_chi2 = l1data[idx_cal].nsc_chi2
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, cal_chi2, 1D2*cal_null, PSYM=8, COLOR=100
        ERRPLOT, cal_chi2, 1D2*(cal_null-cal_nerr), 1D2*(cal_null+cal_nerr), COLOR=100  
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        sci_chi2 = l1data[idx_cal].nsc_chi2
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, cal_chi2, 1D2*sci_null, PSYM=8, COLOR=200
        ERRPLOT, cal_chi2, 1D2*(sci_null-sci_nerr), 1D2*(sci_null+sci_nerr), COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
      
      ; b. Plot astrophysical nulls vs mean phase
      phavg     = l1data.nsc_phavg/!DPI*wav_uniq[0]*1D+6*0.5
      ph_range  = MAX(phavg)-MIN(phavg)
      plotname  =  tf_path + date_lng + '_NAS-VS-PHAVG_APER' + STRING(aper_rad, FORMAT='(I0)') + '.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=17.8, YSIZE=14.7, /TIMES
      xrange = [MIN(phavg)-0.05*ph_range,MAX(phavg)+0.05*ph_range]
      yrange = [MIN(null)-0.75*null_range,MAX(null)+0.75*null_range]*1D2
      PLOT, [0], [0], XTITLE='Best-fit mean OPD offset [um]', YTITLE='Best-fit null depth [%]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
            XRANGE=xrange, YRANGE=yrange
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        cal_phavg = [l1data[idx_cal].nsc_phavg/!DPI*wav_uniq[0]*1D+6*0.5]
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, cal_phavg, 1D2*cal_null, PSYM=8, COLOR=100
        ERRPLOT, cal_phavg, 1D2*(cal_null-cal_nerr), 1D2*(cal_null+cal_nerr), COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        sci_phavg = [l1data[idx_sci].nsc_phavg/!DPI*wav_uniq[0]*1D+6*0.5]
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, sci_phavg, 1D2*sci_null, PSYM=8, COLOR=200
        ERRPLOT, sci_phavg, 1D2*(sci_null-sci_nerr), 1D2*(sci_null+sci_nerr), COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
      
      ; c. Plot astrophysical nulls vs phase RMS
      phrms     = l1data.nsc_phrms/!DPI*wav_uniq[0]*1D+6*0.5
      ph2_range = MAX(phrms)-MIN(phrms)
      plotname  =  tf_path + date_lng + '_NAS-VS-PHRMS_APER' + STRING(aper_rad, FORMAT='(I0)') + '.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=17.8, YSIZE=14.7, /TIMES
      xrange = [MIN(phrms)-0.05*ph_range,MAX(phrms)+0.05*ph_range]
      yrange = [MIN(null)-0.75*null_range,MAX(null)+0.75*null_range]*1D2
      PLOT, [0], [0], XTITLE='Best-fit OPD jitter [um]', YTITLE='Best-fit null depth [%]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
        XRANGE=xrange, YRANGE=yrange
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        cal_phrms = [l1data[idx_cal].nsc_phrms/!DPI*wav_uniq[0]*1D+6*0.5]
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, cal_phrms, 1D2*cal_null, PSYM=8, COLOR=100
        ERRPLOT, cal_phrms, 1D2*(cal_null-cal_nerr), 1D2*(cal_null+cal_nerr), COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        sci_phrms = [l1data[idx_sci].nsc_phrms/!DPI*wav_uniq[0]*1D+6*0.5]
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, sci_phrms, 1D2*sci_null, PSYM=8, COLOR=200
        ERRPLOT, sci_phrms, 1D2*(sci_null-sci_nerr), 1D2*(sci_null+sci_nerr), COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
      
      ; d. Plot mean phase vs time
      ph_range  = MAX(phavg)-MIN(phavg)
      plotname  =  tf_path + date_lng + '_PHAVG_APER' + STRING(aper_rad, FORMAT='(I0)') + '.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
      xrange = [min_time,max_time]*fac
      yrange = [0.,MAX(phavg)+0.05*ph_range]
      PLOT, [0], [0], XTITLE='UT time', YTITLE='Best-fit mean OPD offset [um]', TITLE=title_day, XSTYLE=1, YSTYLE=9, $
            XRANGE=xrange, YRANGE=yrange, POSITION=[0.12,0.12,0.88,0.90]
      ; Add right vertical axis
      AXIS, YTITLE ='Best-fit mean phase offset [rad]', /YAXIS, YSTYLE=1, YRANGE=yrange*!DPI/(wav_uniq[0]*1D+6*0.5)
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        cal_phavg_err = [l1data[idx_cal].nsc_phavg_err/!DPI*wav_uniq[0]*1D+6*0.5]
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, caltime*fac, cal_phavg, PSYM=8, COLOR=100
        ERRPLOT, caltime*fac, cal_phavg-cal_phavg_err, cal_phavg+cal_phavg_err, COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
         sci_phavg_err = [l1data[idx_sci].nsc_phavg_err/!DPI*wav_uniq[0]*1D+6*0.5]
         USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
         OPLOT, scitime*fac, sci_phavg, PSYM=8, COLOR=200
         ERRPLOT, scitime*fac, sci_phavg-sci_phavg_err, sci_phavg+sci_phavg_err, COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
      
      ; e. Plot rms phase vs time
      phrms_err = l1data.nsc_phavg_err/!DPI*wav_uniq[0]*1D+6*0.5
      ph_range  = MAX(phrms+phrms_err)-MIN(phrms-phrms_err)
      plotname  =  tf_path + date_lng + '_PHRMS_APER' + STRING(aper_rad, FORMAT='(I0)') + '.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
      xrange = [min_time,max_time]*fac
      yrange = [0.,MAX(phrms+phrms_err)+0.05*ph_range]
      PLOT, [0], [0], XTITLE='UT time', YTITLE='Best-fit OPD jitter [!9m!3m]', TITLE=title_day, XSTYLE=1, YSTYLE=9, $
            XRANGE=xrange, YRANGE=yrange, POSITION=[0.12,0.12,0.88,0.90]
      ; Add right vertical axis
      AXIS, YTITLE ='Best-fit phase jitter [rad]', /YAXIS, YSTYLE=1, YRANGE=yrange*!DPI/(wav_uniq[0]*1D+6*0.5)
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        cal_phrms_err = [l1data[idx_cal].nsc_phavg_err/!DPI*wav_uniq[0]*1D+6*0.5]
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, caltime*fac, cal_phrms, PSYM=8, COLOR=100
        ERRPLOT, caltime*fac, cal_phrms-cal_phrms_err, cal_phrms+cal_phrms_err, COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        sci_phrms_err = [l1data[idx_sci].nsc_phavg_err/!DPI*wav_uniq[0]*1D+6*0.5]
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, scitime*fac, sci_phrms, PSYM=8, COLOR=200
        ERRPLOT, scitime*fac, sci_phrms-sci_phrms_err, sci_phrms+sci_phrms_err, COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
      
      ; f. Plot chi2  vs time
      plotname  =  tf_path + date_lng + '_CHI2_APER' + STRING(aper_rad, FORMAT='(I0)') + '.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
      yrange = [(MIN(cal_chi2<sci_chi2)-0.05*chi2_range)>0.1,MAX(cal_chi2>sci_chi2)+0.05*chi2_range]
      PLOT, [0], [0], XTITLE='UT time', YTITLE='!9c!3!e2!ir!3', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
        XRANGE=xrange, YRANGE=yrange
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, caltime*fac, cal_chi2, PSYM=8, COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, scitime*fac, sci_chi2, PSYM=8, COLOR=200
      ENDIF
      chi2_mean = MEAN([cal_chi2,sci_chi2])
      OPLOT, [-24,24], [chi2_mean, chi2_mean], LINESTYLE=1, THICK=3
      DEVICE, /CLOSE & END_PS
      
      ; g. Plot nas error vs sigma phi
      plotname  =  tf_path + date_lng + '_NASERR-VS-PHIDPHI_APER' + STRING(aper_rad, FORMAT='(I0)') + '.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
      PLOT, [0], [0], XTITLE='(Best-fit phase jitter !9s!if!3!n)/(Best-fit mean phase !9m!if!3!n)', YTITLE='Error on best-fit null depth [%]', TITLE=title_day, XSTYLE=1, YSTYLE=1,  $
            XRANGE=[0.,5.0], YRANGE=[0.,1D2*(MAX(sci_nerr)>MAX(cal_nerr))+0.1]
      ; Overplot the estimated null values for calibrators measurements
      IF n_scidat GE 1 THEN BEGIN
       phidphi_c = (cal_phrms)/cal_phavg;cal_phrms*SQRT(cal_phavg)
       USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
       OPLOT, phidphi_c, 1D2*cal_nerr, PSYM=8, COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        phidphi_s = (sci_phrms)/sci_phavg;sci_phrms*SQRT(sci_phavg)
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, phidphi_s, 1D2*sci_nerr, PSYM=8, COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
      
      ; h. Plot measured smttau vs phase
      plotname  =  tf_path + date_lng + '_PHRMS-vs_PWV.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
      xrange = [0.,5.0]
      yrange = [0.,MAX(phrms)+0.5*ph_range]
      PLOT, [0], [0], XTITLE='PWV [mm]', YTITLE='Measured OPD jitter [um]', TITLE=title_day, XSTYLE=1, YSTYLE=9, $
        XRANGE=xrange, YRANGE=yrange, POSITION=[0.12,0.12,0.88,0.90]
      ; Add right vertical axis
      AXIS, YTITLE ='Measured phase jitter [rad]', /YAXIS, YSTYLE=1, YRANGE=yrange*!DPI/(wav_uniq[0]*1D+6*0.5)
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, smttau[idx_cal], cal_phrms, PSYM=8, COLOR=100
        ERRPLOT,  smttau[idx_cal], cal_phrms-cal_phrms_err, cal_phrms+cal_phrms_err, COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, smttau[idx_sci], sci_phrms, PSYM=8, COLOR=200
        ERRPLOT,  smttau[idx_sci], sci_phrms-sci_phrms_err, sci_phrms+sci_phrms_err, COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
      
      ; i. Plot raw null RMS vs measured error
      plotname  =  tf_path + date_lng + '_PHERR-vs_NULLRMS.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
      xrange = [0.,MAX(null_rms)]
      yrange = [0.,MAX(phrms)+0.5*ph_range]
      PLOT, [0], [0], XTITLE='(Raw null RMS)/(background RMS relative to peak)', YTITLE='Measured OPD jitter [um]', TITLE=title_day, XSTYLE=1, YSTYLE=9, $
        XRANGE=xrange, YRANGE=yrange, POSITION=[0.12,0.12,0.88,0.90]
      ; Add right vertical axis
      AXIS, YTITLE ='Measured phase jitter [rad]', /YAXIS, YSTYLE=1, YRANGE=yrange*!DPI/(wav_uniq[0]*1D+6*0.5)
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, null_rms[idx_cal], cal_phrms, PSYM=8, COLOR=100
        ERRPLOT,  null_rms[idx_cal], cal_phrms-cal_phrms_err, cal_phrms+cal_phrms_err, COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, null_rms[idx_sci], sci_phrms, PSYM=8, COLOR=200
        ERRPLOT,  null_rms[idx_sci], sci_phrms-sci_phrms_err, sci_phrms+sci_phrms_err, COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
      
      ; j. Plot raw null RMS vs measured error
      plotname  =  tf_path + date_lng + '_NULLAVG-vs_NULLRMS.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
      xrange = [0.,MAX(null_rms)]
      yrange = 2*MAX(sci_null)*[-1.,1.]
      PLOT, [0], [0], XTITLE='(Raw null RMS)/(background RMS relative to peak)', YTITLE='Measured null', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
                XRANGE=xrange, YRANGE=yrange, POSITION=[0.12,0.12,0.88,0.90]
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, null_rms[idx_cal], cal_null, PSYM=8, COLOR=100
        ERRPLOT,  null_rms[idx_cal], cal_null-cal_nerr, cal_null+cal_nerr, COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, null_rms[idx_sci], sci_null, PSYM=8, COLOR=200
        ERRPLOT,  null_rms[idx_sci], sci_null-sci_nerr, sci_null+sci_nerr, COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
      
      ; k. Plot measured error vs error expected from background noise
      plotname  =  tf_path + date_lng + '_NULLERR-vs_BCKGERR.eps'
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
      xrange = [0.,MAX(null_err_phot)]*1D2
      yrange = [0, 2*MAX(sci_nerr)]*1D2
      PLOT, [0], [0], XTITLE='Photometric error on background relative to peak [%]', YTITLE='Null errror [%]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
        XRANGE=xrange, YRANGE=yrange, POSITION=[0.12,0.12,0.88,0.90]
      OPLOT, INDGEN(10), INDGEN(10)
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, null_err_phot[idx_cal]*1D2, cal_nerr*1D2, PSYM=8, COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, null_err_phot[idx_sci]*1D2, sci_nerr*1D2, PSYM=8, COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
      
      ; k. Plot measured error vs error expected from background noise
      plotname  =  tf_path + date_lng + '_EXCESS-PHOT-ERR-vs-time.eps'
      excess_cal= cal_nerr/null_err_phot[idx_cal]
      excess_sci= sci_nerr/null_err_phot[idx_sci]
      fac        = 24D
      min_time   = 0.98*MIN(time)
      max_time   = 1.02*MAX(time)
      xrange     = [min_time,max_time]*fac
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
      yrange = [0, 1.5*MAX(excess_sci)]
      PLOT, [0], [0], XTITLE='UT time', YTITLE='(Null error) / (Photometric error relative to peak)', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
        XRANGE=xrange, YRANGE=yrange, POSITION=[0.12,0.12,0.88,0.90]
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, caltime*fac, excess_cal, PSYM=8, COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, scitime*fac, excess_sci, PSYM=8, COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
      
      ; k. Plot measured error vs error expected from background noise
      plotname  =  tf_path + date_lng + '_EXCESS_PHASE-ERR-vs-time.eps'
      excess_cal= cal_nerr/null_err_phase[idx_cal]
      excess_sci= sci_nerr/null_err_phase[idx_sci]
      fac        = 24D
      min_time   = 0.98*MIN(time)
      max_time   = 1.02*MAX(time)
      xrange     = [min_time,max_time]*fac
      PREP_PS, /BOLD & LOADCT, 12, /SILENT
      DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
      yrange = [0, 1.5*MAX(excess_sci)]
      PLOT, [0], [0], XTITLE='UT time', YTITLE='(Null error) / (Error from phase variations)', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
        XRANGE=xrange, YRANGE=yrange, POSITION=[0.12,0.12,0.88,0.90]
      ; Overplot the estimated null values for calibrators measurements
      IF n_caldat GE 1 THEN BEGIN
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, caltime*fac, excess_cal, PSYM=8, COLOR=100
      ENDIF
      IF n_scidat GE 1 THEN BEGIN
        USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
        OPLOT, scitime*fac, excess_sci, PSYM=8, COLOR=200
      ENDIF
      DEVICE, /CLOSE & END_PS
    ENDIF
    
    ; 6. Plot fpc piston vs time
    ; --------------------------
    
    ph_range  = MAX(fpcp)-MIN(fpcp)
    plotname  =  tf_path + date_lng + '_PCP.eps'
    PREP_PS, /BOLD & LOADCT, 12, /SILENT
    DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
    xrange = [min_time,max_time]*fac
    yrange = [0.,MAX(fpcp)+0.25*ph_range]
    PLOT, [0], [0], XTITLE='UT time', YTITLE='PHASECam OPD jitter [um]', TITLE=title_day, XSTYLE=1, YSTYLE=9, $
      XRANGE=xrange, YRANGE=yrange, POSITION=[0.12,0.12,0.88,0.90]
    ; Add right vertical axis
    AXIS, YTITLE ='PHASECam phase jitter [rad]', /YAXIS, YSTYLE=1, YRANGE=yrange*!DPI/(wav_uniq[0]*1D+6*0.5)
    ; Overplot the estimated null values for calibrators measurements
    IF n_caldat GE 1 THEN BEGIN
      cal_fpcp = [phasec[idx_cal].fpc_pists*plc_wav[idx_cal]/360.]
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, caltime*fac, cal_fpcp, PSYM=8, COLOR=100
      ;ERRPLOT, caltime*fac, cal_fpcp-cal_fpcp_err, cal_fpcp+cal_fpcp_err, COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      sci_fpcp = [phasec[idx_sci].fpc_pists*plc_wav[idx_sci]/360.]
      USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
      OPLOT, scitime*fac, sci_fpcp, PSYM=8, COLOR=200
      ;ERRPLOT, scitime*fac, sci_fpcp-sci_fpcp_err, sci_fpcp+sci_fpcp_err, COLOR=200
    ENDIF
    DEVICE, /CLOSE & END_PS
  
    ; 7. Plot measured RMS phase (over the last DITms)  vs time
    ; --------------------------------------------------------
    
    ph_range  = MAX(phstd)-MIN(phstd)
    plotname  = tf_path + date_lng + '_PHSTD.eps'
    PREP_PS, /BOLD & LOADCT, 12, /SILENT
    DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
    yrange = [0.,MAX(phstd)+2.0*ph_range]
    PLOT, [0], [0], XTITLE='UT time', YTITLE='Measured OPD jitter [um]', TITLE=title_day, XSTYLE=1, YSTYLE=9, $
      XRANGE=xrange, YRANGE=yrange, POSITION=[0.12,0.12,0.88,0.90]
    ; Add right vertical axis
    AXIS, YTITLE ='Measured phase jitter [rad]', /YAXIS, YSTYLE=1, YRANGE=yrange*!DPI/(wav_uniq[0]*1D+6*0.5)
    ; Overplot the estimated null values for calibrators measurements
    IF n_caldat GE 1 THEN BEGIN
      cal_phstd     = [phasec[idx_cal].pcphstd*plc_wav[idx_cal]/360.]          ; Convert to um
      cal_phstd_err = [phasec[idx_cal].pcphstd_err*plc_wav[idx_cal]/360.]           ; Convert to um
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, caltime*fac, cal_phstd, PSYM=8, COLOR=100
      ERRPLOT, caltime*fac, cal_phstd-cal_phstd_err, cal_phstd+cal_phstd_err, COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      sci_phstd     = [phasec[idx_sci].pcphstd*plc_wav[idx_sci]/360.]           ; Convert to um
      sci_phstd_err = [phasec[idx_sci].pcphstd_err*plc_wav[idx_sci]/360.]           ; Convert to um
      USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
      OPLOT, scitime*fac, sci_phstd, PSYM=8, COLOR=200
      ERRPLOT, scitime*fac, sci_phstd-sci_phstd_err, sci_phstd+sci_phstd_err, COLOR=200
    ENDIF
    DEVICE, /CLOSE & END_PS
    
    phavg     = phasec.pcphmean*plc_wav/360.
    ph_range  = (MAX(phavg)-MIN(phavg))
    plotname  = tf_path + date_lng + '_PHAVG.eps'
    PREP_PS, /BOLD & LOADCT, 12, /SILENT
    DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
    yrange = [MIN(phavg)-0.5*ph_range,MAX(phavg)+0.5*ph_range]
    PLOT, [0], [0], XTITLE='UT time', YTITLE='Measured mean OPD [um]', TITLE=title_day, XSTYLE=1, YSTYLE=9, $
      XRANGE=xrange, YRANGE=yrange, POSITION=[0.12,0.12,0.88,0.90]
    ; Add right vertical axis
    AXIS, YTITLE ='Measured mean phase [rad]', /YAXIS, YSTYLE=1, YRANGE=yrange*!DPI/(wav_uniq[0]*1D+6*0.5)
    ; Overplot the estimated null values for calibrators measurements
    IF n_caldat GE 1 THEN BEGIN
      cal_phavg     = [phasec[idx_cal].pcphmean*plc_wav[idx_cal]*1D6/360.]           ; Convert to um
      cal_phavg_err = [phasec[idx_cal].pcphmean_err*plc_wav[idx_cal]*1D6/360.]        ; Convert to um
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, caltime*fac, cal_phavg, PSYM=8, COLOR=100
      ERRPLOT, caltime*fac, cal_phavg-cal_phavg_err, cal_phavg+cal_phavg_err, COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      sci_phavg     = [phasec[idx_sci].pcphmean*plc_wav[idx_sci]*1D6/360.]           ; Convert to um
      sci_phavg_err = [phasec[idx_sci].pcphmean_err*plc_wav[idx_sci]*1D6/360.]      ; Convert to um
      USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
      OPLOT, scitime*fac, sci_phavg, PSYM=8, COLOR=200
      ERRPLOT, scitime*fac, sci_phavg-sci_phavg_err, sci_phavg+sci_phavg_err, COLOR=200
    ENDIF
    DEVICE, /CLOSE & END_PS
    
    ; 8. Plot measured seeing vs time
    ; -------------------------------
    
    plotname  =  tf_path + date_lng + '_SEEING.eps'
    PREP_PS, /BOLD & LOADCT, 12, /SILENT
    DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
    yrange = [0.,2.0]
    PLOT, [0], [0], XTITLE='UT time', YTITLE='Seeing [arcsec]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
          XRANGE=xrange, YRANGE=yrange
    ; Overplot the estimated null values for calibrators measurements
    IF n_scidat GE 1 THEN BEGIN
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, caltime*fac, seeing[idx_cal], PSYM=8, COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
      OPLOT, scitime*fac, seeing[idx_sci], PSYM=8, COLOR=200
    ENDIF
    DEVICE, /CLOSE & END_PS
    
    ; 9. Plot measured smttau vs time
    ; -------------------------------
    
    plotname  =  tf_path + date_lng + '_PWV.eps'
    PREP_PS, /BOLD & LOADCT, 12, /SILENT
    DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
    yrange = [0.,10.0]
    PLOT, [0], [0], XTITLE='UT time', YTITLE='PWV [mm]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
      XRANGE=xrange, YRANGE=yrange
    ; Overplot the estimated null values for calibrators measurements
    IF n_scidat GE 1 THEN BEGIN
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, caltime*fac, smttau[idx_cal], PSYM=8, COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
      OPLOT, scitime*fac, smttau[idx_sci], PSYM=8, COLOR=200
    ENDIF
    DEVICE, /CLOSE & END_PS
    
    ; 9. Plot measured windspeed vs time
    ; ----------------------------------
    
    plotname  =  tf_path + date_lng + '_WIND.eps'
    PREP_PS, /BOLD & LOADCT, 12, /SILENT
    DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
    yrange = [0.,25.0]
    PLOT, [0], [0], XTITLE='UT time', YTITLE='Wind speed [m/s]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
      XRANGE=xrange, YRANGE=yrange;, POSITION=[0.12,0.12,0.88,0.90]
    ; Overplot the estimated null values for calibrators measurements
    IF n_caldat GE 1 THEN BEGIN
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, caltime*fac, wind[idx_cal], PSYM=8, COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
      OPLOT, scitime*fac, wind[idx_sci], PSYM=8, COLOR=200
    ENDIF
    DEVICE, /CLOSE & END_PS
    
    ; 10. Plot position on detector
    ; -----------------------------
    plotname  =  tf_path + date_lng + '_BEAM-POS.eps'
    PREP_PS, /BOLD & LOADCT, 12, /SILENT
    DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
    xrange = [MIN(xcen)-10,MAX(xcen)+10]
    yrange = [MIN(ycen)-10,MAX(ycen)+10]
    PLOT, [0], [0], XTITLE='Initial x position on detector [pix]', YTITLE='Initial y position on detector [pix]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
      XRANGE=xrange, YRANGE=yrange
    ; Overplot the estimated null values for calibrators measurements
    IF n_caldat GE 1 THEN BEGIN
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, xcen[idx_cal], ycen[idx_cal], PSYM=8, COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
      OPLOT, xcen[idx_sci], ycen[idx_sci], PSYM=8, COLOR=200
    ENDIF
    DEVICE, /CLOSE & END_PS
    
    ; 11. Plot position on detector (vs time)
    ; -----------------------------
    fac       = 24D
    min_time  = 0.98*MIN(time)
    max_time  = 1.02*MAX(time)
    xrange    = [min_time,max_time]*fac
    yrange    = [MIN(xcen)-10,MAX(xcen)+10]
    plotname  =  tf_path + date_lng + '_BEAM-XPOS-TIME.eps'
    PREP_PS, /BOLD & LOADCT, 12, /SILENT
    DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
    PLOT, [0], [0], XTITLE='UT hour', YTITLE='Initial x position on detector [pix]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
      XRANGE=xrange, YRANGE=yrange
    ; Overplot the estimated null values for calibrators measurements
    IF n_caldat GE 1 THEN BEGIN
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, caltime*fac, xcen[idx_cal], PSYM=8, COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
      OPLOT, scitime*fac, xcen[idx_sci], PSYM=8, COLOR=200
    ENDIF
    DEVICE, /CLOSE & END_PS
    
    yrange    = [MIN(ycen)-10,MAX(ycen)+10]
    plotname  =  tf_path + date_lng + '_BEAM-YPOS-TIME.eps'
    PREP_PS, /BOLD & LOADCT, 12, /SILENT
    DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
    PLOT, [0], [0], XTITLE='UT hour', YTITLE='Initial x position on detector [pix]', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
      XRANGE=xrange, YRANGE=yrange
    ; Overplot the estimated null values for calibrators measurements
    IF n_caldat GE 1 THEN BEGIN
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, caltime*fac, ycen[idx_cal], PSYM=8, COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
      OPLOT, scitime*fac, ycen[idx_sci], PSYM=8, COLOR=200
    ENDIF
    DEVICE, /CLOSE & END_PS
    
    ; 12. Plot position on detector (vs null)
    ; -----------------------------
  
    xrange     = [MIN(xcen)-10,MAX(xcen)+10]
    null_range = MAX(null)-MIN(null)
    yrange     = [MIN(null)-0.75*null_range,MAX(null)+0.75*null_range]*1D2
    plotname   =  tf_path + date_lng + '_BEAM-XPOS-NULL.eps'
    PREP_PS, /BOLD & LOADCT, 12, /SILENT
    DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
    PLOT, [0], [0], XTITLE='Initial x position on detector [pix]', YTITLE='Null', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
      XRANGE=xrange, YRANGE=yrange
    ; Overplot the estimated null values for calibrators measurements
    IF n_caldat GE 1 THEN BEGIN
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, xcen[idx_cal], 100.*cal_null, PSYM=8, COLOR=100
      ERRPLOT, xcen[idx_cal], 100.*(cal_null-cal_nerr), 100.*(cal_null+cal_nerr), COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
      OPLOT, xcen[idx_sci], 100.*sci_null, PSYM=8, COLOR=200
      ERRPLOT, xcen[idx_sci], 100.*(sci_null-sci_nerr), 100.*(sci_null+sci_nerr), COLOR=200
    ENDIF
    DEVICE, /CLOSE & END_PS
  
    xrange     = [MIN(ycen)-10,MAX(ycen)+10]
    plotname   =  tf_path + date_lng + '_BEAM-YPOS-NULL.eps'
    PREP_PS, /BOLD & LOADCT, 12, /SILENT
    DEVICE, FILEN=plotname, /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7, /TIMES
    PLOT, [0], [0], XTITLE='Initial y position on detector [pix]', YTITLE='Null', TITLE=title_day, XSTYLE=1, YSTYLE=1, $
      XRANGE=xrange, YRANGE=yrange;, POSITION=[0.1,0.1,0.90,0.92]
    ; Overplot the estimated null values for calibrators measurements
    IF n_caldat GE 1 THEN BEGIN
      USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
      OPLOT, ycen[idx_cal], 100.*cal_null, PSYM=8, COLOR=100
      ERRPLOT, ycen[idx_cal], 100.*(cal_null-cal_nerr), 100.*(cal_null+cal_nerr), COLOR=100
    ENDIF
    IF n_scidat GE 1 THEN BEGIN
      USERSYM, [1,0,-1,0,1]*1.0, [0,1,0,-1,0]*1.0, THICK=1.5, /FILL
      OPLOT, ycen[idx_sci], 100.*sci_null, PSYM=8, COLOR=200
      ERRPLOT, ycen[idx_sci], 100.*(sci_null-sci_nerr), 100.*(sci_null+sci_nerr), COLOR=200
    ENDIF
    DEVICE, /CLOSE & END_PS
  ENDIF


; Close file
IF KEYWORD_SET(log_file) THEN BEGIN
  PRINTF,lun, ' '
  PRINTF,lun, ' '
  CLOSE, lun
  FREE_LUN, lun
ENDIF

bck_irad = 0
bck_orad =0
; Write summary file
IF KEYWORD_SET(log_file) AND NOT KEYWORD_SET(RUNBIAS) THEN BEGIN
  ; Create log file
  log_file =  l2_dir + 'reduction_history.txt'
  IF NOT FILE_TEST(log_file) THEN BEGIN
    OPENW, lun, log_file, /GET_LUN, WIDTH=800, /APPEND
    PRINTF, lun, 'Column signification'
    PRINTF, lun, '1: L1 file version'
    PRINTF, lun, '2: Statistical error on TF [%]'
    PRINTF, lun, '3: Systematic error on TF [%]'
    PRINTF, lun, '4: Background subtraction mode'
    PRINTF, lun, '5: Null mode (0: best x%, 1: mode, and 2: NSC)'
    PRINTF, lun, '6: Null correction based on high-frequency phase jitter'
    PRINTF, lun, '7: Multiplier on number of histogram bins (nbins_fac*sqrt(Np))'
    PRINTF, lun, '8: Bin size for NSC reduction (0: constant, 1: variable).'
    PRINTF, lun, '9: Number of bootstrap samples'
    PRINTF, lun, '10: 3-element vector containing the number of elements in the chi2 cube along the null, mean phase and rms phase directions respectively.'
    PRINTF, lun, '11: Minimum number of occurences per bin for the fit'
    PRINTF, lun, '12: Range of acceptable null (before NSC, in %)'
    PRINTF, lun, '13: Radius of the photmetric aperture (in pixels)'
    PRINTF, lun, '14: Inner radius of background region (in pixels)'
    PRINTF, lun, '15: Outer radius of background region (in pixels, 0 to use the closer vertical limit of the channel)'
    PRINTF, lun, '16: L1 DRS version ID number'
    PRINTF, lun, '17: L2 DRS version ID number'
    PRINTF, lun, '18: Execution date'
    PRINTF, lun, ' 1 ', '  ', '   2  ','  ','   3  ','  ',' 4','  ','5','  ','6','  ','7','  ','8','  ','  9 ','  ','    10    ',' ','11','  ','    12     ','  ','13','  ','14','  ','15','  ',' 16 ','  ',' 17 ','  ','        18     '
    PRINTF, lun, '--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
  ENDIF ELSE OPENW, lun, log_file, /GET_LUN, WIDTH=800, /APPEND
  PRINTF,lun, STRING(file_ver, FORMAT='(I03)'), '  ', STRING(MEAN(tf_estat)*1D2, FORMAT='(F6.4)'), '  ', STRING(MEAN(tf_esyst)*1D+2, FORMAT='(F6.4)'), '  ', STRING(bck_mode, FORMAT='(I0)'), '  ', STRING(null_mod, FORMAT='(I0)'), '  ', STRING(null_cor, FORMAT='(I0)'), '  ', $
              STRING(nbin_fac, FORMAT='(I0)'), '  ', STRING(nsc_bins, FORMAT='(I0)'), '  ', STRING(n_btstrp, FORMAT='(I04)'), '  ', nsc_cube, '  ', STRING(nsc_omin, FORMAT='(I02)'), '  ', null_lim, '  ', STRING(aper_rad, FORMAT='(I02)'), '  ', STRING(bck_irad, FORMAT='(I02)'), '  ', STRING(bck_orad, FORMAT='(I02)'), '  ', $
              STRING(FXPAR(header, 'DRS_VERS'), FORMAT='(F4.2)'), '  ', STRING(drs_version, FORMAT='(F4.2)'), '  ', SYSTIME()
  CLOSE, lun
  FREE_LUN, lun
ENDIF
END
