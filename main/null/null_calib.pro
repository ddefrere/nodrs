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
;   Version 6.1,  24-MAY-2024, DD: Update file permission
;   Version 6.2,  27-JAN-2025, DD: Corrected syntax of bckg_bias error term (used only for bckg_mode greater than 1)

pro NULL_CALIB, date, cfg_file, calpob = calpob, log_file = log_file, no_inset = no_inset, remove_id = remove_id, remove_ob = remove_ob, runbias = runbias, info = info, plot = plot, version = version
  compile_opt idl2

  ; Calibration routine version, to be updated manually
  drs_version = 5.7
  drs_date = '10-AUG-2020'

  ; Start actual code
  on_error, 0

  ; Keyword sanity check
  if not keyword_set(info) then info = 1
  if not keyword_set(no_inset) then no_inset = 0

  ; Define running and plotting parameters
  n_lam = 20 ; Number of wavelength bins within the bandwidth (used for broadband computation)
  n_time = 1000 ; Number of points in the interpolated TF curve
  charsize = 1.2 ; Character size in plots
  charthick = 3.0 ; Character thickness in plots

  ; DEFINE GLOBAL VARIABLES
  ; ***********************

  ; Astronomical and physical constants
  GET_PRM, prm

  ; Obtain the definition of the configuration
  GET_CNF, cnf, instrum = 'nomic'

  ; Read config file with reduction parameters
  GET_DRS, drs, 'nodrs/cfg/' + cfg_file
  if keyword_set(runbias) then begin
    if keyword_set(calpob) then begin
      id_name = 'OB-BIAS'
      drs.cal_mode = 1
    endif else id_name = 'PT-BIAS'
  endif else begin
    if keyword_set(calpob) then begin
      id_name = 'OB'
      drs.cal_mode = 1
    endif else id_name = 'PT'
  endelse

  ; Recover the IDL running path
  DECLARE_PATH, pth, instrum = 'nomic'

  ; COMPUTE FILE PATH AND READ DATA (+ basic heaser information)
  ; *******************************

  ; Derive long version of the date
  date_lng = '20' + strmid(date, 0, 2) + '-' + strmid(date, 2, 2) + '-' + strmid(date, 4, 2)
  if max(drs.null_rad) ne 0 and not strmatch(drs.dir_label, '*_APR') then drs.dir_label = drs.dir_label + '_APR' ; If aperture radius is set, redefine dir_label
  l1fits_path = pth.l1Fits_path + date_lng + drs.dir_label + pth.sep
  if not file_test(l1fits_path) then begin
    message, 'No L1 data for ' + date_lng, /continue
    RETURN
  endif

  ; Output directory
  l2_dir = pth.l2Fits_path + pth.sep + date_lng + drs.dir_label + pth.sep
  if not file_test(l2_dir) then file_mkdir, l2_dir
  spawn, 'chmod 775 ' + l2_dir

  ; File version
  if keyword_set(version) then vers = '_v' + string(version, format = '(I0)') else vers = ''

  ; Init log file
  if keyword_set(log_file) then begin
    log_file = l2_dir + date_lng + vers + '.txt'
    openw, lun, log_file, /get_lun, width = 800, /append
    printf, lun, ' '
    printf, lun, '******************************************************************************************************* '
    printf, lun, ' '
    printf, lun, 'NULL_CALIB.pro version ' + string(drs_version, format = '(F3.1)') + ' -- ' + drs_date + ' -- Denis Defrère - Steward Observatory (denis@lbti.org)'
    printf, lun, 'Data calibrated on ' + systime()
    printf, lun, ' '
    printf, lun, '******************************************************************************************************* '
    printf, lun, ' '
  endif else lun = -1

  ; Search for L1 summary file and read it
  l1files = file_search(l1fits_path, 'UT' + date_lng + vers + '.fits', count = n_files)
  if not file_test(l1files[0]) then begin
    message, 'No L1 data for ' + date_lng, /continue
    RETURN
  endif else begin
    ; Search for the right aperture radius
    header = HEADFITS(l1files[0], /silent) ; Read main header
    aper_rad = FXPAR(header, 'NULL_RAD', /nocontinue)
    if aper_rad ne drs.null_rad then begin
      ; Find the latest file with the right aperture
      l1files = file_search(l1fits_path, 'UT' + date_lng + '_v*.fits', count = n_files)
      m_time = (file_info(l1files)).mtime
      l1files = l1files[reverse(sort(m_time))]
      for i_f = 0, n_files - 1 do begin
        file_id = i_f
        if FXPAR(HEADFITS(l1files[file_id], /silent), 'NULL_RAD', /nocontinue) eq drs.null_rad then i_f = n_files
      endfor
    endif else file_id = 0
  endelse
  l1data = MRDFITS(l1files[file_id], 1, hdr_col, /silent) ; Read null data
  detect = MRDFITS(l1files[file_id], 2, /silent) ; Read detector information data
  phasec = MRDFITS(l1files[file_id], 5, /silent) ; Read PHASECam data
  weather = MRDFITS(l1files[file_id], 6, /silent) ; Read weather data
  header = HEADFITS(l1files[file_id], /silent) ; Read main header

  ; Remove duplicated OBs (temporary fix)
  idx_ok = REM_DUP(l1data.ob_id)
  if n_elements(idx_ok) gt 0 then begin
    l1data = l1data[idx_ok]
    phasec = phasec[idx_ok]
    detect = detect[idx_ok]
    weather = weather[idx_ok]
  endif else message, 'No OB to calibrate'

  ; If REMOVE_OB is set, remove the corresponding OBs from the data
  n_data = n_elements(l1data.mjd_obs)
  if keyword_set(remove_ob) then begin
    MATCH, l1data.ob_id, remove_ob, idx_out, idx_tmp, count = n_out
    if n_out gt 0 then begin
      l1data[idx_out].ob_id = -1
      idx_ok = where(l1data.ob_id ne -1, n_data)
      if n_data gt 0 then begin
        l1data = l1data[idx_ok]
        phasec = phasec[idx_ok]
        detect = detect[idx_ok]
        weather = weather[idx_ok]
      endif else message, 'No OB to calibrate'
    endif
  endif

  ; If REMOVEIDB is set, remove the corresponding OBs from the data
  n_data = n_elements(l1data.mjd_obs)
  ; l1data.file_id = INDGEN(n_data)
  if keyword_set(remove_id) then begin
    MATCH, l1data.file_id, remove_id, idx_out, idx_tmp, count = n_out
    if n_out gt 0 then begin
      l1data[idx_out].file_id = -1
      idx_ok = where(l1data.file_id ne -1, n_data)
      if n_data gt 0 then begin
        l1data = l1data[idx_ok]
        phasec = phasec[idx_ok]
        detect = detect[idx_ok]
        weather = weather[idx_ok]
      endif else message, 'No OB to calibrate'
    endif
  endif

  ; Remove high chi2 if requested
  n_data = n_elements(l1data.mjd_obs)
  if drs.chi2_lim and FXPAR(header, 'NUL_MODE') eq 2 then begin
    idx_ok = where(l1data.nsc_chi2 lt drs.chi2_lim, n_data)
    if n_data gt 0 then begin
      l1data = l1data[idx_ok]
      phasec = phasec[idx_ok]
      detect = detect[idx_ok]
      weather = weather[idx_ok]
    endif
  endif

  ; Make sure data is in chronological order
  idx_srt = sort(l1data.mjd_obs)
  l1data = l1data[idx_srt]
  phasec = phasec[idx_srt]
  detect = detect[idx_srt]
  weather = weather[idx_srt]

  ; Read relevant header information
  bck_mode = FXPAR(header, 'BCK_MODE', /nocontinue)
  bfl_mode = FXPAR(header, 'BFL_MODE', /nocontinue)
  file_ver = FXPAR(header, 'FILE_VER', /nocontinue)
  n_btstrp = FXPAR(header, 'N_BTSTRP', /nocontinue)
  nbin_fac = FXPAR(header, 'NBIN_FAC', /nocontinue)
  nsc_bins = FXPAR(header, 'NSC_BINS', /nocontinue)
  nsc_cube = FXPAR(header, 'NSC_CUBE', /nocontinue)
  nsc_mode = FXPAR(header, 'NSC_MODE', /nocontinue)
  nsc_omin = FXPAR(header, 'NSC_OMIN', /nocontinue)
  null_cor = FXPAR(header, 'NULL_COR', /nocontinue)
  null_lim = FXPAR(header, 'NULL_LIM', /nocontinue)
  null_mod = FXPAR(header, 'NUL_MODE', /nocontinue)
  aper_rad = FXPAR(header, 'NULL_RAD', /nocontinue)
  fra_mode = FXPAR(header, 'FRA_MODE', /nocontinue)
  img_mode = FXPAR(header, 'IMG_MODE', /nocontinue)
  fit_mode = FXPAR(header, 'FIT_MODE', /nocontinue)
  flx_mode = FXPAR(header, 'FLX_MODE', /nocontinue)

  ; Read main HDU structure into variables
  objname = strcompress(strlowcase(l1data.objname), /remove_all)
  pid = l1data.pid
  ob_id = l1data.ob_id
  flag = l1data.flag
  time = l1data.mjd_obs
  wav = l1data.wav_eff
  bdw = l1data.bandwidth
  nfr_ob = l1data.nfr_ob
  nfr_rej = l1data.nfr_rej
  utc = l1data.lbt_utc
  lst = l1data.lbt_lst
  ra = l1data.lbt_ra
  dec = l1data.lbt_dec
  alt = l1data.lbt_alt
  az = l1data.lbt_az
  para = l1data.lbt_para

  ; Null or bias
  if not keyword_set(runbias) then begin
    null = l1data.null_meas
    null_err = l1data.null_meas_err
  endif else begin
    null = l1data.bckg_bias
    null_err = l1data.bckg_bias_rms / sqrt(nfr_ob) ; RMS is per frame and we want the error on the mean
    if max(null) eq 0 and max(null_err) eq 0 then begin
      message, 'No background bias information.', /continue
      RETURN
    endif
  endelse

  ; Extract individual error terms
  null_rms = l1data.null_meas_rms / l1data.bckg_cnod_rms
  err_phot = abs(l1data.null_meas_avg) * (l1data.phot_err / l1data.phot_avg) ; error on NULL due to photometric uncertainty (not included in NSC and correlated between all OBs of the same pointing!)
  err_bckg = l1data.bckg_cnod_rms / sqrt(l1data.nfr_rej + l1data.nfr_ob) ; error on background subtraction from complementary nod (this is only an approximation. Go to the L1 log to see the right results)
  err_trm = sqrt(err_phot ^ 2 + err_bckg ^ 2)
  err_nsc = sqrt((l1data.null_meas_err ^ 2 - err_trm ^ 2) > 0) ; error from NSC not saved but one can comptue it from the total error

  ; Compute expected error terms
  if bck_mode ge 0 then null_err_phot = sqrt(err_bckg ^ 2 + l1data.bckg_bias_rms ^ 2) else null_err_phot = sqrt(2) * err_bckg
  null_err_phase = sqrt(4 * (l1data.nsc_phavg) ^ 2 * (l1data.nsc_phrms) ^ 2 + 2 * (l1data.nsc_phrms) ^ 4) / 4

  ; Backwards compatibility
  if TAG_EXIST(l1data, 'aper_rad') then aper = l1data.aper_rad else aper = replicate(FXPAR(header, 'APER_RAD', /nocontinue), n_data)
  if TAG_EXIST(l1data, 'bck_irad') then birad = l1data.bck_irad else birad = replicate(FXPAR(header, 'BCK_IRAD', /nocontinue), n_data)
  if TAG_EXIST(l1data, 'bck_orad') then borad = l1data.bck_orad else borad = replicate(FXPAR(header, 'BCK_ORAD', /nocontinue), n_data)
  if not TAG_EXIST(l1data, 'null_offset') then struct_add_field, l1data, 'NULL_OFFSET', fltarr(n_data)

  ; Compute intensity mismatch
  int_err = 2.0 * abs(l1data.photdx_avg - l1data.photsx_avg) / (l1data.photdx_avg + l1data.photsx_avg)

  ; Read relevant detector information
  xcen = detect.xcen
  ycen = detect.ycen
  dit = detect.int_time

  ; Read relevant weather information
  seeing = weather.seeing
  smttau = weather.smttau
  wind = weather.windspd

  ; Read PHASECam information
  if TAG_EXIST(phasec, 'plc_wav') then plc_wav = phasec.plc_wav * 1d+6 else plc_wav = replicate(2.2, n_data) ; Backward compatibility
  if TAG_EXIST(phasec, 'dith_per') then dith_per = phasec.dith_per else dith_per = replicate(0, n_data) ; Backward compatibility
  fpcp = phasec.fpc_pists * plc_wav / 360. ; Convert to um
  phavg = phasec.pcphmean * plc_wav / 360. ; Convert to um
  phstd = phasec.pcphstd * plc_wav / 360. ; Convert to um
  phavg_err = phasec.pcphmean_err * plc_wav / 360. ; Convert to um
  phstd_err = phasec.pcphstd_err * plc_wav / 360. ; Convert to um

  ; Compute hour angles for each OB
  ha = dblarr(n_data)
  for i = 0, n_data - 1 do begin
    GET_TGT, objname[i], tgt, database = pth.input_path + drs.database
    ha[i] = HOUR_ANGLE(time[i], tgt.ra, tgt.dec, lat = TEN(32, 42, 05.87), lon = -TEN(109, 52, 18.87), altitude = 3267d0) * 12. / 180d0 ; in hours
  endfor

  ; Compute UV coordinates
  u_coord = cnf.base / (206265. * wav) * cos(para * !dtor)
  v_coord = cnf.base / (206265. * wav) * sin(-para * !dtor) ; u positive towards East

  ; Assign pointing ID (successive OB on the same target)
  ; Ensure backward compatibility (pointing ID was not present in the file before)
  if not TAG_EXIST(l1data, 'pt_id') then begin
    pt_id = intarr(n_data)
    for i = 1, n_data - 1 do if strmatch(objname[i], objname[i - 1], /fold_case) eq 0 then pt_id[i] = pt_id[i - 1] + 1 else pt_id[i] = pt_id[i - 1]
  endif else pt_id = l1data.pt_id
  pt_uniq = pt_id[uniq(pt_id, sort(pt_id))]
  n_pt = n_elements(pt_uniq)

  ; Derive the number of nod positions assuming that each nod is located in a different channel (in the Y direction)
  ; If nod are uncorrelated and split_nod is not set, force the number of nods to 1
  if drs.nod_cor eq 1 or drs.split_nod ne 0 then nod_pos = fix(ycen / cnf.y_chan) else nod_pos = intarr(n_data)
  nod_uniq = nod_pos[uniq(nod_pos, sort(nod_pos))]
  n_nod = n_elements(nod_uniq)

  ; Convert time to hour
  DAYCNV, mean(time) + 2400000.5d0, yr, mn, day, hr
  title_day = string(yr, format = '(I0)') + '/' + string(mn, format = '(I0)') + '/' + string(day, format = '(I0)')
  JDCNV, yr, mn, day, 0., t0day
  t0day = t0day - 2400000.5d0
  time = time - t0day

  ; Define arrays
  tf_wb = null ; Wideband transfer function
  tf_wb_estat = null_err ; Statistical error on wideband TF
  tf_wb_etot = null_err ; Total error on wideband TF
  tf_wb_esyst = fltarr(n_data) ; Systematic error on wideband TF
  tgt_flx = fltarr(n_data) ; Will contain the target flux in Jy
  ob_out = intarr(n_data) ; Will contain OBs to be removed

  ; Print some info to screen and log file
  if info gt 0 then begin
    print, '==== REDUCTION INFO ===='
    print, ''
    print, 'File version number     : ' + string(file_ver, format = '(I0)')
    print, 'Input aperture radius   : ' + string(aper_rad, format = '(I0)') + '  (0 for EEID + FWHM)'
    ; PRINT, 'Background inner radius : ' + STRING(bck_irad, FORMAT='(I0)')
    ; PRINT, 'Background outer radius : ' + STRING(bck_orad, FORMAT='(I0)')
    print, 'Image combination mode  : ' + string(img_mode, format = '(I0)') + '  (0: for automatic)'
    print, 'Background nod mode     : ' + string(bck_mode, format = '(I0)') + '  (0: none, 1 nod pairs, 2 closest n frames, 3 dedicated, 3 chopping)'
    print, 'Background floor mode   : ' + string(bfl_mode, format = '(I0)') + '  (0: 5-sigma clipped, 1 median)'
    print, 'Frame mode for flux     : ' + string(fra_mode, format = '(I0)') + '  (0: mean-background subtracted, 1: processed raw, 2: pca-background subtracted)'
    print, 'Flux mode               : ' + string(flx_mode, format = '(I0)') + '  (0: aper. phot., 1: weigh. aper phot., 2: PSF-fitting)'
    print, 'Centroid mode           : ' + string(fit_mode, format = '(I0)') + '  (0: none, 1: CNTRD., 2 GCNTRD, 3: Gaussian, 4: Lorentzan, 5: Moffat)'
    print, 'Null estimator          : ' + string(drs.null_est, format = '(I0)')
    print, 'Null error mode         : ' + string(drs.err_mode, format = '(I0)')
    print, 'Scatter-rejection sigma : ' + string(drs.sig_sca, format = '(I0)')
    print, 'Null correction mode    : ' + string(null_cor, format = '(I0)') + '  (1: Denis'', 2: Bertrand''s method to subtract HF noise)'
    print, 'NSC cube size           : ' + nsc_cube + '  (0 for automatic)'
    print, 'NSC mode                : ' + string(nsc_mode, format = '(I0)') + '  (0: mode, 1: %best, 2: NSC)'
    print, 'Acceptable null range   : ' + null_lim
    print, 'Number of bin factor    : ' + string(nbin_fac, format = '(I0)') + '  (Bin size for NSC reduction (0: constant, 1: variable)'
    print, 'Number of bootstrap     : ' + string(n_btstrp, format = '(I0)')
    print, 'Calibration method      : ' + string(drs.cal_method, format = '(I0)') + '  (0: interpolation nearest neighbors, 2: interpolation on all OBs)'
    print, 'Degree of polynomial    : ' + string(drs.polydeg, format = '(I0)')
    print, 'Calibration mode        : ' + string(drs.cal_mode, format = '(I0)') + '  (0: per pointing, 1: per OB)'
    print, 'Reduced chi2 limit      : ' + string(drs.chi2_lim, format = '(i0)')
    print, 'TF error mode           : ' + string(drs.err_mode, format = '(I0)') + '  (0: unweighted null dispersion, 1: weighted null dispersion, 2: excess variance)'
    print, 'Nod correlation         : ' + string(drs.nod_cor, format = '(I0)') + '  (1 if nulls appear correlated per nod within a given pointing)'
    print, 'Split TF                : ' + string(drs.split_tf, format = '(I0)') + '  (1 to split the TF when there is a dead time longer than split_time between two null measurements)'
    print, 'Split time              : ' + string(drs.split_time, format = '(F3.1)')
    print, 'Split hour              : ' + string(drs.split_hour, format = '(I0)')
    print, 'Split nod               : ' + string(drs.split_nod, format = '(I0)')
    print, ''
    print, 'Number of OBs           : ' + string(n_data, format = '(I03)')
    print, 'Number of pointings     : ' + string(n_pt, format = '(I03)')
    ; PRINT, 'PLC wavelength          : ' + STRING(plc_wav, FORMAT='(F5.2)')
    print, ''
    print, ''
    ; Column signification
    print, 'Column signification'
    print, 'A: OB identification number'
    print, 'B: Object name'
    print, 'C: Flag (SCI/CAL)'
    print, 'D: UT time (from LBT)'
    print, 'E: Telescope elevation in degrees (from LBT)'
    print, 'F: Parallactic angle in degrees (from LBT)'
    print, 'G: DIMM seeing [arcsec]'
    print, 'H: PWV [mm]'
    print, 'I: DIT [ms]'
    print, 'J: Dithering period [number of frames]'
    print, 'K: Wavelength [microns]'
    print, 'L: Bandwidths [microns]'
    print, 'M: SNR on photometry (DX and SX)'
    print, 'N: Intensity mismatch (in %)'
    print, 'O: Tip/tilt error (DX and SX, in mas)'
    print, 'P: RMS of photometry relative to constructive peak (DX and SX, in %)'
    print, 'Q: RMS of background relative to constructive peak (in %)'
    print, 'R: Background bias relative to constructive peak (in %)'
    print, 'S: Null (in %)'
    print, 'T: Null offset between optimum NSC value and bootstrapped or bayesian value (in %)'
    print, 'U: Null uncertainty due to error on constructive peak (in %)'
    print, 'V: Null uncertainty due to external error on background floor (in %)'
    print, 'W: Null uncertainty from NSC fit (in %)'
    print, 'X: Total null uncertainty (in %)'
    print, 'Y: Expected null uncertainty due to photometric errors (in %)'
    print, 'Z: Expected null uncertainty due to phase variation (i.e., sqrt(4mu^2.sig^2 + 2*sig^4)/4, in %)'
    print, 'a: Best-fit mean PHASE from NSC (in microns)'
    print, 'b: Best-fit PHASE jitter from NSC (in microns)'
    print, 'c: Reduced chi2 from NSC'
    print, 'd: Number of rejected null frames'
    print, 'e: Preserved number of frames'
    print, ' '
    print, ' A ', '  ', '    B   ', '  ', ' C ', '  ', '     D     ', '  ', '   E  ', '  ', '    F  ', '  ', '  G ', '  ', '  H ', '  ', '  I  ', '  ', '   J ', '  ', '  K  ', '  ', '      L     ', '  ', '   M   ', '  ', '     N     ', '  ', '    O     ', '  ', '  P  ', '  ', '  Q  ', '  ', '    R  ', ' ', '   S  ', '  ', '    T ', '  ', '   U  ', ' ', '   V  ', '  ', '   W  ', '   ', '   X  ', ' ', '  Y  ', '  ', '   Z  ', '  ', '   a   ', '  ', '  b  ', '  ', '  c  ', '  ', '  d  ', '  ', '  e  ', '  '
    print, '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    for i = 0, n_data - 1 do print, string(l1data[i].ob_id, format = '(I03)'), '  ', string(objname[i], format = '(A8)'), '  ', string(flag[i], format = '(A3)'), '  ', string(utc[i], format = '(A11)'), '  ', string(alt[i], format = '(F6.2)'), '  ', string(para[i], format = '(F7.2)'), '  ', string(seeing[i], format = '(F4.2)'), '  ', string(smttau[i], format = '(F4.2)'), $
      '  ', string(1d+3 * dit[i], format = '(F5.1)'), '  ', string(dith_per[i], format = '(I04)'), '  ', string(1d+6 * wav[i], format = '(F5.1)'), '  ', string(1d+6 * bdw[i], format = '(F4.1)'), ' | ', string(l1data[i].photdx_snr, format = '(F5.1)'), ' ', string(l1data[i].photsx_snr, format = '(F5.1)'), $
      '  ', string(int_err[i], format = '(F7.4)'), '  ', string(l1data[i].ttdx_rms, format = '(F5.1)'), ' ', string(l1data[i].ttsx_rms, format = '(F5.1)'), '  ', string(100. * l1data[i].photdx_rms / l1data[i].phot_avg, format = '(F4.2)'), '  ', string(100. * l1data[i].photsx_rms / l1data[i].phot_avg, format = '(F4.2)'), $
      '  ', string(100. * l1data[i].bckg_cnod_rms, format = '(F5.2)'), '  ', string(100. * l1data[i].bckg_bias, format = '(F5.2)'), ' || ', string(100. * l1data[i].null_meas, format = '(F5.2)'), '  ', string(100. * l1data[i].null_offset, format = '(F5.2)'), ' || ', string(100. * err_phot[i], format = '(F5.2)'), '  ', string(100. * err_bckg[i], format = '(F5.2)'), '  ', string(100. * err_nsc[i], format = '(F5.2)'), $
      ' | ', string(100. * l1data[i].null_meas_err, format = '(F5.2)'), ' || ', string(100. * null_err_phot[i], format = '(F5.2)'), ' ', string(100. * null_err_phase[i], format = '(F5.2)'), ' ||', string(1d6 * l1data[i].nsc_phavg * wav[i] / (2 * !dpi), format = '(F6.2)'), '  ', string(1d6 * l1data[i].nsc_phrms * wav[i] / (2 * !dpi), format = '(F6.2)'), '  ', string(l1data[i].nsc_chi2, format = '(F6.2)'), $
      '   ', string(l1data[i].nfr_rej, format = '(I04)'), '  ', string(l1data[i].nfr_ob, format = '(I04)')
    print, ''
    if lun gt 0 then begin
      printf, lun, '==== REDUCTION INFO ===='
      printf, lun, ''
      printf, lun, 'File version number     : ' + string(file_ver, format = '(I0)')
      printf, lun, 'Input aperture radius   : ' + string(aper_rad, format = '(I0)') + '  (0 for EEID + FWHM)'
      ; PRINTF, lun, 'Background inner radius : ' + STRING(bck_irad, FORMAT='(I0)')
      ; PRINTF, lun, 'Background outer radius : ' + STRING(bck_orad, FORMAT='(I0)')
      printf, lun, 'Image combination mode  : ' + string(img_mode, format = '(I0)') + '  (0 for automatic)'
      printf, lun, 'Background nod mode     : ' + string(bck_mode, format = '(I0)') + '  (0 none, 1 nod pairs, 2 closest n frames, 3 dedicated, 3 chopping)'
      printf, lun, 'Background floor mode   : ' + string(bfl_mode, format = '(I0)') + '  (0 5-sigma clipped, 1 median)'
      printf, lun, 'Frame mode for flux     : ' + string(fra_mode, format = '(I0)') + '  (0 background subtracted, 1: processed raw)'
      printf, lun, 'Flux mode               : ' + string(flx_mode, format = '(I0)') + '  (0 aper. phot., 1: weigh. aper phot., 2: PSF-fitting)'
      printf, lun, 'Centroid mode           : ' + string(fit_mode, format = '(I0)') + '  (0 none, 1: CNTRD., 2 GCNTRD, 3: Gaussian, 4: Lorentzan, 5: Moffat)'
      printf, lun, 'Null estimator          : ' + string(drs.null_est, format = '(I0)')
      printf, lun, 'Null error mode         : ' + string(drs.err_mode, format = '(I0)')
      printf, lun, 'Scatter-rejection sigma : ' + string(drs.sig_sca, format = '(I0)')
      printf, lun, 'Null correction mode    : ' + string(null_cor, format = '(I0)') + '  (1: Denis'', 2: Bertrand''s method to subtract HF noise)'
      printf, lun, 'NSC cube size           : ' + nsc_cube + '  (0 for automatic)'
      printf, lun, 'NSC mode                : ' + string(nsc_mode, format = '(I0)') + '  (0: mode, 1: %best, 2: NSC)'
      printf, lun, 'Acceptable null range   : ' + null_lim
      printf, lun, 'Number of bin factor    : ' + string(nbin_fac, format = '(I0)') + '  (Bin size for NSC reduction (0: constant, 1: variable)'
      printf, lun, 'Number of bootstrap     : ' + string(n_btstrp, format = '(I0)')
      printf, lun, 'Calibration method      : ' + string(drs.cal_method, format = '(I0)') + '  (0: interpolation nearest neighbors, 2: interpolation on all OBs)'
      printf, lun, 'Degree of polynomial    : ' + string(drs.polydeg, format = '(I0)')
      printf, lun, 'Calibration mode        : ' + string(drs.cal_mode, format = '(I0)') + '  (0: per pointing, 1: per OB)'
      printf, lun, 'Reduced chi2 limit      : ' + string(drs.chi2_lim, format = '(i0)')
      printf, lun, 'TF error mode           : ' + string(drs.err_mode, format = '(I0)') + '  ( 0: unweighted null dispersion, 1: weighted null dispersion, 2: excess variance)'
      printf, lun, 'Nod correlation         : ' + string(drs.nod_cor, format = '(I0)') + '  ( 1 if nulls appear correlated per nod within a given pointing)'
      printf, lun, 'Split TF                : ' + string(drs.split_tf, format = '(I0)') + '  (1 to split the TF when there is a dead time longer than split_time between two null measurements)'
      printf, lun, 'Split time              : ' + string(drs.split_time, format = '(F3.1)')
      printf, lun, 'Split hour              : ' + string(drs.split_hour, format = '(I0)')
      printf, lun, 'Split nod               : ' + string(drs.split_nod, format = '(I0)')
      printf, lun, ''
      printf, lun, 'Number of OBs           : ' + string(n_data, format = '(I03)')
      printf, lun, 'Number of pointings     : ' + string(n_pt, format = '(I03)')
      ; PRINTF, lun, 'PLC wavelength          : ' + STRING(plc_wav, FORMAT='(F5.2)')
      printf, lun, ''
      printf, lun, ''
      printf, lun, '==== OB-BASED INFORMATION ===='
      ; Column signification
      printf, lun, 'Column signification'
      printf, lun, 'A: OB identification number'
      printf, lun, 'B: Object name'
      printf, lun, 'C: Flag (SCI/CAL)'
      printf, lun, 'D: UT time (from LBT)'
      printf, lun, 'E: Telescope elevation in degrees (from LBT)'
      printf, lun, 'F: Parallactic angle in degrees (from LBT)'
      printf, lun, 'G: DIMM seeing [arcsec]'
      printf, lun, 'H: PWV [mm]'
      printf, lun, 'I: DIT [ms]'
      printf, lun, 'J: Dithering period [number of frames]'
      printf, lun, 'K: Wavelength [microns]'
      printf, lun, 'L: Bandwidths [microns]'
      printf, lun, 'M: SNR on photometry (DX and SX)'
      printf, lun, 'N: Intensity mismatch (in %)'
      printf, lun, 'O: Tip/tilt error (DX and SX, in mas)'
      printf, lun, 'P: RMS of photometry relative to constructive peak (DX and SX, in %)'
      printf, lun, 'Q: RMS of background relative to constructive peak (in %)'
      printf, lun, 'R: Background bias relative to constructive peak (in %)'
      printf, lun, 'S: Null (in %)'
      printf, lun, 'T: Null offset between optimum NSC value and bootstrapped or bayesian value (in %)'
      printf, lun, 'U: Null uncertainty due to error on constructive peak (in %)'
      printf, lun, 'V: Null uncertainty due to external error on background floor (in %)'
      printf, lun, 'W: Null uncertainty from NSC fit (in %)'
      printf, lun, 'X: Total null uncertainty (in %)'
      printf, lun, 'Y: Expected null uncertainty due to photometric errors (in %)'
      printf, lun, 'Z: Expected null uncertainty due to phase variation (i.e., sqrt(4mu^2.sig^2 + 2*sig^4)/4, in %)'
      printf, lun, 'a: Best-fit mean PHASE from NSC (in microns)'
      printf, lun, 'b: Best-fit PHASE jitter from NSC (in microns)'
      printf, lun, 'c: Reduced chi2 from NSC'
      printf, lun, 'd: Number of rejected null frames'
      printf, lun, 'e: Preserved number of frames'
      printf, lun, ' '
      printf, lun, ' A ', '  ', '    B   ', '  ', ' C ', '  ', '     D     ', '  ', '   E  ', '  ', '    F  ', '  ', '  G ', '  ', '  H ', '  ', '  I  ', '  ', '   J ', '  ', '  K  ', '  ', '      L     ', '  ', '   M   ', '  ', '     N     ', '  ', '    O     ', '  ', '  P  ', '  ', '  Q  ', '  ', '    R  ', ' ', '   S  ', '  ', '    T ', '  ', '   U  ', ' ', '   V  ', '  ', '   W  ', '   ', '   X  ', ' ', '  Y  ', '  ', '   Z  ', '  ', '   a   ', '  ', '  b  ', '  ', '  c  ', '  ', '  d  ', '  ', '  e  ', '  '
      printf, lun, '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
      for i = 0, n_data - 1 do printf, lun, string(l1data[i].ob_id, format = '(I03)'), '  ', string(objname[i], format = '(A8)'), '  ', string(flag[i], format = '(A3)'), '  ', string(utc[i], format = '(A11)'), '  ', string(alt[i], format = '(F6.2)'), '  ', string(para[i], format = '(F7.2)'), '  ', string(seeing[i], format = '(F4.2)'), '  ', string(smttau[i], format = '(F4.2)'), $
        '  ', string(1d+3 * dit[i], format = '(F5.1)'), '  ', string(dith_per[i], format = '(I04)'), '  ', string(1d+6 * wav[i], format = '(F5.1)'), '  ', string(1d+6 * bdw[i], format = '(F4.1)'), ' | ', string(l1data[i].photdx_snr, format = '(F5.1)'), ' ', string(l1data[i].photsx_snr, format = '(F5.1)'), $
        '  ', string(int_err[i], format = '(F7.4)'), '  ', string(l1data[i].ttdx_rms, format = '(F5.1)'), ' ', string(l1data[i].ttsx_rms, format = '(F5.1)'), '  ', string(100. * l1data[i].photdx_rms / l1data[i].phot_avg, format = '(F4.2)'), '  ', string(100. * l1data[i].photsx_rms / l1data[i].phot_avg, format = '(F4.2)'), $
        '  ', string(100. * l1data[i].bckg_cnod_rms, format = '(F5.2)'), '  ', string(100. * l1data[i].bckg_bias, format = '(F5.2)'), ' || ', string(100. * l1data[i].null_meas, format = '(F5.2)'), '  ', string(100. * l1data[i].null_offset, format = '(F5.2)'), ' || ', string(100. * err_phot[i], format = '(F5.2)'), '  ', string(100. * err_bckg[i], format = '(F5.2)'), '  ', string(100. * err_nsc[i], format = '(F5.2)'), $
        ' | ', string(100. * l1data[i].null_meas_err, format = '(F5.2)'), ' || ', string(100. * null_err_phot[i], format = '(F5.2)'), ' ', string(100. * null_err_phase[i], format = '(F5.2)'), ' ||', string(1d6 * l1data[i].nsc_phavg * wav[i] / (2 * !dpi), format = '(F6.2)'), '  ', string(1d6 * l1data[i].nsc_phrms * wav[i] / (2 * !dpi), format = '(F6.2)'), '  ', string(l1data[i].nsc_chi2, format = '(F6.2)'), $
        '   ', string(l1data[i].nfr_rej, format = '(I04)'), '  ', string(l1data[i].nfr_ob, format = '(I04)')
      printf, lun, ''
    endif
  endif

  ; SUBTRACT GEOMETRIC NULL FROM CALIBRATOR NULL
  ; ********************************************

  ; Find index of calibrator data
  idx_cal = where(strmatch(flag, 'CAL') eq 1, n_caldat)
  if n_caldat le 0 then begin
    print, 'No calibrator data => no calibration'
    if lun gt 0 then printf, lun, 'No calibrator data => no calibration'
    goto, skip_cal
  endif

  ; Derive number of unique wavelength
  cal_wav = wav[idx_cal]
  wav_uniq = cal_wav[uniq(cal_wav, sort(cal_wav))]
  n_wav = n_elements(wav_uniq)

  ; Loop over the wavelengths
  for i_wav = 0, n_wav - 1 do begin
    ; Extract data of this wavelength
    idx_wav = idx_cal[where(cal_wav eq wav_uniq[i_wav])]
    calname = objname[idx_wav]
    calbdw = bdw[idx_wav]

    ; Check whether there is only one bandwidth for this wavelength
    bdw_uniq = calbdw[uniq(calbdw, sort(calbdw))]
    if n_elements(bdw_uniq) gt 1 then message, 'Found two different bandwidths for the same wavelength'

    ; Compute wavelength range
    lambda_chnl = wav_uniq[i_wav] - 0.5 * bdw_uniq[0] + bdw_uniq[0] * (1. + 1. / n_lam) * dindgen(n_lam) / (n_lam - 1)

    ; Extract LBTI transmission profile for this wavalength (N' and 8.7um only now)
    if wav_uniq[i_wav] ge 1.1d-5 and wav_uniq[i_wav] le 1.2d-5 then begin
      READ_TABLE, 'nodrs/input/n-band_thruput.txt', lam_tmp, thruput, first = 4, separator = 'TAB'
      trans = interpol(thruput, lam_tmp, lambda_chnl * 1d6)
    endif else begin
      if wav_uniq[i_wav] ge 8.6d-6 and wav_uniq[i_wav] lt 8.8d-6 then begin
        READ_TABLE, 'nodrs/input/f87_thruput.txt', lam_tmp, thruput, first = 4, separator = 'TAB'
        trans = interpol(thruput, lam_tmp, lambda_chnl * 1d6)
      endif else begin
        if wav_uniq[i_wav] ge 8.8d-6 and wav_uniq[i_wav] le 9.0d-6 then begin
          READ_TABLE, 'nodrs/input/n08909_thruput.txt', lam_tmp, thruput, first = 4, separator = 'TAB'
          trans = interpol(thruput, lam_tmp, lambda_chnl * 1d6)
        endif else begin
          if wav_uniq[i_wav] ge 12d-6 and wav_uniq[i_wav] le 12.5d-6 then begin
            READ_TABLE, 'nodrs/input/n12520_thruput.txt', lam_tmp, thruput, first = 4, separator = 'TAB'
            trans = interpol(thruput, lam_tmp, lambda_chnl * 1d6)
          endif else trans = 1
        endelse
      endelse
    endelse

    ; Compute effective wavelength
    lam_eff = total(trans * lambda_chnl) / total(trans)

    ; Compute unique calibrators for this wavelength
    nam_uniq = calname[uniq(calname, sort(calname))]
    n_cal = n_elements(nam_uniq)

    ; Loop over calibrator data
    for i_cal = 0, n_cal - 1 do begin
      ; Extract data for this calibrator
      idx_cur = idx_wav[where(calname eq nam_uniq[i_cal])]

      ; Retrieve target information
      GET_TGT, nam_uniq[i_cal], star, database = pth.input_path + drs.database
      th_mean = star.ldm
      th_min = star.ldm - star.lde
      th_max = star.ldm + star.lde
      uld = star.uld

      ; Generate stellar spectrum for broadband null computation
      star_flux = BLACKBODY(star.temp, lambda_chnl, standard = 3) ; flux in photons/s/m²/Hz/sr
      flux_int = total(trans * star_flux * !dpi * (0.5 * th_mean * prm.m2R) ^ 2) / total(trans) ; Integrated flux in ph/s/m2/Hz
      tgt_flx[idx_cur] = flux_int * prm.h * prm.c * prm.w2Jy / lam_eff ; Integrated flux in Jansky

      ; Compute the expected null for this star
      norm = total((trans * star_flux) ^ 2)
      null_lam = GEOMETRIC_NULL(th_mean, cnf.base, lambda_chnl, uld = uld)
      null_wb = total((trans * star_flux) ^ 2 * null_lam) / norm
      ; Compute minimum null (with minimum diameter)
      null_min = GEOMETRIC_NULL(th_min, cnf.base, lambda_chnl, uld = uld)
      null_wb_min = total((trans * star_flux) ^ 2 * null_min) / norm
      ; Compute maximum null (with maximum diameter)
      null_max = GEOMETRIC_NULL(th_max, cnf.base, lambda_chnl, uld = uld)
      null_wb_max = total((trans * star_flux) ^ 2 * null_max) / norm
      ; Derive the error bar on the estimated null
      null_err_sys = (null_wb_max - null_wb_min) / 2d0

      ; Derive the wide-band instrumental transfer function and its error bar
      ; Statistical error bar remains unchanged
      tf_wb[idx_cur] -= null_wb ; correct for geometric leakage
      tf_wb_esyst[idx_cur] = null_err_sys ; systematic error bar
      tf_wb_etot[idx_cur] = sqrt(tf_wb_estat[idx_cur] ^ 2 + tf_wb_esyst[idx_cur] ^ 2) ; total error bar
    endfor
  endfor

  ; COMPUTE MEAN NULL AND ERROR PER POINTING/NOD
  ; ********************************************

  ; Number of unique pointings
  if info gt 0 then begin
    print, '==== POINTING-BASED NULLS (NO TIME GATING!) ===='
    print, '   '
    print, '  1. Instrumental null floor per pointing (pointing, name, flag, null[%], w_null[%], stat error[%], disp error[%], w_disp error[%], excess error [%])  '
    print, '   '
    if lun gt 0 then begin
      printf, lun, '==== POINTING-BASED NULLS (NO TIME GATING!) ===='
      printf, lun, '   '
      printf, lun, '  1. Instrumental null floor per pointing (pointing, name, flag, null[%], w_null[%], stat error[%], disp error[%], w_disp error[%], excess error [%])'
      printf, lun, '   '
    endif
  endif

  ; Compute weigthed instrumental null per pointing and per nod
  flag_pt = strarr(n_pt)
  rms_exc = fltarr(n_pt)
  time_avg_pt = fltarr(n_pt, n_nod)
  null_avg_pt = fltarr(n_pt, n_nod)
  null_err_pt = fltarr(n_pt, n_nod)
  null_err_stat_pt = fltarr(n_pt, n_nod)
  null_err_disp_pt = fltarr(n_pt, n_nod)
  null_err_wdisp_pt = fltarr(n_pt, n_nod)
  for i_pt = 0, n_pt - 1 do begin
    ; Derive OB of this pointing and parse stats results
    idx_pt = where(pt_id eq pt_uniq[i_pt], n_ob)
    AVGERR, tf_wb[idx_pt], tf_wb_estat[idx_pt], avg_unw, avg_wei, err_stat, err_sca_unw, err_sca_wei, exc_rms, kappa = drs.sig_sca
    flag_pt[i_pt] = flag[idx_pt[0]]
    rms_exc[i_pt] = exc_rms
    ; Now loop over the nods
    for i_nod = 0, n_nod - 1 do begin
      ; Derive nod of this pointing
      idx_npt = idx_pt[where(nod_pos[idx_pt] eq nod_uniq[i_nod], n_tmp)]
      ; Compute mean time
      time_avg_pt[i_pt, i_nod] = mean(time[idx_npt])
      ; Compute average value and error bars
      AVGERR, tf_wb[idx_npt], tf_wb_estat[idx_npt], avg_unw, avg_wei, err_stat, err_sca_unw, err_sca_wei, sys_ecc, kappa = drs.sig_sca
      if drs.null_est eq 0 then null_est = avg_unw else null_est = avg_wei
      null_avg_pt[i_pt, i_nod] = null_est
      null_err_stat_pt[i_pt, i_nod] = err_stat
      null_err_disp_pt[i_pt, i_nod] = err_sca_unw
      null_err_wdisp_pt[i_pt, i_nod] = err_sca_wei
      ; Systematic error
      case (abs(drs.err_mode) mod 3) of
        0: rms_sys = null_err_disp_pt[i_pt, i_nod]
        1: rms_sys = null_err_wdisp_pt[i_pt, i_nod]
        2: rms_sys = sys_ecc
        else: message, 'Undefined error mode'
      endcase
      ; Compute total error (use maximum of the two if greater than 2)
      if drs.err_mode ge 3 then err_tot = null_err_stat_pt[i_pt, i_nod] > rms_sys else err_tot = sqrt(null_err_stat_pt[i_pt, i_nod] ^ 2 + rms_sys ^ 2)
      null_err_pt[i_pt, i_nod] = sqrt(err_tot ^ 2 + mean(tf_wb_esyst[idx_npt]) ^ 2) ; in most cases, the systematic error due to the uncertainty of the stellar diameter is the same for all points
      ; Print info to screen
      if info gt 0 then begin
        print, string(i_pt, format = '(I02)'), '  ', string(objname[idx_npt[0]], format = '(A8)'), '  ', flag[idx_npt[0]], '  ', string(1d2 * avg_unw, format = '(F7.4)'), '  ', string(1d2 * avg_wei, format = '(F7.4)'), '  ', string(1d2 * null_err_stat_pt[i_pt, i_nod], format = '(F6.4)'), '  ', string(1d2 * err_sca_unw, format = '(F6.4)'), '  ', string(1d2 * err_sca_wei, format = '(F6.4)'), '  ', string(1d2 * sys_ecc, format = '(F6.4)')
        if lun gt 0 then printf, lun, string(i_pt, format = '(I02)'), '  ', string(objname[idx_npt[0]], format = '(A8)'), '  ', flag[idx_npt[0]], '  ', string(1d2 * avg_unw, format = '(F7.4)'), '  ', string(1d2 * avg_wei, format = '(F7.4)'), '  ', string(1d2 * null_err_stat_pt[i_pt, i_nod], format = '(F6.4)'), '  ', string(1d2 * err_sca_unw, format = '(F6.4)'), '  ', string(1d2 * err_sca_wei, format = '(F6.4)'), '  ', string(1d2 * sys_ecc, format = '(F6.4)')
      endif
    endfor
  endfor

  ; Merge nods if requested
  if drs.split_nod eq 0 and n_nod gt 1 then begin
    if info gt 0 then begin
      print, '   '
      print, '  2. Instrumental null floor per pointing (pointing, null[%], stat error[%], disp error[%], w_disp error[%], excess error [%], tot error [%])'
      if lun gt 0 then printf, lun, '   '
      if lun gt 0 then printf, lun, '  2. Instrumental null floor per pointing (pointing, null[%], stat error[%], disp error[%], w_disp error[%], excess error [%], tot error[%])  '
    endif
    time_avg_pt = mean(time_avg_pt, dimension = 2)
    null_err_stat_pt = mean(null_err_stat_pt, dimension = 2) / sqrt(n_nod)
    null_avg_pt0 = null_avg_pt
    null_err_pt0 = null_err_pt
    null_err_disp_pt0 = null_err_disp_pt
    null_err_wdisp_pt0 = null_err_wdisp_pt
    null_avg_pt = fltarr(n_pt)
    null_err_pt = fltarr(n_pt)
    null_err_disp_pt = fltarr(n_pt)
    null_err_wdisp_pt = fltarr(n_pt)
    for i_pt = 0, n_pt - 1 do begin
      ; Compute dispersion errors
      AVGSDV, null_avg_pt0[i_pt, *], avg_tmp1, rms_disp, rmsm_disp, kappa = drs.sig_sca
      AVGSDV, null_avg_pt0[i_pt, *], avg_tmp2, rms_wdisp, rmsm_wdisp, kappa = drs.sig_sca, weight = 1. / null_err_pt0[i_pt, *] ^ 2
      if drs.null_est eq 0 then null_est = avg_tmp1 else null_est = avg_tmp2
      null_avg_pt[i_pt] = null_est
      null_err_disp_pt[i_pt] = rmsm_disp
      null_err_wdisp_pt[i_pt] = rmsm_wdisp
      ; Systematic error. RMS_SYS is added quadratically to the statiscal error, unless err_mode is greater than 2 (in that case, the maximum of the 2 is used)
      case (abs(drs.err_mode) mod 3) of
        0: rms_sys = null_err_disp_pt[i_pt]
        1: rms_sys = null_err_wdisp_pt[i_pt]
        2: rms_sys = rms_exc[i_pt]
        else: message, 'Undefined error mode'
      endcase
      ; Compute total error (use maximum of the two if greater than 2)
      if drs.err_mode ge 3 then null_err_pt[i_pt] = null_err_stat_pt[i_pt] > rms_sys else null_err_pt[i_pt] = sqrt(null_err_stat_pt[i_pt] ^ 2 + rms_sys ^ 2)
      if info gt 0 then begin
        print, string(i_pt, format = '(I02)'), '  ', string(1d2 * null_avg_pt[i_pt], format = '(F7.4)'), '  ', string(1d2 * null_err_stat_pt[i_pt], format = '(F6.4)'), '  ', string(1d2 * null_err_disp_pt[i_pt], format = '(F6.4)'), '  ', string(1d2 * null_err_wdisp_pt[i_pt], format = '(F6.4)'), '  ', string(1d2 * rms_exc[i_pt], format = '(F6.4)'), '  ', string(1d2 * null_err_pt[i_pt], format = '(F6.4)')
        if lun gt 0 then printf, lun, string(i_pt, format = '(I02)'), '  ', string(1d2 * null_avg_pt[i_pt], format = '(F7.4)'), '  ', string(1d2 * null_err_stat_pt[i_pt], format = '(F6.4)'), '  ', string(1d2 * null_err_disp_pt[i_pt], format = '(F6.4)'), '  ', string(1d2 * null_err_wdisp_pt[i_pt], format = '(F6.4)'), '  ', string(1d2 * rms_exc[i_pt], format = '(F6.4)'), '  ', string(1d2 * null_err_pt[i_pt], format = '(F6.4)')
      endif
    endfor
  endif

  ; Extract CAL pointings (for later)
  idx_cal_pt = where(strmatch(flag_pt, 'CAL') eq 1, n_caldat_pt)
  caltime = time[idx_cal]
  caltime_pt = time_avg_pt[idx_cal_pt, *]

  ; READ SCIENCE DATA
  ; *****************

  ; Jumping point if no calibrator data
  skip_cal:

  ; Search for science target fits files in the input folder
  idx_sci = where(strmatch(flag, 'SCI') eq 1, n_scidat)
  if n_scidat gt 0 then scidata = l1data[idx_sci] else goto, skip_interp
  idx_sci_pt = where(strmatch(flag_pt, 'SCI') eq 1, n_scidat_pt)

  ; Extract SCI times
  scitime = time[idx_sci]
  scitime_pt = time_avg_pt[idx_sci_pt, *]
  sci_pt = pt_uniq[idx_sci_pt]

  ; Compute number of nod positions to consider (otherwise, it uses nod_pos as defined above)
  if drs.split_nod eq 0 then nod_pos = intarr(n_data)
  nod_uniq = nod_pos[uniq(nod_pos, sort(nod_pos))]
  n_nod = n_elements(nod_uniq)

  ; Define output arrays
  if drs.cal_mode eq 1 then n_out = n_scidat else n_out = n_scidat_pt
  if drs.cal_mode eq 1 then n_bis = 1 else n_bis = n_nod
  tf_sci = dblarr(n_out, n_bis)
  tf_estat = tf_sci
  tf_esyst = tf_sci
  tf_etot = tf_sci ; will contain the estimate of the TF at the science times (and error bars)
  null_sci = dblarr(n_out, n_bis)
  null_sci_estat = null_sci
  null_sci_esyst = null_sci
  null_sci_etot = null_sci ; will contain the final calibrated null and corresponding error bars

  ; INTERPOLATE THE TF OVER THE NIGHT
  ; *********************************

  ; Define the number of sub-TFs based on maximum dead time allowed
  ; This is badly coded...each element of idx_split gives the beginning of the new sub TF
  if drs.split_tf ne 0 and drs.split_time ne 0 then begin
    idx_split = where(abs(l1data[1 : *].mjd_obs - l1data[0 : (n_data - 1)].mjd_obs) * 24. gt drs.split_time, n_tf)
    if n_tf gt 0 then idx_split = [0, idx_split + 1, n_data] else idx_split = [0, n_data]
  endif else idx_split = [0, n_data]

  ; Add sub-TFs based on user-defined hours at which the TF has to be split
  if drs.split_hour[0] ne 0 then begin
    idx_tmp = where((l1data.mjd_obs - t0day) * 24 - drs.split_hour gt 0, n_tf)
    idx_tmp = idx_tmp[0]
    if n_tf gt 0 then begin
      idx_split = [idx_split, idx_tmp]
      idx_split = idx_split[sort(idx_split)]
    endif
  endif

  ; Now split if aper radius changes (when it's EEID dependend)
  for i = 1, n_elements(aper) - 1 do if aper[i] ne aper[i - 1] then idx_split = [idx_split, i]
  idx_split = idx_split[uniq(idx_split, sort(idx_split))]

  ; Also split if aper radius changes
  ; Inner radius
  for i = 1, n_elements(birad) - 1 do if birad[i] ne birad[i - 1] then idx_split = [idx_split, i]
  idx_split = idx_split[uniq(idx_split, sort(idx_split))]
  ; Outer radius (no, don't use it anymore because stars are likely not on the same spot of the detector anyway if 'borad' is different)
  ; FOR i = 1, N_ELEMENTS(borad)-1 DO IF borad[i] NE borad[i-1] THEN idx_split = [idx_split, i]
  ; idx_split = idx_split[UNIQ(idx_split,SORT(idx_split))]

  ; Derive number of sub-TFs
  n_tf = n_elements(idx_split) - 1

  ; Prepare array for interpolated TF
  tf_time = fltarr(n_tf, n_time)
  tf_plot = tf_time
  tf_plot_err = fltarr(n_tf)
  tf_plot_pt = fltarr(n_tf, n_time, n_nod)
  tf_plot_pt_err = fltarr(n_tf, n_nod)

  ; Loop over the TFs
  for i_tf = 0, n_tf - 1 do begin
    ; Extract data of this sub-TF
    if i_tf ne n_tf - 1 then idx_tf = where(time ge time[idx_split[i_tf]] and time lt time[idx_split[i_tf + 1]], n_intf) $
    else idx_tf = where(time ge time[idx_split[i_tf]], n_intf)
    flag_tf = flag[idx_tf]
    idx_cal_tf = idx_tf[where(strmatch(flag_tf, 'CAL') eq 1, n_cal_tf)]
    idx_sci_tf = idx_tf[where(strmatch(flag_tf, 'SCI') eq 1, n_sci_tf)]

    ; Skip this TF if no science and no cal frames
    if n_sci_tf gt 0 and n_cal_tf gt 0 then begin
      ; Number of unique calibrators
      cal_name = objname[idx_cal_tf]
      cal_uniq = cal_name[uniq(cal_name, sort(cal_name))]
      n_cal = n_elements(cal_uniq)

      ; Define plot range
      min_time = 0.98 * min(time[idx_tf])
      max_time = 1.02 * max(time[idx_tf])
      plottime = dindgen(n_time) / (double(n_time) - 1d0) * (max_time - min_time) + min_time
      tf_time[i_tf, *] = plottime

      ; Remove bad points (see rules within AVGERR routine)
      AVGERR, tf_wb[idx_cal_tf], tf_wb_estat[idx_cal_tf], avg_null, avg_null_we, rms_stat_null, disp_m_null, disp_m_null_we, rms_exc, kappa = drs.sig_sca, idx_out = idx_out, n_out = n_out
      if n_out gt 0 then ob_out[idx_cal_tf[idx_out]] = 1
      idx_in = idx_cal_tf[where(ob_out[idx_cal_tf] ne 1)]
      n_cal_tf_f = n_elements(null_cal)

      ; Calibation by linear interpolation between neighbours (REMOVING bad points!)
      if drs.cal_method eq 0 then begin
        ; To be implemented
        message, 'Not yet implemented'
      endif else begin ; Calibation by polynomial interpolation on all calibrator data
        ; Interpolate the data at data and plot times
        polydeg = drs.polydeg
        if polydeg gt n_cal then message, 'Degree of polynomial fit should not exceed the number of calibrator data points.'
        polyfit = poly_fit(time[idx_in], tf_wb[idx_in], polydeg, /double, measure_errors = tf_wb_estat[idx_in] ^ drs.null_est, chisq = chi2, covar = covar, sigma = polysig, yfit = yfit, yband = yband, yerror = yerror) ; include only statistical error bars
        tf_data = 0d
        tf_tmp = 0d
        for ipol = 0, polydeg do begin
          tf_data += polyfit[ipol] * time[idx_tf] ^ ipol
          tf_tmp += polyfit[ipol] * plottime ^ ipol
        endfor
        tf_plot[i_tf, *] = tf_tmp
      endelse

      ; Calibrate null
      null_c = tf_wb[idx_tf] - tf_data

      ; Extract calibrator data
      null_cal = null_c[where(strmatch(flag_tf, 'CAL') eq 1, ntmp)]
      null_cal_err = tf_wb_estat[idx_cal_tf]

      ; Compute mean and errors on TF: OB BASED and also flag points too far from TF
      AVGERR, null_cal, null_cal_err, avg_null, avg_null_we, rms_stat_null, disp_m_null, disp_m_null_we, rms_exc, kappa = drs.sig_sca, idx_out = idx_out, n_out = n_out
      if n_out gt 0 then ob_out[idx_cal_tf[idx_out]] = 1
      n_cal_tf_f = n_elements(null_cal)

      ; Systematic error defined by the err_mode input parameter
      ; This term will be added to the statistcal error
      case (abs(drs.err_mode) mod 3) of
        0: rms_sys = disp_m_null
        1: rms_sys = disp_m_null_we
        2: rms_sys = rms_exc
        else: message, 'Undefined error mode'
      endcase

      ; Compute the systematic error on TF estimate (due to uncertainty on stellar diameter)
      ; Loop over the diffferent calibrators (systematic error is fully correlated between different OBs on the same calibrator)
      null_err_diam = 0d
      for ical = 0, n_cal - 1 do null_err_diam += (total(tf_wb_esyst[idx_cal_tf[where(cal_name eq cal_uniq[ical], n_c)]]) / n_c) ^ 2
      rms_sys_tot = sqrt(rms_sys ^ 2 + null_err_diam ^ 2)

      ; Compute total error (use maximum of the two if greater than 2)
      if drs.err_mode ge 3 then null_err_tot = rms_stat_null > rms_sys_tot else null_err_tot = sqrt(rms_stat_null ^ 2 + rms_sys_tot ^ 2)
      tf_plot_err[i_tf] = null_err_tot

      ; Ok, now save the calibrated null and null floor in output arrays if CAL_MODE is 1
      if drs.cal_mode eq 1 then begin
        idx_tmp = where(strmatch(flag_tf, 'SCI') eq 1, ntmp)
        for i_s = 0, n_sci_tf - 1 do begin
          idx_s = where(scitime eq time[idx_sci_tf[i_s]], n_ok)
          if n_ok le 0 then message, 'Problem with science time computation!!'
          null_sci[idx_s] = null_c[idx_tmp[i_s]]
          tf_sci[idx_s] = tf_data[i_s]
          tf_estat[idx_s] = rms_stat_null
          tf_esyst[idx_s] = rms_sys_tot
          tf_etot[idx_s] = null_err_tot
          null_sci_estat[idx_s] = sqrt(tf_wb_estat[idx_sci_tf[i_s]] ^ 2 + rms_stat_null ^ 2)
          null_sci_esyst[idx_s] = rms_sys_tot
          null_sci_etot[idx_s] = sqrt(null_sci_estat[idx_s] ^ 2 + null_sci_esyst[idx_s] ^ 2)
        endfor
      endif

      ; Print diagnostic info to screen and log
      if info gt 0 then begin
        print, ' '
        print, '==== SUB-TF COMPUTATION (OB-based approach) ==== '
        print, '================================================ '
        print, ' '
        print, '==== PART ' + string(i_tf + 1, format = '(I0)') + ' ==== '
        print, '  Statistical error on TF [%] :'
        print, '    - based on individual error bars    : ', string(rms_stat_null * 1d+2, format = '(F6.4)')
        print, '    - standard deviation on the mean    : ', string(disp_m_null_we * 1d+2, format = '(F6.4)') + ' (' + string(disp_m_null * 1d+2, format = '(F6.4)') + ' if unweighted)'
        ; PRINT, '    - quadratic sum of the 2            : ', STRING(null_err_tot*1D+2, FORMAT='(F6.4)')
        print, '  Systematic error on TF [%]            : '
        print, '    - due to diameter uncertainty       : ', string(null_err_diam * 1d+2, format = '(F6.4)')
        print, '    - sqrt(excess variance)             : ', string(rms_exc * 1d+2, format = '(F6.4)')
        print, '  Total error on TF [%]                 : ', string(null_err_tot * 1d+2, format = '(F6.4)')
        print, '  Number of initial CAL OBs             : ', string(n_cal_tf, format = '(I0)')
        print, '  Number of preserved CAL OBs           : ', string(n_cal_tf_f, format = '(I0)')
        print, '  Number of initial CAL frames          : ', string(total(nfr_ob[idx_cal_tf]), format = '(I0)') ; + ' (' + STRING(TOTAL(sci_fro), FORMAT='(I0)') + ',' + STRING(TOTAL(cal_fro), FORMAT='(I0)') +')'
        print, ' '
        if lun gt 0 then begin
          printf, lun, ' '
          printf, lun, '==== SUB-TF COMPUTATION (OB-based approach) ==== '
          printf, lun, '================================================ '
          printf, lun, ' '
          printf, lun, '==== PART ' + string(i_tf + 1, format = '(I0)') + ' ==== '
          printf, lun, '  Statistical error on TF [%] :'
          printf, lun, '    - based on individual error bars    : ', string(rms_stat_null * 1d+2, format = '(F6.4)')
          printf, lun, '    - standard deviation on the mean    : ', string(disp_m_null_we * 1d+2, format = '(F6.4)') + ' (' + string(disp_m_null * 1d+2, format = '(F6.4)') + ' if unweighted)'
          ; PRINTF, lun,  '    - quadratic sum of the 2            : ', STRING(null_err_tot*1D+2, FORMAT='(F6.4)')
          printf, lun, '  Systematic error on TF [%]            : '
          printf, lun, '    - due to diameter uncertainty       : ', string(null_err_diam * 1d+2, format = '(F6.4)')
          printf, lun, '    - sqrt(excess variance)             : ', string(rms_exc * 1d+2, format = '(F6.4)')
          printf, lun, '  Total error on TF [%]                 : ', string(null_err_tot * 1d+2, format = '(F6.4)')
          printf, lun, '  Number of initial CAL OBs             : ', string(n_cal_tf, format = '(I0)')
          printf, lun, '  Number of preserved CAL OBs           : ', string(n_cal_tf_f, format = '(I0)')
          printf, lun, '  Number of initial CAL frames          : ', string(total(nfr_ob[idx_cal_tf]), format = '(I0)') ; + ' (' + STRING(TOTAL(sci_fro), FORMAT='(I0)') + ',' + STRING(TOTAL(cal_fro), FORMAT='(I0)') +')'
          printf, lun, ' '
        endif
      endif

      ; Now do the pointing approach (i.e., calibrate based on the pointing values
      ; **************************************************************************
      ;
      ; Find pointings in this TF
      pt_id_tf = pt_id[idx_tf]
      nod_pos_tf = nod_pos[idx_tf]
      pt_uniq_tf = pt_id_tf[uniq(pt_id_tf, sort(pt_id_tf))]
      n_pt_tf = n_elements(pt_uniq_tf)
      idx_pt_tf = value_locate(pt_uniq, pt_uniq_tf)
      flag_pt_tf = flag_pt[idx_pt_tf]

      ; Prepare arrays
      time_avg_pt_tf = dblarr(n_pt_tf, n_nod)
      null_avg_pt_tf = time_avg_pt_tf
      null_err_pt_tf = time_avg_pt_tf

      ; Recompute pointing-based value because they might be different from the ones computed above (because of specific TF splitting)
      for ip = 0, n_pt_tf - 1 do begin
        ; Derive nod position for this TF
        idx_p = where(pt_id_tf eq pt_uniq_tf[ip])
        nod_pt_tf = nod_pos_tf[idx_p]
        nod_uniq_tf = nod_pt_tf[uniq(nod_pt_tf, sort(nod_pt_tf))]
        n_nod = n_elements(nod_uniq_tf)
        for i_nod = 0, n_nod - 1 do begin
          ; Extract nod data and parse to TF arrays (also ignore points to far from the TF)
          idx_pn = idx_tf[where(pt_id_tf eq pt_uniq_tf[ip] and nod_pos_tf eq nod_uniq_tf[i_nod] and ob_out[idx_tf] ne 1)]
          AVGERR, tf_wb[idx_pn], tf_wb_estat[idx_pn], avg_null, avg_null_we, rms_stat_null, disp_m_null, disp_m_null_we, rms_exc, kappa = drs.sig_sca, idx_out = idx_out, n_out = n_out
          time_avg_pt_tf[ip, i_nod] = mean(time[idx_pn])
          if drs.null_est eq 0 then null_est = avg_null else null_est = avg_null_we
          null_avg_pt_tf[ip, i_nod] = null_est
          if n_out gt 0 then ob_out[idx_pn[idx_out]] = 1

          ; Systematic error (according to err_mode definition)
          case (abs(drs.err_mode) mod 3) of
            0: rms_sys = disp_m_null
            1: rms_sys = disp_m_null_we
            2: rms_sys = rms_exc
            else: message, 'Undefined error mode'
          endcase
          ; Compute total error as sum of the two or just the maximum one (keep error due to stellar diameter in both case)
          if drs.err_mode ge 3 then null_err_pt_tf[ip, i_nod] = sqrt(rms_stat_null ^ 2 + mean(tf_wb_esyst[idx_pn]) ^ 2) > sqrt(mean(tf_wb_esyst[idx_pn]) ^ 2 + rms_sys ^ 2) else null_err_pt_tf[ip, i_nod] = sqrt(rms_stat_null ^ 2 + mean(tf_wb_esyst[idx_pn]) ^ 2 + rms_sys ^ 2)
          ; Parse to global pointing arrays
          idx_pc = where(pt_uniq eq pt_uniq_tf[ip])
          time_avg_pt[idx_pc, i_nod] = time_avg_pt_tf[ip, i_nod]
          null_avg_pt[idx_pc, i_nod] = null_avg_pt_tf[ip, i_nod]
          null_err_pt[idx_pc, i_nod] = null_err_pt_tf[ip, i_nod]
        endfor
      endfor

      ; Find CAL and SCI pointings
      idx_cal_pt_tf = where(strmatch(flag_pt_tf, 'CAL') eq 1, n_cal_pt_tf)
      idx_sci_pt_tf = where(strmatch(flag_pt_tf, 'SCI') eq 1, n_sci_pt_tf)

      ; Now loop over nods
      for i_nod = 0, n_nod - 1 do begin
        ; Extract value for this nod
        time_avg_pt_tf_cal = time_avg_pt_tf[idx_cal_pt_tf, i_nod]
        null_avg_pt_tf_cal = null_avg_pt_tf[idx_cal_pt_tf, i_nod]
        null_err_pt_tf_cal = null_err_pt_tf[idx_cal_pt_tf, i_nod]

        ; Do interpolation on pointing nulls
        if drs.cal_method eq 0 then begin
          ; To be implemented
          message, 'Not yet implemented'
        endif else begin ; Calibation by polynomial interpolation on all calibrator data
          ; Interpolate the data at data and plot times
          polydeg = drs.polydeg
          if polydeg gt n_cal then message, 'Degree of polynomial fit should not exceed the number of calibrator data points.'
          polyfit = poly_fit(time_avg_pt_tf_cal, null_avg_pt_tf_cal, polydeg, /double, measure_errors = null_err_pt_tf_cal ^ drs.null_est, chisq = chi2, covar = covar, sigma = polysig, yfit = yfit, yband = yband, yerror = yerror) ; include only statistical error bars
          tf_data = 0d
          tf_tmp = 0d
          for ipol = 0, polydeg do begin
            tf_data += polyfit[ipol] * time_avg_pt_tf[idx_pt_tf, i_nod] ^ ipol
            tf_tmp += polyfit[ipol] * plottime ^ ipol
          endfor
          tf_plot_pt[i_tf, *, i_nod] = tf_tmp
        endelse

        ; Calibrate null
        null_c_pt = null_avg_pt_tf[*, i_nod] - tf_data

        ; Extract calibrator data
        null_cal_pt = null_c_pt[where(strmatch(flag_pt_tf, 'CAL') eq 1, ntmp)]
        null_cal_err_pt = null_err_pt_tf_cal

        ; Compute propagated statistical error on TF: POINTING BASED and renove large residues
        AVGERR, null_cal_pt, null_cal_err_pt, avg_null, avg_null_we, rms_stat_null_pt, disp_m_null_pt, disp_m_null_pt_we, rms_exc, kappa = drs.sig_sca, idx_out = idx_out, n_out = n_out

        ; Definition of systematic errors
        case (abs(drs.err_mode) mod 3) of
          0: rms_sys_pt = disp_m_null_pt
          1: rms_sys_pt = disp_m_null_pt_we
          2: rms_sys_pt = rms_exc
          else: message, 'Undefined error mode'
        endcase

        ; Systematic error is the same as the one per OB
        rms_sys_tot_pt = sqrt(rms_sys_pt ^ 2 + null_err_diam ^ 2)

        ; Compute total error with systematic errors
        if drs.err_mode ge 3 then null_err_tot_pt = sqrt(rms_stat_null_pt ^ 2 + null_err_diam ^ 2) > sqrt(rms_sys_pt ^ 2 + null_err_diam ^ 2) else null_err_tot_pt = sqrt(rms_stat_null_pt ^ 2 + rms_sys_tot_pt ^ 2)
        tf_plot_pt_err[i_tf, i_nod] = null_err_tot_pt

        ; Print diagnostic info to screen and log
        if info gt 0 then begin
          print, ' '
          print, '==== SUB-TF COMPUTATION (pointing-based approach) ==== '
          print, '====================================================== '
          print, ' '
          print, '==== PART ' + string(i_tf + 1, format = '(I0)') + ' -- NOD ' + string(i_nod + 1, format = '(I0)') + ' ==== '
          print, '  Statistical error on TF [%] :'
          print, '    - based on individual error bars    : ', string(rms_stat_null_pt * 1d+2, format = '(F6.4)')
          print, '    - standard deviation on the mean    : ', string(disp_m_null_pt_we * 1d+2, format = '(F6.4)') + ' (' + string(disp_m_null_pt * 1d+2, format = '(F6.4)') + ' if unweighted)'
          ; PRINT, '    - quadratic sum of the 2            : ', STRING(null_err_tot_pt*1D+2, FORMAT='(F6.4)')
          print, '  Systematic error on TF [%]            : '
          print, '    - due to diameter uncertainty       : ', string(null_err_diam * 1d+2, format = '(F6.4)')
          print, '    - sqrt(excess variance)             : ', string(rms_sys_pt * 1d+2, format = '(F6.4)')
          print, '  Total error on TF [%]                 : ', string(null_err_tot_pt * 1d+2, format = '(F6.4)')
          print, '  Number of CAL pointings               : ', string(n_cal_pt_tf, format = '(I0)')
          print, ' '
          if lun gt 0 then begin
            printf, lun, ' '
            printf, lun, '==== SUB-TF COMPUTATION (pointing-based approach) ==== '
            printf, lun, '====================================================== '
            printf, lun, ' '
            printf, lun, '==== PART ' + string(i_tf + 1, format = '(I0)') + ' -- NOD ' + string(i_nod + 1, format = '(I0)') + ' ==== '
            printf, lun, '  Statistical error on TF [%] :'
            printf, lun, '    - based on individual error bars    : ', string(rms_stat_null_pt * 1d+2, format = '(F6.4)')
            printf, lun, '    - standard deviation on the mean    : ', string(disp_m_null_pt_we * 1d+2, format = '(F6.4)') + ' (' + string(disp_m_null_pt * 1d+2, format = '(F6.4)') + ' if unweighted)'
            ; PRINTF, lun, '    - quadratic sum of the 2            : ', STRING(null_err_tot_pt*1D+2, FORMAT='(F6.4)')
            printf, lun, '  Systematic error on TF [%]            : '
            printf, lun, '    - due to diameter uncertainty       : ', string(null_err_diam * 1d+2, format = '(F6.4)')
            printf, lun, '    - sqrt(excess variance)             : ', string(rms_sys_pt * 1d+2, format = '(F6.4)')
            printf, lun, '  Total error on TF [%]                 : ', string(null_err_tot_pt * 1d+2, format = '(F6.4)')
            printf, lun, '  Number of CAL pointings               : ', string(n_cal_pt_tf, format = '(I0)')
            printf, lun, ' '
          endif
        endif

        ; Ok, now save the calibrated null and null floor in output arrays if CAL_MODE is 1
        if drs.cal_mode eq 0 then begin
          idx_tmp = where(strmatch(flag_pt_tf, 'SCI') eq 1)
          for i_s = 0, n_sci_pt_tf - 1 do begin
            idx_s = where(pt_uniq_tf[idx_tmp[i_s]] eq sci_pt, n_ok)
            if n_ok le 0 then message, 'Problem with science pointing association!!'
            tf_sci[idx_s, i_nod] = tf_data[idx_tmp[i_s]]
            tf_estat[idx_s, i_nod] = rms_stat_null_pt
            tf_esyst[idx_s, i_nod] = rms_sys_tot_pt
            tf_etot[idx_s, i_nod] = null_err_tot_pt
            null_sci[idx_s, i_nod] = null_c_pt[idx_tmp[i_s]]
            null_sci_estat[idx_s, i_nod] = sqrt(null_err_pt_tf[idx_tmp[i_s], i_nod] ^ 2 + rms_stat_null_pt ^ 2)
            null_sci_esyst[idx_s, i_nod] = rms_sys_tot_pt
            null_sci_etot[idx_s, i_nod] = sqrt(null_sci_estat[idx_s] ^ 2 + null_sci_esyst[idx_s] ^ 2)
          endfor
        endif
      endfor
    endif
  endfor

  ; PREPARE OUTPUT DATA (ALSO USED FOR PLOTS)
  ; *****************************************
  if drs.cal_mode eq 1 then begin
    idx_ok = where(null_sci ne 0, n_scidat) ; 0 can happen when some SCI OB are not associated to any CAL OBs (e.g., if taken too far away from CAL OBs)
    sci_pid = long(pid[idx_sci[idx_ok]])
    sci_name = objname[idx_sci[idx_ok]]
    sci_time = double(scidata[idx_ok].mjd_obs)
    sci_utc = utc[idx_sci[idx_ok]]
    sci_lst = lst[idx_sci[idx_ok]]
    sci_ra = ra[idx_sci[idx_ok]]
    sci_dec = dec[idx_sci[idx_ok]]
    sci_alt = float(alt[idx_sci[idx_ok]])
    sci_az = float(az[idx_sci[idx_ok]])
    sci_para = float(para[idx_sci[idx_ok]])
    sci_ha = float(ha[idx_sci[idx_ok]])
    sci_ucoord = float(u_coord[idx_sci[idx_ok]])
    sci_vcoord = float(v_coord[idx_sci[idx_ok]])
    sci_wav = float(wav[idx_sci[idx_ok]])
    sci_bdw = float(bdw[idx_sci[idx_ok]])
    sci_aper = fix(aper[idx_sci[idx_ok]])
    sci_birad = fix(birad[idx_sci[idx_ok]])
    sci_borad = fix(borad[idx_sci[idx_ok]])
    sci_dit = float(dit[idx_sci[idx_ok]])
    sci_fro = long(nfr_ob[idx_sci[idx_ok]])
    sci_frr = long(nfr_rej[idx_sci[idx_ok]])
    sci_seeing = float(seeing[idx_sci[idx_ok]])
    sci_smttau = float(smttau[idx_sci[idx_ok]])
    sci_fpcp = float(fpcp[idx_sci[idx_ok]])
    qua_flag = 2 + intarr(n_scidat) ; Quality flag (manual for now)
    null_meas = double(null[idx_sci[idx_ok]])
    null_meas_err = double(null_err[idx_sci[idx_ok]])
    null_flo = double(tf_sci[idx_ok])
    null_flo_sta = double(tf_estat[idx_ok])
    null_flo_sys = double(tf_esyst[idx_ok])
    null_flo_err = double(tf_etot[idx_ok])
    null_cal = double(null_sci[idx_ok])
    null_cal_sta = double(null_sci_estat[idx_ok])
    null_cal_sys = double(null_sci_esyst[idx_ok])
    null_cal_err = double(null_sci_etot[idx_ok])
  endif else begin
    sci_pid = lonarr(n_scidat_pt, n_nod)
    sci_name = strarr(n_scidat_pt, n_nod)
    sci_time = dblarr(n_scidat_pt, n_nod)
    sci_utc = strarr(n_scidat_pt, n_nod)
    sci_lst = strarr(n_scidat_pt, n_nod)
    sci_ra = strarr(n_scidat_pt, n_nod)
    sci_dec = strarr(n_scidat_pt, n_nod)
    sci_alt = fltarr(n_scidat_pt, n_nod)
    sci_az = fltarr(n_scidat_pt, n_nod)
    sci_para = fltarr(n_scidat_pt, n_nod)
    sci_ha = fltarr(n_scidat_pt, n_nod)
    sci_wav = fltarr(n_scidat_pt, n_nod)
    sci_bdw = fltarr(n_scidat_pt, n_nod)
    sci_aper = intarr(n_scidat_pt, n_nod)
    sci_birad = intarr(n_scidat_pt, n_nod)
    sci_borad = intarr(n_scidat_pt, n_nod)
    sci_dit = fltarr(n_scidat_pt, n_nod)
    sci_fro = lonarr(n_scidat_pt, n_nod)
    sci_frr = lonarr(n_scidat_pt, n_nod)
    sci_fpcp = fltarr(n_scidat_pt, n_nod)
    sci_seeing = fltarr(n_scidat_pt, n_nod)
    sci_smttau = fltarr(n_scidat_pt, n_nod)
    sci_ucoord = fltarr(n_scidat_pt, n_nod)
    sci_vcoord = fltarr(n_scidat_pt, n_nod)
    qua_flag = 2 + intarr(n_scidat_pt)
    for i_nod = 0, n_nod - 1 do begin
      for i_pt = 0, n_scidat_pt - 1 do begin
        idx_pt = where(pt_id eq pt_uniq[idx_sci_pt[i_pt]] and nod_pos eq nod_uniq[i_nod], n_inpt)
        sci_pid[i_pt, i_nod] = long(pid[idx_pt[0]])
        sci_name[i_pt, i_nod] = objname[idx_pt[0]]
        sci_time[i_pt, i_nod] = double(mean(l1data[idx_pt[0]].mjd_obs))
        sci_utc[i_pt, i_nod] = utc[idx_pt[n_inpt / 2]]
        sci_lst[i_pt, i_nod] = lst[idx_pt[n_inpt / 2]]
        sci_ra[i_pt, i_nod] = ra[idx_pt[n_inpt / 2]]
        sci_dec[i_pt, i_nod] = dec[idx_pt[n_inpt / 2]]
        sci_alt[i_pt, i_nod] = float(mean(alt[idx_pt]))
        sci_az[i_pt, i_nod] = float(mean(az[idx_pt]))
        sci_para[i_pt, i_nod] = float(mean(para[idx_pt]))
        sci_ha[i_pt, i_nod] = float(mean(ha[idx_pt]))
        sci_ucoord[i_pt, i_nod] = float(mean(u_coord[idx_pt]))
        sci_vcoord[i_pt, i_nod] = float(mean(v_coord[idx_pt]))
        sci_wav[i_pt, i_nod] = float(wav[idx_pt[0]])
        sci_bdw[i_pt, i_nod] = float(bdw[idx_pt[0]])
        sci_aper[i_pt, i_nod] = fix(aper[idx_pt[0]])
        sci_birad[i_pt, i_nod] = fix(birad[idx_pt[0]])
        sci_borad[i_pt, i_nod] = fix(borad[idx_pt[0]])
        sci_dit[i_pt, i_nod] = float(dit[idx_pt[0]])
        sci_fro[i_pt, i_nod] = long(total(nfr_ob[idx_pt]))
        sci_frr[i_pt, i_nod] = long(total(nfr_rej[idx_pt]))
        sci_seeing[i_pt, i_nod] = float(mean(seeing[idx_pt]))
        sci_smttau[i_pt, i_nod] = float(mean(smttau[idx_pt]))
        sci_fpcp[i_pt, i_nod] = float(mean(fpcp[idx_pt]))
      endfor
    endfor
    null_meas = double(reform(null_avg_pt[idx_sci_pt, *], 1, n_scidat_pt * n_nod))
    null_meas_err = double(reform(null_err_pt[idx_sci_pt, *], 1, n_scidat_pt * n_nod))
    sci_pid = reform(sci_pid, 1, n_scidat_pt * n_nod)
    sci_name = reform(sci_name, 1, n_scidat_pt * n_nod)
    sci_time = reform(sci_time, 1, n_scidat_pt * n_nod)
    sci_utc = reform(sci_utc, 1, n_scidat_pt * n_nod)
    sci_lst = reform(sci_lst, 1, n_scidat_pt * n_nod)
    sci_ra = reform(sci_ra, 1, n_scidat_pt * n_nod)
    sci_dec = reform(sci_dec, 1, n_scidat_pt * n_nod)
    sci_alt = reform(sci_alt, 1, n_scidat_pt * n_nod)
    sci_az = reform(sci_az, 1, n_scidat_pt * n_nod)
    sci_para = reform(sci_para, 1, n_scidat_pt * n_nod)
    sci_ha = reform(sci_ha, 1, n_scidat_pt * n_nod)
    sci_wav = reform(sci_wav, n_scidat_pt * n_nod)
    sci_bdw = reform(sci_bdw, 1, n_scidat_pt * n_nod)
    sci_dit = reform(sci_dit, 1, n_scidat_pt * n_nod)
    sci_fro = reform(sci_fro, 1, n_scidat_pt * n_nod)
    sci_frr = reform(sci_frr, 1, n_scidat_pt * n_nod)
    sci_seeing = reform(sci_seeing, 1, n_scidat_pt * n_nod)
    sci_smttau = reform(sci_smttau, 1, n_scidat_pt * n_nod)
    sci_fpcp = reform(sci_fpcp, 1, n_scidat_pt * n_nod)
    null_flo = double(reform(tf_sci, 1, n_scidat_pt * n_nod))
    null_flo_sta = double(reform(tf_estat, 1, n_scidat_pt * n_nod))
    null_flo_sys = double(reform(tf_esyst, 1, n_scidat_pt * n_nod))
    null_flo_err = double(reform(tf_etot, 1, n_scidat_pt * n_nod))
    null_cal = double(reform(null_sci, 1, n_scidat_pt * n_nod))
    null_cal_sta = double(reform(null_sci_estat, 1, n_scidat_pt * n_nod))
    null_cal_sys = double(reform(null_sci_esyst, 1, n_scidat_pt * n_nod))
    null_cal_err = double(reform(null_sci_etot, 1, n_scidat_pt * n_nod))
  endelse

  ; ESTIMATE TARGET FLUX BASED ON CALIBRATOR FLUX
  ; *********************************************

  if n_scidat gt 0 and n_caldat gt 1 then begin
    ; Individual photometry
    cal_adu_dx_avg = l1data[idx_cal].photdx_avg
    cal_adu_dx_rms = l1data[idx_cal].photdx_rms
    cal_adu_sx_avg = l1data[idx_cal].photsx_avg
    cal_adu_sx_rms = l1data[idx_cal].photsx_rms
    sci_adu_dx_avg = l1data[idx_sci].photdx_avg
    sci_adu_dx_rms = l1data[idx_sci].photdx_rms
    sci_adu_sx_avg = l1data[idx_sci].photsx_avg
    sci_adu_sx_rms = l1data[idx_sci].photsx_rms
    ; Combined photometry
    cal_adu_avg = 2 * (cal_adu_dx_avg + cal_adu_sx_avg)
    cal_adu_rms = sqrt(cal_adu_dx_rms ^ 2 + cal_adu_sx_rms ^ 2)
    sci_adu_avg = 2 * (sci_adu_dx_avg + sci_adu_sx_avg)
    sci_adu_rms = sqrt(sci_adu_dx_rms ^ 2 + sci_adu_sx_rms ^ 2)
    ; Compute estimated Jy to ADU conversion factor based on estimated calibrator flux
    jy2adu_avg = tgt_flx[idx_cal] / cal_adu_avg
    jy2adu_rms = tgt_flx[idx_cal] / cal_adu_avg ^ 2 * cal_adu_rms
    jy2adu_estat = 1. / sqrt(total(1. / jy2adu_rms ^ 2))
    ; Compute mean over all measurements
    AVGSDV, jy2adu_avg, jy2adu, jy2adu_err, jy2adu_err_mean, kappa = 5
    ; Number of different science objects
    sci_uniq = sci_name[uniq(sci_name, sort(sci_name))]
    n_sci = n_elements(sci_uniq)
    if info gt 0 then begin
      print, '==== PHOTOMETRY INFO ===='
      print, 'Jy to ADU conversion               : ', string(jy2adu, format = '(F8.6)')
      print, 'Science target estimated flux [Jy] : '
      for i = 0, n_sci - 1 do begin
        idx_tmp = where(sci_name eq sci_uniq[i])
        print, '    - ' + sci_uniq[i] + '         : ' + string(mean(jy2adu * sci_adu_avg[idx_tmp]), format = '(F5.1)')
      endfor
      print, ' '
    endif
  endif

  if info gt 0 and FXPAR(header, 'NUL_MODE') eq 2 then begin
    print, '==== NSC-RELATED INFO ===='
    print, 'Mean phase offset [rad]               : ', string(mean(l1data.nsc_phavg), format = '(F8.6)')
    print, 'Mean error on mean phase offset [rad] : ', string(mean(l1data.nsc_phavg_err), format = '(F8.6)')
    print, 'Phase jitter [rad]                    : ', string(mean(l1data.nsc_phrms), format = '(F8.6)')
    print, 'Mean error on phase jitter [rad]      : ', string(mean(l1data.nsc_phrms_err), format = '(F8.6)')
    print, ' '
    if lun gt 0 then begin
      printf, lun, '==== NSC-RELATED INFO ===='
      printf, lun, 'Mean phase offset for calibrator [rad]               : ', string(mean(l1data.nsc_phavg), format = '(F8.6)')
      printf, lun, 'Mean error on mean phase offset for calibrator [rad] : ', string(mean(l1data.nsc_phavg_err), format = '(F8.6)')
      printf, lun, 'Phase jitter for calibrator [rad]                    : ', string(mean(l1data.nsc_phrms), format = '(F8.6)')
      printf, lun, 'Mean error on phase jitter for calibrator [rad]      : ', string(mean(l1data.nsc_phrms_err), format = '(F8.6)')
      printf, lun, ' '
    endif
  endif

  ; SAVE L2 FITS FILE
  ; *****************

  if not keyword_set(NO_SAVE) then begin
    ; Save results to an oifits file
    outfile = l2_dir + 'UT' + date_lng + '_calib_' + id_name + '.fits'

    ; Create new FITS file with primary HDU
    FXHMAKE, hdr, /init, /extend, 0
    FXADDPAR, hdr, 'COMMENT', 'This is a FITS file of calibrated (level 2) data taken with LBTI/' + strtrim(FXPAR(header, 'INSTRUME')) + '.'
    FXADDPAR, hdr, 'COMMENT', 'The raw data has been reduced using the nodrs software version ' + string(FXPAR(header, 'DRS_VERS'), format = '(F3.1)') + '.'
    FXADDPAR, hdr, 'COMMENT', 'Contact: lbtisupp@ipac.caltech.edu'
    FXADDPAR, hdr, 'COMMENT', 'Date and DRS information', before = 'DATE_OBS'
    FXADDPAR, hdr, 'DATE_OBS', date_lng, 'Date of observation'
    FXADDPAR, hdr, 'TELESCOP', 'LBT', 'Telescope'
    FXADDPAR, hdr, 'INSTRUME', FXPAR(header, 'INSTRUME'), 'Instrument'
    FXADDPAR, hdr, 'DRS_VERS', FXPAR(header, 'DRS_VERS'), 'Version of the DRS'
    FXADDPAR, hdr, 'DRS_DATE', FXPAR(header, 'DRS_DATE'), 'DRS version date'
    FXADDPAR, hdr, 'DATE_RED', FXPAR(header, 'DATE_RED'), 'Date of data reduction'
    FXADDPAR, hdr, 'DATE_CAL', systime(), 'Date of data calibration'
    FXADDPAR, hdr, 'COMMENT', 'Image calibration parameters', after = 'DATE_CAL'
    FXADDPAR, hdr, 'BCK_MODE', FXPAR(header, 'BCK_MODE'), '0=nod pairs, 1=0, 2=adjacent nods, 3=closest frames''
    FXADDPAR, hdr, 'N_BIN', FXPAR(header, 'N_BIN'), 'Image re-binning size'
    FXADDPAR, hdr, 'NO_BPM', FXPAR(header, 'NO_BPM'), 'Bad pixel correction off if 1 '
    FXADDPAR, hdr, 'NO_DARK', FXPAR(header, 'NO_DARK'), 'Dark subtraction off if 1 '
    FXADDPAR, hdr, 'NO_FLAT', FXPAR(header, 'NO_FLAT'), 'Flat correction off if 1'
    FXADDPAR, hdr, 'NOD_FRQ', FXPAR(header, 'NOD_FRQ'), 'Mean nodding period [s]'
    FXADDPAR, hdr, 'N_FRBCK', FXPAR(header, 'N_FRBCK'), 'Number of preserved frames for background subtraction'
    FXADDPAR, hdr, 'COMMENT', 'Flux computation parameters', after = 'N_FRBCK'
    FXADDPAR, hdr, 'FLX_MODE', FXPAR(header, 'FLX_MODE'), '0=aperture phot, 1=PSF fitting'
    FXADDPAR, hdr, 'FIT_MODE', FXPAR(header, 'FIT_MODE'), '0=no centroid, 1-3=Gaussian fit, 4=Lorentzian fit, 5=Moffat fit'
    FXADDPAR, hdr, 'COMMENT', 'Null computation parameters', after = 'FIT_MODE'
    FXADDPAR, hdr, 'OB_MODE', FXPAR(header, 'OB_MODE'), '0=all, 1=nod splitted, 2=user defined'
    FXADDPAR, hdr, 'NUL_MODE', FXPAR(header, 'NUL_MODE'), '0=mode, 1=%best, 2=statistical'
    FXADDPAR, hdr, 'R_FRNULL', FXPAR(header, 'R_FRNULL'), 'Ratio of preserved frames for null computation'
    FXADDPAR, hdr, 'NULL_COR', FXPAR(header, 'NULL_COR'), 'Subtract the null estimated from high-frequency phase noise'
    FXADDPAR, hdr, 'N_BTSTRP', FXPAR(header, 'N_BTSTRP'), 'Number of bootstrap samples to compute the error bar'
    FXADDPAR, hdr, 'NSC_CUBE', nsc_cube, 'Cube size for NSC reduction (null, phase offset, phase rms)'
    FXADDPAR, hdr, 'NSC_BFAC', nbin_fac, 'Multiplier on number of histogram bins (nbins_fac*sqrt(Np))'
    FXADDPAR, hdr, 'NSC_BINS', nsc_bins, 'Bin size for NSC reduction (0: constant, 1: variable)'
    FXADDPAR, hdr, 'NSC_OMIN', nsc_omin, 'Minimum number of occurences per bin for the fit'
    FXADDPAR, hdr, 'NULL_COR', null_cor, 'Correct for high frequency phase noise to each null measurement'
    FXADDPAR, hdr, 'NULL_LIM', null_lim, 'Acceptable raw null range (before NSC)'
    FXADDPAR, hdr, 'COMMENT', 'Null calibration parameters', after = 'NULL_LIM'
    FXADDPAR, hdr, 'CAL_METH', fix(drs.cal_method), '0=linear interpolation between two nearest neighbors, 1=global polynomial interpolation'
    FXADDPAR, hdr, 'CAL_MODE', fix(drs.cal_mode), '0=calibrate per pointing, 1=calibrate per OB'
    FXADDPAR, hdr, 'CHI2_LIM', fix(drs.cal_mode), 'Limit on acceptable chi2 (only with NSC)'
    FXADDPAR, hdr, 'NULL_EST', fix(drs.null_est), '0=unweighted, 1=weighted average of OBs'
    FXADDPAR, hdr, 'POLYDEG', fix(drs.polydeg), 'Degree of the polynomial for method CAL_MODE=1'
    FXADDPAR, hdr, 'SIG_SCA', fix(drs.sig_sca), 'Number of sigmas used for OB-based scatter OB removal'
    FXADDPAR, hdr, 'SPLT_NOD', fix(drs.split_nod), '0=calibrate both nod positions at once, 1=calibrate separately'
    FXADDPAR, hdr, 'SPLT_TF', fix(drs.split_tf), '1=split the TF where dead time > time_split'
    FXADDPAR, hdr, 'TIME_SPL', drs.split_time, 'Maximum time between 2 nulls for a single TF [hours]'
    FXWRITE, outfile, hdr

    ; Create header for the main table
    n_row = n_elements(sci_wav)
    FXBHMAKE, hdr, n_row, /init, extver = 1, 'RESULTS_SUMMARY', 'results of data calibration'

    ; Init column number
    n_col = 35
    col = lindgen(n_col) + 1l

    ; Fill extension header with column names
    FXBADDCOL, 1l, hdr, long(0), 'PID', 'Project ID number'
    FXBADDCOL, 2l, hdr, 'aaaaaaaaaaaaaaaaaaaa', 'OBJNAME', 'Object name'
    FXBADDCOL, 3l, hdr, double(0), 'MJD_OBS', 'Modified Julian Date of observation'
    FXBADDCOL, 4l, hdr, '00:00:00.00', 'LBT_UTC', 'UTC from observatory'
    FXBADDCOL, 5l, hdr, '00:00:00.00', 'LBT_LST', 'LST from observatory'
    FXBADDCOL, 6l, hdr, '+00:00:00.000', 'LBT_RA', 'RA from observatory'
    FXBADDCOL, 7l, hdr, '+00:00:00.000', 'LBT_DEC', 'DEC from observatory'
    FXBADDCOL, 8l, hdr, float(0), 'LBT_ALT', 'ALT from observatory'
    FXBADDCOL, 9l, hdr, float(0), 'LBT_AZ', 'AZ from observatory'
    FXBADDCOL, 10l, hdr, float(0), 'LBT_PARA', 'Parralactic angle from obs.'
    FXBADDCOL, 11l, hdr, float(0), 'HR_ANGLE', 'Hour angle from obs.'
    FXBADDCOL, 12l, hdr, float(0), 'U_COORD', 'u coordinate', tunit = 'cycles/arcsec'
    FXBADDCOL, 13l, hdr, float(0), 'V_COORD', 'v coordinate', tunit = 'cycles/arcsec'
    FXBADDCOL, 14l, hdr, float(0), 'INT_TIME', 'Integration time', tunit = 's'
    FXBADDCOL, 15l, hdr, float(0), 'WAV_EFF', 'Central wavelength', tunit = 'm'
    FXBADDCOL, 16l, hdr, float(0), 'BANDWIDTH', 'Bandwidth', tunit = 'm'
    FXBADDCOL, 17l, hdr, fix(0), 'APER_RAD', 'Radius for aperture photometry', tunit = 'pix'
    FXBADDCOL, 18l, hdr, fix(0), 'BCK_IRAD', 'Inner radius for background computation', tunit = 'pix'
    FXBADDCOL, 19l, hdr, fix(0), 'BCK_ORAD', 'Outer radius for background computation', tunit = 'pix'
    FXBADDCOL, 20l, hdr, double(0), 'NULL_MEAS', 'Measured raw null'
    FXBADDCOL, 21l, hdr, double(0), 'NULL_MEAS_ERR', 'Error on the measured raw null'
    FXBADDCOL, 22l, hdr, double(0), 'NULL_FLO', 'Null floor at SCI position'
    FXBADDCOL, 23l, hdr, double(0), 'NULL_FLO_STA', 'Statistical error on null floor'
    FXBADDCOL, 24l, hdr, double(0), 'NULL_FLO_SYS', 'Systematic error on null floor'
    FXBADDCOL, 25l, hdr, double(0), 'NULL_FLO_ERR', 'Total error on null floor'
    FXBADDCOL, 26l, hdr, double(0), 'NULL_CAL', 'Calibrated null'
    FXBADDCOL, 27l, hdr, double(0), 'NULL_CAL_STA', 'Statistical error on the calibrated null'
    FXBADDCOL, 28l, hdr, double(0), 'NULL_CAL_SYS', 'Systematic error on the calibrated null'
    FXBADDCOL, 29l, hdr, double(0), 'NULL_CAL_ERR', 'Total error on the calibrated null'
    FXBADDCOL, 30l, hdr, long(0), 'NFR_OB', 'Number of null frames'
    FXBADDCOL, 31l, hdr, long(0), 'NFR_REJ', 'Number of rejected null frames'
    FXBADDCOL, 32l, hdr, float(0), 'SEEING', 'Median seeing', tunit = 'arcsec'
    FXBADDCOL, 33l, hdr, float(0), 'SMTTAU', 'Precipitable water vapor estimate', tunit = 'mm'
    FXBADDCOL, 34l, hdr, float(0), 'FPC_PISTS', 'FPC RMS piston', tunit = 'um'
    FXBADDCOL, 35l, hdr, fix(0), 'QUALITY_FLG', 'Quality flag'

    ; Write extension header to FITS file
    FXBCREATE, unit, outfile, hdr
    FXBWRITM, unit, col, sci_pid, sci_name, sci_time, sci_utc, sci_lst, sci_ra, sci_dec, sci_alt, sci_az, sci_para, sci_ha, sci_ucoord, sci_vcoord, sci_dit, sci_wav, sci_bdw, $
      sci_aper, sci_birad, sci_borad, null_meas, null_meas_err, null_flo, null_flo_sta, null_flo_sys, null_flo_err, null_cal, null_cal_sta, null_cal_sys, null_cal_err, sci_fro, $
      sci_frr, sci_seeing, sci_smttau, sci_fpcp, qua_flag
    FXBFINISH, unit

    ; Copy file to new version
    file_copy, outfile, l2_dir + 'UT' + date_lng + '_calib_' + id_name + '_aper' + string(aper_rad, format = '(I0)') + '_L1v' + string(file_ver, format = '(I0)') + '.fits', /overwrite
  endif

  ; Jumping point if no science data
  skip_interp:
  if n_scidat eq 0 then begin
    message, 'WARNING, no science target in the input files', /continue
    RETURN
  endif

  ; Now start plotting....lots and lots of plots
  ; ********************************************

  if keyword_set(runbias) then bias_tag = '_BIAS' else bias_tag = ''
  ; 1. Plot the uv plane
  ; --------------------

  ; charthick = 4.0
  ; charsize  = 1.3
  uv_path = pth.result_path + 'uv' + pth.sep
  if not file_test(uv_path) then file_mkdir, uv_path
  spawn, 'chmod 775 ' + uv_path

  ; Number of unique calibrators
  sci_name = objname[idx_sci]
  sci_uc = u_coord[idx_sci]
  sci_vc = v_coord[idx_sci]
  sci_uniq = sci_name[uniq(sci_name, sort(sci_name))]
  n_scitgt = n_elements(sci_uniq)

  for i_st = 0, n_scitgt - 1 do begin
    sci_tgt = sci_uniq[i_st]
    idx_tgt = where(sci_name eq sci_tgt, n_uv)
    u_tmp = sci_uc[idx_tgt]
    v_tmp = sci_vc[idx_tgt]

    PREP_PS, /bold
    device, filename = uv_path + date_lng + '_' + sci_tgt + '_uv_coord.eps', /encapsulate, /color, xsize = 15.4, ysize = 14.7, /times
    xrange = -[-10, 10]
    yrange = -xrange
    PLOT, [0, 0], [0, 0], xtitle = 'u [cycles/arcsec]', ytitle = 'v [cycles/arcsec]', title = 'uv coordinates', $ ; ',$ ;STRTRIM(sci_name[0],2), $
      xrange = xrange, yrange = yrange, xstyle = 1, ystyle = 1, /nodata
    ; PLOTS, [MIN(xrange), MAX(xrange)], [0,0], LINESTYLE=2
    ; PLOTS, [0,0], [MIN(yrange), MAX(yrange)], LINESTYLE=2

    ; Then plot the good data
    loadct, 13, /silent
    for i_uv = 0, n_uv - 1 do oplot, [-u_tmp[i_uv], u_tmp[i_uv]], [-v_tmp[i_uv], v_tmp[i_uv]], psym = 1, color = 250

    ; Overplot the orientation of the main disk (+ 90 because pa is given for the photosphere)
    ; disk_ori = (114) * !Dpi / 180D0
    ; x = xrange[0] + (-xrange[0] + xrange[1]) * DINDGEN(1D+2)/(1D+2-1)
    ; y = TAN(disk_ori) * x
    ; OPLOT, y, x, LINESTYLE=2, THICK = 4          ; PA counted from North to (y,x) and not (x,y)

    ; Draw E-N small axis at the bootm right of the figure
    ARROW, 0.90 * xrange[1], 0.90 * yrange[0], 0.65 * xrange[1], 0.90 * yrange[0], thick = 5.5, head_indent = 0.1, /data
    ARROW, 0.90 * xrange[1], 0.90 * yrange[0], 0.90 * xrange[1], 0.65 * yrange[0], thick = 5.5, head_indent = 0.1, /data
    xyouts, 0.63 * xrange[1], 0.87 * yrange[0], 'E', charthick = charthick, charsize = charsize
    xyouts, 0.85 * xrange[1], 0.64 * yrange[0], 'N', charthick = charthick, charsize = charsize
    ; XYOUTS, 1.2, -2.4, "Moerchen midplane", ORIENTATION=70
    device, /close
    END_PS
  endfor

  ; 2. Plot data and transfer function (OB based)
  ; ---------------------------------------------

  ; Plot path (create a specific path if DIR_LABEL is set)
  ; The label must not constain _APR, which is treated differently
  if strmatch(drs.dir_label, '*_APR') then label = strmid(drs.dir_label, 0, strpos(drs.dir_label, '_APR')) $
  else label = drs.dir_label
  ; If label is not empty, add label to the TF path
  if strlen(strcompress(label, /remove_all)) eq 0 then tf_path = pth.result_path + 'TF' + pth.sep + date_lng + pth.sep $
  else tf_path = pth.result_path + 'TF' + pth.sep + date_lng + pth.sep + label + pth.sep

  if not file_test(tf_path) then file_mkdir, tf_path
  spawn, 'chmod 775 ' + tf_path
  plotname = tf_path + date_lng + '_TF-OB_APER' + string(aper_rad, format = '(I0)') + bias_tag + '.eps'
  ; Init plot
  PREP_PS, /bold
  device, filen = plotname, /encaps, /color, xsize = 17.8, ysize = 13.2, /times
  if keyword_set(no_inset) then begin
    !p.multi = [0, 1, 1]
    xtitle = 'UT hour'
  endif else begin
    !p.multi = [0, 1, 2]
    xtickname = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']
    xtitle = ' '
    postop = [0.151, 0.28, 0.954, 0.932]
    posbot1 = [0.151, 0.10, 0.954, 0.28]
  endelse
  loadct, 12, /silent
  fac = 24d
  min_time = 0.98 * min(time)
  max_time = 1.02 * max(time)
  null_plot = [-1, 5]
  xrange = [min_time, max_time] * fac
  null_range = max(null) - min(null)
  yrange = [min(null) - 0.10 * null_range, max(null) + 0.10 * null_range] * 1d2
  yrange[0] = yrange[0] < null_plot[0]
  yrange[1] = yrange[1] > null_plot[1]
  PLOT, [0], [0], xtitle = xtitle, ytitle = 'Measured null depth per OB [%]', title = title_day, xstyle = 1, ystyle = 1, $
    position = postop, xtickname = xtickname, xrange = xrange, yrange = yrange
  ; Overplot the estimated null values for calibrators measurements
  if n_caldat gt 0 then begin
    usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.7, thick = 1.5
    oplot, caltime * fac, 100. * tf_wb[idx_cal], psym = 8, color = 100
    errplot, caltime * fac, 100. * (tf_wb[idx_cal] - tf_wb_etot[idx_cal]), 100. * (tf_wb[idx_cal] + tf_wb_etot[idx_cal]), color = 100
    cal_null = null[idx_cal]
    cal_nerr = null_err[idx_cal]
    usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.7, thick = 1.5, /fill
    oplot, caltime * fac, 100. * cal_null, psym = 8, color = 100
    errplot, caltime * fac, 100. * (cal_null - cal_nerr), 100. * (cal_null + cal_nerr), color = 100
  endif
  ; Overplot the estimated null values for science measurements
  if n_scidat ge 1 then begin
    sci_null = null[idx_sci]
    sci_nerr = null_err[idx_sci]
    usersym, [1, 0, -1, 0, 1], [0, 1, 0, -1, 0], thick = 1.5, /fill
    oplot, scitime * fac, 100. * sci_null, psym = 8, color = 200
    errplot, scitime * fac, 100. * (sci_null - sci_nerr), 100. * (sci_null + sci_nerr), color = 200
  endif

  ; Overplot the estimated T2 values at science measurements (and confidence intervals)
  loadct, 0, /silent
  if drs.cal_method eq 1 and n_scidat ne 0 then begin
    for i_tf = 0, n_tf - 1 do begin
      oplot, reform(tf_time[i_tf, *]) * fac, 100. * reform(tf_plot[i_tf, *]), linestyle = 0, thick = 4
      oplot, reform(tf_time[i_tf, *]) * fac, 100. * (reform(tf_plot[i_tf, *]) - tf_plot_err[i_tf]), linestyle = 1, thick = 3.5
      oplot, reform(tf_time[i_tf, *]) * fac, 100. * (reform(tf_plot[i_tf, *]) + tf_plot_err[i_tf]), linestyle = 1, thick = 3.5
    endfor
  endif

  ; Plot background bias in the bottom inset if requested
  no_inset = 1
  if not keyword_set(no_inset) then begin
    yrange = [-0.004, 0.004] * 1d2 + avg_bias * 1d2
    PLOT, [0], [0], xtitle = 'UT hour', ytitle = 'Bckg null [%]', title = ' ', xstyle = 1, ystyle = 1, $
      position = posbot1, xrange = xrange, yrange = yrange, yminor = 3
    loadct, 12, /silent
    ; Overplot the estimated null values for calibrators measurements
    if n_caldat gt 0 then begin
      usersym, [1, 1, -1, -1, 1] * 0.5, [1, -1, -1, 1, 1] * 0.5, thick = 1.5, /fill
      oplot, caltime * fac, 1d2 * bias[idx_cal], psym = 8, color = 100
      errplot, caltime * fac, 1d2 * (bias[idx_cal] - bias_err[idx_cal]), 1d2 * (bias[idx_cal] + bias_err[idx_cal]), color = 100
    endif
    if n_scidat gt 0 then begin
      usersym, [1, 0, -1, 0, 1] * 0.5, [0, 1, 0, -1, 0] * 0.5, thick = 1.5, /fill
      oplot, scitime * fac, 1d2 * bias[idx_sci], psym = 8, color = 200
      errplot, scitime * fac, 1d2 * (bias[idx_sci] - bias_err[idx_sci]), 1d2 * (bias[idx_sci] + bias_err[idx_sci]), color = 200
    endif
    AVGSDV, bias, avg_bias_we, rms_bias, rms_bias_mean, kappa = 5, weight = 1. / bias_err ^ 2
    oplot, [0, 24], [avg_bias_we, avg_bias_we] * 1d+2, linestyle = 0, thick = 4
    oplot, [0, 24], [avg_bias_we + rms_bias_mean, avg_bias_we + rms_bias_mean] * 1d2, linestyle = 1, thick = 4
    oplot, [0, 24], [avg_bias_we - rms_bias_mean, avg_bias_we - rms_bias_mean] * 1d2, linestyle = 1, thick = 4
  endif

  ; Finish plot
  device, /close
  END_PS
  !p.multi = [0, 1, 1]

  ; 2bis. Plot filtered data and transfer function (OB based)
  ; --------------------------------------------------------

  ; Filter data
  idx_ok = where(ob_out ne 1)
  time_f = time[idx_ok]
  null_f = null[idx_ok]
  null_err_f = null_err[idx_ok]
  tf_wb_f = tf_wb[idx_ok]
  tf_wb_etot_f = tf_wb_etot[idx_ok]
  flag_f = flag[idx_ok]
  idx_cal_f = where(strmatch(flag_f, 'CAL') eq 1, n_caldat_f)
  idx_sci_f = where(strmatch(flag_f, 'SCI') eq 1, n_scidat_f)
  caltime_f = time_f[idx_cal_f]
  scitime_f = time_f[idx_sci_f]

  ; Plot path
  plotname = tf_path + date_lng + '_TF-OB_APER' + string(aper_rad, format = '(I0)') + bias_tag + '_FILT.eps'
  ; Init plot
  PREP_PS, /bold
  device, filen = plotname, /encaps, /color, xsize = 17.8, ysize = 13.2, /times
  if keyword_set(no_inset) then begin
    !p.multi = [0, 1, 1]
    xtitle = 'UT hour'
  endif else begin
    !p.multi = [0, 1, 2]
    xtickname = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']
    xtitle = ' '
    postop = [0.151, 0.28, 0.954, 0.932]
    posbot1 = [0.151, 0.10, 0.954, 0.28]
  endelse
  loadct, 12, /silent
  PLOT, [0], [0], xtitle = xtitle, ytitle = 'Measured null depth per OB [%]', title = title_day, xstyle = 1, ystyle = 1, $
    position = postop, xtickname = xtickname, xrange = xrange, yrange = yrange
  ; Overplot the estimated null values for calibrators measurements
  if n_caldat_f gt 0 then begin
    usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.7, thick = 1.5
    oplot, caltime_f * fac, 100. * tf_wb_f[idx_cal_f], psym = 8, color = 100
    errplot, caltime_f * fac, 100. * (tf_wb_f[idx_cal_f] - tf_wb_etot_f[idx_cal_f]), 100. * (tf_wb_f[idx_cal_f] + tf_wb_etot_f[idx_cal_f]), color = 100
    cal_null = null_f[idx_cal_f]
    cal_nerr = null_err_f[idx_cal_f]
    usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.7, thick = 1.5, /fill
    oplot, caltime_f * fac, 100. * cal_null, psym = 8, color = 100
    errplot, caltime_f * fac, 100. * (cal_null - cal_nerr), 100. * (cal_null + cal_nerr), color = 100
  endif
  ; Overplot the estimated null values for science measurements
  if n_scidat_f ge 1 then begin
    sci_null = null_f[idx_sci_f]
    sci_nerr = null_err_f[idx_sci_f]
    usersym, [1, 0, -1, 0, 1], [0, 1, 0, -1, 0], thick = 1.5, /fill
    oplot, scitime_f * fac, 100. * sci_null, psym = 8, color = 200
    errplot, scitime_f * fac, 100. * (sci_null - sci_nerr), 100. * (sci_null + sci_nerr), color = 200
  endif

  ; Overplot the estimated T2 values at science measurements (and confidence intervals)
  loadct, 0, /silent
  if drs.cal_method eq 1 and n_scidat ne 0 then begin
    for i_tf = 0, n_tf - 1 do begin
      oplot, reform(tf_time[i_tf, *]) * fac, 100. * reform(tf_plot[i_tf, *]), linestyle = 0, thick = 4
      oplot, reform(tf_time[i_tf, *]) * fac, 100. * (reform(tf_plot[i_tf, *]) - tf_plot_err[i_tf]), linestyle = 1, thick = 3.5
      oplot, reform(tf_time[i_tf, *]) * fac, 100. * (reform(tf_plot[i_tf, *]) + tf_plot_err[i_tf]), linestyle = 1, thick = 3.5
    endfor
  endif

  ; Plot background bias in the bottom inset if requested
  no_inset = 1 ; Depreciated
  if not keyword_set(no_inset) then begin
    yrange = [-0.004, 0.004] * 1d2 + avg_bias * 1d2
    PLOT, [0], [0], xtitle = 'UT hour', ytitle = 'Bckg null [%]', title = ' ', xstyle = 1, ystyle = 1, $
      position = posbot1, xrange = xrange, yrange = yrange, yminor = 3
    loadct, 12, /silent
    ; Overplot the estimated null values for calibrators measurements
    if n_caldat gt 0 then begin
      usersym, [1, 1, -1, -1, 1] * 0.5, [1, -1, -1, 1, 1] * 0.5, thick = 1.5, /fill
      oplot, caltime_f * fac, 1d2 * bias_f[idx_cal_f], psym = 8, color = 100
      errplot, caltime_f * fac, 1d2 * (bias_f[idx_cal_f] - bias_err_f[idx_cal_f]), 1d2 * (bias_f[idx_cal_f] + bias_err_f[idx_cal_f]), color = 100
    endif
    if n_scidat gt 0 then begin
      usersym, [1, 0, -1, 0, 1] * 0.5, [0, 1, 0, -1, 0] * 0.5, thick = 1.5, /fill
      oplot, scitime_f * fac, 1d2 * bias[idx_sci_f], psym = 8, color = 200
      errplot, scitime_f * fac, 1d2 * (bias[idx_sci_f] - bias_err_f[idx_sci_f]), 1d2 * (bias_f[idx_sci_f] + bias_err_f[idx_sci_f]), color = 200
    endif
    AVGSDV, bias, avg_bias_we, rms_bias, rms_bias_mean, kappa = 5, weight = 1. / bias_err ^ 2
    oplot, [0, 24], [avg_bias_we, avg_bias_we] * 1d+2, linestyle = 0, thick = 4
    oplot, [0, 24], [avg_bias_we + rms_bias_mean, avg_bias_we + rms_bias_mean] * 1d2, linestyle = 1, thick = 4
    oplot, [0, 24], [avg_bias_we - rms_bias_mean, avg_bias_we - rms_bias_mean] * 1d2, linestyle = 1, thick = 4
  endif

  ; Finish plot
  device, /close
  END_PS
  !p.multi = [0, 1, 1]

  ; 3. Plot data and transfer function (POINTING based)
  ; ---------------------------------------------------

  plotname = tf_path + date_lng + '_TF-PTG_APER' + string(aper_rad, format = '(I0)') + bias_tag + '.eps'
  PREP_PS, /bold
  loadct, 12, /silent
  device, filen = plotname, /encaps, /color, xsize = 17.8, ysize = 13.2, /times
  xrange = [min_time, max_time] * fac
  null_range = max(null_avg_pt) - min(null_avg_pt)
  yrange = [min(null_avg_pt) - 0.75 * null_range, max(null_avg_pt) + 0.75 * null_range] * 1d2
  yrange[0] = yrange[0] < null_plot[0]
  yrange[1] = yrange[1] > null_plot[1]
  PLOT, [0], [0], xtitle = 'UT hour', ytitle = 'Measured null depth per pointing [%]', title = title_day, xstyle = 1, ystyle = 1, $
    xrange = xrange, yrange = yrange
  ; Overplot the estimated null values for calibrators measurements
  for i_nod = 0, n_nod - 1 do begin
    if n_caldat ge 1 then begin
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, caltime_pt[*, i_nod] * fac, null_avg_pt[idx_cal_pt, i_nod] * 1d2, psym = 8, color = 100
      errplot, caltime_pt[*, i_nod] * fac, (null_avg_pt[idx_cal_pt, i_nod] + null_err_pt[idx_cal_pt, i_nod]) * 1d2, (null_avg_pt[idx_cal_pt, i_nod] - null_err_pt[idx_cal_pt, i_nod]) * 1d2, color = 100
    endif
    if n_scidat ge 1 then begin
      usersym, [1, 0, -1, 0, 1], [0, 1, 0, -1, 0], thick = 1.5, /fill
      oplot, scitime_pt[*, i_nod] * fac, null_avg_pt[idx_sci_pt, i_nod] * 1d2, psym = 8, color = 200
      errplot, scitime_pt[*, i_nod] * fac, (null_avg_pt[idx_sci_pt, i_nod] + null_err_pt[idx_sci_pt, i_nod]) * 1d2, (null_avg_pt[idx_sci_pt, i_nod] - null_err_pt[idx_sci_pt, i_nod]) * 1d2, color = 200
    endif
  endfor
  ; Overplot the estimated T2 values at science measurements (and confidence intervals)
  loadct, 0, /silent
  for i_nod = 0, n_nod - 1 do begin
    if drs.cal_method eq 1 then begin
      for i_tf = 0, n_tf - 1 do begin
        oplot, reform(tf_time[i_tf, *]) * fac, 100. * reform(tf_plot_pt[i_tf, *, i_nod]), linestyle = 0, thick = 4
        oplot, reform(tf_time[i_tf, *]) * fac, 100. * (reform(tf_plot_pt[i_tf, *, i_nod]) - tf_plot_pt_err[i_tf, i_nod]), linestyle = 1, thick = 3.5
        oplot, reform(tf_time[i_tf, *]) * fac, 100. * (reform(tf_plot_pt[i_tf, *, i_nod]) + tf_plot_pt_err[i_tf, i_nod]), linestyle = 1, thick = 3.5
      endfor
    endif
  endfor
  device, /close
  END_PS

  ; Now only if not RUNBIAS
  if not keyword_set(runbias) then begin
    ; 4. Plot the photometries and throughput (only computed if both SCI and CAL data)
    ; --------------------------------------------------------------------------------

    if n_scidat gt 0 and n_caldat gt 1 then begin
      plotname = tf_path + date_lng + '_PHOT_APER' + string(aper_rad, format = '(I0)') + '.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 17.8, ysize = 14.7, /times
      xrange = [min_time, max_time] * fac
      yrange = [0., 1.5 * max(cal_adu_sx_avg > cal_adu_dx_avg > sci_adu_sx_avg > sci_adu_dx_avg)]
      PLOT, [0], [0], xtitle = 'UT time', ytitle = 'Measured photometry [ADU]', title = title_day, xstyle = 1, ystyle = 1, $
        xrange = xrange, yrange = yrange
      ; Overplot the estimated null values for calibrators measurements
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, caltime * fac, cal_adu_dx_avg, psym = 8, color = 100
      errplot, caltime * fac, cal_adu_dx_avg + cal_adu_dx_rms, cal_adu_dx_avg - cal_adu_dx_rms, color = 100
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5
      oplot, caltime * fac, cal_adu_sx_avg, psym = 8, color = 100
      errplot, caltime * fac, cal_adu_sx_avg + cal_adu_sx_rms, cal_adu_sx_avg - cal_adu_sx_rms, color = 100
      usersym, [1, 0, -1, 0, 1], [0, 1, 0, -1, 0], thick = 1.5, /fill
      oplot, scitime * fac, [sci_adu_dx_avg], psym = 8, color = 200
      errplot, scitime * fac, [sci_adu_dx_avg + sci_adu_dx_rms], [sci_adu_dx_avg - sci_adu_dx_rms], color = 200
      usersym, [1, 0, -1, 0, 1], [0, 1, 0, -1, 0], thick = 1.5
      oplot, scitime * fac, [sci_adu_sx_avg], psym = 8, color = 200
      errplot, scitime * fac, [sci_adu_sx_avg + sci_adu_sx_rms], [sci_adu_sx_avg - sci_adu_sx_rms], color = 200
      device, /close
      END_PS

      plotname = tf_path + date_lng + '_TRANS_APER' + string(aper_rad, format = '(I0)') + '.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 17.8, ysize = 14.7, /times
      xrange = [min_time, max_time] * fac
      yrange = [0., 1.5 * max(jy2adu_avg)]
      PLOT, [0], [0], xtitle = 'UT time', ytitle = 'Jy to ADU', title = title_day, xstyle = 1, ystyle = 1, $
        xrange = xrange, yrange = yrange, xthick = xthick, ythick = ythick, charthick = charthick, charsize = charsize
      ; Overplot the estimated null values for calibrators measurements
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, caltime * fac, jy2adu_avg, psym = 8, color = 100
      errplot, caltime * fac, jy2adu_avg + jy2adu_rms, jy2adu_avg - jy2adu_rms, color = 100
      device, /close
      END_PS
    endif

    ; 5. Plot NSC results
    ; -------------------

    if FXPAR(header, 'NUL_MODE') eq 2 then begin
      ; a. Plot astrophysical nulls vs chi2
      chi2 = l1data.nsc_chi2
      chi2_range = max(chi2) - min(chi2)
      plotname = tf_path + date_lng + '_NAS-VS-CHI2_APER' + string(aper_rad, format = '(I0)') + '.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 17.8, ysize = 14.7, /times
      xrange = [(min(chi2) - 0.05 * chi2_range) > 0.1, max(chi2) + 0.05 * chi2_range]
      null_range = max(null) - min(null)
      yrange = [min(null) - 0.75 * null_range, max(null) + 0.75 * null_range] * 1d2
      PLOT, [0], [0], xtitle = '!9c!3!e2!ir!3', ytitle = 'Best-fit null depth [%]', title = title_day, xstyle = 1, ystyle = 1, $
        xrange = xrange, yrange = yrange, /xlog
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        cal_chi2 = l1data[idx_cal].nsc_chi2
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, cal_chi2, 1d2 * cal_null, psym = 8, color = 100
        errplot, cal_chi2, 1d2 * (cal_null - cal_nerr), 1d2 * (cal_null + cal_nerr), color = 100
      endif
      if n_scidat ge 1 then begin
        sci_chi2 = l1data[idx_cal].nsc_chi2
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, cal_chi2, 1d2 * sci_null, psym = 8, color = 200
        errplot, cal_chi2, 1d2 * (sci_null - sci_nerr), 1d2 * (sci_null + sci_nerr), color = 200
      endif
      device, /close
      END_PS

      ; b. Plot astrophysical nulls vs mean phase
      phavg = l1data.nsc_phavg / !dpi * wav_uniq[0] * 1d+6 * 0.5
      ph_range = max(phavg) - min(phavg)
      plotname = tf_path + date_lng + '_NAS-VS-PHAVG_APER' + string(aper_rad, format = '(I0)') + '.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 17.8, ysize = 14.7, /times
      xrange = [min(phavg) - 0.05 * ph_range, max(phavg) + 0.05 * ph_range]
      yrange = [min(null) - 0.75 * null_range, max(null) + 0.75 * null_range] * 1d2
      PLOT, [0], [0], xtitle = 'Best-fit mean OPD offset [um]', ytitle = 'Best-fit null depth [%]', title = title_day, xstyle = 1, ystyle = 1, $
        xrange = xrange, yrange = yrange
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        cal_phavg = [l1data[idx_cal].nsc_phavg / !dpi * wav_uniq[0] * 1d+6 * 0.5]
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, cal_phavg, 1d2 * cal_null, psym = 8, color = 100
        errplot, cal_phavg, 1d2 * (cal_null - cal_nerr), 1d2 * (cal_null + cal_nerr), color = 100
      endif
      if n_scidat ge 1 then begin
        sci_phavg = [l1data[idx_sci].nsc_phavg / !dpi * wav_uniq[0] * 1d+6 * 0.5]
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, sci_phavg, 1d2 * sci_null, psym = 8, color = 200
        errplot, sci_phavg, 1d2 * (sci_null - sci_nerr), 1d2 * (sci_null + sci_nerr), color = 200
      endif
      device, /close
      END_PS

      ; c. Plot astrophysical nulls vs phase RMS
      phrms = l1data.nsc_phrms / !dpi * wav_uniq[0] * 1d+6 * 0.5
      ph2_range = max(phrms) - min(phrms)
      plotname = tf_path + date_lng + '_NAS-VS-PHRMS_APER' + string(aper_rad, format = '(I0)') + '.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 17.8, ysize = 14.7, /times
      xrange = [min(phrms) - 0.05 * ph_range, max(phrms) + 0.05 * ph_range]
      yrange = [min(null) - 0.75 * null_range, max(null) + 0.75 * null_range] * 1d2
      PLOT, [0], [0], xtitle = 'Best-fit OPD jitter [um]', ytitle = 'Best-fit null depth [%]', title = title_day, xstyle = 1, ystyle = 1, $
        xrange = xrange, yrange = yrange
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        cal_phrms = [l1data[idx_cal].nsc_phrms / !dpi * wav_uniq[0] * 1d+6 * 0.5]
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, cal_phrms, 1d2 * cal_null, psym = 8, color = 100
        errplot, cal_phrms, 1d2 * (cal_null - cal_nerr), 1d2 * (cal_null + cal_nerr), color = 100
      endif
      if n_scidat ge 1 then begin
        sci_phrms = [l1data[idx_sci].nsc_phrms / !dpi * wav_uniq[0] * 1d+6 * 0.5]
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, sci_phrms, 1d2 * sci_null, psym = 8, color = 200
        errplot, sci_phrms, 1d2 * (sci_null - sci_nerr), 1d2 * (sci_null + sci_nerr), color = 200
      endif
      device, /close
      END_PS

      ; d. Plot mean phase vs time
      ph_range = max(phavg) - min(phavg)
      plotname = tf_path + date_lng + '_PHAVG_APER' + string(aper_rad, format = '(I0)') + '.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
      xrange = [min_time, max_time] * fac
      yrange = [0., max(phavg) + 0.05 * ph_range]
      PLOT, [0], [0], xtitle = 'UT time', ytitle = 'Best-fit mean OPD offset [um]', title = title_day, xstyle = 1, ystyle = 9, $
        xrange = xrange, yrange = yrange, position = [0.12, 0.12, 0.88, 0.90]
      ; Add right vertical axis
      AXIS, ytitle = 'Best-fit mean phase offset [rad]', /yaxis, ystyle = 1, yrange = yrange * !dpi / (wav_uniq[0] * 1d+6 * 0.5)
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        cal_phavg_err = [l1data[idx_cal].nsc_phavg_err / !dpi * wav_uniq[0] * 1d+6 * 0.5]
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, caltime * fac, cal_phavg, psym = 8, color = 100
        errplot, caltime * fac, cal_phavg - cal_phavg_err, cal_phavg + cal_phavg_err, color = 100
      endif
      if n_scidat ge 1 then begin
        sci_phavg_err = [l1data[idx_sci].nsc_phavg_err / !dpi * wav_uniq[0] * 1d+6 * 0.5]
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, scitime * fac, sci_phavg, psym = 8, color = 200
        errplot, scitime * fac, sci_phavg - sci_phavg_err, sci_phavg + sci_phavg_err, color = 200
      endif
      device, /close
      END_PS

      ; e. Plot rms phase vs time
      phrms_err = l1data.nsc_phavg_err / !dpi * wav_uniq[0] * 1d+6 * 0.5
      ph_range = max(phrms + phrms_err) - min(phrms - phrms_err)
      plotname = tf_path + date_lng + '_PHRMS_APER' + string(aper_rad, format = '(I0)') + '.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
      xrange = [min_time, max_time] * fac
      yrange = [0., max(phrms + phrms_err) + 0.05 * ph_range]
      PLOT, [0], [0], xtitle = 'UT time', ytitle = 'Best-fit OPD jitter [!9m!3m]', title = title_day, xstyle = 1, ystyle = 9, $
        xrange = xrange, yrange = yrange, position = [0.12, 0.12, 0.88, 0.90]
      ; Add right vertical axis
      AXIS, ytitle = 'Best-fit phase jitter [rad]', /yaxis, ystyle = 1, yrange = yrange * !dpi / (wav_uniq[0] * 1d+6 * 0.5)
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        cal_phrms_err = [l1data[idx_cal].nsc_phavg_err / !dpi * wav_uniq[0] * 1d+6 * 0.5]
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, caltime * fac, cal_phrms, psym = 8, color = 100
        errplot, caltime * fac, cal_phrms - cal_phrms_err, cal_phrms + cal_phrms_err, color = 100
      endif
      if n_scidat ge 1 then begin
        sci_phrms_err = [l1data[idx_sci].nsc_phavg_err / !dpi * wav_uniq[0] * 1d+6 * 0.5]
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, scitime * fac, sci_phrms, psym = 8, color = 200
        errplot, scitime * fac, sci_phrms - sci_phrms_err, sci_phrms + sci_phrms_err, color = 200
      endif
      device, /close
      END_PS

      ; f. Plot chi2  vs time
      plotname = tf_path + date_lng + '_CHI2_APER' + string(aper_rad, format = '(I0)') + '.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
      yrange = [(min(cal_chi2 < sci_chi2) - 0.05 * chi2_range) > 0.1, max(cal_chi2 > sci_chi2) + 0.05 * chi2_range]
      PLOT, [0], [0], xtitle = 'UT time', ytitle = '!9c!3!e2!ir!3', title = title_day, xstyle = 1, ystyle = 1, $
        xrange = xrange, yrange = yrange
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, caltime * fac, cal_chi2, psym = 8, color = 100
      endif
      if n_scidat ge 1 then begin
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, scitime * fac, sci_chi2, psym = 8, color = 200
      endif
      chi2_mean = mean([cal_chi2, sci_chi2])
      oplot, [-24, 24], [chi2_mean, chi2_mean], linestyle = 1, thick = 3
      device, /close
      END_PS

      ; g. Plot nas error vs sigma phi
      plotname = tf_path + date_lng + '_NASERR-VS-PHIDPHI_APER' + string(aper_rad, format = '(I0)') + '.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
      PLOT, [0], [0], xtitle = '(Best-fit phase jitter !9s!if!3!n)/(Best-fit mean phase !9m!if!3!n)', ytitle = 'Error on best-fit null depth [%]', title = title_day, xstyle = 1, ystyle = 1, $
        xrange = [0., 5.0], yrange = [0., 1d2 * (max(sci_nerr) > max(cal_nerr)) + 0.1]
      ; Overplot the estimated null values for calibrators measurements
      if n_scidat ge 1 then begin
        phidphi_c = (cal_phrms) / cal_phavg ; cal_phrms*SQRT(cal_phavg)
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, phidphi_c, 1d2 * cal_nerr, psym = 8, color = 100
      endif
      if n_scidat ge 1 then begin
        phidphi_s = (sci_phrms) / sci_phavg ; sci_phrms*SQRT(sci_phavg)
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, phidphi_s, 1d2 * sci_nerr, psym = 8, color = 200
      endif
      device, /close
      END_PS

      ; h. Plot measured smttau vs phase
      plotname = tf_path + date_lng + '_PHRMS-vs_PWV.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
      xrange = [0., 5.0]
      yrange = [0., max(phrms) + 0.5 * ph_range]
      PLOT, [0], [0], xtitle = 'PWV [mm]', ytitle = 'Measured OPD jitter [um]', title = title_day, xstyle = 1, ystyle = 9, $
        xrange = xrange, yrange = yrange, position = [0.12, 0.12, 0.88, 0.90]
      ; Add right vertical axis
      AXIS, ytitle = 'Measured phase jitter [rad]', /yaxis, ystyle = 1, yrange = yrange * !dpi / (wav_uniq[0] * 1d+6 * 0.5)
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, smttau[idx_cal], cal_phrms, psym = 8, color = 100
        errplot, smttau[idx_cal], cal_phrms - cal_phrms_err, cal_phrms + cal_phrms_err, color = 100
      endif
      if n_scidat ge 1 then begin
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, smttau[idx_sci], sci_phrms, psym = 8, color = 200
        errplot, smttau[idx_sci], sci_phrms - sci_phrms_err, sci_phrms + sci_phrms_err, color = 200
      endif
      device, /close
      END_PS

      ; i. Plot raw null RMS vs measured error
      plotname = tf_path + date_lng + '_PHERR-vs_NULLRMS.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
      xrange = [0., max(null_rms)]
      yrange = [0., max(phrms) + 0.5 * ph_range]
      PLOT, [0], [0], xtitle = '(Raw null RMS)/(background RMS relative to peak)', ytitle = 'Measured OPD jitter [um]', title = title_day, xstyle = 1, ystyle = 9, $
        xrange = xrange, yrange = yrange, position = [0.12, 0.12, 0.88, 0.90]
      ; Add right vertical axis
      AXIS, ytitle = 'Measured phase jitter [rad]', /yaxis, ystyle = 1, yrange = yrange * !dpi / (wav_uniq[0] * 1d+6 * 0.5)
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, null_rms[idx_cal], cal_phrms, psym = 8, color = 100
        errplot, null_rms[idx_cal], cal_phrms - cal_phrms_err, cal_phrms + cal_phrms_err, color = 100
      endif
      if n_scidat ge 1 then begin
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, null_rms[idx_sci], sci_phrms, psym = 8, color = 200
        errplot, null_rms[idx_sci], sci_phrms - sci_phrms_err, sci_phrms + sci_phrms_err, color = 200
      endif
      device, /close
      END_PS

      ; j. Plot raw null RMS vs measured error
      plotname = tf_path + date_lng + '_NULLAVG-vs_NULLRMS.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
      xrange = [0., max(null_rms)]
      yrange = 2 * max(sci_null) * [-1., 1.]
      PLOT, [0], [0], xtitle = '(Raw null RMS)/(background RMS relative to peak)', ytitle = 'Measured null', title = title_day, xstyle = 1, ystyle = 1, $
        xrange = xrange, yrange = yrange, position = [0.12, 0.12, 0.88, 0.90]
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, null_rms[idx_cal], cal_null, psym = 8, color = 100
        errplot, null_rms[idx_cal], cal_null - cal_nerr, cal_null + cal_nerr, color = 100
      endif
      if n_scidat ge 1 then begin
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, null_rms[idx_sci], sci_null, psym = 8, color = 200
        errplot, null_rms[idx_sci], sci_null - sci_nerr, sci_null + sci_nerr, color = 200
      endif
      device, /close
      END_PS

      ; k. Plot measured error vs error expected from background noise
      plotname = tf_path + date_lng + '_NULLERR-vs_BCKGERR.eps'
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
      xrange = [0., max(null_err_phot)] * 1d2
      yrange = [0, 2 * max(sci_nerr)] * 1d2
      PLOT, [0], [0], xtitle = 'Photometric error on background relative to peak [%]', ytitle = 'Null errror [%]', title = title_day, xstyle = 1, ystyle = 1, $
        xrange = xrange, yrange = yrange, position = [0.12, 0.12, 0.88, 0.90]
      oplot, indgen(10), indgen(10)
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, null_err_phot[idx_cal] * 1d2, cal_nerr * 1d2, psym = 8, color = 100
      endif
      if n_scidat ge 1 then begin
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, null_err_phot[idx_sci] * 1d2, sci_nerr * 1d2, psym = 8, color = 200
      endif
      device, /close
      END_PS

      ; k. Plot measured error vs error expected from background noise
      plotname = tf_path + date_lng + '_EXCESS-PHOT-ERR-vs-time.eps'
      excess_cal = cal_nerr / null_err_phot[idx_cal]
      excess_sci = sci_nerr / null_err_phot[idx_sci]
      fac = 24d
      min_time = 0.98 * min(time)
      max_time = 1.02 * max(time)
      xrange = [min_time, max_time] * fac
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
      yrange = [0, 1.5 * max(excess_sci)]
      PLOT, [0], [0], xtitle = 'UT time', ytitle = '(Null error) / (Photometric error relative to peak)', title = title_day, xstyle = 1, ystyle = 1, $
        xrange = xrange, yrange = yrange, position = [0.12, 0.12, 0.88, 0.90]
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, caltime * fac, excess_cal, psym = 8, color = 100
      endif
      if n_scidat ge 1 then begin
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, scitime * fac, excess_sci, psym = 8, color = 200
      endif
      device, /close
      END_PS

      ; k. Plot measured error vs error expected from background noise
      plotname = tf_path + date_lng + '_EXCESS_PHASE-ERR-vs-time.eps'
      excess_cal = cal_nerr / null_err_phase[idx_cal]
      excess_sci = sci_nerr / null_err_phase[idx_sci]
      fac = 24d
      min_time = 0.98 * min(time)
      max_time = 1.02 * max(time)
      xrange = [min_time, max_time] * fac
      PREP_PS, /bold
      loadct, 12, /silent
      device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
      yrange = [0, 1.5 * max(excess_sci)]
      PLOT, [0], [0], xtitle = 'UT time', ytitle = '(Null error) / (Error from phase variations)', title = title_day, xstyle = 1, ystyle = 1, $
        xrange = xrange, yrange = yrange, position = [0.12, 0.12, 0.88, 0.90]
      ; Overplot the estimated null values for calibrators measurements
      if n_caldat ge 1 then begin
        usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
        oplot, caltime * fac, excess_cal, psym = 8, color = 100
      endif
      if n_scidat ge 1 then begin
        usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
        oplot, scitime * fac, excess_sci, psym = 8, color = 200
      endif
      device, /close
      END_PS
    endif

    ; 6. Plot fpc piston vs time
    ; --------------------------

    ph_range = max(fpcp) - min(fpcp)
    plotname = tf_path + date_lng + '_PCP.eps'
    PREP_PS, /bold
    loadct, 12, /silent
    device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
    xrange = [min_time, max_time] * fac
    yrange = [0., max(fpcp) + 0.25 * ph_range]
    PLOT, [0], [0], xtitle = 'UT time', ytitle = 'PHASECam OPD jitter [um]', title = title_day, xstyle = 1, ystyle = 9, $
      xrange = xrange, yrange = yrange, position = [0.12, 0.12, 0.88, 0.90]
    ; Add right vertical axis
    AXIS, ytitle = 'PHASECam phase jitter [rad]', /yaxis, ystyle = 1, yrange = yrange * !dpi / (wav_uniq[0] * 1d+6 * 0.5)
    ; Overplot the estimated null values for calibrators measurements
    if n_caldat ge 1 then begin
      cal_fpcp = [phasec[idx_cal].fpc_pists * plc_wav[idx_cal] / 360.]
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, caltime * fac, cal_fpcp, psym = 8, color = 100
      ; ERRPLOT, caltime*fac, cal_fpcp-cal_fpcp_err, cal_fpcp+cal_fpcp_err, COLOR=100
    endif
    if n_scidat ge 1 then begin
      sci_fpcp = [phasec[idx_sci].fpc_pists * plc_wav[idx_sci] / 360.]
      usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
      oplot, scitime * fac, sci_fpcp, psym = 8, color = 200
      ; ERRPLOT, scitime*fac, sci_fpcp-sci_fpcp_err, sci_fpcp+sci_fpcp_err, COLOR=200
    endif
    device, /close
    END_PS

    ; 7. Plot measured RMS phase (over the last DITms)  vs time
    ; --------------------------------------------------------

    ph_range = max(phstd) - min(phstd)
    plotname = tf_path + date_lng + '_PHSTD.eps'
    PREP_PS, /bold
    loadct, 12, /silent
    device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
    yrange = [0., max(phstd) + 2.0 * ph_range]
    PLOT, [0], [0], xtitle = 'UT time', ytitle = 'Measured OPD jitter [um]', title = title_day, xstyle = 1, ystyle = 9, $
      xrange = xrange, yrange = yrange, position = [0.12, 0.12, 0.88, 0.90]
    ; Add right vertical axis
    AXIS, ytitle = 'Measured phase jitter [rad]', /yaxis, ystyle = 1, yrange = yrange * !dpi / (wav_uniq[0] * 1d+6 * 0.5)
    ; Overplot the estimated null values for calibrators measurements
    if n_caldat ge 1 then begin
      cal_phstd = [phasec[idx_cal].pcphstd * plc_wav[idx_cal] / 360.] ; Convert to um
      cal_phstd_err = [phasec[idx_cal].pcphstd_err * plc_wav[idx_cal] / 360.] ; Convert to um
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, caltime * fac, cal_phstd, psym = 8, color = 100
      errplot, caltime * fac, cal_phstd - cal_phstd_err, cal_phstd + cal_phstd_err, color = 100
    endif
    if n_scidat ge 1 then begin
      sci_phstd = [phasec[idx_sci].pcphstd * plc_wav[idx_sci] / 360.] ; Convert to um
      sci_phstd_err = [phasec[idx_sci].pcphstd_err * plc_wav[idx_sci] / 360.] ; Convert to um
      usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
      oplot, scitime * fac, sci_phstd, psym = 8, color = 200
      errplot, scitime * fac, sci_phstd - sci_phstd_err, sci_phstd + sci_phstd_err, color = 200
    endif
    device, /close
    END_PS

    phavg = phasec.pcphmean * plc_wav / 360.
    ph_range = (max(phavg) - min(phavg))
    plotname = tf_path + date_lng + '_PHAVG.eps'
    PREP_PS, /bold
    loadct, 12, /silent
    device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
    yrange = [min(phavg) - 0.5 * ph_range, max(phavg) + 0.5 * ph_range]
    PLOT, [0], [0], xtitle = 'UT time', ytitle = 'Measured mean OPD [um]', title = title_day, xstyle = 1, ystyle = 9, $
      xrange = xrange, yrange = yrange, position = [0.12, 0.12, 0.88, 0.90]
    ; Add right vertical axis
    AXIS, ytitle = 'Measured mean phase [rad]', /yaxis, ystyle = 1, yrange = yrange * !dpi / (wav_uniq[0] * 1d+6 * 0.5)
    ; Overplot the estimated null values for calibrators measurements
    if n_caldat ge 1 then begin
      cal_phavg = [phasec[idx_cal].pcphmean * plc_wav[idx_cal] * 1d6 / 360.] ; Convert to um
      cal_phavg_err = [phasec[idx_cal].pcphmean_err * plc_wav[idx_cal] * 1d6 / 360.] ; Convert to um
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, caltime * fac, cal_phavg, psym = 8, color = 100
      errplot, caltime * fac, cal_phavg - cal_phavg_err, cal_phavg + cal_phavg_err, color = 100
    endif
    if n_scidat ge 1 then begin
      sci_phavg = [phasec[idx_sci].pcphmean * plc_wav[idx_sci] * 1d6 / 360.] ; Convert to um
      sci_phavg_err = [phasec[idx_sci].pcphmean_err * plc_wav[idx_sci] * 1d6 / 360.] ; Convert to um
      usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
      oplot, scitime * fac, sci_phavg, psym = 8, color = 200
      errplot, scitime * fac, sci_phavg - sci_phavg_err, sci_phavg + sci_phavg_err, color = 200
    endif
    device, /close
    END_PS

    ; 8. Plot measured seeing vs time
    ; -------------------------------

    plotname = tf_path + date_lng + '_SEEING.eps'
    PREP_PS, /bold
    loadct, 12, /silent
    device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
    yrange = [0., 2.0]
    PLOT, [0], [0], xtitle = 'UT time', ytitle = 'Seeing [arcsec]', title = title_day, xstyle = 1, ystyle = 1, $
      xrange = xrange, yrange = yrange
    ; Overplot the estimated null values for calibrators measurements
    if n_scidat ge 1 then begin
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, caltime * fac, seeing[idx_cal], psym = 8, color = 100
    endif
    if n_scidat ge 1 then begin
      usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
      oplot, scitime * fac, seeing[idx_sci], psym = 8, color = 200
    endif
    device, /close
    END_PS

    ; 9. Plot measured smttau vs time
    ; -------------------------------

    plotname = tf_path + date_lng + '_PWV.eps'
    PREP_PS, /bold
    loadct, 12, /silent
    device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
    yrange = [0., 10.0]
    PLOT, [0], [0], xtitle = 'UT time', ytitle = 'PWV [mm]', title = title_day, xstyle = 1, ystyle = 1, $
      xrange = xrange, yrange = yrange
    ; Overplot the estimated null values for calibrators measurements
    if n_scidat ge 1 then begin
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, caltime * fac, smttau[idx_cal], psym = 8, color = 100
    endif
    if n_scidat ge 1 then begin
      usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
      oplot, scitime * fac, smttau[idx_sci], psym = 8, color = 200
    endif
    device, /close
    END_PS

    ; 9. Plot measured windspeed vs time
    ; ----------------------------------

    plotname = tf_path + date_lng + '_WIND.eps'
    PREP_PS, /bold
    loadct, 12, /silent
    device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
    yrange = [0., 25.0]
    PLOT, [0], [0], xtitle = 'UT time', ytitle = 'Wind speed [m/s]', title = title_day, xstyle = 1, ystyle = 1, $
      xrange = xrange, yrange = yrange ; , POSITION=[0.12,0.12,0.88,0.90]
    ; Overplot the estimated null values for calibrators measurements
    if n_caldat ge 1 then begin
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, caltime * fac, wind[idx_cal], psym = 8, color = 100
    endif
    if n_scidat ge 1 then begin
      usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
      oplot, scitime * fac, wind[idx_sci], psym = 8, color = 200
    endif
    device, /close
    END_PS

    ; 10. Plot position on detector
    ; -----------------------------
    plotname = tf_path + date_lng + '_BEAM-POS.eps'
    PREP_PS, /bold
    loadct, 12, /silent
    device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
    xrange = [min(xcen) - 10, max(xcen) + 10]
    yrange = [min(ycen) - 10, max(ycen) + 10]
    PLOT, [0], [0], xtitle = 'Initial x position on detector [pix]', ytitle = 'Initial y position on detector [pix]', title = title_day, xstyle = 1, ystyle = 1, $
      xrange = xrange, yrange = yrange
    ; Overplot the estimated null values for calibrators measurements
    if n_caldat ge 1 then begin
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, xcen[idx_cal], ycen[idx_cal], psym = 8, color = 100
    endif
    if n_scidat ge 1 then begin
      usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
      oplot, xcen[idx_sci], ycen[idx_sci], psym = 8, color = 200
    endif
    device, /close
    END_PS

    ; 11. Plot position on detector (vs time)
    ; -----------------------------
    fac = 24d
    min_time = 0.98 * min(time)
    max_time = 1.02 * max(time)
    xrange = [min_time, max_time] * fac
    yrange = [min(xcen) - 10, max(xcen) + 10]
    plotname = tf_path + date_lng + '_BEAM-XPOS-TIME.eps'
    PREP_PS, /bold
    loadct, 12, /silent
    device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
    PLOT, [0], [0], xtitle = 'UT hour', ytitle = 'Initial x position on detector [pix]', title = title_day, xstyle = 1, ystyle = 1, $
      xrange = xrange, yrange = yrange
    ; Overplot the estimated null values for calibrators measurements
    if n_caldat ge 1 then begin
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, caltime * fac, xcen[idx_cal], psym = 8, color = 100
    endif
    if n_scidat ge 1 then begin
      usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
      oplot, scitime * fac, xcen[idx_sci], psym = 8, color = 200
    endif
    device, /close
    END_PS

    yrange = [min(ycen) - 10, max(ycen) + 10]
    plotname = tf_path + date_lng + '_BEAM-YPOS-TIME.eps'
    PREP_PS, /bold
    loadct, 12, /silent
    device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
    PLOT, [0], [0], xtitle = 'UT hour', ytitle = 'Initial x position on detector [pix]', title = title_day, xstyle = 1, ystyle = 1, $
      xrange = xrange, yrange = yrange
    ; Overplot the estimated null values for calibrators measurements
    if n_caldat ge 1 then begin
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, caltime * fac, ycen[idx_cal], psym = 8, color = 100
    endif
    if n_scidat ge 1 then begin
      usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
      oplot, scitime * fac, ycen[idx_sci], psym = 8, color = 200
    endif
    device, /close
    END_PS

    ; 12. Plot position on detector (vs null)
    ; -----------------------------

    xrange = [min(xcen) - 10, max(xcen) + 10]
    null_range = max(null) - min(null)
    yrange = [min(null) - 0.75 * null_range, max(null) + 0.75 * null_range] * 1d2
    plotname = tf_path + date_lng + '_BEAM-XPOS-NULL.eps'
    PREP_PS, /bold
    loadct, 12, /silent
    device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
    PLOT, [0], [0], xtitle = 'Initial x position on detector [pix]', ytitle = 'Null', title = title_day, xstyle = 1, ystyle = 1, $
      xrange = xrange, yrange = yrange
    ; Overplot the estimated null values for calibrators measurements
    if n_caldat ge 1 then begin
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, xcen[idx_cal], 100. * cal_null, psym = 8, color = 100
      errplot, xcen[idx_cal], 100. * (cal_null - cal_nerr), 100. * (cal_null + cal_nerr), color = 100
    endif
    if n_scidat ge 1 then begin
      usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
      oplot, xcen[idx_sci], 100. * sci_null, psym = 8, color = 200
      errplot, xcen[idx_sci], 100. * (sci_null - sci_nerr), 100. * (sci_null + sci_nerr), color = 200
    endif
    device, /close
    END_PS

    xrange = [min(ycen) - 10, max(ycen) + 10]
    plotname = tf_path + date_lng + '_BEAM-YPOS-NULL.eps'
    PREP_PS, /bold
    loadct, 12, /silent
    device, filen = plotname, /encaps, /color, xsize = 19.5, ysize = 14.7, /times
    PLOT, [0], [0], xtitle = 'Initial y position on detector [pix]', ytitle = 'Null', title = title_day, xstyle = 1, ystyle = 1, $
      xrange = xrange, yrange = yrange ; , POSITION=[0.1,0.1,0.90,0.92]
    ; Overplot the estimated null values for calibrators measurements
    if n_caldat ge 1 then begin
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.8, thick = 1.5, /fill
      oplot, ycen[idx_cal], 100. * cal_null, psym = 8, color = 100
      errplot, ycen[idx_cal], 100. * (cal_null - cal_nerr), 100. * (cal_null + cal_nerr), color = 100
    endif
    if n_scidat ge 1 then begin
      usersym, [1, 0, -1, 0, 1] * 1.0, [0, 1, 0, -1, 0] * 1.0, thick = 1.5, /fill
      oplot, ycen[idx_sci], 100. * sci_null, psym = 8, color = 200
      errplot, ycen[idx_sci], 100. * (sci_null - sci_nerr), 100. * (sci_null + sci_nerr), color = 200
    endif
    device, /close
    END_PS
  endif

  ; Close file
  if keyword_set(log_file) then begin
    printf, lun, ' '
    printf, lun, ' '
    close, lun
    free_lun, lun
  endif

  bck_irad = 0
  bck_orad = 0
  ; Write summary file
  if keyword_set(log_file) and not keyword_set(runbias) then begin
    ; Create log file
    log_file = l2_dir + 'reduction_history.txt'
    if not file_test(log_file) then begin
      openw, lun, log_file, /get_lun, width = 800, /append
      printf, lun, 'Column signification'
      printf, lun, '1: L1 file version'
      printf, lun, '2: Statistical error on TF [%]'
      printf, lun, '3: Systematic error on TF [%]'
      printf, lun, '4: Background subtraction mode'
      printf, lun, '5: Null mode (0: best x%, 1: mode, and 2: NSC)'
      printf, lun, '6: Null correction based on high-frequency phase jitter'
      printf, lun, '7: Multiplier on number of histogram bins (nbins_fac*sqrt(Np))'
      printf, lun, '8: Bin size for NSC reduction (0: constant, 1: variable).'
      printf, lun, '9: Number of bootstrap samples'
      printf, lun, '10: 3-element vector containing the number of elements in the chi2 cube along the null, mean phase and rms phase directions respectively.'
      printf, lun, '11: Minimum number of occurences per bin for the fit'
      printf, lun, '12: Range of acceptable null (before NSC, in %)'
      printf, lun, '13: Radius of the photmetric aperture (in pixels)'
      printf, lun, '14: Inner radius of background region (in pixels)'
      printf, lun, '15: Outer radius of background region (in pixels, 0 to use the closer vertical limit of the channel)'
      printf, lun, '16: L1 DRS version ID number'
      printf, lun, '17: L2 DRS version ID number'
      printf, lun, '18: Execution date'
      printf, lun, ' 1 ', '  ', '   2  ', '  ', '   3  ', '  ', ' 4', '  ', '5', '  ', '6', '  ', '7', '  ', '8', '  ', '  9 ', '  ', '    10    ', ' ', '11', '  ', '    12     ', '  ', '13', '  ', '14', '  ', '15', '  ', ' 16 ', '  ', ' 17 ', '  ', '        18     '
      printf, lun, '--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    endif else openw, lun, log_file, /get_lun, width = 800, /append
    printf, lun, string(file_ver, format = '(I03)'), '  ', string(mean(tf_estat) * 1d2, format = '(F6.4)'), '  ', string(mean(tf_esyst) * 1d+2, format = '(F6.4)'), '  ', string(bck_mode, format = '(I0)'), '  ', string(null_mod, format = '(I0)'), '  ', string(null_cor, format = '(I0)'), '  ', $
      string(nbin_fac, format = '(I0)'), '  ', string(nsc_bins, format = '(I0)'), '  ', string(n_btstrp, format = '(I04)'), '  ', nsc_cube, '  ', string(nsc_omin, format = '(I02)'), '  ', null_lim, '  ', string(aper_rad, format = '(I02)'), '  ', string(bck_irad, format = '(I02)'), '  ', string(bck_orad, format = '(I02)'), '  ', $
      string(FXPAR(header, 'DRS_VERS'), format = '(F4.2)'), '  ', string(drs_version, format = '(F4.2)'), '  ', systime()
    close, lun
    free_lun, lun
  endif
end
