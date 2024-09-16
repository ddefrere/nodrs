;+
; NAME: LBTI_DRS
;
; PURPOSE
;   Main procedure to reduce LBTI data.
;
; MANDATORY INPUTS
;   date          :  String with the UT date to be reduced based on the format 'yymmdd' (e.g., '130524')
;   cfg_file      :  String with the name of the config file with the reduction parameters
;
; OPTIONAL INPUT KEYWORDS
;   BAD_IDX       :  Vector with the file number of frames to be removed
;   BCKG_IDX      :  Two-element vector with the lower and the upper file numbers of the background files (e.g., [10,19]). Generally used when BCKG_MODE = 4.
;   DARK_IDX      :  Two-element vector with the lower and the upper file numbers of the dark files (e.g., [0,9])
;   DATA_IDX      :  Two-element vector with the lower and the upper file numbers of the data files (e.g., [20,999])
;   FLAT_IDX      :  Two-element vector with the lower and the upper file numbers of the flat files (e.g., [10,19])
;   NOD_IDX       :  2x(number of nods) array with the lower and the upper file numbers of each nod position (e.g., [[120,239],[240,299],[300,399]])
;   OB_IDX        :  Two-element vector with the lower and upper OB ID numbers to process during the null computation
;   DATA_PATH     :  String vector pointing to the the path of scientific data (superseed date)
;   BCKG_PATH     :  String vector pointing to the the path of background data (superseed bckg_idx keyword)
;   DARK_PATH     :  String vector pointing to the the path of dark data (superseed dark_idx keyword)
;   FLAT_PATH     :  String vector pointing to the the path of flat data (superseed flat_idx keyword)
;
; RUNNING INPUT KEYWORDS
;   MASTERLOG     :  Set this keyword to force the creation of new masterlog file
;   RENEW         :  Set this keyword to force the code to create a new file (e.g., darks, flats, l0_fits, l1_fits, ...).
;                 :  Otherwise, the code doesn't process the data if the appropriate intermediate file already exists (with the same reduction parameters)
;   SKIP_ADI      :  Set to skip ADI processing (superseed the value in the cnfig file)
;   SKIP_FLX      :  Set to skip flux computation (files restored from disk, superseed the value in the config file)
;   SKIP_NULL     :  Set to skip null computation (superseed the value in the cnfig file)
;   SKIP_RED      :  Set to skip image reduction (files restored from disk, superseed the value in the config file)
;   LOG_FILE      :  Set this keyword to save the results in an external log file
;   NO_MULTI      :  Set this keyword to turn off multi-threading
;   NO_SAVE       :  Set this keyword to turn off data saving
;   PLOT          :  Set this keyword to plot the data
;   VERBOSE       :  Define the level of information printed to screen:
;                       - 0: completely silent execution
;                       - 1: minimum level of information
;                       - 2: nominal level of information
;                       - 3: debugging use (print for instance the computation time of each step in the data reduction)
;
; PRE-REQUISITES
;   1. Up-to-date Astrolib (http://idlastro.gsfc.nasa.gov/)
;   2. Astrolib now requires the Coyote library :-/ (http://www.idlcoyote.com/programs/zip_files/coyoteprograms.zip )
;   3. Several Gb of free disk storage to reduce a typical HOSTS night (10000 256x256 frames uses 5 Gb)
;
; CALLING SEQUENCE
;   LBTI_DRS, '132706', 'hosts.cfg', /VERBOSE
;   or LBTI_DRS, '132706', 'hosts.cfg', DATA_IDX=[100,8000], /VERBOSE
;
; MODIFICATION HISTORY:
;   Version 1.0, 17-APR-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 02-MAY-2013, DD: Enabled GPU detection and added keyword 'NO_GPU' to prevent GPU use
;   Version 1.2, 16-MAY-2013, DD: added filter wheel and '*_idx' keywords
;   Version 1.3, 27-JUN-2013, DD: removed keyword 'ONESIDE' and added 'LAMBDA' + 'BANDWIDTH'
;   Version 1.4, 30-JUN-2013, DD: added keyword 'SAVE' and modified call to 'lbti_img2flx.pro' according to new version
;   Version 1.5, 07-JUL-2013, DD: moved parsing of keyword information into 'lbti_readdata.pro'
;   Version 1.6, 29-JUL-2013, DD: added keyword 'N_CLIP'
;   Version 1.7, 01-AUG-2013, DD: added keyword 'BCK_CEN' allowing an user-defined position for the background region
;   Version 1.8, 08-AUG-2013, DD: Included path definition in common 'GLOBAL' + added keyword 'LMIRCAM'
;   Version 1.9, 13-AUG-2013, DD: Improved memory usage
;   Version 2.0, 19-SEP-2013, DD: Added keyword OB_MODE
;   Version 2.1, 10-OCT-2013, DD: Implemented new formalism for L1 files
;   Version 2.2, 14-OCT-2013, DD: Implemented image binning factor (n_bin)
;   Version 2.3, 02-NOV-2013, DD: Implemented the possibility to skip image calibration (+replaced L1_FILES keyword by SKIP_RED)
;   Version 2.4, 03-NOV-2013, DD: Added keyword 'NO_CENTER'
;   Version 2.5, 21-NOV-2013, DD: Renamed keyword METHOD to FLX_METHOD and added keyword FIT_METHOD
;   Version 2.6, 24-NOV-2013, DD: Added keyword OFFSET
;   Version 2.7, 13-JAN-2014, DD: Implemented masterlog file and incremental read of the L0 data files
;   Version 2.8, 14-JAN-2014, DD: Modified the use of the OBSTYPE keyword and implemented XCEN/YCEN as two-element vectors
;   Version 2.9, 15-JAN-2014, DD: Added keyword MASTERLOG to force the creation of a new masterlog file
;   Version 3.0, 07-FEB-2014, DD: Removed obsolote filter keywords and improved speed of image calibration
;   Version 3.1, 23-FEB-2014, DD: Splitted background subtraction and frame centering
;   Version 3.2, 13-MAR-2014, DD: Replaced NO_OVERLAP by OVERLAP
;   Version 3.3, 13-APR-2014, DD: Implemented incremental reduction for null computation
;   Version 3.4, 23-MAY-2014, DD: Added keyword BCKG_MODE and SKIP_FLX
;   Version 3.5, 25-MAY-2014, DD: Added keyword NULL_MODE and NULL_RATIO
;   Version 4.0, 27-MAY-2014, DD: Now compute flux per pointing to avoid producing too large data cube when reducing all data at once
;   Version 4.1, 29-MAY-2014, DD: Added keyword N_FRBCK
;   Version 4.2, 17-SEP-2014, DD: Added keyword NO_FLAT and NO_DARK + removed NO_GPU
;   Version 5.0, 18-SEP-2014, DD: Completely new implementation of intermediate L0 files.
;   Version 5.1, 24-SEP-2014, DD: Added keyword ADI_MODE, SKIP_ADI and SKIP_NULL
;   Version 5.2, 26-SEP-2014, DD: Added keyword N_TRANS
;   Version 5.3, 14-OCT-2014, DD: Added keywords VISI_RATIO, RIN_INIT, SIG_FLX, SIG_POS, SIG_PL, SIG_VIS, and N_COADD
;   Version 5.4, 21-OCT-2014, DD: Renamed NULL_RATIO by KEEP_RATIO for more general use cases (now including Fizeau lucky imaging)
;   Version 5.5, 22-OCT-2014, DD: Now save nods independently
;   Version 5.6, 28-OCT-2014, DD: Added keyword PSF_FILE
;   Version 5.7, 29-OCT-2014, DD: Added keyword BIAS_ESTIM
;   Version 5.8, 01-NOV-2014, DD: Modfified call to LBTI_IMGBCK to use new background subtraction strategy
;   Version 5.9, 05-NOV-2014, DD: Calibrated L0 data now saved in LBTI_IMGBCK + added log file for calibrated L0 files
;   Version 6.0, 11-NOV-2014, DD: Added PHASECam setpoint in the definition of instrumental configurations
;   Version 6.1, 17-NOV-2014, DD: Added keyword NULL_LIM
;   Version 6.2, 21-NOV-2014, DD: Added keyword N_BTSTRP
;   Version 6.3, 12-JAN-2015, DD: Added NO_DARK and NO_FLAT to drs common
;   Version 6.4, 13-JAN-2015, DD: Added keyword NO_BPM
;   Version 6.5, 03-FEB-2015, DD: Corrected nodding number attribution for partial data directories
;   Version 6.6, 06-FEB-2015, DD: Added keyword NULL_COR
;   Version 6.7, 11-FEB-2015, DD: Moved datalog creation to its own procedure
;   Version 6.8, 16-FEB-2015, DD: Moved most keywords to config files
;   Version 6.9, 13-MAR-2015, DD: Implemented multi-threading for null computation
;   Version 7.0, 04-APR-2015, DD: Now remove open-loop frames right from the start (not even read) --- ORR version
;   Version 7.1, 25-APR-2015, DD: Now reduce the raw images nod per nod instead by nod pairs. That's slower but minimize the effective nodding period (=> ELFN mitigation)
;   Version 7.2, 05-MAY-2015, DD: Now OB ID properly assigned when multiple config ID
;   Version 7.3, 05-JUN-2015, DD: Added frame registration function
;   Version 7.4, 16-SEP-2015, DD: Now also save pure background intermediate L0 files (obstype=3)
;   Version 7.5, 02-OCT-2015, DD: Added keyword OB_IDX
;   Version 7.6, 29-OCT-2015, DD: Now read the config ID from the masterlog file
;   Version 7.7, 15-NOV-2015, DD: Updated call to GET_CNF and DECLARE_PATH
;   Version 7.8, 27-NOV-2015, DD: Now use the elevation rather than the time with BCKG_MODE=3
;   Version 7.9, 02-DEC-2015, DD: Cleaned up code (version released to NexSCI)
;   Version 8.0, 22-DEC-2015, DD: Included call to LBTI_READMASTERLOG + modified call to LBTI_READDATA + improved memory usage
;   Version 8.1, 26-MAR-2016, DD: Corrected bug with bad pixel map computation + added keyword RENEW + modified call to LBTI_READDATA
;   Version 8.2, 28-APR-2016, DD: Improved bad pixel map and flat computation
;   Version 8.3, 26-MAY-2016, DD: Added diagnsotic plots
;   Version 8.4, 03-JUL-2016, DD: Implemented use of central value for background selection
;   Version 8.5, 06-JUL-2016, DD: Now avoid reading multiple times the same files when doing the background subtraction
;   Version 8.6, 11-JUL-2016, DD: Now properly assign background OBs in all cases
;   Version 8.7, 03-OCT-2016, DD: Improved backwards compatibility
;   Version 8.8, 11-FEB-2017, DD: Improved speed and added option for frame mode
;   Version 8.9, 04-APR-2017, DD: Added call to GET_TGT
;   Version 9.0, 06-APR-2017, DD: Added call to LBTI_FIXDATA
;   Version 9.1, 06-JUN-2017, DD: Corrected multiple bugs in file association
;   Version 9.2, 18-JUL-2018, DD: Now always use master bad pixel map rather than the one from the same pointing
;   Version 9.3, 31-JAN-2020, DD: Cleaned some keyword definition
;   Version 9.4, 10-AUG-2020, DD: Now do not limit the aperture size in the X direction when drs.sky_col is set to 1
;   Version 9.5, 15-OCT-2023, DD: Updated for FRA_MODE=2 (i.e., PCA background subtraction)
;   Version 9.6, 04-SEP-2024, DD: Added FRA_MODE=3 (mean) and FRA_MODE=4 (median). FRA_MODE=0 is now legacy
;   Version 9.7, 16-SEP-2024, DD: Activated PCA-background frames for PCA too

pro LBTI_DRS, date, cfg_file, $ ; Mandatory inputs (date and config file)
  bad_idx = bad_idx, bckg_idx = bckg_idx, dark_idx = dark_idx, data_idx = data_idx, flat_idx = flat_idx, nod_idx = nod_idx, ob_idx = ob_idx, $ ; Optional inputs (file index, superseed keywords)
  data_path = data_path, bckg_path = bckg_path, dark_path = dark_path, flat_path = flat_path, $ ; Optional inputs (file path, superseed "*_idx" keywords)
  masterlog = masterlog, renew = renew, skip_adi = skip_adi, skip_flx = skip_flx, skip_null = skip_null, skip_red = skip_red, $ ; Running keywords
  log_file = log_file, no_multi = no_multi, no_save = no_save, plot = plot, verbose = verbose ; Control outputs
  compile_opt idl2

  ; DEFINE GLOBAL PARAMETERS
  ; ************************
  ; Global parameters are grouped together into structures, defined hereafter.
  ; The advantages of using structures for global variables is:
  ; a. They are easily modified or extended, via the routine in which they are defined, without the need to update the common blocks
  ; b. Once defined, array length and type of tags are protected, only their value may be modified in the code
  ; c. It is difficult to modify the value inadvertently, since local variables and global variables are clearly distinguishable:
  ; after definition a global variable is accessible as: 'structure_name.variable_name'
  ; The different structures are:
  ; - prm:  invariable astronomical and physical parameters
  ; - cnf:  parameters defining the instrumental configuration
  ; - wav:  parameters pertaining to the propagation and detection of light waves in the array
  ; - tgt:  parameters defining the target star and planetary system components
  ; - pth:  parameters defining the data and result paths
  ; - drs:  parameters related to data reduction (extracted from the config file)
  ; - log:  information contained in the masterlog

  common GLOBAL, prm, cnf, wav, tgt, pth, drs, log
  on_error, 0

  ; DEFINE GLOBAL VARIABLES
  ; ***********************

  ; Astronomical and physical constants
  GET_PRM, prm

  ; Read config file with reduction parameters
  GET_DRS, drs, 'nodrs/cfg/' + cfg_file

  ; Obtain instrumental paramaters
  GET_CNF, cnf, instrum = drs.instrum

  ; Recover the IDL running path
  DECLARE_PATH, pth, instrum = drs.instrum, path_file = drs.path_file

  ; KEYWORD AND INPUT SANITY CHECK
  ; ******************************

  ; Critical parameters
  if strlen(date) ne 6 and not keyword_set(data_path) then message, 'Invalid input date: must be composed of 6 digits'
  if not keyword_set(data_path) then data_path = pth.root_data + date + pth.sep
  pth.data_path = data_path

  ; Keywords
  if keyword_set(skip_adi) then drs.skip_adi = 1
  if keyword_set(skip_flx) then drs.skip_flx = 1
  if keyword_set(skip_null) then drs.skip_null = 1
  if keyword_set(skip_red) then drs.skip_red = 1
  if not keyword_set(plot) then plot = 0
  if keyword_set(verbose) then info = verbose $
  else info = 0

  ; Reformat input date for later
  date_obs = '20' + strmid(date, 0, 2) + '-' + strmid(date, 2, 2) + '-' + strmid(date, 4, 2)

  ; Parse additional info to drs structure
  drs = create_struct(drs, 'VERSION', 9.6, 'DATE', '04-SEP-2024', 'DATE_OBS', date_obs)

  ; INITIALIZE LOG AND TERMINAL OUTPUT
  ; **********************************

  ; First output is a version number of the current code
  if info gt 0 then begin
    print, ' '
    print, 'Nodrs - Version ' + string(drs.version, format = '(F3.1)') + ' -- ' + drs.date + ' -- Denis Defrère - Steward Observatory (denis@lbti.org)'
    print, 'Data processed on ' + systime()
    print, ' '
  endif

  ; Print also the version number in the log file (one file per date in which the results of different reduction are appended)
  if keyword_set(log_file) then begin
    ; Log directory will be the same as the one with the L1 files
    if max(drs.aper_rad) ne 0 then tmp_label = '_APR' else tmp_label = ''
    log_dir = pth.l1Fits_path + pth.sep + date_obs + drs.dir_label + tmp_label + pth.sep
    if not file_test(log_dir) then file_mkdir, log_dir
    ; Create log file
    log_file = log_dir + date_obs + '.txt'
    openw, lun, log_file, /get_lun, width = 800, /append
    printf, lun, ' '
    printf, lun, '******************************************************************************************************* '
    printf, lun, ' '
    printf, lun, 'Nodrs - Version ' + string(drs.version, format = '(F3.1)') + ' -- ' + drs.date + ' -- Denis Defrère - Steward Observatory (denis@lbti.org)'
    printf, lun, 'Data processed on ' + systime()
    printf, lun, ' '
    printf, lun, '******************************************************************************************************* '
    printf, lun, ' '
  endif else lun = -1

  ; READ OR CREATE MASTERLOG FILE
  ; *****************************

  ; Skipped if processing nulling data right away
  if drs.skip_adi eq 0 or drs.skip_red eq 0 then begin
    ; --- Create masterlog file if necessary
    mlog_file = data_path + 'masterlog.dat'
    if not file_test(mlog_file) or keyword_set(masterlog) then begin
      if file_test(data_path + 'config_id.dat') then file_delete, data_path + 'config_id.dat' ; Delete config ID file to tell LBTI_MASTERLOG to create new one
      if file_test(data_path + 'nod_id.dat') then file_delete, data_path + 'nod_id.dat' ; Delete nod ID file to tell LBTI_MASTERLOG to create new one
      print, '  Creating new masterlog file. This may take a while...'
      LBTI_MASTERLOG, data_path, /idl
    endif else if file_test(mlog_file) then LBTI_MASTERLOG, data_path, /append, /idl

    ; --- Run consistancy check on the masterlog (now done separately)
    ; LBTI_FIXDATA, data_path

    ; --- Read masterlog file
    LBTI_READMASTERLOG, mlog_file, log, bad_idx = bad_idx, bckg_idx = bckg_idx, dark_idx = dark_idx, flat_idx = flat_idx, nod_idx = nod_idx

    ; --- Flag out bad frames
    idx_bad = where(log.flag eq 'B', n_b)
    if n_b gt 0 then log.bad_id[idx_bad] = 1

    ; --- Flag out transition frames
    if drs.n_trans then begin
      idx_trans = where(log.nod_id ne shift(log.nod_id, 1))
      for i_trans = 0, drs.n_trans - 1 do log.bad_id[idx_trans + i_trans] = 1
    endif

    ; --- Flag open-loop frames to be removed (AO and phase loops where appropriate. Keep background frames, i.e. obstype = 3)
    if drs.skip_open then begin
      idx_open = where(log.flag eq 'O' and log.obstype ne 3, n_o)
      if n_o gt 0 then log.bad_id[idx_open] = 1
    endif

    ; --- Plot diagnostic info
    if plot gt 0 then begin
      result_path = pth.result_path + pth.sep + 'diagnostic' + pth.sep + 'background' + pth.sep + date_obs + pth.sep
      if not file_test(result_path) then file_mkdir, result_path
      PLOTALL, log.file_id, log.lbtalt, 0, name = result_path + date + '-elevation_vs_file-id', xtitle = 'Target file number', ytitle = 'Telescope elevation [deg]', /bold, /no_fft, /no_histo, /eps
      PLOTALL, log.file_id, log.cv, 0, name = result_path + date + '-cv_vs_file-id', xtitle = 'Target file number', ytitle = 'Central value [ADU]', /bold, /no_fft, /no_histo, /eps
    endif
  endif

  ; 1a. READ AND SAVE CALBRATION IMAGES (BPM, DARK, and FLAT)
  ; *********************************************************

  ; --- Skip calibration if requested and restore the calibrated L0 file from disk
  t0 = systime(1)
  if drs.skip_red eq 1 then begin
    if info gt 0 then print, '1. Restoring calibrated L0 images from disk.'
    goto, skip_imgcal
  endif else if info gt 0 then print, '1. Reading and calibrating L0 data.'
  if info gt 0 then print, ' '

  ; --- Create master BPM, DRK, FLAT, and LIN files
  if drs.no_bpm ne 1 or drs.no_dark ne 1 or drs.no_flat ne 1 then LBTI_IMGCAL, date_obs, dark_path = dark_path, flat_path = flat_path, info = info

  ; 1b. EXTRACT DATA IN RANGE AND DEFINE METRICS
  ; ********************************************

  ; --- Extract useful data and define the number of background pairs (datatype = 3 for backward compatibility)
  file_id = log.file_id
  bad_id = log.bad_id
  datatype = log.datatype
  if keyword_set(data_idx) then idx_data = where(file_id ge data_idx[0] and file_id le data_idx[1] and bad_id ne 1 and ((datatype eq 0 or datatype eq 3)), n_data) $
  else idx_data = where(bad_id ne 1 and ((datatype eq 0 or datatype eq 3)), n_data)
  if n_data le 0 then message, 'No data found in DATA_IDX range!'

  ; --- Deal with BCKG_IDX
  if keyword_set(bckg_idx) then begin
    idx_bckg = where(file_id ge bckg_idx[0] and file_id le bckg_idx[1], n_bckg)
    if n_bckg gt 0 then idx_data = [idx_data, idx_bckg] else print, 'No background frames in BCKG_IDX range!'
  endif

  ; --- Read useful information
  file_id = file_id[idx_data]
  datatype = datatype[idx_data]
  objname = (log.obj_name)[idx_data]
  lbt_alt = (log.lbtalt)[idx_data]
  pt_id = (log.pt_id)[idx_data]
  cfg_id = (log.cfg_id)[idx_data]
  nod_id = (log.nod_id)[idx_data]
  obstype = (log.obstype)[idx_data]
  time_obs = (log.time_obs)[idx_data]
  cv = (log.cv)[idx_data]
  flag = (log.flag)[idx_data]
  n_xpix = (log.n_xpix)[idx_data] ; used to find the corresponding DARK AND FLAT
  n_ypix = (log.n_ypix)[idx_data] ; used to find the corresponding DARK AND FLAT
  log = {file_id: file_id, nod_id: nod_id, cfg_id: cfg_id, pt_id: pt_id, cv: cv} ; Redefine log with minimum required information to save memory

  ; --- Compute number of different nods
  nod_uniq = nod_id[uniq(nod_id, sort(nod_id))]
  idx_pos = where(nod_uniq ge 0, n_pos)
  if n_pos gt 0 then nod_uniq = nod_uniq[idx_pos] else message, 'No valid nods!'
  n_nod = n_elements(nod_uniq)

  ; --- Derive current MJD
  hr = double(strmid(time_obs, 0, 2)) + double(strmid(time_obs, 3, 2)) / 60d + double(strmid(time_obs, 6, 2)) / 3.6d+3 + double(strmid(time_obs, 9, 3)) / 3.6d+6
  JDCNV, fix(strmid(date_obs, 0, 4)), fix(strmid(date_obs, 5, 2)), fix(strmid(date_obs, 8, 2)), hr, jd
  mjd_obs = jd - 2400000.5d ; mjd

  ; --- Remove two NODs if data are being reduced in real-time (use 30 minutes as decision criterion)
  if (max(mjd_obs) + 0.5 / 24) gt (systime(/julian, /utc) - 2400000.5d) then n_nod -= 2

  ; Print WARNING if uneven number of nods when bckg_mode = 1
  if n_nod mod 2 ne 0 and abs(drs.bckg_mode) eq 1 then begin
    message, 'Uneven number of nods not compatible with bckg_mode = 1. Removing the last nod.', /continue
    n_nod -= 1
  endif

  ; --- Define metric used to select the background frames (compute mjd_obs if bckg_sel is not 1)
  case drs.bckg_sel of
    0: metric = mjd_obs
    1: metric = lbt_alt
    2: metric = cv
    else: message, 'Undefined background selection mode (BCKG_SEL)'
  endcase

  ; 2. READ AND REDUCE SCIENCE IMAGES
  ; *********************************

  ; --- Print info to screen
  if info gt 0 then print, ' '
  if info gt 0 then print, '   Processing ' + string(n_data, format = '(I0)') + ' data files. This may take a while...'

  ; --- Loop over the nods (have to proceed per nod here to avoid reading and copying large data arrays)
  t0 = systime(1)
  for i_nod = 0, n_nod - 1 do begin
    ; If file already exists, skip
    case drs.img_mode of ; Image combination mode (0: median, 1: mean, 2: resistant mean)
      0: sub_dir = 'med'
      1: sub_dir = 'mean'
      2: sub_dir = 'resm'
      else: message, 'Unknown IMG_MODE'
    endcase
    if not file_test(pth.l0Fits_path + date_obs + pth.sep + sub_dir + pth.sep + '*_N' + string(nod_uniq[i_nod], format = '(I03)') + '*IMG.fits') or keyword_set(renew) then begin
      ; Init time stamp
      tn = systime(1)

      ; Derive file index of frames in this nod
      idx_nod = where(nod_id eq nod_uniq[i_nod], n)
      if n lt drs.min_frnod then begin
        message, 'No enough frames in nod ' + string(i_nod, format = '(I0)'), /continue
        goto, skip_thisnod
      endif

      ; Derive current pointing ID number and frame size (needed to find dark and flat frames)
      pt_cur = pt_id[idx_nod[0]]
      nx = n_xpix[idx_nod[0]]
      ny = n_ypix[idx_nod[0]]

      ; Print info to screen
      if info gt 0 then print, '   Reducing nod ' + string(nod_uniq[i_nod], format = '(I0)') + ' (' + string(i_nod + 1, format = '(I0)') + '/' + $
        string(n_nod, format = '(I0)') + ', ' + string(n, format = '(I0)') + ' frames)'

      ; Derive obstype. If 3 (pure background), this nod can be skipped unless it is associated with a nulling nod (because it is used by the NSC)
      ; In order to decide whether to reduce this nod or not, we look for nulling frames in the curring pointing.
      obstype_nod = obstype[idx_nod[0]] ; obstype of the current nod
      tmp = where(obstype[where(pt_id eq pt_cur)] eq 2, n_null_ob) ; n_null_ob will contain the number of null OBs in this pointing
      if obstype_nod eq 3 and n_null_ob eq 0 then goto, skip_thisnod ; skip if background OB and no associated NULL OB in the same poiting

      ; Now derive file index of frames in the background nod(s).
      ; The background nods must belong to the same pointing and have the same obstype (or an obstype of 3).
      if obstype_nod ne 3 then begin
        case abs(drs.bckg_mode) of
          1: idx_bck = where(pt_id eq pt_cur and nod_id eq nod_uniq[i_nod + (-1) ^ i_nod] and ((obstype_nod eq obstype or obstype eq 3)), n_bck) ; Background subtraction by nod pairs
          2: idx_bck = where(pt_id eq pt_cur and ((nod_id eq nod_uniq[i_nod] - 1 or nod_id eq nod_uniq[i_nod] + 1)) and ((obstype_nod eq obstype or obstype eq 3)), n_bck) ; Background subtraction using adjacent nods
          3: idx_bck = where(pt_id eq pt_cur and obstype eq 3, n_bck) ; Background subtraction using dedicated background frames
          else: begin
            idx_bck = where(pt_id eq pt_cur and nod_id eq nod_uniq[i_nod] + (-1) ^ i_nod, n_bck)
            if n_bck le 0 then idx_bck = where(pt_id eq pt_cur and nod_id eq nod_uniq[i_nod] - (-1) ^ i_nod, n_bck)
          end
        endcase
        ; If no adjacent background frames, then look within pointing
        if n_bck le 0 then idx_bck = where(pt_id eq pt_cur and obstype eq 3, n_bck)
        ; If still no background frames, then look everywhere
        if n_bck le 0 then idx_bck = where(obstype eq 3, n_bck)
      endif else idx_bck = where(pt_id eq pt_cur and ((nod_id eq nod_uniq[i_nod] - 1 or nod_id eq nod_uniq[i_nod] + 1)), n_bck)

      ; Remove outliers based on central value (5 sigma RMS). Not necessary
      ; cv_nod  = cv[idx_nod]
      ; cv_bck  = cv[idx_bck]
      ; AVGSDV, [cv_bck, cv_nod], cv_avg, cv_rms
      ; idx_bad = WHERE(ABS(cv_nod-cv_avg) GT 5*cv_rms, COMPLEMENT=idx_ok, n_bad) & IF n_bad GT 0 THEN idx_nod = idx_nod[idx_ok] & n     = N_ELEMENTS(idx_nod)
      ; idx_bad = WHERE(ABS(cv_bck-cv_avg) GT 5*cv_rms, COMPLEMENT=idx_ok, n_bad) & IF n_bad GT 0 THEN idx_bck = idx_bck[idx_ok] & n_bck = N_ELEMENTS(idx_bck)

      ; Keep only background frames within maximum allowed time
      if drs.max_time gt 0 then begin
        idx_time = where(abs(mjd_obs[idx_bck] - mjd_obs[idx_nod[0]]) lt drs.max_time / 24. / 60. or abs(mjd_obs[idx_bck] - mjd_obs[idx_nod[-1]]) lt drs.max_time / 24. / 60., n_bck)
        if n_bck gt drs.min_frnod then idx_bck = idx_bck[idx_time]
      endif

      ; Keep only the closest frames based on the chosen metric (e.g., time, elevation, central value)
      if drs.n_frbck ne 0 and drs.n_frbck lt n_bck then begin
        if drs.n_frbck eq -1 then n_bck = (n < n_bck) else n_bck = drs.n_frbck ; if -1, keep the same number of frames as in the current nod
        idx_min = sort(abs(metric[idx_bck] - mean(metric[idx_nod]))) ; we compare to the mean value of the current NOD
        idx_bck = idx_bck[idx_min[0 : n_bck - 1]]
      endif

      ; Skip if no background frames
      if n_bck lt drs.min_frnod then begin
        message, 'No background frames for nod ' + string(i_nod, format = '(I0)'), /continue
        goto, skip_thisnod
      endif

      ; Concatenate file ID to be passed to LBTI_READDATA (and sort them)
      ; PRINT, MIN(idx_nod), MAX(idx_nod), N_ELEMENTS(idx_nod)
      ; PRINT, MIN(idx_bck), MAX(idx_bck), N_ELEMENTS(idx_bck)
      idx_all = [idx_nod, idx_bck]
      fid_all = file_id[idx_all]
      fid_all = fid_all[sort(fid_all)]
      min_id = min(fid_all)
      max_id = max(fid_all)

      ; Now look if some files have already been read in the previous iteration (and remove them from the files to be read in this iteration)
      ; FILE ID must match, not file index!
      if keyword_set(hdr_data) then begin
        MATCH, fid_all, hdr_data.file_id, idx_ok, idx_prev, count = n_prev ; hdr_data.file_id is defined below and comes from the previous iteration
        if n_prev gt 0 then begin
          ; File to be read
          fid_all[idx_ok] = -1
          idx_read = where(fid_all ne -1, n_read)
          ; IF info GE 3 THEN PRINT, '    Frames already read in previous iteration', n_prev, n_read
          if n_read gt 0 then begin
            fid_all = fid_all[idx_read]
            ; Retrieve frames already read
            ; img_data and hdr_data are defined below and come from the previous iteration
            img_prev = img_data[*, *, idx_prev]
            hdr_prev = hdr_data[idx_prev]
            ; Undefine arrays to save memory
            UNDEFINE, img_data, hdr_data
          endif
        endif else n_read = 1
      endif else n_read = 1

      ; Read OB's data (only if new frames to read)
      if n_read gt 0 then begin
        if info gt 0 then print, '    - reading L0 data.'
        img_data = LBTI_READDATA(data_path, crop = drs.pre_crop, data_idx = fid_all, hdr_data = hdr_data, mlog_data = log, info = info, plot = plot)

        ; Reduce L0 images (dark subtraction, flat fielding, wavelentgh calibration, pickup noise, etc...)
        if info gt 0 then print, '    - calibrating L0 images.'

        ; --- Read last dark and flat files (assumes that we used only a single size throught the night)
        if drs.no_bpm ne 1 then LBTI_READCALIMG, pth.bpm_path, date_obs, [nx, ny], -999, bpm, hdr_bpm, crop = drs.pre_crop else if info gt 0 then print, '      - no bad pixel correction' ; -999 will tell this routine to use the master bad pixel map
        if drs.no_bpm ne 1 and (size(bpm))[0] le 1 then goto, skip_thisnod
        if drs.no_dark ne 1 then LBTI_READCALIMG, pth.dark_path, date_obs, [nx, ny], pt_cur, drk, hdr_drk, crop = drs.pre_crop else if info gt 0 then print, '      - no dark subraction'
        if drs.no_dark ne 1 and (size(drk))[0] le 1 then goto, skip_thisnod
        if drs.no_flat ne 1 then LBTI_READCALIMG, pth.flat_path, date_obs, [nx, ny], pt_cur, flt, hdr_flt, crop = drs.pre_crop else if info gt 0 then print, '      - data not flat fielded'
        if drs.no_flat ne 1 and (size(flt))[0] le 1 then goto, skip_thisnod

        ; --- Flat field each frame of this nod
        img_data = LBTI_IMGRED(temporary(img_data), temporary(hdr_data), img_bpm = bpm, img_drk = drk, img_flt = flt, hdr_drk = hdr_drk, hdr_flt = hdr_flt, $ ; Images and correspondind header
          hdr_data = hdr_data, $ ; Output keywords
          log_file = log_file, info = info, plot = plot) ; Input keywords
        if (size(img_data))[0] eq 0 then goto, skip_thisnod

        ; Undefine arrays to save memory
        UNDEFINE, bpm, drk, flt, hdr_bpm, hdr_drk, hdr_flt

        ; Now concatenate files read this iteration with files already read before
        if keyword_set(hdr_prev) gt 0 then begin
          if n_prev gt 0 then begin ; n_prev not defined if i_nod = 0
            ; Concatenate
            img_data = [[[img_data]], [[img_prev]]]
            hdr_data = [hdr_data, hdr_prev]
            ; Sort (better have data in chronolical order for later)
            idx_srt = sort(hdr_data.file_id)
            hdr_data = temporary(hdr_data[idx_srt])
            img_data = temporary(img_data[*, *, idx_srt])
          endif
        endif
      endif

      ; Only keep data in the useful range
      idx_ok = where(hdr_data.file_id ge min_id and hdr_data.file_id le max_id)
      img_data = img_data[*, *, idx_ok]
      hdr_data = hdr_data[idx_ok]

      ; Assign background nod flag
      idx_bck = where(hdr_data.nod_id ne nod_uniq[i_nod], complement = idx_nod)
      hdr_data[idx_nod].bck_nod = 0
      hdr_data[idx_bck].bck_nod = 1

      ; Perform background subtraction and save data
      if info gt 0 then print, '    - performing background subtraction.'
      LBTI_IMGBCK, img_data, hdr_data, bckg_path = bckg_path, log_file = log_file, info = info, no_save = no_save, plot = plot ; Input keywords

      ; Jump point of problem with nod computation
      skip_thisnod:

      ; Print time to process this nod
      if info ge 3 then print, '   Time to process this nod : ', systime(1) - tn
    endif else print, '   Nod ' + string(nod_uniq[i_nod], format = '(I0)') + ' already exists'
  endfor

  ; Print info to screen
  t1 = systime(1)
  if info gt 2 then print, 'Time to perform L0 data calibration :', t1 - t0
  if info gt 0 then print, ' '

  ; Jumping point if restoring calibrated L0 files
  skip_imgcal:

  ; 3. FLUX COMPUTATION (also needed for frame centering and cropping before ADI processing)
  ; *******************

  ; --- Skip flux computation if requested and restore flux from disk
  t2 = systime(1)
  if drs.skip_flx eq 1 and drs.skip_vis eq 1 then begin
    if info gt 0 then begin
      if drs.skip_flx eq 1 then print, '2. Restoring flux measurements from disk.' $
      else if drs.skip_vis eq 1 then print, '2. Restoring visibility measurements from disk.'
    endif
    goto, skip_flx
  endif
  if info gt 0 then print, ' '

  ; --- Restore LO log file and read useful data
  datalog = pth.l0Fits_path + date_obs + pth.sep + 'datalog.sav'
  if not file_test(datalog) then begin
    message, 'No data log file found!', /continue
    RETURN
  endif else restore, datalog

  ; --- Extract relevant information for all data (used to computed OB id of photometric and background frames)
  objname = data_r.objname
  nod_id = data_r.nod_id
  nod_id0 = nod_id
  cfg_id = data_r.cfg_id
  cfg_id0 = cfg_id
  pt_id = data_r.pt_id
  pt_id0 = pt_id
  obstype = data_r.obstype
  xcen_sx = data_r.xcen_sx
  xcen_dx = data_r.xcen_dx
  ycen_sx = data_r.ycen_sx
  ycen_dx = data_r.ycen_dx
  min_fid = data_r.min_fid
  max_fid = data_r.max_fid
  flag = data_r.flag

  time_obs = data_r.ut_time
  ut_time = double(strmid(time_obs, 0, 2)) + double(strmid(time_obs, 3, 2)) / 60d + double(strmid(time_obs, 6, 2)) / 3.6d+3 + double(strmid(time_obs, 9, 3)) / 3.6d+6

  ; --- Read target information (or query online catalogs if not present)
  tgt_uniq = objname[uniq(objname, sort(objname))]
  GET_TGT, tgt_uniq, tgt, database = pth.input_path + drs.database
  struct_add_field, tgt, 'calfor', 'tbd'
  struct_add_field, tgt, 'maxapr', 1d3

  ; --- Assign calibrators to their science object (needed because of EEID computation).
  ; --- Also assign the maximum aperture radius possible for each star.
  pt_uniq = pt_id[uniq(pt_id, sort(pt_id))]
  n_pt = n_elements(pt_uniq)
  flag_pt = strarr(n_pt)
  tgt_pt = flag_pt
  max_apr = intarr(n_pt)
  for i = 0, n_pt - 1 do begin
    idx = where(pt_id eq pt_uniq[i])
    flag_pt[i] = flag[idx[0]]
    tgt_pt[i] = strlowcase(strcompress(objname[idx[0]], /remove_all))
    ; in sub-frame mode, this doesn't capture the distance to the horizontal center but no easy solution for that at this point...
    ; this is also only valid for SX, which encodes the position of the nulling beam. PHOTOMETRIC beam should be close
    idx = where(pt_id eq pt_uniq[i] and ycen_sx ne 0 and xcen_sx ne 0)
    if drs.sky_col eq 0 then max_apr[i] = min(xcen_sx[idx] mod cnf.x_chan) < min(ycen_sx[idx] mod cnf.y_chan) < min(cnf.x_chan - (xcen_sx[idx] mod cnf.x_chan)) < min(cnf.y_chan - (ycen_sx[idx] mod cnf.y_chan)) $
    else max_apr[i] = min(ycen_sx[idx] mod cnf.y_chan) < min(cnf.y_chan - (ycen_sx[idx] mod cnf.y_chan)) ; look only in the Y direction!
    if flag_pt[i] eq 'SCI' then begin
      idx_tgt = where(tgt.name eq tgt_pt[i])
      if tgt[idx_tgt].maxapr gt max_apr[i] then tgt[idx_tgt].maxapr = max_apr[i]
    endif
  endfor
  ; Pre-defined list (should go in a separate file eventually)
  READ_TABLE, 'nodrs/input/calsci_pairs.txt', date_list, ut_min, ut_max, cal_list, sci_list, first = 2, string = [3, 4], separator = ';'

  ; Loop over the pointings
  for i = 0, n_pt - 1 do begin
    ; Current target
    idx_tgt = where(tgt.name eq tgt_pt[i])
    ra_tgt = tgt[idx_tgt].ra
    de_tgt = tgt[idx_tgt].dec
    ; If in pre-defined list, use it!
    idx_ok = where(strlowcase(tgt_pt[i]) eq strlowcase(cal_list) and date eq date_list and ut_time[i] ge ut_min and ut_time[i] le ut_max, n_ok)
    if n_ok eq 1 then tgt[idx_tgt].calfor = strlowcase(sci_list[idx_ok]) else begin
      ; If CAL, find nearest SCIs
      if flag_pt[i] eq 'CAL' then begin
        idx_low = max(where(pt_uniq lt pt_uniq[i] and flag_pt eq 'SCI', n1))
        if n1 gt 0 then tgt_low = tgt_pt[idx_low] else tgt_low = 'none' ; Nearest previous SCI pointings
        idx_up = min(where(pt_uniq gt pt_uniq[i] and flag_pt eq 'SCI', n2))
        if n2 gt 0 then tgt_up = tgt_pt[idx_up] else tgt_up = 'none' ; Nearest following SCI pointings
        ; Current target
        if n1 ne 0 or n2 ne 0 then begin
          ; If not the same, find the one with the closest RA
          if tgt_low ne tgt_up then begin
            idx = where(tgt.name eq tgt_low, n_low)
            if n_low ge 1 then ra_low = tgt[idx].ra else ra_low = 1d9
            if n_low ge 1 then de_low = tgt[idx].dec else de_low = 1d9
            idx = where(tgt.name eq tgt_up, n_up)
            if n_up ge 1 then ra_up = tgt[idx].ra else ra_up = 1d9
            if n_up ge 1 then de_up = tgt[idx].dec else de_up = 1d9
            if ((ra_tgt - ra_low) ^ 2 + (de_tgt - de_low) ^ 2) lt ((ra_tgt - ra_up) ^ 2 + (de_tgt - de_up) ^ 2) then tgt[idx_tgt].calfor = tgt_low else tgt[idx_tgt].calfor = tgt_up
          endif else tgt[idx_tgt].calfor = tgt_up
        endif else begin
          ; Only needed if aper_rad is 0
          if drs.aper_rad eq 0 then begin
            read, ttt, prompt = 'Enter science object for ' + tgt_pt[i] + ' : '
            tgt[idx_tgt].calfor = ttt
          endif else tgt[idx_tgt].calfor = 'UND'
        endelse
      endif
    endelse
  endfor

  ; --- Consitency check
  ; Make sure there is no missing object name
  ; idx_missing = WHERE(obj_name EQ '?' OR obj_name EQ ' ', n_missing)
  ; IF n_missing GE 1 THEN MESSAGE, 'Problem with object name in files : ' + STRING(file_id[idx_missing], FORMAT='(I0)') + ' (' + obj_name[idx_missing] +')'

  ; --- Assign unique ID number to each line of the file
  n_data = n_elements(objname)
  file_id = indgen(n_data) + 1 ; I don't want 0!!!!

  ; --- Assign OB id number to coherent data (so don't consider photometric and background frames in the OB id count)
  idx_coh = where(obstype eq 1 or obstype eq 2, n_coh)
  if n_coh gt 0 then begin
    ob_id = intarr(n_data) - 1 ; Photometry and background frames have an OB id of -1.
    case drs.ob_mode of
      0: dat_id = [transpose(nod_id[idx_coh]), transpose(pt_id[idx_coh]), transpose(cfg_id[idx_coh])]
      1: message, 'OB_MODE=1 not yet implemented'
      2: message, 'OB_MODE=2 not yet implemented'
      3: message, 'OB_MODE=3 not yet implemented (one OB per chop position)'
    endcase
    ob_id[idx_coh] = ATTRIBUTE_ID(dat_id) + 1 ; I don't want 0!!!!
  endif else ob_id = nod_id
  ob_id0 = ob_id

  ; --- Assign background and photometric file ID for each coherent OB (only needed for null OB in principle)
  bck_id = intarr(n_data) ; will contain the corresponding background file ID
  pho_fid = intarr(n_data) ; will contain the corresponding photometric file ID
  nod_pos = fix(ycen_sx / cnf.y_chan) ; this assumes that each NOD ends up in a different vertical channel (which is true for nulling)
  for i_d = 0, n_data - 1 do begin
    ; 1. Compute the corredponding background OB. It must be of the same pointing (and same config ID) but different nod position (+ avoid photometric frames)
    idx_bck = where(pt_id eq pt_id[i_d] and cfg_id eq cfg_id[i_d] and ob_id ne ob_id[i_d] and obstype ne 0 and (nod_pos ne nod_pos[i_d] or obstype eq 3), n_bck)
    if n_bck gt 0 then begin
      file_bck = file_id[idx_bck] ; list of suitable files
      idx_min = where(abs(file_id[i_d] - file_bck) eq min(abs(file_id[i_d] - file_bck))) ; closest ones
      if n_elements(idx_min) gt 1 then begin ; if two, keep the one with more frames
        tmp = max(max_fid[idx_bck[idx_min]] - min_fid[idx_bck[idx_min]], i_max)
      endif else i_max = 0
      bck_id[i_d] = file_bck[idx_min[i_max]] ; retain only one
    endif else bck_id[i_d] = -1 ; no associated background file

    ; 2. Compute the associated photometric file ID. Ideally, it must be of the same pointing but look in different pointings if not found
    idx_pho = where(pt_id eq pt_id[i_d] and cfg_id eq cfg_id[i_d] and obstype eq 0 and xcen_sx ne 0 and ycen_sx ne 0, n_pho)
    if n_pho le 0 then begin
      print, ' No corresponding photometric file for nod ' + string(nod_id[i_d], format = '(I0)') + '. Looking in other pointings.'
      if lun gt 0 then printf, lun, ' No corresponding photometric file for nod ' + string(nod_id[i_d], format = '(I0)') + '. Looking in other pointings.'
      idx_pho = where(objname eq objname[i_d] and cfg_id eq cfg_id[i_d] and obstype eq 0 and xcen_sx ne 0 and ycen_sx ne 0, n_pho)
    endif
    if n_pho gt 0 then begin
      fid_pho = file_id[idx_pho] ; list of suitable file IDs
      idx_min = where(abs(file_id[i_d] - fid_pho) eq min(abs(file_id[i_d] - fid_pho))) ; closest one
      if n_elements(idx_min) gt 1 then begin ; if two, keep the one with the more frames
        tmp = max(max_fid[idx_pho[idx_min]] - min_fid[idx_pho[idx_min]], i_max)
      endif else i_max = 0
      pho_fid[i_d] = fid_pho[idx_min[i_max]] ; retain only one
    endif else begin
      if info gt 0 then begin
        print, ' No corresponding photometric file for nod ' + string(nod_id[i_d], format = '(I0)')
        if lun gt 0 then printf, lun, ' No corresponding photometric file for nod ' + string(nod_id[i_d], format = '(I0)')
      endif
      pho_fid[i_d] = -1 ; no associated photometric NOD
    endelse
  endfor

  ; -- Extract data in DATA_IDX range
  if keyword_set(data_idx) then begin
    idx_data = where(data_r.min_fid ge data_idx[0] and data_r.max_fid le data_idx[1], n_data)
    if n_data gt 0 then begin
      data_r = data_r[idx_data]
      cfg_id0 = cfg_id[idx_data] ; only used for print
      pt_id0 = pt_id[idx_data] ; only used for print
      nod_id0 = nod_id[idx_data] ; only used for print
      ob_id0 = ob_id[idx_data] ; only used for print
    endif else begin
      message, 'No data in the DATA_IDX range', /continue
      RETURN
    endelse
  endif else n_data = n_elements(obj_name)

  ; -- Derive number of unique nod IDs
  nod_uniq = nod_id0[uniq(nod_id0, sort(nod_id0))]
  n_nod = n_elements(nod_uniq) > 1

  ; -- Print info to screen and to the log
  if info gt 0 then begin
    ; Extract data in DATA_IDX range
    wav_eff = data_r.lam_cen
    dit = data_r.int_time
    tgtname = data_r.objname
    ; Print info to screen
    if info gt 2 then begin
      ; Compute the number of different parameters in the input data sequence
      tgt_uniq = tgtname[uniq(tgtname, sort(tgtname))]
      n_tgt = n_elements(tgt_uniq) > 1 ; Distinct object name
      lam_uniq = wav_eff[uniq(wav_eff, sort(wav_eff))]
      n_lam = n_elements(lam_uniq) > 1 ; Distinct integration time
      int_uniq = dit[uniq(dit, sort(dit))]
      n_int = n_elements(int_uniq) > 1 ; Distinct central wavelengths
      pt_uniq = pt_id0[uniq(pt_id0, sort(pt_id0))]
      n_pt = n_elements(pt_uniq) > 1 ; Distinct pointing #ID
      ob_uniq = ob_id0[uniq(ob_id0, sort(ob_id0))]
      n_ob = n_elements(ob_uniq) > 1 ; Distinct OB #ID
      ; Print info to screen
      print, '2. Performing flux/visibility computation over the following parameters:'
      print, ' - date                    : ', strtrim(date_obs)
      print, ' - target name             : ', strtrim(tgt_uniq)
      print, ' - central wavelength [um] : ', string(1d+6 * lam_uniq, format = '(F5.2)')
      print, ' - integration time [ms]   : ', string(1d+3 * int_uniq, format = '(I0)')
      print, ' - number of pointings     : ', string(n_pt, format = '(I0)')
      print, ' - number of nods          : ', string(n_nod, format = '(I0)')
      print, ' - number of OBs           : ', string(n_ob, format = '(I0)')
      print, ' '
      ; Print also the version number in the log file
      if lun gt 0 then begin
        printf, lun, '2. Performing flux/visibility computation over the following parameters:'
        printf, lun, ' - date                    : ', strtrim(date_obs)
        printf, lun, ' - target name             : ', strtrim(tgt_uniq)
        printf, lun, ' - central wavelength [um] : ', string(1d+6 * lam_uniq, format = '(F3.1)')
        printf, lun, ' - integration time [ms]   : ', string(1d+3 * int_uniq, format = '(I0)')
        printf, lun, ' - number of pointings     : ', string(n_pt, format = '(I0)')
        printf, lun, ' - number of nods          : ', string(n_nod, format = '(I0)')
        printf, lun, ' - number of OBs           : ', string(n_ob, format = '(I0)')
        printf, lun, ' '
      endif
    endif

    if drs.flx_mode ne 2 then begin
      print, 'Now computing flux by APERTURE PHOTOMETRY on ' + string(n_data, format = '(I0)') + ' frames.'
      print, '      - Photometric aperture radius [pix] : ' + string(drs.aper_rad, format = '(I0)') + ' (0: 2xEEID)'
      print, '      - Inner background radius [pix]     : ' + string(drs.bck_irad, format = '(I0)') + ' (0: aperture radius)'
      print, '      - Outer background radius [pix]     : ' + string(drs.bck_orad, format = '(I0)') + ' (0: channel edge)'
    endif else print, 'Now computing flux by PSF-FITTING on ' + string(n_data, format = '(I0)') + ' frames.'

    ; Initiate pointing number display
    print, ' '
    print, format = '("Processing nod number ", $)'
  endif
  if lun gt 0 then begin
    if drs.flx_mode ne 2 then begin
      printf, lun, ' '
      printf, lun, 'Flux computed by APERTURE PHOTOMETRY on ' + string(n_data, format = '(I0)') + ' frames:'
      printf, lun, '      - Photometric aperture radius [pix] : ' + string(drs.aper_rad, format = '(I0)') + ' (0: 2xEEID)'
      printf, lun, '      - Inner background radius [pix]     : ' + string(drs.bck_irad, format = '(I0)') + ' (0: aperture radius)'
      printf, lun, '      - Outer background radius [pix]     : ' + string(drs.bck_orad, format = '(I0)') + ' (0: channel edge)'
      printf, lun, ' '
    endif else begin
      printf, lun, ' '
      printf, lun, 'Flux computed by PSF-FITTING on ' + string(n_data, format = '(I0)') + 'frames.'
      printf, lun, ' '
    endelse
  endif

  ; --- If aperture radius is set, save L1 files in a separate directory (no extra label if EEID)
  if max(drs.aper_rad) ne 0 then drs.dir_label = drs.dir_label + '_APR' ;+ STRING(MIN(drs.aper_rad[WHERE(drs.aper_rad GT 0)]), FORMAT='(I0)')

  ; --- Perform flux computation (aperture photometry or PSF fitting)
  for i_nod = 0, n_nod - 1 do begin
    ; Print info to screen
    if info gt 0 then print, nod_uniq[i_nod], format = '($, (x, I0))'

    ; Find most recent files of this nod
    idx_nod = where(nod_id eq nod_uniq[i_nod], n_cfg)
    if n_cfg le 0 then goto, skip_nod

    ; Loop over the files
    for i_cfg = 0, n_cfg - 1 do begin
      ; Derive file ID
      idx_cfg = idx_nod[i_cfg]
      fid_cur = file_id[idx_cfg]
      ob_cur = ob_id[idx_cfg]
      obstype_cur = obstype[idx_cfg]

      ; Check whether this OB has already been reduced
      if obstype_cur eq 0 then tag = 'PHOT'
      if obstype_cur eq 1 then tag = 'FIZ'
      if obstype_cur eq 2 then tag = 'NULL'
      if obstype_cur ge 3 then tag = 'DOANYWAY' ; IF OBSTYPE GE 3, the flux has to be recomputed. I don't remember why but it doesn't hurt to repeat it..

      ; Skip if already exists
      if not file_test(pth.l1Fits_path + drs.date_obs + drs.dir_label + pth.sep + '*_ID' + string(ob_cur, format = '(I03)') + '*' + tag + '*.fits') or keyword_set(renew) then begin
        ; Read raw L0 file (use mean-subtracted, pca-subtracted or raw-subtracted frames)
        case drs.fra_mode of
          0: label = 'bckg' ; Legacy value used in the survey papers and including open-loop frames
          1: label = 'raw'
          2: label = 'pca'
          3: label = 'mean'
          4: label = 'med'
          else: message, 'Undefined frame selection mode (FRA_MODE)'
        endcase
        ; if obstype_cur eq 0 then label = 'bckg' ; 0CT 2023, updated for PCA frames. I don't remember why this bit is here. I leave it here for backward compatibility?

        file_nod = file_search(pth.l0Fits_path + date_obs + pth.sep + label + pth.sep + '*_N' + string(nod_uniq[i_nod], format = '(I03)') + '*' + string(min_fid[idx_cfg], format = '(I06)') + '-' + $
          string(max_fid[idx_cfg], format = '(I06)') + '_IMG.fits', count = n0)
        if n0 gt 0 then begin
          ; If more than 1 file, only keep the latest file and delete the other one
          if n0 gt 1 then begin
            finfo = file_info(file_nod)
            idx = where(finfo.mtime eq max(finfo.mtime), complement = idxd)
            file_delete, file_nod[idxd]
            file_nod = file_nod[idx]
          endif

          ; Read image
          img_data = LBTI_READL0RED(file_nod, hdr_data = hdr_data, info = info)

          ; Current pointing, nod, OB and file ID
          pt_cur = pt_id[idx_cfg]
          cfg_cur = cfg_id[idx_cfg]
          nod_cur = nod_id[idx_cfg]

          ; Parse beam positions
          hdr_data.data[*].xcen_sx = xcen_sx[idx_cfg]
          hdr_data.data[*].ycen_sx = ycen_sx[idx_cfg]
          hdr_data.data[*].xcen_dx = xcen_dx[idx_cfg]
          hdr_data.data[*].ycen_dx = ycen_dx[idx_cfg]

          ; Reset auxiliary positions (used by LBTI_ING2FLX to know at which additional positions to compute the flux)
          xcen_bck = 0
          ycen_bck = 0

          ; Assign OB id to frames and additional beam positions for background
          ; Also assign OB number for all OBs that use this NOD for background/photometry
          case hdr_data.header.obstype of
            0: begin ; Photometric frames (if it exists, take OB ID of associated coherent frames, else take its own asociated ob_id)
              idx_pho = where(pho_fid eq fid_cur and ob_id ne -1, n_ok)
              if n_ok gt 0 then ob_in = ob_id[idx_pho] $ ; contain all OBs that use this photometric file
              else ob_in = -1 ; not used, -1 to skip flux computation below
            end
            1: begin ; Visibility frames
              hdr_data.data[*].ob_id = ob_cur
              idx_bck = where(bck_id eq fid_cur and ob_id ne -1, n_ok)
              if n_ok gt 0 then begin
                ob_in = ob_id[idx_bck] ; contain all OBs that use the current OB as background
                xcen_bck = xcen_sx[idx_bck] ; contain corresponding X positions
                ycen_bck = ycen_sx[idx_bck] ; contain corresponding Y positions
              endif else ob_in = 0 ; not used
            end
            2: begin ; Nulling frames
              hdr_data.data[*].ob_id = ob_cur
              idx_bck = where(pt_id eq pt_cur and cfg_id eq cfg_cur and ob_id ne ob_cur and obstype ne 0 and nod_pos ne nod_pos[idx_cfg], n_bck) ; will return the file ID of the opposite NODs
              if n_bck ge 1 then begin
                nod_bck = nod_id[idx_bck]
                idx_bck = idx_bck[where(abs(nod_cur - nod_bck) eq min(abs(nod_cur - nod_bck)), n_ok)] ; keep only the closest one(s)
              endif else idx_bck = where(bck_id eq fid_cur and ob_id ne -1, n_ok) ; old approach
              if n_ok gt 0 then begin
                ob_in = ob_id[idx_bck] ; contain all OBs that use the current OB as background
                xcen_bck = xcen_sx[idx_bck] ; contain corresponding X positions
                ycen_bck = ycen_sx[idx_bck] ; contain corresponding Y positions
              endif else ob_in = 0 ; not used
            end
            3: begin ; Background frames (only used if assciated to NULL nod)
              idx_bck = where(bck_id eq fid_cur and ob_id ne -1, n_ok)
              if n_ok gt 0 then begin
                ob_in = ob_id[idx_bck] ; contain all OBs that use the current OB as background
                xcen_bck = xcen_sx[idx_bck] ; contain corresponding X positions
                ycen_bck = ycen_sx[idx_bck] ; contain corresponding X positions
              endif else ob_in = -1 ; not used, -1 to skip flux computation below
            end
            else: goto, skip_nod
          endcase
          if info gt 2 then print, ' - PT, NOD, OB | BCK_OB:', pt_cur, nod_cur, ob_cur, ' | ', ob_in

          ; Compute flux/visibility for all frames in this nod and save data if allowed (don't do background frames not used in conjunction with null data)
          if not drs.skip_flx and max(ob_in) ne -1 then LBTI_IMG2FLX, img_data, hdr_data, log_file = lun, info = info, plot = plot, no_save = no_save, ob_in = ob_in, xcen = xcen_bck, ycen = ycen_bck
          if not drs.skip_vis and max(ob_in) ne -1 then LBTI_IMG2VIS, img_data, hdr_data, log_file = lun, info = info, plot = plot, no_save = no_save, ob_in = ob_in, xcen = xcen_bck, ycen = ycen_bck
        endif else print, ' No file for this nod in sub-dir ' + label
      endif ; ELSE PRINT, ' File ' + STRING(fid_cur, FORMAT='(I0)') + '  already exists'
    endfor
    ; Jumping point if no file for this nod ID
    skip_nod:
  endfor

  t3 = systime(1)
  if info gt 0 then print, ' '
  if info gt 2 then print, 'Time to perform flux computation :', t3 - t2

  ; Jumping point if SKIP_FLX is set
  skip_flx:

  ; 4. ADI PROCESSING (from processed L1 frames)
  ; ******************

  ; --- Skip ADI processing if requested
  t1 = systime(1)
  if drs.skip_adi then begin
    if info gt 0 then print, '3. Skipping ADI processing'
    goto, skip_adi
  endif

  ; --- Print info to screen
  if info gt 0 then begin
    print, ' '
    print, '3. Performing ADI processing:'
    print, format = '("Processing nod number ", $)'
  endif

  ; --- Restore LO log file and read useful data
  datalog = pth.l0Fits_path + date_obs + pth.sep + 'datalog.sav'
  if file_test(datalog) then restore, datalog else message, 'No data log file found!'
  objname = data_r.objname
  lam_cen = data_r.lam_cen
  obj_uniq = strcompress(objname[uniq(objname, sort(objname))], /remove_all)
  n_obj = n_elements(obj_uniq) > 1 ; Distinct object name

  ; --- Loop over the objects
  for i_o = 0, n_obj - 1 do begin
    ; --- Derive number of wavelengths (need to be improved)
    tgt_name = obj_uniq[i_o]
    print, tgt_name
    idx_obj = where(objname eq tgt_name)
    lam_cur = lam_cen[idx_obj]
    lam_uniq = lam_cur[uniq(lam_cur, sort(lam_cur))]
    n_lam = n_elements(lam_uniq)
    ; --- Loop over wavelengths
    for i_lam = 0, n_lam - 1 do begin
      ; --- Merge all files
      search_string = '*' + tgt_name + '*' + string(1d+6 * lam_uniq[i_lam], format = '(I0)')
      file = file_search(pth.l1Fits_path + date_obs + drs.dir_label + pth.sep, search_string + '*_IMG.fits', count = n0)
      if n0 lt 1 then goto, skip_file else if not drs.skip_merge then LBTI_MERGEL1FILES, file
      ; --- Frame selection and cropping for all files
      file = file_search(pth.l1Fits_path + date_obs + drs.dir_label + pth.sep, search_string + '*IMG_ALL.fits', count = n0)
      data_file = file_search(pth.l1Fits_path + date_obs + drs.dir_label + pth.sep, search_string + '*DATA_ALL.fits', count = nd0)
      if n0 lt 1 then goto, skip_file else if not drs.skip_sel then LBTI_IMGSEL, file, data_file, info = info, plot = plot ; ,/MEDIAN
      ; --- Derotate (PCA or ADI)
      file = file_search(pth.l1Fits_path + date_obs + drs.dir_label + pth.sep, search_string + '*IMG_SEL.fits', count = n0)
      data_file = file_search(pth.l1Fits_path + date_obs + drs.dir_label + pth.sep, search_string + '*DATA_SEL.fits', count = nd0)
      if n0 lt 1 then goto, skip_file else LBTI_IMGDEROT, file, data_file, plot = plot, verbose = info, /median
      ; --- Skipping point
      skip_file:
    endfor
  endfor

  t2 = systime(1)
  if info gt 2 then print, 'Time to perform ADI processing     :', t2 - t1

  ; Jumping point if SKIP_NULL is set
  skip_adi:

  ; 5. NULL COMPUTATION
  ; *******************

  ; --- Compute null (if nulling data)
  t3 = systime(1)
  if drs.skip_null then begin
    if info gt 0 then print, '4. Skipping null computation'
    goto, skip_null
  endif
  if info gt 0 then print, ' '

  if info gt 0 then print, '4. Now converting flux to null.'
  LBTI_FLX2NULL, date, ob_idx = ob_idx, log_file = lun, info = info, plot = plot, no_multi = no_multi, no_save = no_save, renew = renew
  t4 = systime(1)
  if info gt 2 then print, 'Time to perform null computation     :', t4 - t3

  ; Jumping point if SKIP_NULL is set
  skip_null:

  ; --- Print reduction time to the log file and close it
  if info gt 2 then print, 'Total time to perform data reduction :', systime(1) - t0
  if lun gt 0 then begin
    printf, lun, ' '
    printf, lun, 'Total time to perform data reduction   :', systime(1) - t0
    printf, lun, ' '
    close, lun
    free_lun, lun
  endif
end