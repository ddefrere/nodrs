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

PRO LBTI_DRS, date, cfg_file, $                                                                                                               ; Mandatory inputs (date and config file)
              BAD_IDX=bad_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, DATA_IDX=data_idx,  FLAT_IDX=flat_idx, NOD_IDX=nod_idx, OB_IDX=ob_idx, $ ; Optional inputs (file index, superseed keywords)
              DATA_PATH=data_path, BCKG_PATH=bckg_path, DARK_PATH=dark_path, FLAT_PATH=flat_path, $                                           ; Optional inputs (file path, superseed "*_idx" keywords)
              MASTERLOG=masterlog, RENEW=renew, SKIP_ADI=skip_adi, SKIP_FLX=skip_flx, SKIP_NULL=skip_null, SKIP_RED=skip_red,$                ; Running keywords
              LOG_FILE=log_file, NO_MULTI=no_multi, NO_SAVE=no_save, PLOT=plot, VERBOSE=verbose                                               ; Control outputs

  
; DEFINE GLOBAL PARAMETERS
; ************************
; Global parameters are grouped together into structures, defined hereafter.
; The advantages of using structures for global variables is:
; a. They are easily modified or extended, via the routine in which they are defined, without the need to update the common blocks
; b. Once defined, array length and type of tags are protected, only their value may be modified in the code
; c. It is difficult to modify the value inadvertently, since local variables and global variables are clearly distinguishable:
;    after definition a global variable is accessible as: 'structure_name.variable_name'
; The different structures are:
; - prm:  invariable astronomical and physical parameters
; - cnf:  parameters defining the instrumental configuration
; - wav:  parameters pertaining to the propagation and detection of light waves in the array
; - tgt:  parameters defining the target star and planetary system components
; - pth:  parameters defining the data and result paths
; - drs:  parameters related to data reduction (extracted from the config file)
; - log:  information contained in the masterlog

COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs, log
ON_ERROR, 0

; DEFINE GLOBAL VARIABLES
; ***********************

; Astronomical and physical constants
GET_PRM, prm

; Read config file with reduction parameters
GET_DRS, drs, 'nodrs/cfg/' + cfg_file

; Obtain instrumental paramaters
GET_CNF, cnf, INSTRUM=drs.instrum

; Recover the IDL running path
DECLARE_PATH, pth, INSTRUM=drs.instrum, PATH_FILE=drs.path_file


; KEYWORD AND INPUT SANITY CHECK 
; ******************************

; Critical parameters
IF STRLEN(date) NE 6 AND NOT KEYWORD_SET(DATA_PATH) THEN MESSAGE, 'Invalid input date: must be composed of 6 digits'
IF NOT KEYWORD_SET(DATA_PATH) THEN data_path = pth.root_data + date + pth.sep 
pth.data_path = data_path

; Keywords
IF KEYWORD_SET(SKIP_ADI)  THEN drs.skip_adi  = 1
IF KEYWORD_SET(SKIP_FLX)  THEN drs.skip_flx  = 1
IF KEYWORD_SET(SKIP_NULL) THEN drs.skip_null = 1
IF KEYWORD_SET(SKIP_RED)  THEN drs.skip_red  = 1
IF NOT KEYWORD_SET(PLOT)  THEN plot          = 0
IF KEYWORD_SET(VERBOSE)   THEN info          = verbose $
                          ELSE info          = 0

; Reformat input date for later
date_obs = '20' + STRMID(date, 0, 2) + '-' + STRMID(date, 2, 2) + '-' + STRMID(date, 4, 2)

; Parse additional info to drs structure
drs = CREATE_STRUCT(drs, 'VERSION', 9.5, 'DATE', '15-OCT-2023', 'DATE_OBS', date_obs)


; INITIALIZE LOG AND TERMINAL OUTPUT
; **********************************

; First output is a version number of the current code
IF info GT 0 THEN BEGIN
  PRINT,' '
  PRINT,'Nodrs - Version ' + STRING(drs.version, FORMAT='(F3.1)') +  ' -- ' + drs.date + ' -- Denis Defrère - Steward Observatory (denis@lbti.org)'
  PRINT,'Data processed on '+ SYSTIME()
  PRINT,' '
ENDIF

; Print also the version number in the log file (one file per date in which the results of different reduction are appended)
IF KEYWORD_SET(log_file) THEN BEGIN
  ; Log directory will be the same as the one with the L1 files
  IF MAX(drs.aper_rad) NE 0 THEN tmp_label = '_APR' ELSE tmp_label = ''
  log_dir  = pth.l1fits_path + pth.sep + date_obs + drs.dir_label + tmp_label + pth.sep
  IF NOT FILE_TEST(log_dir) THEN FILE_MKDIR, log_dir
  ; Create log file
  log_file =  log_dir + date_obs + '.txt'
  OPENW, lun, log_file, /GET_LUN, WIDTH=800, /APPEND
  PRINTF,lun, ' '
  PRINTF,lun, '******************************************************************************************************* '
  PRINTF,lun, ' '
  PRINTF,lun, 'Nodrs - Version ' + STRING(drs.version, FORMAT='(F3.1)') +  ' -- ' + drs.date + ' -- Denis Defrère - Steward Observatory (denis@lbti.org)'
  PRINTF,lun, 'Data processed on '+ SYSTIME()
  PRINTF,lun, ' '
  PRINTF,lun, '******************************************************************************************************* '
  PRINTF,lun, ' '
ENDIF ELSE lun = -1


; READ OR CREATE MASTERLOG FILE
; *****************************

; Skipped if processing nulling data right away
IF drs.skip_adi EQ 0 OR drs.skip_red EQ 0 THEN BEGIN
  
  ; --- Create masterlog file if necessary
  mlog_file = data_path + 'masterlog.dat'
  IF NOT FILE_TEST(mlog_file) OR KEYWORD_SET(MASTERLOG) THEN BEGIN
    IF FILE_TEST(data_path + 'config_id.dat') THEN FILE_DELETE, data_path + 'config_id.dat'  ; Delete config ID file to tell LBTI_MASTERLOG to create new one
    IF FILE_TEST(data_path + 'nod_id.dat')    THEN FILE_DELETE, data_path + 'nod_id.dat'     ; Delete nod ID file to tell LBTI_MASTERLOG to create new one
    PRINT, '  Creating new masterlog file. This may take a while...'
    LBTI_MASTERLOG, data_path, /IDL
  ENDIF ELSE IF FILE_TEST(mlog_file) THEN LBTI_MASTERLOG, data_path, /APPEND, /IDL
  
  ; --- Run consistancy check on the masterlog (now done separately)
  ;  LBTI_FIXDATA, data_path
  
  ; --- Read masterlog file
  LBTI_READMASTERLOG, mlog_file, log, BAD_IDX=bad_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx
     
  ; --- Flag out bad frames
  idx_bad = WHERE(log.flag EQ 'B', n_b)
  IF n_b GT 0 THEN log.bad_id[idx_bad] = 1
  
  ; --- Flag out transition frames
  IF drs.n_trans THEN BEGIN
    idx_trans = WHERE(log.nod_id NE SHIFT(log.nod_id, 1))
    FOR i_trans = 0, drs.n_trans-1 DO log.bad_id[idx_trans+i_trans] = 1
  ENDIF
  
  ; --- Flag open-loop frames to be removed (AO and phase loops where appropriate. Keep background frames, i.e. obstype = 3)
  IF drs.skip_open THEN BEGIN
    idx_open = WHERE(log.flag EQ 'O' AND log.obstype NE 3, n_o)
    IF n_o GT 0 THEN log.bad_id[idx_open] = 1
  ENDIF
  
  ; --- Plot diagnostic info
  IF plot GT 0 THEN BEGIN
    result_path = pth.result_path + pth.sep + 'diagnostic' + pth.sep + 'background' + pth.sep + date_obs + pth.sep
    IF NOT FILE_TEST(result_path) THEN FILE_MKDIR, result_path
    PLOTALL, log.file_id, log.lbtalt, 0, NAME=result_path + date + '-elevation_vs_file-id', XTITLE='Target file number', YTITLE='Telescope elevation [deg]', /BOLD, /NO_FFT, /NO_HISTO, /EPS
    PLOTALL, log.file_id, log.cv, 0, NAME=result_path + date + '-cv_vs_file-id', XTITLE='Target file number', YTITLE='Central value [ADU]', /BOLD, /NO_FFT, /NO_HISTO, /EPS
  ENDIF
ENDIF


; 1a. READ AND SAVE CALBRATION IMAGES (BPM, DARK, and FLAT)
; *********************************************************

; --- Skip calibration if requested and restore the calibrated L0 file from disk
t0 = SYSTIME(1)
IF drs.skip_red EQ 1 THEN BEGIN
  IF info GT 0 THEN PRINT, '1. Restoring calibrated L0 images from disk.'
  GOTO, skip_imgcal
ENDIF ELSE IF info GT 0 THEN PRINT, '1. Reading and calibrating L0 data.'
IF info GT 0 THEN PRINT, ' '

; --- Create master BPM, DRK, FLAT, and LIN files
IF drs.no_bpm NE 1 OR drs.no_dark NE 1 OR drs.no_flat NE 1 THEN LBTI_IMGCAL, date_obs, DARK_PATH=dark_path, FLAT_PATH=flat_path, INFO=info

; 1b. EXTRACT DATA IN RANGE AND DEFINE METRICS
; ********************************************

; --- Extract useful data and define the number of background pairs (datatype = 3 for backward compatibility)
file_id  = log.file_id
bad_id   = log.bad_id
datatype = log.datatype
IF KEYWORD_SET(DATA_IDX) THEN idx_data = WHERE(file_id GE data_idx[0] AND file_id LE data_idx[1] AND bad_id NE 1 AND ((datatype EQ 0 OR datatype EQ 3)), n_data) $
                         ELSE idx_data = WHERE(bad_id NE 1 AND ((datatype EQ 0 OR datatype EQ 3)), n_data)
IF n_data LE 0 THEN MESSAGE, 'No data found in DATA_IDX range!'

; --- Deal with BCKG_IDX
IF KEYWORD_SET(BCKG_IDX) THEN BEGIN
  idx_bckg = WHERE(file_id GE bckg_idx[0] AND file_id LE bckg_idx[1], n_bckg)
  IF n_bckg GT 0 THEN idx_data = [idx_data, idx_bckg] ELSE PRINT, 'No background frames in BCKG_IDX range!'
ENDIF               

; --- Read useful information 
file_id   = file_id[idx_data]
datatype  = datatype[idx_data]
objname   = (log.obj_name)[idx_data]
lbt_alt   = (log.lbtalt)[idx_data]
pt_id     = (log.pt_id)[idx_data]
cfg_id    = (log.cfg_id)[idx_data]
nod_id    = (log.nod_id)[idx_data]
obstype   = (log.obstype)[idx_data]
time_obs  = (log.time_obs)[idx_data]
cv        = (log.cv)[idx_data]
flag      = (log.flag)[idx_data]
n_xpix    = (log.n_xpix)[idx_data]   ; used to find the corresponding DARK AND FLAT
n_ypix    = (log.n_ypix)[idx_data]   ; used to find the corresponding DARK AND FLAT
log = {FILE_ID: file_id, NOD_ID: nod_id, CFG_ID: cfg_id, PT_ID: pt_id, CV: cv} ; Redefine log with minimum required information to save memory

; --- Compute number of different nods 
nod_uniq = nod_id[UNIQ(nod_id,  SORT(nod_id))]
idx_pos  = WHERE(nod_uniq GE 0, n_pos)
IF n_pos GT 0 THEN nod_uniq = nod_uniq[idx_pos] ELSE MESSAGE, 'No valid nods!'
n_nod    = N_ELEMENTS(nod_uniq)

; --- Derive current MJD
hr = DOUBLE(STRMID(time_obs, 0, 2)) + DOUBLE(STRMID(time_obs, 3, 2))/60D + DOUBLE(STRMID(time_obs, 6, 2))/3.6D+3 + DOUBLE(STRMID(time_obs, 9, 3))/3.6D+6
JDCNV, FIX(STRMID(date_obs, 0, 4)), FIX(STRMID(date_obs, 5, 2)), FIX(STRMID(date_obs, 8, 2)), hr, jd
mjd_obs = jd - 2400000.5D                 ; mjd

; --- Remove two NODs if data are being reduced in real-time (use 30 minutes as decision criterion)
IF (MAX(mjd_obs)+0.5/24) GT (SYSTIME(/JULIAN, /UTC) - 2400000.5D) THEN n_nod -= 2

; Print WARNING if uneven number of nods when bckg_mode = 1
IF n_nod MOD 2 NE 0 AND ABS(drs.bckg_mode) EQ 1 THEN BEGIN
  MESSAGE, "Uneven number of nods not compatible with bckg_mode = 1. Removing the last nod.", /CONTINUE
  n_nod -= 1
ENDIF

; --- Define metric used to select the background frames (compute mjd_obs if bckg_sel is not 1)
CASE drs.bckg_sel OF 
  0: metric = mjd_obs
  1: metric = lbt_alt
  2: metric = cv
  ELSE : MESSAGE, 'Undefined background selection mode (BCKG_SEL)' 
ENDCASE

; 2. READ AND REDUCE SCIENCE IMAGES
; *********************************
 
; --- Print info to screen
IF info GT 0 THEN PRINT, ' '
IF info GT 0 THEN PRINT, '   Processing ' + STRING(n_data, FORMAT='(I0)') + ' data files. This may take a while...'

; --- Loop over the nods (have to proceed per nod here to avoid reading and copying large data arrays)
t0 = SYSTIME(1)
FOR i_nod = 0, n_nod-1 DO BEGIN
  ; If file already exists, skip
  IF NOT FILE_TEST(pth.l0fits_path + date_obs + pth.sep + 'bckg' + pth.sep + '*_N' + STRING(nod_uniq[i_nod], FORMAT='(I03)') + '*IMG.fits') OR KEYWORD_SET(RENEW) THEN BEGIN
    ; Init time stamp
    tn = SYSTIME(1)
    
    ; Derive file index of frames in this nod
    idx_nod = WHERE(nod_id EQ nod_uniq[i_nod], n)
    IF n LT drs.min_frnod THEN BEGIN
      MESSAGE, 'No enough frames in nod ' + STRING(i_nod, FORMAT='(I0)'), /CONTINUE
      GOTO, skip_thisnod
    ENDIF
    
    ; Derive current pointing ID number and frame size (needed to find dark and flat frames)
    pt_cur = pt_id[idx_nod[0]]
    nx     = n_xpix[idx_nod[0]]
    ny     = n_ypix[idx_nod[0]]
  
    ; Print info to screen
    IF info GT 0 THEN PRINT, '   Reducing nod ' + STRING(nod_uniq[i_nod], FORMAT='(I0)') + ' (' + STRING(i_nod+1, FORMAT='(I0)') + '/' + $
                             STRING(n_nod, FORMAT='(I0)') + ', ' + STRING(n, FORMAT='(I0)') + ' frames)' 
    
    ; Derive obstype. If 3 (pure background), this nod can be skipped unless it is associated with a nulling nod (because it is used by the NSC)
    ; In order to decide whether to reduce this nod or not, we look for nulling frames in the curring pointing.
    obstype_nod = obstype[idx_nod[0]]                                     ; obstype of the current nod
    tmp         = WHERE(obstype[WHERE(pt_id EQ pt_cur)] EQ 2, n_null_ob)  ; n_null_ob will contain the number of null OBs in this pointing
    IF obstype_nod EQ 3 AND n_null_ob EQ 0 THEN GOTO, skip_thisnod        ; skip if background OB and no associated NULL OB in the same poiting
    
    ; Now derive file index of frames in the background nod(s). 
    ; The background nods must belong to the same pointing and have the same obstype (or an obstype of 3).
    IF obstype_nod NE 3 THEN BEGIN
      CASE ABS(drs.bckg_mode) OF
        1: idx_bck = WHERE(pt_id EQ pt_cur AND nod_id EQ nod_uniq[i_nod+(-1)^i_nod] AND ((obstype_nod EQ obstype OR obstype EQ 3)), n_bck)                            ; Background subtraction by nod pairs
        2: idx_bck = WHERE(pt_id EQ pt_cur AND ((nod_id EQ nod_uniq[i_nod]-1 OR nod_id EQ nod_uniq[i_nod]+1)) AND ((obstype_nod EQ obstype OR obstype EQ 3)), n_bck)  ; Background subtraction using adjacent nods
        3: idx_bck = WHERE(pt_id EQ pt_cur AND obstype EQ 3, n_bck)                                                                                                   ; Background subtraction using dedicated background frames
        ELSE: BEGIN
           idx_bck = WHERE(pt_id EQ pt_cur AND nod_id EQ nod_uniq[i_nod]+(-1)^i_nod, n_bck)
           IF n_bck LE 0 THEN idx_bck = WHERE(pt_id EQ pt_cur AND nod_id EQ nod_uniq[i_nod]-(-1)^i_nod, n_bck)
        END
      ENDCASE
      ; If no adjacent background frames, then look within pointing
      IF n_bck LE 0 THEN idx_bck = WHERE(pt_id EQ pt_cur AND obstype EQ 3, n_bck)
      ; If still no background frames, then look everywhere
      IF n_bck LE 0 THEN idx_bck = WHERE(obstype EQ 3, n_bck) 
    ENDIF ELSE idx_bck = WHERE(pt_id EQ pt_cur AND ((nod_id EQ nod_uniq[i_nod]-1 OR nod_id EQ nod_uniq[i_nod]+1)), n_bck)    
    
    ; Remove outliers based on central value (5 sigma RMS). Not necessary
    ; cv_nod  = cv[idx_nod]
    ; cv_bck  = cv[idx_bck]
    ; AVGSDV, [cv_bck, cv_nod], cv_avg, cv_rms
    ; idx_bad = WHERE(ABS(cv_nod-cv_avg) GT 5*cv_rms, COMPLEMENT=idx_ok, n_bad) & IF n_bad GT 0 THEN idx_nod = idx_nod[idx_ok] & n     = N_ELEMENTS(idx_nod)
    ; idx_bad = WHERE(ABS(cv_bck-cv_avg) GT 5*cv_rms, COMPLEMENT=idx_ok, n_bad) & IF n_bad GT 0 THEN idx_bck = idx_bck[idx_ok] & n_bck = N_ELEMENTS(idx_bck)
    
    ; Keep only background frames within maximum allowed time
    IF drs.max_time GT 0 THEN BEGIN
      idx_time = WHERE(ABS(mjd_obs[idx_bck]-mjd_obs[idx_nod[0]]) LT drs.max_time/24./60. OR ABS(mjd_obs[idx_bck]-mjd_obs[idx_nod[-1]]) LT drs.max_time/24./60., n_bck)
      IF n_bck GT drs.min_frnod THEN idx_bck = idx_bck[idx_time]
    ENDIF
    
    ; Keep only the closest frames based on the chosen metric (e.g., time, elevation, central value)
    IF drs.n_frbck NE 0 AND drs.n_frbck LT n_bck THEN BEGIN
      IF drs.n_frbck  EQ -1 THEN n_bck = (n<n_bck) ELSE n_bck = drs.n_frbck         ; if -1, keep the same number of frames as in the current nod
      idx_min = SORT(ABS(metric[idx_bck]-MEAN(metric[idx_nod])))                    ; we compare to the mean value of the current NOD
      idx_bck = idx_bck[idx_min[0:n_bck-1]]
    ENDIF   
      
    ; Skip if no background frames   
    IF n_bck LT drs.min_frnod THEN BEGIN
      MESSAGE, 'No background frames for nod ' + STRING(i_nod, FORMAT='(I0)'), /CONTINUE
      GOTO, skip_thisnod
    ENDIF
    
    ; Concatenate file ID to be passed to LBTI_READDATA (and sort them)
    ;PRINT, MIN(idx_nod), MAX(idx_nod), N_ELEMENTS(idx_nod)
    ;PRINT, MIN(idx_bck), MAX(idx_bck), N_ELEMENTS(idx_bck)
    idx_all = [idx_nod,idx_bck]
    fid_all = file_id[idx_all]
    fid_all = fid_all[SORT(fid_all)]
    min_id  = MIN(fid_all)
    max_id  = MAX(fid_all)
    
    ; Now look if some files have already been read in the previous iteration (and remove them from the files to be read in this iteration)
    ; FILE ID must match, not file index!
    IF KEYWORD_SET(hdr_data) THEN BEGIN
      MATCH, fid_all, hdr_data.file_id, idx_ok, idx_prev, COUNT = n_prev  ; hdr_data.file_id is defined below and comes from the previous iteration
      IF n_prev GT 0 THEN BEGIN
        ; File to be read
        fid_all[idx_ok] = -1
        idx_read = WHERE(fid_all NE -1, n_read)
        ;IF info GE 3 THEN PRINT, '    Frames already read in previous iteration', n_prev, n_read
        IF n_read GT 0 THEN BEGIN
          fid_all = fid_all[idx_read]
          ; Retrieve frames already read 
          ; img_data and hdr_data are defined below and come from the previous iteration
          img_prev = img_data[*,*,idx_prev]
          hdr_prev = hdr_data[idx_prev]
          ; Undefine arrays to save memory
          UNDEFINE, img_data, hdr_data
        ENDIF         
      ENDIF ELSE n_read = 1
    ENDIF ELSE n_read = 1
    
    ; Read OB's data (only if new frames to read)
    IF n_read GT 0 THEN BEGIN
      IF info GT 0 THEN PRINT, '    - reading L0 data.'
      img_data = LBTI_READDATA(data_path, CROP=drs.pre_crop, DATA_IDX=fid_all, HDR_DATA=hdr_data, MLOG_DATA=log, INFO=info, PLOT=plot) 
                                                     
      ; Reduce L0 images (dark subtraction, flat fielding, wavelentgh calibration, pickup noise, etc...)
      IF info GT 0 THEN PRINT, '    - calibrating L0 images.'
      
      ; --- Read last dark and flat files (assumes that we used only a single size throught the night)
      IF drs.no_bpm  NE 1 THEN LBTI_READCALIMG, pth.bpm_path,  date_obs, [nx,ny], -999, bpm, hdr_bpm, CROP=drs.pre_crop ELSE IF info GT 0 THEN PRINT, '      - no bad pixel correction'  ; -999 will tell this routine to use the master bad pixel map
      IF drs.no_bpm  NE 1 AND (SIZE(bpm))[0] LE 1 THEN GOTO, skip_thisnod
      IF drs.no_dark NE 1 THEN LBTI_READCALIMG, pth.dark_path, date_obs, [nx,ny], pt_cur, drk, hdr_drk, CROP=drs.pre_crop ELSE IF info GT 0 THEN PRINT, '      - no dark subraction'
      IF drs.no_dark NE 1 AND (SIZE(drk))[0] LE 1 THEN GOTO, skip_thisnod
      IF drs.no_flat NE 1 THEN LBTI_READCALIMG, pth.flat_path, date_obs, [nx,ny], pt_cur, flt, hdr_flt, CROP=drs.pre_crop ELSE IF info GT 0 THEN PRINT, '      - data not flat fielded'
      IF drs.no_flat NE 1 AND (SIZE(flt))[0] LE 1 THEN GOTO, skip_thisnod
            
      ; --- Flat field each frame of this nod
      img_data = LBTI_IMGRED(TEMPORARY(img_data), TEMPORARY(hdr_data), IMG_BPM=bpm, IMG_DRK=drk, IMG_FLT=flt, HDR_DRK=hdr_drk, HDR_FLT=hdr_flt, $             ; Images and correspondind header        
                             HDR_DATA=hdr_data, $                                                                                                             ; Output keywords
                             LOG_FILE=log_file, INFO=info, PLOT=plot)                                                                                         ; Input keywords
      IF (SIZE(img_data))[0] EQ 0 THEN GOTO, skip_thisnod
      
      ; Undefine arrays to save memory
      UNDEFINE, bpm, drk, flt, hdr_bpm, hdr_drk, hdr_flt
      
      ; Now concatenate files read this iteration with files already read before
      IF KEYWORD_SET(hdr_prev) GT 0 THEN BEGIN
        IF n_prev GT 0 THEN BEGIN ; n_prev not defined if i_nod = 0
          ; Concatenate
          img_data = [[[img_data]],[[img_prev]]]
          hdr_data = [hdr_data,hdr_prev]
          ; Sort (better have data in chronolical order for later)
          idx_srt  = SORT(hdr_data.file_id)
          hdr_data = TEMPORARY(hdr_data[idx_srt])
          img_data = TEMPORARY(img_data[*,*,idx_srt])
        ENDIF
      ENDIF
    ENDIF
    
    ; Only keep data in the useful range
    idx_ok   = WHERE(hdr_data.file_id GE min_id AND hdr_data.file_id LE max_id)
    img_data = img_data[*,*,idx_ok]
    hdr_data = hdr_data[idx_ok]
    
    ; Assign background nod flag
    idx_bck = WHERE(hdr_data.nod_id NE nod_uniq[i_nod], COMPLEMENT=idx_nod)
    hdr_data[idx_nod].bck_nod = 0
    hdr_data[idx_bck].bck_nod = 1
      
    ; Perform background subtraction and save data         
    IF info GT 0 THEN PRINT, '    - performing background subtraction.'
    LBTI_IMGBCK, img_data, hdr_data, BCKG_PATH=bckg_path, LOG_FILE=log_file, INFO=info, NO_SAVE=no_save, PLOT=plot                               ; Input keywords
                                                                             
    ; Jump point of problem with nod computation
    skip_thisnod:
    
    ; Print time to process this nod
    IF info GE 3 THEN PRINT, '   Time to process this nod : ', SYSTIME(1)-tn
  ENDIF ELSE PRINT, '   Nod ' + STRING(nod_uniq[i_nod], FORMAT='(I0)') + ' already exists'
ENDFOR

; Print info to screen
t1 = SYSTIME(1)
IF info GT 2 THEN PRINT, 'Time to perform L0 data calibration :', t1-t0
IF info GT 0 THEN PRINT, ' '

; Jumping point if restoring calibrated L0 files
skip_imgcal:


; 3. FLUX COMPUTATION (also needed for frame centering and cropping before ADI processing)
; *******************

; --- Skip flux computation if requested and restore flux from disk
t2 = SYSTIME(1)
IF drs.skip_flx EQ 1 AND drs.skip_vis EQ 1 THEN BEGIN
  IF info GT 0 THEN BEGIN
    IF drs.skip_flx EQ 1 THEN PRINT, '2. Restoring flux measurements from disk.' $
                         ELSE IF drs.skip_vis EQ 1 THEN PRINT, '2. Restoring visibility measurements from disk.'
  ENDIF
  GOTO, skip_flx
ENDIF
IF info GT 0 THEN PRINT, ' '

; --- Restore LO log file and read useful data
datalog = pth.l0fits_path + date_obs + pth.sep + 'datalog.sav'
IF NOT FILE_TEST(datalog) THEN BEGIN
  MESSAGE, 'No data log file found!', /CONTINUE
  RETURN
ENDIF ELSE RESTORE, datalog

; --- Extract relevant information for all data (used to computed OB id of photometric and background frames)
objname = data_r.objname
nod_id  = data_r.nod_id   & nod_id0 = nod_id
cfg_id  = data_r.cfg_id   & cfg_id0 = cfg_id
pt_id   = data_r.pt_id    & pt_id0  = pt_id
obstype = data_r.obstype
xcen_sx = data_r.xcen_sx  & xcen_dx = data_r.xcen_dx
ycen_sx = data_r.ycen_sx  & ycen_dx = data_r.ycen_dx
min_fid = data_r.min_fid
max_fid = data_r.max_fid
flag    = data_r.flag

time_obs = data_r.ut_time
ut_time  = DOUBLE(STRMID(time_obs, 0, 2)) + DOUBLE(STRMID(time_obs, 3, 2))/60D + DOUBLE(STRMID(time_obs, 6, 2))/3.6D+3 + DOUBLE(STRMID(time_obs, 9, 3))/3.6D+6

; --- Read target information (or query online catalogs if not present)
tgt_uniq = objname[UNIQ(objname,  SORT(objname))]
GET_TGT, tgt_uniq, tgt, DATABASE= pth.input_path + drs.database
struct_add_field, tgt, 'calfor', 'tbd'
struct_add_field, tgt, 'maxapr', 1D3

; --- Assign calibrators to their science object (needed because of EEID computation).
; --- Also assign the maximum aperture radius possible for each star.
pt_uniq = pt_id[UNIQ(pt_id,  SORT(pt_id))]
n_pt    = N_ELEMENTS(pt_uniq)
flag_pt = STRARR(n_pt)
tgt_pt  = flag_pt
max_apr = INTARR(n_pt)
FOR i = 0, n_pt-1 DO BEGIN
  idx        = WHERE(pt_id EQ pt_uniq[i])
  flag_pt[i] = flag[idx[0]]
  tgt_pt[i]  = STRLOWCASE(STRCOMPRESS(objname[idx[0]], /REMOVE_ALL)) 
   ; in sub-frame mode, this doesn't capture the distance to the horizontal center but no easy solution for that at this point...  
   ; this is also only valid for SX, which encodes the position of the nulling beam. PHOTOMETRIC beam should be close
  idx        = WHERE(pt_id EQ pt_uniq[i] AND ycen_sx NE 0 AND xcen_sx NE 0)
  IF drs.sky_col EQ 0 THEN  max_apr[i] = MIN(xcen_sx[idx] MOD cnf.X_CHAN) < MIN(ycen_sx[idx] MOD cnf.Y_CHAN) < MIN(cnf.X_CHAN - (xcen_sx[idx] MOD cnf.X_CHAN)) < MIN(cnf.Y_CHAN - (ycen_sx[idx] MOD cnf.Y_CHAN)) $
                      ELSE  max_apr[i] = MIN(ycen_sx[idx] MOD cnf.Y_CHAN) < MIN(cnf.Y_CHAN - (ycen_sx[idx] MOD cnf.Y_CHAN))  ; look only in the Y direction!
  IF flag_pt[i] EQ 'SCI' THEN BEGIN
    idx_tgt = WHERE(tgt.name EQ tgt_pt[i])
    IF tgt[idx_tgt].maxapr GT max_apr[i] THEN tgt[idx_tgt].maxapr = max_apr[i]
  ENDIF
ENDFOR 
; Pre-defined list (should go in a separate file eventually)
READ_TABLE, 'nodrs/input/calsci_pairs.txt', date_list, ut_min, ut_max, cal_list, sci_list, FIRST=2, STRING=[3,4], SEPARATOR=';'

; Loop over the pointings
FOR i = 0, n_pt-1 DO BEGIN
  ; Current target
  idx_tgt = WHERE(tgt.name EQ tgt_pt[i])
  ra_tgt  = tgt[idx_tgt].ra
  de_tgt  = tgt[idx_tgt].dec
  ; If in pre-defined list, use it!
  idx_ok = WHERE(STRLOWCASE(tgt_pt[i]) EQ STRLOWCASE(cal_list) AND date EQ date_list AND ut_time[i] GE ut_min AND ut_time[i] LE ut_max, n_ok)
  IF n_ok EQ 1 THEN tgt[idx_tgt].calfor = STRLOWCASE(sci_list[idx_ok]) ELSE BEGIN
    ; If CAL, find nearest SCIs
    IF flag_pt[i] EQ 'CAL' THEN BEGIN
      idx_low = MAX(WHERE(pt_uniq LT pt_uniq[i] AND flag_pt EQ 'SCI', n1)) & IF n1 GT 0 THEN tgt_low = tgt_pt[idx_low] ELSE tgt_low = 'none'  ; Nearest previous SCI pointings
      idx_up  = MIN(WHERE(pt_uniq GT pt_uniq[i] AND flag_pt EQ 'SCI', n2)) & IF n2 GT 0 THEN tgt_up  = tgt_pt[idx_up]  ELSE tgt_up  = 'none'  ; Nearest following SCI pointings
      ; Current target
      IF n1 NE 0 OR n2 NE 0 THEN BEGIN        
        ; If not the same, find the one with the closest RA
        IF tgt_low NE tgt_up THEN BEGIN
          idx = WHERE(tgt.name EQ tgt_low, n_low)
          IF n_low GE 1 THEN ra_low  = tgt[idx].ra  ELSE ra_low = 1D9
          IF n_low GE 1 THEN de_low  = tgt[idx].dec ELSE de_low = 1D9
          idx = WHERE(tgt.name EQ tgt_up, n_up)
          IF n_up GE 1 THEN ra_up = tgt[idx].ra  ELSE ra_up = 1D9
          IF n_up GE 1 THEN de_up = tgt[idx].dec ELSE de_up = 1D9
          IF ((ra_tgt-ra_low)^2+(de_tgt-de_low)^2) LT ((ra_tgt-ra_up)^2+(de_tgt-de_up)^2) THEN tgt[idx_tgt].calfor = tgt_low ELSE tgt[idx_tgt].calfor = tgt_up
        ENDIF ELSE tgt[idx_tgt].calfor = tgt_up
      ENDIF ELSE BEGIN
        ; Only needed if aper_rad is 0
        IF drs.aper_rad EQ 0 THEN BEGIN
          READ, ttt, PROMPT='Enter science object for ' + tgt_pt[i] + ' : '
          tgt[idx_tgt].calfor = ttt
        ENDIF ELSE tgt[idx_tgt].calfor = 'UND'
      ENDELSE
    ENDIF
   ENDELSE
ENDFOR

; --- Consitency check
; Make sure there is no missing object name
;idx_missing = WHERE(obj_name EQ '?' OR obj_name EQ ' ', n_missing)
;IF n_missing GE 1 THEN MESSAGE, 'Problem with object name in files : ' + STRING(file_id[idx_missing], FORMAT='(I0)') + ' (' + obj_name[idx_missing] +')'

; --- Assign unique ID number to each line of the file
n_data  = N_ELEMENTS(objname)
file_id = INDGEN(n_data) + 1  ; I don't want 0!!!!

; --- Assign OB id number to coherent data (so don't consider photometric and background frames in the OB id count)
idx_coh = WHERE(obstype EQ 1 OR obstype EQ 2, n_coh)
IF n_coh GT 0 THEN BEGIN
  ob_id = INTARR(n_data) - 1 ; Photometry and background frames have an OB id of -1.
  CASE drs.ob_mode OF    
    0: dat_id = [TRANSPOSE(nod_id[idx_coh]),TRANSPOSE(pt_id[idx_coh]),TRANSPOSE(cfg_id[idx_coh])]
    1: MESSAGE, 'OB_MODE=1 not yet implemented'
    2: MESSAGE, 'OB_MODE=2 not yet implemented'
    3: MESSAGE, 'OB_MODE=3 not yet implemented (one OB per chop position)'
  ENDCASE
  ob_id[idx_coh] = ATTRIBUTE_ID(dat_id) + 1 ; I don't want 0!!!!
ENDIF ELSE ob_id = nod_id
ob_id0 = ob_id

; --- Assign background and photometric file ID for each coherent OB (only needed for null OB in principle)
bck_id  = INTARR(n_data)           ; will contain the corresponding background file ID
pho_fid = INTARR(n_data)           ; will contain the corresponding photometric file ID
nod_pos = FIX(ycen_sx/cnf.y_chan)  ; this assumes that each NOD ends up in a different vertical channel (which is true for nulling)
FOR i_d = 0, n_data-1 DO BEGIN
  
  ; 1. Compute the corredponding background OB. It must be of the same pointing (and same config ID) but different nod position (+ avoid photometric frames)
  idx_bck = WHERE(pt_id EQ pt_id[i_d] AND cfg_id EQ cfg_id[i_d] AND ob_id NE ob_id[i_d] AND obstype NE 0 AND (nod_pos NE nod_pos[i_d] OR obstype EQ 3), n_bck)
  IF n_bck GT 0 THEN BEGIN
    file_bck = file_id[idx_bck]                                                     ; list of suitable files
    idx_min  = WHERE(ABS(file_id[i_d]-file_bck) EQ MIN(ABS(file_id[i_d]-file_bck))) ; closest ones
    IF N_ELEMENTS(idx_min) GT 1 THEN BEGIN                                          ; if two, keep the one with more frames
      tmp = MAX(max_fid[idx_bck[idx_min]]-min_fid[idx_bck[idx_min]], i_max)
    ENDIF ELSE i_max = 0
    bck_id[i_d] = file_bck[idx_min[i_max]]              ; retain only one
  ENDIF ELSE bck_id[i_d] = -1                           ; no associated background file
  
  ; 2. Compute the associated photometric file ID. Ideally, it must be of the same pointing but look in different pointings if not found
  idx_pho = WHERE(pt_id EQ pt_id[i_d] AND cfg_id EQ cfg_id[i_d] AND obstype EQ 0 AND xcen_sx NE 0 AND ycen_sx NE 0, n_pho)
  IF n_pho LE 0 THEN BEGIN
    PRINT, ' No corresponding photometric file for nod ' + STRING(nod_id[i_d], FORMAT='(I0)') + '. Looking in other pointings.'
    IF lun GT 0 THEN PRINTF, lun, ' No corresponding photometric file for nod ' + STRING(nod_id[i_d], FORMAT='(I0)') + '. Looking in other pointings.'
    idx_pho = WHERE(objname EQ objname[i_d] AND cfg_id EQ cfg_id[i_d] AND obstype EQ 0 AND xcen_sx NE 0 AND ycen_sx NE 0, n_pho)
  ENDIF 
  IF n_pho GT 0 THEN BEGIN
    fid_pho = file_id[idx_pho]                                                      ; list of suitable file IDs
    idx_min  = WHERE(ABS(file_id[i_d]-fid_pho) EQ MIN(ABS(file_id[i_d]-fid_pho))) ; closest one
    IF N_ELEMENTS(idx_min) GT 1 THEN BEGIN                                          ; if two, keep the one with the more frames
      tmp = MAX(max_fid[idx_pho[idx_min]]-min_fid[idx_pho[idx_min]], i_max)
    ENDIF ELSE i_max = 0
    pho_fid[i_d] = fid_pho[idx_min[i_max]]            ; retain only one
  ENDIF ELSE BEGIN
    IF info GT 0 THEN BEGIN
      PRINT, ' No corresponding photometric file for nod ' + STRING(nod_id[i_d], FORMAT='(I0)')
      IF lun GT 0 THEN PRINTF, lun, ' No corresponding photometric file for nod ' + STRING(nod_id[i_d], FORMAT='(I0)')
    ENDIF
    pho_fid[i_d] = -1 ; no associated photometric NOD
  ENDELSE  
ENDFOR

; -- Extract data in DATA_IDX range
IF KEYWORD_SET(DATA_IDX) THEN BEGIN
  idx_data = WHERE(data_r.min_fid GE data_idx[0] AND data_r.max_fid LE data_idx[1], n_data)
  IF n_data GT 0 THEN BEGIN
    data_r  = data_r[idx_data] 
    cfg_id0 = cfg_id[idx_data] ; only used for print
    pt_id0  = pt_id[idx_data]  ; only used for print
    nod_id0 = nod_id[idx_data] ; only used for print
    ob_id0  = ob_id[idx_data]  ; only used for print
  ENDIF ELSE BEGIN
    MESSAGE, 'No data in the DATA_IDX range', /CONTINUE
    RETURN
  ENDELSE
ENDIF ELSE n_data = N_ELEMENTS(obj_name)

; -- Derive number of unique nod IDs
nod_uniq = nod_id0[UNIQ(nod_id0,  SORT(nod_id0))]  
n_nod    = N_ELEMENTS(nod_uniq)  > 1 

; -- Print info to screen and to the log
IF info GT 0 THEN BEGIN
  ; Extract data in DATA_IDX range
  wav_eff = data_r.lam_cen
  dit     = data_r.int_time
  tgtname = data_r.objname
  ; Print info to screen
  IF info GT 2 THEN BEGIN
    ; Compute the number of different parameters in the input data sequence
    tgt_uniq  = tgtname[UNIQ(tgtname,  SORT(tgtname))]  & n_tgt  = N_ELEMENTS(tgt_uniq)  > 1  ; Distinct object name
    lam_uniq  = wav_eff[UNIQ(wav_eff, SORT(wav_eff))]   & n_lam  = N_ELEMENTS(lam_uniq)  > 1  ; Distinct integration time
    int_uniq  = dit[UNIQ(dit, SORT(dit))]               & n_int  = N_ELEMENTS(int_uniq)  > 1  ; Distinct central wavelengths
    pt_uniq   = pt_id0[UNIQ(pt_id0,  SORT(pt_id0))]     & n_pt   = N_ELEMENTS(pt_uniq)   > 1  ; Distinct pointing #ID
    ob_uniq   = ob_id0[UNIQ(ob_id0,  SORT(ob_id0))]     & n_ob   = N_ELEMENTS(ob_uniq)   > 1  ; Distinct OB #ID
    ; Print info to screen
    PRINT, '2. Performing flux/visibility computation over the following parameters:'
    PRINT, ' - date                    : ', STRTRIM(date_obs)
    PRINT, ' - target name             : ', STRTRIM(tgt_uniq)
    PRINT, ' - central wavelength [um] : ', STRING(1D+6*lam_uniq, FORMAT="(F5.2)")
    PRINT, ' - integration time [ms]   : ', STRING(1D+3*int_uniq, FORMAT="(I0)")
    PRINT, ' - number of pointings     : ', STRING(n_pt, FORMAT="(I0)")
    PRINT, ' - number of nods          : ', STRING(n_nod, FORMAT="(I0)")
    PRINT, ' - number of OBs           : ', STRING(n_ob, FORMAT="(I0)")
    PRINT, ' '
    ; Print also the version number in the log file
    IF lun GT 0 THEN BEGIN
      PRINTF, lun, '2. Performing flux/visibility computation over the following parameters:'
      PRINTF, lun, ' - date                    : ', STRTRIM(date_obs)
      PRINTF, lun, ' - target name             : ', STRTRIM(tgt_uniq)
      PRINTF, lun, ' - central wavelength [um] : ', STRING(1D+6*lam_uniq, FORMAT="(F3.1)")
      PRINTF, lun, ' - integration time [ms]   : ', STRING(1D+3*int_uniq, FORMAT="(I0)")
      PRINTF, lun, ' - number of pointings     : ', STRING(n_pt, FORMAT="(I0)")
      PRINTF, lun, ' - number of nods          : ', STRING(n_nod, FORMAT="(I0)")
      PRINTF, lun, ' - number of OBs           : ', STRING(n_ob, FORMAT="(I0)")
      PRINTF, lun, ' '
    ENDIF
  ENDIF
  
  IF drs.flx_mode NE 2 THEN BEGIN
    PRINT, 'Now computing flux by APERTURE PHOTOMETRY on ' + STRING(n_data, FORMAT='(I0)') + ' frames.'
    PRINT,'      - Photometric aperture radius [pix] : ' + STRING(drs.aper_rad, FORMAT='(I0)') + ' (0: 2xEEID)'
    PRINT,'      - Inner background radius [pix]     : ' + STRING(drs.bck_irad, FORMAT='(I0)') + ' (0: aperture radius)'
    PRINT,'      - Outer background radius [pix]     : ' + STRING(drs.bck_orad, FORMAT='(I0)') + ' (0: channel edge)'
  ENDIF ELSE PRINT, 'Now computing flux by PSF-FITTING on ' + STRING(n_data, FORMAT='(I0)') + ' frames.'

  ; Initiate pointing number display
  PRINT, ' '
  PRINT, FORMAT='("Processing nod number ", $)'
ENDIF
IF lun GT 0 THEN BEGIN
  IF drs.flx_mode NE 2 THEN BEGIN
    PRINTF,lun,' '
    PRINTF,lun,'Flux computed by APERTURE PHOTOMETRY on ' + STRING(n_data, FORMAT='(I0)') + ' frames:'
    PRINTF,lun,'      - Photometric aperture radius [pix] : ' + STRING(drs.aper_rad, FORMAT='(I0)') + ' (0: 2xEEID)'
    PRINTF,lun,'      - Inner background radius [pix]     : ' + STRING(drs.bck_irad, FORMAT='(I0)') + ' (0: aperture radius)'
    PRINTF,lun,'      - Outer background radius [pix]     : ' + STRING(drs.bck_orad, FORMAT='(I0)') + ' (0: channel edge)'
    PRINTF,lun,' '
  ENDIF ELSE BEGIN
    PRINTF,lun,' '
    PRINTF,lun,'Flux computed by PSF-FITTING on ' + STRING(n_data, FORMAT='(I0)') + 'frames.'
    PRINTF,lun,' '
  ENDELSE
ENDIF

; --- If aperture radius is set, save L1 files in a separate directory (no extra label if EEID)
IF MAX(drs.aper_rad) NE 0 THEN drs.dir_label = drs.dir_label + '_APR' ;+ STRING(MIN(drs.aper_rad[WHERE(drs.aper_rad GT 0)]), FORMAT='(I0)')

; --- Perform flux computation (aperture photometry or PSF fitting)
FOR i_nod = 0, n_nod-1 DO BEGIN
  
  ; Print info to screen
  IF info GT 0 THEN PRINT, nod_uniq[i_nod], format='($, (x, I0))'
  
  ; Find most recent files of this nod
  idx_nod  = WHERE(nod_id EQ nod_uniq[i_nod], n_cfg)
  IF n_cfg LE 0 THEN GOTO, skip_nod
  
  ; Loop over the files
  FOR i_cfg = 0, n_cfg-1 DO BEGIN
    ; Derive file ID
    idx_cfg = idx_nod[i_cfg]
    fid_cur = file_id[idx_cfg]
    ob_cur  = ob_id[idx_cfg]
    obstype_cur = obstype[idx_cfg]
    
    ; Check whether this OB has already been reduced
    IF obstype_cur EQ 0 THEN tag = 'PHOT'
    IF obstype_cur EQ 1 THEN tag = 'FIZ'
    IF obstype_cur EQ 2 THEN tag = 'NULL'
    IF obstype_cur GE 3 THEN tag = 'DOANYWAY'  ; IF OBSTYPE GE 3, the flux has to be recomputed. I don't remember why but it doesn't hurt to repeat it..
    
    ; Skip if already exists
    IF NOT FILE_TEST(pth.l1fits_path + drs.date_obs + drs.dir_label + pth.sep + '*_ID' + STRING(ob_cur, FORMAT='(I03)') + '*' + tag + '*.fits') OR KEYWORD_SET(RENEW) THEN BEGIN
      
      ; Read raw L0 file (use mean-subtracted, pca-subtracted or raw-subtracted frames)
      CASE drs.fra_mode OF 
        0: label = 'bckg'
        1: label = 'raw'
        2: label = 'pca'
        ELSE : MESSAGE, 'Undefined frame selection mode (FRA_MODE)' 
      ENDCASE
      IF obstype_cur EQ 0 THEN label = 'bckg' ; 0CT 2023, updated for PCA frames. I don't remember why this bit is here. I leave it here for backward compatibility? 


      file_nod = FILE_SEARCH(pth.l0fits_path + date_obs + pth.sep + label + pth.sep + '*_N' + STRING(nod_uniq[i_nod], FORMAT='(I03)') + '*' + STRING(min_fid[idx_cfg], FORMAT='(I06)') + '-' +  $
                             STRING(max_fid[idx_cfg], FORMAT='(I06)') + '_IMG.fits', COUNT=n0)
      IF n0 GT 0 THEN BEGIN
        
        ; If more than 1 file, only keep the latest file and delete the other one
        IF n0 GT 1 THEN BEGIN
          finfo    = FILE_INFO(file_nod)
          idx      = WHERE(finfo.mtime EQ MAX(finfo.mtime), COMPLEMENT=idxd)
          FILE_DELETE, file_nod[idxd]
          file_nod = file_nod[idx]
        ENDIF
        
        ; Read image
        img_data = LBTI_READL0RED(file_nod, HDR_DATA=hdr_data, INFO=info)
        
        ; Current pointing, nod, OB and file ID
        pt_cur  = pt_id[idx_cfg]
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
        CASE hdr_data.header.obstype OF 
          0: BEGIN ; Photometric frames (if it exists, take OB ID of associated coherent frames, else take its own asociated ob_id)
               idx_pho  = WHERE(pho_fid EQ fid_cur AND ob_id NE -1, n_ok)
               IF n_ok GT 0 THEN ob_in = ob_id[idx_pho] $      ; contain all OBs that use this photometric file
                            ELSE ob_in = -1                    ; not used, -1 to skip flux computation below
             END
          1: BEGIN ; Visibility frames
               hdr_data.data[*].ob_id = ob_cur
               idx_bck = WHERE(bck_id EQ fid_cur AND ob_id NE -1, n_ok)
               IF n_ok GT 0 THEN BEGIN
                 ob_in    = ob_id[idx_bck]   ; contain all OBs that use the current OB as background
                 xcen_bck = xcen_sx[idx_bck] ; contain corresponding X positions
                 ycen_bck = ycen_sx[idx_bck] ; contain corresponding Y positions
               ENDIF ELSE ob_in = 0          ; not used
             END
          2: BEGIN ; Nulling frames 
               hdr_data.data[*].ob_id = ob_cur               
               idx_bck = WHERE(pt_id EQ pt_cur AND cfg_id EQ cfg_cur AND ob_id NE ob_cur AND obstype NE 0 AND nod_pos NE nod_pos[idx_cfg], n_bck)  ; will return the file ID of the opposite NODs
               IF n_bck GE 1 THEN BEGIN
                 nod_bck  = nod_id[idx_bck]                                                     
                 idx_bck  = idx_bck[WHERE(ABS(nod_cur-nod_bck) EQ MIN(ABS(nod_cur-nod_bck)), n_ok)]  ; keep only the closest one(s)
               ENDIF ELSE idx_bck = WHERE(bck_id EQ fid_cur AND ob_id NE -1, n_ok) ; old approach               
               IF n_ok GT 0 THEN BEGIN
                 ob_in    = ob_id[idx_bck]    ; contain all OBs that use the current OB as background
                 xcen_bck = xcen_sx[idx_bck]  ; contain corresponding X positions
                 ycen_bck = ycen_sx[idx_bck]  ; contain corresponding Y positions
               ENDIF ELSE ob_in = 0           ; not used
             END
          3: BEGIN ; Background frames (only used if assciated to NULL nod)
               idx_bck = WHERE(bck_id EQ fid_cur AND ob_id NE -1, n_ok)
               IF n_ok GT 0 THEN BEGIN
                 ob_in    = ob_id[idx_bck]     ; contain all OBs that use the current OB as background
                 xcen_bck = xcen_sx[idx_bck]   ; contain corresponding X positions
                 ycen_bck = ycen_sx[idx_bck]   ; contain corresponding X positions
               ENDIF ELSE ob_in = -1           ; not used, -1 to skip flux computation below
             END
           ELSE: GOTO, skip_nod
        ENDCASE
        IF info GT 2 THEN PRINT, ' - PT, NOD, OB | BCK_OB:', pt_cur, nod_cur, ob_cur, ' | ', ob_in
        
        ; Compute flux/visibility for all frames in this nod and save data if allowed (don't do background frames not used in conjunction with null data)
        IF NOT drs.skip_flx AND MAX(ob_in) NE -1 THEN LBTI_IMG2FLX, img_data, hdr_data, LOG_FILE=lun, INFO=info, PLOT=plot, NO_SAVE=no_save, OB_IN=ob_in, XCEN=xcen_bck, YCEN=ycen_bck
        IF NOT drs.skip_vis AND MAX(ob_in) NE -1 THEN LBTI_IMG2VIS, img_data, hdr_data, LOG_FILE=lun, INFO=info, PLOT=plot, NO_SAVE=no_save, OB_IN=ob_in, XCEN=xcen_bck, YCEN=ycen_bck
      ENDIF
      
    ENDIF ;ELSE PRINT, ' File ' + STRING(fid_cur, FORMAT='(I0)') + '  already exists'
  ENDFOR
  ; Jumping point if no file for this nod ID
  skip_nod:
ENDFOR

t3 = SYSTIME(1)
IF info GT 0 THEN  PRINT, ' '
IF info GT 2 THEN  PRINT, 'Time to perform flux computation :', t3-t2 
  
; Jumping point if SKIP_FLX is set
skip_flx: 

 ; 4. ADI PROCESSING (from processed L1 frames)
; ******************

; --- Skip ADI processing if requested
t1 = SYSTIME(1)
IF drs.skip_adi THEN BEGIN
  IF info GT 0 THEN PRINT, '3. Skipping ADI processing'
  GOTO, skip_adi
ENDIF

; --- Print info to screen
IF info GT 0 THEN BEGIN
  PRINT, ' '
  PRINT, '3. Performing ADI processing:'
  PRINT, FORMAT='("Processing nod number ", $)'
ENDIF

; --- Restore LO log file and read useful data
datalog = pth.l0fits_path + date_obs + pth.sep + 'datalog.sav'
IF FILE_TEST(datalog) THEN RESTORE, datalog ELSE MESSAGE, 'No data log file found!'
objname  = data_r.objname
lam_cen  = data_r.lam_cen
obj_uniq = STRCOMPRESS(objname[UNIQ(objname,  SORT(objname))], /REMOVE_ALL)
n_obj    = N_ELEMENTS(obj_uniq)  > 1  ; Distinct object name

; --- Loop over the objects
FOR i_o = 0, n_obj-1 DO BEGIN
  ; --- Derive number of wavelengths (need to be improved)
  tgt_name = obj_uniq[i_o]
  PRINT,  tgt_name
  idx_obj  = WHERE(objname EQ tgt_name)
  lam_cur  = lam_cen[idx_obj]
  lam_uniq = lam_cur[UNIQ(lam_cur,  SORT(lam_cur))]
  n_lam    = N_ELEMENTS(lam_uniq)
  ; --- Loop over wavelengths
  FOR i_lam = 0, n_lam-1 DO BEGIN  
    ; --- Merge all files 
    search_string = '*' + tgt_name + '*' + STRING(1D+6*lam_uniq[i_lam], FORMAT='(I0)') 
    file = FILE_SEARCH(pth.l1fits_path + date_obs + drs.dir_label + pth.sep, search_string + '*_IMG.fits', COUNT=n0)
    IF n0 LT 1 THEN GOTO, skip_file ELSE IF NOT drs.skip_merge THEN LBTI_MERGEL1FILES, file
    ; --- Frame selection and cropping for all files
    file      = FILE_SEARCH(pth.l1fits_path + date_obs + drs.dir_label + pth.sep, search_string + '*IMG_ALL.fits', COUNT=n0)
    data_file = FILE_SEARCH(pth.l1fits_path + date_obs + drs.dir_label + pth.sep, search_string + '*DATA_ALL.fits', COUNT=nd0)
    IF n0 LT 1 THEN GOTO, skip_file ELSE IF NOT drs.skip_sel THEN LBTI_IMGSEL, file, data_file, INFO=info, PLOT=plot;,/MEDIAN
    ; --- Derotate (PCA or ADI)
    file      = FILE_SEARCH(pth.l1fits_path + date_obs + drs.dir_label + pth.sep, search_string + '*IMG_SEL.fits', COUNT=n0)
    data_file = FILE_SEARCH(pth.l1fits_path + date_obs + drs.dir_label + pth.sep, search_string + '*DATA_SEL.fits', COUNT=nd0)
    IF n0 LT 1 THEN GOTO, skip_file ELSE LBTI_IMGDEROT, file, data_file, PLOT=plot, VERBOSE=info, /MEDIAN
    ; --- Skipping point
    skip_file:
  ENDFOR
ENDFOR

t2 = SYSTIME(1)
IF info GT 2 THEN  PRINT, 'Time to perform ADI processing     :', t2-t1

; Jumping point if SKIP_NULL is set
skip_adi:

; 5. NULL COMPUTATION
; *******************

; --- Compute null (if nulling data)
t3 = SYSTIME(1)
IF drs.skip_null THEN BEGIN
  IF info GT 0 THEN PRINT, '4. Skipping null computation'
  GOTO, skip_null
ENDIF
IF info GT 0 THEN PRINT, ' '

IF info GT 0 THEN PRINT, '4. Now converting flux to null.'
LBTI_FLX2NULL, date, OB_IDX=ob_idx, LOG_FILE=lun, INFO=info, PLOT=plot, NO_MULTI=no_multi, NO_SAVE=no_save, RENEW=renew
t4 = SYSTIME(1)
IF info GT 2 THEN  PRINT, 'Time to perform null computation     :', t4-t3

; Jumping point if SKIP_NULL is set
skip_null:

; --- Print reduction time to the log file and close it
IF info GT 2 THEN  PRINT, 'Total time to perform data reduction :', SYSTIME(1)-t0
IF lun GT 0 THEN BEGIN
  PRINTF,lun, ' '
  PRINTF,lun, 'Total time to perform data reduction   :', SYSTIME(1)-t0
  PRINTF,lun, ' '
  CLOSE, lun
  FREE_LUN, lun
ENDIF

END
