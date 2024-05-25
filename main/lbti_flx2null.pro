;+
; NAME: LBTI_FLX2NULL
; 
; PURPOSE:
;   This function returns null estimates based on an input flux sequence. This function must work for any input data in a given night, i.e. different
;   targets, instrumental setups, wavelengths, etc... The implicit assumption is that each observing block (OB) contains consistent data[0].
;
; INPUTS:
;   date          :  Date to reduce (format 'yymmdd')
;
; KEYWORDS:
;   OB_IDX        :  Two-element vector with the lower and upper OB ID numbers to process during the null computation
;   LOG_FILE      :
;   INFO          :  Define the level of information printed to screen:
;                       - 0: completely silent execution;
;                       - 1: minimum level of information;
;                       - 2: nominal level of information;
;                       - 3: debugging use.
;   NO_MULTI      :  Set this keyword to turn off multi-threading
;   NO_SAVE       :  Set this keyword to turn-off on data saving
;   PLOT          :  Set this keyword to plot the data to eps files
;   RENEW         :  Set to recompute the OBs already processed
;
; MODIFICATION HISTORY:
;   Version 1.0,  01-OCT-2012, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (adapted from former routine 'nomic_null.pro')
;   Version 1.1,  15-MAY-2012, DD: added new keywords to the output 'fits' files
;   Version 1.2,  30-JUN-2013, DD: added keyword 'SAVE'
;   Version 1.3,  19-SEP-2013, DD: implemented loop over the OBs and improved output files following discussion with Rafael Millan-Gabet
;   Version 1.4,  14-OCT-2013, DD: implemented keyword N_FROB and LOG_FILE
;   Version 1.5,  05-NOV-2013, DD: new way to save L1 files (one per line of the summary table)
;   Version 1.6,  03-MAR-2014, DD: photometry of other OBs is used in case it's missing in a given OB
;   Version 1.7,  08-MAR-2014, DD: now keeping only closed-loop frames for null computation
;   Version 1.8,  28-MAR-2014, DD: added background (+ corresponding plots) and OPD RMS to output files
;   Version 1.9,  09-APR-2014, DD: adapted for new formalism of BIN_NULLDATA.pro
;   Version 2.0,  25-MAY-2014, DD: converted to procedure and removed a bunch of keywords
;   Version 2.1,  29-MAY-2014, DD: added background bias correction
;   Version 2.2,  29-JUL-2014, DD: added error on the photometry to the error on the null
;   Version 2.3,  17-NOV-2014, DD: now reject all null measurements outside drs.null_lim
;   Version 2.4,  19-NOV-2014, DD: added NSC results to output log
;   Version 2.5,  22-NOV-2014, DD: now properly reject OBs with no photometric frames (only when using the NSC pipeline)
;   Version 2.6,  28-NOV-2014, DD: added path to NSC plots
;   Version 2.7,  13-DEC-2014, DD: added errors on NSD-derived terms to output data
;   Version 2.8,  08-JAN-2015, DD: now adopt constructive signal to normalize the null (rather than the total photometry)
;   Version 2.9,  12-JAN-2015, DD: adapted for OBs with multiple integration time
;   Version 3.0,  23-JAN-2015, DD: added error on phasecam telemetry to output data
;   Version 3.1,  06-FEB-2015, DD: added null correction based on high-frequency phase jitter
;   Version 3.2,  10-FEB-2015, DD: added new diagnostic plots
;   Version 3.3,  16-FEB-2015, DD: removed keyword MOVIE
;   Version 3.4,  26-FEB-2015, DD: added raw null RMS to output structure
;   Version 3.5,  14-MAR-2015, DD: added median-based selection of the null
;   Version 3.6,  31-JUL-2015, DD: corrected distribution of photometric data for NSC reduction
;   Version 3.7,  20-OCT-2015, DD: added quality flag
;   Version 3.8,  02-DEC-2015, DD: now moved data defintion to separate routine
;   Version 3.9,  04-DEC-2015, DD: added computation of effective wavelength for PHASEcam (used by NSC)
;   Version 4.0,  23-DEC-2015, DD: Added precision and image mode to output data + corrected FIT_MODE
;   Version 4.1,  03-FEB-2016, DD: Now save intermediate files for each OB (and restore it if reduction paramaters have not changed)
;   Version 4.2,  05-FEB-2016, DD: Added pointing ID
;   Version 4.3,  12-FEB-2016, DD: Corrected minor bug in plots
;   Version 4.4,  16-FEB-2016, DD: Corrected minor bug for NULL_MODE=1 + added plot for H-K phase
;   Version 4.5,  03-APR-2016, DD: Added keyword RENEW
;   Version 4.6,  06-JUL-2016, DD: Added NULL_RANGE to CALL_NSC call
;   Version 4.7,  26-OCT-2016, DD: Now plot also raw nulls (before filtering) + add call to REMOVE_NULLJUMP.pro
;   Version 4.8,  28-OCT-2016, DD: Moved a bunch of plotting stuff to LBTI_IMG2FLX.pro + cleaned code
;   Version 4.9,  29-OCT-2016, DD: Improved output text
;   Version 5.0,  08-NOV-2016, DD: Added null mode to output data + removed error on photometry from the null error per OB because it's correlated between all OBs
;   Version 5.1,  15-NOV-2016, DD: Added PCPHMCOS and PCPHMSIN to call to CALL_NSC.pro
;   Version 5.2,  18-NOV-2016, DD: Added column to output information (null bias)
;   Version 5.3,  31-JAN-2017, DD: Added dither pattern
;   Version 5.4,  17-FEB-2017, DD: Adapted for new CALL_NSC output formalism
;   Version 5.5,  27-MAR-2017, DD: Added the frame mode
;   Version 5.6,  04-APR-2017, DD: Now reads target information from header!
;   Version 5.7,  04-AUG-2017, DD: Updated for new formalism of REMOVE_NULLJUMP.pro
;   Version 5.8,  07-FEB-2018, DD: Added new plots
;   Version 5.8,  15-SEP-2023, DD: Prevent the use of PCPMCOS and PCPHMSIN for 230523 and 230525
;   Version 5.9,  15-OCT-2023, DD: Updated for FRA_MODE=2 (i.e., PCA background subtraction)
;   Version 6.0,  20-FEB-2024, DD: Added FRA_MODE to OB-restoration checking parameters
;   Version 6.1,  03-MAY-2024, DD: Now saved processed L1 files
;   Version 6.2,  24-MAY-2024, DD: Update file permission

PRO LBTI_FLX2NULL, date, OB_IDX=ob_idx, INFO=info, LOG_FILE=log_file, NO_MULTI=no_multi, NO_SAVE=no_save, PLOT=plot, RENEW=renew         

; Global parameter
COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs
ON_ERROR, 0

; Check keywords
IF NOT KEYWORD_SET(INFO)     THEN info = 0
IF NOT KEYWORD_SET(PLOT)     THEN plot = 0 ELSE IF plot GE 2 THEN display = 1
IF NOT KEYWORD_SET(LOG_FILE) THEN lun  = -1 ELSE lun = log_file

; Define outlier rejection limit
sig_out = 5

; If aperture radius is set, redefine dir_label
IF MAX(drs.null_rad) NE 0 AND NOT STRMATCH(drs.dir_label, '*_APR') THEN drs.dir_label = drs.dir_label + '_APR'

; Retrieve data files
data_path  = pth.l1fits_path + drs.date_obs + drs.dir_label + pth.sep
null_files = FILE_SEARCH(data_path,'*NULL.fits', COUNT=n_files)

; Create directories if non existent
sav_fil_path = pth.l1fits_path + pth.sep + drs.date_obs + drs.dir_label + pth.sep + 'filtered' + pth.sep
IF NOT FILE_TEST(sav_fil_path) THEN FILE_MKDIR, sav_fil_path & SPAWN, 'chmod 775 ' + sav_fil_path
sav_int_path = pth.l1fits_path + pth.sep + drs.date_obs + drs.dir_label + pth.sep + 'intermediate' + pth.sep
IF NOT FILE_TEST(sav_int_path) THEN FILE_MKDIR, sav_int_path & SPAWN, 'chmod 775 ' + sav_int_path

; Return if no nulling files
IF n_files EQ 0 THEN BEGIN
  MESSAGE, 'No NULL files to process', /CONTINUE
  RETURN
ENDIF

;Print info to screen
IF info GT 0 THEN BEGIN
  ; Basic info
  PRINT, ''
  PRINT, 'Now processing ' + STRING(n_files, FORMAT='(I0)') + ' OBs'
  IF NOT KEYWORD_SET(no_multi) AND drs.null_mode EQ 2 THEN PRINT, 'Multi-threading on ' + STRING(2^(FLOOR(ALOG(!CPU.HW_NCPU-1)/ALOG(2))), FORMAT='(I0)') + ' cores'
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
  PRINT, 'J: Wavelength [microns]'
  PRINT, 'K: Bandwidths [microns]'
  PRINT, 'L: SNR on photometry (DX and SX)' 
  PRINT, 'M: Intensity mismatch (in %)'  
  PRINT, 'N: Tip/tilt error (DX and SX, in mas)'
  PRINT, 'O: RMS of photometry relative to constructive peak (DX and SX, in %)'  
  PRINT, 'P: RMS of background relative to constructive peak (in %)'
  PRINT, 'Q: Uncorrected background bias relative to constructive peak (in %)'  
  PRINT, 'R: Estimated (from empty region) background bias relative to constructive peak (in %)'  
  PRINT, 'S: Null (in %)'
  PRINT, 'T: Null offset between optimum NSC value and bootstrapped or bayesian value (in %)'  
  PRINT, 'U: Null uncertainty due to error on constructive peak (in %)'
  PRINT, 'V: Null uncertainty due to external error on background floor (in %)'
  PRINT, 'W: Null uncertainty from NSC fit (in %)'  
  PRINT, 'X: Total null uncertainty (in %)'  
  PRINT, 'Y: Expected null uncertainty due to photometric errors (in %)'
  PRINT, 'Z: Expected null uncertainty due to phase variation (i.e., sqrt(4mu^2.sig^2 + 2*sig^4)/4, in %) TO BE CHECKED'  
  PRINT, 'a: Best-fit mean PHASE from NSC (in microns)'
  PRINT, 'b: Best-fit PHASE jitter from NSC (in microns)'
  PRINT, 'c: Best-fit background RMS multiplication factor from NSC'
  PRINT, 'd: Reduced chi2 from NSC'
  PRINT, 'e: Number of rejected null frames'
  PRINT, 'f: Preserved number of frames'  
  PRINT, ' '
  PRINT, ' A ', '  ', '    B   ','  ',' C ','  ','     D     ','  ','   E  ','  ','    F  ','  ','  G ','  ','  H ','  ','  I  ','  ','   J ','  ','  K  ','  ','      L     ','  ','   M   ','  ','     N     ','  ','    O     ','  ','  P  ','  ','  Q  ','  ','  R  ','  ','    S  ',' ','   T  ','  ','    U ','  ','   V  ',' ','   W  ','  ','   X  ','   ','   Y  ','  ','  Z  ','  ','   a  ','  ','   b  ','  ','  c  ','  ','  d  ','  ','  e  ','  ', '  f  '
  PRINT, '---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
  IF lun GT 0 THEN BEGIN
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
    PRINTF, lun, 'J: Wavelength [microns]'
    PRINTF, lun, 'K: Bandwidths [microns]'
    PRINTF, lun, 'L: SNR on photometry (DX and SX)' 
    PRINTF, lun, 'M: Intensity mismatch (in %)'  
    PRINTF, lun, 'N: Tip/tilt error (DX and SX, in mas)'
    PRINTF, lun, 'O: RMS of photometry relative to constructive peak (DX and SX, in %)'  
    PRINTF, lun, 'P: RMS of background relative to constructive peak (in %)'
    PRINTF, lun, 'Q: Uncorrected background bias relative to constructive peak (in %)'  
    PRINTF, lun, 'R: Estimated (from empty region) background bias relative to constructive peak (in %)'  
    PRINTF, lun, 'S: Null (in %)'
    PRINTF, lun, 'T: Null offset between optimum NSC value and bootstrapped or bayesian value (in %)'   
    PRINTF, lun, 'U: Null uncertainty due to error on constructive peak (in %)'
    PRINTF, lun, 'V: Null uncertainty due to external error on background floor (in %)'
    PRINTF, lun, 'W: Null uncertainty from NSC fit (in %)'  
    PRINTF, lun, 'X: Total null uncertainty (in %)'  
    PRINTF, lun, 'Y: Expected null uncertainty due to photometric errors (in %)'
    PRINTF, lun, 'Z: Expected null uncertainty due to phase variation (i.e., sqrt(4mu^2.sig^2 + 2*sig^4)/4, in %) TO BE CHECKED'  
    PRINTF, lun, 'a: Best-fit mean PHASE from NSC (in microns)'
    PRINTF, lun, 'b: Best-fit PHASE jitter from NSC (in microns)'
    PRINTF, lun, 'c: Best-fit background RMS multiplication factor from NSC'
    PRINTF, lun, 'd: Reduced chi2 from NSC'
    PRINTF, lun, 'e: Number of rejected null frames'
    PRINTF, lun, 'f: Preserved number of frames'  
    PRINTF, lun, ' '
    PRINTF, lun, ' A ', '  ', '    B   ','  ',' C ','  ','     D     ','  ','   E  ','  ','    F  ','  ','  G ','  ','  H ','  ','  I  ','  ','   J ','  ','  K  ','  ','      L     ','  ','   M   ','  ','     N     ','  ','    O     ','  ','  P  ','  ','    O     ','  ','  P  ','  ','  Q  ','  ','  R  ','  ','    S  ',' ','   T  ','  ','    U ','  ','   V  ',' ','   W  ','  ','   X  ','   ','   Y  ','  ','  Z  ','  ','   a  ','  ','   b  ','  ','  c  ','  ','  d  ','  ','  e  ','  ', '  f  '
    PRINTF, lun, '---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
  ENDIF
ENDIF

; Initialize output data structure
DEFINE_NULLDATA, data

; Read PHASECam tranmission file 
READ_TABLE, 'nodrs/input/phasecam_K.txt', lam_k, trans, FIRST=0, SEPARATOR=' '
idx_k = WHERE(lam_k GT 1.95 and lam_k LT 2.45, n_lam)
lam_k = lam_k[idx_k]
trans = trans[idx_k]
        
; Prepare plot paths
plot_path = pth.result_path + 'null' + pth.sep + drs.date_obs + drs.dir_label + pth.sep               ; Main path for NULL plots
diag_path = pth.result_path + 'diagnostic' + pth.sep + 'phasecam' + pth.sep + drs.date_obs + pth.sep  ; Main path for DIAGNOSTIC plots
IF NOT FILE_TEST(plot_path) THEN FILE_MKDIR, plot_path & SPAWN, 'chmod 775 ' + plot_path        ; Create directory if it does not exist
IF NOT FILE_TEST(diag_path) THEN FILE_MKDIR, diag_path & SPAWN, 'chmod 775 ' + diag_path        ; Create directory if it does not exist

; Replicate output structure to the number of maximum lines
data    = REPLICATE(data, n_files)

; Jump point to restart reduction using a different photomety (only used when null_mode = 1)
RESTART:

; Loop over the OBs
FOR i_f = 0, n_files-1 DO BEGIN
  
    ; 1. Read null file
    ; *****************
    header    = HEADFITS(null_files[i_f])
    data_null = MRDFITS(null_files[i_f], 1, /SILENT)
    
    ; Make sure this OB is in the user-defined range
    ob_id     = FXPAR(header, 'OB_ID', DATATYPE=0, /NOCONTINUE)
    IF KEYWORD_SET(OB_IDX) THEN IF ob_id LT ob_idx[0] OR ob_id GT ob_idx[1] THEN GOTO, SKIP_OB
            
    ; Read useful information
    objname  = STRCOMPRESS(STRTRIM(FXPAR(header, 'OBJNAME', /NOCONTINUE)), /REMOVE_ALL)
    exptime  = FXPAR(header, 'EXPTIME', /NOCONTINUE)
    lam_cen  = FXPAR(header, 'WAVELENG', /NOCONTINUE)
    bdwdth   = FXPAR(header, 'BANDWIDT', /NOCONTINUE)
    bck_mode = FXPAR(header, 'BCK_MODE', /NOCONTINUE)
    bfl_mode = FXPAR(header, 'BFL_MODE', /NOCONTINUE)
    fit_mode = FXPAR(header, 'FIT_MODE', /NOCONTINUE)
    flx_mode = FXPAR(header, 'FLX_MODE', /NOCONTINUE)
    fra_mode = FXPAR(header, 'FRA_MODE', /NOCONTINUE)
    img_mode = FXPAR(header, 'IMG_MODE', /NOCONTINUE)
    ob_mode  = FXPAR(header, 'OB_MODE', /NOCONTINUE)
    pixscale = FXPAR(header, 'PIXSCALE', /NOCONTINUE) 
    precisio = FXPAR(header, 'PRECISIO', /NOCONTINUE)
    pt_id    = FXPAR(header, 'PT_ID', /NOCONTINUE)
    flag     = FXPAR(header, 'FLAG', /NOCONTINUE)
    n_rej    = FXPAR(header, 'N_REJ', /NOCONTINUE) > 0
    pcspec   = FXPAR(header, 'TGT_SPEC', /NOCONTINUE) 
    calfor   = FXPAR(header, 'TGT_CAL4', /NOCONTINUE) 
           
    ; Skip OB if not enough frames
    n_null0 = N_ELEMENTS(data_null.flx_tot[0])
    IF n_null0 LT drs.min_fr THEN BEGIN
      IF info GT 0 THEN PRINT, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ', 'Not enough null measurements in OB (only ' + STRING(n_null0, FORMAT='(I0)') + ')'
      IF lun GT 0 THEN PRINTF, lun,  STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ', 'Not enough null measurements in OB (only ' + STRING(n_null0, FORMAT='(I0)') + ')'
      GOTO, SKIP_OB
    ENDIF
    
    ; Compute effective wavelength for PHASECam (read target spectrum and instrument throughput)
    pick   = READ_PICKLES(pth.pickles_path + pcspec, LAM_MIN=1.95D-6, LAM_MAX=2.45D-6, NLAMBDA=n_lam, STANDARD=3)
    spec   = INTERPOL(pick[1,*], pick[0,*], lam_k*1D-6)
    pcwav  = 1./(TOTAL((trans*spec)^2*1./lam_k)/TOTAL((trans*spec)^2))*1D-6

    ; Compute EEID to determine the aperture size (use the same for SCI and CAL)
    IF drs.null_rad EQ 0 THEN BEGIN
      IF STRTRIM(flag, 2) EQ 'CAL' THEN GET_TGT, calfor, tgt_sci, DATABASE= pth.input_path + drs.database ELSE GET_TGT, objname, tgt_sci, DATABASE= pth.input_path + drs.database
      eeid_as = COMPUTE_EEID(0.5*tgt_sci.ldm, tgt_sci.dist, tgt_sci.temp)
      aper_rad = ROUND(eeid_as/pixscale) > 8 
    ENDIF ELSE aper_rad = drs.null_rad
    
    ; Find closest aperture radius. If drs.aper_rad is 0, it means that we have to use the EEID.
    ; Use the drs.aper_rad otherwise
    n_apr    = N_ELEMENTS(data_null[0].flx_tot)
    apr_all  = FXPAR(header, 'APERRAD' + STRING(0, FORMAT='(I0)'))
    FOR i_r = 1, n_apr-1 DO apr_all = [apr_all, FXPAR(header, 'APERRAD' + STRING(i_r, FORMAT='(I0)'))]
    tmp      = MIN(ABS(apr_all-aper_rad), i_aper) 
    apr_rad  = apr_all[i_aper]
    napr_pix = !dpi*apr_rad^2    ; Number of pixels in photometric aperture and in annulus
    bck_irad = FXPAR(header, 'BCKIRAD' + STRING(i_aper, FORMAT='(I0)'))
    bck_orad = FXPAR(header, 'BCKORAD' + STRING(i_aper, FORMAT='(I0)'))
    nsky_pix = FXPAR(header, 'NSKYPIX' + STRING(i_aper, FORMAT='(I0)'))  
    IF !err EQ -1 THEN nsky_pix = !dpi*(bck_orad^2-bck_irad^2)           ; this was used for a test but not needed after all
  
    ; Restore OB if already reduced with the same parameters (and saved)
    ; Not working for null_mode = 1 because the reduction neeeds to know the photometry of the faintest star and it is only computed later
    ; But not a big deal because null_mode = 1 is not as slow as NSC
    IF drs.null_mode NE 1 AND NOT KEYWORD_SET(RENEW) THEN BEGIN
      files_ob = FILE_SEARCH(data_path,'null_ob' + STRING(ob_id, FORMAT='(I03)') + '*.sav', COUNT=n_files_ob)
      ; Loop over files of this OB until matching paramaters are found (or not)
      FOR i_fob = 0, n_files_ob-1 DO BEGIN
        ; Restore old file
        RESTORE, files_ob[i_fob]
        ; IF no BFL_MODE, consider as an obsolote file and keep looking
        IF TAG_EXIST(data_ob, 'bfl_mode') THEN BEGIN
          ; Compare to running/current parameters
          IF drs_ob.null_mode     EQ drs.null_mode     AND drs_ob.keep_ratio  EQ drs.keep_ratio  AND drs_ob.min_fr      EQ drs.min_fr      AND drs_ob.n_btstrp    EQ drs.n_btstrp    AND drs_ob.nsc_bfac    EQ drs.nsc_bfac    AND $
             drs_ob.nsc_bins      EQ drs.nsc_bins      AND drs_ob.nsc_cube[0] EQ drs.nsc_cube[0] AND drs_ob.nsc_cube[1] EQ drs.nsc_cube[1] AND drs_ob.nsc_cube[2] EQ drs.nsc_cube[2] AND drs_ob.nsc_cube[3] EQ drs.nsc_cube[3] AND $
             drs_ob.nsc_cube[4]   EQ drs.nsc_cube[4]   AND drs_ob.nsc_omin    EQ drs.nsc_omin    AND drs_ob.null_cor    EQ drs.null_cor    AND drs_ob.null_rad    EQ drs.null_rad    AND drs_ob.null_lim[0] EQ drs.null_lim[0] AND $
             drs_ob.null_lim[1]   EQ drs.null_lim[1]   AND drs_ob.null_range[0] EQ drs.null_range[0] AND drs_ob.null_range[1] EQ drs.null_range[1] AND drs.n_frob EQ drs.n_frob      AND $ ; end of drs parameters
             data_ob.objname      EQ objname       AND data_ob.lam_cen    EQ lam_cen         AND data_ob.bck_mode   EQ bck_mode        AND data_ob.bfl_mode   EQ bfl_mode        AND data_ob.fit_mode   EQ fit_mode AND data_ob.fra_mode   EQ fra_mode        AND data_ob.img_mode EQ img_mode     AND $
             data_ob.ob_mode      EQ ob_mode       AND data_ob.pixscale   EQ pixscale        AND data_ob.precision  EQ precisio        AND data_ob.bck_irad   EQ bck_irad        AND data_ob.bck_orad EQ bck_orad THEN BEGIN                                                                              ; end of paramaters on L1 images
            data[i_f] = data_ob
            GOTO, RESTORED_OB
          ENDIF
        ENDIF
      ENDFOR
    ENDIF
    
    ; Prepare plot and log paths (define plot path anyway because used later for logs)
    plot_path_ob  = plot_path + 'aper-' + STRING(apr_rad, FORMAT='(I0)') + 'pix' + pth.sep + 'ob' + STRING(ob_id, FORMAT='(I0)') + pth.sep
    plot_path_nsc = plot_path_ob + 'nsc' + pth.sep
    IF NOT FILE_TEST(plot_path_ob) THEN FILE_MKDIR, plot_path_ob & SPAWN, 'chmod 775 ' + plot_path_ob
    IF NOT FILE_TEST(plot_path_nsc) THEN FILE_MKDIR, plot_path_nsc & SPAWN, 'chmod 775 ' + plot_path_nsc
    
    ; Plot results
    plotnull = plot_path_ob + drs.date_obs + '_OB' + STRING(ob_id, FORMAT='(I03)') + '_' + objname + '_dit-' + STRING(1000.*exptime, FORMAT='(I0)') + 'ms_wav-' + STRING(1D+6*lam_cen, FORMAT='(I0)') + 'um'
    plotdiag = diag_path + drs.date_obs + '_OB' + STRING(ob_id, FORMAT='(I03)') + '_' + objname + '_dit-' + STRING(1000.*exptime, FORMAT='(I0)') + 'ms_wav-' + STRING(1D+6*lam_cen, FORMAT='(I0)') + 'um'

    ; 2. Read and process photometry files
    ; ************************************
    
    data_phot1 = LBTI_READL1DATA(data_path, ob_id, 'PHOT1', APER=apr_rad, IDX_APER=i_aper1, HEADER=hdr_phot1, /FILTER)
    data_phot2 = LBTI_READL1DATA(data_path, ob_id, 'PHOT2', APER=apr_rad, IDX_APER=i_aper2, HEADER=hdr_phot2, /FILTER)    
    IF (SIZE(data_phot1))[2] EQ 8 AND (SIZE(data_phot2))[2] EQ 8 THEN BEGIN
      phot_tot_phot1 = data_phot1.flx_tot[i_aper1] & phot_err_phot1 = data_phot1.flx_err[i_aper1] & n_phot1 = N_ELEMENTS(data_phot1.mjd_obs)
      phot_tot_phot2 = data_phot2.flx_tot[i_aper2] & phot_err_phot2 = data_phot2.flx_err[i_aper2] & n_phot2 = N_ELEMENTS(data_phot2.mjd_obs)
      IF n_phot1 EQ 1 THEN BEGIN
        phot1          = phot_tot_phot1
        phot1_err      = phot_err_phot1
        phot1_err_mean = phot_err_phot1
      ENDIF ELSE AVGSDV, phot_tot_phot1, phot1, phot1_err, phot1_err_mean, KAPPA=3
      IF n_phot2 EQ 1 THEN BEGIN
        phot2          = phot_tot_phot2
        phot2_err      = phot_err_phot2
        phot2_err_mean = phot_err_phot2
      ENDIF ELSE AVGSDV, phot_tot_phot2, phot2, phot2_err, phot2_err_mean, KAPPA=3
      ; Compute rms tipt/tilt
      tt1 = SQRT((data_phot1.xcen)^2 + (data_phot1.ycen)^2)*pixscale*1D3 & AVGSDV, tt1, avg, rms_tt1
      tt2 = SQRT((data_phot2.xcen)^2 + (data_phot2.ycen)^2)*pixscale*1D3 & AVGSDV, tt2, avg, rms_tt2      
      ; Adjust error bar if negative background mode
      ;IF bck_mode GE 0 THEN BEGIN
      ;  phot1_err = SQRT(phot1_err^2+phot1_err^2/      1+FXPAR(header1, 'N_FRBCK', /NOCONTINUE)/N_ELEMENTS(phot_tot_phot1))
      ;  phot2_err = SQRT(1+FXPAR(header2, 'N_FRBCK', /NOCONTINUE)/N_ELEMENTS(phot_tot_phot2))
      ;ENDIF
      IF n_phot1 GT 1 AND n_phot2 GT 1 THEN BEGIN 
        time1 = REFORM(TRANSPOSE(data_phot1.mjd_obs-MIN(data_phot1.mjd_obs))*24.*60.*60.)      ; convert MJD to ellapsed seconds
        time2 = REFORM(TRANSPOSE(data_phot2.mjd_obs-MIN(data_phot2.mjd_obs))*24.*60.*60.)      ; convert MJD to ellapsed seconds
        PLOTALL, time1, REFORM(phot_tot_phot1), REFORM(phot_err_phot1), NAME=plotnull, TAG='PHOT1', XTITLE='Elapsed time [s]', YTITLE='Photometric measurements [ADU]', TITLE=' ', /BOLD, /EPS
        PLOTALL, time2, REFORM(phot_tot_phot2), REFORM(phot_err_phot2), NAME=plotnull, TAG='PHOT2', XTITLE='Elapsed time [s]', YTITLE='Photometric measurements [ADU]', TITLE=' ', /BOLD, /EPS
      ENDIF
    ENDIF ELSE BEGIN
       ; Skip OB if NSC reduction. Use max interferometric flux otherwise.
       IF drs.null_mode EQ 2 THEN BEGIN
         IF info GT 0 THEN PRINT, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ', 'No photometric frames'
         IF lun  GT 0 THEN PRINTF, lun, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ', 'No photometric frames'
         GOTO, SKIP_OB 
       ENDIF ELSE PRINT, 'WARNING: no photometric frames for OB ' + STRING(i_f, FORMAT='(I0)') + ' (using maximum inteferometric flux as photometry)' 
       phot1     = 0.25*MAX(data_null.flx_tot[i_aper], i_max) & phot2 = phot1
       phot1_err = data_null[i_max].flx_err[i_aper]           & phot2_err = phot1_err
       phot1_err_mean = phot1_err                             & phot2_err_mean = phot2_err
       rms_tt1 = 0 & rms_tt2 = 0
    ENDELSE
    int_err  = 2.0*ABS(phot1-phot2)/(phot1+phot2)
    phot_tot = phot1+phot2+2*SQRT(ABS(phot1*phot2)) 
    phot_err = 2.0*SQRT(phot1_err_mean^2+phot2_err_mean^2)  ; first order estimate dN = N*phot_err/phot_tot
    
    ; If using the best x% null technique, restart the reduction using the faintest CAL to correct for the background bias
    IF drs.null_mode EQ 1 THEN BEGIN
      IF N_ELEMENTS(phot_faintcal) EQ 0 THEN phot_faintcal = phot_tot 
      IF phot_tot LT phot_faintcal THEN BEGIN
        phot_faintcal = phot_tot
        IF info GT 0 THEN PRINT, 'WARNING: restarting null computation using the fainter calibrator for background bias correction'
        IF lun GT 0 THEN PRINTF, lun, 'WARNING: restarting null computation using the fainter calibrator for background bias correction'
        GOTO, RESTART
      ENDIF
    ENDIF

    ; Save filtered L1 file
    IF NOT KEYWORD_SET(no_save) THEN LBTI_SAVEL1FLX_1APER, data_phot1, hdr_phot1, i_aper1, OB_ID=ob_id, TAG='PHOT1'
    IF NOT KEYWORD_SET(no_save) THEN LBTI_SAVEL1FLX_1APER, data_phot2, hdr_phot2, i_aper2, OB_ID=ob_id, TAG='PHOT2'
    
    ; 3. Read and process background file
    ; ***********************************
    
    ; Process background files only if the background mode is negative (which means that a constant frame has been subtraced from each frame of the current nod and
    ; background files are obtained at the same position as the the beam). If the background mode is greater than 0, use background mesured at a different position 
    ; on the detector (so, stored directly in the NULL files). 
    data_bckg = LBTI_READL1DATA(data_path, ob_id, 'BCKG', APER=apr_rad, IDX_APER=id_aper, HEADER=hdr_bckg, /FILTER)
    IF (SIZE(data_bckg))[2] EQ 8 AND bck_mode LT 0 THEN BEGIN
      bck_tot_phot  = data_bckg.flx_tot[id_aper]   ; Photometric aperture flux in background nod
      bck_err_phot  = data_bckg.flx_err[id_aper]   ; Corresponding error
      bck_bck_flr   = data_bckg.bck_tot[id_aper]   ; background floor measured in complementary nod
      bck_tot_phot2 = data_bckg.flx_tot2[id_aper]  ; Photometric aperture flux in background nod
      bck_bck_flr2  = data_bckg.bck_tot2[id_aper]  ; background floor measured in complementary nod
      ;bck_bck_eflr  = data_bckg.bck_err[id_aper]  ; error on background floor
      nod_bckg      = data_bckg.nod_id             ; NOD ID of each background flux
      ; Plot data of requested
      IF N_ELEMENTS(bck_tot_phot) GT 1 THEN BEGIN
        time = REFORM(TRANSPOSE(data_bckg.mjd_obs-MIN(data_bckg.mjd_obs))*24.*60.*60.)     ; convert MJD to ellapsed seconds
        PLOTALL, time, REFORM(bck_tot_phot), REFORM(bck_err_phot), NAME=plotnull, TAG='BCKG', XTITLE='Elapsed time [s]', YTITLE='Background measurements [ADU]', TITLE='', /BOLD, /EPS
        PLOTALL, time, REFORM(bck_bck_flr),  0., NAME=plotnull, TAG='FLOOR-BCKG', XTITLE='Elapsed time [s]', YTITLE='Background floor [ADU]', TITLE='', /BOLD, /EPS
        PLOTALL, REFORM(bck_bck_flr), REFORM(bck_tot_phot), 0., NAME=plotnull, TAG='BIAS-vs-FLOOR', XTITLE='Background floor [ADU]', YTITLE='Background bias [ADU]', TITLE='', /BOLD, /EPS, /SCATTER, /NO_FFT
        IF MAX(ABS(bck_bck_flr2)) NE 0 THEN PLOTALL, REFORM(bck_bck_flr2), REFORM(bck_tot_phot2), 0., NAME=plotnull, TAG='BIAS2-vs-FLOOR2', XTITLE='Background floor [ADU]', YTITLE='Background bias [ADU]', TITLE='', /BOLD, /EPS, /SCATTER, /NO_FFT
      ENDIF
      ; REGION 1 (with the star). 
      ; Compute uncorrected background BIAS
      IF fra_mode EQ 1 THEN AVGSDV, bck_tot_phot, bck_bias_unc, junk1, junk2, KAPPA=5 ELSE bck_bias_unc = 0.  ; If frame mode of 0, then the bias is corrected by image subtraction and we have no info here...
    ENDIF ELSE BEGIN
      coeff        = [0,0]
      nod_bckg     = 0
      bck_avg      = 0
      bck_rms      = 0
      bck_rms_mean = 0   ; this is the external error, which is included in the NSC error if bck_mode is positive 
      bck_bias_unc = 0
      ; Skip OB if no BCKG file
      IF bck_mode LT 0 THEN BEGIN
        IF info GT 0 THEN PRINT, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ', 'No background frames found (skipped)'
        IF lun GT 0  THEN PRINTF, lun, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ', 'No background frames found (skipped)'
        GOTO, SKIP_OB
      ENDIF
    ENDELSE

    ; Save filtered L1 file
    IF NOT KEYWORD_SET(no_save) THEN LBTI_SAVEL1FLX_1APER, data_bckg, hdr_bckg, id_aper, OB_ID=ob_id, TAG='BCKG'
   
    ; 4. Process and filter nulls
    ; ***************************
    
    ; Plot diagnostic infos about the PWV loop (+log file for Phil, only of H-band data)
    ; Plot null sequence before any filtering
    time          = REFORM(TRANSPOSE(data_null.mjd_obs-MIN(data_null.mjd_obs))*24.*60.*60.)      ; convert MJD to ellapsed seconds
    null_tot_phot = REFORM(data_null.flx_tot[i_aper])-bck_bias_unc                               ; correct for background BIAS
    PLOTALL, time, 100.*REFORM(data_null.flx_tot[i_aper])/phot_tot, 100.*REFORM(data_null.flx_err[i_aper])/phot_tot, NAME=plotnull, TAG='NULL', XTITLE='Elapsed time [s]', YTITLE='Instantaneous null depth [%]', TITLE='', /YLOG, YRANGE=[0.1,100.0], /BOLD, /EPS, /KERNEL
    ; Plot PHASECam diagnostic INFO
    IF MAX(data_null.pcphmean2) NE 0 THEN BEGIN      
      ; Extract data        
      pwv_pha = data_null.pcphmean2 - data_null.pcphmean
      PLOTALL, pwv_pha, 100.*null_tot_phot/phot_tot, 0., NAME=plotdiag, TAG='PWV-NULL', XTITLE='H-K phase [deg]', YTITLE='Instantaneous null depth [%]', TITLE='', /BOLD, /EPS, /NO_FFT, /SCATTER
      PLOTALL, time, pwv_pha, 0, NAME=plotdiag, TAG='PWV', XTITLE='Elapsed time [s]', YTITLE='H-K phase [deg]', TITLE='',  /BOLD, /EPS
      ; Temporary hack
      PREP_PS, /BOLD & LOADCT, 0, /SILENT
      DEVICE, FILENAME = plotdiag + '_PWV-vs-NULL.eps', /ENCAPS, /COLOR, XSIZE=17.8, YSIZE=14.7, /TIMES
      nulll = 100.*null_tot_phot/phot_tot
      xrange = [MIN(time),MAX(time)]
      yrange = [MIN(nulll),MAX(nulll)]
      PLOT, [0], [0], XTITLE='Elapsed time [s]', YTITLE='Null depth [%]', TITLE=title_day, XSTYLE=1, YSTYLE=9, $
            XRANGE=xrange, YRANGE=yrange, XTHICK=xthick, YTHICK=ythick, CHARTHICK=charthick, CHARSIZE=charsize
      AXIS, YTITLE ='H-K phase [deg]', /YAXIS, YSTYLE=1, YRANGE=[MIN(pwv_pha), MAX(pwv_pha)]
      LOADCT, 13, /SILENT
      OPLOT, time, nulll, COLOR=100
      OPLOT, time, (pwv_pha-MIN(pwv_pha))*(MAX(nulll)-MIN(nulll))/(MAX(pwv_pha)-MIN(pwv_pha)), COLOR=250
      DEVICE, /CLOSE & END_PS

      ; Print file for Phil
      log_file2 = diag_path + date + '_OB' + STRING(i_f, FORMAT='(I0)') + '.txt'
      OPENW, lun2, log_file2, /GET_LUN, WIDTH=800
      PRINTF,lun2, 't0=' + STRING(data_bckg[0].mjd_obs, FORMAT='(F16.8)')
      PRINTF,lun2, 'time[s];null[%];PHAB1[deg];PHAB2[deg]
      FOR i=0, N_ELEMENTS(nulll)-1 DO PRINTF,lun2, time[i],';', nulll[i],';',data_null[i].pcphmean,';',data_null[i].pcphmean2
      CLOSE, lun2
      FREE_LUN, lun2
    ENDIF 
    
    ; Now filter the null data.
    ; ------------------------
     
    ; First, remove fringe jump. 18% is considerd as a jump, keep lowest sequence + those with a mean within 5% of the lowest one(and enough frames!) 
    ; Only for data after July 2014 (otherwise, CG only and this does not apply)
    IF MEAN(data_null.mjd_obs) GT 56839 THEN BEGIN
      tmp_null = data_null.flx_tot[i_aper]/phot_tot
      tmp_null = REMOVE_NULLJUMP(tmp_null, 0.18, 0.05, 20, IDX_OUT=idx_null) 
      n_null   = N_ELEMENTS(idx_null)
      IF n_null LT drs.min_fr THEN BEGIN
        IF info GT 0 THEN PRINT, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived fringe jump filtering (skipped)'
        IF lun GT 0  THEN PRINTF, lun,  STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived fringe jump filtering  (skipped)'
        GOTO, SKIP_OB
      ENDIF
    ENDIF ELSE idx_null = LINDGEN(N_ELEMENTS(data_null.mjd_obs))

    ; Compute mean fringe SNR and keep only the best one (low SNR usualy means bad phase setpoint)
    MEANCLIP, data_null[idx_null].pcmsnr, avg_snr, rms_snr, CLIPSIG=sig_out
    idx_null = idx_null[WHERE(data_null[idx_null].pcmsnr GE (avg_snr-sig_out*rms_snr), n_null)]
    IF n_null LT drs.min_fr THEN BEGIN
       IF info GT 0 THEN PRINT, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived the SNR cut (skipped)'
       IF lun GT 0  THEN PRINTF, lun,  STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived the SNR cut  (skipped)'
       GOTO, SKIP_OB
    ENDIF
    ; Compute mean phase noise and reject high noise
    MEANCLIP, data_null[idx_null].pcphstd, avg_phstd, rms_phstd, CLIPSIG=sig_out
    idx_null = idx_null[WHERE(ABS(data_null[idx_null].pcphstd-avg_phstd) LE sig_out*rms_phstd, n_null)]
    IF n_null LT drs.min_fr THEN BEGIN
      IF info GT 0 THEN PRINT, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived the phase noise cut (skipped)'
      IF lun GT 0  THEN PRINTF, lun, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived the phase noise cut  (skipped)'
      GOTO, SKIP_OB
    ENDIF
    ; Compute mean phase and reject frame to far from setpoint
    MEANCLIP, data_null[idx_null].pcphmean, avg_phmean, rms_phmean, CLIPSIG=sig_out
    idx_null = idx_null[WHERE(ABS(data_null[idx_null].pcphmean-avg_phmean) LE sig_out*rms_phmean, n_null)]
    IF n_null LT drs.min_fr THEN BEGIN
      IF info GT 0 THEN PRINT, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived the mean phase cut (skipped)'
      IF lun GT 0  THEN PRINTF, lun, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived the mean phase cut  (skipped)'
      GOTO, SKIP_OB
    ENDIF
    ; Only keep the frame in null_lim range
    IF N_ELEMENTS(drs.null_lim) EQ 2 THEN BEGIN
      IF drs.null_lim[0] NE 0 OR drs.null_lim[1] NE 0 THEN BEGIN
        flx_tmp  = data_null[idx_null].flx_tot[i_aper] - bck_bias_unc ; correct for background bias
        idx_null = idx_null[WHERE(flx_tmp/phot_tot GE drs.null_lim[0] AND flx_tmp/phot_tot LE drs.null_lim[1], n_null, /NULL)]
        IF n_null LT drs.min_fr THEN BEGIN
          IF info GT 0 THEN PRINT, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') in the NULL_LIM range (skipped)'
          IF lun GT 0  THEN PRINTF, lun,  STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') in the NULL_LIM range (skipped)'
          GOTO, SKIP_OB
        ENDIF
      ENDIF
    ENDIF
    ; Remove frames more than 5 sigma away from the median
    flx_tmp = data_null[idx_null].flx_tot[i_aper]
    AVGSDV, flx_tmp, avg_tmp, rms_tmp, rms_tmp_m, KAPPA=sig_out
    idx_null = idx_null[WHERE(ABS(flx_tmp-MEDIAN(flx_tmp)) LT sig_out*rms_tmp, n_null, /NULL)] 
    IF n_null LT drs.min_fr THEN BEGIN
      IF info GT 0 THEN PRINT, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived after the median cut (skipped)'
      IF lun GT 0  THEN PRINTF, lun, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived after the median cut (skipped)'
      GOTO, SKIP_OB
    ENDIF
    ; Use background in surrounding regions to remove outliers
    flx_tmp = data_null[idx_null].bck_tot[i_aper]
    AVGSDV, flx_tmp, avg_tmp, rms_tmp, rms_tmp_m, KAPPA=sig_out
    idx_null = idx_null[WHERE(ABS(flx_tmp-avg_tmp) LE sig_out*rms_tmp, n_null, /NULL)]
    IF n_null LT drs.min_fr THEN BEGIN
      IF info GT 0 THEN PRINT, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived after the background check (skipped)'
      IF lun GT 0  THEN PRINTF, lun, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(objname, FORMAT='(A8)'), '  ', STRING(flag, FORMAT='(A3)'), '  ','Not enough null frames (' + STRING(n_null, FORMAT='(I0)') + ') survived after the background check (skipped)'
      GOTO, SKIP_OB
    ENDIF    
    ; Extract data
    data_null = data_null[idx_null]
    ; Extract data for chosen aperture radius
    null_tot_phot  = REFORM(data_null.flx_tot[i_aper])  & null_err_phot  = REFORM(data_null.flx_err[i_aper])     ; flux measurements at null and corresponding errors
    null_tot_bckg  = REFORM(data_null.bck_tot[i_aper])  & null_err_bckg  = REFORM(data_null.bck_err[i_aper])     ; background measurements (per pixel) and corresponding errors (in the background region surrounding the photometric aperture) 
    null_tot_phot2 = REFORM(data_null.flx_tot2[i_aper]) & null_err_phot2 = REFORM(data_null.flx_err2[i_aper])    ; flux measurements at null and corresponding errors
    null_tot_bckg2 = REFORM(data_null.bck_tot2[i_aper]) & null_err_bckg2 = REFORM(data_null.bck_err2[i_aper])    ; background measurements (per pixel) and corresponding errors (in the background region surrounding the photometric aperture) 
    
    ; Save filtered L1 file
    IF NOT KEYWORD_SET(no_save) THEN LBTI_SAVEL1FLX_1APER, data_null, header, i_aper, OB_ID=ob_id, TAG='NULL'

    ; Debias null and off-axis measurements by using measurements in the complementary NOD (only if FRA_MODE of 1)
    ; ----------------------------------------------------------------------------------
    
    ; First, determine degree of POLYNOMIAL fit by comparing the range of background floors in the complentary NOD(s) and the science nod
    ; The rule of thumb here is to use a constant if the difference between the mean background floors in each NOD is less than twice the standard deviation of background floors in teh complementary NOD(s)
    ; Use a degree of 2 otherwise 
    IF fra_mode EQ 1 THEN BEGIN
      AVGSDV, bck_bck_flr, bck_bck_avg, bck_bck_rms, bck_bck_mrms, KAPPA=5                       ; compute mean and RMS of background floors in complementary nods
      IF (MEAN(null_tot_bckg)-bck_bck_avg) LE 2*bck_bck_rms THEN poly_deg = 1 ELSE poly_deg = 0  ; compute polynomial degree
      coeff = POLY_FIT(bck_bck_flr, bck_tot_phot, poly_deg)                                      ; fit polynomial based on complementary nods
      IF MAX(ABS(bck_tot_phot2)) GT 0 THEN coeff2 = POLY_FIT(bck_bck_flr2, bck_tot_phot2, poly_deg) ELSE coeff2 = REPLICATE(0, poly_deg+1)
      ; Now debias background measurements in complementary NODS (both regions of the detector)
      FOR i = 0, poly_deg DO BEGIN
        bck_tot_phot   -= coeff[i]*bck_bck_flr^i     ; FLux in complementary nod (in BCKG file)
        null_tot_phot  -= coeff[i]*null_tot_bckg^i   ; Flux on source
        null_tot_phot2 -= coeff2[i]*null_tot_bckg2^i ; Flux in nearby empty region (in NULL file)
      ENDFOR  
    ENDIF 
    
    ; Now compute NSC-related values and PLOT results
    ; -----------------------------------------------
    
    ; Remove mean value from background for NSC (if from different NODs, do it per NOD!)
    nod_bckg_uniq = nod_bckg[UNIQ(nod_bckg,  SORT(nod_bckg))]
    n_nod_bckg    = N_ELEMENTS(nod_bckg_uniq)
    FOR ib = 0, n_nod_bckg-1 DO BEGIN
      idx_nod_bckg = WHERE(nod_bckg EQ nod_bckg_uniq[ib])
      AVGSDV, bck_tot_phot[idx_nod_bckg], bck_avg, tmp, tmp2, KAPPA=5
      bck_tot_phot[idx_nod_bckg] = bck_tot_phot[idx_nod_bckg] - bck_avg   ; remove mean value for NSC fitting
    ENDFOR
    ; Plot corrected background
    IF N_ELEMENTS(bck_tot_phot) GT 1 THEN BEGIN
      time = REFORM(TRANSPOSE(data_bckg.mjd_obs-MIN(data_bckg.mjd_obs))*24.*60.*60.)     ; convert MJD to ellapsed seconds
      PLOTALL, time, REFORM(bck_tot_phot), REFORM(bck_err_phot), NAME=plotnull, TAG='BCKG-COR', XTITLE='Elapsed time [s]', YTITLE='Background measurements [ADU]', TITLE='', /BOLD, /EPS
    ENDIF
    ; Now compute RMS over full sequence
    AVGSDV, bck_tot_phot, tmp, bck_rms, bck_rms_mean;, KAPPA=5
    ; Derive minimum possible astro NULL (defined arbitrary here as -1 sigma the background variation relative to peak)
    nas_min = -1.*bck_rms/phot_tot
    ; Compute raw null RMS and AVG (avter debias)
    AVGSDV, null_tot_phot/phot_tot, null_avg, null_rms
    ;Degrade null (used to simulate fainter stars)
    ;IF STRTRIM(FXPAR(header, 'FLAG')) EQ 'SCI' THEN BEGIN
    ;  null_tot_phot += bck_tot_phot*SQRT(7./2.-1)
    ;  bck_tot_phot *= SQRT(7./2.)
    ;ENDIF
    ; Compute mean background bias
    IF MAX(null_tot_phot2) NE 0 THEN BEGIN
      idx_ok = WHERE(null_err_phot2 NE 0, n_ok)
      null_tot_bckg2 = null_tot_bckg2[idx_ok]
      null_tot_phot2 = null_tot_phot2[idx_ok]
      null_err_phot2 = null_err_phot2[idx_ok]
      AVGSDV, null_tot_phot2, bck_bias_cor, rms_bias, rms_bias_mean, KAPPA=sig_out, WEIGHT=(1/null_err_phot2)^2
      null_bckg = [bck_bias_unc, bck_bias_cor, rms_bias, rms_bias_mean, 0]/phot_tot   ; Uncorrected background bias (same position), corrected background bias (other position), RMS at other position, and RMS on the mean
    ENDIF ELSE BEGIN
      null_bckg = [bck_bias_unc, 0., 0., 0]/phot_tot 
    ENDELSE
    ; Compute mean background level (per pixel)
    AVGSDV, null_tot_bckg, avg_bck, rms_bck, KAPPA=sig_out, WEIGHT=(1./null_err_bckg)^2
    bck_rms_est = SQRT(napr_pix)*SQRT(nsky_pix)*rms_bck     ; estimated RMS of total background in photometric aperture (not used after all)
    ; Plot data of requested  
    time = REFORM(TRANSPOSE(data_null.mjd_obs-MIN(data_null.mjd_obs))*24.*60.*60.)      ; convert MJD to ellapsed seconds
    PLOTALL, time, 100.*null_tot_phot/phot_tot, 100.*null_err_phot/phot_tot, NAME=plotnull, TAG='NULLOK', XTITLE='Elapsed time [s]', YTITLE='Instantaneous null depth [%]', TITLE='', /YLOG, YRANGE=[0.1,100.0], /BOLD, /EPS, /KERNEL 
    PLOTALL, time, null_tot_bckg, null_err_bckg, NAME=plotnull, TAG='FLOOR-NULL', XTITLE='Elapsed time [s]', YTITLE='Background floor [ADU]', TITLE='', /BOLD,  /EPS
    PLOTALL, null_tot_bckg, null_tot_phot, 0, NAME=plotnull, TAG='NULL-VS-FLOOR', XTITLE='Background floor [ADU]', YTITLE='Background bias [ADU]', TITLE='', /BOLD, /EPS, /SCATTER, /NO_FFT
    ;PLOTALL, time, fake_bckg, 0, NAME=plotnull, TAG='FAKE-BCKG', XTITLE='Elapsed time [s]', YTITLE='Estimated backgroud in photometric aperture [ADU]', TITLE='', /BOLD,  /EPS
    ;PLOTALL, data_null.pcphmean, 100.*null_tot_phot/phot_tot, 0., NAME=plotnull, TAG='PCPHMEAN2', XTITLE='Phase [deg]', YTITLE='Instantaneous null depth [%]', TITLE='', /BOLD, /EPS, /NO_FFT, /SCATTER
    ;PLOTALL, data_null.pcphstd, 100.*null_tot_phot/phot_tot, 0., NAME=plotnull, TAG='PCPHSTD', XTITLE='Phase jitter [deg]', YTITLE='Instantaneous null depth [%]', TITLE='', /BOLD, /EPS, /NO_FFT, /SCATTER
    ;PLOTALL, data_null.fpcazm, 100.*null_tot_phot/phot_tot, 0., NAME=plotnull, TAG='FPCAZM', XTITLE='Azimuth [mas]', YTITLE='Instantaneous null depth [%]', TITLE='', /BOLD, /EPS, /NO_FFT, /SCATTER
    ;PLOTALL, data_null.fpcelm, 100.*null_tot_phot/phot_tot, 0., NAME=plotnull, TAG='FPCELM', XTITLE='Elevation [mas]', YTITLE='Instantaneous null depth [%]', TITLE='', /BOLD, /EPS, /NO_FFT, /SCATTER
    IF MAX(data_null.pcphstd)  NE 0 THEN PLOTALL, time, data_null.pcphstd, 0., NAME=plotnull, TAG='PCPHSTD', XTITLE='Elapsed time [s]', YTITLE='Phase jitter wittin DIT [deg]', TITLE='', /BOLD, /EPS;, /NO_FFT
    IF MAX(data_null.pcphmean) NE 0 THEN PLOTALL, time, data_null.pcphmean, 0., NAME=plotnull, TAG='PCPHMEAN', XTITLE='Elapsed time [s]', YTITLE='Mean phase per DIT [deg]', TITLE='', /BOLD, /EPS;, /NO_FFT
    IF MAX(null_tot_phot2)     NE 0 THEN PLOTALL, time[idx_ok], 100.*null_tot_phot2/phot_tot, 100.*null_err_phot2/phot_tot, NAME=plotnull, TAG='BIAS', XTITLE='Elapsed time [s]', YTITLE='Background null [%]', TITLE='', /BOLD, /EPS
    IF MAX(null_tot_phot2)     NE 0 THEN PLOTALL, null_tot_bckg2, null_tot_phot2, 0., NAME=plotnull, TAG='NULL2-VS-FLOOR2', XTITLE='Background floor [ADU]', YTITLE='Background bias [ADU]', TITLE='', /BOLD, /EPS, /SCATTER, /NO_FFT
      
    ; 4. Compute PHASECam telemetry AVG and RMS per OBs
    ; *************************************************
    
    ; Average phase over the last (x) ms
    AVGSDV, data_null.pcphmean, pcphmean_avg, pcphmean_rms, pcphmean_err, KAPPA=5
    ;IF MEAN(data_null.pcphmean-data_null.pcplsp) GT 5*pcphmean_rms THEN BEGIN
     ;MESSAGE, "Double-check this OB: measured phase doesn't match setpoint ( " + STRING(MEAN(data_null.pcphmean), FORMAT='(I0)') + " degrees instead of " + STRING(MEAN(data_null.pcplsp), FORMAT='(I0)') + ')', /CONTINUE
     ;IF lun GT 0 THEN PRINTF, lun, "Double-check this OB: measured phase doesn't match setpoint ( " + STRING(MEAN(data_null.pcphmean), FORMAT='(I0)') + " degrees instead of " + STRING(MEAN(data_null.pcplsp), FORMAT='(I0)') + ')'
    ;ENDIF 
             
    ; RMS phase over the last (x) ms (remove outliers)
    AVGSDV, data_null.pcphstd, pcphstd_avg, pcphstd_rms, pcphstd_err, KAPPA=5
    
    ; Scale to NOMIC (WARNING: hardcoded for now. These are the values used by Elwood's to compute PCPHMCOS and PCPHMSIN. We should use the same here!!)
    ; WARNING 2: pcphmean should include the dither pattern to compute correctly pcphmcos and pcphmsin below
    pcphmean = data_null.pcphmean*!dtor*2.2/11.1     
    
    ; Compute effective wavelength for PHASECam (read target spectrum and instrument throughput)
    pcwav = FXPAR(header, 'PCWAV', /NOCONTINUE)
    IF !err EQ -1 THEN BEGIN
      READ_TABLE, 'nodrs/input/phasecam_K.txt', lam_k, trans, FIRST=0, SEPARATOR=' '
      idx_k = WHERE(lam_k GT 1.95 and lam_k LT 2.45, n_lam)
      lam_k = lam_k[idx_k]
      trans = trans[idx_k]    
      GET_TGT, objname, tgt, DATABASE= pth.input_path + drs.database
      pick  = READ_PICKLES(pth.pickles_path + tgt.spectrum, LAM_MIN=1.95D-6, LAM_MAX=2.45D-6, NLAMBDA=n_lam, STANDARD=3)
      spec  = INTERPOL(pick[1,*], pick[0,*], lam_k*1D-6)
      pcwav = 1./(TOTAL((trans*spec)^2*1./lam_k)/TOTAL((trans*spec)^2))*1D-6
    ENDIF
    
    ; Subtract estimated null floor if requested (convert measured K-band phase in degrees to N-band phase in radians)
    ; For NSC, this information is used in the fit rather than subtracted from each null measurements
    phi_rms_dit = FLTARR(n_null) & phi_avg_dit = FLTARR(n_null) & pcphmcos = FLTARR(n_null) & pcphmsin = FLTARR(n_null)
    IF drs.null_mode EQ 2 THEN BEGIN 
      IF drs.null_cor EQ 1 THEN null_tot_phot -= phot_tot*(REFORM(data_null.pcphstd)*!dtor*pcwav/lam_cen)^2/2
      IF drs.null_cor GE 2 THEN BEGIN
        phi_rms_dit = REFORM(data_null.pcphstd)*!dtor*pcwav/lam_cen                                                            ; convert to radian and scale to NOMIC's wavelength
        IF TAG_EXIST(data_null, 'pcphmcos') AND date NE '230523' AND date NE '230525' THEN pcphmcos = data_null.pcphmcos*COS(pcphmean)+data_null.pcphmsin*SIN(pcphmean)  ; mean cos(phase-<phase>) over DIT + backward compatibility. For 230523 and 230525, PCPHMCOS and PCPHMSIN are wrong.
        IF TAG_EXIST(data_null, 'pcphmsin') AND date NE '230523' AND date NE '230525' THEN pcphmsin = data_null.pcphmsin*COS(pcphmean)-data_null.pcphmcos*SIN(pcphmean)  ; mean sin(phase-<phase>) over DIT + backward compatibility. For 230523 and 230525, PCPHMCOS and PCPHMSIN are wrong.
        IF drs.null_cor EQ 3 THEN MESSAGE, 'NULL_COR of 3 is depreciated'; phi_avg_dit = REFORM(data_null.pcphmean-pcphmean_avg)*!dtor*pcwav/lam_cen;)^3/3              ; Mean phase over DIT (not used anymnore)
      ENDIF
    ENDIF ELSE IF drs.null_cor GT 0 THEN null_tot_phot -= phot_tot*(REFORM(data_null.pcphstd)*!dtor*pcwav/lam_cen)^2/2
    
    ; Compute phase RMS at NOMIC's frequency (used as minimum phase jitter for NSC)
    pcphmean = (data_null.pcphmean-data_null.spdthpos)*!dtor*2.2/11.1 
    AVGSDV, pcphmean, pcphmean_avg, pcphmean_rms, pcphmean_err, KAPPA=5
    phi_rms_min = 0.5*pcphmean_rms ; use 0.5 to account for possible PHASECam noise (to be checked). Not reliable right now because of fringe jumps!
    
    ; Read dither pattern and scale it to NOMIC's wavelength
    spdthpos = data_null.spdthpos*!dtor*pcwav/lam_cen        ; phase in RADIAN at NOMIC's wavelength
    IF MAX(ABS(spdthpos)) GT 0 THEN PLOTALL, time, spdthpos, 0., NAME=plotnull, TAG='SPDTHPOS', XTITLE='Elapsed time [s]', YTITLE='OPD dither pattern [rad]', TITLE='', /BOLD, /EPS;, /NO_FFT
    IF MAX(ABS(pcphmcos)) GT 0 THEN PLOTALL, time, pcphmcos, 0., NAME=plotnull, TAG='PCPHMCOS', XTITLE='Elapsed time [s]', YTITLE='mean COS(PHASE) over DIT', TITLE='', /BOLD, /EPS, /NO_FFT
    IF MAX(ABS(pcphmsin)) GT 0 THEN PLOTALL, time, pcphmsin, 0., NAME=plotnull, TAG='PCPHMSIN', XTITLE='Elapsed time [s]', YTITLE='mean SIN(PHASE) over DIT', TITLE='', /BOLD, /EPS, /NO_FFT
  
    
    ; 5. Prepare sequence for NSC fit (not used anymore)
    ; *******************************
    
    ; Compute the true photometric variations 
    IF drs.null_mode EQ 2 THEN BEGIN
      ; Number of points to create the synthetic sequences
      n_gen = 1D5                                    
      ; Compute residual RMS (after subtraction of background contribution)
      IF bck_mode GE 0 THEN bck_rms_tmp = rms_bias ELSE bck_rms_tmp = bck_rms
      IF phot1_err GT bck_rms_tmp THEN phot1_nsc_rms = SQRT(phot1_err^2-bck_rms_tmp^2) ELSE phot1_nsc_rms = 0
      IF phot2_err GT bck_rms_tmp THEN phot2_nsc_rms = SQRT(phot2_err^2-bck_rms_tmp^2) ELSE phot2_nsc_rms = 0
      ; Create synthetic photometric sequences (Gaussian assumption) 
      ; WARNING. RANDOMN updates the seed!
      phot_tot_phot1 = phot1 + RANDOMN(noseed1, n_gen)*phot1_nsc_rms 
      phot_tot_phot2 = phot2 + RANDOMN(noseed2, n_gen)*phot2_nsc_rms 
      ;bck_tot_phot   = RANDOMN(noseed, n_gen)*bck_rms_tmp
    ENDIF
    
    ; 6. Compute null per OB
    ; **********************
    ratio = 1
    mode  = 0
    ; Statistical reduction if null mode is 2
    IF drs.null_mode NE 2 THEN BEGIN
      ; Classical reduction
      CASE drs.null_mode OF
        0: mode  = 1
        1: BEGIN
          ; Add noise to account for different brightness (manual at this point)
          ratio = drs.keep_ratio
          x_cor = SQRT((phot_tot/phot_faintcal)^2-1)
          IF bck_mode GE 0 THEN null_tot_phot += (null_tot_phot2-bck_bias_cor)*x_cor ELSE null_tot_phot += (null_tot_bckg-avg_bck)*x_cor ; Subtract mean value to only account for photon noise
        END
        ELSE: PRINT, 'Unknown null mode'
      ENDCASE
      ; Compute null (for star) and parse output strucuture similarly to NSC
      null_res = BIN_NULLDATA(null_tot_phot, null_err_phot, RATIO=ratio, MODE=mode)/phot_tot
      data[i_f].null_star.null_avg.nas = DOUBLE(null_res[0])
      data[i_f].null_star.null_avg.nas_err_low = DOUBLE(null_res[1])
      data[i_f].null_star.null_avg.nas_err_sup = DOUBLE(null_res[1])
    ENDIF ELSE BEGIN
      ; Statistical reduction. Depending on the background mode, we use either the background measured in a nearby empty region of the detector (bck_mode > 0)
      ; or use the same position but in the complemntary nod (bck_mode < 0).
      IF bck_mode GE 0 THEN dark = null_tot_phot2-bck_bias_cor ELSE dark = bck_tot_phot 
      arraydata = {DARK2: dark, A: phot_tot_phot1, B: phot_tot_phot2, AB: null_tot_phot, NAS_MIN: nas_min, PHI_RMS_MIN: phi_rms_min, PHI_AVG_T: phi_avg_dit, PHI_VAR_T: phi_rms_dit^2, PCPHMCOS_T: pcphmcos, PCPHMSIN_T: pcphmsin, SPDTHPOS: spdthpos}
      data[i_f].null_star = CALL_NSC(arraydata, CUBE_SIZE=drs.nsc_cube, N_BOOTSTRAP=drs.n_btstrp, NBINS_FACTOR=drs.nsc_bfac, NO_MULTI=no_multi, NULL_RANGE=drs.null_range,  FILE_ID=ob_id, WAV_EFF=lam_cen, BANDWIDTH=bdwdth,$
                                     OCC_MIN=drs.nsc_omin, VARIABLE=drs.nsc_bins, INFO=info, FILE_PATH=plot_path_nsc, /EPS, /PLOT)
    ENDELSE
        
    ; Save data in output structure
    ; *****************************
    
    ; Header information
    data[i_f].instrum  = STRTRIM(FXPAR(header, 'INSTRUME'))
    data[i_f].objname  = objname
    data[i_f].pid      = FXPAR(header, 'PID')
    data[i_f].flag     = STRTRIM(FXPAR(header, 'FLAG'))
    data[i_f].mjd_obs  = MEAN(data_null.mjd_obs, /DOUBLE)
    data[i_f].lbt_utc  = STRTRIM(FXPAR(header, 'LBT_UTC'))        ; at start of sequence
    data[i_f].lbt_lst  = STRTRIM(FXPAR(header, 'LBT_LST'))        ; at start of sequence
    data[i_f].lbt_ra   = STRTRIM(FXPAR(header, 'LBT_RA'))         ; at start of sequence
    data[i_f].lbt_dec  = STRTRIM(FXPAR(header, 'LBT_DEC'))        ; at start of sequence
    data[i_f].lbt_para = FXPAR(header, 'LBT_PARA')                ; mean of sequence
    data[i_f].lbt_alt  = FXPAR(header, 'LBT_ALT')                 ; mean of sequence
    data[i_f].lbt_az   = FXPAR(header, 'LBT_AZ')                  ; mean of sequence
    ;data[i_f].lbt_ha   = 0.                                      ; not yet in LBT header (temporary implementation)
    
    ; Detector parameters
    data[i_f].int_time = exptime
    data[i_f].acq_time = FXPAR(header, 'ACQTIME')
    data[i_f].lam_cen  = lam_cen
    data[i_f].bandwidth= bdwdth
    data[i_f].n_xpix   = FXPAR(header, 'N_XPIX')
    data[i_f].n_ypix   = FXPAR(header, 'N_YPIX')
    data[i_f].pixscale = pixscale
    data[i_f].xcen     = MEAN(data_null.xcen)
    data[i_f].ycen     = MEAN(data_null.ycen)
    data[i_f].smplmode = FXPAR(header, 'SMPLMODE')
    data[i_f].detmode  = FXPAR(header, 'DETMODE')
    data[i_f].pagain   = FXPAR(header, 'PAGAIN')
    data[i_f].pabandw  = FXPAR(header, 'PABANDW')
    data[i_f].detbias  = FXPAR(header, 'DETBIAS')
    data[i_f].eperadu  = FXPAR(header, 'EPERADU')
    data[i_f].n_coadd  = FXPAR(header, 'N_COADD')
    
    ; Filter keywords
    data[i_f].ubc_dxsp = FXPAR(header, 'UBC_DXSP')
    data[i_f].ubc_sxsp = FXPAR(header, 'UBC_SXSP')
    data[i_f].nic_fstp = FXPAR(header, 'NIC_FSTP')
    data[i_f].nic_beam = FXPAR(header, 'NIC_BEAM')
    data[i_f].lmir_fw1 = FXPAR(header, 'LMIR_FW1')
    data[i_f].lmir_fw2 = FXPAR(header, 'LMIR_FW2')
    data[i_f].lmir_fw3 = FXPAR(header, 'LMIR_FW3')
    data[i_f].lmir_fw4 = FXPAR(header, 'LMIR_FW4')
    data[i_f].nom_fw1  = FXPAR(header, 'NOM_FW1')
    data[i_f].nom_fw2  = FXPAR(header, 'NOM_FW2')
    data[i_f].nom_apw  = FXPAR(header, 'NOM_APW')
    data[i_f].pha_fw1  = FXPAR(header, 'PHA_FW1')
    data[i_f].pha_fw2  = FXPAR(header, 'PHA_FW2')
    data[i_f].pha_img  = FXPAR(header, 'PHA_IMG')
    data[i_f].nil_dic  = FXPAR(header, 'NIL_DIC')

    ; AO keywords
    data[i_f].daomode  = FXPAR(header, 'DAOMODE')
    data[i_f].daostrhl = FXPAR(header, 'DAOSTRHL')
    data[i_f].dcmodes  = FXPAR(header, 'DCMODES')
    data[i_f].dloopon  = FXPAR(header, 'DLOOPON')
    data[i_f].dloopgn  = FXPAR(header, 'DLOOPGN')
    data[i_f].dwfscfrq = FXPAR(header, 'DWFSCFRQ')
    data[i_f].dwfscbin = FXPAR(header, 'DWFSCBIN')
    data[i_f].saomode  = FXPAR(header, 'SAOMODE')
    data[i_f].saostrhl = FXPAR(header, 'SAOSTRHL')
    data[i_f].scmodes  = FXPAR(header, 'SCMODES')
    data[i_f].sloopon  = FXPAR(header, 'SLOOPON')
    data[i_f].sloopgn  = FXPAR(header, 'SLOOPGN')
    data[i_f].swfscfrq = FXPAR(header, 'SWFSCFRQ')
    data[i_f].swfscbin = FXPAR(header, 'SWFSCBIN')

    ; Phasecam keywords
    data[i_f].pcclosed   = FXPAR(header, 'PCCLOSED')
    data[i_f].plc_wav    = pcwav
    data[i_f].plc_spec   = pcspec
    data[i_f].pcplsp     = MEAN(data_null.pcplsp)
    data[i_f].pctipsp    = MEAN(data_null.pctipsp)
    data[i_f].pctltsp    = MEAN(data_null.pctltsp)
    data[i_f].spc_pist   = FXPAR(header, 'SPCPIST')
    data[i_f].spc_az     = FXPAR(header, 'SPCAZ')
    data[i_f].spc_el     = FXPAR(header, 'SPCEL')
    data[i_f].fpc_pistm  = MEAN(data_null.fpcpistm)
    data[i_f].fpc_pists  = MEAN(data_null.fpcpists)
    data[i_f].fpc_azm    = MEAN(data_null.fpcazm)
    data[i_f].fpc_azs    = MEAN(data_null.fpcazs)
    data[i_f].fpc_elm    = MEAN(data_null.fpcelm)
    data[i_f].fpc_els    = MEAN(data_null.fpcels)
    data[i_f].pcphmean   = pcphmean_avg & data[i_f].pcphmean_err = pcphmean_err
    data[i_f].pcphstd    = pcphstd_avg  & data[i_f].pcphstd_err  = pcphstd_err
    data[i_f].pcmsnr     = MEAN(data_null.pcmsnr)   
    
    ; Weather parameters
    data[i_f].seeing   = FXPAR(header, 'SEEING')
    data[i_f].smttau   = FXPAR(header, 'SMTTAU')
    data[i_f].lbttemp  = FXPAR(header, 'LBTTEMP')
    data[i_f].winddir  = FXPAR(header, 'WINDDIR')
    data[i_f].windspd  = FXPAR(header, 'WINDSPD')
    
    ; Nulls
    data[i_f].null_avg  = null_avg    ; Diagnostic info
    data[i_f].null_rms  = null_rms    ; Diagnostic info
           
    ; Photometry
    data[i_f].photdx_avg = phot1
    data[i_f].photsx_avg = phot2
    data[i_f].photdx_rms = phot1_err
    data[i_f].photsx_rms = phot2_err
    data[i_f].photdx_snr = phot1/phot1_err
    data[i_f].photsx_snr = phot2/phot2_err
    data[i_f].phot_avg   = phot_tot
    data[i_f].phot_err   = phot_err
    data[i_f].int_err    = int_err
    
    ; Background
    data[i_f].bckg_avg   = bck_avg/phot_tot       ; AVG background per photometric aperture   
    data[i_f].bckg_rms   = bck_rms/phot_tot       ; RMS background per photometric aperture
    data[i_f].bckg_err   = bck_rms_mean/phot_tot  ; RMS background per photometric aperture 
    data[i_f].bias_unc   = null_bckg[0]           ; Uncorrected background bias per photometric aperture 
    data[i_f].bias_cor   = null_bckg[1]           ; Empty-region background bias per photometric aperture 
    data[i_f].bias_err   = null_bckg[2]           ; Empty-region background bias error per photometric aperture 
    data[i_f].bckg_meas  = avg_bck                ; AVG background per pixel (in the background region surrounding the photometric aperture)          
    data[i_f].bckg_emeas = rms_bck                ; RMS background per pixel (in the background region surrounding the photometric aperture) 
        
    ; Reduction paramaters
    data[i_f].file_id   = ob_id
    data[i_f].ob_id     = ob_id                ; currently OB id and file if are the same
    data[i_f].pt_id     = pt_id               
    data[i_f].nod_id    = data_null[0].nod_id
    data[i_f].nfr_in    = n_null0
    data[i_f].nfr_rej   = n_null0-n_null;+n_rej  ; now n_rej includes the backgrounds frames obtained between nods...so no
    data[i_f].nfr_ob    = FLOOR(ratio*n_null)
    data[i_f].nod_frq   = FXPAR(header, 'NOD_FRQ')
    data[i_f].n_frbck   = FXPAR(header, 'N_FRBCK')    
    data[i_f].rms_opd   = SQRT(data[i_f].null_err)*(0.5*data[i_f].lam_cen)  ; first order approximation
    data[i_f].rms_tt    = [rms_tt1,rms_tt2]
    data[i_f].aper_rad  = apr_rad
    data[i_f].bck_irad  = bck_irad
    data[i_f].bck_orad  = bck_orad
    data[i_f].fit_mode  = fit_mode
    data[i_f].flx_mode  = flx_mode
    data[i_f].fra_mode  = fra_mode
    data[i_f].img_mode  = img_mode
    data[i_f].ob_mode   = ob_mode
    data[i_f].bck_mode  = bck_mode
    data[i_f].bfl_mode  = bfl_mode
    data[i_f].precision = precisio
    
    ; Compute quality flag. Current defintion follwoing discussion with RMG is 1: no suitable for science, 2: poor quality due to seeing or PWV, 3: is good quality.
    ; Per my experience, I use the follwoing definition:
    ;  - 1 : seeing is above 1.5 or PWV above 6mm or more than 50% of frames rejected
    ;  - 2 : seeing is above 1.2 or PWV above 4.5mm or more than 25% of frames rejected
    ;  - 3 : otherwise
    qua_flg = 3
    IF data[i_f].seeing GT 1.2 OR data[i_f].smttau GT 4.5 OR data[i_f].nfr_rej/n_null0 GT 0.25 THEN qua_flg = 2
    IF data[i_f].seeing GT 1.5 OR data[i_f].smttau GT 6   OR data[i_f].nfr_rej/n_null0 GT 0.50 THEN qua_flg = 1    
    data[i_f].qua_flg = qua_flg
    
    ; Save OB file
    IF NOT KEYWORD_SET(no_save) and drs.null_mode NE 1 THEN BEGIN          
      ; Output file (new file if not /RENEW)
      ob_file0 = sav_int_path + 'null_ob' + STRING(ob_id, FORMAT='(I03)') + '.sav'
      ; If exist, archiv
      IF FILE_TEST(ob_file0) THEN BEGIN
        i_red   = 1
        ob_file = sav_int_path + 'null_ob' + STRING(ob_id, FORMAT='(I03)') + '_v1.sav'
        WHILE FILE_TEST(ob_file) EQ 1 DO BEGIN
          ob_file = sav_int_path + 'null_ob' + STRING(ob_id, FORMAT='(I03)') + '_v' + STRING(++i_red, FORMAT='(I0)') + '.sav'
        ENDWHILE
        FILE_MOVE, ob_file0, ob_file
      ENDIF ELSE i_red = 0
      ; If /RENEW, erase the old one (and use previous ob_file0 defined above if doesn't exist)
      IF KEYWORD_SET(RENEW) THEN BEGIN
        files_ob = FILE_SEARCH(data_path,'null_ob' + STRING(ob_id, FORMAT='(I03)') + '*.sav', COUNT=n_files_ob)
        ; Loop over files of this OB until matching paramaters are found (or not)
        FOR i_fob = 0, n_files_ob-1 DO BEGIN
          ; Restore old file
          RESTORE, files_ob[i_fob]
          IF TAG_EXIST(data_ob, 'bfl_mode') THEN BEGIN
            ; Compare to running/current parameters
            IF drs_ob.null_mode EQ drs.null_mode AND drs_ob.keep_ratio  EQ drs.keep_ratio  AND drs_ob.min_fr     EQ drs.min_fr      AND drs_ob.n_btstrp    EQ drs.n_btstrp    AND drs_ob.nsc_bfac  EQ drs.nsc_bfac AND $
              drs_ob.nsc_bins   EQ drs.nsc_bins  AND drs_ob.nsc_cube[0] EQ drs.nsc_cube[0] AND drs_ob.nsc_cube[1] EQ drs.nsc_cube[1] AND drs_ob.nsc_cube[2] EQ drs.nsc_cube[2] AND drs_ob.nsc_omin  EQ drs.nsc_omin AND $
              drs_ob.null_cor   EQ drs.null_cor  AND drs_ob.null_rad    EQ drs.null_rad    AND drs_ob.null_lim[0] EQ drs.null_lim[0] AND drs_ob.null_lim[1] EQ drs.null_lim[1] AND drs_ob.null_range[0] EQ drs.null_range[0] AND $
              drs_ob.null_range[1] EQ drs.null_range[1] AND drs.n_frob EQ drs.n_frob   AND $ ; end of drs parameters
              data_ob.objname  EQ objname       AND data_ob.lam_cen    EQ lam_cen         AND data_ob.bck_mode   EQ bck_mode        AND data_ob.bfl_mode   EQ bfl_mode        AND data_ob.fit_mode   EQ fit_mode        AND data_ob.img_mode EQ img_mode     AND $
              data_ob.ob_mode  EQ ob_mode       AND data_ob.pixscale   EQ pixscale        AND data_ob.precision  EQ precisio        AND data_ob.bck_irad   EQ bck_irad        AND data_ob.bck_orad EQ bck_orad THEN BEGIN                                                                              ; end of paramaters on L1 images
              ob_file0 = files_ob[i_fob]
            ENDIF
          ENDIF
        ENDFOR
      ENDIF      
      ; Save file
      data_ob = data[i_f]
      drs_ob  = drs
      SAVE, drs_ob, data_ob, FILENAME=ob_file0
    ENDIF
    
    ; Jump point if restored OB
    RESTORED_OB:
    
    ; Use NSC mode results if NSC_MODE is 1 (only one is saved in the L1 summary file but they can be reproduced without recomputing everything)
    ; Add error on photometry (not included in NSC) and error on background subtraction
    err_phot = ABS(data[i_f].null_avg)*(data[i_f].phot_err/data[i_f].phot_avg)  ; error on photometry not included in NSC (correlated between all OBs of the same pointing!)
    err_bckg = data[i_f].bckg_err                                          ; error on background subtraction from complementary nod
    err_trm  = SQRT(err_phot^2+err_bckg^2)    
    IF drs.null_mode EQ 2 THEN BEGIN
      CASE drs.nsc_mode OF
        0: null_data = data[i_f].null_star.null_opt   ; use optimum value
        1: null_data = data[i_f].null_star.null_bay   ; use bayesian value
        2: MESSAGE, 'Bootstrap computation is depreciated. Use NSC_MODE=1 for best results.' ;null_data = data[i_f].null_star.null_avg   ; use average of bootsrap samples
        3: MESSAGE, 'Bootstrap computation is depreciated. Use NSC_MODE=1 for best results.' ;null_data = data[i_f].null_star.null_mod   ; use mode of bootstrap sampls
        ELSE: MESSAGE, 'Unknown NSC mode. Must be 0, 1, 2, or 3.'
      ENDCASE
    ENDIF ELSE null_data = data[i_f].null_star.null_avg
    ; Now parse to null results
    err_nsc             = null_data.nas_err_sup > null_data.nas_err_low                          ; Take MAX of error LOW and error SUP
    data[i_f].null_meas = null_data.nas   & data[i_f].null_err      = SQRT(err_nsc^2+err_trm^2)  ; Add error on photometry (not included in NSC) and error on background subtraction
    data[i_f].nsc_phavg = null_data.mu    & data[i_f].nsc_phavg_err = null_data.mu_err_sup > null_data.mu_err_low 
    data[i_f].nsc_phrms = null_data.sig   & data[i_f].nsc_phrms_err = null_data.sig_err_sup > null_data.sig_err_sup
    data[i_f].nsc_kdrk  = null_data.kdrk  & data[i_f].nsc_kdrk_err  = null_data.kdrk_err_sup > null_data.kdrk_err_sup
    data[i_f].nsc_chi2  = null_data.chi2  & data[i_f].nsc_chi2_err  = null_data.chi2_err_sup
    data[i_f].null_snr  = data[i_f].null_meas/data[i_f].null_err
    data[i_f].null_offset = data[i_f].null_meas - data[i_f].null_star.null_opt.nas
    
    ; Compute expected error terms
    IF data[i_f].bck_mode GE 0 THEN null_err_phot = SQRT(err_bckg^2 + data[i_f].bckg_ebias^2) ELSE null_err_phot = SQRT(2)*err_bckg
    null_err_phase = SQRT(4*(data[i_f].nsc_phavg)^2*(data[i_f].nsc_phrms)^2 + 2*(data[i_f].nsc_phrms)^4)/4
    
    ; Print onfo
    IF info GT 0 THEN BEGIN
      PRINT, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(data[i_f].objname, FORMAT='(A8)'), '  ', STRING(data[i_f].flag, FORMAT='(A3)'), '  ', STRING(data[i_f].lbt_utc, FORMAT='(A11)'), '  ', STRING(data[i_f].lbt_alt, FORMAT='(F6.2)'), '  ', STRING(data[i_f].lbt_para, FORMAT='(F7.2)'), '  ', STRING(data[i_f].seeing ,FORMAT='(F4.2)'), '  ', STRING(data[i_f].smttau ,FORMAT='(F4.2)'), $
        '  ', STRING(1D+3*data[i_f].int_time, FORMAT='(F5.1)'), '  ', STRING(1D+6*data[i_f].lam_cen, FORMAT='(F5.1)'), '  ', STRING(1D+6*data[i_f].bandwidth, FORMAT='(F4.1)'), ' || ', STRING(data[i_f].photdx_snr, FORMAT='(F5.1)'), ' ', STRING(data[i_f].photsx_snr, FORMAT='(F5.1)'),$
        '  ', STRING(data[i_f].int_err, FORMAT='(F7.4)'), '  ', STRING(data[i_f].rms_tt[0], FORMAT='(F5.1)'), ' ', STRING(data[i_f].rms_tt[1], FORMAT='(F5.1)'), '  ', STRING(100.*data[i_f].photdx_rms/data[i_f].phot_avg, FORMAT='(F4.2)'), '  ', STRING(100.*data[i_f].photsx_rms/data[i_f].phot_avg, FORMAT='(F4.2)'), $
        '  ', STRING(100.*data[i_f].bckg_rms, FORMAT='(F5.2)'), '  ', STRING(100.*data[i_f].bias_unc, FORMAT='(F5.1)'), '  ', STRING(100.*data[i_f].bias_cor, FORMAT='(F5.2)'), ' || ', STRING(100.*data[i_f].null_meas, FORMAT='(F5.2)'), '  ', STRING(100.*data[i_f].null_offset, FORMAT='(F5.2)'), ' || ', STRING(100.*err_phot, FORMAT='(F5.2)'), '  ', STRING(100.*err_bckg, FORMAT='(F5.2)'), '  ', STRING(100.*err_nsc, FORMAT='(F5.2)'), $
        ' | ', STRING(100.*data[i_f].null_err, FORMAT='(F5.2)'), ' || ',  STRING(100.*null_err_phot, FORMAT='(F5.2)'), ' ',  STRING(100.*null_err_phase, FORMAT='(F5.2)'), ' ||', STRING(1D6*data[i_f].nsc_phavg*lam_cen/(2*!dpi), FORMAT='(F6.2)'),'  ', STRING(1D6*data[i_f].nsc_phrms*lam_cen/(2*!dpi), FORMAT='(F5.2)'),'  ', STRING(data[i_f].nsc_kdrk, FORMAT='(F4.2)'), '  ', STRING(data[i_f].nsc_chi2, FORMAT='(F6.2)'), '   ', STRING(data[i_f].nfr_rej, FORMAT='(I04)'), '  ', STRING(data[i_f].nfr_ob, FORMAT='(I04)')
      IF lun GT 0 THEN $
        PRINTF, lun, STRING(ob_id, FORMAT='(I03)'), '  ', STRING(data[i_f].objname, FORMAT='(A8)'), '  ', STRING(data[i_f].flag, FORMAT='(A3)'), '  ', STRING(data[i_f].lbt_utc, FORMAT='(A11)'), '  ', STRING(data[i_f].lbt_alt, FORMAT='(F6.2)'), '  ', STRING(data[i_f].lbt_para, FORMAT='(F7.2)'), '  ', STRING(data[i_f].seeing ,FORMAT='(F4.2)'), '  ', STRING(data[i_f].smttau ,FORMAT='(F4.2)'), $
                      '  ', STRING(1D+3*data[i_f].int_time, FORMAT='(F5.1)'), '  ', STRING(1D+6*data[i_f].lam_cen, FORMAT='(F5.1)'), '  ', STRING(1D+6*data[i_f].bandwidth, FORMAT='(F4.1)'), ' || ', STRING(data[i_f].photdx_snr, FORMAT='(F5.1)'), ' ', STRING(data[i_f].photsx_snr, FORMAT='(F5.1)'),$
                      '  ', STRING(data[i_f].int_err, FORMAT='(F7.4)'), '  ', STRING(data[i_f].rms_tt[0], FORMAT='(F5.1)'), ' ', STRING(data[i_f].rms_tt[1], FORMAT='(F5.1)'), '  ', STRING(100.*data[i_f].photdx_rms/data[i_f].phot_avg, FORMAT='(F4.2)'), '  ', STRING(100.*data[i_f].photsx_rms/data[i_f].phot_avg, FORMAT='(F4.2)'), $
                      '  ', STRING(100.*data[i_f].bckg_rms, FORMAT='(F5.2)'), '  ', STRING(100.*data[i_f].bias_unc, FORMAT='(F5.1)'), '  ', STRING(100.*data[i_f].bias_cor, FORMAT='(F5.2)'), ' || ', STRING(100.*data[i_f].null_meas, FORMAT='(F5.2)'), '  ', STRING(100.*data[i_f].null_offset, FORMAT='(F5.2)'), ' || ', STRING(100.*err_phot, FORMAT='(F5.2)'), '  ', STRING(100.*err_bckg, FORMAT='(F5.2)'), '  ', STRING(100.*err_nsc, FORMAT='(F5.2)'), $
                      ' | ', STRING(100.*data[i_f].null_err, FORMAT='(F5.2)'), ' || ',  STRING(100.*null_err_phot, FORMAT='(F5.2)'), ' ',  STRING(100.*null_err_phase, FORMAT='(F5.2)'), ' ||', STRING(1D6*data[i_f].nsc_phavg*lam_cen/(2*!dpi), FORMAT='(F6.2)'),'  ', STRING(1D6*data[i_f].nsc_phrms*lam_cen/(2*!dpi), FORMAT='(F5.2)'), '  ', STRING(data[i_f].nsc_kdrk, FORMAT='(F4.2)'),'  ', STRING(data[i_f].nsc_chi2, FORMAT='(F6.2)'), '   ', STRING(data[i_f].nfr_rej, FORMAT='(I04)'), '  ', STRING(data[i_f].nfr_ob, FORMAT='(I04)')
    ENDIF
   
    ; Jump point if no nulling data in this OB
    SKIP_OB:
ENDFOR

; Remove skipped frames from data
idx_ok = WHERE(data.mjd_obs NE 0, n_ok)
IF n_ok GT 0 THEN data = data[idx_ok]

; Save L1 file
IF NOT KEYWORD_SET(no_save) THEN BEGIN  
  ; Create directories if non existent
  sav_path = pth.l1fits_path + pth.sep + drs.date_obs + drs.dir_label + pth.sep 
  IF NOT FILE_TEST(sav_path) THEN FILE_MKDIR, sav_path & SPAWN, 'chmod 775 ' + sav_path        ; Create directory if it does not exist
  
  ; If exist, archive it
  outfile = sav_path + 'UT' + drs.date_obs + '.fits'
  savfile = sav_path + 'UT' + drs.date_obs + '.sav'
  IF FILE_TEST(outfile) THEN BEGIN
    i_red = 1
    exist = 1
    WHILE exist EQ 1 DO BEGIN
      oldfile = sav_path + 'UT' + drs.date_obs + '_v' + STRING(i_red, FORMAT='(I0)') + '.fits'
      IF NOT FILE_TEST(oldfile) THEN exist = 0 ELSE i_red = i_red + 1
    ENDWHILE
    IF FILE_TEST(outfile) THEN FILE_MOVE, outfile, sav_path + 'UT' + drs.date_obs + '_v' + STRING(i_red, FORMAT='(I0)') + '.fits'
    IF FILE_TEST(savfile) THEN FILE_MOVE, savfile, sav_path + 'UT' + drs.date_obs + '_v' + STRING(i_red, FORMAT='(I0)') + '.sav'
  ENDIF ELSE i_red = 0
  
  ; Save data table
  LBTI_SAVEL1SUM, data, OUTFILE=outfile, FILE_VERSION=i_red+1
  
  ; Also save data as an IDL file
  SAVE, data, FILENAME= savfile
ENDIF
END