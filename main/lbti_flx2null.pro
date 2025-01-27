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
;   Version 6.2,  27-AUG-2024, DD: Moved up the computation of the mean background and save the corrected background as filtered file
;   Version 6.3,  27-JAN-2025, DD: Corrected syntax of bckg_bias error term (used only for bckg_mode greater than 1)

pro LBTI_FLX2NULL, date, ob_idx = ob_idx, info = info, log_file = log_file, no_multi = no_multi, no_save = no_save, plot = plot, renew = renew
  compile_opt idl2

  ; Global parameter
  common GLOBAL, prm, cnf, wav, tgt, pth, drs
  on_error, 0

  ; Check keywords
  if not keyword_set(info) then info = 0
  if not keyword_set(plot) then plot = 0 else if plot ge 2 then display = 1
  if not keyword_set(log_file) then lun = -1 else lun = log_file

  ; Define outlier rejection limit
  sig_out = 5

  ; If aperture radius is set, redefine dir_label
  if max(drs.null_rad) ne 0 and not strmatch(drs.dir_label, '*_APR') then drs.dir_label = drs.dir_label + '_APR'

  ; Retrieve data files
  data_path = pth.l1Fits_path + drs.date_obs + drs.dir_label + pth.sep
  null_files = file_search(data_path, '*NULL.fits', count = n_files)

  ; Create directories if non existent
  sav_fil_path = pth.l1Fits_path + pth.sep + drs.date_obs + drs.dir_label + pth.sep + 'filtered' + pth.sep
  if not file_test(sav_fil_path) then file_mkdir, sav_fil_path
  spawn, 'chmod -R 775 ' + sav_fil_path
  sav_int_path = pth.l1Fits_path + pth.sep + drs.date_obs + drs.dir_label + pth.sep + 'intermediate' + pth.sep
  if not file_test(sav_int_path) then file_mkdir, sav_int_path
  spawn, 'chmod -R 775 ' + sav_int_path

  ; Return if no nulling files
  if n_files eq 0 then begin
    message, 'No NULL files to process', /continue
    RETURN
  endif

  ; Print info to screen
  if info gt 0 then begin
    ; Basic info
    print, ''
    print, 'Now processing ' + string(n_files, format = '(I0)') + ' OBs'
    if not keyword_set(no_multi) and drs.null_mode eq 2 then print, 'Multi-threading on ' + string(2 ^ (floor(alog(!cpu.hw_ncpu - 1) / alog(2))), format = '(I0)') + ' cores'
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
    print, 'J: Wavelength [microns]'
    print, 'K: Bandwidths [microns]'
    print, 'L: SNR on photometry (DX and SX)'
    print, 'M: Intensity mismatch (in %)'
    print, 'N: Tip/tilt error (DX and SX, in mas)'
    print, 'O: RMS of photometry relative to constructive peak (DX and SX, in %)'
    print, 'P: RMS of background relative to constructive peak (in %)'
    print, 'Q: Uncorrected background bias relative to constructive peak (in %)'
    print, 'R: Estimated (from empty region) background bias relative to constructive peak (in %)'
    print, 'S: Null (in %)'
    print, 'T: Null offset between optimum NSC value and bootstrapped or bayesian value (in %)'
    print, 'U: Null uncertainty due to error on constructive peak (in %)'
    print, 'V: Null uncertainty due to external error on background floor (in %)'
    print, 'W: Null uncertainty from NSC fit (in %)'
    print, 'X: Total null uncertainty (in %)'
    print, 'Y: Expected null uncertainty due to photometric errors (in %)'
    print, 'Z: Expected null uncertainty due to phase variation (i.e., sqrt(4mu^2.sig^2 + 2*sig^4)/4, in %) TO BE CHECKED'
    print, 'a: Best-fit mean PHASE from NSC (in microns)'
    print, 'b: Best-fit PHASE jitter from NSC (in microns)'
    print, 'c: Best-fit background RMS multiplication factor from NSC'
    print, 'd: Reduced chi2 from NSC'
    print, 'e: Number of rejected null frames'
    print, 'f: Preserved number of frames'
    print, ' '
    print, ' A ', '  ', '    B   ', '  ', ' C ', '  ', '     D     ', '  ', '   E  ', '  ', '    F  ', '  ', '  G ', '  ', '  H ', '  ', '  I  ', '  ', '   J ', '  ', '  K  ', '  ', '      L     ', '  ', '   M   ', '  ', '     N     ', '  ', '    O     ', '  ', '  P  ', '  ', '  Q  ', '  ', '  R  ', '  ', '    S  ', ' ', '   T  ', '  ', '    U ', '  ', '   V  ', ' ', '   W  ', '  ', '   X  ', '   ', '   Y  ', '  ', '  Z  ', '  ', '   a  ', '  ', '   b  ', '  ', '  c  ', '  ', '  d  ', '  ', '  e  ', '  ', '  f  '
    print, '---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    if lun gt 0 then begin
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
      printf, lun, 'J: Wavelength [microns]'
      printf, lun, 'K: Bandwidths [microns]'
      printf, lun, 'L: SNR on photometry (DX and SX)'
      printf, lun, 'M: Intensity mismatch (in %)'
      printf, lun, 'N: Tip/tilt error (DX and SX, in mas)'
      printf, lun, 'O: RMS of photometry relative to constructive peak (DX and SX, in %)'
      printf, lun, 'P: RMS of background relative to constructive peak (in %)'
      printf, lun, 'Q: Uncorrected background bias relative to constructive peak (in %)'
      printf, lun, 'R: Estimated (from empty region) background bias relative to constructive peak (in %)'
      printf, lun, 'S: Null (in %)'
      printf, lun, 'T: Null offset between optimum NSC value and bootstrapped or bayesian value (in %)'
      printf, lun, 'U: Null uncertainty due to error on constructive peak (in %)'
      printf, lun, 'V: Null uncertainty due to external error on background floor (in %)'
      printf, lun, 'W: Null uncertainty from NSC fit (in %)'
      printf, lun, 'X: Total null uncertainty (in %)'
      printf, lun, 'Y: Expected null uncertainty due to photometric errors (in %)'
      printf, lun, 'Z: Expected null uncertainty due to phase variation (i.e., sqrt(4mu^2.sig^2 + 2*sig^4)/4, in %) TO BE CHECKED'
      printf, lun, 'a: Best-fit mean PHASE from NSC (in microns)'
      printf, lun, 'b: Best-fit PHASE jitter from NSC (in microns)'
      printf, lun, 'c: Best-fit background RMS multiplication factor from NSC'
      printf, lun, 'd: Reduced chi2 from NSC'
      printf, lun, 'e: Number of rejected null frames'
      printf, lun, 'f: Preserved number of frames'
      printf, lun, ' '
      printf, lun, ' A ', '  ', '    B   ', '  ', ' C ', '  ', '     D     ', '  ', '   E  ', '  ', '    F  ', '  ', '  G ', '  ', '  H ', '  ', '  I  ', '  ', '   J ', '  ', '  K  ', '  ', '      L     ', '  ', '   M   ', '  ', '     N     ', '  ', '    O     ', '  ', '  P  ', '  ', '  Q  ', '  ', '  R  ', '  ', '    S  ', ' ', '   T  ', '  ', '    U ', '  ', '   V  ', ' ', '   W  ', '  ', '   X  ', '   ', '   Y  ', '  ', '  Z  ', '  ', '   a  ', '  ', '   b  ', '  ', '  c  ', '  ', '  d  ', '  ', '  e  ', '  ', '  f  '
      printf, lun, '---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    endif
  endif

  ; Initialize output data structure
  DEFINE_NULLDATA, data

  ; Read PHASECam tranmission file
  READ_TABLE, 'nodrs/input/phasecam_K.txt', lam_k, trans, first = 0, separator = ' '
  idx_k = where(lam_k gt 1.95 and lam_k lt 2.45, n_lam)
  lam_k = lam_k[idx_k]
  trans = trans[idx_k]

  ; Prepare plot paths
  plot_path = pth.result_path + 'null' + pth.sep + drs.date_obs + drs.dir_label + pth.sep ; Main path for NULL plots
  diag_path = pth.result_path + 'diagnostic' + pth.sep + 'phasecam' + pth.sep + drs.date_obs + pth.sep ; Main path for DIAGNOSTIC plots
  if not file_test(plot_path) then file_mkdir, plot_path
  spawn, 'chmod 775 ' + plot_path ; Create directory if it does not exist
  if not file_test(diag_path) then file_mkdir, diag_path
  spawn, 'chmod 775 ' + diag_path ; Create directory if it does not exist

  ; Replicate output structure to the number of maximum lines
  data = replicate(data, n_files)

  ; Jump point to restart reduction using a different photomety (only used when null_mode = 1)
  restart:

  ; Loop over the OBs
  for i_f = 0, n_files - 1 do begin
    ; 1. Read null file
    ; *****************
    header = HEADFITS(null_files[i_f])
    data_null = MRDFITS(null_files[i_f], 1, /silent)

    ; Make sure this OB is in the user-defined range
    ob_id = FXPAR(header, 'OB_ID', datatype = 0, /nocontinue)
    if keyword_set(ob_idx) then if ob_id lt ob_idx[0] or ob_id gt ob_idx[1] then goto, skip_ob

    ; Read useful information
    objname = strcompress(strtrim(FXPAR(header, 'OBJNAME', /nocontinue)), /remove_all)
    exptime = FXPAR(header, 'EXPTIME', /nocontinue)
    lam_cen = FXPAR(header, 'WAVELENG', /nocontinue)
    bdwdth = FXPAR(header, 'BANDWIDT', /nocontinue)
    bck_mode = FXPAR(header, 'BCK_MODE', /nocontinue)
    bfl_mode = FXPAR(header, 'BFL_MODE', /nocontinue)
    fit_mode = FXPAR(header, 'FIT_MODE', /nocontinue)
    flx_mode = FXPAR(header, 'FLX_MODE', /nocontinue)
    fra_mode = FXPAR(header, 'FRA_MODE', /nocontinue)
    img_mode = FXPAR(header, 'IMG_MODE', /nocontinue)
    ob_mode = FXPAR(header, 'OB_MODE', /nocontinue)
    pixscale = FXPAR(header, 'PIXSCALE', /nocontinue)
    precisio = FXPAR(header, 'PRECISIO', /nocontinue)
    pt_id = FXPAR(header, 'PT_ID', /nocontinue)
    flag = FXPAR(header, 'FLAG', /nocontinue)
    n_rej = FXPAR(header, 'N_REJ', /nocontinue) > 0
    pcspec = FXPAR(header, 'TGT_SPEC', /nocontinue)
    calfor = FXPAR(header, 'TGT_CAL4', /nocontinue)

    ; Skip OB if not enough frames
    n_null0 = n_elements(data_null.flx_tot[0])
    if n_null0 lt drs.min_fr then begin
      if info gt 0 then print, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null measurements in OB (only ' + string(n_null0, format = '(I0)') + ')'
      if lun gt 0 then printf, lun, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null measurements in OB (only ' + string(n_null0, format = '(I0)') + ')'
      goto, skip_ob
    endif

    ; Compute effective wavelength for PHASECam (read target spectrum and instrument throughput)
    pick = READ_PICKLES(pth.pickles_path + pcspec, lam_min = 1.95d-6, lam_max = 2.45d-6, nlambda = n_lam, standard = 3)
    spec = interpol(pick[1, *], pick[0, *], lam_k * 1d-6)
    pcwav = 1. / (total((trans * spec) ^ 2 * 1. / lam_k) / total((trans * spec) ^ 2)) * 1d-6

    ; Compute EEID to determine the aperture size (use the same for SCI and CAL)
    if drs.null_rad eq 0 then begin
      if strtrim(flag, 2) eq 'CAL' then GET_TGT, calfor, tgt_sci, database = pth.input_path + drs.database else GET_TGT, objname, tgt_sci, database = pth.input_path + drs.database
      eeid_as = COMPUTE_EEID(0.5 * tgt_sci.ldm, tgt_sci.dist, tgt_sci.temp)
      aper_rad = round(eeid_as / pixscale) > 8
    endif else aper_rad = drs.null_rad

    ; Find closest aperture radius. If drs.aper_rad is 0, it means that we have to use the EEID.
    ; Use the drs.aper_rad otherwise
    n_apr = n_elements(data_null[0].flx_tot)
    apr_all = FXPAR(header, 'APERRAD' + string(0, format = '(I0)'))
    for i_r = 1, n_apr - 1 do apr_all = [apr_all, FXPAR(header, 'APERRAD' + string(i_r, format = '(I0)'))]
    tmp = min(abs(apr_all - aper_rad), i_aper)
    apr_rad = apr_all[i_aper]
    napr_pix = !dpi * apr_rad ^ 2 ; Number of pixels in photometric aperture and in annulus
    bck_irad = FXPAR(header, 'BCKIRAD' + string(i_aper, format = '(I0)'))
    bck_orad = FXPAR(header, 'BCKORAD' + string(i_aper, format = '(I0)'))
    nsky_pix = FXPAR(header, 'NSKYPIX' + string(i_aper, format = '(I0)'))
    if !err eq -1 then nsky_pix = !dpi * (bck_orad ^ 2 - bck_irad ^ 2) ; this was used for a test but not needed after all

    ; Restore OB if already reduced with the same parameters (and saved)
    ; Not working for null_mode = 1 because the reduction neeeds to know the photometry of the faintest star and it is only computed later
    ; But not a big deal because null_mode = 1 is not as slow as NSC
    if drs.null_mode ne 1 and not keyword_set(renew) then begin
      files_ob = file_search(data_path, 'null_ob' + string(ob_id, format = '(I03)') + '*.sav', count = n_files_ob)
      ; Loop over files of this OB until matching paramaters are found (or not)
      for i_fob = 0, n_files_ob - 1 do begin
        ; Restore old file
        restore, files_ob[i_fob]
        ; IF no BFL_MODE, consider as an obsolote file and keep looking
        if TAG_EXIST(data_ob, 'bfl_mode') then begin
          ; Compare to running/current parameters
          if drs_ob.null_mode eq drs.null_mode and drs_ob.keep_ratio eq drs.keep_ratio and drs_ob.min_fr eq drs.min_fr and drs_ob.n_btstrp eq drs.n_btstrp and drs_ob.nsc_bfac eq drs.nsc_bfac and $
            drs_ob.nsc_bins eq drs.nsc_bins and drs_ob.nsc_cube[0] eq drs.nsc_cube[0] and drs_ob.nsc_cube[1] eq drs.nsc_cube[1] and drs_ob.nsc_cube[2] eq drs.nsc_cube[2] and drs_ob.nsc_cube[3] eq drs.nsc_cube[3] and $
            drs_ob.nsc_cube[4] eq drs.nsc_cube[4] and drs_ob.nsc_omin eq drs.nsc_omin and drs_ob.null_cor eq drs.null_cor and drs_ob.null_rad eq drs.null_rad and drs_ob.null_lim[0] eq drs.null_lim[0] and $
            drs_ob.null_lim[1] eq drs.null_lim[1] and drs_ob.null_range[0] eq drs.null_range[0] and drs_ob.null_range[1] eq drs.null_range[1] and drs.n_frob eq drs.n_frob and $ ; end of drs parameters
            data_ob.objname eq objname and data_ob.lam_cen eq lam_cen and data_ob.bck_mode eq bck_mode and data_ob.bfl_mode eq bfl_mode and data_ob.fit_mode eq fit_mode and data_ob.fra_mode eq fra_mode and data_ob.img_mode eq img_mode and $
            data_ob.ob_mode eq ob_mode and data_ob.pixscale eq pixscale and data_ob.precision eq precisio and data_ob.bck_irad eq bck_irad and data_ob.bck_orad eq bck_orad then begin ; end of paramaters on L1 images
            data[i_f] = data_ob
            goto, restored_ob
          endif
        endif
      endfor
    endif

    ; Prepare plot and log paths (define plot path anyway because used later for logs)
    plot_path_ob = plot_path + 'aper-' + string(apr_rad, format = '(I0)') + 'pix' + pth.sep + 'ob' + string(ob_id, format = '(I0)') + pth.sep
    plot_path_nsc = plot_path_ob + 'nsc' + pth.sep
    if not file_test(plot_path_ob) then file_mkdir, plot_path_ob
    spawn, 'chmod 775 ' + plot_path_ob
    if not file_test(plot_path_nsc) then file_mkdir, plot_path_nsc
    spawn, 'chmod 775 ' + plot_path_nsc

    ; Plot results
    plotnull = plot_path_ob + drs.date_obs + '_OB' + string(ob_id, format = '(I03)') + '_' + objname + '_dit-' + string(1000. * exptime, format = '(I0)') + 'ms_wav-' + string(1d+6 * lam_cen, format = '(I0)') + 'um'
    plotdiag = diag_path + drs.date_obs + '_OB' + string(ob_id, format = '(I03)') + '_' + objname + '_dit-' + string(1000. * exptime, format = '(I0)') + 'ms_wav-' + string(1d+6 * lam_cen, format = '(I0)') + 'um'

    ; 2. Read and process photometry files
    ; ************************************

    data_phot1 = LBTI_READL1DATA(data_path, ob_id, 'PHOT1', aper = apr_rad, idx_aper = i_aper1, header = hdr_phot1, /filter)
    data_phot2 = LBTI_READL1DATA(data_path, ob_id, 'PHOT2', aper = apr_rad, idx_aper = i_aper2, header = hdr_phot2, /filter)
    if (size(data_phot1))[2] eq 8 and (size(data_phot2))[2] eq 8 then begin
      phot_tot_phot1 = data_phot1.flx_tot[i_aper1]
      phot_err_phot1 = data_phot1.flx_err[i_aper1]
      n_phot1 = n_elements(data_phot1.mjd_obs)
      phot_tot_phot2 = data_phot2.flx_tot[i_aper2]
      phot_err_phot2 = data_phot2.flx_err[i_aper2]
      n_phot2 = n_elements(data_phot2.mjd_obs)
      if n_phot1 eq 1 then begin
        phot1 = phot_tot_phot1
        phot1_err = phot_err_phot1
        phot1_err_mean = phot_err_phot1
      endif else AVGSDV, phot_tot_phot1, phot1, phot1_err, phot1_err_mean, kappa = 3
      if n_phot2 eq 1 then begin
        phot2 = phot_tot_phot2
        phot2_err = phot_err_phot2
        phot2_err_mean = phot_err_phot2
      endif else AVGSDV, phot_tot_phot2, phot2, phot2_err, phot2_err_mean, kappa = 3
      ; Compute rms tipt/tilt
      tt1 = sqrt((data_phot1.xcen) ^ 2 + (data_phot1.ycen) ^ 2) * pixscale * 1d3
      AVGSDV, tt1, avg, rms_tt1
      tt2 = sqrt((data_phot2.xcen) ^ 2 + (data_phot2.ycen) ^ 2) * pixscale * 1d3
      AVGSDV, tt2, avg, rms_tt2
      ; Adjust error bar if negative background mode
      ; IF bck_mode GE 0 THEN BEGIN
      ; phot1_err = SQRT(phot1_err^2+phot1_err^2/      1+FXPAR(header1, 'N_FRBCK', /NOCONTINUE)/N_ELEMENTS(phot_tot_phot1))
      ; phot2_err = SQRT(1+FXPAR(header2, 'N_FRBCK', /NOCONTINUE)/N_ELEMENTS(phot_tot_phot2))
      ; ENDIF
      if n_phot1 gt 1 and n_phot2 gt 1 then begin
        time1 = reform(transpose(data_phot1.mjd_obs - min(data_phot1.mjd_obs)) * 24. * 60. * 60.) ; convert MJD to ellapsed seconds
        time2 = reform(transpose(data_phot2.mjd_obs - min(data_phot2.mjd_obs)) * 24. * 60. * 60.) ; convert MJD to ellapsed seconds
        PLOTALL, time1, reform(phot_tot_phot1), reform(phot_err_phot1), name = plotnull, tag = 'PHOT1', xtitle = 'Elapsed time [s]', ytitle = 'Photometric measurements [ADU]', title = ' ', /bold, /eps
        PLOTALL, time2, reform(phot_tot_phot2), reform(phot_err_phot2), name = plotnull, tag = 'PHOT2', xtitle = 'Elapsed time [s]', ytitle = 'Photometric measurements [ADU]', title = ' ', /bold, /eps
      endif
    endif else begin
      ; Skip OB if NSC reduction. Use max interferometric flux otherwise.
      if drs.null_mode eq 2 then begin
        if info gt 0 then print, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'No photometric frames'
        if lun gt 0 then printf, lun, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'No photometric frames'
        goto, skip_ob
      endif else print, 'WARNING: no photometric frames for OB ' + string(i_f, format = '(I0)') + ' (using maximum inteferometric flux as photometry)'
      phot1 = 0.25 * max(data_null.flx_tot[i_aper], i_max)
      phot2 = phot1
      phot1_err = data_null[i_max].flx_err[i_aper]
      phot2_err = phot1_err
      phot1_err_mean = phot1_err
      phot2_err_mean = phot2_err
      rms_tt1 = 0
      rms_tt2 = 0
    endelse
    int_err = 2.0 * abs(phot1 - phot2) / (phot1 + phot2)
    phot_tot = phot1 + phot2 + 2 * sqrt(abs(phot1 * phot2))
    phot_err = 2.0 * sqrt(phot1_err_mean ^ 2 + phot2_err_mean ^ 2) ; first order estimate dN = N*phot_err/phot_tot

    ; If using the best x% null technique, restart the reduction using the faintest CAL to correct for the background bias
    if drs.null_mode eq 1 then begin
      if n_elements(phot_faintcal) eq 0 then phot_faintcal = phot_tot
      if phot_tot lt phot_faintcal then begin
        phot_faintcal = phot_tot
        if info gt 0 then print, 'WARNING: restarting null computation using the fainter calibrator for background bias correction'
        if lun gt 0 then printf, lun, 'WARNING: restarting null computation using the fainter calibrator for background bias correction'
        goto, restart
      endif
    endif

    ; Save filtered L1 file
    if not keyword_set(no_save) then LBTI_SAVEL1FLX_1APER, data_phot1, hdr_phot1, i_aper1, ob_id = ob_id, tag = 'PHOT1'
    if not keyword_set(no_save) then LBTI_SAVEL1FLX_1APER, data_phot2, hdr_phot2, i_aper2, ob_id = ob_id, tag = 'PHOT2'

    ; 3. Read and process background file
    ; ***********************************

    ; Process background files only if the background mode is negative (which means that a constant frame has been subtraced from each frame of the current nod and
    ; background files are obtained at the same position as the the beam). If the background mode is greater than 0, use background mesured at a different position
    ; on the detector (so, stored directly in the NULL files).
    data_bckg = LBTI_READL1DATA(data_path, ob_id, 'BCKG', aper = apr_rad, idx_aper = id_aper, header = hdr_bckg, /filter)
    if (size(data_bckg))[2] eq 8 and bck_mode lt 0 then begin
      bck_tot_phot = data_bckg.flx_tot[id_aper] ; Photometric aperture flux in background nod
      bck_err_phot = data_bckg.flx_err[id_aper] ; Corresponding error
      bck_bck_flr = data_bckg.bck_tot[id_aper] ; background floor measured in complementary nod
      bck_tot_phot2 = data_bckg.flx_tot2[id_aper] ; Photometric aperture flux in background nod
      bck_bck_flr2 = data_bckg.bck_tot2[id_aper] ; background floor measured in complementary nod
      ; bck_bck_eflr  = data_bckg.bck_err[id_aper]  ; error on background floor
      nod_bckg = data_bckg.nod_id ; NOD ID of each background flux
      ; Plot data of requested
      if n_elements(bck_tot_phot) gt 1 then begin
        time = reform(transpose(data_bckg.mjd_obs - min(data_bckg.mjd_obs)) * 24. * 60. * 60.) ; convert MJD to ellapsed seconds
        PLOTALL, time, reform(bck_tot_phot), reform(bck_err_phot), name = plotnull, tag = 'BCKG', xtitle = 'Elapsed time [s]', ytitle = 'Background measurements [ADU]', title = '', /bold, /eps
        PLOTALL, time, reform(bck_bck_flr), 0., name = plotnull, tag = 'FLOOR-BCKG', xtitle = 'Elapsed time [s]', ytitle = 'Background floor [ADU]', title = '', /bold, /eps
        PLOTALL, reform(bck_bck_flr), reform(bck_tot_phot), 0., name = plotnull, tag = 'BIAS-vs-FLOOR', xtitle = 'Background floor [ADU]', ytitle = 'Background bias [ADU]', title = '', /bold, /eps, /scatter, /no_fft
        if max(abs(bck_bck_flr2)) ne 0 then PLOTALL, reform(bck_bck_flr2), reform(bck_tot_phot2), 0., name = plotnull, tag = 'BIAS2-vs-FLOOR2', xtitle = 'Background floor [ADU]', ytitle = 'Background bias [ADU]', title = '', /bold, /eps, /scatter, /no_fft
      endif
      ; REGION 1 (with the star).
      ; Compute uncorrected background BIAS
      if fra_mode eq 1 then AVGSDV, bck_tot_phot, bck_bias_unc, junk1, junk2, kappa = 5 else bck_bias_unc = 0. ; If frame mode of 0, then the bias is corrected by image subtraction and we have no info here...
    endif else begin
      coeff = [0, 0]
      nod_bckg = 0
      bck_avg = 0
      bck_rms = 0
      bck_rms_mean = 0 ; this is the external error, which is included in the NSC error if bck_mode is positive
      bck_bias_unc = 0
      ; Skip OB if no BCKG file
      if bck_mode lt 0 then begin
        if info gt 0 then print, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'No background frames found (skipped)'
        if lun gt 0 then printf, lun, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'No background frames found (skipped)'
        goto, skip_ob
      endif
    endelse

    ; Remove mean value from background for NSC (if from different NODs, do it per NOD!)
    nod_bckg_uniq = nod_bckg[uniq(nod_bckg, sort(nod_bckg))]
    n_nod_bckg = n_elements(nod_bckg_uniq)
    for ib = 0, n_nod_bckg - 1 do begin
      idx_nod_bckg = where(nod_bckg eq nod_bckg_uniq[ib])
      AVGSDV, bck_tot_phot[idx_nod_bckg], bck_avg, tmp, tmp2, kappa = 5
      bck_tot_phot[idx_nod_bckg] = bck_tot_phot[idx_nod_bckg] - bck_avg ; remove mean value for NSC fitting
      data_bckg.flx_tot[id_aper] = bck_tot_phot
    endfor

    ; Plot corrected background
    if n_elements(bck_tot_phot) gt 1 then begin
      time = reform(transpose(data_bckg.mjd_obs - min(data_bckg.mjd_obs)) * 24. * 60. * 60.) ; convert MJD to ellapsed seconds
      PLOTALL, time, reform(bck_tot_phot), reform(bck_err_phot), name = plotnull, tag = 'BCKG-COR', xtitle = 'Elapsed time [s]', ytitle = 'Background measurements [ADU]', title = '', /bold, /eps
    endif

    ; Save filtered L1 background file
    if not keyword_set(no_save) then LBTI_SAVEL1FLX_1APER, data_bckg, hdr_bckg, id_aper, ob_id = ob_id, tag = 'BCKG'

    ; 4. Process and filter nulls
    ; ***************************

    ; Plot diagnostic infos about the PWV loop (+log file for Phil, only of H-band data)
    ; Plot null sequence before any filtering
    time = reform(transpose(data_null.mjd_obs - min(data_null.mjd_obs)) * 24. * 60. * 60.) ; convert MJD to ellapsed seconds
    null_tot_phot = reform(data_null.flx_tot[i_aper]) - bck_bias_unc ; correct for background BIAS
    PLOTALL, time, 100. * reform(data_null.flx_tot[i_aper]) / phot_tot, 100. * reform(data_null.flx_err[i_aper]) / phot_tot, name = plotnull, tag = 'NULL', xtitle = 'Elapsed time [s]', ytitle = 'Instantaneous null depth [%]', title = '', /ylog, yrange = [0.1, 100.0], /bold, /eps, /kernel
    ; Plot PHASECam diagnostic INFO
    if max(data_null.pcphmean2) ne 0 then begin
      ; Extract data
      pwv_pha = data_null.pcphmean2 - data_null.pcphmean
      PLOTALL, pwv_pha, 100. * null_tot_phot / phot_tot, 0., name = plotdiag, tag = 'PWV-NULL', xtitle = 'H-K phase [deg]', ytitle = 'Instantaneous null depth [%]', title = '', /bold, /eps, /no_fft, /scatter
      PLOTALL, time, pwv_pha, 0, name = plotdiag, tag = 'PWV', xtitle = 'Elapsed time [s]', ytitle = 'H-K phase [deg]', title = '', /bold, /eps
      ; Temporary hack
      PREP_PS, /bold
      loadct, 0, /silent
      device, filename = plotdiag + '_PWV-vs-NULL.eps', /encaps, /color, xsize = 17.8, ysize = 14.7, /times
      nulll = 100. * null_tot_phot / phot_tot
      xrange = [min(time), max(time)]
      yrange = [min(nulll), max(nulll)]
      PLOT, [0], [0], xtitle = 'Elapsed time [s]', ytitle = 'Null depth [%]', title = title_day, xstyle = 1, ystyle = 9, $
        xrange = xrange, yrange = yrange, xthick = xthick, ythick = ythick, charthick = charthick, charsize = charsize
      AXIS, ytitle = 'H-K phase [deg]', /yaxis, ystyle = 1, yrange = [min(pwv_pha), max(pwv_pha)]
      loadct, 13, /silent
      oplot, time, nulll, color = 100
      oplot, time, (pwv_pha - min(pwv_pha)) * (max(nulll) - min(nulll)) / (max(pwv_pha) - min(pwv_pha)), color = 250
      device, /close
      END_PS

      ; Print file for Phil
      log_file2 = diag_path + date + '_OB' + string(i_f, format = '(I0)') + '.txt'
      openw, lun2, log_file2, /get_lun, width = 800
      printf, lun2, 't0=' + string(data_bckg[0].mjd_obs, format = '(F16.8)')
      printf, lun2, 'time[s];null[%];PHAB1[deg];PHAB2[deg]'
      for i = 0, n_elements(nulll) - 1 do printf, lun2, time[i], ';', nulll[i], ';', data_null[i].pcphmean, ';', data_null[i].pcphmean2
      close, lun2
      free_lun, lun2
    endif

    ; Now filter the null data.
    ; ------------------------

    ; First, remove fringe jump. 18% is considerd as a jump, keep lowest sequence + those with a mean within 5% of the lowest one(and enough frames!)
    ; Only for data after July 2014 (otherwise, CG only and this does not apply)
    if mean(data_null.mjd_obs) gt 56839 then begin
      tmp_null = data_null.flx_tot[i_aper] / phot_tot
      tmp_null = REMOVE_NULLJUMP(tmp_null, 0.18, 0.05, 20, idx_out = idx_null)
      n_null = n_elements(idx_null)
      if n_null lt drs.min_fr then begin
        if info gt 0 then print, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived fringe jump filtering (skipped)'
        if lun gt 0 then printf, lun, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived fringe jump filtering  (skipped)'
        goto, skip_ob
      endif
    endif else idx_null = lindgen(n_elements(data_null.mjd_obs))

    ; Compute mean fringe SNR and keep only the best one (low SNR usualy means bad phase setpoint)
    MEANCLIP, data_null[idx_null].pcmsnr, avg_snr, rms_snr, clipsig = sig_out
    idx_null = idx_null[where(data_null[idx_null].pcmsnr ge (avg_snr - sig_out * rms_snr), n_null)]
    if n_null lt drs.min_fr then begin
      if info gt 0 then print, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived the SNR cut (skipped)'
      if lun gt 0 then printf, lun, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived the SNR cut  (skipped)'
      goto, skip_ob
    endif
    ; Compute mean phase noise and reject high noise
    MEANCLIP, data_null[idx_null].pcphstd, avg_phstd, rms_phstd, clipsig = sig_out
    idx_null = idx_null[where(abs(data_null[idx_null].pcphstd - avg_phstd) le sig_out * rms_phstd, n_null)]
    if n_null lt drs.min_fr then begin
      if info gt 0 then print, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived the phase noise cut (skipped)'
      if lun gt 0 then printf, lun, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived the phase noise cut  (skipped)'
      goto, skip_ob
    endif
    ; Compute mean phase and reject frame to far from setpoint
    MEANCLIP, data_null[idx_null].pcphmean, avg_phmean, rms_phmean, clipsig = sig_out
    idx_null = idx_null[where(abs(data_null[idx_null].pcphmean - avg_phmean) le sig_out * rms_phmean, n_null)]
    if n_null lt drs.min_fr then begin
      if info gt 0 then print, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived the mean phase cut (skipped)'
      if lun gt 0 then printf, lun, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived the mean phase cut  (skipped)'
      goto, skip_ob
    endif
    ; Only keep the frame in null_lim range
    if n_elements(drs.null_lim) eq 2 then begin
      if drs.null_lim[0] ne 0 or drs.null_lim[1] ne 0 then begin
        flx_tmp = data_null[idx_null].flx_tot[i_aper] - bck_bias_unc ; correct for background bias
        idx_null = idx_null[where(flx_tmp / phot_tot ge drs.null_lim[0] and flx_tmp / phot_tot le drs.null_lim[1], n_null, /null)]
        if n_null lt drs.min_fr then begin
          if info gt 0 then print, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') in the NULL_LIM range (skipped)'
          if lun gt 0 then printf, lun, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') in the NULL_LIM range (skipped)'
          goto, skip_ob
        endif
      endif
    endif
    ; Remove frames more than 5 sigma away from the median
    flx_tmp = data_null[idx_null].flx_tot[i_aper]
    AVGSDV, flx_tmp, avg_tmp, rms_tmp, rms_tmp_m, kappa = sig_out
    idx_null = idx_null[where(abs(flx_tmp - median(flx_tmp)) lt sig_out * rms_tmp, n_null, /null)]
    if n_null lt drs.min_fr then begin
      if info gt 0 then print, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived after the median cut (skipped)'
      if lun gt 0 then printf, lun, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived after the median cut (skipped)'
      goto, skip_ob
    endif
    ; Use background in surrounding regions to remove outliers
    flx_tmp = data_null[idx_null].bck_tot[i_aper]
    AVGSDV, flx_tmp, avg_tmp, rms_tmp, rms_tmp_m, kappa = sig_out
    idx_null = idx_null[where(abs(flx_tmp - avg_tmp) le sig_out * rms_tmp, n_null, /null)]
    if n_null lt drs.min_fr then begin
      if info gt 0 then print, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived after the background check (skipped)'
      if lun gt 0 then printf, lun, string(ob_id, format = '(I03)'), '  ', string(objname, format = '(A8)'), '  ', string(flag, format = '(A3)'), '  ', 'Not enough null frames (' + string(n_null, format = '(I0)') + ') survived after the background check (skipped)'
      goto, skip_ob
    endif
    ; Extract data
    data_null = data_null[idx_null]
    ; Extract data for chosen aperture radius
    null_tot_phot = reform(data_null.flx_tot[i_aper])
    null_err_phot = reform(data_null.flx_err[i_aper]) ; flux measurements at null and corresponding errors
    null_tot_bckg = reform(data_null.bck_tot[i_aper])
    null_err_bckg = reform(data_null.bck_err[i_aper]) ; background measurements (per pixel) and corresponding errors (in the background region surrounding the photometric aperture)
    null_tot_phot2 = reform(data_null.flx_tot2[i_aper])
    null_err_phot2 = reform(data_null.flx_err2[i_aper]) ; flux measurements at null and corresponding errors
    null_tot_bckg2 = reform(data_null.bck_tot2[i_aper])
    null_err_bckg2 = reform(data_null.bck_err2[i_aper]) ; background measurements (per pixel) and corresponding errors (in the background region surrounding the photometric aperture)

    ; Save filtered L1 file
    if not keyword_set(no_save) then LBTI_SAVEL1FLX_1APER, data_null, header, i_aper, ob_id = ob_id, tag = 'NULL'

    ; Debias null and off-axis measurements by using measurements in the complementary NOD (only if FRA_MODE of 1)
    ; ----------------------------------------------------------------------------------

    ; First, determine degree of POLYNOMIAL fit by comparing the range of background floors in the complentary NOD(s) and the science nod
    ; The rule of thumb here is to use a constant if the difference between the mean background floors in each NOD is less than twice the standard deviation of background floors in teh complementary NOD(s)
    ; Use a degree of 2 otherwise
    if fra_mode eq 1 then begin
      AVGSDV, bck_bck_flr, bck_bck_avg, bck_bck_rms, bck_bck_mrms, kappa = 5 ; compute mean and RMS of background floors in complementary nods
      if (mean(null_tot_bckg) - bck_bck_avg) le 2 * bck_bck_rms then poly_deg = 1 else poly_deg = 0 ; compute polynomial degree
      coeff = poly_fit(bck_bck_flr, bck_tot_phot, poly_deg) ; fit polynomial based on complementary nods
      if max(abs(bck_tot_phot2)) gt 0 then coeff2 = poly_fit(bck_bck_flr2, bck_tot_phot2, poly_deg) else coeff2 = replicate(0, poly_deg + 1)
      ; Now debias background measurements in complementary NODS (both regions of the detector)
      for i = 0, poly_deg do begin
        bck_tot_phot -= coeff[i] * bck_bck_flr ^ i ; FLux in complementary nod (in BCKG file)
        null_tot_phot -= coeff[i] * null_tot_bckg ^ i ; Flux on source
        null_tot_phot2 -= coeff2[i] * null_tot_bckg2 ^ i ; Flux in nearby empty region (in NULL file)
      endfor
    endif

    ; Now compute NSC-related values and PLOT results
    ; -----------------------------------------------

    ; Now compute RMS over full sequence
    AVGSDV, bck_tot_phot, tmp, bck_rms, bck_rms_mean ; , KAPPA=5
    ; Derive minimum possible astro NULL (defined arbitrary here as -1 sigma the background variation relative to peak)
    nas_min = -1. * bck_rms / phot_tot
    ; Compute raw null RMS and AVG (avter debias)
    AVGSDV, null_tot_phot / phot_tot, null_avg, null_rms
    ; Degrade null (used to simulate fainter stars)
    ; IF STRTRIM(FXPAR(header, 'FLAG')) EQ 'SCI' THEN BEGIN
    ; null_tot_phot += bck_tot_phot*SQRT(7./2.-1)
    ; bck_tot_phot *= SQRT(7./2.)
    ; ENDIF
    ; Compute mean background bias
    if max(null_tot_phot2) ne 0 then begin
      idx_ok = where(null_err_phot2 ne 0, n_ok)
      null_tot_bckg2 = null_tot_bckg2[idx_ok]
      null_tot_phot2 = null_tot_phot2[idx_ok]
      null_err_phot2 = null_err_phot2[idx_ok]
      AVGSDV, null_tot_phot2, bck_bias_cor, rms_bias, rms_bias_mean, kappa = sig_out, weight = (1 / null_err_phot2) ^ 2
      null_bckg = [bck_bias_unc, bck_bias_cor, rms_bias, rms_bias_mean, 0] / phot_tot ; Uncorrected background bias (same position), corrected background bias (other position), RMS at other position, and RMS on the mean
    endif else begin
      null_bckg = [bck_bias_unc, 0., 0., 0] / phot_tot
    endelse
    ; Compute mean background level (per pixel)
    AVGSDV, null_tot_bckg, avg_bck, rms_bck, kappa = sig_out, weight = (1. / null_err_bckg) ^ 2
    bck_rms_est = sqrt(napr_pix) * sqrt(nsky_pix) * rms_bck ; estimated RMS of total background in photometric aperture (not used after all)
    ; Plot data of requested
    time = reform(transpose(data_null.mjd_obs - min(data_null.mjd_obs)) * 24. * 60. * 60.) ; convert MJD to ellapsed seconds
    PLOTALL, time, 100. * null_tot_phot / phot_tot, 100. * null_err_phot / phot_tot, name = plotnull, tag = 'NULLOK', xtitle = 'Elapsed time [s]', ytitle = 'Instantaneous null depth [%]', title = '', /ylog, yrange = [0.1, 100.0], /bold, /eps, /kernel
    PLOTALL, time, null_tot_bckg, null_err_bckg, name = plotnull, tag = 'FLOOR-NULL', xtitle = 'Elapsed time [s]', ytitle = 'Background floor [ADU]', title = '', /bold, /eps
    PLOTALL, null_tot_bckg, null_tot_phot, 0, name = plotnull, tag = 'NULL-VS-FLOOR', xtitle = 'Background floor [ADU]', ytitle = 'Background bias [ADU]', title = '', /bold, /eps, /scatter, /no_fft
    ; PLOTALL, time, fake_bckg, 0, NAME=plotnull, TAG='FAKE-BCKG', XTITLE='Elapsed time [s]', YTITLE='Estimated backgroud in photometric aperture [ADU]', TITLE='', /BOLD,  /EPS
    ; PLOTALL, data_null.pcphmean, 100.*null_tot_phot/phot_tot, 0., NAME=plotnull, TAG='PCPHMEAN2', XTITLE='Phase [deg]', YTITLE='Instantaneous null depth [%]', TITLE='', /BOLD, /EPS, /NO_FFT, /SCATTER
    ; PLOTALL, data_null.pcphstd, 100.*null_tot_phot/phot_tot, 0., NAME=plotnull, TAG='PCPHSTD', XTITLE='Phase jitter [deg]', YTITLE='Instantaneous null depth [%]', TITLE='', /BOLD, /EPS, /NO_FFT, /SCATTER
    ; PLOTALL, data_null.fpcazm, 100.*null_tot_phot/phot_tot, 0., NAME=plotnull, TAG='FPCAZM', XTITLE='Azimuth [mas]', YTITLE='Instantaneous null depth [%]', TITLE='', /BOLD, /EPS, /NO_FFT, /SCATTER
    ; PLOTALL, data_null.fpcelm, 100.*null_tot_phot/phot_tot, 0., NAME=plotnull, TAG='FPCELM', XTITLE='Elevation [mas]', YTITLE='Instantaneous null depth [%]', TITLE='', /BOLD, /EPS, /NO_FFT, /SCATTER
    if max(data_null.pcphstd) ne 0 then PLOTALL, time, data_null.pcphstd, 0., name = plotnull, tag = 'PCPHSTD', xtitle = 'Elapsed time [s]', ytitle = 'Phase jitter wittin DIT [deg]', title = '', /bold, /eps ; , /NO_FFT
    if max(data_null.pcphmean) ne 0 then PLOTALL, time, data_null.pcphmean, 0., name = plotnull, tag = 'PCPHMEAN', xtitle = 'Elapsed time [s]', ytitle = 'Mean phase per DIT [deg]', title = '', /bold, /eps ; , /NO_FFT
    if max(null_tot_phot2) ne 0 then PLOTALL, time[idx_ok], 100. * null_tot_phot2 / phot_tot, 100. * null_err_phot2 / phot_tot, name = plotnull, tag = 'BIAS', xtitle = 'Elapsed time [s]', ytitle = 'Background null [%]', title = '', /bold, /eps
    if max(null_tot_phot2) ne 0 then PLOTALL, null_tot_bckg2, null_tot_phot2, 0., name = plotnull, tag = 'NULL2-VS-FLOOR2', xtitle = 'Background floor [ADU]', ytitle = 'Background bias [ADU]', title = '', /bold, /eps, /scatter, /no_fft

    ; 4. Compute PHASECam telemetry AVG and RMS per OBs
    ; *************************************************

    ; Average phase over the last (x) ms
    AVGSDV, data_null.pcphmean, pcphmean_avg, pcphmean_rms, pcphmean_err, kappa = 5
    ; IF MEAN(data_null.pcphmean-data_null.pcplsp) GT 5*pcphmean_rms THEN BEGIN
    ; MESSAGE, "Double-check this OB: measured phase doesn't match setpoint ( " + STRING(MEAN(data_null.pcphmean), FORMAT='(I0)') + " degrees instead of " + STRING(MEAN(data_null.pcplsp), FORMAT='(I0)') + ')', /CONTINUE
    ; IF lun GT 0 THEN PRINTF, lun, "Double-check this OB: measured phase doesn't match setpoint ( " + STRING(MEAN(data_null.pcphmean), FORMAT='(I0)') + " degrees instead of " + STRING(MEAN(data_null.pcplsp), FORMAT='(I0)') + ')'
    ; ENDIF

    ; RMS phase over the last (x) ms (remove outliers)
    AVGSDV, data_null.pcphstd, pcphstd_avg, pcphstd_rms, pcphstd_err, kappa = 5

    ; Scale to NOMIC (WARNING: hardcoded for now. These are the values used by Elwood's to compute PCPHMCOS and PCPHMSIN. We should use the same here!!)
    ; WARNING 2: pcphmean should include the dither pattern to compute correctly pcphmcos and pcphmsin below
    pcphmean = data_null.pcphmean * !dtor * 2.2 / 11.1

    ; Compute effective wavelength for PHASECam (read target spectrum and instrument throughput)
    pcwav = FXPAR(header, 'PCWAV', /nocontinue)
    if !err eq -1 then begin
      READ_TABLE, 'nodrs/input/phasecam_K.txt', lam_k, trans, first = 0, separator = ' '
      idx_k = where(lam_k gt 1.95 and lam_k lt 2.45, n_lam)
      lam_k = lam_k[idx_k]
      trans = trans[idx_k]
      GET_TGT, objname, tgt, database = pth.input_path + drs.database
      pick = READ_PICKLES(pth.pickles_path + tgt.spectrum, lam_min = 1.95d-6, lam_max = 2.45d-6, nlambda = n_lam, standard = 3)
      spec = interpol(pick[1, *], pick[0, *], lam_k * 1d-6)
      pcwav = 1. / (total((trans * spec) ^ 2 * 1. / lam_k) / total((trans * spec) ^ 2)) * 1d-6
    endif

    ; Subtract estimated null floor if requested (convert measured K-band phase in degrees to N-band phase in radians)
    ; For NSC, this information is used in the fit rather than subtracted from each null measurements
    phi_rms_dit = fltarr(n_null)
    phi_avg_dit = fltarr(n_null)
    pcphmcos = fltarr(n_null)
    pcphmsin = fltarr(n_null)
    if drs.null_mode eq 2 then begin
      if drs.null_cor eq 1 then null_tot_phot -= phot_tot * (reform(data_null.pcphstd) * !dtor * pcwav / lam_cen) ^ 2 / 2
      if drs.null_cor ge 2 then begin
        phi_rms_dit = reform(data_null.pcphstd) * !dtor * pcwav / lam_cen ; convert to radian and scale to NOMIC's wavelength
        if TAG_EXIST(data_null, 'pcphmcos') and date ne '230523' and date ne '230525' then pcphmcos = data_null.pcphmcos * cos(pcphmean) + data_null.pcphmsin * sin(pcphmean) ; mean cos(phase-<phase>) over DIT + backward compatibility. For 230523 and 230525, PCPHMCOS and PCPHMSIN are wrong.
        if TAG_EXIST(data_null, 'pcphmsin') and date ne '230523' and date ne '230525' then pcphmsin = data_null.pcphmsin * cos(pcphmean) - data_null.pcphmcos * sin(pcphmean) ; mean sin(phase-<phase>) over DIT + backward compatibility. For 230523 and 230525, PCPHMCOS and PCPHMSIN are wrong.
        if drs.null_cor eq 3 then message, 'NULL_COR of 3 is depreciated' ; phi_avg_dit = REFORM(data_null.pcphmean-pcphmean_avg)*!dtor*pcwav/lam_cen;)^3/3              ; Mean phase over DIT (not used anymnore)
      endif
    endif else if drs.null_cor gt 0 then null_tot_phot -= phot_tot * (reform(data_null.pcphstd) * !dtor * pcwav / lam_cen) ^ 2 / 2

    ; Compute phase RMS at NOMIC's frequency (used as minimum phase jitter for NSC)
    pcphmean = (data_null.pcphmean - data_null.spdthpos) * !dtor * 2.2 / 11.1
    AVGSDV, pcphmean, pcphmean_avg, pcphmean_rms, pcphmean_err, kappa = 5
    phi_rms_min = 0.5 * pcphmean_rms ; use 0.5 to account for possible PHASECam noise (to be checked). Not reliable right now because of fringe jumps!

    ; Read dither pattern and scale it to NOMIC's wavelength
    spdthpos = data_null.spdthpos * !dtor * pcwav / lam_cen ; phase in RADIAN at NOMIC's wavelength
    if max(abs(spdthpos)) gt 0 then PLOTALL, time, spdthpos, 0., name = plotnull, tag = 'SPDTHPOS', xtitle = 'Elapsed time [s]', ytitle = 'OPD dither pattern [rad]', title = '', /bold, /eps ; , /NO_FFT
    if max(abs(pcphmcos)) gt 0 then PLOTALL, time, pcphmcos, 0., name = plotnull, tag = 'PCPHMCOS', xtitle = 'Elapsed time [s]', ytitle = 'mean COS(PHASE) over DIT', title = '', /bold, /eps, /no_fft
    if max(abs(pcphmsin)) gt 0 then PLOTALL, time, pcphmsin, 0., name = plotnull, tag = 'PCPHMSIN', xtitle = 'Elapsed time [s]', ytitle = 'mean SIN(PHASE) over DIT', title = '', /bold, /eps, /no_fft

    ; 5. Prepare sequence for NSC fit (not used anymore)
    ; *******************************

    ; Compute the true photometric variations
    if drs.null_mode eq 2 then begin
      ; Number of points to create the synthetic sequences
      n_gen = 1d5
      ; Compute residual RMS (after subtraction of background contribution)
      if bck_mode ge 0 then bck_rms_tmp = rms_bias else bck_rms_tmp = bck_rms
      if phot1_err gt bck_rms_tmp then phot1_nsc_rms = sqrt(phot1_err ^ 2 - bck_rms_tmp ^ 2) else phot1_nsc_rms = 0
      if phot2_err gt bck_rms_tmp then phot2_nsc_rms = sqrt(phot2_err ^ 2 - bck_rms_tmp ^ 2) else phot2_nsc_rms = 0
      ; Create synthetic photometric sequences (Gaussian assumption)
      ; WARNING. RANDOMN updates the seed!
      phot_tot_phot1 = phot1 + randomn(noseed1, n_gen) * phot1_nsc_rms
      phot_tot_phot2 = phot2 + randomn(noseed2, n_gen) * phot2_nsc_rms
      ; bck_tot_phot   = RANDOMN(noseed, n_gen)*bck_rms_tmp
    endif

    ; 6. Compute null per OB
    ; **********************
    ratio = 1
    mode = 0
    ; Statistical reduction if null mode is 2
    if drs.null_mode ne 2 then begin
      ; Classical reduction
      case drs.null_mode of
        0: mode = 1
        1: begin
          ; Add noise to account for different brightness (manual at this point)
          ratio = drs.keep_ratio
          x_cor = sqrt((phot_tot / phot_faintcal) ^ 2 - 1)
          if bck_mode ge 0 then null_tot_phot += (null_tot_phot2 - bck_bias_cor) * x_cor else null_tot_phot += (null_tot_bckg - avg_bck) * x_cor ; Subtract mean value to only account for photon noise
        end
        else: print, 'Unknown null mode'
      endcase
      ; Compute null (for star) and parse output strucuture similarly to NSC
      null_res = BIN_NULLDATA(null_tot_phot, null_err_phot, ratio = ratio, mode = mode) / phot_tot
      data[i_f].null_star.null_avg.nas = double(null_res[0])
      data[i_f].null_star.null_avg.nas_err_low = double(null_res[1])
      data[i_f].null_star.null_avg.nas_err_sup = double(null_res[1])
    endif else begin
      ; Statistical reduction. Depending on the background mode, we use either the background measured in a nearby empty region of the detector (bck_mode > 0)
      ; or use the same position but in the complemntary nod (bck_mode < 0).
      if bck_mode ge 0 then dark = null_tot_phot2 - bck_bias_cor else dark = bck_tot_phot
      arraydata = {dark2: dark, a: phot_tot_phot1, b: phot_tot_phot2, ab: null_tot_phot, nas_min: nas_min, phi_rms_min: phi_rms_min, phi_avg_t: phi_avg_dit, phi_var_t: phi_rms_dit ^ 2, pcphmcos_t: pcphmcos, pcphmsin_t: pcphmsin, spdthpos: spdthpos}
      data[i_f].null_star = CALL_NSC(arraydata, cube_size = drs.nsc_cube, n_bootstrap = drs.n_btstrp, nbins_factor = drs.nsc_bfac, no_multi = no_multi, null_range = drs.null_range, file_id = ob_id, wav_eff = lam_cen, bandwidth = bdwdth, $
        occ_min = drs.nsc_omin, variable = drs.nsc_bins, info = info, file_path = plot_path_nsc, /eps, /plot)
    endelse

    ; Save data in output structure
    ; *****************************

    ; Header information
    data[i_f].instrum = strtrim(FXPAR(header, 'INSTRUME'))
    data[i_f].objname = objname
    data[i_f].pid = FXPAR(header, 'PID')
    data[i_f].flag = strtrim(FXPAR(header, 'FLAG'))
    data[i_f].mjd_obs = mean(data_null.mjd_obs, /double)
    data[i_f].lbt_utc = strtrim(FXPAR(header, 'LBT_UTC')) ; at start of sequence
    data[i_f].lbt_lst = strtrim(FXPAR(header, 'LBT_LST')) ; at start of sequence
    data[i_f].lbt_ra = strtrim(FXPAR(header, 'LBT_RA')) ; at start of sequence
    data[i_f].lbt_dec = strtrim(FXPAR(header, 'LBT_DEC')) ; at start of sequence
    data[i_f].lbt_para = FXPAR(header, 'LBT_PARA') ; mean of sequence
    data[i_f].lbt_alt = FXPAR(header, 'LBT_ALT') ; mean of sequence
    data[i_f].lbt_az = FXPAR(header, 'LBT_AZ') ; mean of sequence
    ; data[i_f].lbt_ha   = 0.                                      ; not yet in LBT header (temporary implementation)

    ; Detector parameters
    data[i_f].int_time = exptime
    data[i_f].acq_time = FXPAR(header, 'ACQTIME')
    data[i_f].lam_cen = lam_cen
    data[i_f].bandwidth = bdwdth
    data[i_f].n_xpix = FXPAR(header, 'N_XPIX')
    data[i_f].n_ypix = FXPAR(header, 'N_YPIX')
    data[i_f].pixscale = pixscale
    data[i_f].xcen = mean(data_null.xcen)
    data[i_f].ycen = mean(data_null.ycen)
    data[i_f].smplmode = FXPAR(header, 'SMPLMODE')
    data[i_f].detmode = FXPAR(header, 'DETMODE')
    data[i_f].pagain = FXPAR(header, 'PAGAIN')
    data[i_f].pabandw = FXPAR(header, 'PABANDW')
    data[i_f].detbias = FXPAR(header, 'DETBIAS')
    data[i_f].eperadu = FXPAR(header, 'EPERADU')
    data[i_f].n_coadd = FXPAR(header, 'N_COADD')

    ; Filter keywords
    data[i_f].ubc_dxsp = FXPAR(header, 'UBC_DXSP')
    data[i_f].ubc_sxsp = FXPAR(header, 'UBC_SXSP')
    data[i_f].nic_fstp = FXPAR(header, 'NIC_FSTP')
    data[i_f].nic_beam = FXPAR(header, 'NIC_BEAM')
    data[i_f].lmir_fw1 = FXPAR(header, 'LMIR_FW1')
    data[i_f].lmir_fw2 = FXPAR(header, 'LMIR_FW2')
    data[i_f].lmir_fw3 = FXPAR(header, 'LMIR_FW3')
    data[i_f].lmir_fw4 = FXPAR(header, 'LMIR_FW4')
    data[i_f].nom_fw1 = FXPAR(header, 'NOM_FW1')
    data[i_f].nom_fw2 = FXPAR(header, 'NOM_FW2')
    data[i_f].nom_apw = FXPAR(header, 'NOM_APW')
    data[i_f].pha_fw1 = FXPAR(header, 'PHA_FW1')
    data[i_f].pha_fw2 = FXPAR(header, 'PHA_FW2')
    data[i_f].pha_img = FXPAR(header, 'PHA_IMG')
    data[i_f].nil_dic = FXPAR(header, 'NIL_DIC')

    ; AO keywords
    data[i_f].daomode = FXPAR(header, 'DAOMODE')
    data[i_f].daostrhl = FXPAR(header, 'DAOSTRHL')
    data[i_f].dcmodes = FXPAR(header, 'DCMODES')
    data[i_f].dloopon = FXPAR(header, 'DLOOPON')
    data[i_f].dloopgn = FXPAR(header, 'DLOOPGN')
    data[i_f].dwfscfrq = FXPAR(header, 'DWFSCFRQ')
    data[i_f].dwfscbin = FXPAR(header, 'DWFSCBIN')
    data[i_f].saomode = FXPAR(header, 'SAOMODE')
    data[i_f].saostrhl = FXPAR(header, 'SAOSTRHL')
    data[i_f].scmodes = FXPAR(header, 'SCMODES')
    data[i_f].sloopon = FXPAR(header, 'SLOOPON')
    data[i_f].sloopgn = FXPAR(header, 'SLOOPGN')
    data[i_f].swfscfrq = FXPAR(header, 'SWFSCFRQ')
    data[i_f].swfscbin = FXPAR(header, 'SWFSCBIN')

    ; Phasecam keywords
    data[i_f].pcclosed = FXPAR(header, 'PCCLOSED')
    data[i_f].plc_wav = pcwav
    data[i_f].plc_spec = pcspec
    data[i_f].pcplsp = mean(data_null.pcplsp)
    data[i_f].pctipsp = mean(data_null.pctipsp)
    data[i_f].pctltsp = mean(data_null.pctltsp)
    data[i_f].spc_pist = FXPAR(header, 'SPCPIST')
    data[i_f].spc_az = FXPAR(header, 'SPCAZ')
    data[i_f].spc_el = FXPAR(header, 'SPCEL')
    data[i_f].fpc_pistm = mean(data_null.fpcpistm)
    data[i_f].fpc_pists = mean(data_null.fpcpists)
    data[i_f].fpc_azm = mean(data_null.fpcazm)
    data[i_f].fpc_azs = mean(data_null.fpcazs)
    data[i_f].fpc_elm = mean(data_null.fpcelm)
    data[i_f].fpc_els = mean(data_null.fpcels)
    data[i_f].pcphmean = pcphmean_avg
    data[i_f].pcphmean_err = pcphmean_err
    data[i_f].pcphstd = pcphstd_avg
    data[i_f].pcphstd_err = pcphstd_err
    data[i_f].pcmsnr = mean(data_null.pcmsnr)

    ; Weather parameters
    data[i_f].seeing = FXPAR(header, 'SEEING')
    data[i_f].smttau = FXPAR(header, 'SMTTAU')
    data[i_f].lbttemp = FXPAR(header, 'LBTTEMP')
    data[i_f].winddir = FXPAR(header, 'WINDDIR')
    data[i_f].windspd = FXPAR(header, 'WINDSPD')

    ; Nulls
    data[i_f].null_avg = null_avg ; Diagnostic info
    data[i_f].null_rms = null_rms ; Diagnostic info

    ; Photometry
    data[i_f].photdx_avg = phot1
    data[i_f].photsx_avg = phot2
    data[i_f].photdx_rms = phot1_err
    data[i_f].photsx_rms = phot2_err
    data[i_f].photdx_snr = phot1 / phot1_err
    data[i_f].photsx_snr = phot2 / phot2_err
    data[i_f].phot_avg = phot_tot
    data[i_f].phot_err = phot_err
    data[i_f].int_err = int_err

    ; Background
    data[i_f].bckg_avg = bck_avg / phot_tot ; AVG background per photometric aperture
    data[i_f].bckg_rms = bck_rms / phot_tot ; RMS background per photometric aperture
    data[i_f].bckg_err = bck_rms_mean / phot_tot ; RMS background per photometric aperture
    data[i_f].bias_unc = null_bckg[0] ; Uncorrected background bias per photometric aperture
    data[i_f].bias_cor = null_bckg[1] ; Empty-region background bias per photometric aperture
    data[i_f].bias_err = null_bckg[2] ; Empty-region background bias error per photometric aperture
    data[i_f].bckg_meas = avg_bck ; AVG background per pixel (in the background region surrounding the photometric aperture)
    data[i_f].bckg_emeas = rms_bck ; RMS background per pixel (in the background region surrounding the photometric aperture)

    ; Reduction paramaters
    data[i_f].file_id = ob_id
    data[i_f].ob_id = ob_id ; currently OB id and file if are the same
    data[i_f].pt_id = pt_id
    data[i_f].nod_id = data_null[0].nod_id
    data[i_f].nfr_in = n_null0
    data[i_f].nfr_rej = n_null0 - n_null ;+n_rej  ; now n_rej includes the backgrounds frames obtained between nods...so no
    data[i_f].nfr_ob = floor(ratio * n_null)
    data[i_f].nod_frq = FXPAR(header, 'NOD_FRQ')
    data[i_f].n_frbck = FXPAR(header, 'N_FRBCK')
    data[i_f].rms_opd = sqrt(data[i_f].null_err) * (0.5 * data[i_f].lam_cen) ; first order approximation
    data[i_f].rms_tt = [rms_tt1, rms_tt2]
    data[i_f].aper_rad = apr_rad
    data[i_f].bck_irad = bck_irad
    data[i_f].bck_orad = bck_orad
    data[i_f].fit_mode = fit_mode
    data[i_f].flx_mode = flx_mode
    data[i_f].fra_mode = fra_mode
    data[i_f].img_mode = img_mode
    data[i_f].ob_mode = ob_mode
    data[i_f].bck_mode = bck_mode
    data[i_f].bfl_mode = bfl_mode
    data[i_f].precision = precisio

    ; Compute quality flag. Current defintion follwoing discussion with RMG is 1: no suitable for science, 2: poor quality due to seeing or PWV, 3: is good quality.
    ; Per my experience, I use the follwoing definition:
    ; - 1 : seeing is above 1.5 or PWV above 6mm or more than 50% of frames rejected
    ; - 2 : seeing is above 1.2 or PWV above 4.5mm or more than 25% of frames rejected
    ; - 3 : otherwise
    qua_flg = 3
    if data[i_f].seeing gt 1.2 or data[i_f].smttau gt 4.5 or data[i_f].nfr_rej / n_null0 gt 0.25 then qua_flg = 2
    if data[i_f].seeing gt 1.5 or data[i_f].smttau gt 6 or data[i_f].nfr_rej / n_null0 gt 0.50 then qua_flg = 1
    data[i_f].qua_flg = qua_flg

    ; Save OB file
    if not keyword_set(no_save) and drs.null_mode ne 1 then begin
      ; Output file (new file if not /RENEW)
      ob_file0 = sav_int_path + 'null_ob' + string(ob_id, format = '(I03)') + '.sav'
      ; If exist, archiv
      if file_test(ob_file0) then begin
        i_red = 1
        ob_file = sav_int_path + 'null_ob' + string(ob_id, format = '(I03)') + '_v1.sav'
        while file_test(ob_file) eq 1 do begin
          ob_file = sav_int_path + 'null_ob' + string(ob_id, format = '(I03)') + '_v' + string(++i_red, format = '(I0)') + '.sav'
        endwhile
        file_move, ob_file0, ob_file
      endif else i_red = 0
      ; If /RENEW, erase the old one (and use previous ob_file0 defined above if doesn't exist)
      if keyword_set(renew) then begin
        files_ob = file_search(data_path, 'null_ob' + string(ob_id, format = '(I03)') + '*.sav', count = n_files_ob)
        ; Loop over files of this OB until matching paramaters are found (or not)
        for i_fob = 0, n_files_ob - 1 do begin
          ; Restore old file
          restore, files_ob[i_fob]
          if TAG_EXIST(data_ob, 'bfl_mode') then begin
            ; Compare to running/current parameters
            if drs_ob.null_mode eq drs.null_mode and drs_ob.keep_ratio eq drs.keep_ratio and drs_ob.min_fr eq drs.min_fr and drs_ob.n_btstrp eq drs.n_btstrp and drs_ob.nsc_bfac eq drs.nsc_bfac and $
              drs_ob.nsc_bins eq drs.nsc_bins and drs_ob.nsc_cube[0] eq drs.nsc_cube[0] and drs_ob.nsc_cube[1] eq drs.nsc_cube[1] and drs_ob.nsc_cube[2] eq drs.nsc_cube[2] and drs_ob.nsc_omin eq drs.nsc_omin and $
              drs_ob.null_cor eq drs.null_cor and drs_ob.null_rad eq drs.null_rad and drs_ob.null_lim[0] eq drs.null_lim[0] and drs_ob.null_lim[1] eq drs.null_lim[1] and drs_ob.null_range[0] eq drs.null_range[0] and $
              drs_ob.null_range[1] eq drs.null_range[1] and drs.n_frob eq drs.n_frob and $ ; end of drs parameters
              data_ob.objname eq objname and data_ob.lam_cen eq lam_cen and data_ob.bck_mode eq bck_mode and data_ob.bfl_mode eq bfl_mode and data_ob.fit_mode eq fit_mode and data_ob.img_mode eq img_mode and $
              data_ob.ob_mode eq ob_mode and data_ob.pixscale eq pixscale and data_ob.precision eq precisio and data_ob.bck_irad eq bck_irad and data_ob.bck_orad eq bck_orad then begin ; end of paramaters on L1 images
              ob_file0 = files_ob[i_fob]
            endif
          endif
        endfor
      endif
      ; Save file
      data_ob = data[i_f]
      drs_ob = drs
      save, drs_ob, data_ob, filename = ob_file0
    endif

    ; Jump point if restored OB
    restored_ob:

    ; Use NSC mode results if NSC_MODE is 1 (only one is saved in the L1 summary file but they can be reproduced without recomputing everything)
    ; Add error on photometry (not included in NSC) and error on background subtraction
    err_phot = abs(data[i_f].null_avg) * (data[i_f].phot_err / data[i_f].phot_avg) ; error on photometry not included in NSC (correlated between all OBs of the same pointing!)
    err_bckg = data[i_f].bckg_err ; error on background subtraction from complementary nod
    err_trm = sqrt(err_phot ^ 2 + err_bckg ^ 2)
    if drs.null_mode eq 2 then begin
      case drs.nsc_mode of
        0: null_data = data[i_f].null_star.null_opt ; use optimum value
        1: null_data = data[i_f].null_star.null_bay ; use bayesian value
        2: message, 'Bootstrap computation is depreciated. Use NSC_MODE=1 for best results.' ; null_data = data[i_f].null_star.null_avg   ; use average of bootsrap samples
        3: message, 'Bootstrap computation is depreciated. Use NSC_MODE=1 for best results.' ; null_data = data[i_f].null_star.null_mod   ; use mode of bootstrap sampls
        else: message, 'Unknown NSC mode. Must be 0, 1, 2, or 3.'
      endcase
    endif else null_data = data[i_f].null_star.null_avg

    ; Now parse to null results
    err_nsc = null_data.nas_err_sup > null_data.nas_err_low ; Take MAX of error LOW and error SUP
    data[i_f].null_meas = null_data.nas
    data[i_f].null_err = sqrt(err_nsc ^ 2 + err_trm ^ 2) ; Add error on photometry (not included in NSC) and error on background subtraction
    data[i_f].nsc_phavg = null_data.mu
    data[i_f].nsc_phavg_err = null_data.mu_err_sup > null_data.mu_err_low
    data[i_f].nsc_phrms = null_data.sig
    data[i_f].nsc_phrms_err = null_data.sig_err_sup > null_data.sig_err_sup
    data[i_f].nsc_kdrk = null_data.kdrk
    data[i_f].nsc_kdrk_err = null_data.kdrk_err_sup > null_data.kdrk_err_sup
    data[i_f].nsc_chi2 = null_data.chi2
    data[i_f].nsc_chi2_err = null_data.chi2_err_sup
    data[i_f].null_snr = data[i_f].null_meas / data[i_f].null_err
    data[i_f].null_offset = data[i_f].null_meas - data[i_f].null_star.null_opt.nas

    ; Compute expected error terms
    if data[i_f].bck_mode ge 0 then null_err_phot = sqrt(err_bckg ^ 2 + data[i_f].bias_err ^ 2) else null_err_phot = sqrt(2) * err_bckg
    print, 'BCKG MODE : ' + data[i_f].bck_mode
    null_err_phase = sqrt(4 * (data[i_f].nsc_phavg) ^ 2 * (data[i_f].nsc_phrms) ^ 2 + 2 * (data[i_f].nsc_phrms) ^ 4) / 4

    ; Print onfo
    if info gt 0 then begin
      print, string(ob_id, format = '(I03)'), '  ', string(data[i_f].objname, format = '(A8)'), '  ', string(data[i_f].flag, format = '(A3)'), '  ', string(data[i_f].lbt_utc, format = '(A11)'), '  ', string(data[i_f].lbt_alt, format = '(F6.2)'), '  ', string(data[i_f].lbt_para, format = '(F7.2)'), '  ', string(data[i_f].seeing, format = '(F4.2)'), '  ', string(data[i_f].smttau, format = '(F4.2)'), $
        '  ', string(1d+3 * data[i_f].int_time, format = '(F5.1)'), '  ', string(1d+6 * data[i_f].lam_cen, format = '(F5.1)'), '  ', string(1d+6 * data[i_f].bandwidth, format = '(F4.1)'), ' || ', string(data[i_f].photdx_snr, format = '(F5.1)'), ' ', string(data[i_f].photsx_snr, format = '(F5.1)'), $
        '  ', string(data[i_f].int_err, format = '(F7.4)'), '  ', string(data[i_f].rms_tt[0], format = '(F5.1)'), ' ', string(data[i_f].rms_tt[1], format = '(F5.1)'), '  ', string(100. * data[i_f].photdx_rms / data[i_f].phot_avg, format = '(F4.2)'), '  ', string(100. * data[i_f].photsx_rms / data[i_f].phot_avg, format = '(F4.2)'), $
        '  ', string(100. * data[i_f].bckg_rms, format = '(F5.2)'), '  ', string(100. * data[i_f].bias_unc, format = '(F5.1)'), '  ', string(100. * data[i_f].bias_cor, format = '(F5.2)'), ' || ', string(100. * data[i_f].null_meas, format = '(F5.2)'), '  ', string(100. * data[i_f].null_offset, format = '(F5.2)'), ' || ', string(100. * err_phot, format = '(F5.2)'), '  ', string(100. * err_bckg, format = '(F5.2)'), '  ', string(100. * err_nsc, format = '(F5.2)'), $
        ' | ', string(100. * data[i_f].null_err, format = '(F5.2)'), ' || ', string(100. * null_err_phot, format = '(F5.2)'), ' ', string(100. * null_err_phase, format = '(F5.2)'), ' ||', string(1d6 * data[i_f].nsc_phavg * lam_cen / (2 * !dpi), format = '(F6.2)'), '  ', string(1d6 * data[i_f].nsc_phrms * lam_cen / (2 * !dpi), format = '(F5.2)'), '  ', string(data[i_f].nsc_kdrk, format = '(F4.2)'), '  ', string(data[i_f].nsc_chi2, format = '(F6.2)'), '   ', string(data[i_f].nfr_rej, format = '(I04)'), '  ', string(data[i_f].nfr_ob, format = '(I04)')
      if lun gt 0 then $
        printf, lun, string(ob_id, format = '(I03)'), '  ', string(data[i_f].objname, format = '(A8)'), '  ', string(data[i_f].flag, format = '(A3)'), '  ', string(data[i_f].lbt_utc, format = '(A11)'), '  ', string(data[i_f].lbt_alt, format = '(F6.2)'), '  ', string(data[i_f].lbt_para, format = '(F7.2)'), '  ', string(data[i_f].seeing, format = '(F4.2)'), '  ', string(data[i_f].smttau, format = '(F4.2)'), $
        '  ', string(1d+3 * data[i_f].int_time, format = '(F5.1)'), '  ', string(1d+6 * data[i_f].lam_cen, format = '(F5.1)'), '  ', string(1d+6 * data[i_f].bandwidth, format = '(F4.1)'), ' || ', string(data[i_f].photdx_snr, format = '(F5.1)'), ' ', string(data[i_f].photsx_snr, format = '(F5.1)'), $
        '  ', string(data[i_f].int_err, format = '(F7.4)'), '  ', string(data[i_f].rms_tt[0], format = '(F5.1)'), ' ', string(data[i_f].rms_tt[1], format = '(F5.1)'), '  ', string(100. * data[i_f].photdx_rms / data[i_f].phot_avg, format = '(F4.2)'), '  ', string(100. * data[i_f].photsx_rms / data[i_f].phot_avg, format = '(F4.2)'), $
        '  ', string(100. * data[i_f].bckg_rms, format = '(F5.2)'), '  ', string(100. * data[i_f].bias_unc, format = '(F5.1)'), '  ', string(100. * data[i_f].bias_cor, format = '(F5.2)'), ' || ', string(100. * data[i_f].null_meas, format = '(F5.2)'), '  ', string(100. * data[i_f].null_offset, format = '(F5.2)'), ' || ', string(100. * err_phot, format = '(F5.2)'), '  ', string(100. * err_bckg, format = '(F5.2)'), '  ', string(100. * err_nsc, format = '(F5.2)'), $
        ' | ', string(100. * data[i_f].null_err, format = '(F5.2)'), ' || ', string(100. * null_err_phot, format = '(F5.2)'), ' ', string(100. * null_err_phase, format = '(F5.2)'), ' ||', string(1d6 * data[i_f].nsc_phavg * lam_cen / (2 * !dpi), format = '(F6.2)'), '  ', string(1d6 * data[i_f].nsc_phrms * lam_cen / (2 * !dpi), format = '(F5.2)'), '  ', string(data[i_f].nsc_kdrk, format = '(F4.2)'), '  ', string(data[i_f].nsc_chi2, format = '(F6.2)'), '   ', string(data[i_f].nfr_rej, format = '(I04)'), '  ', string(data[i_f].nfr_ob, format = '(I04)')
    endif

    ; Jump point if no nulling data in this OB
    skip_ob:
  endfor

  ; Remove skipped frames from data
  idx_ok = where(data.mjd_obs ne 0, n_ok)
  if n_ok gt 0 then data = data[idx_ok]

  ; Save L1 file
  if not keyword_set(no_save) then begin
    ; Create directories if non existent
    sav_path = pth.l1Fits_path + pth.sep + drs.date_obs + drs.dir_label + pth.sep
    if not file_test(sav_path) then file_mkdir, sav_path
    spawn, 'chmod -R 775 ' + sav_path ; Create directory if it does not exist

    ; If exist, archive it
    outfile = sav_path + 'UT' + drs.date_obs + '.fits'
    savfile = sav_path + 'UT' + drs.date_obs + '.sav'
    if file_test(outfile) then begin
      i_red = 1
      exist = 1
      while exist eq 1 do begin
        oldfile = sav_path + 'UT' + drs.date_obs + '_v' + string(i_red, format = '(I0)') + '.fits'
        if not file_test(oldfile) then exist = 0 else i_red = i_red + 1
      endwhile
      if file_test(outfile) then file_move, outfile, sav_path + 'UT' + drs.date_obs + '_v' + string(i_red, format = '(I0)') + '.fits'
      if file_test(savfile) then file_move, savfile, sav_path + 'UT' + drs.date_obs + '_v' + string(i_red, format = '(I0)') + '.sav'
    endif else i_red = 0

    ; Save data table
    LBTI_SAVEL1SUM, data, outfile = outfile, file_version = i_red + 1

    ; Also save data as an IDL file
    save, data, filename = savfile
  endif
end
