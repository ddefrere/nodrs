;+
; NAME: LBTI_SAVEL1SUM
; 
; PURPOSE:
;   This procedure saves reduced images in level 1 FITS files.
;
; INPUTS:
;   data          :
;
; KEYWORDS:
;   FILE_VERSION  :  File version number
;   OUTFILE       :  Name of the output file (superseed the automatic file name)
;
;
; MODIFICATION HISTORY:
;   Version 1.0, 02-NOV-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (adapted from former routine 'remove_bck.pro')
;   Version 1.1, 12-NOV-2013, DD: implemented fixed column format following discussions with Mhiseh and Rafael.
;   Version 1.2, 28-MAR-2014, DD: added background and OPD RMS to output files
;   Version 1.3, 31-MAR-2014, DD: added nod ID number to output files
;   Version 1.3, 26-MAY-2014, DD: added new phasecam outputs
;   Version 1.3, 31-MAY-2014, DD: added background bias
;   Version 1.4, 12-AUG-2014, DD: updated header after revision with RMG
;   Version 1.5, 28-OCT-2014, DD: now properly update PID
;   Version 1.6, 28-NOV-2014, DD: added results from the NSC reduction
;   Version 1.7, 23-JAN-2015, DD: added error on phasecam telemetry to output data
;   Version 1.8, 23-FEB-2015, DD: added average and rms of the raw null
;   Version 1.9, 13-MAR-2015, DD: added NSC keywords to output file header
;   Version 2.0, 05-APR-2015, DD: cleaned header for major archiv release
;   Version 2.1, 15-APR-2015, DD: corrected BCK_MODE (now from data rather than config file!)
;   Version 2.2, 07-MAY-2015, DD: removed duplicate columns
;   Version 2.3, 06-AUG-2015, DD: added more NSC parameter information to header
;   Version 2.4, 07-OCT-2015, DD: added plc_wav to PHASECam infos
;   Version 2.5, 20-OCT-2015, DD: added quality flag
;   Version 2.6, 03-DEC-2015, DD: added plc_spec to output structure
;   Version 2.7, 04-DEC-2015, DD: moved APER_RAD, BCK_IRAD, and BCK_ORAD from the header to columns
;   Version 2.8, 23-DEC-2015, DD: Added precision and image mode to output data
;   Version 2.9, 05-FEB-2016, DD: Added pointing ID to output structure
;   Version 3.0, 08-NOV-2016, DD: Added error on photometry

PRO LBTI_SAVEL1SUM, data, FILE_VERSION=file_version, OUTFILE=outfile
   
; Define operational parameters
COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs, log
                                              
; Define output file name
IF NOT KEYWORD_SET(OUTFILE) THEN BEGIN
  ; Create directory
  sav_path = pth.l1fits_path + pth.sep + drs.date_obs + drs.dir_label + pth.sep
  IF NOT FILE_TEST(sav_path) THEN FILE_MKDIR, sav_path
  ; Create file name
  outfile  = sav_path + 'UT' + drs.date_obs + '_SUM.fits'
ENDIF
                
; Create new FITS file with primary HDU
FXHMAKE,  hdr, /INIT, /EXTEND, 0
FXADDPAR, hdr, "COMMENT", "This is a FITS file of uncalibrated (level 1) data taken with LBTI/" + STRTRIM(data[0].instrum[0]) + "."
FXADDPAR, hdr, "COMMENT", "The raw data has been reduced using the nodrs software version " +  STRING(drs.version, FORMAT='(F3.1)') + "." 
FXADDPAR, hdr, "COMMENT", "Contact: lbtisupp@ipac.caltech.edu."
FXADDPAR, hdr, 'DATE_OBS',  drs.date_obs,                        'Date of observation'
FXADDPAR, hdr, 'TELESCOP', 'LBT',                                'Telescope'
FXADDPAR, hdr, 'INSTRUME',  data[0].instrum[0],                  'Instrument'
FXADDPAR, hdr, 'DRS_VERS',  drs.version,                         'DRS version ID'
FXADDPAR, hdr, 'DRS_DATE',  drs.date,                            'DRS version date'
FXADDPAR, hdr, 'FILE_VER',  file_version,                        'L1 file version number'
FXADDPAR, hdr, 'DATE_RED',  systime(),                           'Date of data reduction'
FXADDPAR, hdr, 'OB_MODE',   FIX(data[0].ob_mode[0]),             '0: all, 1: nod splitted, 2: user defined'
FXADDPAR, hdr, 'BCK_MODE',  FIX(data[0].bck_mode[0]),            '0: none, 1: nod pairs, 2: closest n frames, 3: dedicated, 3 chopping'
FXADDPAR, hdr, 'BFL_MODE',  FIX(data[0].bfl_mode[0]),            'Background floor mode (0: mean, 1: median)'
FXADDPAR, hdr, 'FLX_MODE',  FIX(data[0].flx_mode[0]),            '0: aperture phot, 2: weigh. aper. phot, 3: PSF fitting'
FXADDPAR, hdr, 'FIT_MODE',  FIX(data[0].fit_mode[0]),            '0: no centroid, 1-3=Gaussian fit, 4=Lorentzian fit, 5=Moffat fit'
FXADDPAR, hdr, 'FRA_MODE',  FIX(data[0].fra_mode[0]),            '0: flux from background subtracted images, 1: flux from raw images'
FXADDPAR, hdr, 'IMG_MODE',  FIX(data[0].img_mode[0]),            'Image combination mode (median, mean, mean+)'
FXADDPAR, hdr, 'NUL_MODE',  FIX(drs.null_mode),                  '0=mode, 1=%best, 2=statistical'
FXADDPAR, hdr, 'NULL_COR',  FIX(drs.null_cor),                   'Subtract the null estimated from high-frequency phase noise'
null_lim = '[' + STRING(1D2*drs.null_lim[0], FORMAT='(F4.1)') + ',' + STRING(1D2*drs.null_lim[1], FORMAT='(F4.1)') + ']' 
FXADDPAR, hdr, 'NULL_LIM',  null_lim,                            'Acceptable raw null range (before NSC and in %) '
FXADDPAR, hdr, 'NBIN_FAC' , FIX(drs.nsc_bfac),                   'Mulitplier on number of histogram bins (nbins_fac*sqrt(Np))'
FXADDPAR, hdr, 'N_BTSTRP',  FIX(drs.n_btstrp),                   'Number of bootstrap samples to compute the error bar'
nsc_cube = STRING(drs.nsc_cube[0], FORMAT='(I0)') + 'x' + STRING(drs.nsc_cube[1], FORMAT='(I0)') + 'x' + STRING(drs.nsc_cube[2], FORMAT='(I0)') + 'x' + STRING(drs.nsc_cube[3], FORMAT='(I0)')  + 'x' + STRING(drs.nsc_cube[4], FORMAT='(I0)') 
FXADDPAR, hdr, 'NSC_CUBE',  nsc_cube,                            'Cube size for NSC reduction (null, phase offset, phase rms, OPD seed, background factor)'
FXADDPAR, hdr, 'NSC_BINS' , FIX(drs.nsc_bins),                   'Bin size for NSC reduction (0: constant, 1: variable)'
FXADDPAR, hdr, 'NSC_MODE' , FIX(drs.nsc_mode),                   'NSC null mode (0: mean, 1: mode)'
FXADDPAR, hdr, 'NSC_OMIN' , FIX(drs.nsc_omin),                   'Minimum number of occurences per bin for the fit'
FXADDPAR, hdr, 'NOD_FRQ',   data[0].nod_frq[0],                  '[s]    Mean nodding period'
FXADDPAR, hdr, 'N_FRBCK',   FIX(data[0].n_frbck[0]),             'Number of preserved frames for background subtraction'
FXADDPAR, hdr, 'NULL_RAD',  FIX(drs.null_rad),                   'Photometric radius parameter (0 for EEID)'
FXADDPAR, hdr, 'PRECISIO',  FIX(data[0].precision[0]),           'Image processing precision (0: float, 1: double)'
FXADDPAR, hdr, 'R_FRNULL',  drs.keep_ratio,                      'Ratio of preserved frames for null computation'
FXWRITE, outfile, hdr

; Create header for the main table
n_row = N_ELEMENTS(data.ob_id)
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'RESULTS_SUMMARY', 'results of the data reduction'

; Init column number
n_col = 50
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, LONG(0)               , 'PID',          'Project ID number'  
FXBADDCOL, 2L, hdr, LONG(0)               , 'FILE_ID',      'File identification number'
FXBADDCOL, 3L, hdr, LONG(0)               , 'OB_ID',        'OB identification number'
FXBADDCOL, 4L, hdr, LONG(0)               , 'NOD_ID',       'Nod identification number'
FXBADDCOL, 5L, hdr, LONG(0)               , 'PT_ID',        'Pointing identification number'
FXBADDCOL, 6L, hdr, 'aaaaaaaaaaaaaaaaaaaa', 'OBJNAME',      'Object name'
FXBADDCOL, 7L, hdr, 'SCI'                 , 'FLAG',         'SCI/CAL'
FXBADDCOL, 8L, hdr, DOUBLE(0)             , 'MJD_OBS',      'Modified Julian Date of observation'
FXBADDCOL, 9L, hdr, '00:00:00.00'         , 'LBT_UTC',      'UTC from observatory'
FXBADDCOL, 10L, hdr, '00:00:00.00'        , 'LBT_LST',      'LST from observatory'
FXBADDCOL, 11L, hdr, '+00:00:00.000'      , 'LBT_RA',       'RA from observatory'
FXBADDCOL, 12L, hdr, '+00:00:00.000'      , 'LBT_DEC',      'DEC from observatory'
FXBADDCOL, 13L, hdr, FLOAT(0)             , 'LBT_ALT',      'ALT from observatory'
FXBADDCOL, 14L, hdr, FLOAT(0)             , 'LBT_AZ',       'AZ from observatory'
FXBADDCOL, 15L, hdr, FLOAT(0)             , 'LBT_PARA',     'Parralactic angle from obs.'
FXBADDCOL, 16L, hdr, FLOAT(0)             , 'WAV_EFF',      'Central wavelength', TUNIT='m'
FXBADDCOL, 17L, hdr, FLOAT(0)             , 'BANDWIDTH',    'Bandwidth', TUNIT='m'
FXBADDCOL, 18L, hdr, FIX(0)               , 'APER_RAD',     'Radius for aperture photometry', TUNIT='pix'
FXBADDCOL, 19L, hdr, FIX(0)               , 'BCK_IRAD',     'Inner radius for background computation', TUNIT='pix'
FXBADDCOL, 20L, hdr, FIX(0)               , 'BCK_ORAD',     'Outer radius for background computation', TUNIT='pix'
FXBADDCOL, 21L, hdr, FLOAT(0)             , 'PHOT_AVG',     'Mean combined photometry', TUNIT='ADU'
FXBADDCOL, 22L, hdr, FLOAT(0)             , 'PHOT_ERR',     'Error on mean combined photometry', TUNIT='ADU'
FXBADDCOL, 23L, hdr, FLOAT(0)             , 'PHOTDX_AVG',   'Mean DX photometry', TUNIT='ADU'
FXBADDCOL, 24L, hdr, FLOAT(0)             , 'PHOTDX_RMS',   'RMS DX photometry', TUNIT='ADU'
FXBADDCOL, 25L, hdr, FLOAT(0)             , 'PHOTDX_SNR',   'SNR on DX photometry'
FXBADDCOL, 26L, hdr, FLOAT(0)             , 'PHOTSX_AVG',   'Mean SX photometry', TUNIT='ADU'
FXBADDCOL, 27L, hdr, FLOAT(0)             , 'PHOTSX_RMS',   'RMS SX photometry', TUNIT='ADU' 
FXBADDCOL, 28L, hdr, FLOAT(0)             , 'PHOTSX_SNR',   'SNR on SX photometry'
FXBADDCOL, 29L, hdr, FLOAT(0)             , 'TTDX_RMS',     'RMS DX tip/tilt ', TUNIT='mas' 
FXBADDCOL, 30L, hdr, FLOAT(0)             , 'TTSX_RMS',     'RMS SX tip/tilt', TUNIT='mas' 
FXBADDCOL, 31L, hdr, FLOAT(0)             , 'BCKG_MEAS',    'Mean background per pixel', TUNIT='ADU'
FXBADDCOL, 32L, hdr, FLOAT(0)             , 'BCKG_MEAS_RMS','RMS background per pixel', TUNIT='ADU'
FXBADDCOL, 33L, hdr, FLOAT(0)             , 'BCKG_BIAS',    'Empty-region background per photometric aperture (relative to star)'
FXBADDCOL, 34L, hdr, FLOAT(0)             , 'BCKG_BIAS_RMS','Empty-region background RMS per photometric aperture (relative to star)'
FXBADDCOL, 35L, hdr, FLOAT(0)             , 'BCKG_CNOD_RMS','Background RMS per photometric aperture in complementary nod (relative to star)'
FXBADDCOL, 36L, hdr, DOUBLE(0)            , 'NULL_MEAS',    'Measured raw null'
FXBADDCOL, 37L, hdr, DOUBLE(0)            , 'NULL_MEAS_ERR','Error on the measured raw null'
FXBADDCOL, 38L, hdr, DOUBLE(0)            , 'NULL_OFFSET',  'Measured null offset (NSC only)'
FXBADDCOL, 39L, hdr, FLOAT(0)             , 'NULL_MEAS_SNR','SNR on the null'  
FXBADDCOL, 40L, hdr, FLOAT(0)             , 'NULL_MEAS_AVG','Average of the raw null (non filtered)'  
FXBADDCOL, 41L, hdr, FLOAT(0)             , 'NULL_MEAS_RMS','Standard deviation of the raw null (non filtered)'  
FXBADDCOL, 42L, hdr, FLOAT(0)             , 'NSC_CHI2',     'Chi2 (NSC reduction only)'
FXBADDCOL, 43L, hdr, FLOAT(0)             , 'NSC_CHI2_ERR', 'Error on chi2 (NSC reduction only)'
FXBADDCOL, 44L, hdr, FLOAT(0)             , 'NSC_PHAVG',    'Best-fit phase AVG (NSC reduction only)', TUNIT='rad'  
FXBADDCOL, 45L, hdr, FLOAT(0)             , 'NSC_PHAVG_ERR','Error on best-fit phase AVG', TUNIT='rad'
FXBADDCOL, 46L, hdr, FLOAT(0)             , 'NSC_PHRMS',    'Best-fit phase RMS (NSC reduction only)', TUNIT='rad' 
FXBADDCOL, 47L, hdr, FLOAT(0)             , 'NSC_PHRMS_ERR','Error on best-fit phase RMS', TUNIT='rad'   
FXBADDCOL, 48L, hdr, LONG(0)              , 'NFR_OB',       'Number of null frames'
FXBADDCOL, 49L, hdr, LONG(0)              , 'NFR_REJ',      'Number of rejected null frames' 
FXBADDCOL, 50L, hdr, FIX(0)               , 'QUALITY_FLG',  'Quality flag'

; Write extension header to FITS file
FXBCREATE, unit, outfile, hdr
FXBWRITM,  unit, col, LONG(data.pid), LONG(data.file_id), LONG(data.ob_id), LONG(data.nod_id), LONG(data.pt_id), data.objname, data.flag, DOUBLE(data.mjd_obs), data.lbt_utc, data.lbt_lst, data.lbt_ra, data.lbt_dec, FLOAT(data.lbt_alt), FLOAT(data.lbt_az), FLOAT(data.lbt_para), FLOAT(data.lam_cen), FLOAT(data.bandwidth), $
                      FIX(data.aper_rad), FIX(data.bck_irad), FIX(data.bck_orad), FLOAT(data.phot_avg), FLOAT(data.phot_err), FLOAT(data.photdx_avg), FLOAT(data.photdx_rms), FLOAT(data.photdx_snr), FLOAT(data.photsx_avg), FLOAT(data.photsx_rms), FLOAT(data.photsx_snr), FLOAT(data.rms_tt[0]), FLOAT(data.rms_tt[1]), FLOAT(data.bckg_meas), $
                      FLOAT(data.bckg_emeas), FLOAT(data.bias_cor), FLOAT(data.bias_err), FLOAT(data.bckg_rms), DOUBLE(data.null_meas), DOUBLE(data.null_err), DOUBLE(data.null_offset), FLOAT(data.null_snr), FLOAT(data.null_avg), FLOAT(data.null_rms), FLOAT(data.nsc_chi2), FLOAT(data.nsc_chi2_err), FLOAT(data.nsc_phavg), FLOAT(data.nsc_phavg_err), $
                      FLOAT(data.nsc_phrms), FLOAT(data.nsc_phrms_err), LONG(data.nfr_ob), LONG(data.nfr_rej), FIX(data.qua_flg)
FXBFINISH, unit  

; Create header for the detector table
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'DETECTOR_INFO', 'main parameters of the detector'

; Init column number
n_col = 14
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, LONG(0)               , 'OB_ID',     'OB identification number'
FXBADDCOL, 2L, hdr, LONG(0)               , 'N_XPIX',    'X-size of detector', TUNIT='pix'
FXBADDCOL, 3L, hdr, LONG(0)               , 'N_YPIX',    'Y-size of detector', TUNIT='pix'
FXBADDCOL, 4L, hdr, LONG(0)               , 'XCEN',      'Initial x position on detector', TUNIT='pix'
FXBADDCOL, 5L, hdr, LONG(0)               , 'YCEN',      'Initial y position on detector', TUNIT='pix'
FXBADDCOL, 6L, hdr, FLOAT(0)              , 'PIXSCALE',  'Final pixel scale', TUNIT='arcsec'
FXBADDCOL, 7L, hdr, FLOAT(0)              , 'INT_TIME',  'Integration time', TUNIT='s'
FXBADDCOL, 8L, hdr, FLOAT(0)              , 'ACQ_TIME',  'Acquisition time', TUNIT='s'
FXBADDCOL, 9L, hdr, LONG(0)               , 'N_COADD',   'Number of coaads'
FXBADDCOL, 10L, hdr, LONG(0)              , 'SMPLMODE',  'Sampling mode'
FXBADDCOL, 11L, hdr, LONG(0)              , 'DETMODE',   'Detector mode (0: HG, 1: LG)'
FXBADDCOL, 12L, hdr, LONG(0)              , 'PAGAIN',    'Pre-amp gain'
FXBADDCOL, 13L, hdr, LONG(0)              , 'PABANDW',   'Pre-amp bandwidth'
FXBADDCOL, 14L, hdr, FLOAT(0)             , 'DETBIAS',   'Detector bias', TUNIT='V'
FXBADDCOL, 15L, hdr, LONG(0)              , 'EPERADU',   'Electrons per ADU'

; Write extension header to FITS file
FXBCREATE, unit, outfile, hdr
FXBWRITM,  unit, col, LONG(data.ob_id), LONG(data.n_xpix), LONG(data.n_ypix), LONG(data.xcen), LONG(data.ycen), FLOAT(data.pixscale), FLOAT(data.int_time), FLOAT(data.acq_time), LONG(data.n_coadd), $
                      LONG(data.smplmode), LONG(data.detmode), LONG(data.pagain), LONG(data.pabandw), FLOAT(data.detbias), LONG(data.eperadu)
FXBFINISH, unit

; Create header for the AO table
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'AO_INFO', 'main running parameters of the AO systems'

; Init column number
n_col = 15
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, LONG(0)               , 'OB_ID',     'OB identification number'
FXBADDCOL, 2L, hdr, 'aaaaaaaaaa'          , 'DAOMODE',   'DX AO mode'
FXBADDCOL, 3L, hdr, FLOAT(0)              , 'DAOSTRHL',  'DX Strehl Ratio at V-band', TUNIT='%'
FXBADDCOL, 4L, hdr, FIX(0)                , 'DCMODES',   'DX number of Zernike modes'
FXBADDCOL, 5L, hdr, 'aaaaaaaaaa'          , 'DLOOPON',   'DX AO loop on/off'
FXBADDCOL, 6L, hdr, FLOAT(0)              , 'DLOOPGN',   'DX AO loop gain'
FXBADDCOL, 7L, hdr, FLOAT(0)              , 'DWFSCFRQ',  'DX Wavefront sensor speed', TUNIT='Hz'
FXBADDCOL, 8L, hdr, LONG(0)               , 'DWFSCBIN',  'DX Wavefront CCD binning'
FXBADDCOL, 9L, hdr, 'aaaaaaaaaa'          , 'SAOMODE',   'SX AO mode'
FXBADDCOL, 10L, hdr, FLOAT(0)             , 'SAOSTRHL',  'SX Strehl Ratio at V-band', TUNIT='%'
FXBADDCOL, 11L, hdr, FIX(0)               , 'SCMODES',   'SX number of Zernike modes'
FXBADDCOL, 12L, hdr, 'aaaaaaaaaa'         , 'SLOOPON',   'SX AO loop on/off'
FXBADDCOL, 13L, hdr, FLOAT(0)             , 'SLOOPGN',   'SX AO loop gain'
FXBADDCOL, 14L, hdr, FLOAT(0)             , 'SWFSCFRQ',  'SX Wavefront sensor speed', TUNIT='Hz'
FXBADDCOL, 15L, hdr, LONG(0)              , 'SWFSCBIN',  'SX Wavefront CCD binning'                   
                                                                                         
; Write extension header to FITS file
FXBCREATE, unit, outfile, hdr
FXBWRITM,  unit, col, LONG(data.ob_id), STRING(data.daomode), FLOAT(data.daostrhl), FIX(data.dcmodes), STRING(data.dloopon), FLOAT(data.dloopgn), FLOAT(data.dwfscfrq), LONG(data.dwfscbin), STRING(data.saomode), $
                      FLOAT(data.saostrhl), FIX(data.scmodes), STRING(data.sloopon), FLOAT(data.sloopgn), FLOAT(data.swfscfrq), LONG(data.swfscbin)
FXBFINISH, unit

; Create header for the filter table
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'FILTER_INFO', 'filter wheel position'

; Init column number
n_col = 15
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, LONG(0)               , 'OB_ID',     'OB identification number'
FXBADDCOL, 2L, hdr, 'aaaaaaaaaaaaaaaaaaaa', 'UBC_DXSP',  'UBC DX field stop'
FXBADDCOL, 3L, hdr, 'aaaaaaaaaaaaaaaaaaaa', 'UBC_SXSP',  'UBC SX field stop'
FXBADDCOL, 4L, hdr, 'aaaaaaaaaaaaaaaaaaaa', 'NIC_FSTP',  'NIC Field Stop'
FXBADDCOL, 5L, hdr, 'aaaaaaaaaaaaaaaaaaaa', 'NIC_BEAM',  'NIC Beam Diverter'
FXBADDCOL, 6L, hdr, 'aaaaaaaaaaaaaaaaaaaa', 'LMIR_FW1',  'LMIRCam Filter Wheel 1'
FXBADDCOL, 7L, hdr, 'aaaaaaaaaaaaaaaaaaaa', 'LMIR_FW2',  'LMIRCam Filter Wheel 2'
FXBADDCOL, 8L, hdr, 'aaaaaaaaaaaaaaaaaaaa', 'LMIR_FW3',  'LMIRCam Filter Wheel 3'
FXBADDCOL, 9L, hdr, 'aaaaaaaaaaaaaaaaaaaa', 'LMIR_FW4',  'LMIRCam Filter Wheel 4'
FXBADDCOL, 10L, hdr, 'aaaaaaaaaaaaaaaaaaaa','NOM_FW1',   'NOMIC Filter Wheel 1'
FXBADDCOL, 11L, hdr, 'aaaaaaaaaaaaaaaaaaaa','NOM_FW2',   'NOMIC Filter Wheel 1'
FXBADDCOL, 12L, hdr, 'aaaaaaaaaaaaaaaaaaaa','NOM_APW',   'NOMIC aperture Wheel'
FXBADDCOL, 13L, hdr, 'aaaaaaaaaaaaaaaaaaaa','PHA_FW1',   'PHASECam Filter Wheel 1'
FXBADDCOL, 14L, hdr, 'aaaaaaaaaaaaaaaaaaaa','PHA_FW2',   'PHASECam Filter Wheel 2'
FXBADDCOL, 15L, hdr, 'aaaaaaaaaaaaaaaaaaaa','PHA_IMG',   'PHASECam imaging wheel'

; Write extension header to FITS file
FXBCREATE, unit, outfile, hdr
FXBWRITM,  unit, col, LONG(data.ob_id), data.ubc_dxsp, data.ubc_sxsp, data.nic_fstp, data.nic_beam, data.lmir_fw1, data.lmir_fw2, data.lmir_fw3, data.lmir_fw4, data.nom_fw1, data.nom_fw2, data.nom_apw, data.pha_fw1, data.pha_fw2, data.pha_img;, data.nil_dic  
FXBFINISH, unit

;Create header for the PHASECam table
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'PHASECAM_INFO', 'PLC performance indicators'

; Init column number
n_col = 22
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, LONG(0)                ,  'OB_ID',        'OB identification number'
FXBADDCOL, 2L, hdr, LONG(0)                ,  'PCCLOSED',     'Status of the PLC'
FXBADDCOL, 3L, hdr, FLOAT(0)               ,  'PLC_WAV',      'Effective wavelength', TUNIT='um'
FXBADDCOL, 4L, hdr, '                    ' ,  'PLC_SPEC',     'Spectrum file used for PCWAV'
FXBADDCOL, 5L, hdr, FLOAT(0)               ,  'PCPLSTP',      'Mean pathlength setpoint', TUNIT='deg'
FXBADDCOL, 6L, hdr, FLOAT(0)               ,  'PCTIPSP',      'Mean tip setpoint', TUNIT='mas'
FXBADDCOL, 7L, hdr, FLOAT(0)               ,  'PCTLTSP',      'Mean tilt setpoint', TUNIT='mas'
FXBADDCOL, 8L, hdr, FLOAT(0)               ,  'SPC_PIST',     'SPC piston', TUNIT='um'
FXBADDCOL, 9L, hdr, FLOAT(0)               ,  'SPC_AZ',       'SPC azimuth', TUNIT='arcsec'
FXBADDCOL, 10L, hdr, FLOAT(0)              ,  'SPC_EL',       'SPC elevation', TUNIT='arcsec'
FXBADDCOL, 11L, hdr, FLOAT(0)              ,  'FPC_PISTM',    'FPC average piston', TUNIT='um'
FXBADDCOL, 12L, hdr, FLOAT(0)              ,  'FPC_PISTS',    'FPC RMS piston', TUNIT='um'
FXBADDCOL, 13L, hdr, FLOAT(0)              ,  'FPC_AZM',      'FPC average azimuth', TUNIT='arcsec'
FXBADDCOL, 14L, hdr, FLOAT(0)              ,  'FPC_ELM',      'FPC average elevation', TUNIT='arcsec'
FXBADDCOL, 15L, hdr, FLOAT(0)              ,  'FPC_AZS',      'FPC RMS azimuth', TUNIT='arcsec'
FXBADDCOL, 16L, hdr, FLOAT(0)              ,  'FPC_ELS',      'FPC RMS elevation', TUNIT='arcsec'
FXBADDCOL, 17L, hdr, FLOAT(0)              ,  'PCPHMEAN',     'AVG piston over the last X ms', TUNIT='um'
FXBADDCOL, 18L, hdr, FLOAT(0)              ,  'PCPHMEAN_ERR', 'RMS of AVG piston over the last X ms', TUNIT='um'
FXBADDCOL, 19L, hdr, FLOAT(0)              ,  'PCPHSTD',      'RMS piston over the last X ms', TUNIT='um'
FXBADDCOL, 20L, hdr, FLOAT(0)              ,  'PCPHSTD_ERR',  'RMS of RMS piston over the last X ms', TUNIT='um'
FXBADDCOL, 21L, hdr, FLOAT(0)              ,  'PCMSNR',       'AVG K-band FFT-peak SNR'
FXBADDCOL, 22L, hdr, LONG(0)               , 'DITH_PER',      'Dithering period [in number of frames]'

; Write extension header to FITS file
FXBCREATE, unit, outfile, hdr
FXBWRITM,  unit, col, LONG(data.ob_id), LONG(data.pcclosed), FLOAT(data.plc_wav), data.plc_spec, FLOAT(data.pcplsp), FLOAT(data.pctipsp), FLOAT(data.pctltsp), FLOAT(data.spc_pist), FLOAT(data.spc_az), FLOAT(data.spc_el), $
                      FLOAT(data.fpc_pistm), FLOAT(data.fpc_pists), FLOAT(data.fpc_azm), FLOAT(data.fpc_azs), FLOAT(data.fpc_elm), FLOAT(data.fpc_els), FLOAT(data.pcphmean), FLOAT(data.pcphmean_err), $
                      FLOAT(data.pcphstd), FLOAT(data.pcphstd_err), FLOAT(data.pcmsnr), LONG(data.dith_per)  
FXBFINISH, unit

;Create header for the Wheather table
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'WHEATHER_INFO', 'Wheather information'

; Init column number
n_col = 6
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, LONG(0)               , 'OB_ID',     'OB identification number'
FXBADDCOL, 2L, hdr, FLOAT(0)              , 'SEEING',    'Seeing',  TUNIT='arcsec'
FXBADDCOL, 3L, hdr, FLOAT(0)              , 'SMTTAU',    'Precipitable water vapor estimate',  TUNIT='mm'
FXBADDCOL, 4L, hdr, FLOAT(0)              , 'LBTTEMP',   'Temperature in Dome',  TUNIT='C'
FXBADDCOL, 5L, hdr, FLOAT(0)              , 'WINDDIR',   'Wind direction (degrees E of N)'
FXBADDCOL, 6L, hdr, FLOAT(0)              , 'WINDSPD',   'Wind speed',  TUNIT='m/s'

; Write extension header to FITS file
FXBCREATE, unit, outfile, hdr
FXBWRITM,  unit, col, LONG(data.ob_id), FLOAT(data.seeing), FLOAT(data.smttau), FLOAT(data.lbttemp), FLOAT(data.winddir), FLOAT(data.windspd)
FXBFINISH, unit

END