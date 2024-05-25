;+
; NAME: LBTI_SAVEL1VIS
; 
; PURPOSE:
;   This procedure saves data tables in level 1 FITS files.
;
; INPUTS:
;   vis_in        :  On input, the measured visibilities
;   hdr_in        :  On input, a structure with header information
;   data_in       :  On input, a structure that contains information related to each frame
;
; KEYWORDS:
;   FILE_ID       :  Set this keyword to the file ID number to be printed in the file name (none by default)
;                 :  If it's an array of n elements, n identical files will be created (each with a different file id)
;   OUTFILE       :  Name of the output file (superseed the automatic file name)
;
; MODIFICATION HISTORY:
;   Version 1.0, 02-NOV-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (adapted from former routine 'remove_bck.pro')
;   Version 1.1, 04-APR-2015, DD: Added nodding period
;   Version 1.2, 24-MAY-2024, DD: Update file permission

PRO LBTI_SAVEL1VIS, vis_in, hdr_in, data_in, OUTFILE=outfile, FILE_ID=file_id

; Define operational parameters
COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs, log

; Keyword sanity check
IF NOT KEYWORD_SET(FILE_ID) THEN file_id = [0] 
  
; Derive the name tag for files
CASE ABS(hdr_in.obstype) OF
  0: tag = 'PHOT'
  1: tag = 'COHE'
  2: tag = 'NULL'
  3: tag = 'BCKG'
  ELSE: tag = 'UND'
ENDCASE

; Derive which beam is this
IF hdr_in.beam_id EQ 0 THEN BEGIN
  xcen = data_in.xcen_sx
  ycen = data_in.ycen_sx
  IF ABS(hdr_in.obstype) EQ 0 THEN tag = 'PHOT1'
ENDIF ELSE BEGIN
  xcen = data_in.xcen_dx
  ycen = data_in.ycen_dx
  IF ABS(hdr_in.obstype) EQ 0 THEN tag = 'PHOT2'
ENDELSE

; Number of elements and aperture radii
size_data = SIZE(vis_in.vis)
n         = size_data[1]
n_wav     = N_ELEMENTS(vis_in.wav)

; Nod id
nod_id = INTARR(n) + hdr_in.nod_id
                                  
; Define output file name
IF NOT KEYWORD_SET(OUTFILE) THEN BEGIN
  ; Create directory
  sav_path = pth.l1fits_path + drs.date_obs + pth.sep
  IF NOT FILE_TEST(sav_path) THEN BEGIN FILE_MKDIR, sav_path & SPAWN, 'chmod 775 ' + sav_path
  ; Create file name
  id_string = '_ID' + STRING(file_id[0], FORMAT='(I02)') + '_'
  outfile   = STRCOMPRESS(sav_path + 'UT' + STRTRIM(drs.date_obs, 2) + id_string + STRTRIM(hdr_in.flag, 2) + '_' + STRTRIM(hdr_in.objname, 2) + '_DIT-' + STRING(1D+3*hdr_in.int_time, FORMAT='(I0)') + 'ms_' + $
              STRING(1D+6*hdr_in.lam_cen, FORMAT='(I0)') + 'um_' + tag + '.fits' , /REMOVE_ALL)
ENDIF
            
; Add comment to the header
FXHMAKE,  hdr, /INIT, /EXTEND, 0
FXADDPAR, hdr, "COMMENT", "Level 1 FITS file of uncalibrated data taken with LBTI/" + STRTRIM(hdr_in.instrum) + "."
FXADDPAR, hdr, "COMMENT", "The raw data has been reduced using the nodrs software version " +  STRING(hdr_in.drs_vers, FORMAT='(F3.1)') + "."
FXADDPAR, hdr, "COMMENT", "Contact: lbtisupp@ipac.caltech.edu."

; Add keyword definitions to header
FXADDPAR, hdr, "COMMENT", "Target and observation parameters", BEFORE='DATE_OBS'
FXADDPAR, hdr, 'DATE_OBS',  drs.date_obs,        'date of observation'
FXADDPAR, hdr, 'TELESCOP',  'LBT',               'telescope'
FXADDPAR, hdr, 'INSTRUME',  hdr_in.instrum,      'instrument'
FXADDPAR, hdr, 'PID',       hdr_in.pid,          'Project ID'
FXADDPAR, hdr, 'OBJNAME',   hdr_in.objname,      'target name'
FXADDPAR, hdr, 'OBSTYPE',   hdr_in.obstype,      '0=phot, 1=coh, 2=null, 3=bckg'
FXADDPAR, hdr, 'FLAG',      hdr_in.flag,         'SCI, DARK, FLAT, or CAL'
FXADDPAR, hdr, 'OB_ID',     file_id[0],          'OB identification number'
FXADDPAR, hdr, 'LBT_UTC',   data_in[0].lbt_utc,  'UTC from observatory'
FXADDPAR, hdr, 'LBT_RA',    hdr_in.lbt_ra,       'RA from observatory'
FXADDPAR, hdr, 'LBT_DEC',   hdr_in.lbt_dec,      'DEC from observator'
FXADDPAR, hdr, 'LBT_ALT',   data_in[0].lbt_alt,  'ALT from observatory at start of obs.'
FXADDPAR, hdr, 'LBT_AZ',    data_in[0].lbt_az,   'AZ from observatory at start of obs.'
FXADDPAR, hdr, 'LBT_PARA',  data_in[0].lbt_para, 'PA from observatory at start of obs.'
SXADDPAR, hdr, 'NOD_FRQ',   hdr_in.nod_frq,      'Mean nodding period [s]'
SXADDPAR, hdr, 'N_FRBCK',   FIX(hdr_in.n_frbck), 'Number of frames per background'
FXADDPAR, hdr, 'XCEN',      MEAN(xcen),          'Mean initial x-coordinate on detector [pix]'
FXADDPAR, hdr, 'YCEN',      MEAN(ycen),          'Mean initial y-coordinate on detector [pix]'

; Detector parameter
FXADDPAR, hdr, "COMMENT", "Detector parameters", AFTER='YCEN'
FXADDPAR, hdr, 'N_XPIX',    hdr_in.n_xpix,    'initial number of detector rows [pix]'
FXADDPAR, hdr, 'N_YPIX',    hdr_in.n_ypix,    'initial number of detector columns [pix]'
FXADDPAR, hdr, 'PIXSCALE',  hdr_in.pixscale,  'pixel scale [arcsec]'
FXADDPAR, hdr, 'DETMODE',   hdr_in.detmode,   'detector mode (0:HG, 1:LG)'
FXADDPAR, hdr, 'SMPLMODE',  hdr_in.smplmode,  'sampling mode'
FXADDPAR, hdr, 'EXPTIME',   hdr_in.int_time,  'integration time [s]'
FXADDPAR, hdr, 'ACQTIME',   hdr_in.acq_time,  'acquisition time [s]'
FXADDPAR, hdr, 'N_COADD',   hdr_in.n_coadd,   'number of coadds'
FXADDPAR, hdr, 'PAGAIN',    hdr_in.pagain,    'pre-amp gain'
FXADDPAR, hdr, 'PABANDW',   hdr_in.pabandw,   'pre-amp bandwidth [kHz]'
FXADDPAR, hdr, 'DETBIAS',   hdr_in.detbias,   'detector bias [V]'
FXADDPAR, hdr, 'EPERADU',   hdr_in.eperadu,   'number of electron per ADU'

; Filter information
FXADDPAR, hdr, "COMMENT", "Filter information", AFTER='EPERADU'
FXADDPAR, hdr, 'WAVELENG',  hdr_in.lam_cen,   'central wavelength [m]'
FXADDPAR, hdr, 'BANDWIDT',  hdr_in.bandwidth, 'bandwidth [m]'
FXADDPAR, hdr, 'UBC_DXSP',  hdr_in.ubc_dxsp,  'UBC DX field stop'
FXADDPAR, hdr, 'UBC_SXSP',  hdr_in.ubc_sxsp,  'UBC SX field stop'
FXADDPAR, hdr, 'NIC_FSTP',  hdr_in.nic_fstp,  'NIC field stop'
FXADDPAR, hdr, 'NIC_BEAM',  hdr_in.nic_beam,  'NIC beam diverter (0=out, 1=in)'
FXADDPAR, hdr, 'LMIR_FW1',  hdr_in.lmir_fw1,  'LMIRcam filter wheel 1'
FXADDPAR, hdr, 'LMIR_FW2',  hdr_in.lmir_fw2,  'LMIRcam filter wheel 2'
FXADDPAR, hdr, 'LMIR_FW3',  hdr_in.lmir_fw3,  'LMIRcam filter wheel 3'
FXADDPAR, hdr, 'LMIR_FW4',  hdr_in.lmir_fw4,  'LMIRcam filter wheel 4'
FXADDPAR, hdr, 'NOM_FW1',   hdr_in.nom_fw1,   'NOMIC filter wheel 1'
FXADDPAR, hdr, 'NOM_FW2',   hdr_in.nom_fw2,   'NOMIC filter wheel 2'
FXADDPAR, hdr, 'NOM_APW',   hdr_in.nom_apw,   'NOMIC aperture wheel'
FXADDPAR, hdr, 'PHA_FW1',   hdr_in.pha_fw1,   'PHASECAM filter wheel 1'
FXADDPAR, hdr, 'PHA_FW2',   hdr_in.pha_fw2,   'PHASECAM filter wheel 2'
FXADDPAR, hdr, 'PHA_IMG',   hdr_in.pha_img,   'PHASECAM pupil/star imaging wheel'

; AO keywords
SXADDPAR, hdr, "COMMENT", "AO information", AFTER='PHA_IMG'
SXADDPAR, hdr, 'DAOMODE',  hdr_in.daomode,         'Mode of the DX AO'
SXADDPAR, hdr, 'DAOSTRHL', MEAN(data_in.daostrhl), 'DX strehl ratio at V-band'
SXADDPAR, hdr, 'DCMODES',  hdr_in.dcmodes,         'DX number of Zernike modes'
SXADDPAR, hdr, 'DLOOPON',  data_in[0].dloopon,     'DX AO loop on/off'
SXADDPAR, hdr, 'DLOOPGN',  hdr_in.dloopgn,         'DX AO loop gain'
SXADDPAR, hdr, 'DWFSCFRQ', hdr_in.dwfscfrq,        'DX wavefront sensor speed [Hz]'
SXADDPAR, hdr, 'DWFSCBIN', hdr_in.dwfscbin,        'DX wavefront sensor binning'
SXADDPAR, hdr, 'SAOMODE',  hdr_in.saomode,         'Mode of the SX AO'
SXADDPAR, hdr, 'SAOSTRHL', MEAN(data_in.saostrhl), 'SX strehl ratio at V-band'
SXADDPAR, hdr, 'SCMODES',  hdr_in.scmodes,         'SX number of Zernike modes'
SXADDPAR, hdr, 'SLOOPON',  data_in[0].sloopon,     'SX AO loop on/off'
SXADDPAR, hdr, 'SLOOPGN',  hdr_in.sloopgn,         'SX AO loop gain'
SXADDPAR, hdr, 'SWFSCFRQ', hdr_in.dwfscfrq,        'SX wavefront sensor speed [Hz]'
SXADDPAR, hdr, 'SWFSCBIN', hdr_in.dwfscbin,        'SX wavefront sensor binning'

; Phasecam information
FXADDPAR, hdr, "COMMENT", "PHASECam information", AFTER='DWFSCBIN'
FXADDPAR, hdr, 'PCCLOSED', MEAN(data_in.pcclosed), 'Mean status of the PLC (0=open, 1=close)'
IF TAG_EXIST(data_in, 'fpcpistm') THEN BEGIN
  FXADDPAR, hdr, 'SPCPIST',  MEAN(data_in.spcpist),  'Mean SPC piston [um]'
  FXADDPAR, hdr, 'SPCAZ',    MEAN(data_in.spcaz),    'Mean SPC azimuth [mas]'
  FXADDPAR, hdr, 'SPCEL',    MEAN(data_in.spcel),    'Mean SPC elevation [mas]'
  FXADDPAR, hdr, 'SPCAZ',    MEAN(data_in.spcaz),    'Mean SPC azimuth [mas]'
ENDIF ELSE BEGIN
  FXADDPAR, hdr, 'SPCPIST',  -1,    'Mean SPC piston [um]'
  FXADDPAR, hdr, 'SPCAZ',    -1,    'Mean SPC azimuth [mas]'
  FXADDPAR, hdr, 'SPCEL',    -1,    'Mean SPC elevation [mas]'
  FXADDPAR, hdr, 'SPCAZ',    -1,    'Mean SPC azimuth [mas]'
ENDELSE

; Weather conditions
idx_keep = WHERE(data_in.seeing GT 0, n_keep)
IF n_keep GT 5 THEN BEGIN
  FXADDPAR, hdr, "COMMENT", "Weather information", AFTER='SPCEL'
  FXADDPAR, hdr, 'SEEING',    MEDIAN(data_in.seeing[idx_keep]),    'Median seeing over sequence'
  FXADDPAR, hdr, 'SMTTAU',    MEDIAN(data_in.smttau[idx_keep]),    'Median PWV over sequence'
  FXADDPAR, hdr, 'LBTTEMP',   MEDIAN(data_in.lbttemp[idx_keep]),   'Median dome temperature over sequence [C]'
  FXADDPAR, hdr, 'WINDDIR',   MEDIAN(data_in.winddir[idx_keep]),   'Mean wind direction (east of north) [deg]'
  FXADDPAR, hdr, 'WINDSPD',   MEDIAN(data_in.windspd[idx_keep]),   'Mean wind speed [m/s]'
ENDIF ELSE BEGIN
  FXADDPAR, hdr, "COMMENT", "Weather information", AFTER='SPCEL'
  FXADDPAR, hdr, 'SEEING',    -1,   'Median seeing over sequence'
  FXADDPAR, hdr, 'SMTTAU',    -1,   'Median PWV over sequence'
  FXADDPAR, hdr, 'LBTTEMP',   -1,   'Median dome temperature over sequence [C]'
  FXADDPAR, hdr, 'WINDDIR',   -1,   'Mean wind direction (east of north) [deg]'
  FXADDPAR, hdr, 'WINDSPD',   -1,   'Mean wind speed [m/s]'
ENDELSE

; Data reduction information
FXADDPAR, hdr, "COMMENT", "Data reduction information", AFTER='WINDSPD'
FXADDPAR, hdr, 'DRS_VERS',  hdr_in.drs_vers,       'DRS version for image calibration'
FXADDPAR, hdr, 'DRS_DATE',  hdr_in.drs_date,       'DRS version date for image calibration'
FXADDPAR, hdr, 'DATE_CAL',  hdr_in.date_red,       'Date of image calibration'
FXADDPAR, hdr, 'DATE_FLX',  systime(),             'Date of flux computation'
FXADDPAR, hdr, 'BCK_MODE',  FIX(hdr_in.bck_mode),  'Background subtraction mode'
FXADDPAR, hdr, 'OB_MODE',   FIX(drs.ob_mode),      'OB mode'
FXADDPAR, hdr, 'FLX_MODE',  FIX(drs.flx_mode),     '0=aperture phot, 1=PSF fitting'
FXADDPAR, hdr, 'FIT_MODE',  FIX(drs.fit_mode),     'Centroid fitting technique (0=no centroid)'
FXADDPAR, hdr, 'N_CLIP',    FIX(drs.n_clip),       'Size of the extracted images around [XCEN,YCEN]'
FXADDPAR, hdr, 'N_BIN',     FIX(drs.n_bin),        'Image re-binning size'
FXADDPAR, hdr, 'NO_BPM',    FIX(hdr_in.no_bpm),    'Bad pixel correction off if 1 '
FXADDPAR, hdr, 'NO_DARK',   FIX(hdr_in.no_dark),   'Dark subtraction off if 1 '
FXADDPAR, hdr, 'NO_FLAT',   FIX(hdr_in.no_flat),   'Flat correction off if 1'
FXWRITE, outfile, hdr

; Create header for the main table
n_row = N_ELEMENTS(data_in.file_id)
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'RESULTS_SUMMARY', 'results of the data reduction'

; Init column number
IF TAG_EXIST(data_in, 'fpcpistm') THEN n_col = 20 ELSE n_col = 8
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, 0L,              'FILE_ID',  'Raw file identification number' 
FXBADDCOL, 2L, hdr, 0D,              'MJD_OBS',  'Modified Julian Date of observation'
FXBADDCOL, 3L, hdr, 0,               'NOD_ID',   'Nod identification number'
FXBADDCOL, 4L, hdr, 0,               'CHP_ID',   'Chop identification number'
FXBADDCOL, 5L, hdr, 0.,              'XCEN',     'Star initial X position', TUNIT='pix'
FXBADDCOL, 6L, hdr, 0.,              'YCEN',     'Star initial Y position', TUNIT='pic'
FXBADDCOL, 7L, hdr, FLTARR(n_wav),   'VIS_RAW',  'Measured raw visibility'
FXBADDCOL, 8L, hdr, FLTARR(n_wav),   'PHA_RAW',  'Measured raw phase'
IF TAG_EXIST(data_in, 'fpcpistm') THEN BEGIN
  FXBADDCOL, 9L, hdr, 0.,           'FPCPISTM','FPC average piston', TUNIT='um'
  FXBADDCOL, 10L, hdr, 0.,           'FPCPISTS','FPC RMS piston', TUNIT='um' 
  FXBADDCOL, 11L, hdr, 0.,           'FPCAZM',  'FPC average azimuth', TUNIT='arcsec'
  FXBADDCOL, 12L, hdr, 0.,           'FPCELM',  'FPC average elevation', TUNIT='arcsec'
  FXBADDCOL, 13L, hdr, 0.,           'FPCAZS',  'FPC RMS azimuth', TUNIT='arcsec'
  FXBADDCOL, 14L, hdr, 0.,           'FPCELS',  'FPC RMS elevation', TUNIT='arcsec'
  FXBADDCOL, 15L, hdr, 0.,           'PCPHMEAN','Average measured piston', TUNIT='um' 
  FXBADDCOL, 16L, hdr, 0.,           'PCPHSTD', 'RMS measured piston', TUNIT='um' 
  FXBADDCOL, 17L, hdr, 0.,           'PCPHMCOS','Mean cosinus(phi) over the last X ms'
  FXBADDCOL, 18L, hdr, 0.,           'PCPHMSIN','Mean sinus(phi) over the last X ms'
  FXBADDCOL, 19L, hdr, 0,            'PCFJMPS', 'Number of CG jumps over the last X ms'  
  FXBADDCOL, 20L, hdr, 0.,           'PCMSNR',  'K-band FFT-peak SNR' 
ENDIF

; Write extension header to FITS file
FXBCREATE, unit, outfile, hdr
IF TAG_EXIST(data_in, 'fpcpistm') $
 THEN FXBWRITM,  unit, col, LONG(data_in.file_id), DOUBLE(data_in.mjd_obs), FIX(nod_id), FIX(data_in.chp_id), FLOAT(xcen), FLOAT(ycen), $
                 TRANSPOSE(FLOAT(vis_in.vis)), TRANSPOSE(FLOAT(vis_in.phase)), FLOAT(data_in.fpcpistm), FLOAT(data_in.fpcpists), FLOAT(data_in.fpcazm), FLOAT(data_in.fpcelm), $
                 FLOAT(data_in.fpcazs), FLOAT(data_in.fpcels), FLOAT(data_in.pcphmean), FLOAT(data_in.pcphstd), FLOAT(data_in.pcphmcos), FLOAT(data_in.pcphmsin), FIX(data_in.pcfjmps), FLOAT(data_in.pcmsnr) $
 ELSE FXBWRITM,  unit, col, LONG(data_in.file_id), DOUBLE(data_in.mjd_obs), FIX(nod_id), FIX(data_in.chp_id), FLOAT(xcen), FLOAT(ycen), TRANSPOSE(FLOAT(vis_in.vis)), TRANSPOSE(FLOAT(vis_in.phase))
FXBFINISH, unit

; Copy files if FILE_ID is an array
n_id = N_ELEMENTS(file_id)
FOR i_id = 1, n_id-1 DO BEGIN
  id_string = '_ID' + STRING(file_id[i_id], FORMAT='(I02)') + '_'
  outfile2  = STRCOMPRESS(sav_path + 'UT' + STRTRIM(drs.date_obs, 2) + id_string + STRTRIM(hdr_in.flag, 2) + '_' + STRTRIM(hdr_in.objname, 2) + '_DIT-' + STRING(1D+3*hdr_in.int_time, FORMAT='(I0)') + 'ms_' + $
              STRING(1D+6*hdr_in.lam_cen, FORMAT='(I0)') + 'um_' + tag + '.fits' , /REMOVE_ALL)
  FILE_COPY, outfile, outfile2, /OVERWRITE
ENDFOR
END