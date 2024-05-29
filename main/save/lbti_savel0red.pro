;+
; NAME: LBTI_SAVEL0RED
; 
; PURPOSE:
;   This procedure saves calibrated L0 images and relevant keywords in FITS files.
;
; INPUTS:
;   img_in     :  3-dimension array with the calibrated images (n_xpix, n_ypix, n_frames).
;   hdr_in     :  On input, the structure with the corresponding header information
;
; KEYWORDS:
;   SAVEMEDIAN : Set this keyword to save the median frame as well (for 1000 256x256 frames, this takes 6 seconds! The rest takes ~1sec.)
;   SUB_DIR    : Set this keyword to the sub-directory where the data have be to be saved
;   
; MODIFICATION HISTORY:
;   Version 1.0, 18-SEP-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 27-SEP-2014, DD: added a few keywords to output files
;   Version 1.2, 21-OCT-2014, DD: minor bur corrected
;   Version 1.2, 28-OCT-2014, DD: added porject ID to header
;   Version 1.3, 30-OCT-2014, DD: now save important data in the log
;   Version 1.4, 07-NOV-2014, DD: added slope to output structure
;   Version 1.5, 11-NOV-2014, DD: Added PHASECam setpoint to output header
;   Version 1.6, 12-JAN-2015, DD: Added NO_DARK and NO_FLAT to output header
;   Version 1.7  16-FEB-2015, DD: Added new PHASECam fieds (required by NSC)
;   Version 1.8  17-FEB-2015, DD: Now save datalog file in an external function
;   Version 1.9  20-FEB-2015, DD: Removed slope from the output data
;   Version 2.0  04-APR-2015, DD: Added nodding period and number of frames in background
;   Version 2.1, 06-NOV-2015, DD: Added config ID
;   Version 2.2, 02-DEC-2015, DD: Added pointing ID and background frame selection mode to output header
;   Version 2.3, 03-DEC-2015, DD: Added keyword SAVEMEDIAN + moved PHASECam setpoint to first extension to prepare for PWV tracking
;   Version 2.4, 21-DEC-2015, DD: Added precision and image mode to output header + now use drs.precision for /SAVEMEDIAN
;   Version 2.5, 23-DEC-2015, DD: Added /NO_COPY to MWRFITS to save memory
;   Version 2.6, 21-JAN-2016, DD: Added FWHMX, FWHMY, and SLOPE keywords to output header
;   Version 2.7, 12-FEB-2016, DD: Added PHASECam information for beam 2
;   Version 2.8, 03-JUL-2016, DD: Added central value
;   Version 2.9, 20-OCT-2016, DD: Adapted EPERADU for the new preamp
;   Version 3.0, 26-OCT-2016, DD: Added SLOPE rms columns in output data (and removed unused chop_id one)
;   Version 3.1, 30-JAN-2017, DD: Moved weather information to header
;   Version 3.2, 12-FEB-2017, DD: Added keyword SUB_DIR
;   Version 3.3, 17-APR-2018, DD: Now save the closed-loop median frame too
;   Version 3.4, 24-MAY-2024, DD: Update file permission

PRO LBTI_SAVEL0RED, img_in, hdr_in, SAVEMEDIAN=savemedian, SUB_DIR=sub_dir

; Define operational parameters
COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs, log
       
; Create file header
MKHDR, hdr, img_in

; Add comment to the header
SXADDPAR, hdr, "COMMENT", "Calibrated level 0 FITS file of data taken with LBTI/" + STRTRIM(hdr_in[0].instrum) + "."  
SXADDPAR, hdr, "COMMENT", "The raw data has been reduced using the nodrs software version " +  STRING(drs.version, FORMAT='(F3.1)') + "."   
SXADDPAR, hdr, "COMMENT", "Contact: lbtisupp@ipac.caltech.edu." 

; Add keyword definitions to header   
SXADDPAR, hdr, "COMMENT", "Target and observation parameters", BEFORE='TELESCOP'
SXADDPAR, hdr, 'TELESCOP',  'LBT',                      'Telescope'
SXADDPAR, hdr, 'INSTRUME',  hdr_in[0].instrum[0],       'Instrument'
SXADDPAR, hdr, 'DATE_OBS',  drs.date_obs,               'Date of observation'
SXADDPAR, hdr, 'OBJNAME',   hdr_in[0].objname[0],       'Target name'
SXADDPAR, hdr, 'LBT_RA',    hdr_in[0].lbt_ra[0],        'RA from observatory at start of obs.'
SXADDPAR, hdr, 'LBT_DEC',   hdr_in[0].lbt_dec[0],       'DEC from observatory at start of obs.'
SXADDPAR, hdr, 'LBT_UTC',   hdr_in[0].lbt_utc[0],       'UT time at start of oservations'
SXADDPAR, hdr, 'DATATYPE',  FIX(hdr_in[0].datatype[0]), 'Datatype'
SXADDPAR, hdr, 'OBSTYPE',   FIX(hdr_in[0].obstype[0]),  'Obstype'
SXADDPAR, hdr, 'CFG_ID',    FIX(hdr_in[0].cfg_id[0]),   'Config ID number'
SXADDPAR, hdr, 'NOD_ID',    FIX(hdr_in[0].nod_id[0]),   'Nod ID number'
SXADDPAR, hdr, 'PID',       FIX(hdr_in[0].pid[0]),      'Project ID'
SXADDPAR, hdr, 'PT_ID',     FIX(hdr_in[0].pt_id[0]),    'Pointing ID number'
SXADDPAR, hdr, 'FLAG',      hdr_in[0].flag[0],          'Flag'

; Detector information
SXADDPAR, hdr, "COMMENT", "Detector information", AFTER='FLAG'
SXADDPAR, hdr, 'PIXSCALE',  cnf.pix_size/1D3,       'pixel scale [arcsec]'
SXADDPAR, hdr, 'DETMODE',   hdr_in[0].detmode[0],   'detector mode (0:HG, 1:LG)'
SXADDPAR, hdr, 'SMPLMODE',  hdr_in[0].smplmode[0],  'sampling mode'
SXADDPAR, hdr, 'EXPTIME',   hdr_in[0].int_time[0],  'integration time [s]'
SXADDPAR, hdr, 'ACQTIME',   hdr_in[0].acq_time[0],  'acquisition time [s]'
SXADDPAR, hdr, 'NCOADDS',   hdr_in[0].n_coadd[0],   'number of coadds'
SXADDPAR, hdr, 'MODE',      hdr_in[0].smplmode[0],  'camera readout mode'
SXADDPAR, hdr, 'PAGAIN',    hdr_in[0].pagain[0],    'pre-amp gain'
SXADDPAR, hdr, 'PABANDW',   hdr_in[0].pabandw[0],   'pre-amp bandwidth [kHz] '
SXADDPAR, hdr, 'PACTCDLY',  hdr_in[0].pactcdly[0],  'pre-amp delay ID'
SXADDPAR, hdr, 'DETBIAS',   hdr_in[0].detbias[0],   'detector bias [V]'
IF STRCOMPRESS(hdr_in[0].instrum[0], /REMOVE_ALL) EQ 'nomic' THEN BEGIN
  IF STRCOMPRESS(hdr_in[0].interfac[0], /REMOVE_ALL) EQ 'Preamp' THEN BEGIN
    IF hdr_in[0].detmode[0] EQ 0 THEN eperadu = cnf.adu2e[2] ELSE eperadu = cnf.adu2e[3]
  ENDIF ELSE BEGIN
    IF hdr_in[0].detmode[0] EQ 0 THEN eperadu = cnf.adu2e[0] ELSE eperadu = cnf.adu2e[1]
  ENDELSE
ENDIF ELSE eperadu = hdr_in[0].eperadu[0]
SXADDPAR, hdr, 'EPERADU',   eperadu,   'number of electron per ADU'


; AO information
SXADDPAR, hdr, "COMMENT", "AO information", AFTER='EPERADU'
SXADDPAR, hdr, 'DAOMODE',  hdr_in[0].daomode[0],    'Mode of the DX AO'
SXADDPAR, hdr, 'DCMODES',  hdr_in[0].dcmodes[0],    'DX number of Zernike modes'
SXADDPAR, hdr, 'DLOOPGN',  hdr_in[0].dloopgn[0],    'DX AO loop gain'
SXADDPAR, hdr, 'DWFSCFRQ', hdr_in[0].dwfscfrq[0],   'DX wavefront sensor speed [Hz]'
SXADDPAR, hdr, 'DWFSCBIN', hdr_in[0].dwfscbin[0],   'DX wavefront sensor binning'
SXADDPAR, hdr, 'SAOMODE',  hdr_in[0].saomode[0],    'Mode of the SX AO'
SXADDPAR, hdr, 'SCMODES',  hdr_in[0].scmodes[0],    'SX number of Zernike modes'
SXADDPAR, hdr, 'SLOOPGN',  hdr_in[0].sloopgn[0],    'SX AO loop gain'
SXADDPAR, hdr, 'SWFSCFRQ', hdr_in[0].swfscfrq[0],   'SX wavefront sensor speed [Hz]'
SXADDPAR, hdr, 'SWFSCBIN', hdr_in[0].swfscbin[0],   'SX wavefront sensor binning'

; Filter information
SXADDPAR, hdr, "COMMENT", "Filter information", AFTER='SWFSCBIN'
SXADDPAR, hdr, 'WAVELENG',  hdr_in[0].lam_cen[0],   'central wavelength [m]'
SXADDPAR, hdr, 'BANDWIDT',  hdr_in[0].bandwidt[0],  'bandwidth [m]'
SXADDPAR, hdr, 'UBC_DXSP',  hdr_in[0].ubc_dxsp[0],  'UBC DX field stop'
SXADDPAR, hdr, 'UBC_SXSP',  hdr_in[0].ubc_sxsp[0],  'UBC SX field stop'
SXADDPAR, hdr, 'NIC_FSTP',  hdr_in[0].nic_fstp[0],  'NIC field stop'
SXADDPAR, hdr, 'NIC_BEAM',  hdr_in[0].nic_beam[0],  'NIC beam diverter (0=out, 1=in)'
SXADDPAR, hdr, 'LMIR_FW1',  hdr_in[0].lmir_fw1[0],  'LMIRcam filter wheel 1'
SXADDPAR, hdr, 'LMIR_FW2',  hdr_in[0].lmir_fw2[0],  'LMIRcam filter wheel 2'
SXADDPAR, hdr, 'LMIR_FW3',  hdr_in[0].lmir_fw3[0],  'LMIRcam filter wheel 3'
SXADDPAR, hdr, 'LMIR_FW4',  hdr_in[0].lmir_fw4[0],  'LMIRcam filter wheel 4'
SXADDPAR, hdr, 'NOM_FW1',   hdr_in[0].nom_fw1[0],   'NOMIC filter wheel 1'
SXADDPAR, hdr, 'NOM_FW2',   hdr_in[0].nom_fw2[0],   'NOMIC filter wheel 2'
SXADDPAR, hdr, 'NOM_APW',   hdr_in[0].nom_apw[0],   'NOMIC aperture wheel'
SXADDPAR, hdr, 'PHA_FW1',   hdr_in[0].pha_fw1[0],   'PHASECAM filter wheel 1'
SXADDPAR, hdr, 'PHA_FW2',   hdr_in[0].pha_fw2[0],   'PHASECAM filter wheel 2'
SXADDPAR, hdr, 'PHA_IMG',   hdr_in[0].pha_img[0],   'PHASECAM pupil/star imaging wheel'
;SXADDPAR, hdr, 'NIL_DIC',   hdr_in[0].nil_dic[0],   '[]    NIL Dichoic position'

; Reduction parameter
SXADDPAR, hdr, "COMMENT", "Reduction parameters", AFTER='PHA_IMG'
SXADDPAR, hdr, 'DRS_VERS',  drs.version,              'version of the DRS'
SXADDPAR, hdr, 'DRS_DATE',  drs.date,                 'DRS version date'
SXADDPAR, hdr, 'DATE_RED',  systime(),                'Date of data reduction'
SXADDPAR, hdr, 'BCK_MODE',  FIX(drs.bckg_mode),       'Background subtraction mode'
SXADDPAR, hdr, 'BCK_NOD',   hdr_in[0].bck_nod[0],     'Parity of nod uses for background subraction'
SXADDPAR, hdr, 'BCK_SEL',   FIX(drs.bckg_sel),        'Background frame selection mode'
SXADDPAR, hdr, 'IMG_MODE',  FIX(drs.img_mode),        'Image combination mode'
SXADDPAR, hdr, 'N_FRBCK',   FIX(hdr_in[0].n_frbck[0]),'Number of frames in background'
SXADDPAR, hdr, 'NO_BPM',    FIX(drs.no_bpm),          'Bad pixel correction off if 1 '
SXADDPAR, hdr, 'NO_DARK',   FIX(drs.no_dark),         'Dark subtraction off if 1 '
SXADDPAR, hdr, 'NO_FLAT',   FIX(drs.no_flat),         'Flat correction off if 1'
SXADDPAR, hdr, 'PRECISIO',  FIX(drs.precision),       '0: float, 1: double'

; Weather condition
SXADDPAR, hdr, "COMMENT", "Reduction parameters", AFTER='PRECISIO'
idx_keep = WHERE(hdr_in.seeing GT 0, n_keep)
IF n_keep GT 5 THEN BEGIN
  SXADDPAR, hdr, 'SEEING',    MEDIAN(hdr_in[idx_keep].seeing),  'Median seeing over sequence'
  SXADDPAR, hdr, 'SMTTAU',    MEDIAN(hdr_in[idx_keep].smttau),  'Median precipatable water vapor from SMT'
  SXADDPAR, hdr, 'LBTTEMP',   MEDIAN(hdr_in[idx_keep].lbttemp), 'Median temperature in Dome [C]'
  SXADDPAR, hdr, 'WINDDIR',   MEDIAN(hdr_in[idx_keep].winddir), 'Median wind direction (east of north) [deg]'
  SXADDPAR, hdr, 'WINDSPD',   MEDIAN(hdr_in[idx_keep].windspd), 'Median wind speed [m/s]'
ENDIF ELSE BEGIN
  SXADDPAR, hdr, 'SEEING',    -1,   'Median seeing over sequence'
  SXADDPAR, hdr, 'SMTTAU',    -1,   'Median PWV over sequence'
  SXADDPAR, hdr, 'LBTTEMP',   -1,   'Median dome temperature over sequence [C]'
  SXADDPAR, hdr, 'WINDDIR',   -1,   'Mean wind direction (east of north) [deg]'
  SXADDPAR, hdr, 'WINDSPD',   -1,   'Mean wind speed [m/s]'
ENDELSE

; Reduction parameter
SXADDPAR, hdr, "COMMENT", "Computed parameters", AFTER='WINDSPD'
SXADDPAR, hdr, 'NOD_FRQ',   hdr_in[0].nod_frq[0],       'Nodding period [s]'
SXADDPAR, hdr, 'FWHMXDX',   MEAN(hdr_in.fwhm_x[1]),     'Mean best-fit FWHM in the x direction (DX, pixels)'
SXADDPAR, hdr, 'FWHMYDX',   MEAN(hdr_in.fwhm_y[1]),     'Mean best-fit FWHM in the y direction (DX, pixels)'
SXADDPAR, hdr, 'FWHMXSX',   MEAN(hdr_in.fwhm_x[0]),     'Mean best-fit FWHM in the x direction (SX, pixels)'
SXADDPAR, hdr, 'FWHMYSX',   MEAN(hdr_in.fwhm_y[0]),     'Mean best-fit FWHM in the y direction (SX, pixels)'
SXADDPAR, hdr, 'SLOPEDX',   MEAN(hdr_in.slope[1]),      'Slope of the Moffat profile DX'
SXADDPAR, hdr, 'SLOPESX',   MEAN(hdr_in.slope[0]),      'Slope of the Moffat profile SX'
SXADDPAR, hdr, 'XCEN_DX',   MEAN(hdr_in.xcen[1]),       'Mean x coordinate DX'
SXADDPAR, hdr, 'YCEN_DX',   MEAN(hdr_in.ycen[1]),       'Mean y coordinate DX'
SXADDPAR, hdr, 'XCEN_SX',   MEAN(hdr_in.xcen[0]),       'Mean x coordinate SX'
SXADDPAR, hdr, 'YCEN_SX',   MEAN(hdr_in.ycen[0]),       'Mean y coordinate SX'

; Define file path/name
IF NOT KEYWORD_SET(SUB_DIR) THEN sub_dir = ' '
sav_path = pth.l0fits_path + hdr_in[0].date_obs[0] + pth.sep + sub_dir + pth.sep 
IF NOT FILE_TEST(sav_path) THEN FILE_MKDIR, sav_path & SPAWN, 'chmod -R 775 ' + sav_path
filename = STRCOMPRESS(sav_path + 'c_' + drs.date_obs + '_N' +  STRING(hdr_in[0].nod_id[0], FORMAT='(I03)') + '_' $
           +  hdr_in[0].flag[0] + '_' +  hdr_in[0].objname[0] + '_' + STRING(MIN(hdr_in.file_id), FORMAT='(I06)') + '-' +  $
           STRING(MAX(hdr_in.file_id), FORMAT='(I06)'), /REMOVE_ALL)
  
; Save average image (do it first because of /NO_COPY below)
n_img = N_ELEMENTS(img_in[0,0,*])
IF n_img GT 1 AND KEYWORD_SET(SAVEMEDIAN) THEN BEGIN
  ; Save the median of the whole cube
  IF drs.img_mode EQ 0 THEN img_med = REFORM(MEDIAN(img_in, DIMENSION=3, /EVEN, DOUBLE=drs.precision))
  IF drs.img_mode EQ 1 THEN img_med = REFORM(MEAN(img_in, DIMENSION=3, DOUBLE=drs.precision))
  IF drs.img_mode EQ 2 THEN RESISTANT_MEAN, img_in, 10, img_med, DIMENSION=3, DOUBLE=drs.precision, /SILENT
  IF drs.img_mode GT 2 THEN MESSAGE, "Unknown image mode. Must be 0,1, or 2"
  MKHDR, hdr2, img_med                          ; Make header for the fits file
  WRITEFITS, filename + '_MED.fits', img_med, hdr2
  ; Save the median of only the closed-loop frames
  obstype = hdr_in[0].obstype
  nfr     = N_ELEMENTS(img_in[0,0,*])
  CASE obstype  OF
    ; 0: PHOTOMETRY. We request that both AO loops are closed
    0: idx_keep = WHERE(hdr_in.dloopon EQ 1 AND hdr_in.sloopon EQ 1, nfr)
    ; 1: FIZEAU. We request that both AO loops are closed
    1: idx_keep = WHERE(hdr_in.dloopon EQ 1 AND hdr_in.sloopon EQ 1, nfr)
    ; 2: NULLING. We request that both AO loops and the phase loop are closed
    2: idx_keep = WHERE(hdr_in.dloopon EQ 1 AND hdr_in.sloopon EQ 1 AND hdr_in.pcclosed EQ 1, nfr)
    ; 3: BACKGROUND. No conditions
    3: idx_keep = LINDGEN(nfr)
    ; 4: TEST. No conditions
    4: idx_keep = LINDGEN(nfr)
    ELSE: MESSAGE, 'Unknown OBSTYPE'
  ENDCASE
  IF nfr GT 1 THEN BEGIN
    IF drs.img_mode EQ 0 THEN img_med = REFORM(MEDIAN(img_in[*,*,idx_keep], DIMENSION=3, /EVEN, DOUBLE=drs.precision))
    IF drs.img_mode EQ 1 THEN img_med = REFORM(MEAN(img_in[*,*,idx_keep], DIMENSION=3, DOUBLE=drs.precision))
    IF drs.img_mode EQ 2 THEN RESISTANT_MEAN, img_in[*,*,idx_keep], 10, img_med, DIMENSION=3, DOUBLE=drs.precision, /SILENT
    IF drs.img_mode GT 2 THEN MESSAGE, "Unknown image mode. Must be 0,1, or 2"
    MKHDR, hdr2, img_med                          ; Make header for the fits file
    WRITEFITS, filename + '_MED-CLOSED.fits', img_med, hdr2
  ENDIF
ENDIF

; Save file
MWRFITS, img_in, filename + '_IMG.fits', hdr, /CREATE, /SILENT, /NO_COPY

; Append extension with parameters that change frame by frame
n_row = N_ELEMENTS(hdr_in.mjd_obs)

; Create extension
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'DATA_SERIES', 'parameters that change frame by frame'

; Init column number (FXADDCOL takes 50 colums max!!!)
n_col = 46
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, hdr_in[0].file_id,    'FILE_ID',    'Raw file identification number'
FXBADDCOL, 2L, hdr, hdr_in[0].filename,   'FILENAME',   'Raw file original name'
FXBADDCOL, 3L, hdr, hdr_in[0].mjd_obs,    'MJD_OBS',    'Modified Julian Date of observation'
FXBADDCOL, 4L, hdr, hdr_in[0].lbt_utc,    'LBT_UTC',    'UTC from observatory'
FXBADDCOL, 5L, hdr, hdr_in[0].lbt_lst,    'LBT_LST',    'LST from observatory'
FXBADDCOL, 6L, hdr, hdr_in[0].lbt_alt,    'LBT_ALT',    'ALT from observatory'
FXBADDCOL, 7L, hdr, hdr_in[0].lbt_az,     'LBT_AZ',     'AZ from observatory'
FXBADDCOL, 8L, hdr, hdr_in[0].lbt_para,   'LBT_PARA',   'Parralactic angle from obs.'
FXBADDCOL, 9L, hdr, hdr_in[0].bck_avg[*], 'BCK_AVG',    'Average background in ADU'
FXBADDCOL, 10L, hdr, hdr_in[0].bck_rms[*],'BCK_RMS',    'RMS background in ADU'
FXBADDCOL, 11L, hdr, hdr_in[0].xcen[0],   'XCEN_SX',    'X position of SX beam'
FXBADDCOL, 12L, hdr, hdr_in[0].xcen[1],   'XCEN_DX',    'X position of DX beam'
FXBADDCOL, 13L, hdr, hdr_in[0].ycen[0],   'YCEN_SX',    'Y position of SX beam'
FXBADDCOL, 14L, hdr, hdr_in[0].ycen[1],   'YCEN_DX',    'Y position of DX beam'
FXBADDCOL, 15L, hdr, hdr_in[0].cv,        'CV',         'Central value [ADU]'
; AO information
FXBADDCOL, 16L, hdr, hdr_in[0].daostrhL,   'DAOSTRHL',  'DX strehl ratio at V-band'
FXBADDCOL, 17L, hdr, hdr_in[0].dloopon,    'DLOOPON',   'DX AO loop on/off'
FXBADDCOL, 18L, hdr, hdr_in[0].dsloprms,   'DSLOPRMS',  'DX slope RMS'
FXBADDCOL, 19L, hdr, hdr_in[0].saostrhL,   'SAOSTRHL',  'SX strehl ratio at V-band'
FXBADDCOL, 20L, hdr, hdr_in[0].sloopon,    'SLOOPON',   'SX AO loop on/off'
FXBADDCOL, 21L, hdr, hdr_in[0].ssloprms,   'SSLOPRMS',  'SX slope RMS'
; PLC information
FXBADDCOL, 22L, hdr, hdr_in[0].pcclosed,   'PCCLOSED',  'status of the PLC (0=open, 1=close)'
FXBADDCOL, 23L, hdr, hdr_in[0].pcplsp,     'PCPLSP',    'Pathlentgh setpoint [deg]'
FXBADDCOL, 24L, hdr, hdr_in[0].pctipsp,    'PCTIPSP',   'Tip setpoint [mas]'
FXBADDCOL, 25L, hdr, hdr_in[0].pctltsp,    'PCTLTSP',   'Tilt setpoint [mas]'
FXBADDCOL, 26L, hdr, hdr_in[0].spdthpos,   'SPDTHPOS',  'OPD dithering pattern in degrees'
FXBADDCOL, 27L, hdr, hdr_in[0].spc_pist,   'SPCPIST',   'SPC piston [um]'
FXBADDCOL, 28L, hdr, hdr_in[0].spc_az,     'SPCAZ',     'SPC azimuth [mas]'
FXBADDCOL, 29L, hdr, hdr_in[0].spc_el,     'SPCEL',     'SPC elevation [mas]'
FXBADDCOL, 30L, hdr, hdr_in[0].fpc_pistm,  'FPCPISTM',  'FPC average piston [um]'
FXBADDCOL, 31L, hdr, hdr_in[0].fpc_pists,  'FPCPISTS',  'FPC RMS piston [um]'
FXBADDCOL, 32L, hdr, hdr_in[0].fpc_azm,    'FPCAZM',    'FPC average azimuth [mas]'
FXBADDCOL, 33L, hdr, hdr_in[0].fpc_azs,    'FPCAZS',    'FPC RMS azimuth [mas]'
FXBADDCOL, 34L, hdr, hdr_in[0].fpc_elm,    'FPCELM',    'FPC average elevation [mas]'
FXBADDCOL, 35L, hdr, hdr_in[0].fpc_els,    'FPCELS',    'FPC RMS elevation [mas]'
FXBADDCOL, 36L, hdr, hdr_in[0].pcphmean,   'PCPHMEAN',  'Phase AVG over the last X ms [um]'
FXBADDCOL, 37L, hdr, hdr_in[0].pcphstd,    'PCPHSTD',   'Phase RMS over the last X ms [um]'
FXBADDCOL, 38L, hdr, hdr_in[0].pcphmcos,   'PCPHMCOS',  'Mean cosinus(phi) over the last X ms'
FXBADDCOL, 39L, hdr, hdr_in[0].pcphmsin,   'PCPHMSIN',  'Mean sinus(phi) over the last X ms'
FXBADDCOL, 40L, hdr, hdr_in[0].pcfjmps,    'PCFJMPS',   'Number of CG jumps over the last X ms'
FXBADDCOL, 41L, hdr, hdr_in[0].pcmsnr,     'PCMSNR',    'Mean fringe SNR over the last X ms'
FXBADDCOL, 42L, hdr, hdr_in[0].pcphmean2,  'PCPHMEAN2', 'Phase AVG over the last X ms [um]'
FXBADDCOL, 43L, hdr, hdr_in[0].pcphstd2,   'PCPHSTD2',  'Phase RMS over the last X ms [um]'
FXBADDCOL, 44L, hdr, hdr_in[0].pcphmcos2,  'PCPHMCOS2', 'Mean cosinus(phi) over the last X ms'
FXBADDCOL, 45L, hdr, hdr_in[0].pcphmsin2,  'PCPHMSIN2', 'Mean sinus(phi) over the last X ms'
FXBADDCOL, 46L, hdr, hdr_in[0].pcmsnr2,    'PCMSNR2',   'Mean fringe SNR over the last X ms'

; Write extension header to FITS file
FXBCREATE, unit, filename + '_IMG.fits', hdr
FXBWRITM,  unit, col, hdr_in.file_id, hdr_in.filename, hdr_in.mjd_obs, hdr_in.lbt_utc, hdr_in.lbt_lst, hdr_in.lbt_alt, hdr_in.lbt_az, hdr_in.lbt_para, $
                      hdr_in.bck_avg, hdr_in.bck_rms, hdr_in.xcen[0], hdr_in.xcen[1], hdr_in.ycen[0], hdr_in.ycen[1], hdr_in.cv, $
                      hdr_in.daostrhl, hdr_in.dloopon, hdr_in.dsloprms, hdr_in.saostrhl, hdr_in.sloopon, hdr_in.ssloprms, $
                      hdr_in.pcclosed, hdr_in.pcplsp, hdr_in.pctipsp, hdr_in.pctltsp, hdr_in.spdthpos, hdr_in.spc_pist, $
                      hdr_in.spc_az, hdr_in.spc_el, hdr_in.fpc_pistm, hdr_in.fpc_pists, hdr_in.fpc_azm, hdr_in.fpc_azs, hdr_in.fpc_elm, hdr_in.fpc_els, $
                      hdr_in.pcphmean, hdr_in.pcphstd, hdr_in.pcphmcos, hdr_in.pcphmsin, hdr_in.pcmsnr, hdr_in.pcfjmps, hdr_in.pcphmean2, hdr_in.pcphstd2, $
                      hdr_in.pcphmcos2, hdr_in.pcphmsin2, hdr_in.pcmsnr2
FXBFINISH, unit
             
END