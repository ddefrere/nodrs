;+
; NAME: LBTI_SAVEL1IMG
;
; PURPOSE:
;   This procedure saves reduced images in level 1 FITS files.
;
; INPUTS:
;   img_in        :  3-dimension array with the calibrated images (n_xpix, n_ypix, n_frames).
;   hdr_in        :  On input, a structure with header information
;   data_in       :  On input, a structure that contains information related to each frame
;
; KEYWORDS:
;   FILE_ID       :  Set this keyword to the file ID number to be printed in the file name (none by default)
;                 :  If it's an array of n elements, n identical files will be created (each with a different file id)
;   SPLIT         :  If set, the image cube and the data table are saved in different files
;
; MODIFICATION HISTORY:
;   Version 1.0, 14-OCT-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (adapted from former routine 'remove_bck.pro')
;   Version 1.1, 02-NOV-2013, DD: implemented new structure for output files
;   Version 1.2, 21-NOV-2013, DD: added keyword TAG
;   Version 1.3, 24-MAY-2013, DD: removed keyword OB_MODE
;   Version 1.4, 25-MAY-2013, DD: enabled FILE_ID as an array
;   Version 1.5, 18-AUG-2014, DD: added 'UT' to output file name
;   Version 1.6, 24-SEP-2014, DD: simplified implementation + removed keywords ONEFILE and TAG
;   Version 1.7, 28-OCT-2013, DD: added project ID to output header
;   Version 1.8, 15-JAN-2015, DD: added NO_FLAT, NO_DARK, and NO_BPM to output header
;   Version 1.9, 04-APR-2015, DD: Added nodding period
;   Version 2.0, 04-MAY-2015, DD: Added background nod
;   Version 2.1, 24-AUG-2015, DD: Added setpoint information to output header
;   Version 2.2, 06-OCT-2015, DD: Added PCWAV and PCSPEC to output header
;   Version 2.2, 23-DEC-2015, DD: Added precision and image mode to output header
;   Version 2.3, 24-DEC-2015, DD: Added /NO_COPY to save memory
;   Version 2.4, 02-JUL-2016, DD: Added central value
;   Version 2.5, 04-APR-2017, DD: Added more target information to output header
;   Version 2.6, 26-FEB-2019, DD: Added flx_out
;   Version 2.7, 15-OCT-2023, DD: Updated for FRA_MODE=2 (i.e., PCA background subtraction)
;   Version 2.8, 24-MAY-2024, DD: Minor fix for column 13 (added transpose) + update file permission
;   Version 2.9, 02-SEP-2024, DD: Added parameter SKY_OFF

pro LBTI_SAVEL1IMG, img_in, hdr_in, data_in, flx_out, file_id = file_id, split = split
  compile_opt idl2

  ; Define operational parameters
  common GLOBAL, prm, cnf, wav, tgt, pth, drs, log

  ; Keyword sanity check
  if not keyword_set(file_id) then file_id = [0]

  ; Derive which beam is this
  tag = 'NULL' ; NULL by default unless belwo
  if hdr_in.beam_id eq 0 then begin
    xcen = data_in.xcen_sx
    ycen = data_in.ycen_sx
    slope = data_in.slope_sx
    if abs(hdr_in.obstype) eq 0 then tag = 'PHOT1'
    if abs(hdr_in.obstype) eq 3 then tag = 'BCKG'
  endif else begin
    xcen = data_in.xcen_dx
    ycen = data_in.ycen_dx
    slope = data_in.slope_dx
    if abs(hdr_in.obstype) eq 0 then tag = 'PHOT2' else tag = 'BCKG'
  endelse

  ; Nod id
  nod_id = replicate(hdr_in.nod_id, n_elements(data_in.xcen_dx))
  chp_id = replicate(1, n_elements(data_in.xcen_dx)) ; not used urrently

  ; Create output path if does not exist
  sav_path = pth.l1Fits_path + pth.sep + drs.date_obs + drs.dir_label + pth.sep
  if not file_test(sav_path) then file_mkdir, sav_path
  spawn, 'chmod -R 775 ' + sav_path

  ; Define output file name
  id_string = '_ID' + string(file_id[0], format = '(I03)') + '_'
  outname = strcompress(sav_path + 'UT' + strtrim(drs.date_obs, 2) + id_string + strtrim(hdr_in.flag, 2) + '_' + strtrim(hdr_in.objname, 2) + '_DIT-' + string(1d+3 * hdr_in.int_time, format = '(I0)') + 'ms_' + $
    string(1d+6 * hdr_in.lam_cen, format = '(I0)') + 'um_' + tag, /remove_all)

  ; If file already exists, open and decide what to do
  outfile = outname + '_IMG.fits'
  if file_test(outfile) and tag eq 'BCKG' then begin
    hdr_res = HEADFITS(outfile, /silent)
    bck_nod1 = fix(FXPAR(hdr_res, 'BCK_NOD1', /nocontinue))
    bck_nod2 = fix(FXPAR(hdr_res, 'BCK_NOD2', /nocontinue))
    ; If different background nod, append. Otherwise, erase.
    if bck_nod1 ne hdr_in.nod_id and bck_nod2 ne hdr_in.nod_id then begin
      bck_nod2 = hdr_in.nod_id
      img_res = MRDFITS(outfile, 0, /silent)
      img_in = [[[img_in]], [[img_res]]]
    endif
  endif else begin
    if tag eq 'BCKG' then bck_nod1 = hdr_in.nod_id else bck_nod1 = hdr_in.bck_nod
    bck_nod2 = 0
  endelse

  ; Create file header
  MKHDR, hdr, img_in

  ; Add comment to the header
  SXADDPAR, hdr, 'COMMENT', 'Level 1 FITS file of uncalibrated data taken with LBTI/' + strtrim(hdr_in.instrum) + '.'
  SXADDPAR, hdr, 'COMMENT', 'The raw data has been reduced using the nodrs software version ' + string(hdr_in.drs_vers, format = '(F3.1)') + '.'
  SXADDPAR, hdr, 'COMMENT', 'Contact: lbtisupp@ipac.caltech.edu.'

  ; Add keyword definitions to header
  SXADDPAR, hdr, 'COMMENT', 'Target and observation parameters', before = 'DATE_OBS'
  SXADDPAR, hdr, 'DATE_OBS', drs.date_obs, 'date of observation'
  SXADDPAR, hdr, 'TELESCOP', 'LBT', 'telescope'
  SXADDPAR, hdr, 'INSTRUME', hdr_in.instrum, 'instrument'
  SXADDPAR, hdr, 'OBJNAME', hdr_in.objname, 'target name'
  SXADDPAR, hdr, 'PID', hdr_in.pid, 'Project ID'
  SXADDPAR, hdr, 'OBSTYPE', hdr_in.obstype, '0=phot, 1=coh, 2=null, 3=bckg'
  SXADDPAR, hdr, 'FLAG', hdr_in.flag, 'SCI, DARK, FLAT, or CAL'
  SXADDPAR, hdr, 'CFG_ID', hdr_in.cfg_id, 'Config identification number'
  SXADDPAR, hdr, 'OB_ID', file_id[0], 'OB identification number'
  SXADDPAR, hdr, 'PT_ID', hdr_in.pt_id, 'Pointing identification number'
  SXADDPAR, hdr, 'LBT_UTC', data_in[0].lbt_utc, 'UTC from observatory'
  SXADDPAR, hdr, 'LBT_RA', hdr_in.lbt_ra, 'RA from observatory'
  SXADDPAR, hdr, 'LBT_DEC', hdr_in.lbt_dec, 'DEC from observator'
  SXADDPAR, hdr, 'LBT_ALT', data_in[0].lbt_alt, 'ALT from observatory at start of obs.'
  SXADDPAR, hdr, 'LBT_AZ', data_in[0].lbt_az, 'AZ from observatory at start of obs.'
  SXADDPAR, hdr, 'LBT_PARA', data_in[0].lbt_para, 'PA from observatory at start of obs.'
  idx = where(tgt.name eq strcompress(strlowcase(string(hdr_in.objname)), /remove_all))
  SXADDPAR, hdr, 'TGT_DIST', tgt[idx].dist, 'Target distance [pc]'
  SXADDPAR, hdr, 'TGT_TEFF', tgt[idx].temp, 'Effective temperature [K]'
  SXADDPAR, hdr, 'TGT_SPEC', tgt[idx].spectrum[0], 'Pickles spectrum'
  SXADDPAR, hdr, 'TGT_LDM', tgt[idx].ldm, 'Limb-darkening diameter [mas]'
  SXADDPAR, hdr, 'TGT_LDE', tgt[idx].lde, 'Error on limb-darkening diameter [mas]'
  SXADDPAR, hdr, 'TGT_CAL4', tgt[idx].calfor, 'Calibrator for which SCI'
  SXADDPAR, hdr, 'NOD_FRQ', hdr_in.nod_frq, 'Mean nodding period [s]'
  SXADDPAR, hdr, 'XCEN', mean(xcen), 'Mean initial x-coordinate on detector [pix]'
  SXADDPAR, hdr, 'YCEN', mean(ycen), 'Mean initial y-coordinate on detector [pix]'

  ; Detector parameter
  SXADDPAR, hdr, 'COMMENT', 'Detector parameters', after = 'YCEN'
  SXADDPAR, hdr, 'N_XPIX', hdr_in.n_xpix, 'initial number of detector rows [pix]'
  SXADDPAR, hdr, 'N_YPIX', hdr_in.n_ypix, 'initial number of detector columns [pix]'
  SXADDPAR, hdr, 'PIXSCALE', hdr_in.pixscale, 'pixel scale [arcsec]'
  SXADDPAR, hdr, 'DETMODE', hdr_in.detmode, 'detector mode (0:HG, 1:LG)'
  SXADDPAR, hdr, 'SMPLMODE', hdr_in.smplmode, 'sampling mode'
  SXADDPAR, hdr, 'EXPTIME', hdr_in.int_time, 'integration time [s]'
  SXADDPAR, hdr, 'ACQTIME', hdr_in.acq_time, 'acquisition time [s]'
  SXADDPAR, hdr, 'N_COADD', hdr_in.n_coadd, 'number of coadds'
  SXADDPAR, hdr, 'PAGAIN', hdr_in.pagain, 'pre-amp gain'
  SXADDPAR, hdr, 'PABANDW', hdr_in.pabandw, 'pre-amp bandwidth [kHz]'
  SXADDPAR, hdr, 'DETBIAS', hdr_in.detbias, 'detector bias [V]'
  SXADDPAR, hdr, 'EPERADU', hdr_in.eperadu, 'number of electron per ADU'

  ; Filter information
  SXADDPAR, hdr, 'COMMENT', 'Filter information', after = 'EPERADU'
  SXADDPAR, hdr, 'WAVELENG', hdr_in.lam_cen, 'central wavelength [m]'
  SXADDPAR, hdr, 'BANDWIDT', hdr_in.bandwidth, 'bandwidth [m]'
  SXADDPAR, hdr, 'UBC_DXSP', hdr_in.ubc_dxsp, 'UBC DX field stop'
  SXADDPAR, hdr, 'UBC_SXSP', hdr_in.ubc_sxsp, 'UBC SX field stop'
  SXADDPAR, hdr, 'NIC_FSTP', hdr_in.nic_fstp, 'NIC field stop'
  SXADDPAR, hdr, 'NIC_BEAM', hdr_in.nic_beam, 'NIC beam diverter (0=out, 1=in)'
  SXADDPAR, hdr, 'LMIR_FW1', hdr_in.lmir_fw1, 'LMIRcam filter wheel 1'
  SXADDPAR, hdr, 'LMIR_FW2', hdr_in.lmir_fw2, 'LMIRcam filter wheel 2'
  SXADDPAR, hdr, 'LMIR_FW3', hdr_in.lmir_fw3, 'LMIRcam filter wheel 3'
  SXADDPAR, hdr, 'LMIR_FW4', hdr_in.lmir_fw4, 'LMIRcam filter wheel 4'
  SXADDPAR, hdr, 'NOM_FW1', hdr_in.nom_fw1, 'NOMIC filter wheel 1'
  SXADDPAR, hdr, 'NOM_FW2', hdr_in.nom_fw2, 'NOMIC filter wheel 2'
  SXADDPAR, hdr, 'NOM_APW', hdr_in.nom_apw, 'NOMIC aperture wheel'
  SXADDPAR, hdr, 'PHA_FW1', hdr_in.pha_fw1, 'PHASECAM filter wheel 1'
  SXADDPAR, hdr, 'PHA_FW2', hdr_in.pha_fw2, 'PHASECAM filter wheel 2'
  SXADDPAR, hdr, 'PHA_IMG', hdr_in.pha_img, 'PHASECAM pupil/star imaging wheel'

  ; AO paramaters
  SXADDPAR, hdr, 'COMMENT', 'AO information', after = 'PHA_IMG'
  SXADDPAR, hdr, 'DAOMODE', hdr_in.daomode, 'Mode of the DX AO'
  SXADDPAR, hdr, 'DAOSTRHL', mean(data_in.daostrhl), 'DX strehl ratio at V-band'
  SXADDPAR, hdr, 'DCMODES', hdr_in.dcmodes, 'DX number of Zernike modes'
  SXADDPAR, hdr, 'DLOOPON', data_in[0].dloopon, 'DX AO loop on/off'
  SXADDPAR, hdr, 'DLOOPGN', hdr_in.dloopgn, 'DX AO loop gain'
  SXADDPAR, hdr, 'DWFSCFRQ', hdr_in.dwfscfrq, 'DX wavefront sensor speed [Hz]'
  SXADDPAR, hdr, 'DWFSCBIN', hdr_in.dwfscbin, 'DX wavefront sensor binning'
  SXADDPAR, hdr, 'SAOMODE', hdr_in.saomode, 'Mode of the SX AO'
  SXADDPAR, hdr, 'SAOSTRHL', mean(data_in.saostrhl), 'SX strehl ratio at V-band'
  SXADDPAR, hdr, 'SCMODES', hdr_in.scmodes, 'SX number of Zernike modes'
  SXADDPAR, hdr, 'SLOOPON', data_in[0].sloopon, 'SX AO loop on/off'
  SXADDPAR, hdr, 'SLOOPGN', hdr_in.sloopgn, 'SX AO loop gain'
  SXADDPAR, hdr, 'SWFSCFRQ', hdr_in.dwfscfrq, 'SX wavefront sensor speed [Hz]'
  SXADDPAR, hdr, 'SWFSCBIN', hdr_in.dwfscbin, 'SX wavefront sensor binning'

  ; Phasecam information
  SXADDPAR, hdr, 'COMMENT', 'PHASECam information', after = 'SWFSCBIN'
  SXADDPAR, hdr, 'PCCLOSED', mean(data_in.pcclosed), 'Mean status of the PLC (0=open, 1=close)'
  if mean(data_in.pcclosed) gt 0 then begin
    SXADDPAR, hdr, 'PCPLSP', data_in[0].pcplsp, 'Pathlength setpoint at start [degrees]'
    SXADDPAR, hdr, 'PCTIPSP', data_in[0].pctipsp, 'Tip setpoint at start [mas]'
    SXADDPAR, hdr, 'PCTLTSP', data_in[0].pctltsp, 'Tilt setpoint at start [mas]'
    if TAG_EXIST(data_in, 'fpcpistm') then begin
      SXADDPAR, hdr, 'SPCPIST', mean(data_in.spcpist), 'Mean SPC piston [um]'
      SXADDPAR, hdr, 'SPCAZ', mean(data_in.spcaz), 'Mean SPC azimuth [mas]'
      SXADDPAR, hdr, 'SPCEL', mean(data_in.spcel), 'Mean SPC elevation [mas]'
    endif else begin
      SXADDPAR, hdr, 'SPCPIST', -1, 'Mean SPC piston [um]'
      SXADDPAR, hdr, 'SPCAZ', -1, 'Mean SPC azimuth [mas]'
      SXADDPAR, hdr, 'SPCEL', -1, 'Mean SPC elevation [mas]'
    endelse
  endif

  ; Weather conditions (make sure the telescope server was connected, report -1 if not)
  SXADDPAR, hdr, 'COMMENT', 'Weather information', after = 'SPCEL'
  if n_elements(data_in.cv) gt 1 then med_cv = median(data_in.cv) else med_cv = abs(data_in.cv)
  SXADDPAR, hdr, 'CV', med_cv, 'Central value [ADU]'
  SXADDPAR, hdr, 'SEEING', hdr_in.seeing, 'Median seeing over sequence'
  SXADDPAR, hdr, 'SMTTAU', hdr_in.smttau, 'Median PWV over sequence'
  SXADDPAR, hdr, 'LBTTEMP', hdr_in.lbttemp, 'Median dome temperature over sequence [C]'
  SXADDPAR, hdr, 'WINDDIR', hdr_in.winddir, 'Median wind direction (east of north) [deg]'
  SXADDPAR, hdr, 'WINDSPD', hdr_in.windspd, 'Median wind speed [m/s]'

  ; Data reduction information
  SXADDPAR, hdr, 'COMMENT', 'Data reduction information', after = 'WINDSPD'
  SXADDPAR, hdr, 'DRS_VERS', hdr_in.drs_vers, 'DRS version for image calibration'
  SXADDPAR, hdr, 'DRS_DATE', hdr_in.drs_date, 'DRS version date for image calibration'
  SXADDPAR, hdr, 'DATE_CAL', hdr_in.date_red, 'Date of image calibration'
  SXADDPAR, hdr, 'DATE_FLX', systime(), 'Date of flux computation'
  SXADDPAR, hdr, 'BCK_MODE', fix(hdr_in.bck_mode), 'Background subtraction mode'
  SXADDPAR, hdr, 'BCK_NOD1', bck_nod1, 'Background nod ID1'
  SXADDPAR, hdr, 'BCK_NOD2', bck_nod2, 'Background nod ID2'
  SXADDPAR, hdr, 'BCK_SEL', fix(hdr_in.bck_sel), 'Background frame selection mode'
  SXADDPAR, hdr, 'FIT_MODE', fix(drs.fit_mode), 'Centroid fitting technique (0=no centroid)'
  SXADDPAR, hdr, 'IMG_MODE', fix(hdr_in.img_mode), 'Image combination mode (median, mean, mean+)'
  SXADDPAR, hdr, 'OB_MODE', fix(drs.ob_mode), 'OB mode'
  SXADDPAR, hdr, 'N_BIN', fix(drs.n_bin), 'Image re-binning size'
  SXADDPAR, hdr, 'N_CLIP', fix(drs.n_clip), 'Size of the extracted images around [XCEN,YCEN]'
  SXADDPAR, hdr, 'N_FRBCK', fix(hdr_in.n_frbck), 'Number of frames per background'
  SXADDPAR, hdr, 'NO_BPM', fix(hdr_in.no_bpm), 'Bad pixel correction off if 1 '
  SXADDPAR, hdr, 'NO_DARK', fix(hdr_in.no_dark), 'Dark subtraction off if 1 '
  SXADDPAR, hdr, 'NO_FLAT', fix(hdr_in.no_flat), 'Flat correction off if 1'
  SXADDPAR, hdr, 'PRECISIO', fix(hdr_in.precision), '0: float, 1: double'
  SXADDPAR, hdr, 'SKY_OFF', fix(drs.sky_off), '1 to turn off the sky annulus'
  SXADDPAR, hdr, 'SKY_WEI', fix(drs.sky_weight), '1 if sky pixels are weighted'

  ; Save average image (do it first because of /NO_COPY below)
  n_img = n_elements(img_in[0, 0, *])
  if n_img gt 1 then begin
    ; img_med = REFORM(MEDIAN(img_in, DIMENSION=3))
    RESISTANT_MEAN, img_in, 10, img_med, rms, num_rej, dimension = 3, /silent
    file_med = outname + '_IMG_MED.fits'
    MKHDR, hdr2, img_med ; Make header for the fits file
    WRITEFITS, file_med, img_med, hdr2
  endif

  ; Save image-cube file in the primary extension
  outfile = outname + '_IMG.fits'
  MWRFITS, img_in, outfile, hdr, /create, /silent, /no_copy

  ; Append extension with parameters that change frame by frame
  n_row = n_elements(data_in.mjd_obs)
  if keyword_set(split) then begin
    outfile = outname + '_DATA.fits'
    FXHMAKE, hdr, /init, /extend, 0
    SXADDPAR, hdr, 'COMMENT', 'Level 1 FITS file of uncalibrated data taken with LBTI/' + strtrim(hdr_in.instrum) + '.'
    SXADDPAR, hdr, 'COMMENT', 'The raw data has been reduced using the nodrs software version ' + string(hdr_in.drs_vers, format = '(F3.1)') + '.'
    SXADDPAR, hdr, 'COMMENT', 'Contact: rafael@ipac.caltech.edu.'
    SXADDPAR, hdr, 'DATE_OBS', drs.date_obs, 'date of observation'
    SXADDPAR, hdr, 'TELESCOP', 'LBT', 'telescope'
    SXADDPAR, hdr, 'INSTRUME', hdr_in.instrum, 'instrument'
    SXADDPAR, hdr, 'OBJNAME', hdr_in.objname, 'target name'
    SXADDPAR, hdr, 'DRS_VERS', hdr_in.drs_vers, 'DRS version for image calibration'
    SXADDPAR, hdr, 'DRS_DATE', hdr_in.drs_date, 'DRS version date for image calibration'
    SXADDPAR, hdr, 'DATE_CAL', hdr_in.date_red, 'Date of image calibration'
    SXADDPAR, hdr, 'DATE_FLX', systime(), 'Date of flux computation'
    SXADDPAR, hdr, 'BCK_MODE', hdr_in.bck_mode, 'Background subtraction mode'
    SXADDPAR, hdr, 'OB_MODE', drs.ob_mode, 'OB mode'
    SXADDPAR, hdr, 'FIT_MODE', drs.fit_mode, 'Centroid fitting technique (0=no centroid)'
    SXADDPAR, hdr, 'FRA_MODE', drs.fra_mode, '0: use mean-background subtracted images, 1: use raw images, 2: pca-background subtracted images'
    SXADDPAR, hdr, 'N_XPIX', hdr_in.n_xpix, 'initial number of detector rows [pix]'
    SXADDPAR, hdr, 'N_YPIX', hdr_in.n_ypix, 'initial number of detector columns [pix]'
    SXADDPAR, hdr, 'N_CLIP', drs.n_clip, 'Size of the extracted images around [XCEN,YCEN]'
    SXADDPAR, hdr, 'N_BIN', drs.n_bin, 'Image re-binning size'
    SXADDPAR, hdr, 'XCEN', mean(xcen), 'Mean initial x-coordinate on detector [pix]'
    SXADDPAR, hdr, 'YCEN', mean(ycen), 'Mean initial y-coordinate on detector [pix]'
    FXWRITE, outfile, hdr
  endif

  ; Create extension
  FXBHMAKE, hdr, n_row, /init, extver = 1, 'DATA_SERIES', 'parameters that change frame by frame'

  ; Init column number
  n_col = 13
  col = lindgen(n_col) + 1l

  ; Fill extension header with column names
  FXBADDCOL, 1l, hdr, data_in[0].file_id, 'FILE_ID', 'Raw file identification number'
  FXBADDCOL, 2l, hdr, data_in[0].mjd_obs, 'MJD_OBS', 'Modified Julian Date of observation'
  FXBADDCOL, 3l, hdr, data_in[0].lbt_utc, 'LBT_UTC', 'UTC from observatory'
  FXBADDCOL, 4l, hdr, data_in[0].lbt_lst, 'LBT_LST', 'LST from observatory'
  FXBADDCOL, 5l, hdr, data_in[0].lbt_alt, 'LBT_ALT', 'ALT from observatory'
  FXBADDCOL, 6l, hdr, data_in[0].lbt_az, 'LBT_AZ', 'AZ from observatory'
  FXBADDCOL, 7l, hdr, data_in[0].lbt_para, 'LBT_PARA', 'Parralactic angle from obs.'
  FXBADDCOL, 8l, hdr, nod_id[0], 'NOD_ID', 'Nod identification number'
  FXBADDCOL, 9l, hdr, chp_id[0], 'CHP_ID', 'Chop identification number'
  FXBADDCOL, 10l, hdr, xcen[0], 'XCEN', 'X position of the star'
  FXBADDCOL, 11l, hdr, ycen[0], 'YCEN', 'Y position of the star'
  FXBADDCOL, 12l, hdr, slope[0], 'SLOPE', 'Slope of the fitted Moffat profile'
  FXBADDCOL, 13l, hdr, flx_out.bckg_err[0], 'BCK_ERR', 'Background error'

  ; Write extension header to FITS file
  FXBCREATE, unit, outfile, hdr
  FXBWRITM, unit, col, data_in.file_id, data_in.mjd_obs, data_in.lbt_utc, data_in.lbt_lst, data_in.lbt_alt, data_in.lbt_az, data_in.lbt_para, $
    nod_id, chp_id, xcen, ycen, slope, transpose(reform(flx_out.bckg_err[*, 0]))
  FXBFINISH, unit

  ; Copy files if FILE_ID is an array
  n_id = n_elements(file_id)
  for i_id = 1, n_id - 1 do begin
    id_string = '_ID' + string(file_id[i_id], format = '(I03)') + '_'
    outname = strcompress(sav_path + 'UT' + strtrim(drs.date_obs, 2) + id_string + strtrim(hdr_in.flag, 2) + '_' + strtrim(hdr_in.objname, 2) + '_DIT-' + string(1d+3 * hdr_in.int_time, format = '(I0)') + 'ms_' + $
      string(1d+6 * hdr_in.lam_cen, format = '(I0)') + 'um_' + tag, /remove_all)

    ; Save image-cube file in the primary extension
    outfile2 = outname + '_IMG.fits'
    file_copy, outfile, outfile2, /overwrite
    if n_img gt 1 then begin
      file_med2 = outname + '_IMG_MED.fits'
      file_copy, file_med, file_med2, /overwrite
    endif
  endfor
end