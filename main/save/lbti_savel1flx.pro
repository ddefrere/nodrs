;+
; NAME: LBTI_SAVEL1FLX
;
; PURPOSE:
;   This procedure saves data tables in level 1 FITS files.
;
; INPUTS:
;   flx_in        :  On input, the measured fluxes
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
;   Version 1.1, 25-MAY-2014, DD: removed keyword OB_MODE
;   Version 1.3, 25-MAY-2013, DD: enabled FILE_ID as an array
;   Version 1.4, 26-MAY-2013, DD: minor bug corrections
;   Version 1.5, 29-MAY-2013, DD: added background bias
;   Version 1.6, 30-MAY-2013, DD: transposed output arrays
;   Version 1.7, 22-SEP-2013, DD: updated for new formalism and removed keyword TAG
;   Version 1.8, 28-OCT-2013, DD: added project ID to output header
;   Version 1.9, 15-JAN-2015, DD: added NO_FLAT, NO_DARK, and NO_BPM to output header
;   Version 2.0, 16-FEB-2015, DD: Added new PHASECam fieds (required by NSC)
;   Version 2.1, 04-APR-2015, DD: Added nodding period
;   Version 2.2, 04-MAY-2015, DD: Added background nod
;   Version 2.3  24-AUG-2015, DD: Added setpoint information to output header
;   Version 2.4, 06-OCT-2015, DD: Added PCWAV and PCSPEC to output header
;   Version 2.5, 03-DEC-2015, DD: Now save the individual background radii separately
;   Version 2.6, 07-DEC-2015, DD: Now save BCK_IRAD, and BCK_ORAD seperately for each APER_RAD
;   Version 2.7, 07-FEB-2016, DD: Added PHASECam stats for beam2
;   Version 2.8, 02-JUL-2016, DD: Added central value
;   Version 2.9, 25-OCT-2016, DD: Added slope RMS
;   Version 3.0, 30-JAN-2017, DD: Adapted to new format
;   Version 3.1, 12-FEB-2017, DD: Added concatenation
;   Version 3.1, 19-FEB-2017, DD: Now works if restored data have only one frame
;   Version 3.2, 04-APR-2017, DD: Added more target information to output header
;   Version 3.3, 13-MAY-2017, DD: Added NSKYPIX (number of pixels in background region) to output files
;   Version 3.4, 28-JUL-2017, DD: Implemented background floor mode
;   Version 3.5, 10-AUG-2020, DD: Added SKY_COL to output header
;   Version 3.6, 15-OCT-2023, DD: Updated for FRA_MODE=2 (i.e., PCA background subtraction)
;   Version 3.7, 24-MAY-2024, DD: Update file permission
;   Version 3.8, 02-SEP-2024, DD: Added parameter SKY_OFF
;
pro LBTI_SAVEL1FLX, flx_in, hdr_in, data_in, outfile = outfile, file_id = file_id
  compile_opt idl2

  ; Define operational parameters
  common GLOBAL, prm, cnf, wav, tgt, pth, drs, log

  ; Keyword sanity check
  if not keyword_set(file_id) then file_id = [0]

  ; Derive which beam is this
  ; If beam ID is 0, rely on OBSTYPE to define TAG
  ; If beam ID is 1, tag cannot be NULL (must be PHOT2 or BCKG)
  tag = 'NULL' ; NULL by default unless below
  if hdr_in.beam_id eq 0 then begin
    xcen = data_in.xcen_sx
    ycen = data_in.ycen_sx
    fwhm = 0.5 * ((data_in.fwhmx_sx) + (data_in.fwhmy_sx))
    slrms = data_in.ssloprms
    slope = data_in.slope_sx
    if abs(hdr_in.obstype) eq 0 then tag = 'PHOT1'
    if abs(hdr_in.obstype) eq 3 then tag = 'BCKG' ; + STRING(FIX(hdr_in.bck_nod), FORMAT='(I03)')
  endif else begin
    xcen = data_in.xcen_dx
    ycen = data_in.ycen_dx
    fwhm = 0.5 * ((data_in.fwhmx_dx) + (data_in.fwhmy_dx))
    slrms = data_in.dsloprms
    slope = data_in.slope_dx
    if abs(hdr_in.obstype) eq 0 then tag = 'PHOT2' else tag = 'BCKG' ; + STRING(FIX(hdr_in.bck_nod), FORMAT='(I03)')
  endelse

  ; Number of elements and aperture radii
  n_aper = n_elements(drs.aper_rad)
  nod_id = replicate(hdr_in.nod_id, n_elements(data_in.xcen_dx))
  chp_id = replicate(1, n_elements(data_in.xcen_dx)) ; not used currently

  ; Prepare data
  raw_id = long(data_in.file_id)
  mjd_obs = double(data_in.mjd_obs)
  bckg_meas = float(flx_in.bckg_meas)
  bckg_err = float(flx_in.bckg_err)
  flx_tot = double(flx_in.flx_tot)
  flx_err = double(flx_in.flx_err)
  bckg_meas2 = float(flx_in.bckg_meas2)
  bckg_err2 = float(flx_in.bckg_err2)
  flx_tot2 = double(flx_in.flx_tot2)
  flx_err2 = double(flx_in.flx_err2)
  if TAG_EXIST(data_in, 'fpcpistm') then begin
    pcplsp = float(data_in.pcplsp)
    pctipsp = float(data_in.pctipsp)
    pctltsp = float(data_in.pctltsp)
    spdthpos = float(data_in.spdthpos)
    fpcpistm = float(data_in.fpcpistm)
    fpcpists = float(data_in.fpcpists)
    fpcazm = float(data_in.fpcazm)
    fpcelm = float(data_in.fpcelm)
    fpcazs = float(data_in.fpcazs)
    fpcels = float(data_in.fpcels)
    pcphmean = float(data_in.pcphmean)
    pcphstd = float(data_in.pcphstd)
    pcphmcos = float(data_in.pcphmcos)
    pcphmsin = float(data_in.pcphmsin)
    pcfjmps = fix(data_in.pcfjmps)
    pcmsnr = float(data_in.pcmsnr)
    pcphmean2 = float(data_in.pcphmean2)
    pcphstd2 = float(data_in.pcphstd2)
    pcphmcos2 = float(data_in.pcphmcos2)
    pcphmsin2 = float(data_in.pcphmsin2)
    pcmsnr2 = float(data_in.pcmsnr2)
  endif

  ; Define output file name
  if not keyword_set(outfile) then begin
    ; Create directory
    sav_path = pth.l1Fits_path + drs.date_obs + drs.dir_label + pth.sep
    if not file_test(sav_path) then file_mkdir, sav_path
    spawn, 'chmod -R 775 ' + sav_path
    ; Create file name
    id_string = '_ID' + string(file_id[0], format = '(I03)') + '_'
    outfile = strcompress(sav_path + 'UT' + strtrim(drs.date_obs, 2) + id_string + strtrim(hdr_in.flag, 2) + '_' + strtrim(hdr_in.objname, 2) + '_DIT-' + string(1d+3 * hdr_in.int_time, format = '(I0)') + 'ms_' + $
      string(1d+6 * hdr_in.lam_cen, format = '(I0)') + 'um_' + tag + '.fits', /remove_all)
  endif

  ; If file already exists, open and decide what to do
  if file_test(outfile) and tag eq 'BCKG' then begin
    hdr_res = HEADFITS(outfile, /silent)
    bck_nod1 = fix(FXPAR(hdr_res, 'BCK_NOD1', /nocontinue))
    bck_nod2 = fix(FXPAR(hdr_res, 'BCK_NOD2', /nocontinue))
    ; If different background nod, append. Otherwise, erase.
    if bck_nod1 ne hdr_in.nod_id and bck_nod2 ne hdr_in.nod_id then begin
      bck_nod2 = hdr_in.nod_id
      data_res = MRDFITS(outfile, 1, /silent)
      raw_id = [data_res.file_id, raw_id]
      mjd_obs = [data_res.mjd_obs, mjd_obs]
      nod_id = [data_res.nod_id, nod_id]
      chp_id = [data_res.chp_id, chp_id]
      xcen = [data_res.xcen, xcen]
      ycen = [data_res.ycen, ycen]
      fwhm = [data_res.fwhm, fwhm]
      slrms = [data_res.slrms, slrms]
      if (size(bckg_meas))[0] gt 1 then begin
        bckg_meas = transpose([[data_res.bck_tot], [transpose(bckg_meas)]])
        bckg_err = transpose([[data_res.bck_err], [transpose(bckg_err)]])
        flx_tot = transpose([[data_res.flx_tot], [transpose(flx_tot)]])
        flx_err = transpose([[data_res.flx_err], [transpose(flx_err)]])
        bckg_meas2 = transpose([[data_res.bck_tot2], [transpose(bckg_meas2)]])
        bckg_err2 = transpose([[data_res.bck_err2], [transpose(bckg_err2)]])
        flx_tot2 = transpose([[data_res.flx_tot2], [transpose(flx_tot2)]])
        flx_err2 = transpose([[data_res.flx_err2], [transpose(flx_err2)]])
      endif else begin
        bckg_meas = [data_res.bck_tot, bckg_meas]
        bckg_err = [data_res.bck_err, bckg_err]
        flx_tot = [data_res.flx_tot, flx_tot]
        flx_err = [data_res.flx_err, flx_err]
        bckg_meas2 = [data_res.bck_tot2, bckg_meas2]
        bckg_err2 = [data_res.bck_err2, bckg_err2]
        flx_tot2 = [data_res.flx_tot2, flx_tot2]
        flx_err2 = [data_res.flx_err2, flx_err2]
      endelse
      if TAG_EXIST(data_in, 'fpcpistm') then begin
        pcplsp = [data_res.pcplsp, pcplsp]
        pctipsp = [data_res.pctipsp, pctipsp]
        pctltsp = [data_res.pctltsp, pctltsp]
        spdthpos = [data_res.spdthpos, spdthpos]
        fpcpistm = [data_res.fpcpistm, fpcpistm]
        fpcpists = [data_res.fpcpists, fpcpists]
        fpcazm = [data_res.fpcazm, fpcazm]
        fpcelm = [data_res.fpcelm, fpcelm]
        fpcazs = [data_res.fpcazs, fpcazs]
        fpcels = [data_res.fpcels, fpcels]
        pcphmean = [data_res.pcphmean, pcphmean]
        pcphstd = [data_res.pcphstd, pcphstd]
        pcphmcos = [data_res.pcphmcos, pcphmcos]
        pcphmsin = [data_res.pcphmsin, pcphmsin]
        pcfjmps = [data_res.pcfjmps, pcfjmps]
        pcmsnr = [data_res.pcmsnr, pcmsnr]
        pcphmean2 = [data_res.pcphmean2, pcphmean2]
        pcphstd2 = [data_res.pcphstd2, pcphstd2]
        pcphmcos2 = [data_res.pcphmcos2, pcphmcos2]
        pcphmsin2 = [data_res.pcphmsin2, pcphmsin2]
        pcmsnr2 = [data_res.pcmsnr2, pcmsnr2]
      endif
    endif
  endif else begin
    if tag eq 'BCKG' then bck_nod1 = hdr_in.nod_id else bck_nod1 = hdr_in.bck_nod
    bck_nod2 = 0
  endelse

  ; Add comment to the header
  FXHMAKE, hdr, /init, /extend, 0
  SXADDPAR, hdr, 'COMMENT', 'Level 1 FITS file of uncalibrated data taken with LBTI/' + strtrim(hdr_in.instrum) + '.'
  SXADDPAR, hdr, 'COMMENT', 'The raw data has been reduced using the nodrs software version ' + string(hdr_in.drs_vers, format = '(F3.1)') + '.'
  SXADDPAR, hdr, 'COMMENT', 'Contact: lbtisupp@ipac.caltech.edu.'

  ; Add keyword definitions to header
  SXADDPAR, hdr, 'COMMENT', 'Target and observation parameters', before = 'DATE_OBS'
  SXADDPAR, hdr, 'DATE_OBS', drs.date_obs, 'date of observation'
  SXADDPAR, hdr, 'TELESCOP', 'LBT', 'telescope'
  SXADDPAR, hdr, 'INSTRUME', hdr_in.instrum, 'instrument'
  SXADDPAR, hdr, 'PID', hdr_in.pid, 'Project ID'
  SXADDPAR, hdr, 'OBJNAME', hdr_in.objname, 'target name'
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

  ; AO keywords
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
  if TAG_EXIST(data_in, 'fpcpistm') then begin
    SXADDPAR, hdr, 'SPCPIST', mean(data_in.spcpist), 'Mean SPC piston [um]'
    SXADDPAR, hdr, 'SPCAZ', mean(data_in.spcaz), 'Mean SPC azimuth [mas]'
    SXADDPAR, hdr, 'SPCEL', mean(data_in.spcel), 'Mean SPC elevation [mas]'
    SXADDPAR, hdr, 'SPCAZ', mean(data_in.spcaz), 'Mean SPC azimuth [mas]'
  endif else begin
    SXADDPAR, hdr, 'SPCPIST', -1, 'Mean SPC piston [um]'
    SXADDPAR, hdr, 'SPCAZ', -1, 'Mean SPC azimuth [mas]'
    SXADDPAR, hdr, 'SPCEL', -1, 'Mean SPC elevation [mas]'
    SXADDPAR, hdr, 'SPCAZ', -1, 'Mean SPC azimuth [mas]'
  endelse

  ; Weather conditions
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
  SXADDPAR, hdr, 'BFL_MODE', fix(drs.bfl_mode), 'Background floor mode (0: mean, 1: median)'
  SXADDPAR, hdr, 'FLX_MODE', fix(drs.flx_mode), '0=aperture phot, 1=PSF fitting'
  SXADDPAR, hdr, 'FIT_MODE', fix(drs.fit_mode), 'Centroid fitting technique (0=no centroid)'
  SXADDPAR, hdr, 'FRA_MODE', fix(drs.fra_mode), '0: use mean-background subtracted images, 1: use raw images, 2: pca-background subtracted images'
  SXADDPAR, hdr, 'IMG_MODE', fix(hdr_in.img_mode), 'Image combination mode'
  SXADDPAR, hdr, 'OB_MODE', fix(drs.ob_mode), 'OB mode'
  if drs.flx_mode lt 2 then begin
    for i_a = 0, n_aper - 1 do begin
      i_size = string(i_a, format = '(I0)')
      SXADDPAR, hdr, 'APERRAD' + i_size, fix(flx_in.aper_rad[i_a]), 'Radius ' + i_size + ' for aperture photometry [pix]'
      SXADDPAR, hdr, 'BCKIRAD' + i_size, fix(flx_in.bck_irad[i_a]), 'Inner radius ' + i_size + ' for background computation [pix]'
      SXADDPAR, hdr, 'BCKORAD' + i_size, fix(flx_in.bck_orad[i_a]), 'Outer radius ' + i_size + ' for background computation [pix]'
      SXADDPAR, hdr, 'NSKYPIX' + i_size, fix(flx_in.nsky[i_a]), 'Number of pixels in background region'
    endfor
  endif
  SXADDPAR, hdr, 'BCK_NOD1', bck_nod1, 'Background nod ID1'
  SXADDPAR, hdr, 'BCK_NOD2', bck_nod2, 'Background nod ID2'
  SXADDPAR, hdr, 'BCK_SEL', fix(hdr_in.bck_sel), 'Background frame selection mode'
  SXADDPAR, hdr, 'N_BIN', fix(drs.n_bin), 'Image re-binning size'
  SXADDPAR, hdr, 'N_CLIP', fix(drs.n_clip), 'Size of the extracted images around [XCEN,YCEN]'
  SXADDPAR, hdr, 'N_FRBCK', fix(hdr_in.n_frbck), 'Number of frames per background'
  SXADDPAR, hdr, 'N_REJ', fix(hdr_in.n_rej), 'Number of rejected frames'
  SXADDPAR, hdr, 'NO_BPM', fix(hdr_in.no_bpm), 'Bad pixel correction off if 1 '
  SXADDPAR, hdr, 'NO_DARK', fix(hdr_in.no_dark), 'Dark subtraction off if 1 '
  SXADDPAR, hdr, 'NO_FLAT', fix(hdr_in.no_flat), 'Flat correction off if 1'
  SXADDPAR, hdr, 'PRECISIO', fix(hdr_in.precision), '0: float, 1: double'
  SXADDPAR, hdr, 'SKY_COL', fix(drs.sky_col), '1 to use only pixels in the same column as photometric aperture'
  SXADDPAR, hdr, 'SKY_OFF', fix(drs.sky_off), '1 to turn off the sky annulus'
  SXADDPAR, hdr, 'SKY_WEI', fix(drs.sky_weight), '1 if sky pixels are weighted'
  FXWRITE, outfile, hdr

  ; Create header for the main table
  n_row = n_elements(raw_id)
  FXBHMAKE, hdr, n_row, /init, extver = 1, 'RESULTS_SUMMARY', 'results of the data reduction'

  ; Init column number
  if TAG_EXIST(data_in, 'fpcpistm') then n_col = 37 else n_col = 16
  col = lindgen(n_col) + 1l

  ; Fill extension header with column names
  FXBADDCOL, 1l, hdr, 0l, 'FILE_ID', 'Raw file identification number'
  FXBADDCOL, 2l, hdr, 0d, 'MJD_OBS', 'Modified Julian Date of observation'
  FXBADDCOL, 3l, hdr, 0s, 'NOD_ID', 'Nod identification number'
  FXBADDCOL, 4l, hdr, 0s, 'CHP_ID', 'Chop identification number'
  FXBADDCOL, 5l, hdr, 0., 'XCEN', 'Star initial X position', tunit = 'pix'
  FXBADDCOL, 6l, hdr, 0., 'YCEN', 'Star initial Y position', tunit = 'pix'
  FXBADDCOL, 7l, hdr, 0., 'FWHM', 'Fitted FWHM', tunit = 'pix'
  FXBADDCOL, 8l, hdr, 0., 'SLRMS', 'AO slope RMS'
  FXBADDCOL, 9l, hdr, fltarr(n_aper), 'BCK_TOT', 'Mean background per pixel', tunit = 'ADU'
  FXBADDCOL, 10l, hdr, fltarr(n_aper), 'BCK_ERR', 'RMS background per pixel', tunit = 'ADU'
  FXBADDCOL, 11l, hdr, dblarr(n_aper), 'FLX_TOT', 'Measured flux', tunit = 'ADU'
  FXBADDCOL, 12l, hdr, dblarr(n_aper), 'FLX_ERR', 'Error on the measured flux', tunit = 'ADU'
  FXBADDCOL, 13l, hdr, fltarr(n_aper), 'BCK_TOT2', 'Empty-region mean background per pixel', tunit = 'ADU'
  FXBADDCOL, 14l, hdr, fltarr(n_aper), 'BCK_ERR2', 'Empty-region RMS background per pixel', tunit = 'ADU'
  FXBADDCOL, 15l, hdr, dblarr(n_aper), 'FLX_TOT2', 'Empty-region measured flux', tunit = 'ADU'
  FXBADDCOL, 16l, hdr, dblarr(n_aper), 'FLX_ERR2', 'Empty-region error on the measured flux', tunit = 'ADU'
  if TAG_EXIST(data_in, 'fpcpistm') then begin
    FXBADDCOL, 17l, hdr, 0., 'PCPLSP', 'Pathlength setpoint', tunit = 'deg'
    FXBADDCOL, 18l, hdr, 0., 'PCTIPSP', 'Tip setpoint', tunit = 'mas'
    FXBADDCOL, 19l, hdr, 0., 'PCTLTSP', 'Tilt setpoint', tunit = 'mas'
    FXBADDCOL, 20l, hdr, 0., 'SPDTHPOS', 'Tilt setpoint', tunit = 'mas'
    FXBADDCOL, 21l, hdr, 0., 'FPCPISTM', 'FPC average piston', tunit = 'um'
    FXBADDCOL, 22l, hdr, 0., 'FPCPISTS', 'FPC RMS piston', tunit = 'um'
    FXBADDCOL, 23l, hdr, 0., 'FPCAZM', 'FPC average azimuth', tunit = 'arcsec'
    FXBADDCOL, 24l, hdr, 0., 'FPCELM', 'FPC average elevation', tunit = 'arcsec'
    FXBADDCOL, 25l, hdr, 0., 'FPCAZS', 'FPC RMS azimuth', tunit = 'arcsec'
    FXBADDCOL, 26l, hdr, 0., 'FPCELS', 'FPC RMS elevation', tunit = 'arcsec'
    FXBADDCOL, 27l, hdr, 0., 'PCPHMEAN', 'Average measured piston', tunit = 'um'
    FXBADDCOL, 28l, hdr, 0., 'PCPHSTD', 'RMS measured piston', tunit = 'um'
    FXBADDCOL, 29l, hdr, 0., 'PCPHMCOS', 'Mean cosinus(phi) over the last X ms'
    FXBADDCOL, 30l, hdr, 0., 'PCPHMSIN', 'Mean sinus(phi) over the last X ms'
    FXBADDCOL, 31l, hdr, 0s, 'PCFJMPS', 'Number of CG jumps over the last X ms'
    FXBADDCOL, 32l, hdr, 0., 'PCMSNR', 'K-band FFT-peak SNR'
    FXBADDCOL, 33l, hdr, 0., 'PCPHMEAN2', 'Average measured piston', tunit = 'um'
    FXBADDCOL, 34l, hdr, 0., 'PCPHSTD2', 'RMS measured piston', tunit = 'um'
    FXBADDCOL, 35l, hdr, 0., 'PCPHMCOS2', 'Mean cosinus(phi) over the last X ms'
    FXBADDCOL, 36l, hdr, 0., 'PCPHMSIN2', 'Mean sinus(phi) over the last X ms'
    FXBADDCOL, 37l, hdr, 0., 'PCMSNR2', 'K-band FFT-peak SNR'
  endif

  ; Write extension header to FITS file
  FXBCREATE, unit, outfile, hdr
  if TAG_EXIST(data_in, 'fpcpistm') $
    then FXBWRITM, unit, col, raw_id, mjd_obs, fix(nod_id), fix(chp_id), float(xcen), float(ycen), float(fwhm), float(slrms), $
      transpose(bckg_meas), transpose(bckg_err), transpose(flx_tot), transpose(flx_err), transpose(bckg_meas2), transpose(bckg_err2), $
      transpose(flx_tot2), transpose(flx_err2), pcplsp, pctipsp, pctltsp, spdthpos, fpcpistm, $
      fpcpists, fpcazm, fpcelm, fpcazs, fpcels, pcphmean, pcphstd, pcphmcos, pcphmsin, fix(pcfjmps), pcmsnr, pcphmean2, pcphstd2, $
      pcphmcos2, pcphmsin2, pcmsnr2 $
    else FXBWRITM, unit, col, raw_id, mjd_obs, fix(nod_id), fix(chp_id), float(xcen), float(ycen), float(fwhm), float(slrms), $
      transpose(bckg_meas), transpose(bckg_err), transpose(flx_tot), transpose(flx_err), transpose(bckg_meas2), transpose(bckg_err2), $
      transpose(flx_tot2), transpose(flx_err2)
  FXBFINISH, unit

  ; Copy files if FILE_ID is an array
  n_id = n_elements(file_id)
  for i_id = 1, n_id - 1 do begin
    id_string = '_ID' + string(file_id[i_id], format = '(I03)') + '_'
    outfile2 = strcompress(sav_path + 'UT' + strtrim(drs.date_obs, 2) + id_string + strtrim(hdr_in.flag, 2) + '_' + strtrim(hdr_in.objname, 2) + '_DIT-' + string(1d+3 * hdr_in.int_time, format = '(I0)') + 'ms_' + $
      string(1d+6 * hdr_in.lam_cen, format = '(I0)') + 'um_' + tag + '.fits', /remove_all)
    file_copy, outfile, outfile2, /overwrite
  endfor
end