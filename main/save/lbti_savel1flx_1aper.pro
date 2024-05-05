;+
; NAME: LBTI_SAVEL1FLX_1APER
; 
; PURPOSE:
;   This procedure saves data tables in level 1 FITS files for only one aperture.
;
; INPUTS:
;   data_in       :  On input, the measured fluxes
;   header        :  On input, the header of the fits file
;   aper_id       :  Integer giving the aperture column to be saved
;
; KEYWORDS:
;   OUTFILE       :  Name of the output file (superseed the automatic file name)
;   TAG           :  Optional TAG name to give the ttype of file (i.e., NULL, PHOT1, PHOT2, BCKG)
;
; MODIFICATION HISTORY:
;   Version 1.0, 05-MAY-2024, by Denis Defrère, University of Arizona, denis@lbti.com

PRO LBTI_SAVEL1FLX_1APER, data, header, aper_id, OUTFILE=outfile, TAG=tag

; Define operational parameters
COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs, log

; Read some header information
objname  = STRCOMPRESS(STRTRIM(FXPAR(header, 'OBJNAME', /NOCONTINUE)), /REMOVE_ALL)
obstype  = FXPAR(header, 'OBSTYPE', /NOCONTINUE)
exptime  = FXPAR(header, 'EXPTIME', /NOCONTINUE)
lam_cen  = FXPAR(header, 'WAVELENG', /NOCONTINUE)
flag     = FXPAR(header, 'FLAG', /NOCONTINUE)

; Get data type
IF NOT KEYWORD_SET(TAG) THEN BEGIN
    IF ABS(obstype) EQ 0 THEN tag = 'PHOT1'
    IF ABS(obstype) EQ 2 THEN tag = 'NULL' 
    IF ABS(obstype) EQ 3 THEN tag = 'BCKG' 
ENDIF

; Prepare data
raw_id     = LONG(data.file_id)
mjd_obs    = DOUBLE(data.mjd_obs)
bckg_meas  = FLOAT(data.bck_tot[aper_id])
bckg_err   = FLOAT(data.bck_err[aper_id])
flx_tot    = DOUBLE(data.flx_tot[aper_id])
flx_err    = DOUBLE(data.flx_err[aper_id])
bckg_meas2 = FLOAT(data.bck_tot2[aper_id])
bckg_err2  = FLOAT(data.bck_err2[aper_id])
flx_tot2   = DOUBLE(data.flx_tot2[aper_id])
flx_err2   = DOUBLE(data.flx_err2[aper_id])
IF TAG_EXIST(data, 'fpcpistm') THEN BEGIN
  pcplsp    = FLOAT(data_in.pcplsp)
  pctipsp   = FLOAT(data_in.pctipsp)
  pctltsp   = FLOAT(data_in.pctltsp)
  spdthpos  = FLOAT(data_in.spdthpos)
  fpcpistm  = FLOAT(data_in.fpcpistm)
  fpcpists  = FLOAT(data_in.fpcpists)
  fpcazm    = FLOAT(data_in.fpcazm)
  fpcelm    = FLOAT(data_in.fpcelm)
  fpcazs    = FLOAT(data_in.fpcazs)
  fpcels    = FLOAT(data_in.fpcels)
  pcphmean  = FLOAT(data_in.pcphmean)
  pcphstd   = FLOAT(data_in.pcphstd)
  pcphmcos  = FLOAT(data_in.pcphmcos)
  pcphmsin  = FLOAT(data_in.pcphmsin)
  pcfjmps   = FIX(data_in.pcfjmps)
  pcmsnr    = FLOAT(data_in.pcmsnr)
  pcphmean2 = FLOAT(data_in.pcphmean2)
  pcphstd2  = FLOAT(data_in.pcphstd2)
  pcphmcos2 = FLOAT(data_in.pcphmcos2)
  pcphmsin2 = FLOAT(data_in.pcphmsin2)
  pcmsnr2   = FLOAT(data_in.pcmsnr2) 
ENDIF  

; Define output file name
IF NOT KEYWORD_SET(OUTFILE) THEN BEGIN
  ; Create directory
  sav_path = pth.l1fits_path + drs.date_obs + drs.dir_label + pth.sep + 'filtered'
  IF NOT FILE_TEST(sav_path) THEN FILE_MKDIR, sav_path
  ; Create file name
  id_string = '_ID' + STRING(raw_id[0], FORMAT='(I03)') + '_'
  outfile   = STRCOMPRESS(sav_path + 'UT' + STRTRIM(drs.date_obs, 2) + id_string + STRTRIM(flag, 2) + '_' + STRTRIM(objname, 2) + '_DIT-' + STRING(1D+3*exptime, FORMAT='(I0)') + 'ms_' + $
              STRING(1D+6*lam_cen, FORMAT='(I0)') + 'um_' + '_APER-' + STRING(aper_rad, FORMAT='(I0)') + '_FILT-' + tag + '.fits' , /REMOVE_ALL)
ENDIF

; Now add the right aperture
aper_rad = FXPAR(header, 'APERRAD' + STRING(aper_id, FORMAT='(I0)'), /NOCONTINUE)
bck_irad = FXPAR(header, 'BCKIRAD' + STRING(aper_id, FORMAT='(I0)'), /NOCONTINUE)
bck_orad = FXPAR(header, 'BCKORAD' + STRING(aper_id, FORMAT='(I0)'), /NOCONTINUE)
nskypix  = FXPAR(header, 'NSKYPIX' + STRING(aper_id, FORMAT='(I0)'), /NOCONTINUE)
SXADDPAR, hdr, 'APERRAD', FIX(aper_rad), 'Radius for aperture photometry [pix]'
SXADDPAR, hdr, 'BCKIRAD', FIX(bck_irad), 'Inner radius for background computation [pix]'
SXADDPAR, hdr, 'BCKORAD', FIX(bck_orad), 'Outer radius for background computation [pix]'
SXADDPAR, hdr, 'NSKYPIX', FIX(nskypix),  'Number of pixels in background region'

; REMOVE other aperture info from the header
!err = 0
i_aper = 0
WHILE !err NE -1 DO BEGIN
    SXDELPAR, header, 'APERRAD' + STRING(i_aper, FORMAT='(I0)')
    SXDELPAR, header, 'BCKIRAD' + STRING(i_aper, FORMAT='(I0)')
    SXDELPAR, header, 'BCKORAD' + STRING(i_aper, FORMAT='(I0)')
    SXDELPAR, header, 'NSKYPIX' + STRING(i_aper, FORMAT='(I0)')
    i_aper = i_aper + 1 
    tmp = FXPAR(header, 'APERRAD' + STRING(i_aper, FORMAT='(I0)'), /NOCONTINUE)
ENDWHILE

; Add comment to the header
WRITEFITS, outfile, header

; Create header for the main table
n_row = N_ELEMENTS(raw_id)
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'RESULTS_SUMMARY', 'results of the data reduction'

; Init column number
IF TAG_EXIST(data_in, 'fpcpistm') THEN n_col = 37 ELSE n_col = 16
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, 0L, 'FILE_ID',  'Raw file identification number' 
FXBADDCOL, 2L, hdr, 0D, 'MJD_OBS',  'Modified Julian Date of observation'
FXBADDCOL, 3L, hdr, 0,  'NOD_ID',   'Nod identification number'
FXBADDCOL, 4L, hdr, 0,  'CHP_ID',   'Chop identification number'
FXBADDCOL, 5L, hdr, 0.,  'XCEN',     'Star initial X position', TUNIT='pix'
FXBADDCOL, 6L, hdr, 0.,  'YCEN',     'Star initial Y position', TUNIT='pix'
FXBADDCOL, 7L, hdr, 0.,  'FWHM',     'Fitted FWHM', TUNIT='pix'
FXBADDCOL, 8L, hdr, 0.,  'SLRMS',    'AO slope RMS'
FXBADDCOL, 9L, hdr, 0.,  'BCK_TOT',  'Mean background per pixel', TUNIT='ADU'
FXBADDCOL, 10L, hdr, 0., 'BCK_ERR',  'RMS background per pixel', TUNIT='ADU'
FXBADDCOL, 11L, hdr, 0D, 'FLX_TOT',  'Measured flux', TUNIT='ADU'
FXBADDCOL, 12L, hdr, 0D, 'FLX_ERR',  'Error on the measured flux', TUNIT='ADU'
FXBADDCOL, 13L, hdr, 0., 'BCK_TOT2', 'Empty-region mean background per pixel', TUNIT='ADU'
FXBADDCOL, 14L, hdr, 0., 'BCK_ERR2', 'Empty-region RMS background per pixel', TUNIT='ADU'
FXBADDCOL, 15L, hdr, 0D, 'FLX_TOT2', 'Empty-region measured flux', TUNIT='ADU'
FXBADDCOL, 16L, hdr, 0D, 'FLX_ERR2', 'Empty-region error on the measured flux', TUNIT='ADU'
IF TAG_EXIST(data_in, 'fpcpistm') THEN BEGIN
  FXBADDCOL, 17L, hdr, 0., 'PCPLSP',   'Pathlength setpoint', TUNIT='deg'
  FXBADDCOL, 18L, hdr, 0., 'PCTIPSP',  'Tip setpoint', TUNIT='mas'
  FXBADDCOL, 19L, hdr, 0., 'PCTLTSP',  'Tilt setpoint', TUNIT='mas' 
  FXBADDCOL, 20L, hdr, 0., 'SPDTHPOS', 'Tilt setpoint', TUNIT='mas'  
  FXBADDCOL, 21L, hdr, 0., 'FPCPISTM', 'FPC average piston', TUNIT='um'
  FXBADDCOL, 22L, hdr, 0., 'FPCPISTS', 'FPC RMS piston', TUNIT='um' 
  FXBADDCOL, 23L, hdr, 0., 'FPCAZM',   'FPC average azimuth', TUNIT='arcsec'
  FXBADDCOL, 24L, hdr, 0., 'FPCELM',   'FPC average elevation', TUNIT='arcsec'
  FXBADDCOL, 25L, hdr, 0., 'FPCAZS',   'FPC RMS azimuth', TUNIT='arcsec'
  FXBADDCOL, 26L, hdr, 0., 'FPCELS',   'FPC RMS elevation', TUNIT='arcsec'
  FXBADDCOL, 27L, hdr, 0., 'PCPHMEAN', 'Average measured piston', TUNIT='um' 
  FXBADDCOL, 28L, hdr, 0., 'PCPHSTD',  'RMS measured piston', TUNIT='um' 
  FXBADDCOL, 29L, hdr, 0., 'PCPHMCOS', 'Mean cosinus(phi) over the last X ms'
  FXBADDCOL, 30L, hdr, 0., 'PCPHMSIN', 'Mean sinus(phi) over the last X ms'
  FXBADDCOL, 31L, hdr, 0,  'PCFJMPS',  'Number of CG jumps over the last X ms'  
  FXBADDCOL, 32L, hdr, 0., 'PCMSNR',   'K-band FFT-peak SNR' 
  FXBADDCOL, 33L, hdr, 0., 'PCPHMEAN2','Average measured piston', TUNIT='um'
  FXBADDCOL, 34L, hdr, 0., 'PCPHSTD2', 'RMS measured piston', TUNIT='um'
  FXBADDCOL, 35L, hdr, 0., 'PCPHMCOS2','Mean cosinus(phi) over the last X ms'
  FXBADDCOL, 36L, hdr, 0., 'PCPHMSIN2','Mean sinus(phi) over the last X ms'
  FXBADDCOL, 37L, hdr, 0., 'PCMSNR2',  'K-band FFT-peak SNR'
ENDIF

; Write extension header to FITS file
FXBCREATE, unit, outfile, hdr
IF TAG_EXIST(data_in, 'fpcpistm') $
 THEN FXBWRITM,  unit, col, raw_id, mjd_obs, FIX(nod_id), FIX(chp_id), FLOAT(xcen), FLOAT(ycen), FLOAT(fwhm), FLOAT(slrms), $
                 TRANSPOSE(bckg_meas), TRANSPOSE(bckg_err), TRANSPOSE(flx_tot), TRANSPOSE(flx_err), TRANSPOSE(bckg_meas2), TRANSPOSE(bckg_err2), $
                 TRANSPOSE(flx_tot2), TRANSPOSE(flx_err2), pcplsp, pctipsp, pctltsp, spdthpos, fpcpistm, $
                 fpcpists, fpcazm, fpcelm, fpcazs, fpcels, pcphmean, pcphstd, pcphmcos, pcphmsin, pcfjmps, pcmsnr, pcphmean2, pcphstd2, $
                 pcphmcos2, pcphmsin2, pcmsnr2 $
 ELSE FXBWRITM,  unit, col, raw_id, mjd_obs, FIX(nod_id), FIX(chp_id), FLOAT(xcen), FLOAT(ycen), FLOAT(fwhm), FLOAT(slrms), $
                 TRANSPOSE(bckg_meas), TRANSPOSE(bckg_err), TRANSPOSE(flx_tot), TRANSPOSE(flx_err), TRANSPOSE(bckg_meas2), TRANSPOSE(bckg_err2), $
                 TRANSPOSE(flx_tot2), TRANSPOSE(flx_err2)
FXBFINISH, unit

END