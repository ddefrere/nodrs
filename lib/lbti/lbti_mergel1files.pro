;+
; NAME: LBTI_MERGEL1FILES
; 
; PURPOSE:
;   This procedure merges reduced L1 image cubes
;
; INPUTS:
;   l1files
;
; KEYWORDS:
;
; MODIFICATION HISTORY:
;   Version 1.0, 02-NOV-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (adapted from former routine 'remove_bck.pro')
;   Version 1.1, 04-APR-2015, DD: added keyword FILENAME
;   Version 1.2, 29-NOV-2015, DD: added wavelength to FILENAME

PRO LBTI_MERGEL1FILES, l1files, FILENAME=filename

; Derive number of files
n_files = N_ELEMENTS(l1files)

; Read first file
img0 = MRDFITS(l1files[0], 0, hdr0, /SILENT)
dat0 = MRDFITS(l1files[0], 1, /SILENT)
n_img    = N_ELEMENTS(img0[0,0,*])
date_obs = FXPAR(hdr0, 'DATE_OBS', /NOCONTINUE) 
objname  = STRCOMPRESS(FXPAR(hdr0, 'OBJNAME', /NOCONTINUE), /REMOVE_ALL)
flag     = STRCOMPRESS(FXPAR(hdr0, 'FLAG', /NOCONTINUE), /REMOVE_ALL)
wav      = FXPAR(hdr0, 'WAVELENG', /NOCONTINUE)

; Loop over the remaining files and append
FOR i_f = 1, n_files-1 DO BEGIN
  img_tmp  = MRDFITS(l1files[i_f], 0, hdr1, /SILENT)
  n_img    = N_ELEMENTS(img_tmp[0,0,*])
  img0     = [[[TEMPORARY(img0)]],[[img_tmp]]]
  ; Now data series
  dat_tmp  = MRDFITS(l1files[i_f], 1, /SILENT)
  dat0     = [TEMPORARY(dat0),dat_tmp]
ENDFOR

; Add slope to data if not defined
IF NOT TAG_EXIST(dat0, 'SLOPE') THEN struct_add_field, dat0, 'slope', FLTARR(N_ELEMENTS(dat0.mjd_obs))

; Save file
IF NOT KEYWORD_SET(filename) THEN filename = FILE_DIRNAME(l1files[0]) + '/' + date_obs + '_' + flag + '_' + objname + '_' + STRING(1D+6*wav, FORMAT='(I0)') + 'um'
MWRFITS, img0, filename + '_IMG_ALL.fits', hdr0, /CREATE, /SILENT  

; Append extension with parameters that change frame by frame
n_row = N_ELEMENTS(dat0.mjd_obs)
outfile = filename + '_DATA_ALL.fits'
FXHMAKE,  hdr, /INIT, /EXTEND, 0
FXWRITE, outfile, hdr
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'DATA_SERIES', 'Frame to frame variable parameters'

; Init column number
n_col = 12
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, dat0[0].mjd_obs,    'MJD_OBS',   'Modified Julian Date of observation'
FXBADDCOL, 2L, hdr, dat0[0].lbt_utc,    'LBT_UTC',   'UTC from observatory'
FXBADDCOL, 3L, hdr, dat0[0].lbt_lst,    'LBT_LST',   'LST from observatory'
FXBADDCOL, 4L, hdr, dat0[0].lbt_alt,    'LBT_ALT',   'ALT from observatory'
FXBADDCOL, 5L, hdr, dat0[0].lbt_az,     'LBT_AZ',    'AZ from observatory'
FXBADDCOL, 6L, hdr, dat0[0].lbt_para,   'LBT_PARA',  'Parralactic angle from obs.'
FXBADDCOL, 7L, hdr, dat0[0].file_id,    'FILE_ID',   'File identification number'
FXBADDCOL, 8L, hdr, dat0[0].nod_id,     'NOD_ID',    'Nod identification number'
FXBADDCOL, 9L, hdr, dat0[0].chp_id,     'CHP_ID',    'Chop identification number'
FXBADDCOL, 10L, hdr, dat0[0].xcen,      'XCEN',      'X position of the star'
FXBADDCOL, 11L, hdr, dat0[0].ycen,      'YCEN',      'Y position of the star'
FXBADDCOL, 12L, hdr, dat0[0].slope,     'SLOPE',     'Slope of the fitted Moffat model'

; Write extension header to FITS file
FXBCREATE, unit, outfile, hdr
FXBWRITM,  unit, col, dat0.mjd_obs, dat0.lbt_utc, dat0.lbt_lst, dat0.lbt_alt, dat0.lbt_az, dat0.lbt_para, dat0.file_id, $
                      dat0.nod_id, dat0.chp_id, dat0.xcen, dat0.ycen, dat0.slope
FXBFINISH, unit
END