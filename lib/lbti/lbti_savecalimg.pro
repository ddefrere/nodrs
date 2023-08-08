;+
; NAME: LBTI_SAVECALIMG
; 
; PURPOSE:
;   This procedure saves calibrated L0 images and relevant keywords in FITS files.
;
; INPUTS:
;   img_in    :  3-dimension array with the images to save (n_xpix, n_ypix, n_frames).
;   hdr_in    :  The structure with the corresponding header information
;
; KEYWORDS:
;   PATH      :  Full path where the FITS file will be saved
;   TAG       :  Tag name added to the file name
;
; MODIFICATION HISTORY:
;   Version 1.0, 12-JAN-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu

PRO LBTI_SAVECALIMG, img_in, hdr_in, PATH=path, TAG=tag

; Keyword sanity check
IF NOT KEYWORD_SET(PATH) THEN path = '' ELSE IF NOT FILE_TEST(path) THEN FILE_MKDIR, path
IF NOT KEYWORD_SET(TAG)  THEN tag  = ''

; Create file header
MKHDR, hdr, img_in

; Add comment to the header
SXADDPAR, hdr, "COMMENT", "Master " + tag + " image of data taken with LBTI/" + STRTRIM(hdr_in[0].instrum) + "."  
SXADDPAR, hdr, "COMMENT", "Observation parameters"
SXADDPAR, hdr, 'TELESCOP',  'LBT',                  'telescope'
SXADDPAR, hdr, 'INSTRUME',  hdr_in[0].instrum[0],   'instrument'
SXADDPAR, hdr, 'DATE_OBS',  hdr_in[0].date_obs,     'date of observation'

; Define filename and save
size_img = SIZE(img_in)
n_xpix   = size_img[1]
n_ypix   = size_img[2]
filename = STRCOMPRESS(path + hdr_in[0].date_obs + '_' + STRING(n_xpix, FORMAT='(I0)') + 'x' + STRING(n_ypix, FORMAT='(I0)') + '_' + tag + '.fits', /REMOVE_ALL)
MWRFITS, img_in, filename, hdr, /CREATE, /SILENT

; Append extension with parameters that change frame by frame
n_row = N_ELEMENTS(hdr_in.lam_cen)

; Create extension
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'DATA_SERIES', 'parameters that change frame by frame'

; Init column number
n_col = 10
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, hdr_in[0].cfg_id,       'CFG_ID',     'Config ID'
FXBADDCOL, 2L, hdr, hdr_in[0].lam_cen,      'LAM_CEN',    'Effective wavelentgh'
FXBADDCOL, 3L, hdr, hdr_in[0].bandwidt,     'BANDWIDT',   'Bandwidth'
FXBADDCOL, 4L, hdr, hdr_in[0].smplmode,     'SMPLMODE',   'Sampling mode'
FXBADDCOL, 5L, hdr, hdr_in[0].int_time,     'INT_TIME',   'Integration time [s]'
FXBADDCOL, 6L, hdr, hdr_in[0].n_coadd,      'N_COADD',    'Number of coadds'
FXBADDCOL, 7L, hdr, hdr_in[0].pagain,       'PAGAIN',     'Pre-amp gain'
FXBADDCOL, 8L, hdr, hdr_in[0].pabandw,      'PABANDW',    'Pre-amp bandwidth [Hz]'
FXBADDCOL, 9L, hdr, hdr_in[0].detmode,      'DETMODE',    'Detector mode (0:HG, 1:LG)'
FXBADDCOL, 10L, hdr, hdr_in[0].detbias,      'DETBIAS',    'Detector bias [V]'

; Write extension header to FITS file
FXBCREATE, unit, filename, hdr
FXBWRITM,  unit, col, hdr_in.cfg_id, hdr_in.lam_cen, hdr_in.bandwidt, hdr_in.smplmode, hdr_in.int_time, hdr_in.n_coadd, hdr_in.pagain, hdr_in.pabandw, hdr_in.detmode, hdr_in.detbias
FXBFINISH, unit  
         
END