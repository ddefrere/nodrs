;+
; NAME: LBTI_READCALIMG
; 
; PURPOSE:
;   Read a calibration image file and return the image data cube.
;   Keyword information is returned by the keyword 'HDR_DATA'.
;
; MANDATORY INPUTS:
;   path : the path to the directory with the calibration files
;   date : the date (6-digit format)
;   size : the frame size 
;   pt_id: the pointing ID 
;   
; KEYWORDS
;   CROP    :  If set, crop the frames (4-element vector: x_min, y_min, x_max, y_max)
;
; OUTPUT
;   img_out : a data cube with the images contained in the input FITS file.
;   hdr_out : the corresponding header information
;   
; MODIFICATION HISTORY:
;   Version 1.0, 12-JAN-2015, Denis Defrère, based on a more general routine LBTI_READL0
;   Version 1.1, 27-MAR-2016, DD: added pt_id
;   Version 1.2, 28-MAR-2016, DD: added keyword CROP
;   Version 1.3, 25-MAY-2016, DD: added instrument info to output data
;   Version 1.4, 11-AUG-2018, DD: now look first in neighboring pointings if no flat files found

PRO LBTI_READCALIMG, path, date, size, pt_id, img_out, hdr_out, CROP=crop

; Find the file
files = FILE_SEARCH(path + date + '*' + STRING(size[0], FORMAT='(I0)') + 'x' + STRING(size[1], FORMAT='(I0)') + '*_p' + STRING(pt_id, FORMAT='(I0)') + '.fits', COUNT=n_files)

; If no file found, take the one in the next pointing
IF n_files LE 0 THEN files = FILE_SEARCH(path + date + '*' + STRING(size[0], FORMAT='(I0)') + 'x' + STRING(size[1], FORMAT='(I0)') + '*_p' + STRING(pt_id+1, FORMAT='(I0)') + '.fits', COUNT=n_files)
IF n_files LE 0 THEN files = FILE_SEARCH(path + date + '*' + STRING(size[0], FORMAT='(I0)') + 'x' + STRING(size[1], FORMAT='(I0)') + '*_p' + STRING(pt_id-1, FORMAT='(I0)') + '.fits', COUNT=n_files)

; If still no file found, look for the master one
IF n_files LE 0 THEN files = FILE_SEARCH(path + date + '*' + STRING(size[0], FORMAT='(I0)') + 'x' + STRING(size[1], FORMAT='(I0)') + '*m*.fits', COUNT=n_files)

; If still no file found, take the most recent of other dates or pointings
IF n_files LE 0 THEN BEGIN
  files = FILE_SEARCH(path + '*' + STRING(size[0], FORMAT='(I0)') + 'x' + STRING(size[1], FORMAT='(I0)') + '*.fits', COUNT=n_files)
  IF n_files LE 0 THEN BEGIN
    MESSAGE, 'No matching calibration (DARK, FLAT, or BPM) data!', /CONTINUE
    img_out = -1
    hdr_out = -1
    RETURN 
  ENDIF
ENDIF

; Read the data
img_out = READFITS(files[0], header, /SILENT) 
hdr_out = MRDFITS(files[0], 1, junk, /SILENT) 

; Add instrument info to structure
struct_add_field, hdr_out, 'instrum', FXPAR(header, 'INSTRUM')
struct_add_field, hdr_out, 'date_obs', FXPAR(header, 'DATE_OBS')

; Crop the frame if requested
IF KEYWORD_SET(CROP) THEN img_out = TEMPORARY(img_out[crop[0]:crop[2],crop[1]:crop[3],*])

END
