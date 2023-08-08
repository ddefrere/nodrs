; +
; NAME: PLOT_FIELD
; 
; PURPOSE:
;  Plot input field vs file numer for a given date
;
; INPUTS:
;   date       :  String vector with the path to the fits files with data.
;   field      :  Field name
;
; MODIFICATION HISTORY:
;   Version 1.0,  14-JAN-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu


PRO PLOT_LBTI_FIELD, date, field
  
  ; Recover the IDL running path
  DECLARE_PATH, pth, INSTRUM='NOMIC'

  ; Critical parameters
  IF STRLEN(date) NE 6 AND NOT KEYWORD_SET(DATA_PATH) THEN MESSAGE, 'Invalid input date: must be composed of 6 digits'
  IF NOT KEYWORD_SET(DATA_PATH) THEN data_path = pth.root_data + date + pth.sep
  
  ; Retrieve FITS files in the input directory
  data_files = FILE_SEARCH(data_path,'*.fits', COUNT=n_files) 
  IF n_files LT 1 THEN MESSAGE, 'Input data path empty, skipping it.' 
  
  ; Prepare and loop
  data = DBLARR(n_files)
  FOR i=0, n_files-1 DO BEGIN
    header  = HEADFITS(data_files[i], /SILENT)
    IF !err GT 0 THEN data[i] = FLOAT(SXPAR(header, field, /NOCONTINUE)) 
  ENDFOR
  
  ; Plot
  PRINT, data
  PLOT, data
END
