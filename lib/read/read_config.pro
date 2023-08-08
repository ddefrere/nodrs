;+
; NAME: READ_CONFIG
; 
; PURPOSE:
;   Read a properly formatted config file into a structure
;
; INPUT:
;   Path to the config file
;
; OUTPUTS:
;   A structure
;
; MODIFICATION HISTORY:
;   Version 1.0,  22-MAR-2013, by Denis Defrère, Steward Observatory (ddefrere@email.arizona.edu)

FUNCTION READ_CONFIG, cfg_file
    
  ; Open config file. In case of an error go to the defaults. If that fails find a defaults file.
  IF NOT FILE_TEST(cfg_file) THEN cfg_file=DIALOG_PICKFILE(/READ, TITLE='Locate config file', FILTER='*.cfg')
  OPENR, lun, cfg_file, /GET_LUN, ERROR=err
  IF err THEN BEGIN
    MESSAGE,'Error opening config file '+cfg_file+'; starting file finding dialog',/INFO
    cfg_file=DIALOG_PICKFILE(/READ, TITLE='Locate config file', FILTER='*.cfg')
    OPENR, lun, cfg_file, /GET_LUN, ERROR=err
  ENDIF
  
  ; Read in all the parameters:
  line = ''
  WHILE NOT EOF(lun) DO BEGIN
    ; Read line
    READF, lun, line
    ; Check whether it's a good line
    sym_pos = STRPOS(line, ';')      ; Extract until comment symbol
    str_len = STRLEN(line)           ; Lenght of the line
    ; Skip if comment line or empty
    IF sym_pos EQ 0 OR str_len EQ 0 THEN GOTO, skip_line
    ; Extract relevant part (befre comment if comment symbol is found)
    IF sym_pos GT 0 THEN line = STRMID(line, 0, sym_pos)
    ; Split line and assign values
    str   = STRSPLIT(line, /EXTRACT)
    field = STRTRIM(str[0], 2)
    value = STRTRIM(str[1], 2)
    ; Look for the ' or " character to determine whether it's a STRING.
    ; If it's not a STRING, look for the , character to derive the number of elements.
    IF STRMATCH(value, "*'*") NE 1 AND STRMATCH(value, '*"*') NE 1 THEN BEGIN
      n_el    = N_ELEMENTS(STRSPLIT(value, ','))          ; Number of elements
      IF n_el GT 1 THEN val_for = FLTARR(n_el) ELSE val_for = 0.
      READS,value,val_for                                 ; Convert to formatted number
    ENDIF ELSE val_for = STRMID(value,1,STRLEN(value)-2)  ; Remve leading and trailing " or ' characters
    ; Append to structure
    IF (SIZE(stru))[0] NE 0 THEN stru = CREATE_STRUCT(stru, field, val_for) $
                            ELSE stru = CREATE_STRUCT(field, val_for)
    ; Jump point if comment line
    skip_line:
  ENDWHILE
  
  ; Close config file
  CLOSE, lun
  FREE_LUN, lun

RETURN, stru
END