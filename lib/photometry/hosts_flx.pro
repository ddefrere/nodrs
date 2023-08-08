;+
; NAME 
;   HOSTS_FLX
;   
; DESCRIPTION:
;   This function returns the flux in Jansky of a HOSTS target (11.1 microns).
;
; INPUTS:
;   input_name:    Name of the target
;
; KEYWORDS:
;   IRAS      :   Query IRAS if not found (return -1 otherwise)
;
; MODIFICATION HISTORY:
;   Version 1.0,  25-AUG-2015, Denis Defrère, denis@lbti.org
;   Version 1.1,  25-OCT-2016, DD: added keyword IRAS

FUNCTION HOSTS_FLX, input_name, IRAS=iras

  ; Read input file
  sep       = PATH_SEP()
  host_list = 'nodrs' + sep + 'input' + sep + 'hosts_flx.txt'
  READ_TABLE, host_list, hd, name, flx, STRING_ARRAY=[1], SEPARATOR=';', FIRST=1
  
  ; Remove spaces
  name_in = STRCOMPRESS(STRLOWCASE(input_name), /REMOVE_ALL)
  name    = STRCOMPRESS(STRLOWCASE(name), /REMOVE_ALL)
      
  ; Find input object in list
  idx = WHERE(name EQ name_in OR 'hd' + STRING(hd, FORMAT='(I0)') EQ name_in, n_tgt)
  IF n_tgt EQ 1 THEN flx_out = flx[idx] ELSE flx_out = -1
  
  ; Look in calibrator list if not found
  IF flx_out EQ -1 THEN BEGIN
    host_list = 'nodrs' + sep + 'input' + sep + 'calib_flx.txt'
    READ_TABLE, host_list, hd, name, flx, STRING_ARRAY=[1], SEPARATOR=';', FIRST=1
  
    ; Find input object in list
    idx = WHERE(name EQ name_in OR 'hd' + STRING(hd, FORMAT='(I0)') EQ name_in, n_tgt)
    IF n_tgt EQ 1 THEN flx_out = flx[idx] ELSE flx_out = -1
    
    ; Query IRAS catalog if not found
    IF flx_out EQ -1 AND KEYWORD_SET(IRAS) THEN BEGIN
      n_tgt = 1
      query_iras, input_name, iras12, iras60
      ; Do color correction
      color_c = 1.43       ; temporary, should be teff dependent
      iras12 /= color_c
      ; Extrapolate to N'
      flx_out = 1.2*iras12 ; 1.2 holds for stars between 3500K and 10000K
    ENDIF
  ENDIF
  IF n_tgt LE 0 THEN PRINT, 'Input target does not exist in the data base (' + name_in + ')'  
  IF n_tgt GT 1 THEN PRINT, 'Input target exists with various definitions. Pick one above.'
  
  RETURN, flx_out
END