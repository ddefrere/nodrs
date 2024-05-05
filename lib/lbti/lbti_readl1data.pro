;+
; NAME: LBTI_READL1DATA
; 
; PURPOSE:
;   Read a processed L1 file given an input data path, OB ID, and data type (e.g., NULL, PHOT, BCKG).
;
; INPUTS:
;   data_path   
;   ob_id     
;   flag
;
; INPUT KEYWORDS
;   APER       :   
;   FILTER     :
;   
; OUTPUT
;   A data structure containing the flux measurements and related information (e.g., error)
;
; MODIFICATION HISTORY:
;   Version 1.0,  26-OCT-2016, Denis Defrère
;   Version 1.1,  05-MAY-2024, DD: now return the header

FUNCTION LBTI_READL1DATA, data_path, ob_id, flag, APER=aper, FILTER=filter, IDX_APER=i_aper, HEADER=header

  ; First find file
  file = FILE_SEARCH(data_path,'*ID'+STRING(ob_id,FORMAT='(I03)') + '*' + flag +  '.fits', COUNT=n_files)
  IF n_files LT 1 THEN RETURN, -1 ; not found
  IF n_files GT 1 THEN MESSAGE, 'Found multiple ' + flag + ' files. Keeping only the first one', /CONTINUE
    
  ; Read the data
  data   = MRDFITS(file[0], 1, /SILENT)
  header = HEADFITS(file[0])
  
  ; Check whether the same photometric aperture is available
  IF KEYWORD_SET(APER) THEN BEGIN
    ; Find the right aperture radius
    n_apr   = N_ELEMENTS(data[0].flx_tot)
    apr_all = FXPAR(header, 'APERRAD' + STRING(0, FORMAT='(I0)'))
    FOR i_r = 1, n_apr-1 DO apr_all = [apr_all, FXPAR(header, 'APERRAD' + STRING(i_r, FORMAT='(I0)'))]
    IF MIN(ABS(apr_all-aper), i_aper) NE 0 THEN BEGIN
      MESSAGE, 'No corresponding aperture radius for ' + flag + ' files. Skipping this OB.', /CONTINUE
      RETURN, -1
    ENDIF 
  ENDIF 
    
  ; Now filter the data if requested (only works if one aperture)
  ; Remove obvious outliers from sequence and sanity check (do it twice)
  IF KEYWORD_SET(FILTER) AND KEYWORD_SET(aper) THEN BEGIN
    ; First filter bassed on flux 
    flx = data.flx_tot[i_aper]
    AVGSDV, flx, avg, rms, KAPPA=5 
    idx_ok  = WHERE(ABS(flx-avg) LE 5*rms)
    flx = data[idx_ok].flx_tot[i_aper]
    AVGSDV, flx, avg, rms, KAPPA=5
    idx_ok  = idx_ok[WHERE(ABS(flx-avg) LE 5*rms)]
    ; Now filter based on background
    flx = data[idx_ok].bck_tot[i_aper]
    AVGSDV, flx, avg, rms, KAPPA=5
    idx_ok  = idx_ok[WHERE(ABS(flx-avg) LE 5*rms)]
    flx = data[idx_ok].bck_tot[i_aper]
    AVGSDV, flx, avg, rms, KAPPA=5
    idx_ok = idx_ok[WHERE(ABS(flx-avg) LE 5*rms)]
    data   = data[idx_ok]    
  ENDIF

RETURN, data
END
