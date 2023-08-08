; +
; NAME: LBTI_READMASTERLOG
;
; PURPOSE:
;   Read a masterlog file and return a structure with all the information 
;
; INPUTS:
;   mlog_file  : The masterlog file to read (including full path). 
;   
; OPTIONAL INPUT KEYWORDS
;   BAD_IDX       :  Vector with the file number of frames to be removed
;   BCKG_IDX      :  Two-element vector with the lower and the upper file numbers of the background files (e.g., [10,19]). Generally used when BCKG_MODE = 4.
;   DARK_IDX      :  Two-element vector with the lower and the upper file numbers of the dark files (e.g., [0,9])
;   DATA_IDX      :  Two-element vector with the lower and the upper file numbers of the data files (e.g., [20,999])
;   FLAT_IDX      :  Two-element vector with the lower and the upper file numbers of the flat files (e.g., [10,19])
;   NOD_IDX       :  2x(number of nods) array with the lower and the upper file numbers of each nod position (e.g., [[120,239],[240,299],[300,399]])
;
; OUTPUT
;   mlog_data  : A structure containing the masterlog data 
;
; MODIFICATION HISTORY:
;   Version 1.0,  14-JAN-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  26-MAY-2016, DD: added central value

PRO LBTI_READMASTERLOG, mlog_file, mlog_data, BAD_IDX=bad_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx

  ; --- Read the masterlog file
  READ_TABLE, mlog_file, file_id, obj_name, time_obs, lbtalt, wav_eff, dit, pt_id, cfg_id, nod_id, chp_id, obstype, datatype, cv, n_xpix, n_ypix, flag, FIRST=2, SKIP=[15], STRING_ARRAY=[1,2,15], SEP=';'

  ; --- Prepare bad ID for later
  bad_id = 0*nod_id  ; assume all files are good for now
  
  ; --- Superseed entries based on input keywords
  IF KEYWORD_SET(BAD_IDX) THEN BEGIN
    idx_bad = VALUE_LOCATE(file_id, bad_idx)
    IF MIN(idx_bad) GE 0 THEN bad_id[idx_bad] = 1
  ENDIF
  IF KEYWORD_SET(BCKG_IDX) THEN BEGIN
    idx_bckg = WHERE(file_id GE MIN(bckg_idx) AND file_id LE MAX(bckg_idx) AND bad_id NE 1, n_bckg)
    IF MIN(idx_bckg) GE 0 THEN obstype[idx_bckg] = 3   ; Default value for background frames
  ENDIF
  IF KEYWORD_SET(DARK_IDX) THEN BEGIN
    idx_dark = WHERE(file_id GE MIN(dark_idx) AND file_id LE MAX(dark_idx) AND bad_id NE 1, n_drk)
    IF n_drk GE 1 THEN datatype[idx_dark] = 1
  ENDIF
  IF KEYWORD_SET(FLAT_IDX) THEN BEGIN
    idx_flat = WHERE(file_id GE MIN(flat_idx) AND file_id LE MAX(flat_idx) AND bad_id NE 1, n_flt)
    IF n_flt GE 1 THEN datatype[idx_flat] = 2
  ENDIF
  IF KEYWORD_SET(NOD_IDX) THEN BEGIN
    nod_id[*] = -99
    n_nod     = N_ELEMENTS(nod_idx[0,*])
    FOR i_nod = 0, n_nod-1 DO nod_id[WHERE(file_id GE MIN(nod_idx[0,i_nod]) AND file_id LE MAX(nod_idx[1,i_nod]), /NULL)] = i_nod
  ENDIF
  
  ; --- Create output structure
  mlog_data = {FILE_ID: LONG(file_id), OBJ_NAME: obj_name, TIME_OBS: time_obs, LBTALT: lbtalt, WAV_EFF: wav_eff, DIT: dit, BAD_ID: LONG(bad_id), CFG_ID: FIX(cfg_id),$
               NOD_ID: FIX(nod_id), CHP_ID: FIX(chp_id), PT_ID: FIX(pt_id), OBSTYPE: FIX(obstype), DATATYPE: FIX(datatype), CV: cv, N_XPIX: FIX(n_xpix), N_YPIX: FIX(n_ypix), FLAG: flag}
END