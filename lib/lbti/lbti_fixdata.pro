;+
; NAME: LBTI_FIXDATA
;
; PURPOSE:
;   This procedure runs basic consistency checks on the masterlog file of a given date (and fixes it if necessary). In particular, it checks that:
;     - OBJNAME is the same within a given pointing;
;     - DARKs are associated to the right pointing;
;     - Flagg as BAD background frames followed b background frames of a different nod;
;     - OBSTYPE is the same within a given nod;
;     - Check that data are recorded in chronological order
;
; MADATORY INPUT
;   data_path     :  String vector with the path to the fits files with data.
;
; KEYWORDS
;   LMIRCAM       : Set for LMIRCam
;   FIX_ORDER     : Set to fix the data file order (chronologically)
;
; MODIFICATION HISTORY:
;   Version 1.0,  08-JUL-2013, by Denis Defrère, University of Arizona, denis@lbti.org
;   Version 1.1,  04-OCT-2016, DD: updated for new formalism
;   Version 1.2,  07-APR-2017, DD: now backup initial log files

PRO LBTI_FIXDATA, data_path, FIX_ORDER=fix_order

  ; Retrieve FITS files in the input directory
  data_files = FILE_SEARCH(data_path,'*.fits') & n_files = N_ELEMENTS(data_files)
  IF n_files LT 1 THEN MESSAGE, 'Input data path empty'

  ; --- Read masterlog file (create it if necessary)
  mlog_file = data_path + 'masterlog.dat'
  IF NOT FILE_TEST(mlog_file) THEN PRINT, 'Creating new masterlog file. This may take a while...'
  IF NOT FILE_TEST(mlog_file) THEN LBTI_MASTERLOG, data_path, /IDL
  LBTI_READMASTERLOG, mlog_file, mlog_data
  
  ; --- Consitancy check
  IF n_files NE N_ELEMENTS(mlog_data.cv) THEN MESSAGE, 'Masterlog not matching the data directory', /CONTINUE

  ; --- Initiate time variable
  t0 = SYSTIME(1)
  
  ; --- Extract useful data
  obj_name = mlog_data.obj_name
  file_id  = mlog_data.file_id
  time_obs = mlog_data.time_obs
  nod_id   = mlog_data.nod_id
  pt_id    = mlog_data.pt_id
  datatype = mlog_data.datatype
  obstype  = mlog_data.obstype 
  flag     = mlog_data.flag
  date     = STRMID(data_path, STRLEN(data_path)-7, 6 )
    
  ; --- Convert time_obs to hour
  IF KEYWORD_SET(FIX_ORDER) THEN BEGIN
    hr  = DOUBLE(STRMID(time_obs, 0, 2)) + DOUBLE(STRMID(time_obs, 3, 2))/60D + DOUBLE(STRMID(time_obs, 6, 2))/3.6D+3 + DOUBLE(STRMID(time_obs, 9, 3))/3.6D+6
    
    ; --- Sort by increasing time
    idx = SORT(hr)
  
    ; --- Check that data is recorded in chronological order
    IF MAX(ABS(idx-LINDGEN(n_files))) NE 0 THEN BEGIN
       ; --- Create new data path
      new_path = pth.root_data + date + '_fix' + pth.sep
      IF NOT FILE_TEST(new_path) THEN FILE_MKDIR, new_path
  
      ; --- Loop over the files, rename them, and save them by increasing time
      FOR i_f = 0, n_files-1 DO BEGIN
        img_file = READFITS(data_files[idx[i_f]], hdr, /SILENT)
  
        ; Derive new file name
        filename = 'n_' + date + '_' + STRING(i_f, FORMAT='(I06)') + '.fits'
  
        ; Parse it to header
        SXADDPAR, hdr, 'FILENAME',  filename
  
        ; Save fits
        WRITEFITS, new_path + filename, img_file, hdr
  
        ;FILE_DELETE, data_files[idx[i_f]]
      ENDFOR    
    ENDIF
  ENDIF
    
  ; --- Init variable
  n_mod = 0
  
  ; --- Find bad OBJNAME (can be 0 sometimes when NOMIC dies)
  idx_bad = WHERE(VALID_NUM(obj_name), n_bad, COMPLEMENT=idx_ok)
  IF n_bad GT 0 AND n_bad NE n_files THEN BEGIN
    ; Because LBTI_FITSUTIL can be very slow on large directory, find first consecutive BAD file id
    idx_gr = [WHERE(file_id[idx_bad] - SHIFT(file_id[idx_bad], 1) NE 1, n_group), n_bad]
    FOR ig = 0, n_group-1 DO BEGIN
      idx_cur   = idx_bad[idx_gr[ig] + LINDGEN(idx_gr[ig+1]-idx_gr[ig])]
      hdr_cur   = HEADFITS(data_path + '*' + STRING(file_id[idx_cur[0]], FORMAT='(I06)') + '*.fits')
      ra_cur    = TEN(SXPAR(hdr_cur, 'LBT_RA'))
      fid_prev  = file_id[idx_ok[MAX(WHERE(idx_ok LT idx_cur[0], n_prev))]]
      fid_aft   = file_id[idx_ok[MIN(WHERE(idx_ok GT idx_cur[0], n_aft))]]
      ; If good files before and after, replace
      IF n_prev GE 1 AND n_aft GE 1 THEN BEGIN
        hdr_prev  = HEADFITS(data_path + '*' + STRING(fid_prev, FORMAT='(I06)') + '*.fits')
        ra_prev   = TEN(SXPAR(hdr_prev, 'LBT_RA'))
        hdr_aft   = HEADFITS(data_path + '*' + STRING(fid_aft, FORMAT='(I06)') + '*.fits')
        ra_aft    = TEN(SXPAR(hdr_aft, 'LBT_RA'))
        IF ABS(ra_aft-ra_cur) LE ABS(ra_prev-ra_cur) THEN name = SXPAR(hdr_aft, 'OBJNAME') ELSE name = SXPAR(hdr_prev, 'OBJNAME')
        LBTI_FITSUTIL, date, [MIN(file_id[idx_cur]),MAX(file_id[idx_cur])], OBJNAME=name
        n_mod += N_ELEMENTS(idx_cur)
      ENDIF
    ENDFOR
  ENDIF
    
  ; --- Replace OBJNAME if same pointing
  data1 = pt_id
  data2 = obj_name  
  data1_un = data1[UNIQ(data1,  SORT(data1))]
  n_data1  = N_ELEMENTS(pt_uniq)
  FOR i = 0, n_data1-1 DO BEGIN 
    d1_idx = WHERE(data1 EQ data1_un[i], n)
    f1_fid = file_id[d1_idx]
    d2_tmp = data2[d1_idx]
    data2_un = d2_tmp[UNIQ(d2_tmp, SORT(d2_tmp))]
    n_data2  = N_ELEMENTS(data2_un)
    ; If more than one object name for this pointing, find the most frequent one
    IF n_data2 GT 1 THEN BEGIN
      n0 = 0
      FOR j = 0, n_data2-1 DO BEGIN
        idx_tmp = WHERE(d2_tmp EQ data2_un[j], n)
        IF n GT n0 THEN data2_max = data2_un[j]
      ENDFOR
      ; Now replace the bad entries
      idx_bad  = WHERE(d2_tmp NE data2_max, n_bad)
      fid_bad = f1_fid[idx_bad]
      FOR i_b = 0, n_bad-1 DO LBTI_FITSUTIL, date, [fid_bad[i_b],fid_bad[i_b]], OBJNAME=data2_max
      ; Increase counter
      n_mod += n_bad
    ENDIF
  ENDFOR 
  
  ; --- Assign correct pointing to DARKs (based on RA and DEC)
  idx_drk = WHERE(datatype EQ 1, n_drk)
  IF n_drk GT 1 THEN BEGIN
    nod_drk = nod_id[idx_drk]
    nod_un  = nod_drk[UNIQ(nod_drk, SORT(nod_drk))]
    n_seq   = N_ELEMENTS(nod_un)
    FOR i = 0, n_seq - 1 DO BEGIN
      idx_cur = WHERE(nod_id EQ nod_un[i] AND datatype EQ 1, n_cur)
      IF n_cur GT 1 THEN BEGIN
        file_id_prev = file_id[MIN(idx_cur)] - 1 & idx_prev = WHERE(file_id_prev EQ file_id, n_prev)
        file_id_aft  = file_id[MAX(idx_cur)] + 1 & idx_aft  = WHERE(file_id_aft EQ file_id, n_aft)
        IF n_prev GE 1 AND n_aft GE 1 THEN BEGIN                            ; If one file before and after the dark sequence, continue
          IF datatype[idx_prev] NE 1 AND datatype[idx_aft] NE 1 THEN BEGIN  ; If no other darks in adjacent nods, continue
            IF obj_name[idx_prev] NE obj_name[idx_aft] THEN BEGIN           ; If different OBJECTNAME continue and compare LBT_RA and LBT_DEC
              header    = HEADFITS(data_path + '*' + STRING(file_id_prev, FORMAT='(I06)') + '*.fits')
              ra_prev   = TEN(SXPAR(header, 'LBT_RA'))
              flag_prev = SXPAR(header, 'FLAG')
              header    = HEADFITS(data_path + '*' + STRING(file_id_aft, FORMAT='(I06)') + '*.fits')
              ra_aft    = TEN(SXPAR(header, 'LBT_RA'))
              flag_aft  = SXPAR(header, 'FLAG')
              header    = HEADFITS(data_path + '*' + STRING(file_id[MEAN(idx_cur)], FORMAT='(I06)') + '*.fits')
              ra_cur    = TEN(SXPAR(header, 'LBT_RA'))
              IF ra_cur LT 24 THEN BEGIN  ; Only for valid RA
                IF ABS(ra_prev-ra_cur) LE ABS(ra_aft-ra_cur) THEN BEGIN
                  obj_tmp  = obj_name[idx_prev] 
                  flag_tmp = flag_prev
                ENDIF ELSE BEGIN
                  obj_tmp  = obj_name[idx_aft]
                  flag_tmp = flag_aft
                ENDELSE
                ; Replace if new and valid name
                IF obj_tmp NE obj_name[idx_cur] AND NOT VALID_NUM(obj_tmp) THEN BEGIN
                  LBTI_FITSUTIL, date, [MIN(file_id[idx_cur]),MAX(file_id[idx_cur])], OBJNAME=obj_tmp, FLAG=flag_tmp
                  n_mod += N_ELEMENTS(idx_cur)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDFOR
  ENDIF
  
  ; --- Now look for successive BACKGROUND nods (i.e., OBSTYPE = 3), merge them and flag the first ones as BAD!
  ; Ignore frames already flaggd as BAD or this will mess up everything!!
  idx_bckg = WHERE(obstype EQ 3 AND flag NE 'B' AND datatype NE 1 AND datatype NE 2, n_bckg)
  IF n_bckg GE 1 THEN BEGIN
    nod_bckg = nod_id[idx_bckg]
    nod_un   = nod_bckg[UNIQ(nod_bckg, SORT(nod_bckg))]
    n_seq   = N_ELEMENTS(nod_un)
    FOR i = 1, n_seq - 1 DO BEGIN
      ; If different nod IDs, continue
      IF nod_un[i] EQ nod_un[i-1]+1 THEN BEGIN
        idx_cur = WHERE(nod_id EQ nod_un[i])
        idx_pre = WHERE(nod_id EQ nod_un[i-1])
        ; If same OBJNAME and DATATYPE, continue
        IF obj_name[idx_cur[0]] EQ obj_name[idx_pre[0]] THEN BEGIN
          LBTI_FITSUTIL, date, [MIN(file_id[idx_pre]),MAX(file_id[idx_pre])], FLAG='BAD' 
          n_mod += N_ELEMENTS(idx_pre)
        ENDIF
      ENDIF     
    ENDFOR
  ENDIF
    
  ; Create new masterlog if modifications
  If n_mod GT 0 THEN BEGIN
    IF NOT FILE_TEST(data_path + 'masterlog_bu.dat') THEN FILE_COPY, data_path + 'masterlog.dat', data_path + 'masterlog_bu.dat'
    IF NOT FILE_TEST(data_path + 'config_id_bu.dat') THEN FILE_COPY, data_path + 'config_id.dat', data_path + 'config_id_bu.dat'
    IF NOT FILE_TEST(data_path + 'nod_id_bu.dat')    THEN FILE_COPY, data_path + 'nod_id.dat', data_path + 'nod_id_bu.dat'
    LBTI_MASTERLOG, data_path, /IDL
  ENDIF
 
  ; --- Print info to screen
  PRINT, 'Number of modified files', n_mod
  PRINT, 'Time to fix header and create new masterlog [s]:', SYSTIME(1)-t0
END