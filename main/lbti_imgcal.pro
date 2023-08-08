;+
; NAME: LBTI_CALIMG
;
; PURPOSE
;   Main procedure to create calibration images for a given date (BPM, DARKS, FLATS, and LINEARITY)
;
; MANDATORY INPUTS
;   date          :  String with the UT date to be reduced based on the format 'yymmdd' (e.g., '130524')
;   
; OPTIONAL INPUT KEYWORDS
;   DARK_PATH     :  String vector pointing to the the path of dark data (superseed dark_idx keyword)
;   FLAT_PATH     :  String vector pointing to the the path of flat data (superseed flat_idx keyword)
;   
; MODIFICATION HISTORY:
;   Version 1.0, 17-APR-2017, by Denis Defrère, University of Arizona, denis@lbti.org
;   Version 1.1, 29-JUL-2018, DD: now include the bad pixel map of each pointing in the master bad pixel map!

PRO LBTI_IMGCAL, date, DARK_PATH=dark_path, FLAT_PATH=flat_path, INFO=info

  COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs, log
  ON_ERROR, 0
  
; If one of the file is not set, re-do all of them
IF NOT FILE_TEST(pth.dark_path + date + '*m_dark.fits') OR NOT FILE_TEST(pth.bpm_path + date + '*m_bpm.fits') OR NOT FILE_TEST(pth.flat_path + date + '*m_flat.fits')  THEN BEGIN
  ; --- Create master dark file if doesn't already exists (DARKs first because used for BPM later)
  IF info GT 0 THEN PRINT, '   Computing new master dark file.'
  IF KEYWORD_SET(DARK_PATH) THEN BEGIN
    IF info GT 0 THEN PRINT, '   Processing external dark frames.'
    med_drk = LBTI_READDATA(dark_path, HDR_DATA=hdr_drk, INFO=info, PLOT=plot, /MEAN)
  ENDIF ELSE BEGIN
    idx_drk = WHERE(log.datatype EQ 1 AND log.cv LT cnf.max_adu[1], n_drk) ; AND log.bad_id NE 1
    IF N_ELEMENTS(idx_drk) GT 1 THEN med_drk = LBTI_READDATA(pth.data_path, DATA_IDX=idx_drk, HDR_DATA=hdr_drk, MLOG_DATA=log, INFO=info, PLOT=plot, /MEAN) $
    ELSE PRINT, '   !! No dark frames found in data directory !!'
  ENDELSE
  IF (SIZE(med_drk))[0] NE 0 THEN LBTI_SAVECALIMG, med_drk, hdr_drk, PATH=pth.dark_path, TAG='m_dark' ELSE MESSAGE, 'No files found to compute master dark', /CONTINUE
  
  ; --- Now create BPM
  IF info GT 0 THEN PRINT, '   Computing new master bad pixel map.'
  ; First compute hot bad pixel map (high dark current pixels, from DARKS)
  idx_drk = WHERE(log.datatype EQ 1 AND log.cv LT cnf.max_adu[1], n_drk) ;AND log.bad_id NE 1
  IF n_drk GE 1 THEN BEGIN
    bpm_img = LBTI_READDATA(pth.data_path, DATA_IDX=idx_drk, HDR_DATA=m_hdr_bpm, MLOG_DATA=log, INFO=info, PLOT=plot, /SCALE)
    bpm_hot = BPM_HOT(bpm_img, INFO=info)
  ENDIF ELSE bpm_hot = 1
  ; Then compute cold bad pixel map (no flux response, from FLATS)
  idx_flt = WHERE(log.datatype EQ 2 AND log.cv LT cnf.max_adu[1], n_flt) ;  AND log.bad_id NE 1
  IF n_flt GT 1 THEN BEGIN
    ; Read flat files
    bpm_img = LBTI_READDATA(pth.data_path, DATA_IDX=idx_flt, HDR_DATA=m_hdr_bpm, MLOG_DATA=log, INFO=info, PLOT=plot, /MEDIAN, /SCALE) ; /SCALE because RANGE below is for one coadd
    ; Subtract corresponding DARKs
    IF (SIZE(med_drk))[0] NE 0 THEN BEGIN
      n_bpm = N_ELEMENTS(bpm_img[0,0,*])
      FOR i=0, n_bpm-1 DO BEGIN
        idx = WHERE(m_hdr_bpm[i].cfg_id EQ hdr_drk.cfg_id, n_drk)
        IF n_drk GT 0 THEN bpm_img[*,*,i] = bpm_img[*,*,i] - med_drk[*,*,idx]
      ENDFOR
    ENDIF
    bpm_cld = BPM_COLD(bpm_img, INFO=info)
  ENDIF ELSE bpm_cld = 1
  ; Merge bad pixel maps
  m_bpm_map = bpm_cld * bpm_hot
  IF (SIZE(m_bpm_map))[0] NE 0 THEN LBTI_SAVECALIMG, m_bpm_map, m_hdr_bpm, PATH=pth.bpm_path, TAG='m_bpm' ELSE MESSAGE, 'No files found to compute master bad pixel map', /CONTINUE
  
  ; --- Now FLAT
  IF info GT 0 THEN PRINT, '   Computing new master flat file.'
  IF KEYWORD_SET(FLAT_PATH) THEN BEGIN
    IF info GT 0 THEN PRINT, '   Processing external flat frames...'
    med_flt = LBTI_READDATA(flat_path, HDR_DATA=hdr_flt, INFO=info, /MEAN)
    med_flt = LBTI_IMGFLAT(TEMPORARY(med_flt), HDR_FLT=hdr_flt, MED_DRK=med_drk, HDR_DRK=hdr_drk)
  ENDIF ELSE BEGIN
    idx_flt = WHERE(log.datatype EQ 2 AND log.cv LT cnf.max_adu[1], n_flt) ;AND log.bad_id NE 1
    IF N_ELEMENTS(idx_flt) GT 1 THEN BEGIN
      med_flt = LBTI_READDATA(pth.data_path, DATA_IDX=idx_flt, HDR_DATA=hdr_flt, MLOG_DATA=log, INFO=info, PLOT=plot, /MEAN)
      med_flt = LBTI_IMGFLAT(TEMPORARY(med_flt), HDR_FLT=hdr_flt, IMG_DRK=med_drk, HDR_DRK=hdr_drk)
    ENDIF ELSE PRINT, '   !! No flat frames found in data directory !!'
  ENDELSE
  IF (SIZE(med_flt))[0] NE 0 THEN LBTI_SAVECALIMG, med_flt, hdr_flt, PATH=pth.flat_path, TAG='m_flat' ELSE MESSAGE, 'No files found to compute master flat', /CONTINUE
ENDIF 


; --- Compute bpm, dark, and flat calibration images for each pointing (nulling only)
IF info GT 0 THEN PRINT, '   Computing flat files (per pointing).'
; Loop over the pointings
pt_id    = log.pt_id
bad_id   = log.bad_id
datatype = log.datatype
obstype  = log.obstype
n_xpix   = log.n_xpix
n_ypix   = log.n_ypix
pt_uniq  = pt_id[UNIQ(pt_id,  SORT(pt_id))]
n_pt     = N_ELEMENTS(pt_uniq)
FOR i_pt = 0, n_pt-1 DO BEGIN

  ; First, make sure there are flat frames in this pointing (obstype = 3 + make sure it's not a dark with a wrong obstype)
  ; If they are and the file doesn't already exists, compute a flat image for this pointing
  idx_flt = WHERE(pt_id EQ pt_uniq[i_pt] AND bad_id NE 1 AND obstype EQ 3 AND datatype NE 1, n_flt)
  IF NOT FILE_TEST(pth.flat_path + date + '*_p' + STRING(pt_uniq[i_pt], FORMAT='(I0)') + '.fits') OR NOT FILE_TEST(pth.bpm_path + date + '*_p' + STRING(pt_uniq[i_pt], FORMAT='(I0)') + '.fits') THEN BEGIN  
      ; Find dark frames and compute dark for this pointing (if not found, master dark is used)
      idx_drk = WHERE(pt_id EQ pt_uniq[i_pt] AND bad_id NE 1 AND datatype EQ 1, n_drk)
      IF n_drk GT 0 THEN BEGIN
        med_drk = LBTI_READDATA(pth.data_path, DATA_IDX=idx_drk, HDR_DATA=hdr_drk, MLOG_DATA=log, INFO=info, PLOT=plot, /MEAN)
        IF (SIZE(med_drk))[0] NE 0 THEN LBTI_SAVECALIMG, med_drk, hdr_drk, PATH=pth.dark_path, TAG='dark_p' + STRING(pt_uniq[i_pt], FORMAT='(I0)')
      ENDIF ELSE BEGIN
        IF info GT 0 THEN PRINT, ' No dark frames in pointing ' + STRING(pt_uniq[i_pt], FORMAT='(I0)') + '. Using master dark instead.'
        LBTI_READCALIMG, pth.dark_path, date, [n_xpix[idx_flt[0]],n_ypix[idx_flt[0]]], -1, med_drk, hdr_drk   ; -1 to read master dark
        IF (SIZE(med_drk))[0] EQ 0 THEN MESSAGE, 'No master dark found! Data cannot be flat fielded.', /CONTINUE
      ENDELSE
      
      ; Compute bad pixel map
      ; First compute hot bad pixel map (high dark current pixels, from DARKS)
      ; If not found, read the master bad pixel map as hot bad pixel map
      idx_drk = WHERE(pt_id EQ pt_uniq[i_pt] AND datatype EQ 1 AND bad_id NE 1, n_drk)
      IF n_drk GE 1 THEN BEGIN
        bpm_img = LBTI_READDATA(pth.data_path, DATA_IDX=idx_drk, HDR_DATA=hdr_bpm, MLOG_DATA=log, INFO=info, PLOT=plot, /SCALE)
        bpm_hot = BPM_HOT(bpm_img, INFO=info)
      ENDIF ELSE LBTI_READCALIMG, pth.bpm_path, date, [n_xpix[idx_flt[0]],n_ypix[idx_flt[0]]], -1, bpm_hot, hdr_bpm
      ; Then compute cold bad pixel map (no flux response, from FLATS)
      IF n_flt GE 1 THEN BEGIN
        bpm_img = LBTI_READDATA(pth.data_path, DATA_IDX=idx_flt, HDR_DATA=hdr_bpm, MLOG_DATA=log, INFO=info, PLOT=plot, /MEDIAN, /SCALE) ; /SCALE because RANGE below is for one coadd
        ; Subtract corresponding DARKs
        IF (SIZE(med_drk))[0] NE 0 THEN BEGIN
          n_bpm = N_ELEMENTS(bpm_img[0,0,*])
          FOR i=0, n_bpm-1 DO BEGIN
            idx = WHERE(hdr_bpm[i].cfg_id EQ hdr_drk.cfg_id, n_drk)
            IF n_drk GT 0 THEN bpm_img[*,*,i] = bpm_img[*,*,i] - med_drk[*,*,idx]
          ENDFOR
        ENDIF
        bpm_cld = BPM_COLD(bpm_img, INFO=info)
      ENDIF ELSE bpm_cld = 1
      ; Merge bad pixel maps
      bpm_map = bpm_cld * bpm_hot
      IF (SIZE(bpm_map))[0] NE 0 THEN BEGIN
        LBTI_SAVECALIMG, bpm_map, hdr_bpm, PATH=pth.bpm_path, TAG='bpm_p' + STRING(pt_uniq[i_pt], FORMAT='(I0)')
        IF (SIZE(bpm_map))[1] EQ (SIZE(m_bpm_map))[1] THEN m_bpm_map *= bpm_map
      ENDIF ELSE MESSAGE, 'No bad pixels found for this pointing. Skipped.', /CONTINUE
         
      ; Now compute the flat image
      IF (SIZE(med_drk))[0] NE 0 AND n_flt GE 1 THEN BEGIN
        med_flt = LBTI_READDATA(pth.data_path, DATA_IDX=idx_flt, HDR_DATA=hdr_flt, MLOG_DATA=log, INFO=info, PLOT=plot, /MEAN)
        med_flt = LBTI_IMGFLAT(TEMPORARY(med_flt), HDR_FLT=hdr_flt, IMG_DRK=med_drk, HDR_DRK=hdr_drk)
        IF (SIZE(med_flt))[0] NE 0 THEN LBTI_SAVECALIMG, med_flt, hdr_flt, PATH=pth.flat_path, TAG='flat_p' + STRING(pt_uniq[i_pt], FORMAT='(I0)')
      ENDIF ELSE MESSAGE, ' No dark frame for pointing ' + STRING(pt_uniq[i_pt], FORMAT='(I0)'), /CONTINUE
    ENDIF ELSE IF info GT 0 THEN PRINT, ' Calibration frames for pointing ' + STRING(pt_uniq[i_pt], FORMAT='(I0)') + ' already exist'
ENDFOR

; Save master bad pixel map
IF (SIZE(m_bpm_map))[0] NE 0 THEN LBTI_SAVECALIMG, m_bpm_map, m_hdr_bpm, PATH=pth.bpm_path, TAG='m_bpm'

END