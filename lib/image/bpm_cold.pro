;+
; NAME: BPM_COLD
; 
; PURPOSE:
;   Compute cold bad pixel map based on input images (computed using the spatial variance and a hard-coded threshold of 10 sigma). 
;   If data cube, concatenate bad pixels found in each image of the cube
;
; INPUTS:
;   img_in : the input image cube
;
; KEYWORDS
;   RANGE  : in ADU, range of acceptable values [min, max]
;   INFO   : set to print info to screen
;   PLOT   : obsolete
;
; OUTPUT
;   The computed bad pixel map.
;   WARNING: as of March 2016, the binarity is reversed!, i.e. 0 is a bad pixel (anything else is good pixel)
;   
; MODIFICATION HISTORY:
;   Version 1.0, 17-FEB-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 27-MAR-2016, DD: added keyword range + WARNING: reversed the binarity, i.e. 0 is a bad pixel (anything else is good pixel)
;   Version 1.2, 27-APR-2016, DD: renamed routine to BPM_COLD.pro
;   Version 1.3, 15-NOV-2016, DD: added a second box size

FUNCTION BPM_COLD, img_in, RANGE=range, INFO=info, PLOT=plot
  
  ; Keyword sanity check
  IF NOT KEYWORD_SET(INFO) THEN info = 0
  IF NOT KEYWORD_SET(PLOT) THEN plot = 0
  
  ; Running parameters
  rms_pix1   = 8
  rms_pix2   = 5
  box_width1 = 15
  box_width2 = 5

  ; Number of frames in the input image
  size_img = SIZE(img_in)
  IF size_img[0] EQ 3 THEN n_fr = size_img[3] ELSE n_fr = 1
  
  ; Init output bad pixel map
  ; Make sure its a integer
  bpm_map = 1 + INTARR(size_img[1], size_img[2])
  
  ; Loop over the frames
  bad_pix = 0
  FOR i_fr = 0, n_fr-1 DO BEGIN
    ; Find bad pixels
    img_cur  = REFORM(img_in[*,*,i_fr])                                                   ; Reform image
    img_rms  = SQRT(IMAGE_VARIANCE(img_cur, box_width1, MEAN_IM=img_avg, /NEIGHBOURHOOD))  ; image of standard deviation
    bad_pix1 = WHERE(FINITE(img_rms) EQ 0 OR ABS(img_cur-img_avg)/img_rms GT rms_pix1, n_bad, COMPLEMENT=good_pix)
    IF n_bad GT 0 THEN bad_pix = [bad_pix, bad_pix1]
    
    ; Replace bad pixels and do it a second time
    IF n_bad GT 0 THEN img_cur[bad_pix1] = img_avg[bad_pix1]
    img_rms  = SQRT(IMAGE_VARIANCE(img_cur, box_width1, MEAN_IM=img_avg, /NEIGHBOURHOOD))  ; image of standard deviation 
    bad_pix2 = WHERE(FINITE(img_rms) EQ 0 OR ABS(img_cur-img_avg)/img_rms GT rms_pix1, n_bad, COMPLEMENT=good_pix)
    IF n_bad GT 0 THEN bad_pix = [bad_pix, bad_pix2]
    
    ; Find bad pixels
    IF n_bad GT 0 THEN img_cur[bad_pix2] = img_avg[bad_pix2]
    img_rms  = SQRT(IMAGE_VARIANCE(img_cur, box_width2, MEAN_IM=img_avg, /NEIGHBOURHOOD))  ; image of standard deviation
    bad_pix3 = WHERE(FINITE(img_rms) EQ 0 OR ABS(img_cur-img_avg)/img_rms GT rms_pix2, n_bad, COMPLEMENT=good_pix)
    IF n_bad GT 0 THEN bad_pix = [bad_pix, bad_pix3]

    ; Replace bad pixels and do it a second time
    IF n_bad GT 0 THEN img_cur[bad_pix3] = img_avg[bad_pix3]
    img_rms  = SQRT(IMAGE_VARIANCE(img_cur, box_width2, MEAN_IM=img_avg, /NEIGHBOURHOOD))  ; image of standard deviation
    bad_pix4 = WHERE(FINITE(img_rms) EQ 0 OR ABS(img_cur-img_avg)/img_rms GT rms_pix2, n_bad, COMPLEMENT=good_pix)
    IF n_bad GT 0 THEN bad_pix = [bad_pix, bad_pix4]
     
    ; Find value out of range
    IF KEYWORD_SET(RANGE) THEN BEGIN
      IF n_bad GT 0 THEN img_cur[bad_pix4] = img_avg[bad_pix4]
      bad_pix5  = WHERE(img_cur LT range[0] AND img_cur GT range[1], n_bad)
      IF n_bad GT 0 THEN bad_pix = [bad_pix, bad_pix5]
    ENDIF
    
    ; Remove duplicate
    bad_pix = bad_pix[UNIQ(bad_pix, SORT(bad_pix))]
  ENDFOR
  
  ; Parse bad pixels to output image
  bpm_map[bad_pix] = 0
  
  ; Print info to screen
  IF info GT 1 THEN PRINT, '  Number of cold bad pixels: ', N_ELEMENTS(bad_pix)
  
  RETURN, bpm_map
END

; TEST HARNESS
PRO TEST_BPM
  ; Read images
  test_img = READFITS('/Volumes/LaCie/data/LBTI/nomic/160326/n_160326_002255.fits', /SILENT)
  bpm_img  = READFITS('/Volumes/LaCie/nodrs/results/nomic/bpm/2016-03-26_256x256_bpm1.fits', /SILENT)
  drk_img  = READFITS('/Volumes/LaCie/nodrs/results/nomic/dark/2016-03-26_256x256_dark1.fits', /SILENT)
  flt_img  = READFITS('/Volumes/LaCie/nodrs/results/nomic/flat/2016-03-26_256x256_flat1.fits', /SILENT)
  ; Now flat field the test image
  img = (test_img-drk_img)*flt_img
  FIXPIX, img, bpm_img, img_out, NPIX=8, /SILENT
  ; Plot
  ;LOADCT, 0, /SILENT
  PLOTXY, test_img 
  PLOTXY, flt_img
  PLOTXY, img 
  PRINT, MIN(img), MAX(img), MEAN(img)
  PLOTXY, img_out
  PRINT, MIN(img_out),  MAX(img_out), MEAN(img_out)
END

; TEST HARNESS
PRO TEST_BPM2
  ; Read images
  data_path = '/Volumes/nodrs/data/LBTI/nomic/161114/'
  idx_flt   = [14356,15335]
  bpm_img = LBTI_READDATA(data_path, DATA_IDX=idx_flt, HDR_DATA=hdr_bpm, MLOG_DATA=mlog_data, INFO=info, PLOT=plot, /MEDIAN, /SCALE) ; /SCALE because RANGE below is for one coadd
  bpm_cld = BPM_COLD(bpm_img, INFO=info)
  PLOTXY, bpm_cld, NWINDOW=1
  PLOTXY, bpm_img, NWINDOW=2
  idx_drk   = [1000,1999]
  bpm_img = LBTI_READDATA(data_path, DATA_IDX=idx_drk, HDR_DATA=hdr_bpm, MLOG_DATA=mlog_data, INFO=info, PLOT=plot, /SCALE)
  bpm_hot = BPM_HOT(bpm_img, INFO=info)
  PLOTXY, bpm_hot, NWINDOW=3
END


; Plot if requested
;  IF KEYWORD_SET(plot) THEN BEGIN
;    ; Scale the image for plotting purposes
;    drk_rms[WHERE(drk_rms GT 5.*MEDIAN(drk_rms), /NULL)] = 5.*MEDIAN(drk_rms)
;
;    ; Create PATH if does not exist
;    sep      = pth.sep
;    bpm_path = pth.result_path + sep + 'diagnostic' + sep + 'bpm' + sep + med_hdr[i_med].date_obs + sep   ; Path to the BPM folder
;    rnm_path = pth.result_path + sep + 'diagnostic' + sep + 'rnm' + sep + med_hdr[i_med].date_obs + sep   ; Path to the BPM folder
;    IF NOT FILE_TEST(bpm_path) THEN FILE_MKDIR, bpm_path                                                  ; Create directory if it does not exist
;    IF NOT FILE_TEST(rnm_path) THEN FILE_MKDIR, rnm_path                                                  ; Create directory if it does not exist
;
;    ; Initiate plot
;    fit     = 20./1720.
;    LOADCT, 3, /SILENT
;
;    ; Plot rms image
;    filerms = rnm_path + med_hdr[i_med].date_obs + '_' + STRTRIM(med_hdr[i_med].instrum) + '_DIT-' + STRING(1D+3*med_hdr[i_med].int_time, FORMAT='(I0)') + '_smplm-' + $
;      STRING(med_hdr[i_med].smplmode, FORMAT='(I0)') + '_coad-' + STRING(med_hdr[i_med].n_coadd, FORMAT='(I0)') + $
;      '_pag-' + STRING(med_hdr[i_med].pagain, FORMAT='(I0)') + '_pab-' + STRING(med_hdr[i_med].pabandw, FORMAT='(I0)') + '_bias-' + $
;      STRING(med_11[i_med].detbias, FORMAT='(I0)') + '_detm-' + STRING(med_hdr[i_med].detmode, FORMAT='(I0)') + '_RMS.eps'
;    PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME = filerms
;    PLOTXY, drk_rms,  /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='', XTITLE='pixel', YTITLE='pixel', GRID=0, $
;      CHARSIZE=1.0, CHARTHICK=1.5, THICK=2.5, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=[250,100,1050,900]*fit;, INSET_UR='a.'
;    !X.TITLE = ' '
;    !Y.TITLE = ' '
;    COLORBAR, COLOR=color, DIVISIONS=10, POSITION=[0.80,0.1,0.86,0.9], MAXRANGE=MAX(drk_rms), MINRANGE=MIN(drk_rms), CHARTICK=1.5, FORMAT='(F6.2)', /VERTICAL, /RIGHT
;    LOADCT, 0, /SILENT
;    PLOTXY, /FIN
;
;    ; Compute and plot bad pixel map
;    bpm_img          = 0.*drk_img
;    bpm_img[bad_pix] = 1000.
;    filebmp = bpm_path + med_hdr[i_med].date_obs + '_' + STRTRIM(med_hdr[i_med].instrum) + '_DIT-' + STRING(1D+3*med_hdr[i_med].int_time, FORMAT='(I0)') + '_smplm-' + $
;      STRING(med_hdr[i_med].smplmode, FORMAT='(I0)') + '_coad-' + STRING(med_hdr[i_med].n_coadd, FORMAT='(I0)') + $
;      '_pag-' + STRING(med_hdr[i_med].pagain, FORMAT='(I0)') + '_pab-' + STRING(med_hdr[i_med].pabandw, FORMAT='(I0)') + '_bias-' + $
;      STRING(med_hdr[i_med].detbias, FORMAT='(I0)') + '_detm-' + STRING(med_hdr[i_med].detmode, FORMAT='(I0)') + '_BPM.eps'
;    PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME = filebmp
;    PLOTXY, bpm_img, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='',  XTITLE='pixel', YTITLE='pixel', GRID=0, $
;      CHARSIZE=1.0, CHARTHICK=1.5, THICK=2.5, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=[250,100,1050,900]*fit;, INSET_UR='a.'
;    ;!X.TITLE = ' '
;    ;!Y.TITLE = ' '
;    ;COLORBAR, COLOR=color, DIVISIONS=10, POSITION=[0.82,0.1,0.88,0.9], MAXRANGE=MAX(bpm_img), MINRANGE=MIN(bpm_img), /VERTICAL, /RIGHT
;    LOADCT, 0, /SILENT
;    PLOTXY, /FIN
;  ENDIF