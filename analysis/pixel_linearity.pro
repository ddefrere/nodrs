; +
; NAME: PIXEL_LINEARITY
;   
; PURPOSE:
;   Basic routine to plot a linearity curve.
;
; INPUTS:
;   date : the date to analyze (e.g., '150208')
;
; OPTIONAL INPUT KEYWORDS
;   DARK_IDX   :  Two-element vector with the lower and the upper file numbers of the dark files (e.g., [0,9])
;   DATA_IDX   :  Two-element vector with the lower and the upper file numbers of the data files (e.g., [20,999])
;   FLAT_IDX   :  Two-element vector with the lower and the upper file numbers of the flat files (e.g., [10,19])
;   LMIRCAM    :  Set if LMIRCam data
;   RANDOM     :  Set to pick random pixels (look in a box near 64,64 by default)
;   
; MODIFICATION HISTORY:
;   Version 1.0,  23-JAN-2015, Denis Defrère, University of Arizona

PRO PIXEL_LINEARITY, date, DATA_IDX=data_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, LMIRCAM=lmircam, RANDOM=random

; Recover the IDL running path
IF KEYWORD_SET(LMIRCAM) THEN instrum = 'lmircam' ELSE instrum = 'nomic'
DECLARE_PATH, pth, INSTRUM=instrum

; Read data
data_path = pth.root_data + date + pth.sep
img_data = LBTI_READDATA(data_path, DATA_IDX=data_idx, HDR_DATA=hdr_data, INFO=info, PLOT=plot)

; Assign a different config ID for each integration time
hdr_data.cfg_id = FIX(1D3*hdr_data.int_time)

; Find the science frames
idx_drk = WHERE(hdr_data.datatype EQ 1, n_drk, COMPLEMENT=idx_sci)
IF MIN(idx_sci) NE -1 THEN n_sci = N_ELEMENTS(idx_sci) ELSE MESSAGE, 'No science frames to calibrate!'

; Create output header and image array
hdr_sci = TEMPORARY(hdr_data[idx_sci])
img_sci = TEMPORARY(img_data[*,*,idx_sci])

; Create a median image per integration time
img_med = LBTI_IMGMED(img_sci, HDR_DATA=hdr_sci) 
img_cor = img_med
n_xpix  = N_ELEMENTS(img_med[*,0,0])
n_ypix  = N_ELEMENTS(img_med[0,*,0])
n_med   = N_ELEMENTS(img_med[0,0,*])

; Create output header and image array
IF n_drk GT 1 OR KEYWORD_SET(DARK_IDX) THEN BEGIN
  IF NOT KEYWORD_SET(DARK_IDX) THEN BEGIN
    hdr_drk = TEMPORARY(hdr_data[idx_drk])
    img_drk = TEMPORARY(img_data[*,*,idx_drk])
  ENDIF ELSE img_drk = LBTI_READDATA(data_path, DATA_IDX=dark_idx, HDR_DATA=hdr_drk)
  hdr_drk.cfg_id = FIX(1D3*hdr_drk.int_time)      ; Assign a different config ID for each integration time
  med_drk = LBTI_IMGMED(img_drk, HDR_DATA=hdr_drk) 
  ; Compute correction
  FOR i_med = 0, n_med-1 DO img_cor[0,0,i_med] = img_med[*,*,i_med]-med_drk[*,*,i_med]
ENDIF 

; Perform flat-fielding if FLAT_IDX is set
IF KEYWORD_SET(FLAT_IDX) THEN BEGIN
  ; Compute flat image
  img_flt = LBTI_READDATA(data_path, DATA_IDX=flat_idx, HDR_DATA=hdr_flt)
  hdr_flt.cfg_id = FIX(1D3*hdr_flt.int_time)      ; Assign a different config ID for each integration time
  med_flt = LBTI_IMGMED(img_flt, HDR_DATA=hdr_flt)
  flt_cor = LBTI_IMGFLAT(med_flt, HDR_FLT=hdr_flt, IMG_DRK=med_drk, HDR_DRK=hdr_drk)
  ; Compute correction
  FOR i_med = 0, n_med-1 DO img_cor[0,0,i_med] = img_cor[*,*,i_med]*flt_cor
ENDIF

; Create output folder and plot the linearity curves
lin_path = pth.result_path + pth.sep + 'diagnostic' + pth.sep + 'linearity'           ; Path to the linearity folder
IF NOT FILE_TEST(lin_path) THEN FILE_MKDIR, lin_path                          ; Create directory if it does not exist

; Prepare plot
fit     = 20./1720.
LOADCT, 3, /SILENT
xrange = [MIN(hdr_data.int_time),MAX(hdr_data.int_time)]*1D3
IF KEYWORD_SET(LMIRCAM) THEN yrange = [0.,16000.] ELSE yrange = [0.,20000.]
n_r    = 100
n_xpix = N_ELEMENTS(img_data[*,0,0])
n_ypix = N_ELEMENTS(img_data[0,*,0])
 
; Now plot
PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME = lin_path + pth.sep + hdr_data[0].date_obs + '_PIX_LINEAR.eps'
PLOTXY, [0], [0], /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='', XTITLE='Integration time [ms]', YTITLE='Flux [ADU]', GRID=0, $
        CHARSIZE=1.0, CHARTHICK=1.5, THICK=2.5, XSTYLE=1, YSTYLE=1, XTHICK=4, YTHICK=4, /NOERASE, WINDOW=[250,100,1350,900]*fit;, INSET_UR='a.'
FOR i_r = 0, n_r-1 DO BEGIN
  IF KEYWORD_SET(RANDOM) THEN BEGIN
    rand_x = ROUND(RANDOMU(seed)*n_xpix)-1
    rand_y = ROUND(RANDOMU(seed)*n_ypix)-1
  ENDIF ELSE BEGIN
    rand_x = 59 + (i_r MOD 10)
    rand_y = 59 + i_r / 10 
  ENDELSE
  PLOTXY, [hdr_sci.int_time*1D3], [REFORM(img_cor[rand_x,rand_y,*])], LINESTYLE=0, /ADD
  PRINT, i_r, [hdr_sci.int_time*1D3], [REFORM(img_cor[rand_x,rand_y,*])]
ENDFOR
;PLOTXY, [0,1D5], [12750,12750], LINESTYLE=1, /ADD
LOADCT, 0, /SILENT
PLOTXY, /FIN

; Plot median curve
edge    = 10 ; edge pixels to ignore
flx_med = DBLARR(n_med)
drk_med = DBLARR(n_med)
FOR i = 0, n_med-1 DO BEGIN
  flx_med[i] = MEDIAN(img_med[edge:(n_xpix-edge),edge:(n_ypix-edge),i])
  drk_med[i] = MEDIAN(med_drk[edge:(n_xpix-edge),edge:(n_ypix-edge),i])
ENDFOR
PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME = lin_path + pth.sep + hdr_data[0].date_obs + '_MED_LINEAR.eps'
PLOTXY, [0], [0], /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='', XTITLE='Integration time [ms]', YTITLE='Flux [ADU]', GRID=0, $
  CHARSIZE=1.0, CHARTHICK=1.5, THICK=2.5, XTHICK=4, YTHICK=4, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=[250,100,1350,900]*fit;, INSET_UR='a.'
PLOTXY, [hdr_sci.int_time*1D3], flx_med, LINESTYLE=0, THICK=4, /ADD
PLOTXY, [hdr_sci.int_time*1D3], drk_med, LINESTYLE=1, THICK=4, /ADD
LOADCT, 0, /SILENT
PLOTXY, /FIN
END