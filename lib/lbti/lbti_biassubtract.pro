;+
; NAME: LBTI_BIASSUBTRACT
; 
; PURPOSE:
;   This procedure performs bias correction for NOMIC and LMIRCAM.
;
; INPUT:
;   img_in        :  An image data cube
;
; KEYWORDS
;   LMIRCAM
;   NOMIC
;
; OUTPUT
;   Data cube with the centered input images
;
; MODIFICATION HISTORY:
;   Version 1.0,  08-JUL-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu

FUNCTION LBTI_BIASSUBTRACT, img_in, LMIRCAM=lmircam, NOMIC=nomic

; Derive size of the input image
size_img = SIZE(img_in)   
n_row    = size_img[1]
n_col    = size_img[2]
IF size_img[0] EQ 3 THEN n_fr = size_img[3] ELSE IF size_img[0] EQ 2 THEN n_fr = 1 ELSE MESSAGE, 'Invalid input image'

; LMIRCAM has a global floating bias level that is integration time dependent, and can change enough over time to
; affect the background level in the short (unsaturated) frames. To correct for this, subtract, in each frame, the
; overscan pixels (bottom two rows and top two rows). Do it per channel not per column. Previous attempts to do it 
; per column were unsuccessful.
IF KEYWORD_SET(LMIRCAM) AND n_row EQ 1024 THEN BEGIN
  ; define bias pixels per column
  idx_bias = [INDGEN(2), INDGEN(2)+n_row-2]
  ; loop over the frames and channels
  chan_pix = 64L
  n_chan   = ROUND(n_col/chan_pix)
  FOR i_fr=0, n_fr-1 DO FOR i_chan = 0, n_chan-1 DO BEGIN
    idx_chan                = i_chan*chan_pix+LINDGEN(chan_pix)
    img_in[idx_chan,*,i_fr] = img_in[idx_chan,*,i_fr] - MEAN(img_in[idx_chan,idx_bias,i_fr])
  ENDFOR
ENDIF

; NOMIC. The reference pixels are those the closest from the edge horizontally (TBC).
;IF STRTRIM(STRUPCASE(instrum)) EQ 'NOMIC' THEN BEGIN
;  n_ref = 1
;  ; Do left side of the array
;  idx_bias = n_ref+INDGEN(0.5*n_x-n_ref)
;  FOR i_fr=0, n_fr-1 DO FOR i_lin = 0, n_y-1 DO img_in[idx_bias,i_lin,i_fr] = img_in[idx_bias,i_lin,i_fr] - MEAN(img_in[INDGEN(n_ref),i_lin,i_fr])
;  ; Do right side of the array
;  idx_bias = 0.5*n_x + INDGEN(0.5*n_x-n_ref)
;  FOR i_fr=0, n_fr-1 DO FOR i_lin = 0, n_y-1 DO img_in[idx_bias,i_lin,i_fr] = img_in[idx_bias,i_lin,i_fr] - MEAN(img_in[n_x-n_ref+INDGEN(n_ref),i_lin,i_fr])
;ENDIF

RETURN, img_in
END