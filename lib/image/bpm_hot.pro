;+
; NAME: BPM_HOT
; 
; PURPOSE:
;   Compute hot bad pixel map based on input image cube (dark frames in principle).
;   This code compute the standard deviation of each pixel and flag as bad pixels those with higher than normal read noise (5 sigma threshold).
 ;
; INPUTS:
;   img_in : the input image cube
;
; KEYWORDS
;   RANGE  : in ADU, range of acceptable values [min, max]
;   INFO   : set to print info to screen
;
; OUTPUT
;   The computed bad pixel map (0 is a bad pixel, anything else is good pixel)
;   
; MODIFICATION HISTORY:
;   Version 1.0, 27-APR-2016, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 25-MAY-2016, DD: replaced RESISTANT_MEAN by STDDEV to compute the standard deviation
;   Version 1.2, 23-MAR-2017, DD: replaced MEAN by MEDIAN to compute the average read noise

FUNCTION BPM_HOT, img_in, RANGE=range, INFO=info, PLOT=plot
  
  ; Keyword sanity check
  IF NOT KEYWORD_SET(INFO) THEN info = 0
  
  ; Running parameters
  rms_tre = 5.   ; threshold for bad pixels
  
  ; Number of frames in the input image
  size_img = SIZE(img_in)
  IF size_img[0] EQ 3 THEN n_fr = size_img[3] ELSE BEGIN
    MESSAGE, ' Input image cube should be a cube!', /CONTINUE
    RETURN, 0
  ENDELSE
  
  ; Compute STD per pixel
  sdv = STDDEV(img_in, DIMENSION=3)
   
  ; Compute bad pixels
  ron     = MEDIAN(sdv)
  idx_bad = WHERE(ABS(sdv-ron) GT rms_tre*ron, n_bad)
  
  ; Bad pixel map
  bpm_map = INTARR(size_img[1],size_img[2]) + 1
  IF n_bad GT 1 THEN bpm_map[idx_bad] = 0
   
  ; Print info to screen
  IF info GT 1 THEN PRINT, '  Number of hot bad pixels  : ', n_bad
  
  RETURN, bpm_map
END
