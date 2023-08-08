;+
; NAME: NOMIC_MASK
; 
; PURPOSE:
;   This procedure computes a 1024x1024 mask for the NOMIC detector with 0. between channels and 1 elsewhere.
;   The width of the 0 sritp is given as input (10. by default).
;
; KEYWORDS
;   STRIP_WIDTH      : in pixels, the width of the 0 strip between channels
;
; OUTPUT:
;   det_mask         : the 1024x1024 output mask
;
; MODIFICATION HISTORY:
;   Version 1.0,  16-JAN-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu

PRO NOMIC_MASK, det_mask, STRIP_WIDTH=strip_width, PLOT=plot

n_pix   = 1024.
IF NOT KEYWORD_SET(strip_width) THEN strip_width = 10.

; Initialize mask
det_mask = REPLICATE(1., n_pix, n_pix)

; Add middle strip
FOR i_pix = 0, strip_width-1 DO BEGIN
  i_str             = 0.5*(n_pix-strip_width) + i_pix
  det_mask[i_str,*] = 0. 
  FOR i_ch = 0, 7 DO BEGIN
    i_str             = i_ch*0.125*n_pix-0.5*strip_width + i_pix
    det_mask[*,i_str] = 0. 
  ENDFOR
ENDFOR
END