; +
; NAME: CORRECT_BIAS
;
; PURPOSE:
;   This procedure proceeds in an iterative correction of the high-frequency spatial noise seen on NOMIC.
;   It subtracts from each row its mean, then do the same with the columns, and repeat twice by default.
;
; INPUTS:
;   img_in     :  2-dimension array with the image to process
;
; KEYWORDS:
;   EXCL_WIDTH : Width of the region to exclude around the star (in pixels, see star_pos)
;   LIMIT      : TWo-element vector with the row/column limits to use
;   N_ITER     : Number of iterations
;   STAR_POS   : Position of the star in img_in (assumed at the center if not set)
;
; MODIFICATION HISTORY:
;   Version 1.0, 18-SEP-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu

FUNCTION CORRECT_BIAS, img_in, EXCL_WIDTH=excl_width, LIMIT=limit, N_ITER=n_iter, STAR_POS=star_pos

; Input image info
n_xpix = N_ELEMENTS(img_in[*,0])
n_ypix = N_ELEMENTS(img_in[0,*])

; Keyword sanity check
IF NOT KEYWORD_SET(LIMIT)  THEN limit  = [0, n_xpix-1]
IF NOT KEYWORD_SET(N_ITER) THEN n_iter = 2

; Star position
IF NOT KEYWORD_SET(STAR_POS) THEN BEGIN
  dx = 0.5*n_xpix
  dy = 0.5*n_ypix
ENDIF ELSE BEGIN
  dx = star_pos[0]
  dy = star_pos[1]
ENDELSE

; Compute distance to star
dxsq  = (FINDGEN(n_xpix)-dx)^2
rsq   = FLTARR(n_xpix, n_ypix, /NOZERO)
for ii = 0, n_ypix-1 do rsq[0,ii] = dxsq + (ii-dy)^2

; Distance to bottom left pixel amd from star
x = LINDGEN(n_xpix) 
y = LINDGEN(n_ypix) 
r = sqrt(rsq) - 0.5

; Begin
FOR i_iter = 0, n_iter-1 DO BEGIN
  ; Now remove median of each row to each pixel of the same row
  FOR ir = 0, n_ypix-1 DO BEGIN
    idx_row = WHERE(x GT dx-limit[0] AND x LT dx+limit[1] AND r[*,ir] GE 0.5*excl_width)
    img_in[*,ir] -= MEAN(img_in[idx_row,ir])
  ENDFOR
  
  ; Do the same for the columns
  FOR ic = 0, n_xpix-1 DO BEGIN
    idx_col = WHERE(y GT dy-limit[2] AND y LT dy+limit[3] AND r[ic,*] GE 0.5*excl_width)
    img_in[ic,*] -= MEAN(img_in[ic,idx_col])
  ENDFOR
ENDFOR

RETURN, img_in
END