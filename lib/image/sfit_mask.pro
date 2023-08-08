;+
; NAME: SFIT_MASK
;
; PURPOSE:
;   Fit a surface to the input image, excluding a circular region centered in the middle of the image. 
;   The region exlcuded from the fit is then interpolated using the best-fit coefficients. 
;
; INPUTS:
;   img       :  the input image to fit
;   mask      :  the diameter of the mask (in pixels)
;
; KEYWORDS
;   DEGREE    :  the degree of the polyniomal fit (3 by default)
;
; OUTPUT
;   The best-fit surface
;
; MODIFICATION HISTORY:
;   Version 1.0, 12-APR-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu 
   
FUNCTION SFIT_MASK, img, mask, DEGREE=degree
  
  ; Keyword consistancy check
  IF NOT KEYWORD_SET(DEGREE) THEN degree = 3
  
  ; Read array dimension
  nx = N_ELEMENTS(img[*,0])
  ny = N_ELEMENTS(img[0,*])
  n  = nx*ny
  
  ; Create the corresponding position arrays
  x = LINDGEN(nx)#REPLICATE(1,ny)
  y = REPLICATE(1,nx)#LINDGEN(ny)
  
  ; Create and parse the output array
  data = DBLARR(3,n)
  data[0,*] = REFORM(x, n, 1)
  data[1,*] = REFORM(y, n, 1)
  data[2,*] = REFORM(img, n, 1)  
  
  ; Compute mask
  DIST_CIRCLE, d, [nx,ny]
  ; Extract part to be fitted
  d_vec   = REFORM(d, n, 1)  
  idx_ok  = WHERE(d_vec GT mask)
  data    = data[*,idx_ok]
  ; Fit the surface
  img_fit = SFIT(data, degree, /IRREGULAR, KX=kx)
  ; Reconstruct the part in the mask
  img_fit = DBLARR(nx,ny)
  FOR i=0,degree DO FOR j=0,degree DO img_fit += kx[j,i] * x^i * y^j
  
RETURN, img_fit
END