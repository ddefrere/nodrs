;+
; NAME: EXTRACT_PIXELS
; 
; PURPOSE:
;   This routine extracts a set of pixels from an input image. The extracted pixels lie in a region
;   centered at [xin,yin] and within a radius of siz.
;
; INPUTS:
;  image        :   The 2-D input image
;  xin          :   The x coordinates of the center pixel [in pixels]
;  yin          :   The y coordinates of the center pixel [in pixels]
;  siz          :   The radius in pixels of the region to extract around [xin,yin]
;
; KEYWORDS
;   COMPLEMEMT  :   Extract the complement pixels if set
;   PLOT        :   If set, plot the output extracted pixels
;
; OUTPUT
;   The set of pixels located in the region to extract
;
; MODIFICATION HISTORY:
;   Version 1.0,  11-APR-2009, by Denis Defrère, MPIFR, ddefrere@mpifr.de
;   Version 1.1,  11-NOV-2012, DD, added keyword complement

FUNCTION EXTRACT_PIXELS, image, xin, yin, siz, PLOT=plot, COMPLEMENT=complement

; Number of pixels in the input image
n_pix = N_ELEMENTS(image[*,0])

; Coordinates of each pixel
x = -(n_pix-1)/2D0 + DINDGEN(n_pix)/(n_pix-1) * (n_pix-1) & y=x

; Check whether each pixel is in the region to extract. If not, set the pixel to 0.
result = image
IF NOT KEYWORD_SET(complement) THEN BEGIN 
  FOR ix = 0, n_pix-1 DO BEGIN
    FOR iy = 0, n_pix-1 DO BEGIN
      IF SQRT((x[ix]-xin)^2+(y[iy]-yin)^2) GT siz THEN result[ix,iy]=0.
    ENDFOR
  ENDFOR  
ENDIF ELSE BEGIN
  FOR ix = 0, n_pix-1 DO BEGIN
    FOR iy = 0, n_pix-1 DO BEGIN
      IF SQRT((x[ix]-xin)^2+(y[iy]-yin)^2) LE siz THEN result[ix,iy]=0.
    ENDFOR
  ENDFOR  
ENDELSE

; Extract pixel not equal to 0
idx = WHERE(result NE 0, ni)
IF ni GT 0 THEN output = result[idx] ELSE MESSAGE, 'No pixel extracted or input image empty'

; Plot the flux in the asymmetric image (considering only an annulus at 1 AU)
IF KEYWORD_SET(plot) THEN BEGIN
  LOADCT, 38, /SILENT
  fit    = 1.
  PLOTXY, /INIT, /COLOR, /INV, WINDOW=[0, 0, 500, 400]*fit
  PLOTXY, ABS(result), /NEW, XRANGE=[-n_pix/2,n_pix/2], YRANGE=[-n_pix/2,n_pix/2], TITLE='Output pixels', XTITLE='!5X [mas]', YTITLE='!5Y [mas]', GRID=0, $
     CHARSIZE=1.0, CHARTHICK=1.5, THICK=2.5, XSTYLE=1, YSTYLE=1, /COLOR, /NOERASE, WINDOW=[90,50,390,350]*fit;, INSET_UR='a.'
  PLOTXY, /FIN  
ENDIF

RETURN, output
END 