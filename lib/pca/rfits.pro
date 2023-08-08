;+
; NAME:
;	RFITS
;
; PURPOSE:
;	Read multiple FITS files into an image cube.
;
; EXPLANATION:
;	PHARO files are saved as either 512x512x4 (difference) or 512x512x8
;	(endpoint read) cubes.  In the case of endpoint read cubes, will
;	subtract first read from second, returning the difference.  This 
;	program is intended to replace READFITS.pro for PHARO images.
;
; CALLING SEQUENCE:
;	Result = RFITS(Filenames, Headers, subfr=, exten_no=, /pharo, /silent)
;
; INPUTS:
;       Filenames = String vector containing the names of the FITS files
;		to be read.   If the filename has a ".gz" extension, it
;		will be treated as a gzip compressed file.   If it has
;		a ".Z" extension, it will be treated as a Unix compressed file.
;
; OUTPUTS:
;       Result = FITS data array constructed from designated record.
;                If the specified file was not found, then Result = -1
;
; OPTIONAL OUTPUT:
;       Headers = String array containing the headers from the FITS file.
;
; OPTIONAL INPUT KEYWORDS:
;	SUBFR = Vector defining the subframe to extract from each FITS image.
;		defined in pixels as [x0, x1, y0, y1].
;
;       EXTEN_NO - non-negative scalar integer specifying the FITS extension to
;               read.  For example, specify EXTEN = 1 or /EXTEN to read the
;               first FITS extension.
;
;       /PHARO - Keyword passed to RFITS, which calls READPHARO instead of 
;               READFITS to correctly interpret 512x512x4 PHARO image cubes.
;
;       /SILENT - Normally, READFITS will display the size the array at the
;                 terminal.  The SILENT keyword will suppress this.
;
; PROCEDURES USED:
;	Functions: READFITS(), READPHARO()
;
; MODIFICATION HISTORY:
;	Create 04/04 by A. Bouchez, W.M. Keck Observatory
;	Added /PHARO keyword 06/06, A. Bouchez, Caltech Optical Observatories
;	Modified to use READPHARO 05/07, A. Bouchez, COO
;-
FUNCTION	RFITS, filenames, headers, subfr=subfr, $
  exten_no=exten_no, pharo=pharo, naco=naco, nici=nici, nirc2=nirc2,silent=silent

;;; Read single OR multiple fits files, into an array of images and headers.
;;; Create AHB 04/04

  hd0 = HEADFITS(filenames[0], exten=0, /silent)
subfr=0
  if KEYWORD_SET(naco) then begin
  hd0 = HEADFITS(filenames[0], exten=0, /silent)
  sz = [SXPAR(hd0,'NAXIS1'), SXPAR(hd0,'NAXIS2'), N_ELEMENTS(filenames)]
  itype = SXPAR(hd0,'BITPIX')
  endif

  if KEYWORD_SET(nici) then begin
  hd1 = HEADFITS(filenames[0], exten=1, /silent)
  sz = [SXPAR(hd1,'NAXIS1'), SXPAR(hd1,'NAXIS2'), N_ELEMENTS(filenames)]
  itype = SXPAR(hd1,'BITPIX')	
  endif 
 
  if KEYWORD_SET(nirc2) then begin
  hd0 = HEADFITS(filenames[0], exten=0, /silent)
  sz = [SXPAR(hd0,'NAXIS1'), SXPAR(hd0,'NAXIS2'), N_ELEMENTS(filenames)]
  itype = SXPAR(hd0,'BITPIX')
  endif
  
  if KEYWORD_SET(pharo) then begin
    hd0 = HEADFITS(filenames[0], exten=0, /silent)
    sz[0:1] = sz[0:1] * 2
    if KEYWORD_SET(subfr) then begin
    sz[0:1] = [subfr[1]-subfr[0]+1, subfr[3]-subfr[2]+1] 
    endif else begin 
    subfr = [0, sz[0]-1, 0, sz[1]-1]
  endelse
  endif
  
  case itype of
      8: im = BYTARR(sz[0],sz[1],sz[2],/nozero)
     16: im = INTARR(sz[0],sz[1],sz[2],/nozero)
     32: im = LONARR(sz[0],sz[1],sz[2],/nozero)
     64: im = LON64ARR(sz[0],sz[1],sz[2],/nozero)
    -32: im = FLTARR(sz[0],sz[1],sz[2],/nozero)
    -64: im = DBLARR(sz[0],sz[1],sz[2],/nozero)
  endcase
  headers = STRARR(N_ELEMENTS(hd0), sz[2])

  for n=0,sz[2]-1 do begin
    if KEYWORD_SET(pharo) then begin 
      dat = READPHARO(filenames[n], hd, exten_no=exten_no, silent=silent) 
      im[*,*,n] = dat[subfr[0]:subfr[1],subfr[2]:subfr[3]]
    endif else begin
      dat = READFITS(filenames[n], hd, exten_no=exten_no, silent=silent, /noscale)
      im[*,*,n] = dat;[subfr[0]:subfr[1],subfr[2]:subfr[3]]
    endelse

    if KEYWORD_SET(exten_no) then $
      hd = HEADFITS(filenames[n], /silent)
    if N_ELEMENTS(tmp) ge N_ELEMENTS(hd0) then $
      headers[*,n] = hd[0:N_ELEMENTS(hd0)-1] else $
      headers[0:N_ELEMENTS(hd)-1,n] = hd
  endfor

RETURN,im
END
