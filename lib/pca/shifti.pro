;+
; NAME:
;	SHIFTI
; PURPOSE:
;	Shift a 1-d or 2-d array by fractional pixels using interpolation
; CALLING SEQUENCE:
;	SHIFTI, IMAGE, XSHIFT, YSHIFT [,/CUBIC]
; INPUTS:
;	IMAGE : Vector or 2-d image array to shift.
;	XSHIFT, YSHIFT : Pixel amounts to shift IMAGE (XSHIFT only for vector).
; OPTIONAL INPUT KEYWORDS:
;	CUBIC : If set, uses cubic instead of bilinear interpolation.
; RETURNS:
;	Image shifted by XSHIFT,YSHIFT amount.
; MODIFICATION HISTORY:
;	Written by John Krist, STScI   Nov 1995
;-

function shifti, image, xshift, yshift, CUBIC=CUBIC, MISSING=MISSING, $
	NOMISSING=NOMISSING, MAG=MAG

on_error,2

if (not keyword_set(MAG)) then MAG=1.0
if (not keyword_set(CUBIC)) then CUBIC = 0

s = size(image)
ndim = s(0)
nx = s(1)
x = (findgen(nx) - xshift) / MAG

if ( ndim eq 2 ) then begin
	ny = s(2)
	y = (findgen(ny) - yshift) / MAG
endif

if (keyword_set(NOMISSING)) then begin
	if ( ndim eq 2 ) then begin
		return, interpolate(image,x,y,/GRID,CUBIC=CUBIC)
	endif else begin
		return, interpolate(image,x,/GRID,CUBIC=CUBIC)
	endelse
endif else if (n_elements(MISSING) eq 0) then begin
	if ( ndim eq 2 ) then begin
		return, interpolate(image,x,y,/GRID,CUBIC=CUBIC,MISSING=0)
	endif else begin
		return, interpolate(image,x,/GRID,CUBIC=CUBIC,MISSING=0)
	endelse
endif else begin
	if ( ndim eq 2 ) then begin
		return, interpolate(image,x,y,/GRID,CUBIC=CUBIC,MISSING=MISSING)
	endif else begin
		return, interpolate(image,x,/GRID,CUBIC=CUBIC,MISSING=MISSING)
	endelse
endelse

end
