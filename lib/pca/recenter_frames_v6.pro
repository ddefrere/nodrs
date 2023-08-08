;------------------------------------------------------------------------;
;Subpixel registration of a datacube of images			         ;
;D. Mawet, Jet Propulsion Laboratory - California Institute of Technology;
;------------------------------------------------------------------------;
;v1: June 2, 2006, D. Mawet
;v2: sometimes in 2008, D. Mawet
;v3: 14 September 2010, D. Mawet
;v4: 31 May 2011, D. Mawet: a mess! 
;v5: cleaned up obsolete methods, make it compatible with IDL 7.0 and IDL 8.0
;Latest revision, v4: 31 May 2011, D. Mawet

;Needs the IDL ASTROLIB;
;IMPORTANT NOTE: assumes window size of 500 x 500

pro recenter_frames_v6, obj, frame_init, cx, cy, fwhm, method, use_cursor=use_cursor

sz = 2*fwhm; characteristic size for centering
sh_limit=2*fwhm; fool proof shift limit
rd = 5*fwhm; mask radius (allows neglecting the center for shift-subtract and correlations)
rw = 18 * fwhm ; display FoV

;Use_cursor, manual pre-centering for wild ill-behaved cases
if keyword_set(use_cursor) then begin
obj_z=obj[cx-rw:cx+rw-1,cy-rw:cy+rw-1]
wset, 2
tvscl, congrid((obj_z),500,500)
print, 'Left click on the approximate center'
cursor, x, y, /device
print, 'x, y', x , y 
shx=-rw+float(x)/500*rw*2
shy=-rw+float(y)/500*rw*2
print, 'shx, shy', shx, shy
obj=shifti(obj, -shx,-shy, cubic=-0.5)
endif

CASE method OF

	0: BEGIN 
	obj_tmp=obj
	;Method 0: centroiding CNTRD (derivative search)
	;-----------------------------------------------
	print, 'Using cntrd method'
	cntrd, obj_tmp, cx, cy, xcen, ycen, sz;, /KeepCenter
	shx=cx-xcen & shy=cy-ycen
	print, 'Target Offset (x,y)', shx, shy
	END
	
	1: BEGIN 
	obj_tmp=obj
	;Method 1: centroiding GCNTRD (Gaussian fits to marginal X,Y)
	;------------------------------------------------------------
	print, 'Using gcntrd method'
	gcntrd, obj_tmp, cx, cy, xcen, ycen, sz;, /KeepCenter
	shx=cx-xcen & shy=cy-ycen
	print, 'Target Offset (x,y)', shx, shy
	END

	2: BEGIN
	obj_tmp=obj
	;Method 2: Simple centroiding (center of gravity)
	;------------------------------------------------
	print, 'Using CG centroiding'
	obj_tmp=obj_tmp[cx-sz/2:cx+sz/2-1,cy-sz/2:cy+sz/2-1]
	cm=centroid(obj_tmp)
	shx=sz/2-cm[0] & shy=sz/2-cm[1]
	print, 'Target Offset (x,y)', shx, shy
	END

	3: BEGIN
	obj_tmp=obj
	;Method 3: Shift-subtract (suggested by J. Krist)
	;------------------------------------------------
	print, 'Using shift subtract method'
	range = 4 ;(x pixels: -x/2:+x/2)
	mag = range * 10; mag (10 or 20 or 40)
	if rd eq 0 then mask=1.0 else mask=(shift(dist(rw),rw/2,rw/2) ge rd)
	obj_tmp=obj_tmp[cx-rw/2:cx+rw/2-1,cy-rw/2:cy+rw/2-1]*mask
	frame_init_tmp=frame_init[cx-rw/2:cx+rw/2-1,cy-rw/2:cy+rw/2-1]*mask
	magobj=fltarr(mag,mag)
	shx_t=(findgen(mag)-mag/2)/mag * range
	shy_t=(findgen(mag)-mag/2)/mag * range
		for x=0,mag-1 do begin
			for y=0,mag-1 do begin
				if shx_t[x] eq 0 then shx_t[x]=1e-4
				if shy_t[y] eq 0 then shy_t[y]=1e-4
				obj_tmp_sh=shifti(obj_tmp,shx_t[x],shy_t[y],cubic=-0.5)
				magobj[x,y]=total(abs(obj_tmp_sh-frame_init_tmp)^2d);/(2*sz)^2
			endfor
		endfor
	wset, 2
	tvscl, congrid(magobj,500,500)
	tmp=min(magobj,location)
	ind = ARRAY_INDICES(magobj, location)
	shx = shx_t[ind[0]] & shy = shy_t[ind[1]]
	print, 'Target Offset (x,y)', shx, shy
	END
		
	4: BEGIN
	obj_tmp=obj
	;Method 4: correlation
	;---------------------
	print, 'Using image correlation'
	mag=4; Enter correlation magnification
	if rd eq 0 then mask=1.0 else mask=(shift(dist(rw),rw/2,rw/2) ge rd)
	obj_tmp=obj_tmp[cx-rw/2:cx+rw/2-1,cy-rw/2:cy+rw/2-1]*mask
	frame_init_tmp=frame_init[cx-rw/2:cx+rw/2-1,cy-rw/2:cy+rw/2-1]*mask
	wset, 2
	tvscl, congrid(obj_tmp,500,500)
	wset, 3
	tvscl, congrid(frame_init_tmp,500,500)
	correl_optimize, frame_init_tmp, obj_tmp, shx, shy, MAGNIFICATION=mag, /NUMPIX
	END

	5: BEGIN
	;Method 5: NO recentering
	;------------------------
	print, 'No centering'
	END
ENDCASE


;Shift double-Checking and subpixel shifting
;-------------------------------------------
if abs(shx) gt sh_limit or abs(shy) gt sh_limit then begin
	print, 'Shift outside valide range'
	shx = 0 & shy = 0
endif else begin
	obj=shifti(obj, shx,shy, cubic=-0.5)
endelse
wset,1
tvscl, congrid((obj[cx-rw:cx+rw,cy-rw:cy+rw]^0.25),500,500)>0<10000

END
