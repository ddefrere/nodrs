PRO make_ff_nirc2, nirc2_filter, datadir, procdir, flat=flat, bpix=bpix

	;Make Flat field
	;---------------
	if nirc2_filter eq 'H' then begin 
		if (file_test(datadir+'flat_H.fits') eq 0) then begin	
		print, 'H FLAT and BAD PIXEL maps not detected in current directory, please enter starting and ending frame #:'
		read, st, PROMPT='Enter starting frame # for FLATs: '
		read, nd, PROMPT='Enter ending frame # for FLATs: '
		  ;st = 19			; starting FF frame number
		  ;nd = 36			; ending FF frame number
		  fr = st + findgen(nd-st+1)
		  nirc2_mkflat, flat, bpix, frno=fr, datadir=datadir, procdir=procdir
		endif else begin
		 ;flat = readfits(datadir+'flat_Kp.fits')
		 ;bpix = readfits(datadir+'bpix_Kp.fits')
		  flat = readfits(datadir+'flat_H.fits')
		  bpix = readfits(datadir+'bpix_H.fits')
		endelse
	endif
	if nirc2_filter eq 'Kp' then begin 
		if (file_test(datadir+'flat_Kp.fits') eq 0) then begin	
		print, 'Kp FLAT and BAD PIXEL maps not detected in current directory, please enter starting and ending frame #:'
		read, st, PROMPT='Enter starting frame # for FLATs: '
		read, nd, PROMPT='Enter ending frame # for FLATs: '
		  ;st = 19			; starting FF frame number
		  ;nd = 36			; ending FF frame number
		  fr = st + findgen(nd-st+1)
		  nirc2_mkflat, flat, bpix, frno=fr, datadir=datadir, procdir=procdir
		endif else begin
		  flat = readfits(datadir+'flat_Kp.fits')
		  bpix = readfits(datadir+'bpix_Kp.fits')
		endelse
	endif
	
	if nirc2_filter eq 'Ks' then begin 
		if (file_test(datadir+'flat_Ks.fits') eq 0) then begin	
		print, 'Ks FLAT and BAD PIXEL maps not detected in current directory, please enter starting and ending frame #:'
		read, st, PROMPT='Enter starting frame # for FLATs: '
		read, nd, PROMPT='Enter ending frame # for FLATs: '
		  ;st = 19			; starting FF frame number
		  ;nd = 36			; ending FF frame number
		  fr = st + findgen(nd-st+1)
		  nirc2_mkflat, flat, bpix, frno=fr, datadir=datadir, procdir=procdir
		endif else begin
		  flat = readfits(datadir+'flat_Ks.fits')
		  bpix = readfits(datadir+'bpix_Ks.fits')
		endelse
	endif

	if nirc2_filter eq 'K' then begin 
		if (file_test(datadir+'flat_K.fits') eq 0) then begin	
		print, 'K FLAT and BAD PIXEL maps not detected in current directory, please enter starting and ending frame #:'
		read, st, PROMPT='Enter starting frame # for FLATs: '
		read, nd, PROMPT='Enter ending frame # for FLATs: '
		  ;st = 19			; starting FF frame number
		  ;nd = 36			; ending FF frame number
		  fr = st + findgen(nd-st+1)
		  nirc2_mkflat, flat, bpix, frno=fr, datadir=datadir, procdir=procdir
		endif else begin
		  flat = readfits(datadir+'flat_K.fits')
		  bpix = readfits(datadir+'bpix_K.fits')
		endelse
	endif
	
	if nirc2_filter eq 'Lp' then begin 
		if (file_test(datadir+'flat_Lp.fits') eq 0) then begin	
		print, 'Lp FLAT and BAD PIXEL maps not detected in current directory, please enter starting and ending frame #:'
		read, st, PROMPT='Enter starting frame # for FLATs: '
		read, nd, PROMPT='Enter ending frame # for FLATs: '
		  ;st = 19			; starting FF frame number
		  ;nd = 36			; ending FF frame number
		  fr = st + findgen(nd-st+1)
		  nirc2_mkflat, flat, bpix, frno=fr, datadir=datadir, procdir=procdir
		endif else begin
		  flat = readfits(datadir+'flat_Lp.fits')
		  bpix = readfits(datadir+'bpix_Lp.fits')
		endelse
	endif
	
END