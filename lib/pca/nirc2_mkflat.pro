;+
; NAME:
;	PH_MKFLAT
;
; PURPOSE:
;       Make flat and bad pixel maps from PHARO data
;
; CALLING SEQUENCE:
;       PH_MKFLAT, flat, FRNO=, DATADIR=, PROCDIR=, DATE=
;
; INPUTS:
;	frno=	Vector of frame numbers for one set of flats.  If not given,
;		will prompt for first and last frame numbers.
;
;	datadir= Full path to data directory.  If not specified, determined
;		from date.
;
;	procdir= Full path to processing directory.  If not specified,
;		determined from date.
;
;       date = UT date in 'YYMMDD' format.  User will be prompted for input
;               if datadir and procdir are not defined.
;
; OUTPUTS:
;	flat	flat map (1024x1024 FLOAT array)
;
;	bpix	bad pixel map (1024x1024 BYTE array)
;
; EXAMPLE:
;	IDL> PH_MKFLAT, flat, bpix, fr=436+INDGEN(12), date='060908'
;
; MODIFICATION HISTORY:
;	Original writen 9/10/06, A. Bouchez, Caltech Optical Observatories
;	Cleaned up and updated 5/07 A. Bouchez, COO
;-
PRO 	nirc2_MKFLAT, flat, bpix, frno=frno, datadir=datadir, procdir=procdir

;;; Define constants

  pfmt = '("n",I4.4,".fits")'       ; NIRC2 filename format
  min_gain = 0.7		; Min. allowed gain value (less = bad pixel)
  max_gain = 1.3		; Max. allowed gain value (more = bad pixel)

;;; Query for input parameters

  if not KEYWORD_SET(frno) then begin
    fstr= ''
    READ, 'Please enter first frame number: ', fstr
    lstr = ''
    READ, 'Please enter last frame number: ', lstr
    flno = FIX([fstr, lstr])
    frno = flno[0] + LINDGEN(flno[1]-flno[0]+1)
  endif

;;; Read image stack

  MESSAGE, /info, 'Processing skyflat frames: ' + $
    STRING(MINMAX(frno),f='(I0.0,"-",I0.0)')
  im = rfits(datadir+STRING(frno,f=pfmt), hd_t, /nirc2)
  hd=hd_t[*,0]
;;; Compute flat map

  filt = STRTRIM(SXPAR(hd,'FWINAME'),2)
  
  im_reb = REBIN(im, 1024, 1024)
  flat = (im_reb / MEDIAN(im_reb)) > min_gain < max_gain
  sname = 'flat_' + filt + '.fits'
  SXADDPAR, hd, 'BITPIX', -32
  MESSAGE, /info, '  Writing ' + sname
  WRITEFITS, datadir+sname, flat, hd
 

;;; Compute bad pixel map

  bpix = BYTARR(1024, 1024)
  wbad = WHERE((flat eq min_gain) or (flat eq max_gain))
  if wbad[0] ne -1 then bpix[wbad] = 1
  sname = 'bpix_' + filt + '.fits'
  SXADDPAR, hd, 'BITPIX', 8
  MESSAGE, /info, '  Writing ' + sname
  WRITEFITS, datadir+sname, bpix, hd

  
END
