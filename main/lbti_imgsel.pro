;+
; NAME: LBTI_IMGSEL
; 
; PURPOSE:
;   This procedure selects and coadds the images in the input file
;
; INPUTS:
;   file        : Input FITS file with the image cube to process
;   data_file   : Data file with relevant information about the frames in file (in the first extension of file by default)
;
; KEYWORDS:
;   MEDIAN      : If set, the median of the selected image cube is removed from each image
;   INFO        : Define the level of information printed to screen (see main routine for more info)
;   PLOT        : Set this keyword to plot the data
;
; MODIFICATION HISTORY:
;   Version 1.0, 02-FEB-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (based on LBTI_MERGEL1IMG.pro)
;   Version 1.1, 23-APR-2015: DD, added background selection

PRO LBTI_IMGSEL, file, data_file, MEDIAN=median, INFO=info, PLOT=plot

COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs, log

; Read the image file 
img_in = MRDFITS(file, 0, hdr_in, RANGE=range, /SILENT)
size   = SIZE(img_in)
IF size[0] LT 3 THEN MESSAGE, 'Input image cube should at least have 2 images'
n_xpix = size[1]
n_ypix = size[2]
n_img  = size[3]

; Read header information
pixsize = FXPAR(hdr_in, 'PIXSCALE')

; Read the data file  (or look in the first extension for the data)
IF KEYWORD_SET(DATA_FILE) THEN data_in = MRDFITS(data_file, 1, RANGE=range, /SILENT) ELSE data_in = MRDFITS(file, 1, RANGE=range, /SILENT)
  
; Sanity checks
ncoadd = drs.n_coadd > 1
IF n_img/ncoadd LT 2 THEN ncoadd = 1 ;PRINT, 'WARNING: Not enough images for the specifcied number of coadds'

; Perform background selection
IF N_ELEMENTS(drs.bck_lim) GT 1 THEN BEGIN
  ; Compute flux in a 30x30 region around the beam center
  bckg = DBLARR(n_img)
  FOR i_img = 0, n_img-1 DO bckg[i_img] = MEAN(img_in[*,*,i_img])
  ; Compute good range
  idx_bck = WHERE(bckg GE drs.bck_lim[0] AND bckg LE drs.bck_lim[1], n_bck)
  ; Extract the good data
  img_in  = img_in[*,*,idx_bck]
  data_in = data_in[idx_bck]
  PRINT, 'Number of rejected frames based on background selection : ', n_img-n_bck
  n_img = n_bck
ENDIF ELSE n_bck = 0

; Perform flux selection
IF drs.sig_flx THEN BEGIN
  ; Compute flux in a 30x30 region around the beam center
  n_cen   = 30
  flux    = DBLARR(n_img)
  idx_cen = 0.5*size[1]-0.5*n_cen + INDGEN(n_cen)
  FOR i_img = 0, n_img-1 DO flux[i_img] = TOTAL(img_in[idx_cen,idx_cen,i_img])
  ; Compute mean and rms flux
  AVGSDV, flux, avg, rms
  idx_flx = WHERE(flux GE avg-drs.sig_flx*rms, n_flx)
  ; Extract the good data
  img_in  = img_in[*,*,idx_flx]
  data_in = data_in[idx_flx]
  PRINT, 'Number of rejected frames based on flux selection : ', n_img-n_flx
  n_img = n_flx
ENDIF ELSE n_flx = 0

; Keep only the highest visibility images
IF drs.sig_vis THEN BEGIN
  visi    = DBLARR(n_img)
  FOR i_img = 0, n_img-1 DO BEGIN
    idx_max     = WHERE(REFORM(img_in[*,*,i_img]) EQ MAX(img_in[*,*,i_img]))
    idx         = ARRAY_INDICES(REFORM(img_in[*,*,i_img]), idx_max)
    x_pix       = INDGEN(6)-3+idx[0]
    y_pix       = INDGEN(6)-3+idx[1]
    visi[i_img] = (MAX(img_in[x_pix,y_pix,i_img])-MIN(img_in[x_pix,y_pix,i_img]))/(MAX(img_in[x_pix,y_pix,i_img])+MIN(img_in[x_pix,y_pix,i_img]))
  ENDFOR
  AVGSDV, visi, avg, rms
  idx_vis = WHERE(visi GE MEDIAN(visi)+drs.sig_vis*rms, n_visi)
  PRINT, 'Median visibility :', MEDIAN(visi[idx_vis])
  ; Extract the good data
  img_in  = img_in[*,*,idx_vis]
  data_in = data_in[idx_vis]
  PRINT, 'Number of rejected frames based on visibility selection : ', n_img-n_visi
  n_img = n_visi
ENDIF ELSE n_flx = 0

; Keep only the best centered frames
IF drs.sig_pos THEN BEGIN
  dst     = SQRT((data_in.xcen-MEAN(data_in.xcen))^2+(data_in.ycen-MEAN(data_in.ycen))^2)
  AVGSDV, dst, avg_dst, rms_dst
  idx_dst = WHERE(ABS(dst-avg_dst) LE drs.sig_pos*rms_dst, n_dst)
  img_in  = img_in[*,*,idx_dst]
  data_in = data_in[idx_dst]
  PRINT, 'Number of rejected frames based on stellar position : ', n_img-n_dst
  n_img = n_dst
ENDIF ELSE n_dst = 0

; Keep only the best slopes
IF drs.sig_slo THEN BEGIN
  slope = data_in.slope
  AVGSDV, slope, avg_slope, rms_slope
  idx_slope = WHERE(ABS(slope-avg_slope) LE drs.sig_slo*rms_slope, n_slope)
  IF n_slope NE n_img THEN BEGIN
    img_in  = img_in[*,*,idx_slope]
    data_in = data_in[idx_slope]
  ENDIF
  PRINT, 'Number of rejected frames based on Moffat slope : ', n_img-n_slope
  n_img = n_slope
ENDIF ELSE n_slope = 0

; Keep only the frames which provides the highest correlation (only work for synthetic PSF now)
shft_x = INTARR(n_img)
shft_y = INTARR(n_img)
n_keep = 0
IF drs.psf_file EQ 'synthetic' THEN BEGIN
  IF drs.r_keep LT 1 AND drs.r_keep GT 0 THEN BEGIN
    ; Create and save synthetic PSF image
    wav     = FXPAR(hdr_in, 'WAVELENG')
    bdwth   = FXPAR(hdr_in, 'BANDWIDT')
    img_psf = FIZEAU_PSF(cnf.tel_diam, wav, cnf.base, BANDWIDTH=bdwth, CEN_OBS=cnf.cen_obs, PIX_SIZE=1D3*pixsize, VERBOSE=verbose)
    p=STRPOS(file,'.fits')
    IF p GT 0 THEN psf_file=STRMID(file, 0, p) + '_PSF.fits' ELSE MESSAGE, 'Invalid file format.'
    MKHDR, hdr_psf, img_psf
    MWRFITS, img_psf, psf_file, hdr_psf, /CREATE, /SILENT
    ; Extract a sub-image around beam position and compute correlation
    n_ext   = N_ELEMENTS(img_psf[*,0])
    corel   = DBLARR(n_img)
    pos_max = DBLARR(n_img,2)
    FOR i_img = 0, n_img-1 DO BEGIN
      img_tmp      = EXTRAC(img_in[*,*,i_img], 0.5*n_xpix-0.5*n_ext,0.5*n_ypix-0.5*n_ext, n_ext, n_ext)  ; image should be centered here
      corel_map    = CORREL_IMAGES(img_tmp/MAX(img_tmp), img_psf, XSHIFT=5, YSHIFT=5)
      corel[i_img] = MAX(corel_map, idx_max) > 1D-10   ; Make sure correlation is positive!
      pos_max[i_img,*] = ARRAY_INDICES(corel_map, idx_max)
    ENDFOR
    ; Extract best x% 
    n_keep  = FLOOR(drs.r_keep*n_img)
    idx_cor = N_SMALLEST(1./corel, n_keep)
    idx_cor = idx_cor[SORT(idx_cor)]
    img_in  = TEMPORARY(img_in[*,*,idx_cor])
    data_in = TEMPORARY(data_in[idx_cor])
    pos_max = pos_max[idx_cor,*]
    ; Derive mean position of the best correlation
    shft_x = ROUND(0.5*MEAN(pos_max[*,0]));-1
    shft_y = ROUND(0.5*MEAN(pos_max[*,1]));-1
    PRINT, 'Number of rejected frames based on correlation selection : ', n_img-n_keep
    n_img = n_keep
  ENDIF
ENDIF

; Keep only a given X position on the detector
IF KEYWORD_SET(drs.X_RANGE) THEN BEGIN
  xcen = data_in.xcen
  idx_range = WHERE(xcen GT drs.x_range[0] AND xcen LE drs.x_range[1], n_range)
  IF n_range NE n_img THEN BEGIN
    img_in  = TEMPORARY(img_in[*,*,idx_range])
    data_in = TEMPORARY(data_in[idx_range])
  ENDIF
  PRINT, 'Number of rejected frames based on X position : ', n_img-n_range
  n_img = n_range
ENDIF

; Keep only a given Y position on the detector
IF N_ELEMENTS(drs.Y_RANGE) EQ 2 THEN BEGIN
  ycen = data_in.ycen
  idx_range = WHERE(ycen GT drs.y_range[0] AND ycen LE drs.y_range[1], n_range)
  IF n_range NE n_img THEN BEGIN
    img_in  = TEMPORARY(img_in[*,*,idx_range])
    data_in = TEMPORARY(data_in[idx_range])
  ENDIF
  PRINT, 'Number of rejected frames based on Y position : ', n_img-n_range
  n_img = n_range
ENDIF

; Interpol the 999 that occur when the telescope server dies
time   = data_in.mjd_obs
paral  = data_in.lbt_para
idx999 = WHERE(paral EQ 999, n999, COMPLEMENT=idx_ok)
IF n999 GT 0 THEN paral = INTERPOL(paral[idx_ok], time[idx_ok], time)
PRINT, 'Number of bad parallactic angle : ', n999

; Keep only images in user-defined parallactic range
IF N_ELEMENTS(drs.para_range) GT 1 THEN BEGIN
  par_lim = drs.para_range
  idx_par = WHERE(paral GE par_lim[0] AND paral LE par_lim[1], n_img, /NULL)
  IF n_img LE 0 THEN MESSAGE, 'No image survived the parallactic angle selection!'
  paral   = paral[idx_par]
  data_in = data_in[idx_par]
  img_in  = img_in[*,*,idx_par]
ENDIF

; Bin images per user-defined parallactic angle bin
IF drs.para_bin THEN BEGIN
  par_range = MAX(paral)-MIN(paral)
  n_par     = FLOOR(par_range/drs.para_bin)
  IF n_par GT 0 THEN BEGIN
    img_out   = FLTARR(n_xpix, n_ypix, n_par)
    par_cur   = FLTARR(n_par)
    FOR i_par=0, n_par-1 DO BEGIN
      idx_par = WHERE(paral GE i_par*drs.para_bin + MIN(paral) AND paral LT (i_par+1)*drs.para_bin + MIN(paral), n)
      If n LE 0 THEN GOTO, SKIP_BIN
      par_cur[i_par] = MEAN(paral[idx_par])
      IF n GT 1 THEN RESISTANT_MEAN, img_in[*,*,idx_par], 3, img_tmp, rms, num_rej, DIMENSION=3, /SILENT
      IF n EQ 1 THEN img_tmp = img_in[*,*,idx_par]
      img_out[0,0,i_par] = img_tmp
      SKIP_BIN:
    ENDFOR
    ; Remove empty bins
    idx_ok = WHERE(par_cur NE 0, n_img)
    img_in  = img_out[*,*,idx_ok]
    data_in = data_in[idx_ok]
    data_in.lbt_para = par_cur[idx_ok]
  ENDIF
ENDIF

; Number of images in the output cube
n_out = FLOOR(n_img/ncoadd)

; Init output image cube
img_out = DBLARR(n_xpix, n_ypix, n_out)

; Loop over the images (shift and coadd)
FOR i_img = 0, n_img-1 DO BEGIN
  ; Shift image
  IF MAX(shft_x) GT 0 OR MAX(shft_y) GT 0 THEN img_shft = SHIFTI(img_in[*,*,i_img], shft_x, shft_y) ELSE img_shft = img_in[*,*,i_img]
  IF i_img MOD ncoadd EQ 0 THEN img_tmp = img_shft ELSE img_tmp += img_shft
  ; Save image to output array if ncoadd coadds (doesn't work without the ABS!!!)
  IF ABS(i_img+1) MOD ncoadd EQ 0 THEN img_out[0,0,FLOOR(i_img/ncoadd)] = img_tmp/ncoadd
ENDFOR

; Subtract median if set
img_med = MEDIAN(img_out, DIMENSION=3)
IF KEYWORD_SET(MEDIAN) THEN FOR i = 0, n_out-1 DO img_out[*,*,i]=img_out[*,*,i]-img_med

; Save image file and median image
p=STRPOS(file,'_ALL')
IF p GT 0 THEN BEGIN
   outfile=STRMID(file, 0, p) + '_SEL.fits' 
   vipfile=STRMID(file, 0, p) + '_VIP.fits' 
   medfile=STRMID(file, 0, p) + '_SEL-MED.fits'
ENDIF ELSE MESSAGE, 'Invalid file format.'
FXADDPAR, hdr_in, 'NAXIS3', n_out & MWRFITS, img_out, outfile, hdr_in, /CREATE, /SILENT  
FXADDPAR, hdr_in, 'NAXIS3', n_out & MWRFITS, img_out, vipfile, hdr_in, /CREATE, /SILENT  
FXADDPAR, hdr_in, 'NAXIS3', 1     & MWRFITS, img_med, medfile, hdr_in, /CREATE, /SILENT

; Consistancy check
IF N_ELEMENTS(data_in[*]) NE n_img THEN MESSAGE, 'The length of the data files does not match the size of the image cube'

; Init output data
data_out = REPLICATE(data_in[0], n_out)

; Adjust parallactic angle
idx_neg = WHERE(data_in.lbt_para LT 0, n_neg)
;IF n_neg GT 0 THEN data_in[idx_neg].lbt_para += 360

; Loop over the data
FOR i_out = 0, n_out-1 DO BEGIN
  ; Data index
  idx_data = i_out*LONG(ncoadd) + INDGEN(ncoadd) 
  ; Average parameters or use mid value for strings (TEMPORARY implementation)
  data_out[i_out].mjd_obs = MEAN(data_in[idx_data].mjd_obs)
  data_out[i_out].lbt_utc = data_in[idx_data[0.5*ncoadd]].lbt_utc
  data_out[i_out].lbt_lst = data_in[idx_data[0.5*ncoadd]].lbt_lst
  data_out[i_out].lbt_alt = MEAN(data_in[idx_data].lbt_alt)
  data_out[i_out].lbt_az  = MEAN(data_in[idx_data].lbt_az)
  data_out[i_out].lbt_para = MEAN(data_in[idx_data].lbt_para)
  data_out[i_out].file_id  = MEAN(data_in[idx_data].file_id)
  data_out[i_out].chp_id   = MEAN(data_in[idx_data].chp_id)
  data_out[i_out].xcen     = MEAN(data_in[idx_data].xcen)
  data_out[i_out].ycen     = MEAN(data_in[idx_data].ycen)
  data_out[i_out].slope    = MEAN(data_in[idx_data].slope)
ENDFOR

; Append extension with parameters that change frame by frame
n_row = N_ELEMENTS(data_out.mjd_obs)
FXBHMAKE, hdr, n_row, /INIT, EXTVER=1, 'DATA_SERIES', 'Frame to frame variable parameters'

; Init column number
n_col = 11
col   = LINDGEN(n_col)+1L

; Fill extension header with column names
FXBADDCOL, 1L, hdr, data_out[0].mjd_obs,    'MJD_OBS',   'Modified Julian Date of observation'
FXBADDCOL, 2L, hdr, data_out[0].lbt_utc,    'LBT_UTC',   'UTC from observatory'
FXBADDCOL, 3L, hdr, data_out[0].lbt_lst,    'LBT_LST',   'LST from observatory'
FXBADDCOL, 4L, hdr, data_out[0].lbt_alt,    'LBT_ALT',   'ALT from observatory'
FXBADDCOL, 5L, hdr, data_out[0].lbt_az,     'LBT_AZ',    'AZ from observatory'
FXBADDCOL, 6L, hdr, data_out[0].lbt_para,   'LBT_PARA',  'Parralactic angle from obs.'
FXBADDCOL, 7L, hdr, data_out[0].file_id,    'FILE_ID',   'File identification number'
FXBADDCOL, 8L, hdr, data_out[0].chp_id,     'CHP_ID',    'Chop identification number'
FXBADDCOL, 9L, hdr, data_out[0].xcen,       'XCEN',      'X position of the star'
FXBADDCOL, 10L, hdr, data_out[0].ycen,      'YCEN',      'Y position of the star'
FXBADDCOL, 11L, hdr, data_out[0].slope,     'SLOPE',     'Slope of the fitted Moffat model'

; Write extension header to FITS file
FXBCREATE, unit, outfile, hdr
FXBWRITM,  unit, col, data_out.mjd_obs, data_out.lbt_utc, data_out.lbt_lst, data_out.lbt_alt, data_out.lbt_az, data_out.lbt_para, data_out.file_id, $
                      data_out.chp_id, data_out.xcen, data_out.ycen, data_out.slope
FXBFINISH, unit

; Write angle file for VIP
vipang = STRMID(file, 0, p) + '_PARAL_VIP.fits' 
MWRFITS, data_out.lbt_para, vipang, hdr_in, /CREATE, /SILENT

IF KEYWORD_SET(PLOT) THEN BEGIN
  ; Initiate plot parameters
  fit    = 20./1720.
  inv    = 0
  thick  = 4.0
  xthick = 3.5
  ythick = xthick
  cthick = 3.5
  
  ; Plot Histogram of distance to mean center
  IF drs.sig_pos AND n_dst GT 100 THEN BEGIN
    binsize = (MAX(dst) - MIN(dst))/50.
    hist    = HISTOGRAM(dst, BINSIZE=binsize, MAX=MAX(dst), MIN=MIN(dst))
    n_bin   = N_ELEMENTS(hist)
    bins    = (MIN(dst) + FINDGEN(n_bin)/(n_bin-1.) * (MAX(dst)-MIN(dst)))
    PLOT, [0], [0], XTITLE='Distance to mean center', YTITLE='Number of occurence', XRANGE=[MIN(dst),MAX(dst)], YRANGE=[0.,1.2*MAX(hist)], XSTYLE=1, YSTYLE=1, CHARSIZE=charsize
    OPLOT, bins, hist, COLOR=color
  ENDIF
  ; Plot Histogram of slope
;  IF KEYWORD_SET(SIG_SLO) AND n_slope GT 100 THEN BEGIN
;    binsize = (MAX(slope) - MIN(slope))/50.
;    hist    = HISTOGRAM(slope, BINSIZE=binsize, MAX=MAX(slope), MIN=MIN(slope))
;    n_bin   = N_ELEMENTS(hist)
;    bins    = (MIN(slope) + FINDGEN(n_bin)/(n_bin-1.) * (MAX(slope)-MIN(slope)))
;    PLOT, [0], [0], XTITLE='Moffat profile slope', YTITLE='Number of occurence', XRANGE=[MIN(slope),MAX(slope)], YRANGE=[0.,1.2*MAX(hist)], XSTYLE=1, YSTYLE=1, CHARSIZE=charsize
;    OPLOT, bins, hist, COLOR=color
; ENDIF

  ; Plot measured image (split nods)
  IF KEYWORD_SET(psf_file) THEN BEGIN
    ; Plot name
    plotfile = STRMID(file, 0, STRPOS(file,'.fits'))
    
    nod_id   = data_in.nod_id
    nod_uniq = nod_id[UNIQ(nod_id, SORT(nod_id))]
    n_nod    = N_ELEMENTS(nod_uniq)
    FOR in = 0, n_nod-1 DO BEGIN
      idx_nod = WHERE(nod_id EQ nod_uniq[in])
      img_nod = img_out[*,*,idx_nod]
      img_med = MEDIAN(img_nod, DIMENSION=3)
      xrange  = [-n_xpix,n_ypix]*pixsize*0.5
      yrange  = xrange
      PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1000, 1000]*fit, FILENAME = plotfile + '_NOD' + STRING(nod_uniq[in], FORMAT='(i0)') + '_MED.eps'
      LOADCT, 13, /SILENT
      PLOTXY, img_med, PSYM=10, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='Nod median', XTITLE='Angular separation [arcsec]', YTITLE= 'Angular separation [arcsec]', GRID=0, $
              XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NOERASE, WINDOW=[200,150,950,900]*fit;, INSET_UR='a.'
      ;PLOTXY, [res,res], yrange, /ADD, PSYM=10, COLOR=255, THICK=thick, LINESTYLE=2
      ;PLOTXY, -[res,res], yrange, /ADD, PSYM=10, COLOR=255, THICK=thick, LINESTYLE=2
      PLOTXY, /FIN
      LOADCT, 0, /SILENT
        
      ; Plot cut
      r_psf   = (-FIX(0.5*n_ext) + DINDGEN(n_ext))*pixsize
      r_data  = (-FIX(0.5*n_xpix) + DINDGEN(n_xpix))*pixsize
      img_fit = SFIT(img_med, 3)
      img_med = img_med-img_fit
      xrange  = [-n_ext,n_ext]*pixsize*0.5
      yrange  = [0, 1.2]
      PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1000, 1000]*fit, FILENAME = plotfile  + '_NOD' + STRING(nod_uniq[in], FORMAT='(i0)') + '_CUT.eps'
      PLOTXY, [0], [0], PSYM=10, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='PSF Cross Sections', XTITLE='Angular separation [arcsec]', YTITLE= 'Normalized intensity', GRID=0, $
              XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NOERASE, WINDOW=[200,150,950,900]*fit, /NODATA;, INSET_UR='a.'
      LOADCT, 13, /SILENT
      PLOTXY, r_psf, img_psf[*,0.5*n_ext], /ADD, COLOR=90, THICK=thick
      PLOTXY, r_data, img_med[*,0.5*n_ypix]/MAX(img_med[*,0.5*n_ypix]), /ADD, COLOR=250, THICK=thick
      PLOTXY, /FIN
      LOADCT, 0, /SILENT
    ENDFOR
  ENDIF
ENDIF
END