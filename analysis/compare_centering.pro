PRO COMPARE_CENTERING, date, ob_id

; INPUTS
;   date  : data date (format = 'yymmdd')
;   ob_id : Nod ID of the data 
;   
; KEYWORDS
;
; NOTES
;   fit_method=2 gives poor results. fit_method=5 is the more accurate but very time consuming and not much better then fit_method=1
;

; Define telescope parameters
GET_CNF, cnf, LAMBDA=1.1D-5, INSTRUM='NOMIC'
psf_fwhm = cnf.psf_fwhm

file_nod = '/Users/denis/Desktop/UT2016-10-16_ID011_SCI_alf_Cep_DIT-53ms_11um_PHOT2_IMG.fits'

; Declare path
;DECLARE_PATH, pth, INSTRUM='nomic'

; Read data
;data_path = pth.l1fits_path + date + pth.sep
;file_nod = FILE_SEARCH(data_path, '*_ID' + STRING(ob_id, FORMAT='(I02)') + '*_PHOT1_IMG.fits', COUNT=n0)
img_data = READFITS(file_nod[0], header, /SILENT)

; Generate PSF
n_img  = 100
nxpix  = 200
nypix  = 200
tt_rms = 3
;img_data = FLTARR(nxpix,nypix,n_img)
x_r      = FLTARR(n_img)
y_r      = FLTARR(n_img)
;FOR i = 0, n_img-1 DO BEGIN
;   img_tmp =  PSF_GAUSSIAN(FWHM=cnf.psf_fwhm,npixel=nxpix,ndimen=2,/double)
;   x_r[i]  = tt_rms*RANDOMN(noseed)
;   y_r[i]  = tt_rms*RANDOMN(noseed)
;  img_data[*,*,i] = SHIFTI(img_tmp, x_r[i], y_r[i])
;ENDFOR

; Derive number of images
n_img = N_ELEMENTS(img_data[0,0,*])
nxpix = N_ELEMENTS(img_data[*,0,0])
nypix = N_ELEMENTS(img_data[0,*,0])

; Initial position of the star
xcen0 = 0.5*(nxpix-1)
ycen0 = 0.5*(nypix-1)

; Frame centering
t0 = SYSTIME(1)
fit_method = 1
xcen_all1 = FLTARR(n_img)
ycen_all1 = xcen_all1
FOR i_img = 0, n_img-1 DO BEGIN
  ;PLOTXY, img_data[*,*,i_img]
  img_tmp = PSF_FIT(img_data[*,*,i_img], xcen0, ycen0, psf_fwhm, FIT_METHOD=fit_method, FILE=file, MEDIAN=0, OFFSET=offset, ZOOM_RATIO=1, $
                    XCEN=xcen, YCEN=ycen, SLOPE=slope, INFO=info, PLOT=plot)
  xcen_all1[i_img] = xcen
  ycen_all1[i_img] = ycen
  IF KEYWORD_SET(PLOT) THEN BEGIN
    PLOTXY, /INIT, NWINDOW=1, WINDOW=[0, 0, 620, 620]
    PLOTXY, img_tmp, NWINDOW=1, /NEW, /NOERASE, WINDOW=[80,50,580,550]
    LOADCT, 0, /SILENT
    WAIT, 0.5
  ENDIF
ENDFOR
tt1 = SQRT((xcen_all1-(xcen0+x_r))^2+(ycen_all1-(ycen0+y_r))^2)*cnf.pix_size
AVGSDV, tt1, avg_tt1, rms_tt1, KAPPA=3
t1 = SYSTIME(1)
PRINT, 'fit_method=1 (avg, rms, time) : ', avg_tt1, rms_tt1, t1-t0

; Frame centering
fit_method = 2
xcen_all2 = FLTARR(n_img)
ycen_all2 = xcen_all2
FOR i_img = 0, n_img-1 DO BEGIN
  img_tmp = PSF_FIT(img_data[*,*,i_img], xcen0, ycen0, psf_fwhm, FIT_METHOD=fit_method, FILE=file, MEDIAN=0, OFFSET=offset, ZOOM_RATIO=1, $
                    XCEN=xcen, YCEN=ycen, SLOPE=slope, INFO=info, PLOT=plot)
  xcen_all2[i_img] = xcen
  ycen_all2[i_img] = ycen
  IF KEYWORD_SET(PLOT) THEN BEGIN
    PLOTXY, /INIT, NWINDOW=1, WINDOW=[0, 0, 620, 620]
    PLOTXY, img_tmp, NWINDOW=1, /NEW, /NOERASE, WINDOW=[80,50,580,550]
    LOADCT, 0, /SILENT
    WAIT, 0.5
  ENDIF
ENDFOR
tt2 = SQRT((xcen_all2-(xcen0+x_r))^2+(ycen_all2-(ycen0+y_r))^2)*cnf.pix_size
AVGSDV, tt2, avg_tt2, rms_tt2, KAPPA=3
t2 = SYSTIME(1)
PRINT, 'fit_method=2 (avg, rms, time) : ', avg_tt2, rms_tt2, t2-t1

; Frame centering
fit_method = 3
xcen_all3 = FLTARR(n_img)
ycen_all3 = xcen_all3
fwhm_all3 = xcen_all3
FOR i_img = 0, n_img-1 DO BEGIN
  img_tmp = PSF_FIT(img_data[*,*,i_img], xcen0, ycen0, psf_fwhm, FIT_METHOD=fit_method, FILE=file, MEDIAN=0, OFFSET=offset, ZOOM_RATIO=1, $
                    XCEN=xcen, YCEN=ycen, FWHM_X=fwhm_x, FWHM_Y=fwhm_y, SLOPE=slope, INFO=info, PLOT=plot)
  xcen_all3[i_img] = xcen
  ycen_all3[i_img] = ycen
  fwhm_all3[i_img] = SQRT(fwhm_x^2+fwhm_y^2)
  IF KEYWORD_SET(PLOT) THEN BEGIN
    PLOTXY, /INIT, NWINDOW=1, WINDOW=[0, 0, 620, 620]
    PLOTXY, img_tmp, NWINDOW=1, /NEW, /NOERASE, WINDOW=[80,50,580,550]
    LOADCT, 0, /SILENT
    WAIT, 0.5
  ENDIF
ENDFOR
tt3 = SQRT((xcen_all3-(xcen0+x_r))^2+(ycen_all3-(ycen0+y_r))^2)*cnf.pix_size
AVGSDV, tt3, avg_tt3, rms_tt3, KAPPA=3
AVGSDV, (psf_fwhm-fwhm_all3), avg_fwhm3, rms_fwhm3, KAPPA=3
t3 = SYSTIME(1)
PRINT, 'fit_method=3 (avg, rms, fwhm, time) : ', avg_tt3, rms_tt3, rms_fwhm3, t3-t2

; Frame centering
fit_method = 6
xcen_all5 = FLTARR(n_img)
ycen_all5 = xcen_all5
fwhm_all5 = xcen_all5
FOR i_img = 0, n_img-1 DO BEGIN
  img_tmp = PSF_FIT(img_data[*,*,i_img], xcen0, ycen0, psf_fwhm, FIT_METHOD=fit_method, FILE=file, MEDIAN=0, OFFSET=offset, ZOOM_RATIO=1, $
                    XCEN=xcen, YCEN=ycen, SLOPE=slope, FWHM_X=fwhm_x, FWHM_Y=fwhm_y, INFO=info, PLOT=plot)
  xcen_all5[i_img] = xcen
  ycen_all5[i_img] = ycen
  fwhm_all5[i_img] = SQRT(fwhm_x^2+fwhm_y^2)
  IF KEYWORD_SET(PLOT) THEN BEGIN
    PLOTXY, /INIT, NWINDOW=1, WINDOW=[0, 0, 620, 620]
    PLOTXY, img_tmp, NWINDOW=1, /NEW, /NOERASE, WINDOW=[80,50,580,550]
    LOADCT, 0, /SILENT
    WAIT, 0.5
  ENDIF
ENDFOR
tt5 = SQRT((xcen_all5-(xcen0+x_r))^2+(ycen_all5-(ycen0+y_r))^2)*cnf.pix_size
AVGSDV, tt5, avg_tt5, rms_tt5, KAPPA=3
AVGSDV, (psf_fwhm-fwhm_all5), avg_fwhm5, rms_fwhm5, KAPPA=3
t5 = SYSTIME(1)
PRINT, 'fit_method=5 (avg, rms, fwhm, time) : ', avg_tt5, rms_tt5, rms_fwhm5, t5-t3

; Plot time sequence
x = INDGEN(n_img)
PLOTXY, /INIT, NWINDOW=1, WINDOW=[0, 0, 920, 620]
PLOTXY, [0], [0], NWINDOW=1, /NEW, XRANGE=[0,n_img-1], YRANGE=[-50,50], GRID=0, CHARSIZE=1.0, CHARTHICK=1.5, THICK=2.5, XTITLE='Frame', YTITLE='Tilt [mas]',$
        XSTYLE=1, YSTYLE=1, /NOERASE, /NODATA, WINDOW=[80,50,880,550]
LOADCT, 13, /SILENT
;PLOTXY, x, (tt1-avg_tt1), COLOR=90,  LINESTYLE=0, /ADD
;PLOTXY, x, (tt2-avg_tt2),  COLOR=150, LINESTYLE=0, /ADD
;PLOTXY, x, (tt3-avg_tt3),  COLOR=240, LINESTYLE=0, /ADD
LOADCT, 0, /SILENT

; Now plot the image sequence 
xrange = [-0.5,0.5]*nxpix*cnf.pix_size
yrange = [-0.5,0.5]*nypix*cnf.pix_size
FOR i_img = 0, n_img-1 DO BEGIN
  PLOTXY, /INIT, NWINDOW=2, WINDOW=[0, 0, 620, 620]
  PLOTXY, img_data[*,*,i_img], NWINDOW=2, /NEW, XRANGE=xrange, YRANGE=yrange, GRID=0, CHARSIZE=1.0, CHARTHICK=1.5, THICK=2.5, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=[80,50,580,550]
  ;PLOTXY, xrange[0]+[0.5*nxpix]*cnf.pix_size, yrange[0]+[0.5*nypix]*cnf.pix_size, SYMBOL=5, COLOR=255, /ADD
  LOADCT, 13, /SILENT
  ;PLOTXY, xrange[0]+[xcen_all1[i_img]]*cnf.pix_size, yrange[0]+[ycen_all1[i_img]]*cnf.pix_size, COLOR=90,  SYMBOL=2, /ADD
  ;PLOTXY, xrange[0]+[xcen_all2[i_img]]*cnf.pix_size, yrange[0]+[ycen_all2[i_img]]*cnf.pix_size, COLOR=150, SYMBOL=2, /ADD
  PLOTXY, xrange[0]+[xcen_all3[i_img]]*cnf.pix_size, yrange[0]+[ycen_all3[i_img]]*cnf.pix_size, COLOR=230, SYMBOL=2, /ADD
  ;PLOTXY, xrange[0]+[xcen_all5[i_img]]*cnf.pix_size, yrange[0]+[ycen_all5[i_img]]*cnf.pix_size, COLOR=250, SYMBOL=3, /ADD
  LOADCT, 0, /SILENT
  WAIT, 0.1
ENDFOR

END