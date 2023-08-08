; +
; NAME: OPD_VS_MAG
;   
; PURPOSE:
;   This function plots various nulling metrics versus one another (initially used to plot the OPD jitter versus the stellar magnitude, hence the name).  
;
; INPUTS:
;   date : vector of dates to be reduced (e.g., ['150208','150309'])
;
; MODIFICATION HISTORY:
;   Version 1.0,  23-JAN-2015, Denis Defrère, University of Arizona

PRO OPD_VS_MAG, date

; Start actual code
ON_ERROR, 0

; Keyword sanity check
wav = 2.2

; Define running and plotting paramaters
n_lam     = 20                     ;Number of wavelength bins within the bandwidth
n_time    = 1000                   ;Number of points in the interpolated TF curve
charsize  = 1.2
charthick = 3.0

; DEFINE GLOBAL VARIABLES
; ***********************

; Astronomical and physical constants
GET_PRM, prm

; Obtain the definition of the configuration
GET_CNF, cnf

; Recover the IDL running path
DECLARE_PATH, pth, INSTRUM='nomic'


; PROCESS FILES
; *************

n_date = N_ELEMENTS(date)
FOR i_date=0, n_date-1 DO BEGIN
  ; Derive long version of the date
  date_lng    = '20' + STRMID(date[i_date], 0, 2) + '-' + STRMID(date[i_date], 2, 2) + '-' + STRMID(date[i_date], 4, 2)
  l1fits_path = pth.l1fits_path + date_lng + pth.sep
  IF NOT FILE_TEST(l1fits_path) THEN MESSAGE, 'No L1 data for ' + date_lng
  
  ; Search for L1 file and read it
  IF KEYWORD_SET(VERSION) THEN version = '_v' + STRING(version, FORMAT='(I0)') ELSE version = ''
  l1files = FILE_SEARCH(l1fits_path,'UT' + date_lng + version + '.fits')
  l1data  = MRDFITS(l1files[0], 1, hdr_col, /SILENT)   ; Read null data
  phasec  = MRDFITS(l1files[0], 5, /SILENT)            ; Read PHASECam data

  ; Extract useful data
  tgt_name  = l1data.objname
  int_time  = l1data.int_time
  tgt_alt   = l1data.lbt_alt
  ttdx_rms  = l1data.ttdx_rms
  ttsx_rms  = l1data.ttsx_rms
  phot_dx   = l1data.photdx_avg
  phot_sx   = l1data.photsx_avg
  null_rms  = l1data.NULL_MEAS_ERR
  phstd_avg = phasec.pcphstd*wav/360.
  phstd_err = phasec.pcphstd_err*wav/360.
  
  ; Compute unique objects
  tgt_uniq = tgt_name[UNIQ(tgt_name, SORT(tgt_name))]
  n_tgt    = N_ELEMENTS(tgt_uniq)
  
  ; Loop over target and compute flux
  n_data = N_ELEMENTS(tgt_name)
  magK   = FLTARR(n_data)
  magV   = FLTARR(n_data)
  FOR i_tgt=0, n_tgt-1 DO BEGIN
    ; Get magbitude
    PRINT, 'Retrieving magnitude of ', tgt_uniq[i_tgt]
    querysimbad, tgt_uniq[i_tgt], ra, de, VMAG=Vmag, KMAG=Kmag, /CFA
    IF NOT KEYWORD_SET(KMAG) THEN query_2mass,tgt_uniq[i_tgt],J,H,Kmag,Jerr,Herr,Kserr

    ; Extract data for this object and compute average OPD error
    idx_tgt = WHERE(tgt_name EQ tgt_uniq[i_tgt], n_tmp)
    magV[idx_tgt] = Vmag
    magK[idx_tgt] = Kmag
  ENDFOR
  
  ; Store results
  IF i_date EQ 0 THEN BEGIN
    magV_all = magV
    magK_all = magK
    dit_all  = int_time
    opd_all  = phstd_avg
    ttdx_all = ttdx_rms
    ttsx_all = ttsx_rms
    pdx_all  = phot_dx
    psx_all  = phot_sx
    tel_all  = tgt_alt
    null_all = null_rms
  ENDIF ELSE BEGIN
    magV_all = [magV_all, magV]
    magK_all = [magK_all, magK]
    dit_all  = [dit_all, int_time]
    opd_all  = [opd_all, phstd_avg]
    ttdx_all = [ttdx_all, ttdx_rms]
    ttsx_all = [ttsx_all, ttsx_rms]
    pdx_all  = [pdx_all, phot_dx]
    psx_all  = [psx_all, phot_sx]
    tel_all  = [tel_all, tgt_alt]
    null_all = [null_all,null_rms]
  ENDELSE
ENDFOR
PRINT, magV_all

; Remove O that might happen when no photometric frames
idx_ok = WHERE(ttdx_all GT 0 AND ttsx_all GT 0, n0)
ttdx_all = ttdx_all[idx_ok] & ttsx_all = ttsx_all[idx_ok]
magV_all = magV_all[idx_ok]    
    
; Plot results
xrange = [-2,5]
yrange = [0.,0.8]
PREP_PS, /BOLD & LOADCT, 12, /SILENT
DEVICE, FILEN=pth.result_path + 'phasecam_perfo_OPD-vs-MAG.eps', /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7
PLOT, [0], [0], XTITLE='K-band magnitude', YTITLE='Measured OPD jitter over DIT [um]', TITLE=title, XSTYLE=1, YSTYLE=9, $
      XRANGE=xrange, YRANGE=yrange, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize, POSITION=[0.1,0.1,0.90,0.92]
; Add right vertical axis
AXIS, YTITLE ='Measured OPD jitter (>1/27ms)[rad]', /YAXIS, YSTYLE=1, YRANGE=yrange*!DPI/(wav*0.5)
; Overplot the estimated null values for calibrators measurements
USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
OPLOT, magK_all, opd_all, PSYM=8, COLOR=100
DEVICE, /CLOSE & END_PS

xrange = [30.,MAX(tel_all)>MAX(tel_all)]*1.05
yrange = [0.,0.8]
PREP_PS, /BOLD & LOADCT, 12, /SILENT
DEVICE, FILEN=pth.result_path + 'phasecam_perfo_OPD-vs-ALT.eps', /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7
PLOT, [0], [0], XTITLE='Telescope altitude [deg]', YTITLE='Measured OPD jitter over DIT [um]', TITLE=title, XSTYLE=1, YSTYLE=9, $
  XRANGE=xrange, YRANGE=yrange, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize, POSITION=[0.1,0.1,0.90,0.92]
; Add right vertical axis
AXIS, YTITLE ='Measured OPD jitter (>1/27ms)[rad]', /YAXIS, YSTYLE=1, YRANGE=yrange*!DPI/(wav*0.5)
; Overplot the estimated null values for calibrators measurements
USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
OPLOT, tel_all, opd_all, PSYM=8, COLOR=100
DEVICE, /CLOSE & END_PS

xrange = [-2,5]
yrange = [0.,50]
PREP_PS, /BOLD & LOADCT, 13, /SILENT
DEVICE, FILEN=pth.result_path + 'phasecam_perfo_TT-vs-MAG.eps', /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7
PLOT, [0], [0], XTITLE='V-band magnitude', YTITLE='Measured tip/tilt [mas]', TITLE=title, XSTYLE=1, YSTYLE=1, $
      XRANGE=xrange, YRANGE=yrange, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize, POSITION=[0.1,0.1,0.90,0.92]
; Overplot the estimated null values for calibrators measurements
USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
OPLOT, magV_all, ttdx_all, PSYM=8, COLOR=100
OPLOT, magV_all, ttsx_all, PSYM=8, COLOR=250
DEVICE, /CLOSE & END_PS

xrange = [0,100]
PREP_PS, /BOLD & LOADCT, 13, /SILENT
DEVICE, FILEN=pth.result_path + 'phasecam_perfo_TT-vs-DIT.eps', /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7
PLOT, [0], [0], XTITLE='Integration time [ms]', YTITLE='Measured tip/tilt [mas]', TITLE=title, XSTYLE=1, YSTYLE=1, $
      XRANGE=xrange, YRANGE=yrange, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize, POSITION=[0.1,0.1,0.90,0.92]
; Overplot the estimated null values for calibrators measurements
USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
OPLOT, dit_all*1D+3, ttdx_all, PSYM=8, COLOR=100
OPLOT, dit_all*1D+3, ttsx_all, PSYM=8, COLOR=250
DEVICE, /CLOSE & END_PS

xrange = [30.,MAX(tel_all)>MAX(tel_all)]*1.05
PREP_PS, /BOLD & LOADCT, 13, /SILENT
DEVICE, FILEN=pth.result_path + 'phasecam_perfo_TT-vs-ALT.eps', /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7
PLOT, [0], [0], XTITLE='Telescope altitude [deg]', YTITLE='Measured tip/tilt [mas]', TITLE=title, XSTYLE=1, YSTYLE=1, $
  XRANGE=xrange, YRANGE=yrange, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize, POSITION=[0.1,0.1,0.90,0.92]
; Overplot the estimated null values for calibrators measurements
USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
OPLOT, tel_all, ttdx_all, PSYM=8, COLOR=100
OPLOT, tel_all, ttsx_all, PSYM=8, COLOR=250
DEVICE, /CLOSE & END_PS

yrange = [0.,0.10]
PREP_PS, /BOLD & LOADCT, 13, /SILENT
DEVICE, FILEN=pth.result_path + 'phasecam_perfo_NULL-vs-ALT.eps', /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7
PLOT, [0], [0], XTITLE='Telescope altitude [deg]', YTITLE='Measured null variation', TITLE=title, XSTYLE=1, YSTYLE=1, $
  XRANGE=xrange, YRANGE=yrange, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize, POSITION=[0.1,0.1,0.90,0.92]
; Overplot the estimated null values for calibrators measurements
USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
OPLOT, tel_all, null_all, PSYM=8, COLOR=100
DEVICE, /CLOSE & END_PS

xrange = [0.,MAX(pdx_all)>MAX(psx_all)]*1.05
yrange = [0.,50]
PREP_PS, /BOLD & LOADCT, 13, /SILENT
DEVICE, FILEN=pth.result_path + 'phasecam_perfo_TT-vs-PHOT.eps', /ENCAPS, /COLOR, XSIZE=19.5, YSIZE=14.7
PLOT, [0], [0], XTITLE='Measured photometry [ADU]', YTITLE='Measured tip/tilt [mas]', TITLE=title, XSTYLE=1, YSTYLE=1, $
  XRANGE=xrange, YRANGE=yrange, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize, POSITION=[0.1,0.1,0.90,0.92]
; Overplot the estimated null values for calibrators measurements
USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
OPLOT, pdx_all, ttdx_all, PSYM=8, COLOR=100
OPLOT, psx_all, ttsx_all, PSYM=8, COLOR=250
DEVICE, /CLOSE & END_PS
END
