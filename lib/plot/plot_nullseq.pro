;+
; NAME: PLOT_NULLSEQ
; 
; PURPOSE:
;   Main routine to plot results of a given null sequence (null/flux histograms, null sequence, etc)
;
; INPUTS:
;   pth           :  Structure with paths of result folders
;   phot_seq      :  Constructive flux
;   flx           :  Measured flux
;   flx_err       :  Error on the measured flux
;   img_in        :  Image of the detector stored as
;
; KEYWORDS
;   BCK_FLX       :  Backrgound flux per photometric aperture [ADU]
;   PHOT1_FLX     :
;   PHOT2_FLX     :
;   HDR_DATA      :  Mandatory keyword with header information
;   SPEED         :  Set this keyword to 1 to speed-up the movie computation which will however not represent the true acquisition pace anymore
;   PLOT_FILE     :  Set this keyword to the tgt_name of the eps plot
;   DISPLAY       :  Set this keyword to display the null sequence and histogram on screen
;   MOVIE         :  Set this keyword to produce a movie of the sequence
;
; MODIFICATION HISTORY:
;   Version 1.0,  01-FEB-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  04-MAR-2013, DD: Removed keyword movie
;   Version 1.2,  04-APR-2013, DD: Improved movie quality and information display
;   Version 1.3,  01-MAY-2013, DD: Removed various keywords. All file information is now passed through HDR_DATA
;   Version 1.4,  14-OCT-2013, DD: Adapted for variable screen size
;   Version 1.5,  14-MAR-2014, DD: Aded keyword NO_LOG
;   Version 1.6,  26-MAR-2014, DD: Corrected x-axis of histograms plots
;   Version 1.7,  28-MAR-2014, DD: Added photometry output plots
;   Version 1.8,  28-SEP-2014, DD: Completely new implementation

PRO PLOT_NULLSEQ, pth, phot_seq, flx, flx_err, img_in, BCK_FLX=bck_flx, PHOT1_FLX=phot1_flx, PHOT2_FLX=phot2_flx, HDR_DATA=hdr_data, SPEED=speed, DISPLAY=display, MOVIE=movie, NO_LOG=no_log
   
; Keyword sanity check
IF NOT KEYWORD_SET(speed) THEN speed = 0.

; Get screen size
screen = GET_SCREEN_SIZE(RESOLUTION=resolution)

; Reduction parameters
rec_time  = hdr_data.mjd_obs
date      = hdr_data[0].date_obs                                   ; Date at the beginning of the sequence
tgt_name  = hdr_data[0].objname                                    ; Object name
time_span = REFORM(TRANSPOSE(rec_time-MIN(rec_time))*24.*60.*60.)  ; Convert to sec (rec_time is in MJD)
pix_size  = hdr_data[0].pixscale
aper_rad  = hdr_data[0].aper_rad
bck_irad  = hdr_data[0].bck_irad
bck_orad  = hdr_data[0].bck_orad
int_uniq  = hdr_data[0].int_time
coad_uniq = hdr_data[0].n_coadd
smod_uniq = hdr_data[0].smplmode
pag_uniq  = hdr_data[0].pagain
pab_uniq  = hdr_data[0].pabandw
lam_uniq  = hdr_data[0].lam_cen
i_ob      = hdr_data[0].ob_id
n_img     = N_ELEMENTS(img_in[0,0,*])                              ; Number of image frames

; Define filename (plot + movie)
file_name = date + '_OB-' + STRING(i_ob, FORMAT='(I0)') + '_' + tgt_name + '_null_dit-' + STRING(1000.*int_uniq, FORMAT='(I0)') + 'ms' + $
            '_coa-' + STRING(coad_uniq, FORMAT='(I0)') + '_wav-' + STRING(1D+6*lam_uniq, FORMAT='(I0)') + 'um_mode-' + STRING(smod_uniq, FORMAT='(I0)') + '_pag-' + STRING(pag_uniq, FORMAT='(I0)') + $
            '_pab-' + STRING(pab_uniq, FORMAT='(I0)') 

; Define result paths
; 1. Plot path
plot_path = pth.result_path + 'null'  + pth.sep + date + pth.sep            ; Plot path
IF NOT FILE_TEST(plot_path) THEN FILE_MKDIR, plot_path                      ; Create directory if it does not exist
plot_file = plot_path + file_name
; 2. Movie papth
IF KEYWORD_SET(MOVIE) THEN BEGIN
  mov_path  = pth.result_path + 'movie' + pth.sep + date + pth.sep          ; Plot path
  IF NOT FILE_TEST(mov_path)  THEN FILE_MKDIR, mov_path                     ; Create directory if it does not exist
  mov_file  =  mov_path + file_name + '.mpg'
ENDIF

; Null measurements and corresponding errors
null     = flx/phot_seq 
null_err = flx_err/phot_seq
n_null   = N_ELEMENTS(null)
                                                           
; Initiate plot parameters    
LOADCT, 0, /SILENT  
fit      = 20./1720.
inv      = 0. 
  
; 1. Plot the instantenous null
; *****************************

fr_time   = time_span[1]-time_span[0]
xrange    = [MIN(time_span),MAX(time_span)]
IF MIN(null) LT -10. THEN BEGIN
  yrange = [-10.,100.]
  ylog   = 0
  text_pos2 = -8
  text_pos1 = -4 
ENDIF ELSE BEGIN
  yrange = [0.1, 100.]
  ylog   = 1
  text_pos2 = 0.22
  text_pos1 = 0.15
ENDELSE
PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plot_file + '.eps'
PLOTXY, REFORM(time_span), null, YLOG=ylog, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='N-band nulling with LBTI-NOMIC', XTITLE='Elapsed time [s]', YTITLE='Instantaneous null [%]', GRID=0, $
        XSTYLE=1, YSTYLE=1, XTHICK=3.5, YTHICK=3.5, THICK=3., CHARTHICK=2.0, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
LOADCT, 0, /SILENT
DAYCNV, rec_time[n_img-1] +2400000.5D0, yr, mn, day, hr
hr_floor = FLOOR(hr)
min      = (hr-FLOOR(hr))*60.
sec      = (min-FLOOR(min))*60.
XYOUTS, 8*fr_time, text_pos1, tgt_name + ' on ' + date + ' at ' + STRING(hr_floor, FORMAT='(I2.2)') + ':' + STRING(min, FORMAT='(I2.2)') + ':' + STRING(sec, FORMAT='(I2.2)') + ' UT', CHARSIZE=0.8, CHARTHICK=2.0, COLOR=0
;XYOUTS, 3*fr_time, 166, 'Theoretical astrophysical null: ' + STRING(100*geom_null, FORMAT='(F4.1)') + '%', CHARSIZE=2.0, COLOR=0
;XYOUTS, 3*fr_time, 156, 'Measured instantaneous null: ' + STRING(null[n_bin-1], FORMAT='(F4.1)') + '%', CHARSIZE=2.0, COLOR=0
;XYOUTS, 3*fr_time, 146, 'Instantaneous instrumental null: ' + STRING(null_inst[n_bin-1], FORMAT='(F4.1)') + '%', CHARSIZE=2.0, COLOR=0
XYOUTS, 8*fr_time, text_pos2, 'Best measured null: ' + STRING(MIN(100*null, idx_min), FORMAT='(F4.1)') + '% !9+!3' + STRING(100*null_err[MIN(idx_min)], FORMAT='(F4.1)') + '%',  CHARSIZE=0.8, CHARTHICK=2.0, COLOR=0  
LOADCT, 1, /SILENT
PLOTXY, time_span, 100.*null, /ADD, LINESTYLE=0, COLOR=210, THICK=4
LOADCT, 13, /SILENT
;PLOTXY, [0.,1D+16], 100*[geom_null,geom_null], /ADD, LINESTYLE=2, CHARSIZE=0.5, COLOR=245, THICK=4 
;XYOUTS, 0.33*fr_time*n_bin, -8, 'Theoretical astrophysical null (!7h!3!D*!N of ' + star.name + ' = ' + STRING(star.ldm, FORMAT='(F4.1)') +' mas)', CHARSIZE=0.7, COLOR=245
  
; Second plot: overplot the inset with the detector
n_clip   = N_ELEMENTS(img_in[0,*,0])
LOADCT, 0, /SILENT
yrange = 0.5*[-n_clip,n_clip]*pix_size[0]
xrange = 0.5*[-n_clip,n_clip]*pix_size[0]
!P.TITLE = 'Image at minimum'
!X.TITLE = 'Angular offset [arcsec]'
!Y.TITLE = 'Angular offset [arcsec]'
null_max = MAX(100*null, idx_max)
PLOTXY, REFORM(img_in[*,*,idx_max]), XRANGE=xrange, YRANGE=yrange, TITLE='Image at peak', XTITLE='Angular offset [arcsec]', YTITLE='Angular offset [arcsec]', /NEW, GRID=0, THICK=3.0, XSTYLE=1, YSTYLE=1, CHARSIZE=0.5, /NOERASE, WINDOW=[920,160,1100,340]*fit;, INSET_UR='a.'
IF KEYWORD_SET(APER_RAD) THEN BEGIN
  LOADCT, 13, /SILENT
  TVCIRCLE, aper_rad*pix_size, 0., 0., 90, THICK=2.5, LINESTYLE=0, /DATA
  TVCIRCLE, bck_irad*pix_size, 0., 0., 250, THICK=2.5, LINESTYLE=0, /DATA
  TVCIRCLE, bck_orad*pix_size, 0., 0., 250, THICK=2.5, LINESTYLE=0, /DATA
  LOADCT, 0, /SILENT
ENDIF
;COLORBAR, XTITLE='', YTITLE='', /RIGHT, RANGE=[MIN(img_null[*,*,i_f]),MAX(img_null[*,*,i_f])],  POSITION=[0.77,0.15,0.83,0.85], FORMAT='(F6.3)', CHARTHICK=3.0, THICK=1.0, /VERTICAL
PLOTXY, /FIN
LOADCT, 0, /SILENT  

; Plot instantenous null and creat a movie if requested
; *****************************************************

IF KEYWORD_SET(movie) OR KEYWORD_SET(display) THEN BEGIN
  ; Create temporary path for jpeg plots
  temp_path = plot_path + pth.sep + 'temp' + pth.sep          ; temporary folder to store the images for the movie
  FILE_MKDIR, temp_path                                       ; create directory
      
  ; Compute frame rate if a movie has to be created
  IF speed EQ 0 AND KEYWORD_SET(movie) THEN BEGIN
    fr_rate  = 1/fr_time                   ; [fps]
    IF fr_rate LE 24. THEN BEGIN
      mpg_rate = 24. & fps_id = 2     ; [fps], standard mpeg fps
    ENDIF ELSE BEGIN
      IF fr_rate GT 24. AND fr_rate LE 40. THEN BEGIN
       mpg_rate = 30. & fps_id = 5    ; [fps], NTSC drop frame video frame rate (the default)
      ENDIF ELSE BEGIN
        IF fr_rate GT 40. AND fr_rate LE 55. THEN BEGIN
         mpg_rate = 50. & fps_id = 6    ; [fps], NTSC drop frame video frame rate (the default)
        ENDIF ELSE BEGIN
         mpg_rate = 60. & fps_id = 8    ; [fps], double frame rate NTSC drop frame video 
         IF fr_rate GT 60. THEN PRINT, "WARNING: the frame rate is too fast and will be limited to 60 Hz in the movie"
        ENDELSE
      ENDELSE
    ENDELSE
    ratio = FLOOR(mpg_rate/fr_rate) > 1      ; This is too account for very slow frame rates (#ratio identical images will be produced).
  ENDIF ELSE ratio = 1   
  
  ; Initiate plot
  fit        = 1. ;20./1720.
  inv        = 0
  IF MIN(null) LT -10 OR KEYWORD_SET(NO_LOG) THEN BEGIN
    yrange = [-10.,100.]
    ylog   = 0
  ENDIF ELSE BEGIN
    yrange = [0.1, 100.]
    ylog   = 1
  ENDELSE
  xrange     = [MIN(time_span),MAX(time_span)]
  nwindow    = 4
  sx = 0 & sy = screen[0]-800  ; shift window
  PLOTXY, /INIT, INV=inv, NWINDOW=nwindow, /JPEG, /COLOR, WINDOW=[0+sx, 0+sy, 1200+sx, 800+sy]*fit    
  PLOTXY, [0], [0], YLOG=ylog, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='N-band nulling with LBTI-NOMIC', XTITLE='Elapsed time [s]', YTITLE='Instantaneous null [%]', GRID=0, $
          CHARSIZE=2.5, CHARTHICK=2.5, THICK=3.5, XSTYLE=1, YSTYLE=1, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.' 
       
  ; Loop over the files
  i_c = 0                                ; initiate file number counter
  FOR i_f=0, n_img-1 DO BEGIN  
    ; Extract the null sequence
    null_tmp  = 100.*null[INDGEN(i_f+1)]
    min_null  = MIN(null_tmp[INDGEN(i_f+1)], i_min)
    min_nerr  = 100.*null_err[i_min]
    img_out   = img_in[*,*,i_f]          ; current image
    
    ; Compute elasped time
    elaps_time = time_span[INDGEN(i_f+1)]
    
    ; Put the (0,0) pixel to the maximum value for keep the same scaling between plots
    img_out[0,0] = MAX(img_in)
      
    ; Loop over the various duplicate images
    FOR i_d=0, ratio-1 DO BEGIN
        ; File name
        img_name = 'nomic_img_' + STRING(i_c, FORMAT='(I4.4)') + '.jpg'
  
        ; First plot: the instantenous null
        xrange = [MIN(time_span),MAX(time_span)]
        PLOTXY, [0], [0], YLOG=ylog, /NEW, NWINDOW=nwindow, XRANGE=xrange, YRANGE=yrange, TITLE='N-band nulling with LBTI-NOMIC', XTITLE='Elapsed time [s]', YTITLE='Instantaneous null [%]', GRID=0, $
                 CHARSIZE=2.5, CHARTHICK=2.5, THICK=4.0, XSTYLE=1, YSTYLE=1, /NODATA, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
        PLOTXY, elaps_time, null_tmp, /ADD, NWINDOW=nwindow,  LINESTYLE=0, COLOR=250, THICK=4
        DAYCNV, rec_time[i_f] +2400000.5D0, yr, mn, day, hr
        hr_floor = FLOOR(hr)
        min      = 60.*(hr-FLOOR(hr))   & min_floor = FLOOR(min)  
        sec      = 60.*(min-FLOOR(min)) & sec_floor = FLOOR(sec)
        XYOUTS, 3.5, -3,  tgt_name + ' on ' + date + ' at ' + STRING(hr_floor, FORMAT='(I2.2)') + ':' + STRING(min_floor, FORMAT='(I2.2)') + ':' + STRING(sec_floor, FORMAT='(I2.2)') + ' UT', CHARSIZE=2.0, CHARTHICK=2.2, COLOR=220
        XYOUTS, 3.5, -7, 'Instantaneous null: ' + STRING(100.*null[i_f], FORMAT='(F5.2)') + '% !9+!3 ' + STRING(100.*null_err[i_f], FORMAT='(F4.2)') + '%', CHARSIZE=2.0, CHARTHICK=2.2, COLOR=220
        ;XYOUTS, 3.5, -9, 'Best null: ' + STRING(min_null, FORMAT='(F5.2)') + '% !9+!3 ' + STRING(min_nerr, FORMAT='(F4.2)') + '%', CHARSIZE=2.0, CHARTHICK=2.2, COLOR=220  
  
        ; Second plot: overplot the inset with the detector
        n_clip   = N_ELEMENTS(img_in[0,*,0])
        LOADCT, 0, /SILENT
        yrange2 = 0.5*[-n_clip,n_clip]*pix_size[0]
        xrange2 = 0.5*[-n_clip,n_clip]*pix_size[0]
        PLOTXY, CONGRID(img_out, 256, 256), /NEW, NWINDOW=nwindow, XRANGE=xrange2, YRANGE=yrange2, TITLE='Real-time detector view', XTITLE='Angular offset [arcsec]', YTITLE='Angular offset [arcsec]', GRID=0, THICK=3.0, XSTYLE=1, YSTYLE=1, CHARSIZE=1.4, $
                CHARTHICK=1.5, /NOERASE, WINDOW=[880,440,1100,660]*fit;, INSET_UR='a.'
        ; Overplot circles for PSF and background computation
        IF KEYWORD_SET(APER_RAD) THEN BEGIN
  	      LOADCT, 13, /SILENT
  	      TVCIRCLE, aper_rad*pix_size, 0., 0., 90, THICK=2.5, LINESTYLE=0, /DATA
  	      TVCIRCLE, bck_irad*pix_size, 0., 0., 250, THICK=2.5, LINESTYLE=0, /DATA
  	      TVCIRCLE, bck_orad*pix_size, 0., 0., 250, THICK=2.5, LINESTYLE=0, /DATA
  	      LOADCT, 0, /SILENT   
        ENDIF
        ;COLORBAR, XTITLE='', YTITLE='', /RIGHT, RANGE=[MIN(img_null[*,*,i_f]),MAX(img_null[*,*,i_f])],  POSITION=[0.77,0.15,0.83,0.85], FORMAT='(F6.3)', CHARTHICK=3.0, THICK=1.0, /VERTICAL
            
        ; Save display
        IF KEYWORD_SET(movie) THEN SAVEIMAGE, temp_path + img_name, /JPEG, QUALITY=100., DITHER=DITHER, CUBE=CUBE, /QUIET       
   
        ; Increment file number counter
        i_c = i_c + 1
    ENDFOR
  ENDFOR
  
  ; Close plot
  PLOTXY, /FIN
  LOADCT, 0, /SILENT  
  
  ; Now make the movie
  IF KEYWORD_SET(movie) THEN BEGIN
    make_mpg, prefix = temp_path + 'nomic_img_', suffix='jpg', n_start=0, n_end=i_c-1, frame_rate=fps_id, digits=4, mpeg_file = mov_file
    FILE_DELETE, temp_path, /RECURSIVE                  ; delete temporary directory  
  ENDIF
ENDIF      
END
