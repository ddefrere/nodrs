; +
; NAME: IMG_FINDBEAM
; 
; PURPOSE:
;   This function returns the mean beam positions of all img_in.
;
; INPUT:
;   img_in        :  An image data cube
;   fwhm          :  FWHM (in pixels) to be used in the convolve filter
;
; KEYWORDS
;   AUTO_BEAM     :  Set to 1 to turn off user interaction (beam position is forced to one possible value when a confirmation is needed)
;   CURSOR        :  Set to activate CURSOR mode (the user has to click on the image to mark the position of the beam)
;   EDGE          :  Two element vector giving the X,Y channel size (beams found within FWHM of the edges will be discarded)
;   FIT_MODE      :  Define the centroid fitting technique (set it negative to fit inverse function, e.g. for AGPM data):
;                       - 0: no stellar centroid fitting (report position computed on the median combined cube)
;                       - 1: compute the centroid of a star using a derivative search (using CNTRD function)
;                       - 2: compute the stellar centroid by Gaussian fits to marginal X,Y, sums (using GCNTRD function)
;                       - 3: compute the stellar centroid by Gaussian fit using MPFIT2DPEAK (more time consuming)
;                       - 4: compute the stellar centroid by Lorentzian fit using MPFIT2DPEAK
;                       - 5: compute the stellar centroid by Moffat model fit using MPFIT2DPEAK
;   N_BEAM        :  Number of beams to be found in each image (1 if not set)
;   PIXSIZE       :  Pixel size (in mas)
;   PRECISE       :  If set, improved beam centering by calling PSF_FIT.pro to refine the position found by FIND.pro
;   ROI           :  A 4-element vector defining the region of interest (i.e., where to look for a beam). The format is [x_min,y_min,x_max,y_max].
;   INFO          :  Define the level of information printed to screen:
;                       - 0: completely silent execution;
;                       - 1: minimum level of information;
;                       - 2: nominal level of information;
;                       - 3: debugging use.
;
; OUTPUT
;   Array with the coordinates of the source in the coaaded image (n_beam, 2)
;
; MODIFICATION HISTORY:
;   Version 1.0, 08-JUL-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (based on former routine "lbti_imgcen.pro")
;   Version 1.1, 07-NOV-2014, DD: Added slope to output data
;   Version 1.2, 02-FEB-2015, DD: Now handle single frame "cube"
;   Version 1.3, 23-FEB-2015, DD: Now autmatically handle beams near the edge of the frame
;   Version 1.4, 23-NOV-2015, DD: Replaced RESISTANT_MEAN by faster MEDIAN routine 
;   Version 1.4, 03-DEC-2015, DD: Replaced MEDIAN by RESISTANT_MEAN which is more robust + added call to SIGMA_FILTER and GAUSS_SMOOTH
;   Version 1.5, 21-JAN-2016, DD: Added best-fit FWHM to output structure (+ corrected implemtation of slope)
;   Version 1.6, 23-MAY-2016, DD: Now discard positions closer than 'fwhm' pixels from the edge
;   Version 1.7, 15-JUL-2016, DD: Added CURSOR mode
;   Version 1.8, 03-NOV-2016, DD: Added second call to PSF_FIT to improve robustness
;   Version 1.8, 09-JAN-2017, DD: Now pause if a beam is found to cloose to the edge
;   Version 1.9, 14-FEB-2017, DD: Added keyword AUTO_BEAM
;   Version 2.0, 05-APR-2017, DD: Added keyword EDGE
;   Version 2.1, 05-APR-2017, DD: Improved robustness
;   Version 2.2, 18-MAY-2017, DD: Added keyword PLOT

FUNCTION IMG_FINDBEAM, img_in, fwhm, AUTO_BEAM=auto_beam, CURSOR=cursor, EDGE=edge, FIT_MODE=fit_mode, N_BEAM=n_beam, PIXSIZE=pixsize, QUICK=quick, INFO=info, PLOT=plot                                                                      
                       
; Keyword sanity check
IF NOT KEYWORD_SET(FIT_MODE) THEN fit_mode  = 0
IF NOT KEYWORD_SET(N_BEAM)   THEN n_beam    = 1
IF NOT KEYWORD_SET(PIXSIZE)  THEN pixsize   = 1
IF NOT KEYWORD_SET(PLOT)     THEN auto_beam = 1  ; auto beam if no plot!

; Running paramaters
n_keep    = n_beam           ; number of point sources to keep in the list of those found by FIND.pro (goes as decreasing flux)
n_rms     = 4                ; starting threshold for point source (=n_rms*rms)
i_iter    = 0                ; init number of iterations
max_iter  = 20               ; max number of iterations
IF KEYWORD_SET(CURSOR) THEN flx_thr = 1D16 ELSE flx_thr = 1  ; Flux threshold below which the user is asked to confirm the position

; Convertion factor between FWHM and Gaussian sigma
fwhm2sig = 1/(2*SQRT(2*ALOG(2D)))

; Initiate variables and arrays
n_xpix0   = N_ELEMENTS(img_in[*,0,0])
n_ypix0   = N_ELEMENTS(img_in[0,*,0])
n_img     = N_ELEMENTS(img_in[0,0,*])
beam_pos  = FLTARR(n_beam, 2)
fwhm_fit  = FLTARR(n_beam, 2)
slope_pro = FLTARR(n_beam)

; Compute median image if necessary (also use by display)
IF n_img GT 1 THEN img_tot = MEDIAN(img_in, DIMENSION=3) $; RESISTANT_MEAN, img_in, 5, img_tot, rms, num_rej, DIMENSION=3, /SILENT $ (median is faster)
              ELSE img_tot = img_in

; Set fwhm pixels from each channel edge to mean value
img_tmp = img_tot

; Smooth and clean image to improve the automatic beam finding
;img_tmp = SIGMA_FILTER(TEMPORARY(img_tmp), 20, N_SIGMA=5, /ALL_PIXELS, ITERATE=0, MONITOR=monitor, N_CHANGE=nchange)  ; This helps with very faint signals
img_tmp = SIGMA_FILTER(TEMPORARY(img_tmp), 10, N_SIGMA=4, /ALL_PIXELS, ITERATE=0, MONITOR=monitor, N_CHANGE=nchange)  ; This helps with very faint signals
img_tmp = GAUSS_SMOOTH(TEMPORARY(img_tmp), fwhm2sig*fwhm, WIDTH=(4*fwhm < n_xpix0 < n_ypix0), /EDGE_TRUNCATE)
AVGSDV, img_tmp, avg, rms, KAPPA=3 ; compute RMS
img_tmp -= MIN(img_tmp)

; If EDGE is set, set pixels within fwhm of the edge to MIN values
IF KEYWORD_SET(EDGE) THEN BEGIN
  IF N_ELEMENTS(edge) EQ 2 THEN BEGIN
    edge[0] = edge[0] < n_xpix0
    edge[1] = edge[1] < n_ypix0
    n_xchan = CEIL(n_xpix0/edge[0]) & FOR i = 0, n_xchan-1 DO img_tmp[[i*edge[0]+INDGEN(FIX(fwhm)),(i+1)*edge[0]-INDGEN(FIX(fwhm))-1],*] = 0
    n_ychan = CEIL(n_ypix0/edge[1]) & FOR i = 0, n_ychan-1 DO img_tmp[*,[i*edge[1]+INDGEN(FIX(fwhm)),(i+1)*edge[1]-INDGEN(FIX(fwhm))-1]] = 0
  ENDIF ELSE MESSAGE, 'The edge keyword must have 2 elements!'
ENDIF

; First, basic way to find the beam by searching for the maximum of the array
n_star  = n_beam
img_max = MAX(img_tmp, idx_max)
idxmax  = ARRAY_INDICES(img_tmp, idx_max[0])
xcen0   = [idxmax[0]]
ycen0   = [idxmax[1]]
flux    = img_max
IF n_beam GT 1 THEN BEGIN
  ; First set pixels near first maximum to minimum
  img_tmp[((xcen0-2*fwhm)>0):((xcen0+2*fwhm)<(n_xpix0-1)),((ycen0-2*fwhm)>0):((ycen0+2*fwhm)<(n_ypix0-1))] = MIN(img_tmp)
  img_max= MAX(img_tmp, idx_max)
  idxmax = ARRAY_INDICES(img_tmp, idx_max[0])
  xcen0  = [xcen0,idxmax[0]]
  ycen0  = [ycen0,idxmax[1]]
  flux   = [flux,img_max]
ENDIF

; If it doesn't work, we use the FIND function (it is used by default if n_bean is greater than 1)
; Find the beams by looping over the threshold until at least n_beam beams are found.
hmin     = (n_rms*rms < img_max) + rms ; Threshold for point source detection (rms is subtracted later)
roundlim = [-0.5,0.5]                  ; 2 element vector giving low and high cutoff for the roundness statistic (Default: [-1.0,1.0]).
sharplim = [-0.1,0.5]                  ; 2 element vector giving low and high cutoff for the sharpness statistic (Default: [0.2,1.0])
search_again:                          ; jumping point if the star is there but not found (this will try again with a lower threshold)
IF KEYWORD_SET(do_find) THEN BEGIN
  flx_thr = 0.2
  n_star  = 0                          ; will loop through until this is equal to n_beam 
  WHILE n_star LT n_beam DO BEGIN    
    FIND, img_tmp, xcen0, ycen0, flux, sharp, round, hmin, fwhm, roundlim, sharplim, /SILENT
    ; Remove bad position
    IF N_ELEMENTS(x_bad) GT 0 THEN BEGIN
      n_bad = N_ELEMENTS(x_bad)
      FOR ib = 0, n_bad-1 DO BEGIN
        xcen0 = xcen0[WHERE(xcen0 NE x_bad[ib])]
        ycen0 = ycen0[WHERE(xcen0 NE x_bad[ib])]
        flux  = flux[WHERE(xcen0 NE x_bad[ib])]
      ENDFOR
    ENDIF
    ; Number of valid beams
    IF MAX(xcen0) NE 0 THEN n_star = N_ELEMENTS(xcen0)
    ; If too many beams, decide how many to keep (those above average - RMS)
    IF n_star GT 2 THEN BEGIN
      AVGSDV, flux, avg_f, rms_f
      idx_tmp = WHERE(flux GE avg_f - rms_f, n0) 
      IF n0 EQ n_beam THEN n_keep = n0 ELSE n_keep = N_ELEMENTS(xcen0)
      i_iter += 1
      hmin -= rms  ; decrease threshold for next iteration
    ENDIF ELSE BEGIN
      ; If h_min becomes negative, redefine sharplim/roundlim and re-init hmin
      IF hmin LT 0.1*rms THEN BEGIN
        hmin     = n_rms*rms  
        n_keep   = N_ELEMENTS(xcen0)        
        sharplim *= 1.5
        roundlim *= 1.5
        i_iter   += 1
        IF i_iter GT max_iter THEN BEGIN
          PRINT, 'Star not found after ' + STRING(i_iter, FORMAT='(i0)') + ' iterations. Skipped this nod'
          beam_pos[*] = 0
          GOTO, skip_search
        ENDIF
      ENDIF ELSE hmin -= rms
    ENDELSE
  ENDWHILE
ENDIF

; Activate the FIND function
do_find = 1

; Pick the highest flux star if several founds
IF n_star GT n_keep THEN BEGIN
  flx_srt = REVERSE(SORT(flux))
  idx_max = INDGEN(n_star < n_keep)
  xcen0   = xcen0[flx_srt[idx_max]]
  ycen0   = ycen0[flx_srt[idx_max]]
  n_star  = n_keep
ENDIF

; Plot the image to screen to make sure the good star is found
IF KEYWORD_SET(PLOT) THEN BEGIN
  xrange  = [-0.5,0.5]*n_xpix0*pixsize/1D3
  yrange  = [-0.5,0.5]*n_ypix0*pixsize/1D3
  dim     = GET_SCREEN_SIZE()
  window  = [80, 50, 80 + 0.8*dim[1]*n_xpix0/n_ypix0, 50 + 0.8*dim[1]]
  PLOTXY, /INIT, NWINDOW=30, WINDOW=[0, 0, 120 + 0.8*dim[1]*n_xpix0/n_ypix0, 90 + 0.8*dim[1]]
  LOADCT, 0, /SILENT
  PLOTXY, img_tot, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='Locate star', XTITLE='Angular separation [arcsec]', YTITLE='Angular separation [arcsec]', GRID=0, $
    CHARSIZE=1.0, CHARTHICK=1.5, THICK=2.5, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=window;, INSET_UR='a.'
  LOADCT, 13, /SILENT
  IF xcen0[0] NE 0 AND ycen0[0] NE 0 THEN BEGIN
    PLOTXY, [xrange[0]+xcen0*pixsize/1D3], [yrange[0]+ycen0*pixsize/1D3], COLOR=255, SYMBOL=2, /ADD
    FOR is=0, N_ELEMENTS(xcen0)-1 DO XYOUTS, [xrange[0]+(xcen0[is]+2)*pixsize/1D3], [yrange[0]+(ycen0[is]+2)*pixsize/1D3], STRING(is, FORMAT='(I0)'), COLOR=255
  ENDIF ELSE XYOUTS, 0, 0, 'No beam found!', COLOR=255, CHARSIZE=2.5, CHARTHICK=2.0
ENDIF

; If several beams are found, ask which one to use. If star found very faint, ask whether it's the good one.
; If n_beam beams found, record the position right away or search again if too close from the edge 
IF xcen0[0] NE 0 AND ycen0[0] NE 0 THEN BEGIN
  IF n_star GT n_beam OR MAX(flux) LT flx_thr OR ABS(xcen0[0]-n_xpix0) LT fwhm OR ABS(ycen0[0]-n_ypix0) LT fwhm OR ABS(xcen0[0]) LT fwhm OR ABS(ycen0[0]) LT fwhm THEN BEGIN
    ; Ask user to identify good beams
    idx_star = 0D
    PRINT, 'Found too many objects or too low SNR or too close to the edge. Neeed confirmation.'
    FOR i_beam = 0, n_beam-1 DO BEGIN
      ; FInd the more likely position if AUTO_BEAM is set
      ; If not, ask the user
      IF KEYWORD_SET(AUTO_BEAM) THEN BEGIN
        ; Here we make the assumption that the more likely position is that near the vertical middle of a channel
        idx_star = SORT(SQRT(((xcen0 MOD 128)-64)^2 + ((ycen0 MOD 128)-64)^2))
        beam_pos[i_beam,0] = xcen0[idx_star[i_beam]]
        beam_pos[i_beam,1] = ycen0[idx_star[i_beam]]
      ENDIF ELSE BEGIN
        IF NOT KEYWORD_SET(CURSOR) THEN BEGIN 
          READ, idx_star, PROMPT='Enter the location of beam ' + STRING(i_beam+1, FORMAT='(I0)') + ' (see window,-1 if no star,-2 if star not marked): '  ; Read input from the terminal.
          IF idx_star GE 0 THEN BEGIN
            beam_pos[i_beam,0] = xcen0[idx_star]
            beam_pos[i_beam,1] = ycen0[idx_star]
          ENDIF ELSE BEGIN
            xcen0 = 0 
            IF idx_star EQ -1 THEN GOTO, skip_search
            IF idx_star EQ -2 THEN GOTO, search_again
          ENDELSE
        ENDIF ELSE BEGIN
          PRINT, 'Click on the beam position (middle mouse click)'
          CURSOR, x, y, /DEVICE
          beam_pos[i_beam,0] = (x-window[0])*n_xpix0/(window[2]-window[0])
          beam_pos[i_beam,1] = (y-window[1])*n_ypix0/(window[3]-window[1])
        ENDELSE
      ENDELSE
    ENDFOR
  ENDIF ELSE BEGIN
    ; Check if at least the right number of beams have been found!
    IF n_star EQ n_beam THEN BEGIN
      ; Consistancy check
      IF N_ELEMENTS(xcen0) NE n_beam THEN GOTO, search_again
      ; If n_beam = 2, arrange the beams so that the lower or left beam goes to SX (TO BE IMPROVED)
      ; Here, I compute the nod direction
      idx_srt = SORT(xcen0)
      IF n_beam EQ 2 THEN IF ABS(xcen0[0]-xcen0[1]) LT ABS(ycen0[0]-ycen0[1]) THEN idx_srt = SORT(ycen0)
      ; Parse results
      FOR i_beam = 0, n_beam-1 DO BEGIN
        beam_pos[i_beam,0] = xcen0[idx_srt[i_beam]]
        beam_pos[i_beam,1] = ycen0[idx_srt[i_beam]]
      ENDFOR
    ENDIF ELSE GOTO, search_again
  ENDELSE
ENDIF ELSE GOTO, search_again

; Refine the position of the beam (on the coadded stack) if requested
IF KEYWORD_SET(FIT_MODE) THEN BEGIN
  FOR i_beam = 0, n_beam-1 DO BEGIN
    xcen0   = beam_pos[i_beam,0]
    ycen0   = beam_pos[i_beam,1]
    img_jnk = PSF_FIT(img_tot, xcen0, ycen0, fwhm, FIT_METHOD=fit_mode, FILE=file, MEDIAN=0, OFFSET=offset, ZOOM_RATIO=1, $
                      FWHM_X=xsig, FWHM_Y=ysig, SLOPE=slope, XCEN=xcen, YCEN=ycen, INFO=info, PLOT=plot, /NO_SHIFT)
    IF (SQRT((xcen-xcen0)^2 + (ycen-ycen0)^2) GT 0.25*fwhm OR xsig LT 0.5*fwhm OR ysig LT 0.5*fwhm OR xsig GT 1.9*fwhm OR ysig GT 1.5*fwhm) AND NOT KEYWORD_SET(CURSOR) THEN BEGIN
      ; Try fit_mode of 4, which is more robust but less precise (and slower)
      img_jnk = PSF_FIT(img_tot, xcen0, ycen0, fwhm, FIT_METHOD=4, FILE=file, MEDIAN=0, OFFSET=offset, ZOOM_RATIO=1, $
                        FWHM_X=xsig, FWHM_Y=ysig, SLOPE=slope, XCEN=xcen, YCEN=ycen, INFO=info, PLOT=plot, /NO_SHIFT)
      ; If still not found, search again...
      IF (SQRT((xcen-xcen0)^2 + (ycen-ycen0)^2) GT 0.25*fwhm OR xsig LT 0.5*fwhm OR ysig LT 0.5*fwhm OR xsig GT 1.9*fwhm OR ysig GT 1.5*fwhm) AND NOT KEYWORD_SET(CURSOR) THEN BEGIN   
        IF i_iter GT max_iter THEN BEGIN
          PRINT, 'Star not found after ' + STRING(i_iter, FORMAT='(i0)') + ' iterations. Skipped this nod'
          beam_pos[*] = 0
          GOTO, skip_search
        ENDIF  
        IF N_ELEMENTS(x_bad) GT 0 THEN x_bad = [x_bad, xcen0] ELSE x_bad = xcen0         
        MESSAGE, 'Problem with beam position. Restart.', /CONTINUE
        GOTO, search_again
      ENDIF
    ENDIF
    beam_pos[i_beam,0] = xcen
    beam_pos[i_beam,1] = ycen
    fwhm_fit[i_beam,0] = xsig
    fwhm_fit[i_beam,1] = ysig
    slope_pro[i_beam]  = slope
    IF KEYWORD_SET(PLOT) THEN PLOTXY, [xrange[0]+xcen*pixsize/1D3], [yrange[0]+ycen*pixsize/1D3], COLOR=95, SYMBOL=2, /ADD
  ENDFOR
ENDIF

; End plot
IF KEYWORD_SET(PLOT) THEN BEGIN
  LOADCT, 0, /SILENT
  PLOTXY, /FIN
ENDIF
     
; Jumping point if no star found
skip_search:

; Build output data
data_out = {BEAM_POS: beam_pos, FWHM_FIT: fwhm_fit, SLOPE: slope_pro}

RETURN, data_out
END