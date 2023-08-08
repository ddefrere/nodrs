;+
; NAME: LBTI_IMGRED
; 
; PURPOSE:
;   This function performs standard image reduction procedures (flat, dark, wavelength, ...) and return a data cube of calibrated images.
;
; INPUT:
;   img_all       :  An image data cube
;   hdr_all       :  Corresponding header information
;
; KEYWORDS
;   IMG_BPM       :
;   IMG_DRK       :
;   IMG_FLT       :
;   HDR_DRK       :
;   HDR_FLT       :
;   HDR_DATA      :  On output, a structure with header information
;   LOG_FILE      :  Set this keyword the path of a file where the log will be printed
;   INFO          :  Define the level of information printed to screen:
;                       - 0: completely silent execution;
;                       - 1: minimum level of information;
;                       - 2: nominal level of information;
;                       - 3: debugging use.
;   PLOT          :  Set this keyword to plot the data to eps files
;
; OUTPUT
;   Data cube with the calibrated L0 images
;
; MODIFICATION HISTORY:
;   Version 1.0,  30-APR-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (based on former routine "remove_bck.pro")
;   Version 1.1,  02-MAY-2013, DD: Added keyword NO_GPU
;   Version 1.1,  02-JUN-2013, DD: Added '*_idx' keywords
;   Version 1.2,  08-JUK-2013, DD: Completely new formalism
;   Version 1.3,  02-AUG-2013, DD: Speed-up code by a factor ~100 by avoiding stupid concatenations of large data cubes!
;   Version 1.4,  08-AUG-2013, DD: Included path definition in common 'GLOBAL'
;   Version 1.5,  28-OCT-2013, DD: Improved memory usage
;   Version 1.6,  13-JAN-2014, DD: New keyword formalism based on the incremental version of the pipeline
;   Version 1.7,  08-FEB-2014, DD: Improved speed and move MJD computation to LBTI_READL0.pro
;   Version 1.8,  18-SEP-2014, DD: Removed NO_GPU keyword
;   Version 1.9,  18-FEB-2015, DD: Added bias subtraction
;   Version 2.0,  19-FEB-2015, DD: Added keyword IMG_BPM
;   Version 2.1,  22-DEC-2015, DD: Implemented use of drs.precision
;   Version 2.2,  26-MAR-2016, DD: New way to replace bad pixels
;   Version 2.3,  26-NOV-2017, DD: Improved speed!
;   Version 2.4,  17-APR-2018, DD: Added a limit on the number of bad pixels
;   Version 2.5,  29-JUL-2018, DD: Now save BCK_AVG and BCK_RMS separately for each channel

FUNCTION LBTI_IMGRED, img_all, hdr_all, IMG_BPM=img_bpm, IMG_DRK=img_drk, IMG_FLT=img_flt, HDR_DRK=hdr_drk, HDR_FLT=hdr_flt,$ ; path to data directories              
                      HDR_DATA=hdr_data, $                                                                                    ; Output keywords
                      LOG_FILE=log_file, INFO=info, PLOT=plot                                                                 ; Optional input keywords for output information
                       
; Define operational parameters
COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs
ON_ERROR, 0

; Find the science frames
idx_oth = WHERE(hdr_all.datatype EQ 1 OR hdr_all.datatype EQ 2, n_oth, COMPLEMENT=idx_sci)
IF MIN(idx_sci) EQ -1 THEN BEGIN
  MESSAGE, 'No science frames to calibrate! Skipped', /CONTINUE
  RETURN, -1
ENDIF ELSE n_sci = N_ELEMENTS(idx_sci)

; Create output header and image array 
hdr_data = hdr_all[idx_sci]
IF drs.precision THEN img_all = DOUBLE(img_all[*,*,idx_sci]) $
                 ELSE img_all = FLOAT(img_all[*,*,idx_sci])

; Read number of bad pixels
IF KEYWORD_SET(IMG_BPM) THEN BEGIN
  idx_bpm = WHERE(img_bpm EQ 0, nbp)
  IF nbp LE 0 THEN PRINT, 'WARNING, no bad pixels found in bad pixel map.'
ENDIF

; Loop over the science frames and calibrates
FOR i_s=0, n_sci-1 DO BEGIN
  ; Extract current image
  img_tmp = img_all[*,*,i_s]
  hdr_tmp = hdr_data[i_s]
  
  ; Start image calibration
  ; 1. Remove BIAS (estimated from reference pixels, only for LMIRCam for now)
  ; IF STRTRIM(STRUPCASE(instrum)) EQ 'LMIRCAM' AND n_xpix EQ 1024 AND n_ypix EQ 1024 THEN img_in = LBTI_BIASSUBTRACT(TEMPORARY(img_in), /LMIRCAM)
  
  ; 2. Remove dark frame and bad pixels from the data to estimate the background flux
  IF KEYWORD_SET(IMG_DRK) THEN BEGIN
    ; Select the good dark and subtract it
    idx_drk     = WHERE(ABS(hdr_tmp.int_time[0]-hdr_drk.int_time) EQ MIN(ABS(hdr_tmp.int_time[0]-hdr_drk.int_time)) $
                        AND ABS(hdr_tmp.n_coadd[0]-hdr_drk.n_coadd) EQ MIN(ABS(hdr_tmp.n_coadd[0]-hdr_drk.n_coadd)) AND hdr_tmp.pagain[0] EQ hdr_drk.pagain $
                        AND hdr_tmp.detmode[0] EQ hdr_drk.detmode AND  hdr_tmp.detbias[0] EQ hdr_drk.detbias, na)
    ; If not found, remove the integration time from the list
    IF na LE 0 THEN idx_drk = WHERE(ABS(hdr_tmp.n_coadd[0]-hdr_drk.n_coadd) EQ MIN(ABS(hdr_tmp.n_coadd[0]-hdr_drk.n_coadd)) AND hdr_tmp.pagain[0] EQ hdr_drk.pagain $
                                    AND hdr_tmp.detmode[0] EQ hdr_drk.detmode AND  hdr_tmp.detbias[0] EQ hdr_drk.detbias, na)
    ; Subtract DARK
    IF na GT 0 THEN img_tmp = TEMPORARY(img_tmp) - MEAN(hdr_tmp.n_coadd[0]/hdr_drk.n_coadd[idx_drk])*img_drk[*,*,idx_drk] ELSE MESSAGE, 'WARNING, no dark file found for file ' + hdr_tmp[i_s].filename
  ENDIF
    
  ; 3. Perform flat-fielding
  IF KEYWORD_SET(IMG_FLT) THEN BEGIN
    ; Find the good flat frame
    idx_flt  = WHERE(hdr_tmp.smplmode[0] EQ hdr_flt.smplmode AND ABS(hdr_tmp.int_time[0]-hdr_flt.int_time) EQ MIN(ABS(hdr_tmp.int_time[0]-hdr_flt.int_time)) $
                    AND hdr_tmp.pagain[0] EQ hdr_flt.pagain AND hdr_tmp.detbias[0] EQ hdr_flt.detbias AND hdr_tmp.detmode[0] EQ hdr_flt.detmode, na)
    IF na LT 1 THEN MESSAGE, 'WARNING, no flat file found for file ' + hdr_tmp[i_s].filename
    ; Do the flat field correction
    img_tmp  = TEMPORARY(img_tmp)*img_flt[*,*,idx_flt]
  ENDIF
  
  ; Compute size of image
  n_xpix  = (SIZE(img_tmp))[1]
  n_ypix  = (SIZE(img_tmp))[2]
  
  ; 4. Correct bad pixels
  IF KEYWORD_SET(IMG_BPM) THEN BEGIN
    IF nbp GT 1 AND nbp LT 0.1*n_xpix*n_ypix THEN BEGIN    ; only correct for bqd pixels if less thqn 10% of the pixels qre bqd
      FIXPIX, img_tmp, img_bpm, img_out, NPIX=8, /SILENT
      img_tmp = img_out
    ENDIF ELSE IF nbp GT 0.1*n_xpix*n_ypix AND i_s EQ 0 THEN PRINT, 'WARNING, too many bad pixels!!! Skipped...'
  ENDIF
        
  ; 5. Remove median background estimated from test region (different for each channel!)
  ; Also record estimated background in the header (consider only pixels vertically centered to avoid vignetted regions in full frame mode)
  IF hdr_tmp.instrum EQ 'NOMIC' THEN BEGIN
    pix_off = 5  ; pixel offset from the edges
    box_wdt = 50 ; box width to exlude (centered around the center of the channel)
    n_xchan = n_xpix/2.
    n_ychan = cnf.y_chan
    n_chan  = FIX(n_ypix/n_ychan)
    FOR ic = 0, n_chan-1 DO BEGIN
      ; a. left aprt of NOMIC
      mask = INTARR(n_xpix, n_ypix) ; init mask to 0                     
      mask[pix_off:n_xchan-1-pix_off,ic*n_ychan+pix_off:ic*n_ychan+n_ychan-1-pix_off] = 1 ; set everything further than pix_off from the edges to 1
      mask[(0.5*(n_xchan-box_wdt)):(0.5*(n_xchan+box_wdt)),ic*n_ychan+(0.5*(n_ychan-box_wdt)):ic*n_ychan+(0.5*(n_ychan+box_wdt))] = 0 ; set the central box to 0
      idx = WHERE(mask EQ 1)
      med = MEDIAN(img_tmp[idx])
      rms = STDDEV(img_tmp[idx])
      hdr_data[i_s].bck_avg[ic] = med
      hdr_data[i_s].bck_rms[ic] = rms
      img_tmp[0:n_xchan-1,ic*n_ychan:(ic+1)*n_ychan-1] -= med
      ; b. right part of NOMIC (only if no PRE_CROP)
      IF n_xpix EQ N_ELEMENTS(img_tmp[*,0]) THEN BEGIN
        mask = INTARR(n_xpix, n_ypix) ; init mask to 0                     
        mask[n_xchan+pix_off:2*n_xchan-1-pix_off,ic*n_ychan+pix_off:ic*n_ychan+n_ychan-1-pix_off] = 1 ; set everything further than pix_off from the edges to 1
        mask[n_xchan+(0.5*(n_xchan-box_wdt)):n_xchan+(0.5*(n_xchan+box_wdt)),ic*n_ychan+(0.5*(n_ychan-box_wdt)):ic*n_ychan+(0.5*(n_ychan+box_wdt))] = 0 ; set the central box to 0
        idx = WHERE(mask EQ 1)
        med = MEDIAN(img_tmp[idx])
        rms = STDDEV(img_tmp[idx])
        hdr_data[i_s].bck_avg[n_chan+ic] = med
        hdr_data[i_s].bck_rms[n_chan+ic] = rms
        img_tmp[n_xchan:2*n_xchan-1,ic*n_ychan:(ic+1)*n_ychan-1] -= med
      ENDIF
    ENDFOR
  ENDIF
  
  ; 6. Basic cosmetics
  ;img_tmp = LBTI_CLEANIMG(img_tmp, NOMIC=cnf.nomic, LMIRCAM=cnf.lmircam)
  
  ; 7. Uniformize bad detector regions (now done after background subtraction)
  ;n_pix  = N_ELEMENTS(img_tmp[*,0,0]) 
  ;n_chan = n_pix/128.
  ;FOR i_chan = 1, n_chan-2 DO BEGIN
  ;  FOR i_c = 0, n_pix-1 DO img_tmp[i_chan*128,i_c] = (img_tmp[i_chan*128+1,i_c]+img_tmp[i_chan*128-1,i_c])/2.
  ;ENDFOR
  
  ; 8. Remove remaining bad pixels (now done after background subtraction)
  ;box_width = 10.
  ;Nsigma    = 3.
  ;img_tmp   = SIGMA_FILTER(TEMPORARY(img_tmp), box_width, N_SIGMA=Nsigma, /ALL_PIXELS, ITERATE=iterate, MONITOR=monitor, KEEP_OUTLIERS=keep, RADIUS=radius, $
  ;                         N_CHANGE=nchange, VARIANCE_IMAGE=imvar, DEVIATION_IMAGE=imdev)
    
  ; 9. Remove pickup noise. For now, assume the PSF to be Gaussian (RMS=0.5*FWHM) and cut all frequency beyond 10*(1/RMS).
  ; TO BE BETTER DEFINED AND MODIFIED TO BE USED IN IMAGING MODE
  ;freq_cut = 5./(0.5*cnf.psf_fwhm/cnf.pix_size)
  ;img_null = PICKUP_REMOVAL(img_tmp, FREQ_CUT=freq_cut, PIXSIZE=cnf.pix_size, PERIOD=period, NO_GPU=no_gpu, PLOT=plot)
  
  ; 10. Distortion and wavelength calibration will be implemented here
  ; TBD
  
  ; Record calibrated image
  img_all[0,0,i_s] = img_tmp
ENDFOR

RETURN, img_all
END