;+
; NAME: LBTI_IMGFLAT
; 
; PURPOSE:
;   Compute flat correction images channel per channel based on input flat and dark images.
;
; INPUTS:
;   img_flt :
;
; KEYWORDS
;
; OUTPUT
;   A data cube with the flat images. The header info replaces "HDR_FLT" on output.
;
; MODIFICATION HISTORY:
;   Version 1.0,  17-FEB-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  17-OCT-2013, DD: implemented flat fielding for LMIRCam
;   Version 1.2,  08-FEB-2014, DD: improved speed
;   Version 1.3,  08-FEB-2014, DD: NOMIC and LMIRCAM where swapped!
;   Version 1.4,  22-DEC-2015, DD: Flat was not saved properly for LMIRCam
;   Version 1.5,  27-MAR-2016, DD: Remove call to SIGMA_FILTER, which is nit necessary to define the FLAT!!!
;   Version 1.6,  25-APR-2016, DD: Better handle LMIRCam's reference rows

FUNCTION LBTI_IMGFLAT, img_flt, HDR_FLT=hdr_flt, IMG_DRK=img_drk, HDR_DRK=hdr_drk

; Number of flat images
n_flt  = N_ELEMENTS(img_flt[0,0,*])
n_xpix = N_ELEMENTS(img_flt[*,0,0])
n_ypix = N_ELEMENTS(img_flt[0,*,0])
 
; Loop over each input image
FOR i_flt = 0, n_flt-1 DO BEGIN 
  ; Extract the flat image and header
  hdr_tmp = hdr_flt[i_flt]
  flt_img = REFORM(img_flt[*,*,i_flt])

  ; If set, find best corresponding DARK
  IF KEYWORD_SET(img_drk) AND KEYWORD_SET(hdr_drk) THEN BEGIN
     idx_drk = WHERE(hdr_tmp.n_coadd[0] EQ hdr_drk.n_coadd AND hdr_tmp.pagain[0] EQ hdr_drk.pagain AND hdr_tmp.pabandw[0] EQ hdr_drk.pabandw AND hdr_tmp.detmode[0] EQ hdr_drk.detmode, na)
     IF na GT 1 THEN drk_img = MEDIAN(img_drk[*,*,idx_drk], DIMENSION=3) ELSE IF na EQ 1 THEN drk_img = img_drk[*,*,idx_drk] ELSE BEGIN
      MESSAGE, 'No corresponding dark found for flat-fielding.', /CONTINUE
      RETURN, 0
     ENDELSE
  ENDIF ELSE drk_img = 0.*flt_img
  
  ; Instrument-dependent computation
  IF STRTRIM(STRUPCASE(hdr_tmp.instrum[0])) EQ 'LMIRCAM' THEN BEGIN
    ; Number of channels, pixels per channel, and reference rows
    n_pix   = 64.
    n_ref   = 4    
    n_chan  = n_xpix/n_pix
   
    ; Remove top and bottom 3 pixels (reference pixels in full frame)
    IF n_ypix EQ 1024 THEN idx_col = n_ref + INDGEN(n_ypix-n_ref) ELSE idx_col = INDGEN(n_ypix)
    
    ; Perform flat-fielding column by column
    FOR i_chan=0,n_chan-1 DO BEGIN
      ; Subtract dark for this column
      idx_chan = INDGEN(n_pix) + i_chan*n_pix
      dif_col  = flt_img[idx_chan,*] - drk_img[idx_chan,*]
      
      ; Remove 0 or NaN pixels
      idx0 = WHERE(dif_col EQ 0. OR FINITE(dif_col) NE 1, n0)
      IF n0 GT 0 THEN dif_col[idx0] = MEAN(dif_col[*,idx_col])
            
      ; Compute flat-fielding frame (remove reference pixels)
      flt_img[idx_chan,*]  = MEAN(dif_col[*,idx_col])/dif_col
    ENDFOR
    
    ; Set top and bottom three rows to 1
    IF n_ypix EQ 1024 THEN BEGIN
      idx_ref = [INDGEN(n_ref),1023-INDGEN(n_ref)]
      flt_img[*,idx_ref] = 1
    ENDIF
  ENDIF ELSE BEGIN
    ; Number of channels and pixels per channel
    n_pix   = 128.
    n_chan  = n_ypix/n_pix
    
    ; Loop over the columns
    n_col   = 2
    FOR i_col = 0, n_col-1 DO BEGIN   
      ; Column index
      idx_col = INDGEN(0.5*n_xpix) + i_col*0.5*n_xpix   
      ; Loop over the channels
      FOR i_chan = 0, n_chan-1 DO BEGIN
        ; Column index and channel images
        idx_chan = INDGEN(n_pix) + i_chan*n_pix
        drk_chan = drk_img[*,idx_chan] & drk_chan = drk_chan[idx_col,*]  
        flt_chan = flt_img[*,idx_chan] & flt_chan = flt_chan[idx_col,*]  
            
        ; Subtract dark for this channel
        dif_chan = flt_chan - drk_chan
  
        ; Remove possible 0 or NaN pixels
        idx0 = WHERE(dif_chan EQ 0. OR FINITE(dif_chan) NE 1, n0)
        IF n0 GT 0 THEN dif_chan[idx0] = MEAN(dif_chan)
                 
        ; Compute flat-fielding frame (remove reference pixels)
        flt_img[idx_col[0],idx_chan[0]] = MEAN(dif_chan)/dif_chan
      ENDFOR
    ENDFOR
  ENDELSE
  
  ; Combined channels
  img_flt[0,0,i_flt] = flt_img
ENDFOR

RETURN, img_flt
END