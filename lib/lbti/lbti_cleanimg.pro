;+
; NAME: LBTI_CLEANIMAGE
;
; PURPOSE:
;   This procedure performs various cosmetics operations on the input image.
;
; INPUT:
;   img_in        :  The image to "clean"
;
; KEYWORDS
;   LMIRCAM       : Set for LMIRCam
;   NOMIC         : Set for NOMIC
;
; OUTPUT
;   Data cube with the centered input images
;
; MODIFICATION HISTORY:
;   Version 1.0,  08-JUL-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  02-DEC-2015, DD: removed SIGMA_FILTER from here to speed up code (now done in LBTI_IMG2FLX)
;   Version 1.2,  24-DEC-2015, DD: restore SIGMA_FILTER which proved useful to remove blinking bad pixels

FUNCTION LBTI_CLEANIMG, img_in, NOMIC=nomic, LMIRCAM=lmircam

; Compute detector size
n_xpix = N_ELEMENTS(img_in[*,0])
n_ypix = N_ELEMENTS(img_in[0,*])

; NOMIC
IF KEYWORD_SET(NOMIC) THEN BEGIN
  ; Derive the number of verticql channels (subframes will always be an integer value of channels).
  n_chan = n_ypix/128.
  
  ; First take out blinking pixels at the edges
  ;idx = [0, 1, n_xpix-2, n_xpix-1]
  ;FOR i = 0, N_ELEMENTS(idx)-1 DO BEGIN
  ;  AVGSDV, img_in[idx[i],*], avg, rms, KAPPA=5
  ;  idx_bad = WHERE(ABS(img_in[idx[i],*]-avg) GT 3*rms, n_bad, COMPLEMENT=idx_ok)
  ;  IF n_bad GT 0 THEN img_in[idx[i],idx_bad] = MEAN(img_in[idx[i],idx_ok])
  ;ENDFOR
  
  ; Replace each channel edge by the average of neighbor columns + remove median
  img_in[*,0] = 0.5*(img_in[*,1]+img_in[*,2])
  img_in[*,n_ypix-1] = 0.5*(img_in[*,n_ypix-2]+img_in[*,n_ypix-3])
  FOR i_chan = 1, n_chan-1 DO BEGIN
    img_in[*,i_chan*128] = 0.5*(img_in[*,128*i_chan+1]+img_in[*,128*i_chan+2])
    img_in[*,i_chan*128-1] = 0.5*(img_in[*,128*i_chan-2]+img_in[*,128*i_chan-3])
    ;img_in[0:0.5*n_xpix-1,i_chan*128:(i_chan+1)*128-1] -= MEDIAN(img_in[0:0.5*n_xpix-1,i_chan*128:(i_chan+1)*128-1])
    ;img_in[0.5*n_xpix:n_xpix-1,i_chan*128:(i_chan+1)*128-1] -= MEDIAN(img_in[0.5*n_xpix:n_xpix-1,i_chan*128:(i_chan+1)*128-1])
  ENDFOR
  
  ; Subtract ELFN by removing from pixel the mean of its column (not to be used with nulling)
  ;FOR i_chan = 0, n_chan-1 DO BEGIN
  ;  idx_chan = 128*i_chan + INDGEN(128)
  ;  FOR i_x = 0, n_xpix-1 DO img_in[i_x,idx_chan] = img_in[i_x,idx_chan] - MEAN(img_in[i_x,idx_chan])
  ;ENDFOR
    
  ; Number of horizontal channels
  ; As of March 2016, this is not always 2 because of possible pre crops.
  ; So, only do the first column until I implement a way to know where we are 
  ;n_chan = 2.                     ; always 2  
  ;FOR i_chan = 1, n_chan-1 DO img_in[i_chan*n_xpix*0.5,*] = 0.5*(img_in[n_xpix*0.5*i_chan+1,*]+img_in[n_xpix*0.5*i_chan-1,*])
  ;img_in[n_xpix-1,*] = img_in[n_xpix-2,*]
  
  ; Extreme columns are always problematic
  avg                = MEDIAN(img_in[2+INDGEN(n_xpix-4),*])
  img_in[0,*]        = avg;0.5*(img_in[1,*]+img_in[2,*])
  img_in[1,*]        = avg
  img_in[n_xpix-2,*] = avg
  img_in[n_xpix-1,*] = avg;0.5*(img_in[n_xpix-2,*]+img_in[n_xpix-3,*])
  
  ; Bottom left pixels are corrupted and not always captured by bad pixel map
  img_in[0:1,0:1] = 0.33*(img_in[2,2]+img_in[1,2]+img_in[2,1])
  IF n_ypix GT 132 THEN img_in[0:1,128:129] = 0.33*(img_in[2,129]+img_in[1,130]+img_in[2,130])
  
  ; Subtract background bias by removing from pixel the mean of its own row (excluding the region around the star)
;  n_chan = 2
;  n_ref  = 30
;  FOR i_y = 0, n_ypix-1 DO BEGIN
;    FOR i_chan = 0, n_chan-1 DO BEGIN
;      idx_chan = n_xpix/2.*i_chan + INDGEN(n_xpix/2.)
;      idx_ref  = idx_chan[i_chan*(n_xpix/2.-n_ref)+INDGEN(n_ref)]
;      AVGSDV, img_in[idx_ref,i_y], avg_row, rms_row, KAPPA=5.
;      img_in[idx_chan,i_y] = img_in[idx_chan,i_y] - avg_row
;    ENDFOR
;  ENDFOR
  
  ;  General cosmetics of NOMIC (smooth vignetted parts, depending on the frame size)
  IF n_xpix EQ 1024 THEN BEGIN
    ; Vertical stripes
    idx_out           = [INDGEN(128),INDGEN(478)+545]
    img_in[idx_out,*] = MEAN(img_in)
    ; Horizontal stripes
    idx_out           = [INDGEN(356),INDGEN(198)+825] 
    img_in[*,idx_out] = MEAN(img_in)
  ENDIF
  IF n_xpix EQ 512 AND n_ypix EQ 512 THEN BEGIN
    ; Vertical stripes
    idx_out           = [INDGEN(128),INDGEN(128)+384]
    img_in[*,idx_out] = MEAN(img_in)
    ; Horizontal stripes
    idx_out           = [INDGEN(25),256+INDGEN(256)]
    img_in[idx_out,*] = MEAN(img_in)
  ENDIF
  IF n_xpix EQ 512 AND n_ypix EQ 256 THEN BEGIN
    ; Horizontal stripes
    idx_out           = [INDGEN(20),256+INDGEN(256)]
    img_in[idx_out,*] = MEAN(img_in)
  ENDIF
ENDIF 

; LMIRCAM
IF KEYWORD_SET(LMIRCAM) THEN BEGIN
  
  ; Number of vertical channels
  n_pix   = 64              ; number of pixels per channel
  n_edge  = 4               ; number of row of bias pixels (it's 3 at the bottom and 4 at the top actually, but use 4 here)
  n_chan  = n_xpix/n_pix    ; number of channles
  top_fra = 0.20            ; Top fraction of pixel used fr bias subtraction
  bot_fra = 0.20            ; Bottom fraction of pixel used fr bias subtraction
  lef_fra = 0.005
  rig_fra = 0.25

  ; Removed sigma-clipped median of each column
  FOR i_col = 0, n_xpix-1 DO BEGIN
    ; Keep only top and bottom 25% pixels
    top_edge = n_edge + INDGEN(FLOOR(top_fra*n_ypix))
    bot_edge = n_edge + INDGEN(FLOOR(bot_fra*n_ypix))
    idx_keep = [top_edge,n_ypix-1-bot_edge]
    ; Perform the sigma-clipped median on the edge pixels
    MEANCLIP, img_in[i_col,idx_keep], med, rms, CLIPSIG=3, /DOUBLE
    ; Sub-tract median from channel
    img_in[i_col,*] -= med   
  ENDFOR
  
  ; Now remove sigma-clipped median of each channel
  FOR i_chan = 0, n_chan-1 DO BEGIN
     ; Extract channel pixels
     idx_chan = INDGEN(n_pix) + n_pix*i_chan
     pix_chan = img_in[idx_chan,*]
     ; Keep only top and bottom 25% pixels
     top_edge = n_edge + INDGEN(FLOOR(top_fra*n_ypix))
     bot_edge = n_edge + INDGEN(FLOOR(bot_fra*n_ypix))
     idx_keep = [top_edge,n_ypix-1-bot_edge]
     pix_chan = pix_chan[*,idx_keep]
     ; Perform the sigma-clipped median on the edge pixels
     MEANCLIP, pix_chan, med, rms, CLIPSIG=3, /DOUBLE
     ; Sub-tract median from channel
     img_in[idx_chan,*] -= med
  ENDFOR
  
  ; Removed sigma-clipped median of each line
  FOR i_line = 0, n_ypix-1 DO BEGIN
    ; Keep only top and bottom 25% pixels
    top_edge = n_edge + INDGEN(FLOOR(lef_fra*n_xpix))
    bot_edge = n_edge + INDGEN(FLOOR(rig_fra*n_xpix))
    idx_keep = [top_edge,n_xpix-1-bot_edge]
    ; Perform the sigma-clipped median on the edge pixels
    MEANCLIP, img_in[idx_keep,i_line], med, rms, CLIPSIG=3, /DOUBLE
    ; Sub-tract median from channel
    img_in[*,i_line] -= med
  ENDFOR
    
  ; Deal with bad column and smooth bad/vignetted regions
  IF n_xpix EQ 1024 THEN BEGIN
     img_in[0,*]   = 0.5*(img_in[1,*]+img_in[2,*])  
     img_in[715,*] = 0.5*(img_in[713,*]+img_in[718,*])  
     img_in[716,*] = img_in[715,*] 
  ;  idx_out = [INDGEN(5),1023-INDGEN(5)]
  ;  img_in[*,idx_out] = MEAN(img_in)
  ;  idx_out = [INDGEN(5),1023-INDGEN(5)]
  ;  img_in[idx_out,*] = MEAN(img_in)
  ENDIF
  
  ; Set reference pixels to median value
  IF n_ypix EQ 1024 THEN BEGIN
    img_med = MEDIAN(img_in)
    img_in[*,INDGEN(n_edge)]      = img_med
    img_in[*,1023-INDGEN(n_edge)] = img_med
  ENDIF
ENDIF

; Remove obvious outliers
img_in = SIGMA_FILTER(TEMPORARY(img_in), 10, N_SIGMA=5, /ALL_PIXELS, ITERATE=0, MONITOR=monitor, N_CHANGE=nchange)

RETURN, img_in
END