;+
; NAME: LBTI_IMGMED
; 
; PURPOSE:
;   This function sorts the input image cubes according critical parameters in HDR_DATA and performs the median computation
;   of each sub-image cube. As of November 2015, only the config ID is used to sort the images.
;
; INPUT:
;   img_in       :  An image data cube
;
; KEYWORDS
;   HDR_DATA     :  On input (resp. output), a structure with header information corresponding to img_in (resp. img_out)
;   MEAN         :  If set, return the mean image instead of the median
;   PRECISION    :  If set, compute median/mean image in double precision (float by defaut)
;   SCALE        :  If set, scale the output image to one coadd
;   INFO         :  Define the level of information printed to screen:
;                       - 0: completely silent execution;
;                       - 1: minimum level of information;
;                       - 2: nominal level of information;
;                       - 3: debugging use.
;
; OUTPUT
;   Data cube with median combined images
;
; MODIFICATION HISTORY:
;   Version 1.0,  30-APR-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.0,  14-OCT-2013, DD: added wavelength to the critical parameters
;   Version 1.1,  08-FEB-2014, DD: improved speed
;   Version 1.2,  30-OCT-2014, DD: added keyword MEAN
;   Version 1.3,  13-JAN-2015, DD: removed option BAD_PIX
;   Version 1.4,  10-OCT-2015, DD: added PACTCDLY to relevant paramaters
;   Version 1.5,  19-NOV-2015, DD: now use only the config ID
;   Version 1.6,  02-DEC-2015, DD: removed unused common BLOCK
;   Version 1.7,  22-DEC-2015, DD: improved memory usage + added keyword PRECISION
;   Version 1.8,  20-JAN-2016, DD: corrected stupid bug for MEAN=2!!
;   Version 1.9,  26-MAR-2016, DD: added keyword SCALE
;   Version 2.0,  26-MAR-2016, DD: improved memory usage

FUNCTION LBTI_IMGMED, img_in, HDR_DATA=hdr_data, MEAN=mean, PRECISION=precision, SCALE=scale, INFO=info

; Keyword checks
IF NOT KEYWORD_SET(info) THEN info = 0
IF NOT KEYWORD_SET(plot) THEN plot = 0

; Derive the number of config ID
cfg_uniq = hdr_data(UNIQ(hdr_data.cfg_id, SORT(hdr_data.cfg_id))).cfg_id 
n_cfg    = N_ELEMENTS(cfg_uniq) > 1  ; Distinct config IDs

; Inititate median array and header data
n_fr    = N_ELEMENTS(img_in[0,0,*])
n_xpix  = N_ELEMENTS(img_in[*,0,0])
n_ypix  = N_ELEMENTS(img_in[0,*,0])
IF KEYWORD_SET(precision) THEN med_img = DBLARR(n_xpix, n_ypix, n_cfg, /NOZERO) $
                          ELSE med_img = FLTARR(n_xpix, n_ypix, n_cfg, /NOZERO)
med_hdr = REPLICATE(hdr_data[0], n_cfg)

; Loop over the config IDs
FOR i_cfg = 0, n_cfg-1 DO BEGIN
  ; Select the first frame not already used for a previous 'i_med'
  idx_cfg = WHERE(hdr_data.cfg_id EQ cfg_uniq[i_cfg], n_cur)

  ; Compute median array
  IF n_cur GT 1 THEN BEGIN
    IF KEYWORD_SET(MEAN) THEN BEGIN
       ; RESISTANT_MEAN too damn slow
       IF mean EQ 1 THEN med_img[0,0,i_cfg] = MEAN(img_in[*,*,idx_cfg], DIMENSION=3, DOUBLE=precision) 
       IF mean EQ 2 THEN BEGIN
        RESISTANT_MEAN,img_in[*,*,idx_cfg], 10, avg, DIMENSION=3, DOUBLE=precision, /SILENT
        med_img[*,*,i_cfg] = avg
       ENDIF
       IF mean GT 2 OR mean LT 0 THEN MESSAGE, 'Undefined mean mode. Must be 1 or 2.'
    ENDIF ELSE med_img[*,*,i_cfg] = MEDIAN(img_in[*,*,idx_cfg], DIMENSION=3, DOUBLE=precision, /EVEN)
  ENDIF ELSE med_img[*,*,i_cfg] = img_in[*,*,idx_cfg]
  
  ; Compute median header (be careful, only critical parameters are meaningful from now on)
  med_hdr[i_cfg] = hdr_data[idx_cfg[0]]
  
  ; Now scale if requested
  IF KEYWORD_SET(SCALE) THEN med_img[0,0,i_cfg] = med_img[*,*,i_cfg]/med_hdr[i_cfg].n_coadd
ENDFOR

; Store median results in outputs arrays
hdr_data = med_hdr

RETURN, med_img
END