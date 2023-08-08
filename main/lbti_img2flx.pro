;+
; NAME: LBTI_IMG2FLX
; 
; PURPOSE:
;   This is the main function to compute the stellar flux of LBTI images and save the L1 files.
;
; INPUTS:
;   img_in        :  An input image data cube from which to compute the flux
;   hdr_in        :  Input header corresponding to the images
;
; KEYWORDS
;   LOG_FILE      :  Set this keyword the path of a file where the log will be printed
;   INFO          :  Define the level of information printed to screen:
;                       - 0: completely silent execution;
;                       - 1: minimum level of information;
;                       - 2: nominal level of information;
;                       - 3: debugging use.
;   NO_SAVE       :  Set this keyword to prevent the code to save the L1 files
;   OB_IN         :  Set to save the input image/values also with these OB numbers 
;   PLOT          :  Set this keyword to plot the data to eps files
;   XCEN          :  Vector with X positions at which the flux must be computed
;   YCEN          :  Vector with Y positions at which the flux must be computed
;
; MODIFICATION HISTORY:
;   Version 1.0, 24-MAY-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (adapted from former routine 'remove_bck.pro')
;   Version 1.1, 29-MAY-2014, DD: added background bias computation
;   Version 1.2, 10-SEP-2014, DD: revised the position where the bias is computed
;   Version 1.3, 29-OCT-2014, DD: now background bias computation is an option
;   Version 1.4, 30-OCT-2014, DD: minor bug corrected
;   Version 1.5, 03-FEB-2015, DD: now works with only one frame!
;   Version 1.6, 04-APR-2015, DD: added PH2ADU to IMG2FLUX routine call
;   Version 1.7, 23-AUG-2015, DD: corrected major bug in the way the photometric frames are arranged and retrived
;   Version 1.8, 12-SEP-2015, DD: added SKY_LIM to IMG2FLX call
;   Version 1.9, 15-SEP-2015, DD: added call to CORRECT_BIAS.pro
;   Version 2.0, 06-OCT-2015, DD: Added computation of PHASECam effective wavelength
;   Version 2.1, 02-DEC-2015, DD: Now save BCK_IRAD, and BCK_ORAD seperately for each APER_RAD
;   Version 2.2, 17-DEC-2015, DD: Now subtract the best-fit surface to the image before flux computation (best-fit excluding the star)   
;   Version 2.3, 21-MAR-2016, DD: Improved code for bck_orad = -2 
;   Version 2.4, 27-MAR-2016, DD: Improved way to compute nearby position 
;   Version 2.5, 06-JUL-2016, DD: Don't compute BIAS if not enough space around beam center
;   Version 2.6, 07-JUL-2016, DD: Improve outer background computation
;   Version 2.7, 11-JUL-2016, DD: Added keywords XCEN and YCEN
;   Version 2.8, 19-JUL-2016, DD: Optimized version and reject now bad photometry based on SLOPE
;   Version 2.9, 26-OCT-2016, DD: Added diagnostic plots 
;   Version 3.0, 06-MAR-2017, DD: Now use EEID as aperture radius if drs.aper_rad is 0.
;   Version 3.1, 23-MAR-2017, DD: Added keyword SKY_COL to call to IMG2FLX
;   Version 3.2, 04-APR-2017, DD: Target information is now passed through the 'tgt' common
;   Version 3.3, 28-JUL-2017, DD: Implemented background floor mode
;   Version 3.4, 28-JUL-2018, DD: Now add the measured flux floor before background subtraction to the background floor

PRO LBTI_IMG2FLX, img_in, hdr_in, LOG_FILE=log_file, INFO=info, NO_SAVE=no_save, PLOT=plot, OB_IN=ob_in, XCEN=xcen_in, YCEN=ycen_in

; Define operational parameters
COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs, log
ON_ERROR, 0

; Keyword sanity check
IF NOT KEYWORD_SET(LOG_FILE) THEN lun = -1 ELSE lun = log_file
IF NOT KEYWORD_SET(OB_IN)    THEN ob_in = 0
IF NOT KEYWORD_SET(INFO)     THEN info  = 0

; Sanity checks
IF drs.flx_mode LT 2 THEN BEGIN  
  ; Compute EEID to determine the aperture size (use the same for SCI and CAL)
  tgt_cur = tgt[WHERE(tgt.name EQ STRLOWCASE(STRCOMPRESS(hdr_in.header.objname, /remove_all)), ntgt)]
  IF ntgt LE 0 THEN BEGIN
    MESSAGE, 'Target not found!', /CONTINUE
    RETURN
  ENDIF
  IF STRTRIM(hdr_in.header.flag, 2) EQ 'CAL' THEN tgt_sci = tgt[WHERE(tgt.name EQ tgt_cur.calfor)] ELSE tgt_sci = tgt_cur  ; If CAL, find corresponding SCI target
  ; IF tgt_sci.calfor EQ 'none' THEN tgt_sci = tgt   ; What is this line????
  eeid_as  = COMPUTE_EEID(0.5*tgt_sci.ldm, tgt_sci.dist, tgt_sci.temp)
  eeid_pix = ROUND(eeid_as/hdr_in.header.pixscale + cnf.psf_fwhm) 
  ; Define maximum radii
  ; Maximum aperture radius based on beam position  (remove 5 to leave some space for the background annulus, 3 pixels are excluded at the edges and we want at least a width of 2)
  aper_max     = tgt_sci.maxapr-5
  bck_irad_max = aper_max + 1
  bck_orad_max = aper_max + 3                               
  ; If aper_rad is 0, we used the EEID + FWHM as aperture radius
  IF MAX(drs.aper_rad) EQ 0 THEN aper_rad = eeid_pix  > ROUND(0.625*cnf.psf_fwhm) $ ; Use EEID+ FWHM as aperture radius. The second term is the optimum radius for aperture photometry (0.625*lambda/D)
                            ELSE aper_rad = drs.aper_rad  
  ; Now define inner radius (as close as possible to aper radius to minimize background bias but at least at lam/D).
  IF MAX(drs.bck_irad) EQ 0 THEN bck_irad = ((aper_rad+1) > 0.9*cnf.psf_pix) > eeid_pix ELSE bck_irad = drs.bck_irad
  ; Make sure that the inner radius is at least 1 pixel larger than the parture radius
  IF N_ELEMENTS(bck_irad) EQ 1 THEN bck_irad = REPLICATE(bck_irad, N_ELEMENTS(aper_rad))  
  ; Consitency checks
  ; WARNING! This means that for some bright stars, the background annulus can cover the EEID!!!
  aper_rad = aper_rad < aper_max  
  bck_irad = bck_irad > (aper_rad + 1)    
  bck_irad = bck_irad < bck_irad_max
ENDIF 

; Number of frames and pixels
n_img   = N_ELEMENTS(img_in[0,0,*])
n_ypix  = N_ELEMENTS(img_in[0,*,0])
n_xpix  = N_ELEMENTS(img_in[*,0,0])
n_xchan = n_xpix/2.
n_ychan = cnf.y_chan
n_chan  = FIX(n_ypix/n_ychan)

; Keep only closed-loop frames for photometry and Fizeau 
CASE hdr_in.header.obstype  OF
  ; 0: PHOTOMETRY. We request that both AO loops are closed
  0: idx_keep = WHERE(hdr_in.data.dloopon EQ 1 OR hdr_in.data.sloopon EQ 1, n_img, NCOMPLEMENT=n_rej)
  ; 1: FIZEAU. We request that both AO loops are closed
  1: idx_keep = WHERE(hdr_in.data.dloopon EQ 1 AND hdr_in.data.sloopon EQ 1, n_img, NCOMPLEMENT=n_rej)
  ; 2: NULLING. We keep all frames at this point because they can be used for background (this selection is done later now because we want open-loop frames with an OBSTYPE of 2 for the background bias)
  2: n_rej = 0 ; idx_keep = WHERE(hdr_in.data.dloopon EQ 1 AND hdr_in.data.sloopon EQ 1 AND hdr_in.data.pcclosed EQ 1, n_img, NCOMPLEMENT=n_rej) ; (We request that both AO loops and the phase loop are closed)
  ; 3: BACKGROUND. No conditions
  3: n_rej = 0
  ELSE: MESSAGE, 'Unknown OBSTYPE'
ENDCASE

; Check that enough frames survived and extract
IF n_img GT 0 THEN BEGIN
  IF n_rej GT 0 THEN BEGIN
    img_in  = TEMPORARY(img_in[*,*,idx_keep])
    data_in = hdr_in.data[idx_keep]
  ENDIF ELSE data_in = hdr_in.data
ENDIF ELSE BEGIN
  MESSAGE, 'No closed-loop frames to process for this nod', /CONTINUE
  RETURN
ENDELSE

; If set, replace beam position by those in the config file (will be easier for later)
IF drs.xcen[0] GT 0 AND drs.ycen[0] GT 0 THEN BEGIN
  data_in[*].xcen_sx = drs.xcen[0]
  data_in[*].ycen_sx = drs.ycen[0]
  IF N_ELEMENTS(drs.xcen) EQ 2 AND N_ELEMENTS(drs.ycen) EQ 2 THEN BEGIN
    IF drs.xcen[1] GT 0 AND drs.ycen[1] GT 0 THEN BEGIN
      data_in[*].xcen_dx = drs.xcen[1]
      data_in[*].ycen_dx = drs.ycen[1]
    ENDIF
  ENDIF
ENDIF

; Define array with beam positions for each image
; But first compute the total number of beam positions per image (specified in the file data + those defined as input)
n_beam = 0
IF MAX(data_in[*].xcen_sx) NE 0 OR MAX(data_in[*].ycen_sx) NE 0 THEN n_beam = n_beam + 1 
IF MAX(data_in[*].xcen_dx) NE 0 OR MAX(data_in[*].ycen_dx) NE 0 THEN n_beam = n_beam + 1 
IF KEYWORD_SET(xcen_in) THEN n_auxil = N_ELEMENTS(xcen_in) ELSE n_auxil = 0
n_beam   = n_beam + n_auxil
IF n_beam EQ 0 THEN BEGIN
  MESSAGE, 'No beam found for this nod. Skip', /CONTINUE
  RETURN
ENDIF
idx      = 0
ob_all   = FLTARR(n_beam) 
xcen_all = FLTARR(n_img, n_beam)
ycen_all = FLTARR(n_img, n_beam)
IF MAX(data_in[*].xcen_sx) NE 0 OR MAX(data_in[*].ycen_sx) NE 0 THEN BEGIN
  ob_all[0]     = data_in[0].ob_id
  xcen_all[*,0] = data_in[*].xcen_sx
  ycen_all[*,0] = data_in[*].ycen_sx
  idx           = 1
ENDIF
IF MAX(data_in[*].xcen_dx) NE 0 OR MAX(data_in[*].ycen_dx) NE 0 THEN BEGIN
  ;ob_all[1]     = data_in.data.ob_id  ; never used because photometric frames don't use ob_all
  xcen_all[*,1] = data_in[*].xcen_dx
  ycen_all[*,1] = data_in[*].ycen_dx
  idx           = 2
ENDIF
FOR i=0, n_auxil-1 DO BEGIN
  ob_all[i+idx]     = ob_in[i]
  xcen_all[*,i+idx] = xcen_in[i]
  ycen_all[*,i+idx] = ycen_in[i]
ENDFOR

; If null beam, only do closed-loop frames
IF hdr_in.header.obstype EQ 2 THEN BEGIN
  idx_skip = WHERE(hdr_in.data.dloopon NE 1 OR hdr_in.data.sloopon NE 1 OR hdr_in.data.pcclosed NE 1, n_rej, NCOMPLEMENT=n_ok)
  IF n_rej GE 1 THEN xcen_all[idx_skip,0] = 0  ; this will skip flux computation in the loop below (for NULL beam only!!! These frames will be used for background at other position)
  IF n_ok  LE 0 THEN BEGIN
    MESSAGE, 'No closed-loop frames to process for this nod', /CONTINUE
    RETURN
  ENDIF
ENDIF

; If no valid position, skip
IF MAX(xcen_all) EQ 0 OR MAX(ycen_all) EQ 0 THEN BEGIN
  MESSAGE, 'No valid beam for this nod', /CONTINUE
  RETURN
ENDIF

; Running paramaters
n_bin  = drs.n_bin > 1
n_clip = drs.n_clip < (n_xpix > n_ypix)   ; keep at least the larger size in the clipped image (but not larger than the biggest size of the image)
n_pix  = n_clip/n_bin
IF n_clip MOD n_bin NE 0 THEN MESSAGE, 'Binning factor must divide clipping factor'

; Scale the pixel size according to nbin
hdr_in.header.pixscale *= n_bin

; Compute psf's FWHM
psf_rad  = 1.028*hdr_in.header.lam_cen/cnf.tel_diam      ; [rad], psf fwhm 
psf_fwhm = psf_rad*prm.r2m/(hdr_in.header.pixscale*1D+3) ; [pixels], psf fwhm 

; Initiate arrays
n_aper      = N_ELEMENTS(drs.aper_rad)
str_flx     = DBLARR(n_img, n_beam, n_aper)
str_err     = DBLARR(n_img, n_beam, n_aper)
bck_flx     = DBLARR(n_img, n_beam, n_aper)
bck_err     = DBLARR(n_img, n_beam, n_aper)
str_flx2    = DBLARR(n_img, n_beam, n_aper)
str_err2    = DBLARR(n_img, n_beam, n_aper)
bck_flx2    = DBLARR(n_img, n_beam, n_aper)
bck_err2    = DBLARR(n_img, n_beam, n_aper)
star_pos    = FLTARR(n_img, 2, n_beam)
img_sav     = DBLARR(n_pix, n_pix, n_img, n_beam)

; Loop over the images
FOR i_img = 0, n_img-1 DO BEGIN 
  ; Loop over the beam position
  FOR i_beam = 0, n_beam-1 DO BEGIN
    ; Extract xcen/ycen
    xcen = REFORM(xcen_all[i_img,*])
    ycen = REFORM(ycen_all[i_img,*])
    ; Make sure xcen and xcen are defined
    IF xcen[i_beam] NE 0 AND ycen[i_beam] NE 0 THEN BEGIN
      ; Derive max distance to channel edges (NOMIC only)
      x_down = ABS(xcen[i_beam]-0.) - 2          ; (Avoid the first 2 pixels from the horizontal edge)
      x_up   = ABS(xcen[i_beam]-0.5*n_xpix) - 2  ; (Avoid the first 2 pixels from the middle)
      IF xcen[i_beam] GT 0.5*n_xpix THEN BEGIN
        x_down = x_up
        x_up   = ABS(xcen[i_beam]-n_xpix) - 2
      ENDIF 
      y_min  = (ycen[i_beam] MOD 128)
      y_down = y_min - 3            ; -3 because we exclude 3 rows at the bottom (strange behavior, reference pixels?)
      y_up   = ABS(y_min-128) - 3   ; -3 because we exclude 3 rows at the top (strange behavior, reference pixels?)
      ; Make sure they all are within the extracted image
      x_down = x_down<0.5*n_clip & x_up = x_up<(0.5*n_clip-1)  ; -1 because 0 based
      y_down = y_down<0.5*n_clip & y_up = y_up<(0.5*n_clip-1)  ; 0 based
      ; If aperture photometry, make sure the outer background radius is within one channel
      IF drs.flx_mode NE 2 THEN BEGIN
        CASE drs.bck_orad OF 
            ; 0 means the channel edge
            0: BEGIN
              ; if /sky_weigth, we only care about the distance to the vertical edges, else we look in all 4 directions
              IF drs.sky_weight THEN bck_orad = y_down < y_up ELSE bck_orad = x_down < x_up < y_down < y_up
              ; Consistancy check
              IF drs.bck_irad GT MIN(bck_orad)+1 THEN MESSAGE, 'The outer radius of the background region is two small compared to the inner radius'
            END
            ; same number of pixels in the photometric aperture and the background region (ensure at least two pixels wide)
            -1: bck_orad = (CEIL(SQRT(aper_rad^2+bck_irad^2)) > (bck_irad + 2)) < y_down < y_up
            ; this will tell the flux computation routine to use this whole region for the background computation 
            -2: sky_lim  = [x_down, x_up, y_down, y_up]
            ELSE: bck_orad = drs.bck_orad ;IF drs.bck_irad GT drs.bck_orad+1 THEN MESSAGE, 'The outer radius of the background region is two small compared to the inner radius'
        ENDCASE
        ; If sky_col is set, make sure that bck_orad is at least bck_irad + aper_rad (but remain within channel)
        IF drs.sky_weight THEN bck_orad = ( bck_orad > (bck_irad + aper_rad)) < y_down < y_up
        ;IF drs.bck_irad GT MIN(bck_orad) THEN BEGIN
        ;  PRINT, 'Inner background radius too large for file ' + STRING(data_in[i_img].file_id, FORMAT='(I0)') + '. Using ' + STRING(MIN(bck_orad)-1, FORMAT='(I0)') + ' pixels instead.'
        ;  bck_irad = MIN(bck_orad)-1
        ;ENDIF
      ENDIF      
      ; Crop image around center, which is at 0.5*(npix-1) for even images. Hence, the +0.5.
      ;PRINT, i_img, xcen[i_beam], ycen[i_beam], bck_orad, y_down, y_up
      img_tmp = EXTRAC(REFORM(img_in[*,*,i_img]), ROUND(xcen[i_beam]+0.5)-0.5*n_clip, ROUND(ycen[i_beam]+0.5)-0.5*n_clip, n_clip, n_clip)
      ; Clean image
      img_tmp = SIGMA_FILTER(TEMPORARY(img_tmp), 10, N_SIGMA=5, /ALL_PIXELS, ITERATE=0, MONITOR=monitor, N_CHANGE=nchange)
      ; Mitigate background bias (not needed after all)
      ; img_tmp = CORRECT_BIAS(img_tmp, LIMIT=sky_lim, STAR_POS=0, EXCL_WIDTH=50)
      ; Subtract polynomial to mitigate background bias (not needed either)
      ; img_tmp[FIX(0.5*n_clip-x_down):FIX(0.5*n_clip+x_up),FIX(0.5*n_clip-y_down):FIX(0.5*n_clip+y_up)] -= SFIT_MASK(img_tmp[FIX(0.5*n_clip-x_down):FIX(0.5*n_clip+x_up),FIX(0.5*n_clip-y_down):FIX(0.5*n_clip+y_up)], 35, DEGREE=3)
      ; Rebin the image if N_BIN is set
      IF n_bin GT 1 THEN img_tmp = REBIN(img_tmp,n_pix,n_pix) ;ELSE img_tmp = img_tmp
      ; Beam position in extracted array
      x_off = (xcen[i_beam] - (ROUND(xcen[i_beam] + 0.5) - 0.5))/n_bin            ; X offset w/r to closest 4-pixel corner
      y_off = (ycen[i_beam] - (ROUND(ycen[i_beam] + 0.5) - 0.5))/n_bin            ; Y offset w/r to closest 4-pixel corner
      star_pos[i_img,*,i_beam] = [0.5*(n_pix-1) + x_off, 0.5*(n_pix-1) + y_off]   ; 0.5*(n_pix-1) is the middle 4-pixel corner in sub-image
      ; Status is ok, unless changed below
      status = [0]
      ; If photometric frame, compute beam offset w/r to image center (for tip/tilt diagnostic).
      IF hdr_in.header.obstype EQ 0 AND drs.fit_mode GT 0 THEN BEGIN
        img_tmp = PSF_FIT(img_tmp, star_pos[i_img,0,i_beam], star_pos[i_img,1,i_beam], psf_fwhm, FIT_METHOD=drs.fit_mode, FILE=file, MEDIAN=0, OFFSET=offset, ZOOM_RATIO=1, $
                          XCEN=x0, YCEN=y0, FWHM_X=xsig, FWHM_Y=ysig, SLOPE=slope, STATUS=status, INFO=info, PLOT=plot, NO_SHIFT=drs.skip_regis)
        IF drs.skip_regis THEN star_pos[i_img,*,i_beam] = [x0, y0] ELSE star_pos[i_img,*,i_beam] = 0.5*(n_pix-1)*[1, 1]  ; will recenter the photometric aperture on the star in call of IMG2FLX
        IF status[0] EQ 0 THEN BEGIN
          IF i_beam EQ 0 THEN BEGIN
            data_in[i_img].xcen_sx  = (ROUND(xcen[i_beam] + 0.5) - 0.5) + (x0 - 0.5*(n_pix-1))      ; X position in the initial images
            data_in[i_img].ycen_sx  = (ROUND(ycen[i_beam] + 0.5) - 0.5) + (y0 - 0.5*(n_pix-1))      ; Y position in the initial images
            data_in[i_img].fwhmx_sx = xsig                               ; Fitted FWHM in X direction
            data_in[i_img].fwhmy_sx = ysig                               ; Fitted FWHM in Y direction
            IF drs.fit_mode EQ 5 THEN data_in[i_img].slope_sx = slope  
          ENDIF ELSE BEGIN
            data_in[i_img].xcen_dx  = (ROUND(xcen[i_beam] + 0.5) - 0.5) + (x0 - 0.5*(n_pix-1))      ; X position in the initial images
            data_in[i_img].ycen_dx  = (ROUND(ycen[i_beam] + 0.5) - 0.5) + (y0 - 0.5*(n_pix-1))      ; Y position in the initial images
            data_in[i_img].fwhmx_dx = xsig                           ; Fitted FWHM in X direction
            data_in[i_img].fwhmy_dx = ysig                           ; Fitted FWHM in Y direction
            IF drs.fit_mode EQ 5 THEN data_in[i_img].slope_dx = slope 
          ENDELSE
        ENDIF ELSE star_pos[i_img,*,i_beam] = [0,0]
      ENDIF
      ; Compute flux if status OK
      IF status[0] EQ 0 THEN BEGIN
        ;IF hdr_in[i_img].obstype EQ 2 THEN star_pos = [35,75] ELSE star_pos = 0
        good_range = [-16000.,16000.]*MEAN(hdr_in.header.n_coadd)        ; define good range for pixel values (in ADU)
        IMG2FLX, img_tmp, BCK_METHOD=drs.bfl_mode, FLX_METHOD=drs.flx_mode, APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, BCK_CEN=bck_cen, GOOD_RANGE=good_range, PSF_FWHM=psf_fwhm, PHPADU=hdr_in.header.eperadu/cnf.qe, $  
                 SKY_COL=drs.sky_col, SKY_LIM=sky_lim, SKY_WEIGHT=drs.sky_weight, STAR_POS=star_pos[i_img,*,i_beam], NSKY=n_sky, STR_FLX=s0, STR_ERR=s0_err, BCK_FLX=b0, BCK_ERR=b0_err                                                                                                                                                     ; Outputs keywords
        ; Parse results
        str_flx[i_img,i_beam,*] = s0 & str_err[i_img,i_beam,*] = s0_err
        bck_flx[i_img,i_beam,*] = b0 & bck_err[i_img,i_beam,*] = b0_err
        ; Save image (subtract residual background)
        img_sav[*,*,i_img,i_beam] = img_tmp - b0[0]
        ; Estimate bias if needed (only nulling frames). Take symmetric position
        ; Only work if at least 128 columns
        IF drs.bias_estim GT 0 AND n_xpix GT 128 AND hdr_in.header.obstype EQ 2 THEN BEGIN
          ; Derive new center
          x_tmp = ABS(xcen[i_beam]-n_xpix)
          ; Crop image around new center
          img_tmp = EXTRAC(REFORM(img_in[*,*,i_img]), ROUND(x_tmp+0.5)-0.5*n_clip, ROUND(ycen[i_beam]+0.5)-0.5*n_clip, n_clip, n_clip)
          ; Clean image (step 2 here, already applied once in image reduction. Don't use iterate which is quite time consuming)
          img_tmp = SIGMA_FILTER(img_tmp, 10, N_SIGMA=5, /ALL_PIXELS, ITERATE=0, MONITOR=monitor, N_CHANGE=nchange)
          ; Mitigate background bias (doesn't work well)
          ; img_tmp = CORRECT_BIAS(img_tmp, LIMIT=sky_lim, STAR_POS=0, EXCL_WIDTH=50)
          ; sky_lim = [x_up, x_down, bck_orad, 122-bck_orad]
          IF n_bin GT 1 THEN img_tmp = REBIN(img_tmp,n_pix,n_pix)
          ; Beamn position in extracted array
          x_off   = (x_tmp - (ROUND(x_tmp + 0.5) - 0.5))/n_bin                           ; X offset w/r to closest 4-pixel corner
          pos_tmp = [0.5*(n_pix-1) + x_off, 0.5*(n_pix-1) + y_off]   ; 0.5*(n_pix-1) is the middle 4-pixel corner in sub-image
          IMG2FLX, img_tmp, BCK_METHOD=drs.bfl_mode, FLX_METHOD=drs.flx_mode, APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, BCK_CEN=bck_cen, GOOD_RANGE=good_range,  PSF_FWHM=psf_fwhm, PHPADU=adu2ph, $
                   SKY_COL=drs.sky_col, SKY_LIM=sky_lim, SKY_WEIGHT=drs.sky_weight, STAR_POS=pos_tmp, $   ; Input keywords
                   STR_FLX=s0, STR_ERR=s0_err, BCK_FLX=b0, BCK_ERR=b0_err 
          str_flx2[i_img,i_beam,*] = s0 & str_err2[i_img,i_beam,*] = s0_err
          bck_flx2[i_img,i_beam,*] = b0 & bck_err2[i_img,i_beam,*] = b0_err
        ENDIF
        ; Add background flux stored in the header data and computed before background subtraction and before the median of each channel is subtracted
        IF N_ELEMENTS(data_in[i_img].bck_avg) GT 1 THEN BEGIN
          i_chan = FLOOR(xcen[i_beam]/n_xchan)
          j_chan = FLOOR(ycen[i_beam]/n_ychan)
          bck_flx[i_img,i_beam,*] += data_in[i_img].bck_avg[n_chan*i_chan+j_chan]
          IF drs.bias_estim GT 0 AND n_xpix GT 128 AND hdr_in.header.obstype EQ 2 THEN BEGIN
            i_chan2 = FLOOR(x_tmp/n_xchan)
            j_chan2 = FLOOR(ycen[i_beam]/n_ychan)
            bck_flx2[i_img,i_beam,*] += data_in[i_img].bck_avg[n_chan*i_chan2+j_chan2]
          ENDIF
        ENDIF ELSE bck_flx[i_img,i_beam,*] += data_in[i_img].bck_avg  ; backward compatibility
      ENDIF
      ; Enable file saving below
      save_on = 1      
    ENDIF
  ENDFOR
ENDFOR

; Undefine img_in to save memory
UNDEFINE, img_in

; Additional filtering based on fwhm and slope rms
fwhm_sx  = 0.5*((data_in.fwhmx_sx) + (data_in.fwhmy_sx))
fwhm_dx  = 0.5*((data_in.fwhmx_dx) + (data_in.fwhmy_dx))
slrms_sx = data_in.ssloprms 
slrms_dx = data_in.dsloprms 
loop_sx  = data_in.sloopon
loop_dx  = data_in.dloopon
tip_sx   = data_in.xcen_sx*hdr_in.header.pixscale*1D3 & IF N_ELEMENTS(tip_sx) GT 1 THEN tip_sx -= MEDIAN(tip_sx) 
tlt_sx   = data_in.ycen_sx*hdr_in.header.pixscale*1D3 & IF N_ELEMENTS(tlt_sx) GT 1 THEN tlt_sx -= MEDIAN(tlt_sx)
tt_sx    = SQRT(tip_sx^2+tlt_sx^2) 
tip_dx   = data_in.xcen_dx*hdr_in.header.pixscale*1D3 & IF N_ELEMENTS(tip_dx) GT 1 THEN tip_dx -= MEDIAN(tip_dx)
tlt_dx   = data_in.ycen_dx*hdr_in.header.pixscale*1D3 & IF N_ELEMENTS(tlt_dx) GT 1 THEN tlt_dx -= MEDIAN(tlt_dx)
tt_dx    = SQRT(tip_dx^2+tlt_dx^2) 

; Plot data before filtering (filtered data are plotted during the null computation)
IF KEYWORD_SET(PLOT) AND n_img GT 1 THEN BEGIN
  plot_path = pth.result_path + 'diagnostic' + pth.sep + 'ao'  + pth.sep + drs.date_obs + pth.sep + 'nod' + STRING(hdr_in.header.nod_id, FORMAT='(I0)') + pth.sep
  IF NOT FILE_TEST(plot_path) THEN FILE_MKDIR, plot_path  
  fit = 20./1720. 
  !P.FONT = 0  
  thick  = 5.0
  xthick = 5.0
  ythick = xthick
  cthick = 3.5
  csize  = 1.3  
  time   = REFORM(TRANSPOSE(data_in.mjd_obs-MIN(data_in.mjd_obs))*24.*60.*60.)  ; convert MJD to ellapsed seconds
  ; 1. Plot FWHM vs time
  IF MAX(fwhm_sx) NE 0 OR MAX(fwhm_dx) NE 0 THEN BEGIN
    plotname = plot_path + drs.date_obs + '_NOD' + STRING(hdr_in.header.nod_id, FORMAT='(I03)') + '_FWHM.eps'
    PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plotname
    LOADCT, 0, /SILENT
    PLOTXY, time, fwhm_sx, YLOG=ylog, /NEW, XRANGE=[0, MAX(time)], YRANGE=[0,MAX(fwhm_sx)>MAX(fwhm_dx)], TITLE=title, XTITLE="Elasped time [s]" , YTITLE="FWHM [pix]", GRID=0, $
             XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
    LOADCT, 13, /SILENT
    PLOTXY, time, fwhm_sx, /ADD, COLOR=85, THICK=thick
    PLOTXY, time, fwhm_dx, /ADD, COLOR=255, THICK=thick
    PLOTXY, time, loop_sx, /ADD, COLOR=85, THICK=thick
    PLOTXY, time, loop_dx, /ADD, COLOR=255, THICK=thick
    LOADCT, 0, /SILENT
    PLOTXY, /FIN
  ENDIF
  ; 2. Plot SLOPE rms vs time
  IF MAX(slrms_sx) NE 0 OR MAX(slrms_dx) NE 0 THEN BEGIN
    plotname = plot_path + drs.date_obs + '_NOD' + STRING(hdr_in.header.nod_id, FORMAT='(I03)') + '_SLOPE-RMS.eps'
    PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plotname
    LOADCT, 0, /SILENT
    PLOTXY, time, slrms_dx, YLOG=ylog, /NEW, XRANGE=[0, MAX(time)], YRANGE=[0,MAX(slrms_sx)>MAX(slrms_dx)], TITLE=title, XTITLE="Elasped time [s]" , YTITLE="SLOPE RMS", GRID=0, $
            XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
    LOADCT, 13, /SILENT
    PLOTXY, time, slrms_sx, /ADD, COLOR=85, THICK=thick
    PLOTXY, time, slrms_dx, /ADD, COLOR=255, THICK=thick
    PLOTXY, time, loop_sx, /ADD, COLOR=85, THICK=thick
    PLOTXY, time, loop_dx, /ADD, COLOR=255, THICK=thick
    LOADCT, 0, /SILENT
    PLOTXY, /FIN
  ENDIF
  ; 3. Plot tip/tilt/tt vs time
  FOR ib = 0, n_beam-1 DO BEGIN
      IF ib EQ 0 THEN BEGIN
        tip = tip_sx
        tlt = tlt_sx
      ENDIF ELSE BEGIN
        tip = tip_dx
        tlt = tlt_dx    
      ENDELSE
      tt = SQRT(tip^2+tlt^2) 
      plotname = plot_path + drs.date_obs + '_NOD' + STRING(hdr_in.header.nod_id, FORMAT='(I03)')
      IF N_ELEMENTS(time) LE 10 THEN no_histo = 1
      IF MAX(tip) NE 0 THEN PLOTALL, time, REFORM(tip), 0, NAME=plotname, TAG='TIP' + STRING(ib, FORMAT='(I0)'), XTITLE='Elapsed time [s]', YTITLE='Tip [mas]', TITLE=' ', NO_HISTO=no_histo, /BOLD, /EPS
      IF MAX(tlt) NE 0 THEN PLOTALL, time, REFORM(tip), 0, NAME=plotname, TAG='TLT' + STRING(ib, FORMAT='(I0)'), XTITLE='Elapsed time [s]', YTITLE='Tilt [mas]', TITLE=' ', NO_HISTO=no_histo, /BOLD, /EPS
      IF MAX(tt)  NE 0 THEN PLOTALL, time, REFORM(tt),  0, NAME=plotname, TAG='TT'  + STRING(ib, FORMAT='(I0)'), XTITLE='Elapsed time [s]', YTITLE='Tip/tilt [mas]', TITLE=' ', NO_HISTO=no_histo, /BOLD, /EPS
  ENDFOR
ENDIF

; Clean results and save
; Loop over beam positions
IF NOT KEYWORD_SET(NO_SAVE) AND KEYWORD_SET(SAVE_ON) THEN BEGIN
  FOR i_b=0, n_beam-1 DO BEGIN
    ; Only if beam poistion not 0
    IF MAX(star_pos[*,0,i_b]) THEN BEGIN
      ; Store beam ID
      hdr_in.header.beam_id = i_b  
      ; Remove bad frames (including open-loop nulling frames)
      idx_keep = WHERE(star_pos[*,0,i_b] GT 0 AND star_pos[*,1,i_b] GT 0, n_keep)
      ; Remove obvious outliers in aperture flux (5 sigma).     
      flx_ob = REFORM(str_flx[idx_keep,i_b,n_aper-1])    ; Consider only maximum aperture.
      AVGSDV, flx_ob, avg, rms, KAPPA=3
      idx_keep = idx_keep[WHERE(ABS(flx_ob-avg) LE 5.*rms, /NULL)]
      ; Remove obvious outliers in background floor (5 sigma).
      bck_ob = REFORM(bck_flx[idx_keep,i_b,n_aper-1])  
      AVGSDV, bck_ob, avg, rms, KAPPA=3
      idx_keep = idx_keep[WHERE(ABS(bck_ob-avg) LE 5.*rms, /NULL, n_keep)] 
      ; If photometry, remove also images with bad FWHM and too large tip/tilt
      IF hdr_in.header.obstype EQ 0 THEN BEGIN
        ; FWHM
        IF i_b EQ 0 THEN fwhm = fwhm_sx[idx_keep] ELSE fwhm = fwhm_dx[idx_keep]      
        idx_ok = WHERE(FINITE(fwhm) EQ 1, n_ok) ; Make sure there are no Infinite values left
        IF n_ok LE 0 THEN BEGIN
          MESSAGE, 'All fitted FWHM are infinite!!!', /CONTINUE
          RETURN
        ENDIF ELSE fwhm = fwhm[idx_ok]
        AVGSDV, fwhm, avg, rms, KAPPA=3
        idx_keep = idx_keep[idx_ok[WHERE(ABS(fwhm-avg) LE 3*rms, /NULL)]]  ; 1.2 is just a guess, especially because we don't use this criterium for the NULL data
        ; Tip/tilt
        IF i_b EQ 0 THEN tt = tt_sx[idx_keep] ELSE tt = tt_dx[idx_keep]
        idx_keep = idx_keep[WHERE(tt LT 40., /NULL, n_keep)]             ; Hard-coded limit of 40 mas (~4 pixels on NOMIC)
        ob_save  = ob_in 
        IF n_keep LE 0 THEN BEGIN
          MESSAGE, 'No frames survived data filtering. Skipped.', /CONTINUE
          RETURN
        ENDIF
      ENDIF ELSE ob_save = ob_all[i_b]
      ; Total number of rejected frames
      hdr_in.header.n_rej = n_rej + (n_img-n_keep)
      ; Save flux in structure
      flx_out = {APER_RAD: aper_rad, BCK_IRAD: bck_irad, BCK_ORAD: bck_orad, NSKY: n_sky, BCKG_MEAS: REFORM(bck_flx[idx_keep,i_b,*]), BCKG_ERR: REFORM(bck_err[idx_keep,i_b,*]), FLX_TOT:REFORM(str_flx[idx_keep,i_b,*]), FLX_ERR:REFORM(str_err[idx_keep,i_b,*]),$
                BCKG_MEAS2: REFORM(bck_flx2[idx_keep,i_b,*]) , BCKG_ERR2: REFORM(bck_err2[idx_keep,i_b,*]), FLX_TOT2: REFORM(str_flx2[idx_keep,i_b,*]), FLX_ERR2: REFORM(str_err2[idx_keep,i_b,*])}
      ; Save data if requested
      LBTI_SAVEL1IMG, REFORM(img_sav[*,*,idx_keep,i_b]), hdr_in.header, data_in[idx_keep], flx_out, FILE_ID=ob_save
      LBTI_SAVEL1FLX, flx_out, hdr_in.header, data_in[idx_keep], FILE_ID=ob_save
      ; Save idx for later
      IF i_b EQ 0 THEN idx_keep1 = idx_keep
      IF i_b EQ 1 THEN idx_keep2 = idx_keep
    ENDIF
  ENDFOR
ENDIF

; Display image for inspection
; Destroy window if exists
nwindow = 29
IF !D.WINDOW EQ nwindow THEN WDELETE, nwindow

; ; Plot the image to screen to make sure the good star is found
IF KEYWORD_SET(PLOT) THEN BEGIN
  n_show  = 2 < n_beam ; show max 2 beams
  xrange  = [-0.5,0.5]*n_pix*hdr_in.header.pixscale
  yrange  = [-0.5,0.5]*n_pix*hdr_in.header.pixscale
  dim     = GET_SCREEN_SIZE()
  window  = [0.02*dim[0], dim[1]-0.55*dim[1], 0.05*dim[0]+n_show*0.25*dim[1], dim[1]-0.05*dim[1]]
  sub_img = 0.9*0.45*(window[3]-window[1])
  PLOTXY, /INIT, NWINDOW=nwindow, WINDOW=window
  LOADCT, 0, /SILENT
  FOR i=0, n_show-1 DO BEGIN
    ; Show the 8 pixels case
    tmp = MIN(ABS(aper_rad-8), i_aper)
    ; Median frame (weighted by the flux)
    idx_keep = INDGEN(N_ELEMENTS(str_flx[*,0,0]))
    IF i EQ 0 AND N_ELEMENTS(idx_keep1) GT 1 THEN idx_keep = idx_keep1 ;ELSE idx_keep = N_ELEMENTS(str_flx[*,0,0])
    IF i EQ 1 AND N_ELEMENTS(idx_keep2) GT 1 THEN idx_keep = idx_keep2 ;ELSE 
    n_keep = N_ELEMENTS(idx_keep)
    ;AVGSDV, n_keep*str_flx[idx_keep,i,i_aper]^2*(star_pos[idx_keep,0,i])/TOTAL(str_flx[idx_keep,i,i_aper]^2), xcen, tmp, KAPPA=3
    ;AVGSDV, n_keep*str_flx[idx_keep,i,i_aper]^2*(star_pos[idx_keep,1,i])/TOTAL(str_flx[idx_keep,i,i_aper]^2), ycen, tmp, KAPPA=3
    AVGSDV, (star_pos[idx_keep,0,i]), xcen, tmp, KAPPA=3
    AVGSDV, (star_pos[idx_keep,1,i]), ycen, tmp, KAPPA=3
    IF n_keep GT 1 THEN img_med = MEDIAN(img_sav[*,*,idx_keep,i], DIMENSION=3) ELSE img_med = REFORM(img_sav[*,*,idx_keep,i])
    sub_win = [0.08*(window[2]-window[0])+i*0.5*(window[2]-window[0]), 0.54*(window[3]-window[1]), 0.08*(window[2]-window[0])+i*0.5*(window[2]-window[0])+sub_img, 0.54*(window[3]-window[1])+sub_img]
    PLOTXY, img_med, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='OB median -- beam ' + STRING(i+1, FORMAT='(I0)'), XTITLE='Angular separation [arcsec]', YTITLE='Angular separation [arcsec]', GRID=0, $
      CHARSIZE=1.0, CHARTHICK=1.5, THICK=2.5, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=sub_win;, INSET_UR='a.'
    IF KEYWORD_SET(APER_RAD) THEN BEGIN
      LOADCT, 13, /SILENT
      alpha = findgen(100)/99.*2.*!dpi
      PLOTXY, xrange[0] + [xcen]*hdr_in.header.pixscale, yrange[0] + [ycen]*hdr_in.header.pixscale, COLOR=90, SYMBOL=2, /ADD
      PLOTXY, xrange[0] + [xcen+aper_rad[i_aper]*cos(alpha)]*hdr_in.header.pixscale, yrange[0] + [ycen+aper_rad[i_aper]*sin(alpha)]*hdr_in.header.pixscale, THICK=1.5, LINESTYLE=0, COLOR=90, /ADD
      PLOTXY, xrange[0] + [xcen+bck_irad[i_aper]*cos(alpha)]*hdr_in.header.pixscale, yrange[0] + [ycen+bck_irad[i_aper]*sin(alpha)]*hdr_in.header.pixscale, THICK=1.5, LINESTYLE=0, COLOR=250, /ADD
      PLOTXY, xrange[0] + [xcen+bck_orad[i_aper]*cos(alpha)]*hdr_in.header.pixscale, yrange[0] + [ycen+bck_orad[i_aper]*sin(alpha)]*hdr_in.header.pixscale, THICK=1.5, LINESTYLE=0, COLOR=250, /ADD
      LOADCT, 0, /SILENT
    ENDIF
    ; Best null
    tmp   = MIN(str_flx[idx_keep,i,i_aper], i_min)
    i_min = idx_keep[i_min]
    xcen  = star_pos[i_min,0,i]
    ycen  = star_pos[i_min,1,i]
    img_med = REFORM(img_sav[*,*,i_min,i])
    sub_win = [0.08*(window[2]-window[0])+i*0.5*(window[2]-window[0]), 0.05*(window[3]-window[1]), 0.08*(window[2]-window[0])+i*0.5*(window[2]-window[0])+sub_img, 0.05*(window[3]-window[1])+sub_img]
    PLOTXY, img_med, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='Lowest flux -- beam ' + STRING(i+1, FORMAT='(I0)'), XTITLE='Angular separation [arcsec]', YTITLE='Angular separation [arcsec]', GRID=0, $
      CHARSIZE=1.0, CHARTHICK=1.5, THICK=2.5, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=sub_win;, INSET_UR='a.'
    IF KEYWORD_SET(APER_RAD) THEN BEGIN
      LOADCT, 13, /SILENT
      alpha = findgen(100)/99.*2.*!dpi
      PLOTXY, xrange[0] + [xcen]*hdr_in.header.pixscale, yrange[0] + [ycen]*hdr_in.header.pixscale, COLOR=90, SYMBOL=2, /ADD
      PLOTXY, xrange[0] + [xcen+aper_rad[i_aper]*cos(alpha)]*hdr_in.header.pixscale, yrange[0] + [ycen+aper_rad[i_aper]*sin(alpha)]*hdr_in.header.pixscale, THICK=1.5, LINESTYLE=0, COLOR=90, /ADD
      PLOTXY, xrange[0] + [xcen+bck_irad[i_aper]*cos(alpha)]*hdr_in.header.pixscale, yrange[0] + [ycen+bck_irad[i_aper]*sin(alpha)]*hdr_in.header.pixscale, THICK=1.5, LINESTYLE=0, COLOR=250, /ADD
      PLOTXY, xrange[0] + [xcen+bck_orad[i_aper]*cos(alpha)]*hdr_in.header.pixscale, yrange[0] + [ycen+bck_orad[i_aper]*sin(alpha)]*hdr_in.header.pixscale, THICK=1.5, LINESTYLE=0, COLOR=250, /ADD
      LOADCT, 0, /SILENT
    ENDIF
  ENDFOR
  PLOTXY, /FIN
ENDIF

END
