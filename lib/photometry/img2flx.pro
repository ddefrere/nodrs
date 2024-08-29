;+
; NAME: IMG2FLX
; 
; PURPOSE:
;   This function computes the stellar flux either by aperture photometry or PSF-fitting. Unless STAR_POS is set, the star is assumed to be *centered* in each frame.
;   Output keywords return the stellar flux and the background flux, computed either by aperture photometry or PSF fitting.
;
; INPUTS:
;   img_in        :  An input image to compute the flux
;
; KEYWORDS
;   NSKY          :  On output, number of pixels in the background region (can be a vector)
;   STR_FLX       :  On output, the measured stellar flux
;   STR_ERR       :  On output, the 1-sigma error on the measured flux
;   BCK_FLX       :  On output, the measured background flux after image subtraction if bckg_file is set
;   BCK_ERR       :  On output, the 1-sigma error on the measured background flux after image subtraction if bckg_file is set
;   BCK_METHOD    :  Set this keyword to select the method use for background estimation
;                       - 0: sigma-clipped mean (default)
;                       - 1: median
;   FLX_METHOD    :  Set this keyword to select the method use for flux computation
;                       - 0: aperture photometry (default)
;                       - 1: weighted aperture photometry
;                       - 2: PSF-fitting photometry
;   APER_RAD      :  Set this keyword to the radius of the photmetric aperture (in pixels). It can be an array of aperture radii.
;   BCK_IRAD      :  Set this keyword to the inner radius of the background ring (in pixels). It can be a vector of the same size as APER_RAD.
;   BCK_ORAD      :  Set this keyword to the outer radius of the background ring (in pixels). It can be a vector of the same size as APER_RAD.
;   BCK_CEN       :  Set this keyword to a two-element vector with the center position of the background region (by default concentric region around XCEN,YCEN)
;   GOOD_RANGE    :  Set this keyword to the range of acceptable pixel values
;   PHPADU        :  Set this keyword to the estimated photon per ADU ratio
;   SKY_COL       :  Set this keyword in order to use only pixels in the same column for the sky (+impose same number as in photometric aperture)
;   SKY_LIM       :  Four element vector with the outer limit of the background region in all 4 directions (computed from the star): [x_down, x_up, y_down, y_up].
;   SKY_WEIGHT    :  Set this keyword in order to weight the number of pixels per column in the background region according to the corresponding number of pixels in the photometric aperture. 
;   STAR_POS      :  Set this keyword to a two-element vector with the X and Y coordinates of the star (in pixels)
;   INFO          :  Define the level of information printed to screen:
;                       - 0: completely silent execution;
;                       - 1: minimum level of information;
;                       - 2: nominal level of information;
;                       - 3: debugging use.
;   PLOT          :  Set this keyword to plot the data to eps files
;
; OUTPUT
;   Image cube with the star centered in each image and the residual background removed.
;
; MODIFICATION HISTORY:
;   Version 1.0, 30-APR-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (adapted from former routine 'remove_bck.pro')
;   Version 1.1, 29-MAY-2013, DD: when set to 0, the outer radius of the background ring is now set to the closest distance to the- channel edge
;   Version 1.2, 06-JUL-2013, DD: new keyword definitions
;   Version 1.3, 01-AUG-2013, DD: added keyword 'BCK_CEN' allowing an user-defined position for the background region
;   Version 1.4, 02-AUG-2013, DD: Speed-up code by a factor ~15 by avoiding stupid concatenations of large data cubes!
;   Version 1.5, 08-AUG-2013, DD: Included path definition in common 'GLOBAL'
;   Version 1.6, 13-AUG-2013, DD: Now also save the median-combined image
;   Version 1.7, 28-OCT-2013, DD: Improved memory usage
;   Version 1.8, 21-NOV-2013, DD: Renamed keyword METHOD to FLX_METHOD
;   Version 1.9, 16-DEC-2013, DD: Implemented first version of PSF-fitting photometry (by calling PSF_FIT.PRO)
;   Version 2.0, 21-FEB-2014, DD: Now use only pixels in the same column as the star for background computation
;   Version 2.1, 24-MAY-2014, DD: Added keyword STAR_POS and removed keyword NO_SAVE
;   Version 2.2, 25-MAY-2014, DD: Renamed to IMG2FLX and removed COMMON blocks
;   Version 2.3, 21-DEC-2014, DD: Added keywords XCEN, YCEN, and PSF_FWHM
;   Version 2.4, 04-APR-2015, DD: Added keyword PHPADU
;   Version 2.5, 12-SEP-2015, DD: Added keyword SKY_LIM
;   Version 2.6, 16-SEP-2015, DD: BCK_IRAD and BCK_ORAD can now be a vector
;   Version 2.7, 03-MAY-2016, DD: Corrected bug with FLX_METHOD=1
;   Version 2.8, 10-MAY-2016, DD: Added keyword SKY_WEIGHT
;   Version 2.8, 23-MAR-2017, DD: Added keyword SKY_COL
;   Version 2.9, 23-MAY-2017, DD: Added keyword NSKY
;   Version 3.0, 28-JUL-2017, DD: Implemented background floor mode

PRO IMG2FLX, img_in, BCK_METHOD=bck_method, FLX_METHOD=flx_method, APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, BCK_CEN=bck_cen, GOOD_RANGE=good_range, PHPADU=phpadu, PSF_FWHM=psf_fwhm, $
             SKY_COL=sky_col, SKY_LIM=sky_lim, SKY_WEIGHT=sky_weight, STAR_POS=star_pos, $                                                                                                      ; input keywords
             NSKY=n_sky, STR_FLX=str_flx, STR_ERR=str_err, BCK_FLX=bck_flx, BCK_ERR=bck_err, XCEN=xcen, YCEN=ycen                                                                               ; outputs keywords
      
  ; Keyword check
  IF NOT KEYWORD_SET(FLX_METHOD) THEN bck_method = 0
  IF NOT KEYWORD_SET(FLX_METHOD) THEN flx_method = 0
  IF NOT KEYWORD_SET(PHPADU)     THEN phpadu     = 145./0.4         ; that's the default for nulling in high gain mode
  IF flx_method EQ 0 THEN BEGIN
    IF NOT KEYWORD_SET(APER_RAD)   THEN MESSAGE, 'APER_RAD is a mandatory keyword'
    IF NOT KEYWORD_SET(BCK_IRAD)   THEN MESSAGE, 'BCK_IRAD is a mandatory keyword'
    IF NOT KEYWORD_SET(BCK_ORAD)   THEN MESSAGE, 'BCK_ORAD is a mandatory keyword'
  ENDIF
  
  ; Number of frames and pixels
  n_img  = N_ELEMENTS(img_in[0,0,*])
  n_ypix = N_ELEMENTS(img_in[0,*,0])
  n_xpix = N_ELEMENTS(img_in[*,0,0])
  
  ; Star position (should be centered). In IDL, the center of the image is in the middle of the 4 central pixels
  IF NOT KEYWORD_SET(STAR_POS) THEN xcen = 0.5*(n_xpix-1) ELSE xcen = star_pos[0]
  IF NOT KEYWORD_SET(STAR_POS) THEN ycen = 0.5*(n_ypix-1) ELSE ycen = star_pos[1]

  ; Prepare photometric computation (aperture photometry or PSF-fitting)
  IF flx_method EQ 2 THEN BEGIN
    moffat           = 1
    psf_xrms         = 0.5*psf_fwhm
    psf_yrms         = 0.5*psf_fwhm
    parinfo          = REPLICATE({value:0, fixed:0, limited:[0,0], limits:[0,0]}, 8)
    parinfo[2].fixed = 0 & parinfo[2].limits = [0.8,2.0]*psf_xrms  & parinfo[2].limited = [1,1]
    parinfo[3].fixed = 0 & parinfo[3].limits = [0.8,2.0]*psf_yrms  & parinfo[3].limited = [1,1]
    parinfo[6].fixed = 1 & parinfo[6].limits = [0.,!Dpi]           & parinfo[6].limited = [1,1]
    parinfo[7].fixed = 0 & parinfo[7].limits = [0.,20.]            & parinfo[7].limited = [1,1]
  ENDIF ELSE BEGIN
    ; Initiate arrays
    n_aper  = N_ELEMENTS(aper_rad)
    n_sky   = DBLARR(n_aper)
    str_flx = DBLARR(n_aper)
    str_err = DBLARR(n_aper)
    bck_flx = DBLARR(n_aper)
    bck_err = DBLARR(n_aper)
    ; Sanity checks
    IF N_ELEMENTS(bck_irad) EQ 1 THEN bck_irad = REPLICATE(bck_irad, n_aper) ELSE IF N_ELEMENTS(bck_irad) NE n_aper THEN MESSAGE, 'BCK_IRAD and APER_RAD must have the same size' 
    IF N_ELEMENTS(bck_orad) EQ 1 THEN bck_orad = REPLICATE(bck_orad, n_aper) ELSE IF N_ELEMENTS(bck_orad) NE n_aper THEN MESSAGE, 'BCK_ORAD and APER_RAD must have the same size' 
  ENDELSE
  
  ; Do aperture photometry
  IF flx_method NE 2 THEN BEGIN   
    ; Loop over the image frames
    FOR i_aper = 0, n_aper-1 DO BEGIN
      ; Sanity checks
      bck_irad[i_aper] = (aper_rad[i_aper]+1) > bck_irad[i_aper]
      bck_orad[i_aper] = (bck_irad[i_aper]+1) > bck_orad[i_aper]   
      ; Do aperture photometry
      IF NOT KEYWORD_SET(bck_cen) THEN BEGIN
        ; Because the pixels are correlated per column (at least for NOMIC), don't use APER anymore but define the background region as the pixels directly above and below the the star.
        ;APER, img_in, xcen, ycen, flx_tot0, flx_err0, bck_flx0, bck_err0, phpadu, [aper_rad[i_aper]], [bck_irad,bck_orad], good_range, /FLUX, /SILENT, /EXACT, /MEANBACK
        APER_WEIGHT, img_in, xcen, ycen, flx_tot0, flx_err0, bck_flx0, bck_err0, phpadu, [aper_rad[i_aper]], [bck_irad[i_aper],bck_orad[i_aper]], good_range, CLIPSIG=5, SKYLIM=sky_lim, $
                     NSKY=n_bck, SKY_COL=sky_col, SKY_WEIGHT=sky_weight, WEIGHT=flx_method, FWHM=psf_fwhm, MEANBACK=ABS(bck_method-1),SETSKYVAL=1D-14, /FLUX, /SILENT, /EXACT ;(tests show that /MEANBACK gives the best results)
      ENDIF ELSE BEGIN
        ; Compute the array where each value is its distance to xcen, ycen
        DIST_CIRCLE, dist_map, [n_xpix,n_ypix], xcen, ycen, /DOUBLE
        idx_aper = WHERE(dist_map LE aper_rad, np_aper)
        DIST_CIRCLE, dist_map, [n_xpix,n_ypix], xcen+bck_cen[0], ycen+bck_cen[1], /DOUBLE
        idx_bck  = WHERE(dist_map LE bck_rout, n_bck)
        ; Compute fluxes
        ; Extract input fluxes by aperture photometry and store results. The background and background error are given per pixel, hence one has to multiply by the square root of the aperture area in pixels.
        AVGSDV, img_tmp[idx_bck], bck_flx0, bck_err0
        flx_tot0 = TOTAL(img_tmp[idx_aper]) - TOTAL(img_tmp[idx_bck])
        flx_err0 = bck_err0*SQRT(np_aper)+bck_err0*np_aper/n_bck                     ; estimated error on the total flux (standard deviation + error on the mean value)
      ENDELSE
      ; Store results
      ; If undefined (because of a bad pixel for instance)
      IF FINITE(flx_tot0[0]) NE 0 THEN BEGIN
        n_sky[i_aper]   = n_bck
        str_flx[i_aper] = flx_tot0 & str_err[i_aper] = flx_err0
        bck_flx[i_aper] = bck_flx0 & bck_err[i_aper] = bck_err0   
      ENDIF   
    ENDFOR
  ENDIF ELSE BEGIN
    ; Do weighted-PSF fitting (under progress)
    img_tmp    = PSF_FIT(img_in, xcen, ycen, psf_fwhm, FIT_METHOD=5, FILE=file, PSF_FILE=psf_file, MEDIAN=median, OFFSET=offset, ZOOM_RATIO=5, $
                         BCK_FLX=bck_flx, BCK_ERR=bck_err, STR_FLX=str_flx, STR_ERR=str_err, XCEN=xcen0, YCEN=ycen0, SLOPE=slope, INFO=info, PLOT=plot)
  ENDELSE
END
