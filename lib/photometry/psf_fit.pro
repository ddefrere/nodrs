; --- Function called by MPFIT2DFUN to do weighted-PSF fitting ---
FUNCTION FITPSF, x_coord, y_coord, p
  COMMON PSFFIT, psf_img 
  img_out = p[0] + p[1]*SHIFT(psf_img, p[2], p[3])
  RETURN, img_out
END

; --- Main function starts here ---
FUNCTION PSF_FIT, img_in, xcen0, ycen0, fwhm, FIT_METHOD=fit_method, FILE=file, PSF_FILE=psf_file, MEDIAN=median, NO_SHIFT=no_shift, OFFSET=offset, ZOOM_RATIO=zoom_ratio, $          ; Inputs
                  BCK_FLX=bck_flx, BCK_ERR=bck_err, FWHM_X=xsig, FWHM_Y=ysig, IMG_FIT=img_fit, SLOPE=slope, STATUS=status, STR_FLX=str_flx, STR_ERR=str_err, XCEN=xpos, YCEN=ypos, $  ; Outpus
                  INFO=info, PLOT=plot  

  ; PURPOSE:
  ;   This function computes the centroid of a star located near [xcen0,ycen0] using one of the functions specified by FIT_METHOD 
  ;   and shifts it to [xcen0,ycen0] (if NO_SHIFT is not set). Testing of this routine shows good behavior with FIT_METHOD = 1, 2 
  ;   or 5, with option 5 ~50x longer. FIT_METHOD = 3 and 4 show poor results on high noise images.
  ;
  ; INPUTS:
  ;   img_in     :  Input image cube (ignored if FILE is set)
  ;   xcen0      :  Approximate x position of the star
  ;   ycen0      :  Approximate y position of the star
  ;   fwhm       :  FWHM of the PSF in pixels (ignored if PSF_FILE is set)
  ;
  ; KEYWORDS:
  ;   FIT_METHOD :  Define the centroid fitting technique (set it negative to fit inverse function, e.g. for AGPM data):
  ;                       - 0: no stellar centroid fitting
  ;                       - 1: compute the centroid of a star using a derivative search (using CNTRD function)
  ;                       - 2: compute the stellar centroid by Gaussian fits to marginal X,Y, sums (using GCNTRD function)
  ;                       - 3: compute the stellar centroid by Gaussian fit using MPFIT2DPEAK (more time consuming) 
  ;                       - 4: compute the stellar centroid by Lorentzian fit using MPFIT2DPEAK
  ;                       - 5: compute the stellar centroid by Moffat model fit using MPFIT2DPEAK
  ;                       - 6: compute the stellar centroid by PSF image fit using MPFIT2DFUN
  ;   FILE       :  Set this keyword to the name of the file with the image cube (img_in is ignored in that case)
  ;   PSF_FILE   :  Set this keyword to the name of the file with an image of the PSF (FWHM and FIT_METHOD are ignored in that case)
  ;   MEDIAN     :  If set, the median of the centered image cube will be subtracted from each image
  ;   NO_SHIFT   :  If set, the input images are not shifted to be center around the centroid
  ;   OFFSET     :
  ;   ZOOM_RATIO :
  ;   BCK_FLX    :  On output, the background flux (ignored if FIT_METHOD < 3)
  ;   BCK_ERR    :  On output, the error on the background flux (ignored if FIT_METHOD < 3)
  ;   FWHM_X     :  On output, the best-fit FWHM in the x direction (in pixels, only for FIT_METHOD 3, 4, and 5)
  ;   FWHM_Y     :  On output, the best-fit FWHM in the y direction (in pixels, only for FIT_METHOD 3, 4, and 5)
  ;   IMG_FIT    :  On output, the best-fit PSF image (only if FIT_METHOD = 6)
  ;   SLOPE      :  On output, the slope of the fitted Moffat model (only with fit_method = 5)
  ;   STATUS     :  On output, the status if the fit (0: no problem, 1: failed)
  ;   STR_FLX    :  On output, the stellar flux (ignored if FIT_METHOD < 3)
  ;   STR_ERR    :  On output, the error on the stellar flux (ignored if FIT_METHOD < 3)
  ;   XCEN       :  On output, the fitted x position of the star
  ;   YCEN       :  On output, the fitted y position of the star
  ;   INFO       :  Define the level of information printed to screen:
  ;                       - 0: completely silent execution;
  ;                       - 1: minimum level of information;
  ;                       - 2: nominal level of information;  
  ;   PLOT       :  If set, plot results to screen
  ;   
  ; PRE-REQUISITE:
  ;   MPFIT
  ;   
  ; NOTES:
  ;   fit_method=2 currently gives poor results. fit_method=5 works best but is very time-consuming. fit_method=1 works well too for high-SNR images and much faster.
  ; 
  ; LIMITATION
  ;   As of version 2.2, this routine cannot be used anymore to compute the flux! (only for centroiding)
  ;   
  ; MODIFICATION HISTORY:
  ;   Version 1.0, 02-NOV-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu 
  ;   Version 1.1, 21-NOV-2013, DD: added step size used by MPFIT2DPEAK to improve speed computation
  ;   Version 1.2, 22-NOV-2013, DD: implemented call to shifti.pro used by the AGPM team and provided by O. Absil
  ;   Version 1.3, 23-NOV-2013, DD: implemented the use of negative FIT_METHOD to fit inverse functions (e.g., for AGPM data)
  ;   Version 1.4, 24-NOV-2013, DD: added keywords OFFSET and ZOOM_RATIO
  ;   Version 1.5, 13-DEC-2013, DD: added option 6 for FIT_METHOD and keyword PSF_FILE
  ;   Version 1.6, 16-DEC-2013, DD: added BCK_FLX, BCK_ERR, STR_FLX, and STR_ERR to the output keywords
  ;   Version 1.7, 21-MAR-2014, DD: improved speed for fit_method > 2 and added keyword NO_SHIFT
  ;   Version 1.8, 16-NOV-2014, DD: improved robustness of FIT_METHOD=6 
  ;   Version 1.8, 23-NOV-2015, DD: now compute the FWHM automatically when a PSF_FILE is set
  ;   Version 1.9, 24-NOV-2015, DD: added keyword IMG_FIT
  ;   Version 2.0, 21-JAN-2016, DD: now fit the fwhm too (added keywords FWHM_X and FWHM_Y)
  ;   Version 2.1, 16-OCT-2017, DD: added kewyord status
  ;   Version 2.2, 27-JUL-2017, DD: now convolve the image with a Gaussian before searching for the PSF center (much more robust!!)
   
  ; Sanity checks
  IF NOT KEYWORD_SET(info)       THEN info       = 0
  IF NOT KEYWORD_SET(plot)       THEN plot       = 0
  IF NOT KEYWORD_SET(fit_method) THEN fit_method = 0
  IF NOT KEYWORD_SET(offset)     THEN offset     = [0.,0.]
  IF NOT KEYWORD_SET(zoom_ratio) THEN zoom_ratio = 1.
  
  ; Turn off informational messages
  quiet0 = !quiet
  !quiet = 1
  
  ; If Fit methid is 2, force zoom_ratio of 1 (doesn't work with zoomed image)
  IF fit_method EQ 2 THEN zoom_ratio = 1
  
  ; Convertion factor between FWHM and Gaussian sigma
  fwhm2sig = 1/(2*SQRT(2*ALOG(2D)))
  
  ; Read the file and its extension
  IF KEYWORD_SET(FILE)     THEN img_in = MRDFITS(file, 0, hdr_in, /SILENT)
  IF KEYWORD_SET(PSF_FILE) THEN BEGIN
    psf_in = MRDFITS(psf_file, 0, /SILENT)
    tmp    = GAUSS2DFIT(psf_in, a)
    fwhm_x = a[2]
    fwhm_y = a[3]
    PRINT, 'X and Y FWHM in pixels :', fwhm_x, fwhm_y
    fwhm   = fwhm_x > fwhm_y         ; use maximum of the two
  ENDIF
   
  ; Number of images in the input file
  size   = SIZE(img_in)
  IF size[0] EQ 3 THEN n_img = size[3] ELSE n_img =1
  
  ; Size of extracted image 
  nsub   = 4.*ROUND(fwhm)                                ; Must be at least 4 (initial pixel scale => will be even)
  nsub_z = ROUND(nsub*zoom_ratio)                        ; zoomed
  
  ; Redefine parameters based on the zoom ratio
  fwhm_z  = ROUND(fwhm*zoom_ratio)                        ; Zoomed FWHM
  xcen_z  = 0.5*(nsub_z-1)                                ; center X position after congrid      
  ycen_z  = 0.5*(nsub_z-1)                                ; center Y position after congrid     
  xcen0_z = xcen_z+(xcen0+0.5-ROUND(xcen0))*zoom_ratio    ; center X position after congrid + zoomed offset (this assumes even images)
  ycen0_z = ycen_z+(ycen0+0.5-ROUND(ycen0))*zoom_ratio    ; center Y position after congrid + zoomed offset (this assumes even images)
   
  ; Init arrays
  bck_flx = DBLARR(n_img) & bck_err  = DBLARR(n_img)
  str_flx = DBLARR(n_img) & str_err  = DBLARR(n_img)
  xpos    = DBLARR(n_img) & xpos_err = DBLARR(n_img)
  ypos    = DBLARR(n_img) & ypos_err = DBLARR(n_img)
  xsig    = FLTARR(n_img) & ysig     = FLTARR(n_img)
  slope   = FLTARR(n_img) & status  = INTARR(n_img)
  
  ; Prepare MPFIT fitting for FIT_METHOD = 6
  IF fit_method EQ 7 THEN BEGIN
    ; Init MPFIT2DFUN COMMON
    COMMON PSFFIT, psf_img
    ; Make sure psf_img is properly normalized
    psf_in  = psf_in/MAX(psf_in)
    ;psf_img = psf_img/TOTAL(psf_img)
    ; Congrid PSF image
    psf_dim = SIZE(psf_in)
    psf_img = CONGRID(REFORM(psf_in[ROUND(0.5*psf_dim[1]-0.5*nsub):ROUND(0.5*psf_dim[1]+0.5*nsub-1),ROUND(0.5*psf_dim[2]-0.5*nsub):ROUND(0.5*psf_dim[2]+0.5*nsub-1)]), nsub_z, nsub_z, 1, cubic=-0.5)
    ; Prepare PSF-fitting (4 parmaters here: floor, peak, xcen, and ycen)
    parinfo = REPLICATE({value:0., fixed:0, limited:[0,0], limits:[0.,0.], step:0D}, 4)
    ; Output image
    img_fit = 0*img_in
  ENDIF
   
  ; Loop over the frames
  FOR i_img = 0, n_img-1 DO BEGIN
    ; Extract current image and congrid it around approximate center (necessary to sub-pixel alignement)
    IF zoom_ratio GT 1 THEN img_tmp = CONGRID(EXTRAC(REFORM(img_in[*,*,i_img]), ROUND(xcen0+0.5)-0.5*nsub, ROUND(ycen0+0.5)-0.5*nsub, nsub, nsub), nsub_z, nsub_z, 1, cubic=-0.5, /center) $
                       ELSE img_tmp = EXTRAC(REFORM(img_in[*,*,i_img]), ROUND(xcen0+0.5)-0.5*nsub, ROUND(ycen0+0.5)-0.5*nsub, nsub, nsub)
    
    ; Weight more center pixels 
    r_dist  = SHIFT(DIST(nsub_z,nsub_z), 0.5*nsub_z, 0.5*nsub_z)
    ;img_err = r_dist
    
    ; CONVOLVE image
    img_tmp = GAUSS_SMOOTH(img_tmp, fwhm_z*fwhm2sig, /EDGE_TRUNCATE)
    
    ; Prepare MPFIT fitting for FIT_METHOD = 3, 4, or 5
    IF fit_method EQ 4 OR fit_method EQ 5 OR fit_method EQ 6 THEN BEGIN
      ; Derive image floor (mean of outer region) and peak (just the max)
      img_floor = MEAN(img_tmp[WHERE(r_dist GT 1.5*fwhm_z)])
      img_roof  = MAX(img_tmp)
      img_peak  = img_roof-img_floor
      ; Prepare fitting Moffat profile...the more free parameters, the slower.
      parinfo          = REPLICATE({value:0., fixed:0, limited:[0,0], limits:[0.,0.], step:0D}, 8)
      ; Fixed parameters      
      parinfo[0].fixed = 1 ;& parinfo[0].limits = [img_floor,img_roof] & parinfo[0].limited = [1,1] & parinfo[0].step = (img_roof-img_floor)/1D2
      parinfo[1].fixed = 1 ;& parinfo[1].limits = [0.5,2.]*img_peak    & parinfo[1].limited = [1,1] & parinfo[1].step = 0.1*img_peak
      parinfo[6].fixed = 1 ;& parinfo[6].limits = [0.,!Dpi]            & parinfo[6].limited = [1,1] & parinfo[6].step = 0.01*!Dpi
      ; Parameters to fit
      parinfo[2].fixed = 0 & parinfo[2].limits = [0.2,2.0]*fwhm_z  & parinfo[2].limited = [1,1] & parinfo[2].step = 0.02*fwhm_z
      parinfo[3].fixed = 0 & parinfo[3].limits = [0.2,2.0]*fwhm_z  & parinfo[3].limited = [1,1] & parinfo[3].step = 0.02*fwhm_z
      parinfo[4].fixed = 0 & parinfo[4].limits = [0.,nsub_z-1]     & parinfo[4].limited = [1,1] & parinfo[4].step = 0.1*zoom_ratio  ; will be 0.1 of initial pixel size
      parinfo[5].fixed = 0 & parinfo[5].limits = [0.,nsub_z-1]     & parinfo[5].limited = [1,1] & parinfo[5].step = 0.1*zoom_ratio  ; will be 0.1 of initial pixel size
      parinfo[7].fixed = 0 & parinfo[7].limits = [0.,1D+4]         & parinfo[7].limited = [1,1] & parinfo[7].step = 1.0    
      ; Initial guess
      parinfo[*].value = [img_floor, img_peak, fwhm2sig*fwhm_z, fwhm2sig*fwhm_z, xcen0_z, ycen0_z, 0., 1.]
      estimates        = [img_floor, img_peak, fwhm2sig*fwhm_z, fwhm2sig*fwhm_z, xcen0_z, ycen0_z, 0., 1.]
      ; Remove last item if necessary
      IF fit_method EQ 3 OR fit_method EQ 4 THEN BEGIN
        parinfo   = parinfo[INDGEN(7)]
        estimates = estimates[INDGEN(7)]
      ENDIF
    ENDIF
       
    ; Prepare MPFIT fitting for FIT_METHOD = 7
    IF fit_method EQ 7 THEN BEGIN
      ; Derive image floor and peak
      img_floor = MIN(img_tmp)
      img_roof  = MAX(img_tmp)
      img_peak  = img_roof-img_floor
      ; Parameters to fit
      parinfo[0].fixed = 0 & parinfo[0].limits = [img_floor,img_roof] & parinfo[0].limited = [1,1] & parinfo[0].step = 1;(img_roof-img_floor)/1D2
      parinfo[1].fixed = 0 & parinfo[1].limits = [0.5 ,1.0]*img_peak  & parinfo[1].limited = [1,1] & parinfo[1].step = 1.;0.1*img_peak
      parinfo[2].fixed = 0 & parinfo[2].limits = [-fwhm_z,fwhm_z]     & parinfo[2].limited = [1,1] & parinfo[2].step = 1
      parinfo[3].fixed = 0 & parinfo[3].limits = [-fwhm_z,fwhm_z]     & parinfo[3].limited = [1,1] & parinfo[3].step = 1
      ; Initial guess
      parinfo[*].value = [img_floor+0.01*img_peak, 0.99*img_peak, 0, 0]    ; Cannot use the limits
      estimates        = [img_floor+0.01*img_peak, 0.99*img_peak, 0, 0]    ; Cannot use the limits
    ENDIF
    
    ; If fit method is negative, fit the inverse function
    IF fit_method LT 0 THEN img_tmp = 1./img_tmp
    
    ; Perform the stellar centroid fit
    CASE ABS(fit_method) OF
      1: CNTRD, REFORM(img_tmp), xcen0_z, ycen0_z, x0, y0, fwhm_z , /SILENT
      2: GCNTRD, REFORM(img_tmp), xcen0_z, ycen0_z, x0, y0, fwhm_z, /SILENT;, /KEEPCENTER ;(work better without keepcenter in most cases)
      3: img_fit = GAUSS2DFIT(REFORM(img_tmp), param_out)
      4: img_fit = MPFIT2DPEAK(REFORM(img_tmp), param_out, DOF=dof, ERROR=img_err, PARINFO=parinfo, ESTIMATES=estimates, /GAUSSIAN, $
                               MEASURE_ERRORS=measure_errors, BESTNORM=bestnorm, NFREE=nfree, PERROR=perror, ERRMSG=errmsg)
      5: img_fit = MPFIT2DPEAK(REFORM(img_tmp), param_out, DOF=dof, ERROR=img_err, PARINFO=parinfo, ESTIMATES=estimates, /LORENTZIAN, $
                               MEASURE_ERRORS=measure_errors, BESTNORM=bestnorm, NFREE=nfree, PERROR=perror, ERRMSG=errmsg)
      6: img_fit = MPFIT2DPEAK(REFORM(img_tmp), param_out, DOF=dof, ERROR=img_err, PARINFO=parinfo, ESTIMATES=estimates, /MOFFAT, $
                               MEASURE_ERRORS=measure_errors, BESTNORM=bestnorm, NFREE=nfree, PERROR=perror, ERRMSG=errmsg)
      7: param_out = MPFIT2DFUN('FITPSF', INDGEN(2.*fwhm_z), INDGEN(2.*fwhm_z), REFORM(img_tmp), img_err, estimates, PARINFO=parinfo, DOF=dof, BESTNORM=bestnorm, NFREE=nfree, $
                                PERROR=perror, ERRMSG=errmsg, /QUIET)                          
      ELSE: PRINT, 'Undefined fitting method -- no fitting done'
    ENDCASE
   
    ; Extract results if FIT_METHOD = 3, 4, or 5
    IF fit_method EQ 3 OR fit_method EQ 4 OR fit_method EQ 5 OR fit_method EQ 6 THEN BEGIN
      bck_flx[i_img] = param_out[0] & IF fit_method GT 3 THEN bck_err[i_img] = perror[0]*SQRT(bestnorm/dof)
      str_flx[i_img] = param_out[1] & IF fit_method GT 3 THEN str_err[i_img] = perror[1]*SQRT(bestnorm/dof)
      xsig[i_img]    = SQRT(((param_out[2]/fwhm2sig)^2-fwhm^2) > 0)/zoom_ratio ;& xsig_err = perror[2]*SQRT(bestnorm/dof) ; correct for the GAUSSIAN smoothing
      ysig[i_img]    = SQRT(((param_out[3]/fwhm2sig)^2-fwhm^2) > 0)/zoom_ratio ;& ysig_err = perror[3]*SQRT(bestnorm/dof) ; correct for the GAUSSIAN smoothing   
      x0             = param_out[4] ;& xcen_err = perror[4]*SQRT(bestnorm/dof)
      y0             = param_out[5] ;& ycen_err = perror[5]*SQRT(bestnorm/dof)
      IF fit_method EQ 6 THEN slope[i_img] = param_out[7]
      IF NOT FINITE(ysig^2) OR NOT FINITE(xsig^2) OR NOT FINITE(xsig) OR NOT FINITE(ysig) OR xsig LT 0.5*fwhm_z OR xsig GT 2*fwhm_z OR ysig LT 0.5*fwhm_z OR ysig GT 2*fwhm_z OR ABS(x0-xcen_z) GT fwhm_z OR ABS(y0-ycen_z) GT fwhm_z THEN status[i_img] = 1
    ENDIF
    
    ; Extract results if FIT_METHOD = 6
    IF fit_method EQ 7 THEN BEGIN
      bck_flx[i_img] = param_out[0] & bck_err[i_img] = perror[0]*SQRT(bestnorm/dof)
      str_flx[i_img] = param_out[1] & str_err[i_img] = perror[1]*SQRT(bestnorm/dof)
      x0 = xcen0_z + param_out[2] ;& xcen_err = perror[4]*SQRT(bestnorm/dof)
      y0 = ycen0_z + param_out[3] ;& ycen_err = perror[5]*SQRT(bestnorm/dof)
      img_fit[0,0,i_img] = str_flx[i_img]*psf_in+bck_flx[i_img]
    ENDIF
    
    ; Compute status
    IF x0 LE 0 OR y0 LE 0 THEN status[i_img] = 1
    
    ; Compute shifts to center of the image
    shft_x = (x0-xcen_z)/zoom_ratio + offset[0]
    shft_y = (y0-ycen_z)/zoom_ratio + offset[1]
    
    ; Recompute stellar coordinate in non-zoomed coordinates
    ; Center of image is at ROUND(xcen0), ROUND(ycen0) because of 
    xpos[i_img] = ROUND(xcen0 + 0.5) - 0.5 + shft_x
    ypos[i_img] = ROUND(ycen0 + 0.5) - 0.5 + shft_y

    ; Shift the image to re-center it
    IF NOT KEYWORD_SET(NO_SHIFT) THEN img_in[0,0,i_img] = SHIFTI(img_in[*,*,i_img], -shft_x, -shft_y, /CUBIC)
  ENDFOR
  
  ; Plot Histogram of slopes
  IF fit_method EQ 5 AND n_img GE 100 AND plot GT 1 THEN BEGIN
    binsize = (MAX(slope) - MIN(slope))/50.
    hist    = HISTOGRAM(slope, BINSIZE=binsize, MAX=MAX(slope), MIN=MIN(slope))
    n_bin   = N_ELEMENTS(hist)
    bins    = (MIN(slope) + FINDGEN(n_bin)/(n_bin-1.) * (MAX(slope)-MIN(slope)))
    PLOT, [0], [0], XTITLE='Moffat profile slope', YTITLE='Number of occurence', XRANGE=[MIN(slope),MAX(slope)], YRANGE=[0.,1.2*MAX(hist)], XSTYLE=1, YSTYLE=1, CHARSIZE=charsize
    OPLOT, bins, hist, COLOR=color
  ENDIF
  
  ; Compute and subtract median
  IF KEYWORD_SET(MEDIAN) THEN BEGIN
    img_med = MEDIAN(img_in, DIMENSION=3)
    FOR i_img = 0, n_img-1 DO img_in[*,*,i_img] = img_in[*,*,i_img]-img_med
  ENDIF
  
  ; Save file 
  IF KEYWORD_SET(FILE) THEN BEGIN
    outfile = 'l1img_centered.fits'
    MWRFITS, img_in, outfile, hdr_in, /CREATE, /SILENT
  ENDIF
  
  ; Restore informational messages
  !quiet = quiet0
  
  ; Return centered image
  RETURN, img_in
END
; ---------
PRO TEST_PSFFIT
  ;psf_file = '/Volumes/LaCie/nodrs/results/nomic/l0_fits/2014-11-07/2014-11-07_PSF.fits'
  ;img_file = '/Volumes/LaCie/nodrs/results/nomic/l0_fits/2014-11-07/2014-11-07_hltau_IMG_ALL.fits'
  ;psf_file = '/Volumes/nodrs/nodrs/results/nomic/l1_fits/2015-02-08/2015-02-08_CAL_IMG_SEL-MED.fits'
  ;img_file = '/Volumes/nodrs/nodrs/results/nomic/l1_fits/2015-02-08/2015-02-08_SCI_IMG_SEL.fits'
  ;
  
  ; Test on real image
  img_file = '/Users/denis/Desktop/UT2016-10-16_ID011_SCI_alf_Cep_DIT-53ms_11um_PHOT2_IMG.fits'
  img_in   = READFITS(img_file, /SILENT)
  n_img    = N_ELEMENTS(img_in[0,0,*])
  npix     = (size(img_in))[2]
  fwhm   = 16
  fac    = ( 2.0d* sqrt( 2.0d* aLog(2.0d) ) )
  sig    = fwhm / fac
  tt_rms = 3    ; [pix]
  x0     = 0.5*(npix-1)
  y0     = x0
  
  ; Center image file and fit
  pos_err  = FLTARR(n_img,2)
  fwhm_err = pos_err
  FOR i=0, n_img-1 DO BEGIN
    img_tmp = img_in[*,*,i]
    img_cen = PSF_FIT(img_tmp, 0.5*(npix-1), 0.5*(npix-1), fwhm, FIT_METHOD=3, PSF_FILE=psf_file, ZOOM_RATIO=1, BCK_FLX=bck_flx, BCK_ERR=bck_err, IMG_FIT=img_fit, $
                      STR_FLX=str_flx, STR_ERR=str_err, XCEN=xcen, YCEN=ycen, FWHM_X=xsig, FWHM_Y=ysig)
    img_cen = GAUSS2DFIT(img_tmp, a)
    a[2] *= fac
    a[3] *= fac
    PRINT, 'POS X: ' + STRING(xcen, FORMAT='(F5.2)') + ' ' + STRING(a[4], FORMAT='(F5.2)') + ' POS Y: ' + STRING(ycen, FORMAT='(F5.2)') + ' ' + STRING(a[5], FORMAT='(F5.2)') + $
           ' FWHM X: ' + STRING(xsig, FORMAT='(F5.2)') + ' ' +STRING( a[2], FORMAT='(F5.2)') + ' FWHM X: ' + STRING(ysig, FORMAT='(F5.2)') + ' ' + STRING(a[3], FORMAT='(F5.2)') 
    LOADCT, 0, /silent
    TVSCL, img_tmp;, XRANGE=[0,npix-1], YRANGE=[0,npix-1];, NWINDOW=2
    LOADCT, 13, /silent
    OPLOT, [xcen], [ycen], PSYM=2, COLOR=90;, /ADD
    ;PLOTXY, [a[3]], [a[4]], SYMBOL=3, COLOR=250, /ADD
    WAIT, 0.2
  ENDFOR
  
  ; PSF Gaussian test
  npix   = 400
  n_img  = 10
  fwhm   = 16
  fac    = ( 2.0d* sqrt( 2.0d* aLog(2.0d) ) )
  sig    = fwhm / fac 
  tt_rms = 3    ; [pix]
  img_in = PSF_GAUSSIAN(FWHM=fwhm,npixel=npix,ndimen=2,/double)
  x0     = 0.5*(npix-1)
  y0     = x0
  
  ; Center image file and fit
  pos_err  = FLTARR(n_img,2)
  fwhm_err = pos_err
  FOR i=0, n_img-1 DO BEGIN
    x_r     = tt_rms*RANDOMN(noseed) 
    y_r     = tt_rms*RANDOMN(noseed)
    img_tmp = SHIFTI(img_in[*,*], x_r, y_r)
    img_cen = PSF_FIT(img_tmp, 0.5*(npix-1), 0.5*(npix-1), fwhm, FIT_METHOD=3, PSF_FILE=psf_file, ZOOM_RATIO=1, BCK_FLX=bck_flx, BCK_ERR=bck_err, IMG_FIT=img_fit, $
                      STR_FLX=str_flx, STR_ERR=str_err, XCEN=xcen, YCEN=ycen, FWHM_X=xsig, FWHM_Y=ysig)
    pos_err[i,*]  = [x0+x_r-xcen, y0+y_r-ycen]
    fwhm_err[i,*] = [fwhm-xsig, fwhm-ysig]
    ;tmp           = gauss2Dfit(img_tmp, a)
    ;pos_err[i,*]  = [x0+x_r-a[4], y0+y_r-a[5]]
    ;fwhm_err[i,*] = [sig-a[2], sig-a[3]]
  ENDFOR
  LOADCT, 0, /SILENT
  PLOT, pos_err[*,0]
  LOADCT, 13, /SILENT
  OPLOT, pos_err[*,1], COLOR=255
  PRINT, 'Initial position and RMS : ', x0, y0, tt_rms
  AVGSDV, pos_err[*,0], avg_x, rms_x  & PRINT, 'Retrived position error and RMS : ', avg_x, rms_x
  AVGSDV, pos_err[*,1], avg_y, rms_y  & PRINT, 'Retrived position error and RMS : ', avg_y, rms_y
  AVGSDV, fwhm_err[*,0], avg_x, rms_x & PRINT, 'Retrived FWHM error and RMS : ', avg_x, rms_x 
  AVGSDV, fwhm_err[*,1], avg_y, rms_y & PRINT, 'Retrived FWHM error and RMS : ', avg_y, rms_y 
END