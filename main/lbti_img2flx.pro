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
;   Version 3.5, 27-MAY-2024, DD: Now make sure bck_orad has the same size as bck_irad
;   Version 3.6, 27-AUG-2024, DD: Remove transition frames in the background data
;   Version 3.7, 02-SEP-2024, DD: Added option to disregard the sky annulus
;   Version 3.8, 24-JAN-2025, DD: Corrected bug introduced in version 3.6 for cases where no star was found in nulling mode
;
pro LBTI_IMG2FLX, img_in, hdr_in, log_file = log_file, info = info, no_save = no_save, plot = plot, ob_in = ob_in, xcen = xcen_in, ycen = ycen_in
  compile_opt idl2

  ; Define operational parameters
  common GLOBAL, prm, cnf, wav, tgt, pth, drs, log
  on_error, 0

  ; Keyword sanity check
  if not keyword_set(log_file) then lun = -1 else lun = log_file
  if not keyword_set(ob_in) then ob_in = 0
  if not keyword_set(info) then info = 0

  ; Sanity checks
  if drs.flx_mode lt 2 then begin
    ; Compute EEID to determine the aperture size (use the same for SCI and CAL)
    tgt_cur = tgt[where(tgt.name eq strlowcase(strcompress(hdr_in.header.objname, /remove_all)), ntgt)]
    if ntgt le 0 then begin
      message, 'Target not found!', /continue
      RETURN
    endif
    if strtrim(hdr_in.header.flag, 2) eq 'CAL' then tgt_sci = tgt[where(tgt.name eq tgt_cur.calfor)] else tgt_sci = tgt_cur ; If CAL, find corresponding SCI target
    ; IF tgt_sci.calfor EQ 'none' THEN tgt_sci = tgt   ; What is this line????
    eeid_as = COMPUTE_EEID(0.5 * tgt_sci.ldm, tgt_sci.dist, tgt_sci.temp)
    eeid_pix = round(eeid_as / hdr_in.header.pixscale + cnf.psf_fwhm)
    ; Define maximum radii
    ; Maximum aperture radius based on beam position  (remove 5 to leave some space for the background annulus, 3 pixels are excluded at the edges and we want at least a width of 2)
    aper_max = tgt_sci.maxapr - 5
    bck_irad_max = aper_max + 1
    bck_orad_max = aper_max + 3
    ; If aper_rad is 0, we used the EEID + FWHM as aperture radius
    if max(drs.aper_rad) eq 0 then aper_rad = eeid_pix > round(0.625 * cnf.psf_fwhm) $ ; Use EEID+ FWHM as aperture radius. The second term is the optimum radius for aperture photometry (0.625*lambda/D)
    else aper_rad = drs.aper_rad
    ; Now define inner radius (as close as possible to aper radius to minimize background bias but at least at lam/D).
    if max(drs.bck_irad) eq 0 then bck_irad = ((aper_rad + 1) > 0.9 * cnf.psf_pix) > eeid_pix else bck_irad = drs.bck_irad
    ; Make sure that the inner radius is at least 1 pixel larger than the parture radius
    if n_elements(bck_irad) eq 1 then bck_irad = replicate(bck_irad, n_elements(aper_rad))
    ; Consitency checks
    ; WARNING! This means that for some bright stars, the background annulus can cover the EEID!!!
    aper_rad = aper_rad < aper_max
    bck_irad = bck_irad > (aper_rad + 1)
    bck_irad = bck_irad < bck_irad_max
  endif

  ; Number of frames and pixels
  n_img = n_elements(img_in[0, 0, *])
  n_ypix = n_elements(img_in[0, *, 0])
  n_xpix = n_elements(img_in[*, 0, 0])
  n_xchan = n_xpix / 2.
  n_ychan = cnf.y_chan
  n_chan = fix(n_ypix / n_ychan)

  ; Keep only closed-loop frames for photometry and Fizeau
  case hdr_in.header.obstype of
    ; 0: PHOTOMETRY. We request that both AO loops are closed
    0: idx_keep = where(hdr_in.data.dloopon eq 1 or hdr_in.data.sloopon eq 1, n_img, ncomplement = n_rej)
    ; 1: FIZEAU. We request that both AO loops are closed
    1: idx_keep = where(hdr_in.data.dloopon eq 1 and hdr_in.data.sloopon eq 1, n_img, ncomplement = n_rej)
    ; 2: NULLING. We keep all frames at this point because they can be used for background (this selection is done later now because we want open-loop frames with an OBSTYPE of 2 for the background bias)
    2: n_rej = 0 ; idx_keep = WHERE(hdr_in.data.dloopon EQ 1 AND hdr_in.data.sloopon EQ 1 AND hdr_in.data.pcclosed EQ 1, n_img, NCOMPLEMENT=n_rej) ; (We request that both AO loops and the phase loop are closed)
    ; 3: BACKGROUND. No conditions
    3: n_rej = 0
    else: message, 'Unknown OBSTYPE'
  endcase

  ; Check that enough frames survived and extract
  if n_img gt 0 then begin
    if n_rej gt 0 then begin
      img_in = temporary(img_in[*, *, idx_keep])
      data_in = hdr_in.data[idx_keep]
    endif else data_in = hdr_in.data
  endif else begin
    message, 'No closed-loop frames to process for this nod', /continue
    RETURN
  endelse

  ; If set, replace beam position by those in the config file (will be easier for later)
  if drs.xcen[0] gt 0 and drs.ycen[0] gt 0 then begin
    data_in[*].xcen_sx = drs.xcen[0]
    data_in[*].ycen_sx = drs.ycen[0]
    if n_elements(drs.xcen) eq 2 and n_elements(drs.ycen) eq 2 then begin
      if drs.xcen[1] gt 0 and drs.ycen[1] gt 0 then begin
        data_in[*].xcen_dx = drs.xcen[1]
        data_in[*].ycen_dx = drs.ycen[1]
      endif
    endif
  endif

  ; Define array with beam positions for each image
  ; But first compute the total number of beam positions per image (specified in the file data + those defined as input)
  n_star = 0
  if max(data_in[*].xcen_sx) ne 0 or max(data_in[*].ycen_sx) ne 0 then n_star = n_star + 1
  if max(data_in[*].xcen_dx) ne 0 or max(data_in[*].ycen_dx) ne 0 then n_star = n_star + 1
  if keyword_set(xcen_in) then n_auxil = n_elements(xcen_in) else n_auxil = 0 ; Contain the position of the star in the next nod. We need this to compute the background flux at the same position here
  if n_star eq 0 and n_auxil eq 0 then begin
    message, 'No beam found for this nod. Skip', /continue
    RETURN
  endif
  if hdr_in.header.obstype eq 2 and n_star eq 0 then n_star = 1 ; Make sure there is at least one beam in nulling mode. The auxiliary fluwx has to go in the second column (sometimes, there is no nulling beam but we need this nod for background)
  n_beam = n_star + n_auxil
  ob_all = fltarr(n_beam)
  xcen_all = fltarr(n_img, n_beam)
  ycen_all = fltarr(n_img, n_beam)
  if max(data_in[*].xcen_sx) ne 0 or max(data_in[*].ycen_sx) ne 0 then begin
    ob_all[0] = data_in[0].ob_id
    xcen_all[*, 0] = data_in[*].xcen_sx
    ycen_all[*, 0] = data_in[*].ycen_sx
  endif
  if max(data_in[*].xcen_dx) ne 0 or max(data_in[*].ycen_dx) ne 0 then begin
    ; ob_all[1]     = data_in.data.ob_id  ; never used because photometric frames don't use ob_all
    xcen_all[*, 1] = data_in[*].xcen_dx
    ycen_all[*, 1] = data_in[*].ycen_dx
  endif
  for i = 0, n_auxil - 1 do begin
    ob_all[i + n_star] = ob_in[i]
    xcen_all[*, i + n_star] = xcen_in[i]
    ycen_all[*, i + n_star] = ycen_in[i]
  endfor

  ; If null beam, only do closed-loop frames for the NULL beam. For the background beam, skip the first 5 frames which are often transition frames (when the loop is open)
  if hdr_in.header.obstype eq 2 then begin
    ; Only keep closed-loop frames for nulling beam
    idx_skip = where(hdr_in.data.dloopon ne 1 or hdr_in.data.sloopon ne 1 or hdr_in.data.pcclosed ne 1, n_rej, ncomplement = n_ok)
    if n_rej ge 1 then xcen_all[idx_skip, 0] = 0 ; this will skip flux computation in the loop below (for NULL beam only!!! These frames will be used for background at other position)
    if n_ok le 0 then begin
      message, 'No closed-loop null frames to process for this nod', /continue
      RETURN
    endif
    ; Remove first ten open-loop transition frames for background beam. In recent data, we continue data acquisition when the loop is being closed and, sometimes, the star is still moving to the other quadrant.
    ; This makes sure that the star has completely moved.
    idx10 = indgen(10)
    idx_skip = where(hdr_in[idx10].data.dloopon ne 1 or hdr_in.data[idx10].sloopon ne 1 or hdr_in[idx10].data.pcclosed ne 1, n_rej)
    if n_rej ge 1 then xcen_all[idx_skip, 1] = 0 ; this will skip flux computation in the loop below (for the bckg beam only!!!)
  endif

  ; If no valid position, skip
  if max(xcen_all) eq 0 or max(ycen_all) eq 0 then begin
    message, 'No valid beam for this nod', /continue
    RETURN
  endif

  ; Running paramaters
  n_bin = drs.n_bin > 1
  n_clip = drs.n_clip < (n_xpix > n_ypix) ; keep at least the larger size in the clipped image (but not larger than the biggest size of the image)
  n_pix = n_clip / n_bin
  if n_clip mod n_bin ne 0 then message, 'Binning factor must divide clipping factor'

  ; Scale the pixel size according to nbin
  hdr_in.header.pixscale *= n_bin

  ; Compute psf's FWHM
  psf_rad = 1.028 * hdr_in.header.lam_cen / cnf.tel_diam ; [rad], psf fwhm
  psf_fwhm = psf_rad * prm.r2M / (hdr_in.header.pixscale * 1d+3) ; [pixels], psf fwhm

  ; Initiate arrays
  n_aper = n_elements(drs.aper_rad)
  str_flx = dblarr(n_img, n_beam, n_aper)
  str_err = dblarr(n_img, n_beam, n_aper)
  bck_flx = dblarr(n_img, n_beam, n_aper)
  bck_err = dblarr(n_img, n_beam, n_aper)
  str_flx2 = dblarr(n_img, n_beam, n_aper)
  str_err2 = dblarr(n_img, n_beam, n_aper)
  bck_flx2 = dblarr(n_img, n_beam, n_aper)
  bck_err2 = dblarr(n_img, n_beam, n_aper)
  star_pos = fltarr(n_img, 2, n_beam)
  img_sav = dblarr(n_pix, n_pix, n_img, n_beam)

  ; Loop over the images
  for i_img = 0, n_img - 1 do begin
    ; Loop over the beam position
    for i_beam = 0, n_beam - 1 do begin
      ; Extract xcen/ycen
      xcen = reform(xcen_all[i_img, *])
      ycen = reform(ycen_all[i_img, *])
      ; Make sure xcen and xcen are defined
      if xcen[i_beam] ne 0 and ycen[i_beam] ne 0 then begin
        ; Derive max distance to channel edges (NOMIC only)
        x_down = abs(xcen[i_beam] - 0.) - 2 ; (Avoid the first 2 pixels from the horizontal edge)
        x_up = abs(xcen[i_beam] - 0.5 * n_xpix) - 2 ; (Avoid the first 2 pixels from the middle)
        if xcen[i_beam] gt 0.5 * n_xpix then begin
          x_down = x_up
          x_up = abs(xcen[i_beam] - n_xpix) - 2
        endif
        y_min = (ycen[i_beam] mod 128)
        y_down = y_min - 3 ; -3 because we exclude 3 rows at the bottom (strange behavior, reference pixels?)
        y_up = abs(y_min - 128) - 3 ; -3 because we exclude 3 rows at the top (strange behavior, reference pixels?)
        ; Make sure they all are within the extracted image
        x_down = x_down < 0.5 * n_clip
        x_up = x_up < (0.5 * n_clip - 1) ; -1 because 0 based
        y_down = y_down < 0.5 * n_clip
        y_up = y_up < (0.5 * n_clip - 1) ; 0 based
        ; If aperture photometry, make sure the outer background radius is within one channel
        if drs.flx_mode ne 2 then begin
          case drs.bck_orad of
            ; 0 means the channel edge
            0: begin
              ; if /sky_weigth, we only care about the distance to the vertical edges, else we look in all 4 directions
              if drs.sky_weight then bck_orad = y_down < y_up else bck_orad = x_down < x_up < y_down < y_up
              ; Consistancy check
              if drs.bck_irad gt min(bck_orad) + 1 then message, 'The outer radius of the background region is two small compared to the inner radius'
            end
            ; same number of pixels in the photometric aperture and the background region (ensure at least two pixels wide)
            - 1: bck_orad = (ceil(sqrt(aper_rad ^ 2 + bck_irad ^ 2)) > (bck_irad + 2)) < y_down < y_up
            ; this will tell the flux computation routine to use this whole region for the background computation
            - 2: sky_lim = [x_down, x_up, y_down, y_up]
            else: begin
              bck_orad = drs.bck_orad ; IF drs.bck_irad GT drs.bck_orad+1 THEN MESSAGE, 'The outer radius of the background region is two small compared to the inner radius'
              if n_elements(bck_orad) eq 1 then bck_orad = replicate(bck_orad, n_elements(aper_rad))
            end
          endcase
          ; If sky_col is set, make sure that bck_orad is at least bck_irad + aper_rad (but remain within channel)
          if drs.sky_weight then bck_orad = (bck_orad > (bck_irad + aper_rad)) < y_down < y_up
          ; IF drs.bck_irad GT MIN(bck_orad) THEN BEGIN
          ; PRINT, 'Inner background radius too large for file ' + STRING(data_in[i_img].file_id, FORMAT='(I0)') + '. Using ' + STRING(MIN(bck_orad)-1, FORMAT='(I0)') + ' pixels instead.'
          ; bck_irad = MIN(bck_orad)-1
          ; ENDIF
        endif
        ; Crop image around center, which is at 0.5*(npix-1) for even images. Hence, the +0.5.
        ; PRINT, i_img, xcen[i_beam], ycen[i_beam], bck_orad, y_down, y_up
        img_tmp = extrac(reform(img_in[*, *, i_img]), round(xcen[i_beam] + 0.5) - 0.5 * n_clip, round(ycen[i_beam] + 0.5) - 0.5 * n_clip, n_clip, n_clip)
        ; Clean image
        img_tmp = SIGMA_FILTER(temporary(img_tmp), 10, n_sigma = 5, /all_pixels, iterate = 0, monitor = monitor, n_change = nchange)
        ; Mitigate background bias (not needed after all)
        ; img_tmp = CORRECT_BIAS(img_tmp, LIMIT=sky_lim, STAR_POS=0, EXCL_WIDTH=50)
        ; Subtract polynomial to mitigate background bias (not needed either)
        ; img_tmp[FIX(0.5*n_clip-x_down):FIX(0.5*n_clip+x_up),FIX(0.5*n_clip-y_down):FIX(0.5*n_clip+y_up)] -= SFIT_MASK(img_tmp[FIX(0.5*n_clip-x_down):FIX(0.5*n_clip+x_up),FIX(0.5*n_clip-y_down):FIX(0.5*n_clip+y_up)], 35, DEGREE=3)
        ; Rebin the image if N_BIN is set
        if n_bin gt 1 then img_tmp = rebin(img_tmp, n_pix, n_pix) ; ELSE img_tmp = img_tmp
        ; Beam position in extracted array
        x_off = (xcen[i_beam] - (round(xcen[i_beam] + 0.5) - 0.5)) / n_bin ; X offset w/r to closest 4-pixel corner
        y_off = (ycen[i_beam] - (round(ycen[i_beam] + 0.5) - 0.5)) / n_bin ; Y offset w/r to closest 4-pixel corner
        star_pos[i_img, *, i_beam] = [0.5 * (n_pix - 1) + x_off, 0.5 * (n_pix - 1) + y_off] ; 0.5*(n_pix-1) is the middle 4-pixel corner in sub-image
        ; Status is ok, unless changed below
        status = [0]
        ; If photometric frame, compute beam offset w/r to image center (for tip/tilt diagnostic).
        if hdr_in.header.obstype eq 0 and drs.fit_mode gt 0 then begin
          img_tmp = PSF_FIT(img_tmp, star_pos[i_img, 0, i_beam], star_pos[i_img, 1, i_beam], psf_fwhm, fit_method = drs.fit_mode, file = file, median = 0, offset = offset, zoom_ratio = 1, $
            xcen = x0, ycen = y0, fwhm_x = xsig, fwhm_y = ysig, slope = slope, status = status, info = info, plot = plot, no_shift = drs.skip_regis)
          if drs.skip_regis then star_pos[i_img, *, i_beam] = [x0, y0] else star_pos[i_img, *, i_beam] = 0.5 * (n_pix - 1) * [1, 1] ; will recenter the photometric aperture on the star in call of IMG2FLX
          if status[0] eq 0 then begin
            if i_beam eq 0 then begin
              data_in[i_img].xcen_sx = (round(xcen[i_beam] + 0.5) - 0.5) + (x0 - 0.5 * (n_pix - 1)) ; X position in the initial images
              data_in[i_img].ycen_sx = (round(ycen[i_beam] + 0.5) - 0.5) + (y0 - 0.5 * (n_pix - 1)) ; Y position in the initial images
              data_in[i_img].fwhmx_sx = xsig ; Fitted FWHM in X direction
              data_in[i_img].fwhmy_sx = ysig ; Fitted FWHM in Y direction
              if drs.fit_mode eq 5 then data_in[i_img].slope_sx = slope
            endif else begin
              data_in[i_img].xcen_dx = (round(xcen[i_beam] + 0.5) - 0.5) + (x0 - 0.5 * (n_pix - 1)) ; X position in the initial images
              data_in[i_img].ycen_dx = (round(ycen[i_beam] + 0.5) - 0.5) + (y0 - 0.5 * (n_pix - 1)) ; Y position in the initial images
              data_in[i_img].fwhmx_dx = xsig ; Fitted FWHM in X direction
              data_in[i_img].fwhmy_dx = ysig ; Fitted FWHM in Y direction
              if drs.fit_mode eq 5 then data_in[i_img].slope_dx = slope
            endelse
          endif else star_pos[i_img, *, i_beam] = [0, 0]
        endif
        ; Compute flux if status OK
        if status[0] eq 0 then begin
          ; IF hdr_in[i_img].obstype EQ 2 THEN star_pos = [35,75] ELSE star_pos = 0
          good_range = [-16000., 16000.] * mean(hdr_in.header.n_coadd) ; define good range for pixel values (in ADU)
          IMG2FLX, img_tmp, bck_method = drs.bfl_mode, flx_method = drs.flx_mode, aper_rad = aper_rad, bck_irad = bck_irad, bck_orad = bck_orad, bck_cen = bck_cen, good_range = good_range, psf_fwhm = psf_fwhm, phpadu = hdr_in.header.eperadu / cnf.qe, $
            sky_col = drs.sky_col, sky_lim = sky_lim, sky_off = drs.sky_off, sky_weight = drs.sky_weight, star_pos = star_pos[i_img, *, i_beam], nsky = n_sky, str_flx = s0, str_err = s0_err, bck_flx = b0, bck_err = b0_err ; Outputs keywords
          ; Parse results
          str_flx[i_img, i_beam, *] = s0
          str_err[i_img, i_beam, *] = s0_err
          bck_flx[i_img, i_beam, *] = b0
          bck_err[i_img, i_beam, *] = b0_err
          ; Save image (subtract residual background)
          img_sav[*, *, i_img, i_beam] = img_tmp - b0[0]
          ; Estimate bias if needed (only nulling frames). Take symmetric position
          ; Only work if at least 128 columns
          if drs.bias_estim gt 0 and n_xpix gt 128 and hdr_in.header.obstype eq 2 then begin
            ; Derive new center
            x_tmp = abs(xcen[i_beam] - n_xpix)
            ; Crop image around new center
            img_tmp = extrac(reform(img_in[*, *, i_img]), round(x_tmp + 0.5) - 0.5 * n_clip, round(ycen[i_beam] + 0.5) - 0.5 * n_clip, n_clip, n_clip)
            ; Clean image (step 2 here, already applied once in image reduction. Don't use iterate which is quite time consuming)
            img_tmp = SIGMA_FILTER(img_tmp, 10, n_sigma = 5, /all_pixels, iterate = 0, monitor = monitor, n_change = nchange)
            ; Mitigate background bias (doesn't work well)
            ; img_tmp = CORRECT_BIAS(img_tmp, LIMIT=sky_lim, STAR_POS=0, EXCL_WIDTH=50)
            ; sky_lim = [x_up, x_down, bck_orad, 122-bck_orad]
            if n_bin gt 1 then img_tmp = rebin(img_tmp, n_pix, n_pix)
            ; Beamn position in extracted array
            x_off = (x_tmp - (round(x_tmp + 0.5) - 0.5)) / n_bin ; X offset w/r to closest 4-pixel corner
            pos_tmp = [0.5 * (n_pix - 1) + x_off, 0.5 * (n_pix - 1) + y_off] ; 0.5*(n_pix-1) is the middle 4-pixel corner in sub-image
            IMG2FLX, img_tmp, bck_method = drs.bfl_mode, flx_method = drs.flx_mode, aper_rad = aper_rad, bck_irad = bck_irad, bck_orad = bck_orad, bck_cen = bck_cen, good_range = good_range, psf_fwhm = psf_fwhm, phpadu = adu2ph, $
              sky_col = drs.sky_col, sky_lim = sky_lim, sky_off = drs.sky_off, sky_weight = drs.sky_weight, star_pos = pos_tmp, $ ; Input keywords
              str_flx = s0, str_err = s0_err, bck_flx = b0, bck_err = b0_err
            str_flx2[i_img, i_beam, *] = s0
            str_err2[i_img, i_beam, *] = s0_err
            bck_flx2[i_img, i_beam, *] = b0
            bck_err2[i_img, i_beam, *] = b0_err
          endif
          ; Add background flux stored in the header data and computed before background subtraction and before the median of each channel is subtracted
          if n_elements(data_in[i_img].bck_avg) gt 1 then begin
            i_chan = floor(xcen[i_beam] / n_xchan)
            j_chan = floor(ycen[i_beam] / n_ychan)
            bck_flx[i_img, i_beam, *] += data_in[i_img].bck_avg[n_chan * i_chan + j_chan]
            if drs.bias_estim gt 0 and n_xpix gt 128 and hdr_in.header.obstype eq 2 then begin
              i_chan2 = floor(x_tmp / n_xchan)
              j_chan2 = floor(ycen[i_beam] / n_ychan)
              bck_flx2[i_img, i_beam, *] += data_in[i_img].bck_avg[n_chan * i_chan2 + j_chan2]
            endif
          endif else bck_flx[i_img, i_beam, *] += data_in[i_img].bck_avg ; backward compatibility
        endif
        ; Enable file saving below
        save_on = 1
      endif
    endfor
  endfor

  ; Undefine img_in to save memory
  UNDEFINE, img_in

  ; Additional filtering based on fwhm and slope rms
  fwhm_sx = 0.5 * ((data_in.fwhmx_sx) + (data_in.fwhmy_sx))
  fwhm_dx = 0.5 * ((data_in.fwhmx_dx) + (data_in.fwhmy_dx))
  slrms_sx = data_in.ssloprms
  slrms_dx = data_in.dsloprms
  loop_sx = data_in.sloopon
  loop_dx = data_in.dloopon
  tip_sx = data_in.xcen_sx * hdr_in.header.pixscale * 1d3
  if n_elements(tip_sx) gt 1 then tip_sx -= median(tip_sx)
  tlt_sx = data_in.ycen_sx * hdr_in.header.pixscale * 1d3
  if n_elements(tlt_sx) gt 1 then tlt_sx -= median(tlt_sx)
  tt_sx = sqrt(tip_sx ^ 2 + tlt_sx ^ 2)
  tip_dx = data_in.xcen_dx * hdr_in.header.pixscale * 1d3
  if n_elements(tip_dx) gt 1 then tip_dx -= median(tip_dx)
  tlt_dx = data_in.ycen_dx * hdr_in.header.pixscale * 1d3
  if n_elements(tlt_dx) gt 1 then tlt_dx -= median(tlt_dx)
  tt_dx = sqrt(tip_dx ^ 2 + tlt_dx ^ 2)

  ; Plot data before filtering (filtered data are plotted during the null computation)
  if keyword_set(plot) and n_img gt 1 then begin
    plot_path = pth.result_path + 'diagnostic' + pth.sep + 'ao' + pth.sep + drs.date_obs + pth.sep + 'nod' + string(hdr_in.header.nod_id, format = '(I0)') + pth.sep
    if not file_test(plot_path) then file_mkdir, plot_path
    fit = 20. / 1720.
    !p.font = 0
    thick = 5.0
    xthick = 5.0
    ythick = xthick
    cthick = 3.5
    csize = 1.3
    time = reform(transpose(data_in.mjd_obs - min(data_in.mjd_obs)) * 24. * 60. * 60.) ; convert MJD to ellapsed seconds
    ; 1. Plot FWHM vs time
    if max(fwhm_sx) ne 0 or max(fwhm_dx) ne 0 then begin
      plotname = plot_path + drs.date_obs + '_NOD' + string(hdr_in.header.nod_id, format = '(I03)') + '_FWHM.eps'
      PlotXY, /init, inv = inv, /color, /eps, window = [0, 0, 1200, 800] * fit, filename = plotname
      loadct, 0, /silent
      PlotXY, time, fwhm_sx, ylog = ylog, /new, xrange = [0, max(time)], yrange = [0, max(fwhm_sx) > max(fwhm_dx)], title = title, xtitle = 'Elasped time [s]', ytitle = 'FWHM [pix]', grid = 0, $
        xstyle = 1, ystyle = 1, xthick = xthick, ythick = ythick, thick = thick, charthick = cthick, /nodata, /noerase, window = [150, 100, 1150, 700] * fit ; , INSET_UR='a.'
      loadct, 13, /silent
      PlotXY, time, fwhm_sx, /add, color = 85, thick = thick
      PlotXY, time, fwhm_dx, /add, color = 255, thick = thick
      PlotXY, time, loop_sx, /add, color = 85, thick = thick
      PlotXY, time, loop_dx, /add, color = 255, thick = thick
      loadct, 0, /silent
      PlotXY, /fin
    endif
    ; 2. Plot SLOPE rms vs time
    if max(slrms_sx) ne 0 or max(slrms_dx) ne 0 then begin
      plotname = plot_path + drs.date_obs + '_NOD' + string(hdr_in.header.nod_id, format = '(I03)') + '_SLOPE-RMS.eps'
      PlotXY, /init, inv = inv, /color, /eps, window = [0, 0, 1200, 800] * fit, filename = plotname
      loadct, 0, /silent
      PlotXY, time, slrms_dx, ylog = ylog, /new, xrange = [0, max(time)], yrange = [0, max(slrms_sx) > max(slrms_dx)], title = title, xtitle = 'Elasped time [s]', ytitle = 'SLOPE RMS', grid = 0, $
        xstyle = 1, ystyle = 1, xthick = xthick, ythick = ythick, thick = thick, charthick = cthick, /nodata, /noerase, window = [150, 100, 1150, 700] * fit ; , INSET_UR='a.'
      loadct, 13, /silent
      PlotXY, time, slrms_sx, /add, color = 85, thick = thick
      PlotXY, time, slrms_dx, /add, color = 255, thick = thick
      PlotXY, time, loop_sx, /add, color = 85, thick = thick
      PlotXY, time, loop_dx, /add, color = 255, thick = thick
      loadct, 0, /silent
      PlotXY, /fin
    endif
    ; 3. Plot tip/tilt/tt vs time
    for ib = 0, n_beam - 1 do begin
      if ib eq 0 then begin
        tip = tip_sx
        tlt = tlt_sx
      endif else begin
        tip = tip_dx
        tlt = tlt_dx
      endelse
      tt = sqrt(tip ^ 2 + tlt ^ 2)
      plotname = plot_path + drs.date_obs + '_NOD' + string(hdr_in.header.nod_id, format = '(I03)')
      if n_elements(time) le 10 then no_histo = 1
      if max(tip) ne 0 then PLOTALL, time, reform(tip), 0, name = plotname, tag = 'TIP' + string(ib, format = '(I0)'), xtitle = 'Elapsed time [s]', ytitle = 'Tip [mas]', title = ' ', no_histo = no_histo, /bold, /eps
      if max(tlt) ne 0 then PLOTALL, time, reform(tip), 0, name = plotname, tag = 'TLT' + string(ib, format = '(I0)'), xtitle = 'Elapsed time [s]', ytitle = 'Tilt [mas]', title = ' ', no_histo = no_histo, /bold, /eps
      if max(tt) ne 0 then PLOTALL, time, reform(tt), 0, name = plotname, tag = 'TT' + string(ib, format = '(I0)'), xtitle = 'Elapsed time [s]', ytitle = 'Tip/tilt [mas]', title = ' ', no_histo = no_histo, /bold, /eps
    endfor
  endif

  ; Clean results and save
  ; Loop over beam positions
  if not keyword_set(no_save) and keyword_set(save_on) then begin
    for i_b = 0, n_beam - 1 do begin
      ; Only if beam position not 0
      if max(star_pos[*, 0, i_b]) then begin
        ; Store beam ID
        hdr_in.header.beam_id = i_b
        ; Remove bad frames (including open-loop nulling frames)
        idx_keep = where(star_pos[*, 0, i_b] gt 0 and star_pos[*, 1, i_b] gt 0, n_keep)
        ; Remove obvious outliers in aperture flux (5 sigma).
        flx_ob = reform(str_flx[idx_keep, i_b, n_aper - 1]) ; Consider only maximum aperture.
        AVGSDV, flx_ob, avg, rms, kappa = 3
        idx_keep = idx_keep[where(abs(flx_ob - avg) le 5. * rms, /null)]
        ; Remove obvious outliers in background floor (5 sigma).
        bck_ob = reform(bck_flx[idx_keep, i_b, n_aper - 1])
        AVGSDV, bck_ob, avg, rms, kappa = 3
        idx_keep = idx_keep[where(abs(bck_ob - avg) le 5. * rms, /null, n_keep)]
        ; If photometry, remove also images with bad FWHM and too large tip/tilt
        if hdr_in.header.obstype eq 0 then begin
          ; FWHM
          if i_b eq 0 then fwhm = fwhm_sx[idx_keep] else fwhm = fwhm_dx[idx_keep]
          idx_ok = where(finite(fwhm) eq 1, n_ok) ; Make sure there are no Infinite values left
          if n_ok le 0 then begin
            message, 'All fitted FWHM are infinite!!!', /continue
            RETURN
          endif else fwhm = fwhm[idx_ok]
          AVGSDV, fwhm, avg, rms, kappa = 3
          idx_keep = idx_keep[idx_ok[where(abs(fwhm - avg) le 3 * rms, /null)]] ; 1.2 is just a guess, especially because we don't use this criterium for the NULL data
          ; Tip/tilt
          if i_b eq 0 then tt = tt_sx[idx_keep] else tt = tt_dx[idx_keep]
          idx_keep = idx_keep[where(tt lt 40., /null, n_keep)] ; Hard-coded limit of 40 mas (~4 pixels on NOMIC)
          ob_save = ob_in
          if n_keep le 0 then begin
            message, 'No frames survived data filtering. Skipped.', /continue
            RETURN
          endif
        endif else ob_save = ob_all[i_b]
        ; Total number of rejected frames
        hdr_in.header.n_rej = n_rej + (n_img - n_keep)
        ; Save flux in structure
        flx_out = {aper_rad: aper_rad, bck_irad: bck_irad, bck_orad: bck_orad, nsky: n_sky, bckg_meas: reform(bck_flx[idx_keep, i_b, *]), bckg_err: reform(bck_err[idx_keep, i_b, *]), flx_tot: reform(str_flx[idx_keep, i_b, *]), flx_err: reform(str_err[idx_keep, i_b, *]), $
          bckg_meas2: reform(bck_flx2[idx_keep, i_b, *]), bckg_err2: reform(bck_err2[idx_keep, i_b, *]), flx_tot2: reform(str_flx2[idx_keep, i_b, *]), flx_err2: reform(str_err2[idx_keep, i_b, *])}
        ; Save data if requested
        LBTI_SAVEL1IMG, reform(img_sav[*, *, idx_keep, i_b]), hdr_in.header, data_in[idx_keep], flx_out, file_id = ob_save
        LBTI_SAVEL1FLX, flx_out, hdr_in.header, data_in[idx_keep], file_id = ob_save
        ; Save idx for later
        if i_b eq 0 then idx_keep1 = idx_keep
        if i_b eq 1 then idx_keep2 = idx_keep
      endif
    endfor
  endif

  ; Display image for inspection
  ; Destroy window if exists
  nwindow = 29
  if !d.window eq nwindow then wdelete, nwindow

  ; ; Plot the image to screen to make sure the good star is found
  if keyword_set(plot) then begin
    n_show = 2 < n_beam ; show max 2 beams
    xrange = [-0.5, 0.5] * n_pix * hdr_in.header.pixscale
    yrange = [-0.5, 0.5] * n_pix * hdr_in.header.pixscale
    dim = get_screen_size()
    window = [0.02 * dim[0], dim[1] - 0.55 * dim[1], 0.05 * dim[0] + n_show * 0.25 * dim[1], dim[1] - 0.05 * dim[1]]
    sub_img = 0.9 * 0.45 * (window[3] - window[1])
    PlotXY, /init, nwindow = nwindow, window = window
    loadct, 0, /silent
    for i = 0, n_show - 1 do begin
      ; Show the 8 pixels case
      tmp = min(abs(aper_rad - 8), i_aper)
      ; Median frame (weighted by the flux)
      idx_keep = indgen(n_elements(str_flx[*, 0, 0]))
      if i eq 0 and n_elements(idx_keep1) gt 1 then idx_keep = idx_keep1 ; ELSE idx_keep = N_ELEMENTS(str_flx[*,0,0])
      if i eq 1 and n_elements(idx_keep2) gt 1 then idx_keep = idx_keep2 ; ELSE
      n_keep = n_elements(idx_keep)
      ; AVGSDV, n_keep*str_flx[idx_keep,i,i_aper]^2*(star_pos[idx_keep,0,i])/TOTAL(str_flx[idx_keep,i,i_aper]^2), xcen, tmp, KAPPA=3
      ; AVGSDV, n_keep*str_flx[idx_keep,i,i_aper]^2*(star_pos[idx_keep,1,i])/TOTAL(str_flx[idx_keep,i,i_aper]^2), ycen, tmp, KAPPA=3
      AVGSDV, (star_pos[idx_keep, 0, i]), xcen, tmp, kappa = 3
      AVGSDV, (star_pos[idx_keep, 1, i]), ycen, tmp, kappa = 3
      if n_keep gt 1 then img_med = median(img_sav[*, *, idx_keep, i], dimension = 3) else img_med = reform(img_sav[*, *, idx_keep, i])
      sub_win = [0.08 * (window[2] - window[0]) + i * 0.5 * (window[2] - window[0]), 0.54 * (window[3] - window[1]), 0.08 * (window[2] - window[0]) + i * 0.5 * (window[2] - window[0]) + sub_img, 0.54 * (window[3] - window[1]) + sub_img]
      PlotXY, img_med, /new, xrange = xrange, yrange = yrange, title = 'OB median -- beam ' + string(i + 1, format = '(I0)'), xtitle = 'Angular separation [arcsec]', ytitle = 'Angular separation [arcsec]', grid = 0, $
        charsize = 1.0, charthick = 1.5, thick = 2.5, xstyle = 1, ystyle = 1, /noerase, window = sub_win ; , INSET_UR='a.'
      if keyword_set(aper_rad) then begin
        loadct, 13, /silent
        alpha = findgen(100) / 99. * 2. * !dpi
        PlotXY, xrange[0] + [xcen] * hdr_in.header.pixscale, yrange[0] + [ycen] * hdr_in.header.pixscale, color = 90, symbol = 2, /add
        PlotXY, xrange[0] + [xcen + aper_rad[i_aper] * cos(alpha)] * hdr_in.header.pixscale, yrange[0] + [ycen + aper_rad[i_aper] * sin(alpha)] * hdr_in.header.pixscale, thick = 1.5, linestyle = 0, color = 90, /add
        PlotXY, xrange[0] + [xcen + bck_irad[i_aper] * cos(alpha)] * hdr_in.header.pixscale, yrange[0] + [ycen + bck_irad[i_aper] * sin(alpha)] * hdr_in.header.pixscale, thick = 1.5, linestyle = 0, color = 250, /add
        PlotXY, xrange[0] + [xcen + bck_orad[i_aper] * cos(alpha)] * hdr_in.header.pixscale, yrange[0] + [ycen + bck_orad[i_aper] * sin(alpha)] * hdr_in.header.pixscale, thick = 1.5, linestyle = 0, color = 250, /add
        loadct, 0, /silent
      endif
      ; Best null
      tmp = min(str_flx[idx_keep, i, i_aper], i_min)
      i_min = idx_keep[i_min]
      xcen = star_pos[i_min, 0, i]
      ycen = star_pos[i_min, 1, i]
      img_med = reform(img_sav[*, *, i_min, i])
      sub_win = [0.08 * (window[2] - window[0]) + i * 0.5 * (window[2] - window[0]), 0.05 * (window[3] - window[1]), 0.08 * (window[2] - window[0]) + i * 0.5 * (window[2] - window[0]) + sub_img, 0.05 * (window[3] - window[1]) + sub_img]
      PlotXY, img_med, /new, xrange = xrange, yrange = yrange, title = 'Lowest flux -- beam ' + string(i + 1, format = '(I0)'), xtitle = 'Angular separation [arcsec]', ytitle = 'Angular separation [arcsec]', grid = 0, $
        charsize = 1.0, charthick = 1.5, thick = 2.5, xstyle = 1, ystyle = 1, /noerase, window = sub_win ; , INSET_UR='a.'
      if keyword_set(aper_rad) then begin
        loadct, 13, /silent
        alpha = findgen(100) / 99. * 2. * !dpi
        PlotXY, xrange[0] + [xcen] * hdr_in.header.pixscale, yrange[0] + [ycen] * hdr_in.header.pixscale, color = 90, symbol = 2, /add
        PlotXY, xrange[0] + [xcen + aper_rad[i_aper] * cos(alpha)] * hdr_in.header.pixscale, yrange[0] + [ycen + aper_rad[i_aper] * sin(alpha)] * hdr_in.header.pixscale, thick = 1.5, linestyle = 0, color = 90, /add
        PlotXY, xrange[0] + [xcen + bck_irad[i_aper] * cos(alpha)] * hdr_in.header.pixscale, yrange[0] + [ycen + bck_irad[i_aper] * sin(alpha)] * hdr_in.header.pixscale, thick = 1.5, linestyle = 0, color = 250, /add
        PlotXY, xrange[0] + [xcen + bck_orad[i_aper] * cos(alpha)] * hdr_in.header.pixscale, yrange[0] + [ycen + bck_orad[i_aper] * sin(alpha)] * hdr_in.header.pixscale, thick = 1.5, linestyle = 0, color = 250, /add
        loadct, 0, /silent
      endif
    endfor
    PlotXY, /fin
  endif
end
