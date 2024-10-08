pro aper_weight, image, xc, yc, mags, errap, sky, skyerr, phpadu, apr, skyradii, badpix, $
  setskyval = setskyval, print = print, silent = silent, flux = flux, fwhm = fwhm, $
  exact = exact, nan = nan, readnoise = readnoise, meanback = meanback, $
  clipsig = clipsig, maxiter = maxiter, converge_num = converge_num, $
  minsky = minsky, nsky = nsky, skylim = skylim, sky_col = sky_col, $
  sky_weight = sky_weight, weight = weight
  ;+
  ; NAME:
  ;      APER
  ; PURPOSE:
  ;      Compute concentric "weighted" aperture photometry (adapted from DAOPHOT)
  ;      Based on ASTROLIB routine APER.pro
  ; EXPLANATION:
  ;     APER can compute photometry in several user-specified aperture radii.
  ;     A separate sky value is computed for each source using specified inner
  ;     and outer sky radii.
  ;
  ; CALLING SEQUENCE:
  ;     APER, image, xc, yc, [ mags, errap, sky, skyerr, phpadu, apr, skyrad,
  ;                       badpix, /NAN, /EXACT, /FLUX, PRINT = , /SILENT,
  ;                       /MEANBACK, MINSKY=, SETSKYVAL = ]
  ; INPUTS:
  ;     IMAGE -  input image array
  ;     XC     - vector of x coordinates.
  ;     YC     - vector of y coordinates
  ;
  ; OPTIONAL INPUTS:
  ;     PHPADU - Photons per Analog Digital Units, numeric scalar.  Converts
  ;               the data numbers in IMAGE to photon units.  (APER assumes
  ;               Poisson statistics.)
  ;     APR    - Vector of up to 12 REAL photometry aperture radii.
  ;     SKYRAD - Two element vector giving the inner and outer radii
  ;               to be used for the sky annulus.   Ignored if the SETSKYVAL
  ;              keyword is set.
  ;     BADPIX - Two element vector giving the minimum and maximum value
  ;               of a good pixel.   If badpix is not supplied or if BADPIX[0] is
  ;               equal to BADPIX[1] then it is assumed that there are no bad
  ;               pixels.     Note that fluxes will not be computed for any star
  ;               with a bad pixel within the aperture area, but that bad pixels
  ;               will be simply ignored for the sky computation.    The BADPIX
  ;               parameter is ignored if the /NAN keyword is set.
  ;
  ; OPTIONAL KEYWORD INPUTS:
  ;     CLIPSIG - if /MEANBACK is set, then this is the number of sigma at which
  ;             to clip the background.  Default=3
  ;     CONVERGE_NUM:  if /MEANBACK is set then if the proportion of
  ;           rejected pixels is less than this fraction, the iterations stop.
  ;           Default=0.02, i.e., iteration stops if fewer than 2% of pixels
  ;           excluded.
  ;     /EXACT -  By default, APER counts subpixels, but uses a polygon
  ;             approximation for the intersection of a circular aperture with
  ;             a square pixel (and normalizes the total area of the sum of the
  ;             pixels to exactly match the circular area).   If the /EXACT
  ;             keyword, then the intersection of the circular aperture with a
  ;             square pixel is computed exactly.    The /EXACT keyword is much
  ;             slower and is only needed when small (~2 pixels) apertures are
  ;             used with very undersampled data.
  ;     /FLUX - By default, APER uses a magnitude system where a magnitude of
  ;               25 corresponds to 1 flux unit.   If set, then APER will keep
  ;              results in flux units instead of magnitudes.
  ;     FWHM   - Full width at half maximum of the PSF (used for weighted aperture
  ;              photometry). Only when /WEIGHT
  ;     MAXITER if /MEANBACK is set then this is the ceiling on number of
  ;             clipping iterations of the background.  Default=5
  ;     /MEANBACK - if set, then the background is computed using the 3 sigma
  ;             clipped mean (using meanclip.pro) rather than using the mode
  ;             computed with mmm.pro.    This keyword is useful for the Poisson
  ;             count regime or where contamination is known  to be minimal.
  ;      MINSKY - Integer giving mininum number of sky values to be used with MMM
  ;             APER will not compute a flux if fewer valid sky elements are
  ;               within the sky annulus.   Default = 20.
  ;      NSKY - On output, the number of pixels in the sky region
  ;     /NAN  - If set then APER will check for NAN values in the image.   /NAN
  ;             takes precedence over the BADPIX parameter.   Note that fluxes
  ;             will not be computed for any star with a NAN pixel within the
  ;             aperture area, but that NAN pixels will be simply ignored for
  ;             the sky computation.
  ;     PRINT - if set and non-zero then APER will also write its results to
  ;               a file aper.prt.   One can specify the output file name by
  ;               setting PRINT = 'filename'.
  ;     READNOISE - Scalar giving the read noise (or minimum noise for any
  ;              pixel.   This value is passed to the procedure mmm.pro when
  ;              computing the sky, and is only need for images where
  ;              the noise is low, and pixel values are quantized.
  ;     /SILENT -  If supplied and non-zero then no output is displayed to the
  ;               terminal.
  ;     SETSKYVAL - Use this keyword to force the sky to a specified value
  ;               rather than have APER compute a sky value.    SETSKYVAL
  ;               can either be a scalar specifying the sky value to use for
  ;               all sources, or a 3 element vector specifying the sky value,
  ;               the sigma of the sky value, and the number of elements used
  ;               to compute a sky value.   The 3 element form of SETSKYVAL
  ;               is needed for accurate error budgeting.
  ;     SKYLIM -  Four element vector with the outer limit of the background
  ;               region in all 4 directions (computed from the star): [x_down,
  ;               x_up, y_down, y_up]. Superseed (skyrad[1])
  ;     /SKY_COL    - Set this keyword in order to use only pixels in the same column
  ;                   for the sky (+impose same number as in photometric aperture)
  ;     /SKY_WEIGHT - Set this keyword in order to weight the number of pixels per
  ;                   column in the background region according to the corresponding
  ;                   number of pixels in the photometric aperture (time consuming!).
  ;     /WEIGHT - Set this keyword in order to perform weighted aperture photometry
  ;
  ; OUTPUTS:
  ;     MAGS   -  NAPER by NSTAR array giving the magnitude for each star in
  ;               each aperture.  (NAPER is the number of apertures, and NSTAR
  ;               is the number of stars).   If the /FLUX keyword is not set, then
  ;               a flux of 1 digital unit is assigned a zero point magnitude of
  ;               25.
  ;     ERRAP  -  NAPER by NSTAR array giving error for each star.  If a
  ;               magnitude could not be determined then  ERRAP = 9.99 (if in
  ;                magnitudes) or ERRAP = !VALUES.F_NAN (if /FLUX is set).
  ;     SKY  -    NSTAR element vector giving sky value for each star in
  ;               flux units
  ;     SKYERR -  NSTAR element vector giving error in sky values
  ;
  ; EXAMPLE:
  ;       Determine the flux and error for photometry radii of 3 and 5 pixels
  ;       surrounding the position 234.2,344.3 on an image array, im.   Compute
  ;       the partial pixel area exactly.    Assume that the flux units are in
  ;       Poisson counts, so that PHPADU = 1, and the sky value is already known
  ;       to be 1.3, and that the range [-32767,80000] for bad low and bad high
  ;       pixels
  ;
  ;
  ;       IDL> aper, im, 234.2, 344.3, flux, eflux, sky,skyerr, 1, [3,5], -1, $
  ;            [-32767,80000],/exact, /flux, setsky = 1.3
  ;
  ; PROCEDURES USED:
  ;       GETOPT, MMM, PIXWT(), STRN(), STRNUMBER()
  ; NOTES:
  ;       Reasons that a valid magnitude cannot be computed include the following:
  ;      (1) Star position is too close (within 0.5 pixels) to edge of the frame
  ;      (2) Less than 20 valid pixels available for computing sky
  ;      (3) Modal value of sky could not be computed by the procedure MMM
  ;      (4) *Any* pixel within the aperture radius is a "bad" pixel
  ;      (5) The total computed flux is negative.     In this case the negative
  ;          flux and error are returned.
  ;
  ;
  ;       For the case where the source is fainter than the background, APER will
  ;       return negative fluxes if /FLUX is set, but will otherwise give
  ;       invalid data (since negative fluxes can't be converted to magnitudes)
  ;
  ;       APER was modified in June 2000 in two ways: (1) the /EXACT keyword was
  ;       added (2) the approximation of the intersection of a circular aperture
  ;       with square pixels was improved (i.e. when /EXACT is not used)
  ; REVISON HISTORY:
  ;       Adapted to IDL from DAOPHOT June, 1989   B. Pfarr, STX
  ;       FLUX keyword added                       J. E. Hollis, February, 1996
  ;       SETSKYVAL keyword, increase maxsky       W. Landsman, May 1997
  ;       Work for more than 32767 stars           W. Landsman, August 1997
  ;       Don't abort for insufficient sky pixels  W. Landsman  May 2000
  ;       Added /EXACT keyword                     W. Landsman  June 2000
  ;       Allow SETSKYVAL = 0                      W. Landsman  December 2000
  ;       Set BADPIX[0] = BADPIX[1] to ignore bad pixels W. L.  January 2001
  ;       Fix chk_badpixel problem introduced Jan 01 C. Ishida/W.L. February 2001
  ;       Set bad fluxes and error to NAN if /FLUX is set  W. Landsman Oct. 2001
  ;       Remove restrictions on maximum sky radius W. Landsman  July 2003
  ;       Added /NAN keyword  W. Landsman November 2004
  ;       Set badflux=0 if neither /NAN nor badpix is set  M. Perrin December 2004
  ;       Added READNOISE keyword   W. Landsman January 2005
  ;       Added MEANBACK keyword   W. Landsman October 2005
  ;       Correct typo when /EXACT and multiple apertures used.  W.L. Dec 2005
  ;       Remove VMS-specific code W.L. Sep 2006
  ;       Add additional keywords if /MEANBACK is set W.L  Nov 2006
  ;       Allow negative fluxes if /FLUX is set  W.L.  Mar 2008
  ;       Previous update would crash if first star was out of range W.L. Mar 2008
  ;       Fix floating equality test for bad magnitudes W.L./J.van Eyken Jul 2009
  ;       Added MINSKY keyword W.L. Dec 2011
  ;       Don't ever modify input skyrad variable  W. Landsman Aug 2013
  ;       Avoid integer overflow for very big images W. Landsman/R. Gutermuth   Mar 2016
  ;       Eliminate limit on maximum number of sky pixels W. Landsman  Dec 2016
  ;       Updated for weighted aperture photometry D. Defrere Feb 2017
  ;       Added NSKY to output keywords D. Defrere May 2017
  ;-
  compile_opt idl2
  on_error, 2
  ; Set parameter limits
  ; Smallest number of pixels from which the sky may be determined
  if n_elements(minsky) eq 0 then minsky = 20
  ;
  if n_params() lt 3 then begin ; Enough parameters supplied?
    print, $
      'Syntax - APER, image, xc, yc, [ mags, errap, sky, skyerr, phpadu, apr, '
    print, '             skyrad, badpix, /EXACT, /FLUX, SETSKYVAL = ,PRINT=, ]'
    print, '             /SILENT, /NAN, MINSKY='
    return
  endif

  s = size(image)
  if (s[0] ne 2) then message, $
    'ERROR - Image array (first parameter) must be 2 dimensional'
  ncol = s[1]
  nrow = s[2] ; Number of columns and rows in image array

  silent = keyword_set(silent)

  if ~keyword_set(nan) then begin
    if (n_elements(badpix) ne 2) then begin ; Bad pixel values supplied
      get_badpix:
      ans = ''
      print, 'Enter low and high bad pixel values, [RETURN] for defaults'
      read, 'Low and high bad pixel values [none]: ', ans
      if ans eq '' then badpix = [0, 0] else begin
        badpix = getopt(ans, 'F')
        if (n_elements(badpix) ne 2) then begin
          message, 'Expecting 2 scalar values', /continue
          goto, get_badpix
        endif
      endelse
    endif

    chk_badpix = badpix[0] lt badpix[1] ; Ignore bad pixel checks?
  endif

  if (n_elements(apr) lt 1) then begin ; Read in aperture sizes?
    apr = fltarr(10)
    read, 'Enter first aperture radius: ', ap
    apr[0] = ap
    ap = 'aper'
    for i = 1, 9 do begin
      getap:
      read, 'Enter another aperture radius, [RETURN to terminate]: ', ap
      if ap eq '' then goto, done
      result = strnumber(ap, val)
      if result eq 1 then apr[i] = val else goto, getap
    endfor
    done:
    apr = apr[0 : i - 1]
  endif

  if n_elements(setskyval) gt 0 then begin
    if n_elements(setskyval) eq 1 then setskyval = [setskyval, 0., 1.]
    if n_elements(setskyval) ne 3 then message, $
      'ERROR - Keyword SETSKYVAL must contain 1 or 3 elements'
    skyrad = [0., max(apr) + 1]
  endif else begin
    if n_elements(skyradii) ne 2 then begin
      skyrad = fltarr(2)
      read, 'Enter inner and outer sky radius (pixel units): ', skyrad
    endif else skyrad = float(skyradii)
  endelse

  if (n_elements(phpadu) lt 1) then $
    read, 'Enter scale factor in Photons per Analog per Digital Unit: ', phpadu

  Naper = n_elements(apr) ; Number of apertures
  Nstars = min([n_elements(xc), n_elements(yc)]) ; Number of stars to measure

  ms = strarr(Naper) ; String array to display mag for each aperture
  if keyword_set(flux) then $
    fmt = '(F8.1,1x,A,F7.1)' else $ ; Flux format
    fmt = '(F9.3,A,F5.3)' ; Magnitude format
  fmt2 = '(I5,2F8.2,F7.2,1x,3A,3(/,28x,4A,:))' ; Screen format
  fmt3 = '(I4,5F8.2,1x,6A,2(/,44x,9A,:))' ; Print format

  mags = fltarr(Naper, Nstars)
  errap = mags ; Declare arrays
  sky = fltarr(Nstars)
  skyerr = sky
  area = !dpi * apr * apr ; Area of each aperture

  if keyword_set(exact) then begin
    bigrad = apr + 0.5
    smallrad = apr / sqrt(2) - 0.5
  endif

  if n_elements(setskyval) eq 0 then begin
    rinsq = (skyrad[0] > 0.) ^ 2
    routsq = skyrad[1] ^ 2
  endif

  if keyword_set(print) then begin ; Open output file and write header info?
    if size(print, /tname) ne 'STRING' then file = 'aper.prt' $
    else file = print
    message, 'Results will be written to a file ' + file, /inf
    openw, lun, file, /get_lun
    printf, lun, 'Program: APER: ' + systime(), '   User: ', $
      getenv('USER'), '  Host: ', getenv('HOST')
    for j = 0, Naper - 1 do printf, lun, $
      format = '(a,i2,a,f4.1)', 'Radius of aperture ', j, ' = ', apr[j]
    if n_elements(setskyval) eq 0 then begin
      printf, lun, f = '(/a,f4.1)', 'Inner radius for sky annulus = ', skyrad[0]
      printf, lun, f = '(a,f4.1)', 'Outer radius for sky annulus = ', skyrad[1]
    endif else printf, lun, 'Sky values fixed at ', strtrim(setskyval[0], 2)
    if keyword_set(flux) then begin
      printf, lun, f = '(/a)', $
        'Star   X       Y        Sky   SkySig    SkySkw   Fluxes'
    endif else printf, lun, f = '(/a)', $
      'Star   X       Y        Sky   SkySig    SkySkw   Magnitudes'
  endif
  print = keyword_set(print)

  ; Print header
  if ~silent then begin
    if keyword_set(flux) then begin
      print, format = '(/1X,''Star'',5X,''X'',7X,''Y'',6X,''Sky'',8X,''Fluxes'')'
    endif else print, $
      format = '(/1X,''Star'',5X,''X'',7X,''Y'',6X,''Sky'',8X,''Magnitudes'')'
  endif

  ; Compute the limits of the submatrix.   Do all stars in vector notation.
  if keyword_set(skylim) then begin
    if n_elements(skylim) ne 4 then message, 'skylim must have 4 elements'
    lx = long(xc - skylim[0]) > 0 ; Lower limit X direction
    ux = long(xc + skylim[1]) < (ncol - 1) ; Upper limit X direction
    ly = long(yc - skylim[2]) > 0 ; Lower limit Y direction
    uy = long(yc + skylim[3]) < (nrow - 1) ; ;Upper limit Y direction
  endif else begin
    lx = long(xc - skyrad[1]) > 0 ; Lower limit X direction
    ux = long(xc + skyrad[1]) < (ncol - 1) ; Upper limit X direction
    ly = long(yc - skyrad[1]) > 0 ; Lower limit Y direction
    uy = long(yc + skyrad[1]) < (nrow - 1) ; ;Upper limit Y direction
  endelse
  nx = ux - lx + 1 ; Number of pixels X direction
  ny = uy - ly + 1 ; Number of pixels Y direction
  dx = xc - lx ; X coordinate of star's centroid in subarray
  dy = yc - ly ; Y coordinate of star's centroid in subarray

  edge = (dx - 0.5) < (nx + 0.5 - dx) < (dy - 0.5) < (ny + 0.5 - dy) ; Closest edge to array
  badstar = ((xc lt 0.5) or (xc gt ncol - 1.5) $ ; Stars too close to the edge
  or (yc lt 0.5) or (yc gt nrow - 1.5))
  ;
  badindex = where(badstar, Nbad) ; Any stars outside image
  if (Nbad gt 0) then message, /inf, $
    'WARNING - ' + strn(Nbad) + ' star positions outside image'
  if keyword_set(flux) then begin
    badval = !values.f_nan
    baderr = badval
  endif else begin
    badval = 99.999
    baderr = 9.999
  endelse

  for i = 0l, Nstars - 1 do begin ; Compute magnitudes for each star
    apmag = replicate(badval, Naper)
    magerr = replicate(baderr, Naper)
    skymod = 0.
    skysig = 0.
    skyskw = 0. ; Sky mode sigma and skew
    if badstar[i] then goto, badstar
    error1 = apmag
    error2 = apmag
    error3 = apmag

    rotbuf = image[lx[i] : ux[i], ly[i] : uy[i]] ; Extract subarray from image
    ; RSQ will be an array, the same size as ROTBUF containing the square of
    ; the distance of each pixel to the center pixel.

    dxsq = (findgen(nx[i]) - dx[i]) ^ 2
    rsq = fltarr(nx[i], ny[i], /nozero)
    for ii = 0, ny[i] - 1 do rsq[0, ii] = dxsq + (ii - dy[i]) ^ 2

    if keyword_set(exact) or keyword_set(sky_col) or keyword_set(sky_weight) then begin
      nbox = lindgen(nx[i] * ny[i])
      xx = reform((nbox mod nx[i]), nx[i], ny[i])
      yy = reform((nbox / nx[i]), nx[i], ny[i])
      x1 = abs(xx - dx[i])
      y1 = abs(yy - dy[i])
    endif

    r = sqrt(rsq) - 0.5 ; 2-d array of the radius of each pixel in the subarray

    ; Select pixels within sky annulus, and eliminate pixels falling
    ; below BADLO threshold.  SKYBUF will be 1-d array of sky pixels
    if n_elements(setskyval) eq 0 then begin
      if not keyword_set(sky_weight) then begin
        if keyword_set(sky_col) then skypix = (r ge sqrt(rinsq)) and ((x1 - 0.5) ^ 2 le apr[i] ^ 2) and ((y1 - 0.5) ^ 2 le (sqrt(rinsq) - apr[i] + 2 * sqrt((apr[i] ^ 2 - x1 ^ 2) > 0)) ^ 2) else $ ; if sky_col is set, we want pixel between SQRT(rinsq) and SQRT(rinsq) + 2*SQRT(apr[i]^2-x^2)
          skypix = (rsq ge rinsq) and (rsq le routsq)
        if keyword_set(nan) then skypix = skypix and finite(rotbuf) $
        else if chk_badpix then skypix = skypix and (rotbuf gt badpix[0]) and $
          (rotbuf lt badpix[1])
        sindex = where(skypix, nsky)

        if (nsky lt minsky) then begin ; Sufficient sky pixels?
          if ~silent then $
            message, 'There aren''t enough valid pixels in the sky annulus.', /con
          goto, badstar
        endif
        skybuf = rotbuf[sindex[0 : nsky - 1]]

        if keyword_set(meanback) then $
          meanclip, skybuf, skymod, skysig, $
          clipsig = clipsig, maxiter = maxiter, converge_num = converge_num, /double else $
          mmm, skybuf, skymod, skysig, skyskw, readnoise = readnoise, minsky = minsky

        ; Obtain the mode, standard deviation, and skewness of the peak in the
        ; sky histogram, by calling MMM.

        skyvar = skysig ^ 2 ; Variance of the sky brightness
        sigsq = skyvar / nsky ; Square of standard error of mean sky brightness
      endif else begin
        skycol = fltarr(2 * apr[i] + 1)
        sigcol = skycol
        weicol = skycol
        for ia = 0, 2 * apr[i] do begin
          skypix = (rsq ge rinsq) and ((y1 - 0.5) ^ 2 le routsq) and (xx ge dx[i] - apr[i] + ia - 0.5) and (xx lt dx[i] - apr[i] + ia + 0.5)
          if keyword_set(nan) then skypix = skypix and finite(rotbuf) $
          else if chk_badpix then skypix = skypix and (rotbuf gt badpix[0]) and $
            (rotbuf lt badpix[1])
          sindex = where(skypix, nsky)
          skybuf = rotbuf[sindex[0 : nsky - 1]]
          if keyword_set(meanback) then $
            meanclip, skybuf, skymod, skysig, $
            clipsig = clipsig, maxiter = maxiter, converge_num = converge_num, /double else $
            mmm, skybuf, skymod, skysig, skyskw, readnoise = readnoise
          skycol[ia] = skymod
          sigcol[ia] = skysig
          ; Compute weight as number of pixels in the same column and in the photometric aperture
          colpix = (rsq le apr[i] ^ 2) and (xx ge dx[i] - apr[i] + ia - 0.5) and (xx lt dx[i] - apr[i] + ia + 0.5)
          sindex_tmp = where(colpix, Nwei)
          weicol[ia] = Nwei
        endfor
        ; Compute weithed sky value and standard deviation
        skymod = total(weicol * skycol) / total(weicol)
        skyvar = (total(weicol * sigcol) / total(weicol)) ^ 2 ; Variance of the sky brightness
        sigsq = total(weicol * (skycol - skymod) ^ 2) / total(weicol) / n_elements(skycol) ; Square of standard error of mean sky brightness
      endelse

      ; If the modal sky value could not be determined, then all apertures for this
      ; star are bad

      if (skysig lt 0.0) then goto, badstar

      skysig = skysig < 999.99 ; Don't overload output formats
      skyskw = skyskw > (-99) < 999.9
    endif else begin
      skymod = setskyval[0]
      skysig = setskyval[1]
      nsky = setskyval[2]
      skyvar = skysig ^ 2
      sigsq = skyvar / nsky
      skyskw = 0
    endelse

    ; Compute weights
    if keyword_set(weight) then begin
      if not keyword_set(fwhm) then message, 'fwhm must be set when /weight'
      wei = psf_gaussian(npixel = [nx[i], ny[i]], centroid = [dx - 0.5, dy - 0.5], fwhm = fwhm, /double, /normalize)
    endif

    for k = 0, Naper - 1 do begin ; Find pixels within each aperture

      if (edge[i] ge apr[k]) then begin ; Does aperture extend outside the image?
        if keyword_set(exact) then begin
          mask = fltarr(nx[i], ny[i])
          good = where((x1 lt smallrad[k]) and (y1 lt smallrad[k]), Ngood)
          if Ngood gt 0 then mask[good] = 1.0
          bad = where((x1 gt bigrad[k]) or (y1 gt bigrad[k])) ; Fix 05-Dec-05
          mask[bad] = -1

          gfract = where(mask eq 0.0, Nfract)
          if Nfract gt 0 then mask[gfract] = $
            PIXWT(dx[i], dy[i], apr[k], xx[gfract], yy[gfract]) > 0.0
          thisap = where(mask gt 0.0)
          thisapd = rotbuf[thisap]
          fractn = mask[thisap]
        endif else begin
          ;
          thisap = where(r lt apr[k]) ; Select pixels within radius
          thisapd = rotbuf[thisap]
          thisapr = r[thisap]
          fractn = (apr[k] - thisapr < 1.0 > 0.0) ; Fraction of pixels to count
          full = fractn eq 1.0
          gfull = where(full, Nfull)
          gfract = where(1 - full)
          factor = (area[k] - Nfull) / total(fractn[gfract])
          fractn[gfract] = fractn[gfract] * factor
        endelse

        ; Subtract background
        thisapd -= skymod

        ; Weight
        if keyword_set(weight) then thisapd = total(wei[thisap] * thisapd) / total(wei[thisap] ^ 2) * wei[thisap]

        ; If the pixel is bad, set the total counts in this aperture to a large
        ; negative number
        ;
        if keyword_set(nan) then $
          badflux = min(finite(thisapd)) eq 0 $
        else if chk_badpix then begin
          minthisapd = min(thisapd, max = maxthisapd)
          badflux = (minthisapd le badpix[0]) or (maxthisapd ge badpix[1])
        endif else badflux = 0

        if ~badflux then $
          apmag[k] = total(thisapd * fractn) ; Total over irregular aperture
      endif
    endfor ; k
    if keyword_set(flux) then g = where(finite(apmag), Ng) else $
      g = where(abs(apmag - badval) gt 0.01, Ng)
    if Ng gt 0 then begin
      ; apmag[g] = apmag[g] - skymod*area[g]  ;Subtract sky from the integrated brightnesses #Feb 2017, now remove background above (line 472) for better precision

      ; Add in quadrature 3 sources of error: (1) random noise inside the star
      ; aperture, including readout noise and the degree of contamination by other
      ; stars in the neighborhood, as estimated by the scatter in the sky values
      ; (this standard error increases as the square root of the area of the
      ; aperture); (2) the Poisson statistics of the observed star brightness;
      ; (3) the uncertainty of the mean sky brightness (this standard error
      ; increases directly with the area of the aperture).

      error1[g] = area[g] * skyvar ; Scatter in sky values
      error2[g] = (apmag[g] > 0) / phpadu ; Random photon noise
      error3[g] = sigsq * area[g] ^ 2 ; Uncertainty in mean sky brightness
      magerr[g] = sqrt(error1[g] + error2[g] + error3[g])

      if ~keyword_set(flux) then begin
        good = where(apmag gt 0.0, Ngood) ; Are there any valid integrated fluxes?
        if (Ngood gt 0) then begin ; If YES then compute errors
          magerr[good] = 1.0857 * magerr[good] / apmag[good] ; 1.0857 = log(10)/2.5
          apmag[good] = 25. - 2.5 * alog10(apmag[good])
        endif
      endif
    endif

    badstar:

    ; Print out magnitudes for this star

    for ii = 0, Naper - 1 do $ ; Concatenate mags into a string

      ms[ii] = string(apmag[ii], '+-', magerr[ii], form = fmt)
    if print then printf, lun, $ ; Write results to file?
      form = fmt3, i, xc[i], yc[i], skymod, skysig, skyskw, ms
    if ~silent then print, form = fmt2, $ ; Write results to terminal?
      i, xc[i], yc[i], skymod, ms

    sky[i] = skymod
    skyerr[i] = skysig ; Store in output variable
    mags[0, i] = apmag
    errap[0, i] = magerr
  endfor ; i

  if print then free_lun, lun ; Close output file

  return
end

; TEST ROUTINE
; ------------

pro TEST_APERWEIGHT
  compile_opt idl2
  ; Notes:
  ; 1. constant image shows that /EXACT introduce a bias in the results (both APER and APER_WEIGHT). Now solved (Feb 2017!)
  ; 2. mode option introduce a bias for low SNR
  ; 3.
  ; CONCLUSION: use /MEANBACK

  ; Image and PSF paraeaters
  n_pix = 100
  xcen = 50
  ycen = 40
  fwhm = 16
  rad = fwhm / 2
  bck_irad = rad
  bck_orad = bck_irad + rad
  flx_star = 10000 ; typical flux at null
  sig_noise = 10
  real_bck = 6000
  skylim = 0 ; [62,62,62,62]
  clipsig = 5
  sky_weight = 0
  sky_col = 1

  ; Simple test
  k = lonarr(n_pix, n_pix) + real_bck

  ; Look at a constant image with no noise => should find 0!
  APER, k, xcen, ycen, flx_tot0, flx_err0, bck_flx0, bck_err0, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /silent, /meanback ; , /EXACT
  APER, k, xcen, ycen, flx_tot6, flx_err6, bck_flx6, bck_err6, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /silent ; , /EXACT
  aper_weight, k, xcen, ycen, flx_tot1, flx_err1, bck_flx1, bck_err0, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /exact, /silent, /weight, /meanback, fwhm = fwhm, sky_weight = sky_weight, sky_col = sky_col
  aper_weight, k, xcen, ycen, flx_tot2, flx_err2, bck_flx2, bck_err0, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /exact, /silent, /meanback, sky_weight = sky_weight, sky_col = sky_col
  aper_weight, k, xcen, ycen, flx_tot5, flx_err5, bck_flx5, bck_err0, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /silent, /exact, sky_weight = sky_weight, sky_col = sky_col
  aper_weight, k, xcen, ycen, flx_tot3, flx_err3, bck_flx3, bck_err0, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /silent, sky_weight = sky_weight, sky_col = sky_col
  aper_weight, k, xcen, ycen, flx_tot4, flx_err4, bck_flx4, bck_err0, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /silent, /meanback, sky_weight = sky_weight, sky_col = sky_col ; , CLIPSIG=5
  print, 'Constant image'
  print, '  - Stellar flux (APER MEAN, APER MODE, WEIGHT+EXACT+MEAN, EXACT+MEAN, EXACT+MODE, MODE, MEAN)    : ', [flx_tot0, flx_tot6, flx_tot1, flx_tot2, flx_tot5, flx_tot3, flx_tot4]
  print, '  - Background flux (APER MEAN, APER MODE, WEIGHT+EXACT+MEAN, EXACT+MEAN, EXACT+MODE, MODE, MEAN) : ', [bck_flx0, bck_flx6, bck_flx1, bck_flx2, bck_flx5, bck_flx3, bck_flx4]

  ; Now do simulation with NOISE
  t0 = systime(1)
  ts = systime(/seconds)
  n_img = 5000
  flx0 = fltarr(n_img)
  bck0 = flx0
  flx1 = flx0
  bck1 = flx0
  flx2 = flx0
  bck2 = flx0
  flx3 = flx0
  bck3 = flx0
  flx4 = flx0
  bck4 = flx0
  flx5 = flx0
  bck5 = flx0
  for i_img = 0, n_img - 1 do begin
    shft_x = 0 ; i_img
    shft = 0 ; RANDOMN(seed)
    ; k    = img[*,*,i_img];+ sig_noise*RANDOMN(seed, n_pix, n_pix);+flx_star*PSF_GAUSSIAN(NPIXEL=[n_pix,n_pix], FWHM=fwhm, /double, /NORMALIZE)
    seed = ulong(i_img * (systime(/seconds)))
    k = real_bck + sig_noise * randomn(seed, n_pix, n_pix)
    APER, k, xcen + shft_x, ycen, flx_tot0, flx_err0, bck_flx0, bck_err0, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /silent, /meanback, clipsig = clipsig
    APER, k, xcen + shft_x, ycen, flx_tot1, flx_err1, bck_flx1, bck_err1, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /silent ; , /MEANBACK
    aper_weight, k, xcen + shft_x, ycen + shft, flx_tot2, flx_err0, bck_flx2, bck_err0, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /exact, /silent, /meanback, /weight, fwhm = fwhm, clipsig = clipsig, sky_weight = sky_weight, sky_col = sky_col
    aper_weight, k, xcen + shft_x, ycen + shft, flx_tot3, flx_err0, bck_flx3, bck_err0, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /exact, /silent, /meanback, clipsig = clipsig, sky_weight = sky_weight, sky_col = sky_col
    aper_weight, k, xcen + shft_x, ycen + shft, flx_tot4, flx_err0, bck_flx4, bck_err0, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /silent, /meanback, /weight, fwhm = fwhm, clipsig = clipsig, sky_weight = sky_weight, sky_col = sky_col
    aper_weight, k, xcen + shft_x, ycen + shft, flx_tot5, flx_err0, bck_flx5, bck_err0, 154 / 0.4, [rad], [bck_irad, bck_orad], [-10000, 10000], /flux, /silent, /meanback, clipsig = clipsig, sky_weight = sky_weight, sky_col = sky_col
    flx0[i_img] = flx_tot0
    bck0[i_img] = bck_flx0 - real_bck
    flx1[i_img] = flx_tot1
    bck1[i_img] = bck_flx1 - real_bck
    flx2[i_img] = flx_tot2
    bck2[i_img] = bck_flx2 - real_bck
    flx3[i_img] = flx_tot3
    bck3[i_img] = bck_flx3 - real_bck
    flx4[i_img] = flx_tot4
    bck4[i_img] = bck_flx4 - real_bck
    flx5[i_img] = flx_tot5
    bck5[i_img] = bck_flx5 - real_bck
  endfor
  PLOT, indgen(n_img), flx5
  oplot, indgen(n_img), flx0, linestyle = 1
  AVGSDV, flx0, avg0, rms0
  AVGSDV, bck0, avg_b0, rms_b0
  AVGSDV, flx1, avg1, rms1
  AVGSDV, bck1, avg_b1, rms_b1
  AVGSDV, flx2, avg2, rms2
  AVGSDV, bck2, avg_b2, rms_b2
  AVGSDV, flx3, avg3, rms3
  AVGSDV, bck3, avg_b3, rms_b3
  AVGSDV, flx4, avg4, rms4
  AVGSDV, bck4, avg_b4, rms_b4
  AVGSDV, flx5, avg5, rms5
  AVGSDV, bck5, avg_b5, rms_b5
  print, 'Constant+noise image :'
  print, '   - Expected flux  :', flx_star / 2. ; default for fwhm
  print, '   - Retrieved flux :', avg0, avg1, avg2, avg3, avg4, avg5
  print, '   - Retrieved bck  :', avg_b0, avg_b1, avg_b2, avg_b3, avg_b4, avg_b5
  print, '   - Noise on mean  :', rms0 / sqrt(n_img), rms1 / sqrt(n_img), rms2 / sqrt(n_img), rms3 / sqrt(n_img), rms4 / sqrt(n_img), rms5 / sqrt(n_img)
  print, '   - Noise per point:', rms0, rms1, rms2, rms3, rms4, rms5
  print, '   - Background flux:', avg_b0, avg_b1, avg_b2, avg_b3, avg_b4, avg_b5
  print, '   - Expected noise :', sqrt(!dpi * rad ^ 2) * sqrt(2) ; SQRT(2) for background region
  print, '   - SNR            :', avg0 / rms0, avg1 / rms1, avg2 / rms2, avg3 / rms3, avg4 / rms4, avg5 / rms5
  print, ' Execution time     :', systime(1) - t0
end