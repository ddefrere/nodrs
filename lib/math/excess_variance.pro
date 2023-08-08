FUNCTION excess_variance, X, DX, MEANX=meanx, MSERR=mserr, $
  VARX=varx, DIMENSION=dimension, $
  CHATTER=chatter, ERROR=error

  ; ----------------------------------------------------------
  ;+
  ; NAME:
  ;       EXCESS_VARIANCE
  ;
  ; PURPOSE:
  ;       Compute the 'excess' variance of a sample
  ;
  ; AUTHOR:
  ;       Simon Vaughan (U.Leicester)
  ;
  ; CALLING SEQUENCE:
  ;       xs = EXCESS_VARIANCE(x, dx)
  ;
  ; INPUTS:
  ;       x        - (array) sample of data
  ;       dx       - (array) std devs of data x
  ;
  ; OPTIONAL INPUTS:
  ;       dimension - (integer) dimension of array use (0, 1 or 2)
  ;
  ; OUTPUTS:
  ;       xs       - (vector) the excess variance
  ;
  ; OPTIONAL OUTPUTS:
  ;       meanx    - (vector) mean value of x
  ;       varx     - (vector) total variance of x
  ;       mserr    - (vector) mean square 'error' on x
  ;       error    - (vector) Error on xs variance
  ;
  ; DETAILS:
  ;       Calculates the 'excess' variance, which is the difference
  ;       between the total variance of a sample of estimates and that
  ;       expected based on the 'errors' in the estimates.
  ;
  ;       Consider estimates of x[i], i=0,1,2,...,N-1. Our data are
  ;       estimates x_obs[i] = x[i] + err[i], where err[i] are the
  ;       'errors' - the standard deviations of the estimates about the
  ;       true values. So err[i]^2 are the variances of the
  ;       estimates. If the population variance of x is S_x then the
  ;       total variance of the data will be (in expectation):
  ;
  ;              S_tot = S_x + E[err^2]
  ;
  ;       where E[err^2] is the expectation of the squared error. So
  ;       the intrinsic sample variance is
  ;
  ;              S_x = S_tot - E[err^2]
  ;
  ;       And we can estimate the two terms on the right sides using the
  ;       sample variance of x_obs and the mean square error. The
  ;       difference between these two is the 'excess' variance. In the
  ;       limit of S_x = 0 the total variance is dominated by the error
  ;       terms, and so S_tot ~ MEAN(err^2). In this case the excess
  ;       variance will be randomly distributed about zero (so may be
  ;       negative).
  ;
  ;       We compute the mean of X as MEANX, the total variance as VARX,
  ;       the mean square error as MSERR and the difference as XS.
  ;       All of these are available as output if needed.
  ;
  ;       If the input (X and DX) are one-dimensional arrays (vectors) then
  ;       we compute a single number of excess variance, etc. If the input
  ;       as an array of dimensions [N,M] then we have a choice: compute the
  ;       excess variance from all elements of the array, or compute N
  ;       variances (one for each column vector of length M),or compute M
  ;       variances (one for each row vector of length N). The
  ;       DIMENSION keyword allows for this choice. A value of 0 (default)
  ;       will compute using all elements, producing a single valued output.
  ;       DIMENSION=1 will give an M element vector, with the variances
  ;       computed within each column; DIMENSION=2 will give an N element
  ;       vector with the variance computed within each row. (NB: in IDL
  ;       an [N,M] array is said to have N columns, M rows.)
  ;
  ;       We return an error on the excess variance using the formula of
  ;       Vaughan et al. (2003; MNRAS, 345, 1271), see also Appendix B of
  ;       A. Wilkinson's PhD thesis (2011). This accounts for the effect of
  ;       the "measurement errors" DX on the estimation of excess variance,
  ;       it does not account for the intrinsic randomness in the true
  ;       variance itself (if the underlying process is random).
  ;
  ; EXAMPLE USAGE:
  ;
  ;       N = 1000
  ;       x = RANDOMN(seed, N) * 5 + 100
  ;       sigma = SQRT(ABS(x))
  ;       err = sigma * RANDOMN(seed, N)
  ;       x_obs = x + err
  ;       xs = EXCESS_VARIANCE(x_obs, sigma, mserr=mserr)
  ;       PRINT, xs, mserr
  ;
  ; HISTORY:
  ;       03/06/2012 - v1.1 - added DIMENSION keyword
  ;       19/07/2012 - v1.2 - removed DOUBLE keyword, now always
  ;                            work in double prec.
  ;                            Added more explanation of DIMENSION
  ;                            keyword.
  ;       28/06/2013 - v1.3 - added CHATTER input, NaN checking,
  ;                            added ERROR output
  ;
  ; NOTES:
  ;-
  ; ----------------------------------------------------------

  ; options for compilation (recommended by RSI)

  COMPILE_OPT idl2, HIDDEN

  ; watch out for errors

  ON_ERROR, 0

  ; ---------------------------------------------------------
  ; Check arguments

  N = N_ELEMENTS(x)

  IF (N le 3) THEN BEGIN
    PRINT,'** Not enough data for EXCESS_VARIANCE'
    RETURN, !NULL
  ENDIF

  if NOT KEYWORD_SET(dimension) then dimension=0

  ; define the level out output to screen

  IF NOT KEYWORD_SET(chatter) THEN chatter = 0

  ; ---------------------------------------------------------
  ; Main routine

  ; compute the mean square error MSERR

  N_err = N_ELEMENTS(dx)

  IF (N_err eq 1) THEN BEGIN
    mserr = dx*dx
  ENDIF ELSE BEGIN
    IF (N_err lt N) THEN BEGIN
      PRINT,'** X and DX lengths do not match in EXCESS_VARIANCE'
      RETURN, -1
    ENDIF
    mserr = MEAN(dx*dx, DOUBLE=KEYWORD_SET(double), DIMENSION=dimension, /NAN)
  ENDELSE

  mask = WHERE(dx le 0, count)

  IF (chatter gt 0) THEN BEGIN
    IF (count gt 0) THEN BEGIN
      PRINT,'** Negative or zero values of DX in EXCESS_VARIANCE'
    ENDIF
  ENDIF

  ; compute mean and total variance

  meanx = MEAN(x, /DOUBLE, DIMENSION=dimension, /NAN)
  varx = VARIANCE(x, /DOUBLE, DIMENSION=dimension, /NAN)

  ; compute excess variance

  xs = varx - mserr

  ; [Added 28/06/2013]
  ; estimate error on xs using eqn 11 of Vaughan et al. (2003)
  ; (see also appendix B of A. Wilkinson's PhD thesis, 2011)

  if (dimension eq 0) then N = N_ELEMENTS(x)
  if (dimension gt 1) then N = (SIZE(x))[dimension-1]
  err2 = 2*mserr^2/N + 4*mserr^2*xs^2/N
  error = SQRT(err2)

  ; ---------------------------------------------------------
  ; Return to user

  RETURN, xs

END