FUNCTION chi2_pdf, x, df=df

  ; ----------------------------------------------------------
  ;+
  ; NAME:
  ;       CHI2_PDF
  ;
  ; PURPOSE:
  ;       Compute density for chi-sq distribution
  ;
  ; AUTHOR:
  ;       Simon Vaughan
  ;
  ; CALLING SEQUENCE:
  ;       d = CHI2_PDF(x, df=2)
  ;
  ; INPUTS:
  ;       x       - (float vector) values (x>0) at which to compute density
  ;       df      - (integer) degrees of freedom
  ;
  ;
  ; OUTPUTS:
  ;       d       - (float vector) densities pdf(x)
  ;
  ; OPTIONAL OUTPUTS:
  ;       none
  ;
  ; DETAILS:
  ;       Compute the probability density function (PDF) for the
  ;       chi-square distribution.
  ;
  ;       Example usage
  ;       IDL> x = INDGEN(100, /FLOAT) * 0.1
  ;       IDL> d = CHI2_PDF(x, df=4)
  ;       IDL> PLOT, x, d
  ;
  ; HISTORY:
  ;       18/07/13  - v1.0 - first working version
  ;
  ; USES:
  ;       none
  ;
  ; NOTES:
  ;       none
  ;
  ;-
  ; ----------------------------------------------------------

  ; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

  ; watch out for errors

  ON_ERROR, 3

  ; ----------------------------------------------------------
  ; check arguments

  n = N_ELEMENTS(x)

  ; is the sampling period dt supplied? (default = 1.0)

  if (n eq 0) then begin
    PRINT, '-- Missing X values in CHI2_PDF'
    RETURN, !NULL
  ENDIF

  if NOT KEYWORD_SET(df) THEN df=5
  df = ROUND(df)

  ; ----------------------------------------------------------
  ; main routine

  d = x^(0.5*df-1) * EXP(-0.5*x) / 2.0^(0.5*df) / GAMMA(0.5*df)

  ; ----------------------------------------------------------
  ; return values and end

  RETURN, d

END
