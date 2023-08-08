;+
; NAME: READ_PICKLES
; 
; DESCRIPTION
;   Opens an ASCII file and reads values in the lines containing the spectrum tabulated by Pickles.
;   Lines starting with the character # are discarded.
;
; INPUT
;   spectrum_file: the name of the file containing the Pickles spectrum.
;
; KEYWORDS
;   lam_min  : if set, truncates the output at the specified minimum wavelength [m]
;   lam_max  : if set, truncates the output at the specified maximum wavelength [m]
;   nlambda  : if set, specifies the number of wavelength channels on which the spectrum will be sampled
;   standard : 0 (default) - result in energy per unit wavenumber (proportional to Jy)
;              1           -           energy per unit wavelength (proportional to W/m2/m)
;              2           -           photons per unit wavelength (proportional to ph/s/m)
;              3           -           photons per unit wavenumber (proportional to ph/s/Hz)
;
; OUTPUT
;   spec     : two-columned table containing the wavelengths [m] and the spectrum
;
; MODIFICATION HISTORY
;   Version 1.0, 18-MAR-2008, by Olivier Absil, LAOG, olivier.absil@obs.ujf-grenoble.fr
;   Version 1.1, 02-DEC-2015, DD: cleanup comments

FUNCTION READ_PICKLES, spectrum_file, LAM_MIN=lam_min, LAM_MAX=lam_max, NLAMBDA=nlambda, STANDARD=standard

h=6.6262D-34  ; J s
c=2.9979D+8   ; m/s

; Read the file containing the spectrum of the disk (optical thickness for DISKPIC, number of zodis for ZODIPIC)
OPENR, lun, spectrum_file, /GET_LUN, ERROR=err

IF err THEN BEGIN
  MESSAGE,'Error opening defaults file; starting file finding dialog',/INFO
  deffile=DIALOG_PICKFILE(/READ, TITLE='Locate text file containing the disk spectrum')
  OPENR, lun, deffile, /GET_LUN, ERROR=err3
  IF err3 THEN MESSAGE,'Selected defaults file corrupt'
ENDIF

text = ''
WHILE NOT EOF(lun) DO BEGIN
  READF, lun, text
  IF STRMID(text,0,1) NE '#' THEN BEGIN
    pos = STRPOS(text, ' ')
    pos2 = STRPOS(text, ' ', pos+1)
    IF pos EQ 0 THEN BEGIN
      WHILE pos2 - pos EQ 1 DO BEGIN
        pos = pos2
        pos2 = STRPOS(text, ' ', pos+1)
      ENDWHILE
    ENDIF
    lam = DOUBLE(STRMID(text, 0, pos2))
    pos3 = STRPOS(text, ' ', pos2+1)
    WHILE pos3 - pos2 EQ 1 DO BEGIN
      pos2 = pos3
      pos3 = STRPOS(text, ' ', pos2+1)
    ENDWHILE
    flux = DOUBLE(STRMID(text, pos2, pos3-pos2))
    IF NOT KEYWORD_SET(spec) THEN spec = [lam*1D-10,flux] ELSE spec = [[TEMPORARY(spec)], [lam*1D-10,flux]]
  ENDIF
ENDWHILE
CLOSE,lun & FREE_LUN, lun

IF KEYWORD_SET(lam_min) THEN imin = MIN(WHERE(spec[0,*] GE lam_min)) ELSE imin=0
IF KEYWORD_SET(lam_max) THEN imax = MAX(WHERE(spec[0,*] LE lam_max)) ELSE imax=N_ELEMENTS(spec[0,*])-1
IF (imin EQ -1) OR (imax EQ -1) THEN MESSAGE, 'The input file does not contain data between the two specified wavelengths.'
IF KEYWORD_SET(nlambda) THEN BEGIN
  lambda = DINDGEN(nlambda)/(nlambda-1D0)*(spec[0,imax]-spec[0,imin]) + spec[0,imin]
  spectrum = INTERPOL(REFORM(spec[1,imin:imax]),REFORM(spec[0,imin:imax]),lambda) ; flux per unit wavelength
ENDIF ELSE BEGIN
  lambda = REFORM(spec[0,imin:imax])
  spectrum = REFORM(spec[1,imin:imax]) ; flux per unit wavelength
ENDELSE

IF NOT KEYWORD_SET(standard) THEN standard = 0
CASE standard OF
  0   : spectrum = spectrum  * lambda^2/c ; energy per unit wavenumber
  ;1   : spectrum = spectrum * c / lambda^2
  2   : spectrum = spectrum / (h*c/lambda) ; photons per unit wavelength
  3   : spectrum = spectrum * lambda^3/c^2/h ; photons per unit wavenumber
  ELSE: IF standard NE 1 THEN MESSAGE, 'The STANDARD keyword must be an integer ranging between 0 and 3.'
ENDCASE

; Normalise the output spectrum to avoid small or large numbers
mean_spec = MEAN(spectrum)
IF mean_spec GT 0 THEN spectrum = TEMPORARY(spectrum)/mean_spec

RETURN, [TRANSPOSE(lambda),TRANSPOSE(spectrum)]

END