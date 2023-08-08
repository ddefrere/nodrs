FUNCTION SPTYPE2NUM, type
; Converts spectral type from O1 to M9 to numbers from 1 to 69
CASE STRMID(type, 0, 1) OF
  'O': num = 0
  'B': num = 10
  'A': num = 20
  'F': num = 30
  'G': num = 40
  'K': num = 50
  'M': num = 60
  ELSE: MESSAGE, 'Wrong spectral type'
ENDCASE
subnum = FIX(STRMID(type, 1, 1))
RETURN, num + subnum
END

; ---------------------------

FUNCTION PMAGTABlE, band
; Load the pseudomagnitude table for the appopriate IR band
CASE band OF
  'K': ptable = $
    [[ 5, 0.2432, 0.0443], $
    [ 6, 0.2534, 0.0337], $
    [ 7, 0.2670, 0.0255], $
    [ 8, 0.2830, 0.0192], $
    [ 9, 0.3009, 0.0147], $
    [10, 0.3201, 0.0116], $
    [11, 0.3401, 0.0096], $
    [12, 0.3603, 0.0085], $
    [13, 0.3805, 0.0077], $
    [14, 0.4003, 0.0071], $
    [15, 0.4195, 0.0067], $
    [16, 0.4378, 0.0062], $
    [17, 0.4551, 0.0059], $
    [18, 0.4712, 0.0055], $
    [19, 0.4862, 0.0053], $
    [20, 0.4998, 0.0051], $
    [21, 0.5121, 0.0050], $
    [22, 0.5231, 0.0049], $
    [23, 0.5329, 0.0048], $
    [24, 0.5414, 0.0047], $
    [25, 0.5488, 0.0046], $
    [26, 0.5551, 0.0045], $
    [27, 0.5604, 0.0043], $
    [28, 0.5648, 0.0041], $
    [29, 0.5684, 0.0039], $
    [30, 0.5714, 0.0037], $
    [31, 0.5738, 0.0034], $
    [32, 0.5757, 0.0032], $
    [33, 0.5773, 0.0030], $
    [34, 0.5786, 0.0028], $
    [35, 0.5798, 0.0026], $
    [36, 0.5810, 0.0025], $
    [37, 0.5822, 0.0024], $
    [38, 0.5835, 0.0023], $
    [39, 0.5850, 0.0022], $
    [40, 0.5867, 0.0021], $
    [41, 0.5888, 0.0020], $
    [42, 0.5912, 0.0020], $
    [43, 0.5939, 0.0019], $
    [44, 0.5971, 0.0019], $
    [45, 0.6006, 0.0019], $
    [46, 0.6045, 0.0019], $
    [47, 0.6088, 0.0019], $
    [48, 0.6134, 0.0020], $
    [49, 0.6183, 0.0020], $
    [50, 0.6235, 0.0021], $
    [51, 0.6289, 0.0021], $
    [52, 0.6345, 0.0022], $
    [53, 0.6402, 0.0023], $
    [54, 0.6460, 0.0025], $
    [55, 0.6518, 0.0026], $
    [56, 0.6575, 0.0028], $
    [57, 0.6632, 0.0029], $
    [58, 0.6686, 0.0030], $
    [59, 0.6739, 0.0031], $
    [60, 0.6790, 0.0029], $
    [61, 0.6838, 0.0026], $
    [62, 0.6884, 0.0023], $
    [63, 0.6928, 0.0025], $
    [64, 0.6971, 0.0039], $
    [65, 0.7012, 0.0065], $
    [66, 0.7054, 0.0104]]
  'H': MESSAGE, 'No pseudomag table available for H band yet'
  'J': MESSAGE, 'No pseudomag table available for J band yet'
ENDCASE
RETURN, ptable
END

; ---------------------------

FUNCTION SBR_CHELLI, mv, e_mv, mir, e_mir, type, BAND=band

; DESCRIPTION
;   Uses experimental surface-brightness relations to estimate stellar LD angular diameters.
;   Based on Chelli et al. (2016) -- applicable for all luminosity classes
;
; INPUTS
;   mv   : the V band magnitude
;   e_mv : the error on the V band magnitude
;   mir  : the IR magnitude (can be J, H or K band)
;   e_mir: the error on the IR magnitude (can be J, H or K band)
;   type : a string or array of two strings defining the spectral type or range of spectral types
;
; KEYWORDS
;   band : specifiy the IR band in which the IR magnitude is given. Default = K band
;
; OUTPUT
;   A two-element array containing the LD diameter and its error bar [mas]
;
; REFERENCES
;   Chelli et al. 2016, A&A 589, A112: "Pseudomagnitudes and differential surface brightness: 
;     Application to the apparent diameter of stars"
;
; MODIFICATION HISTORY:
;   Version 1.0, 22-NOV-2016, by Olivier Absil (ULg), absil@astro.ulg.ac.be
  
IF NOT KEYWORD_SET(band) THEN band = 'K'

IF N_ELEMENTS(type) EQ 1 THEN BEGIN
  type0 = type 
ENDIF ELSE IF N_ELEMENTS(type) EQ 2 THEN BEGIN
  type0 = type[0]
  type1 = type[1]
ENDIF ELSE MESSAGE, 'The keyword TYPE can contain only one or two elements'

; Define the extinction coefficient ratios
cv = 1.0
CASE band OF
  'J': cir = 0.28
  'H': cir = 0.17
  'K': cir = 0.12
  ELSE: MESSAGE, 'Input band is not a valid band'
ENDCASE

; Convert spectral type to numbers and define the mean spectral type if a range of sepctral types is defined
num0 = SPTYPE2NUM(type0)
IF num0 LT 5 THEN MESSAGE, 'Spectral type cannot be earlier than O5'
IF num0 GT 66 THEN MESSAGE, 'Spectral type cannot be later than M6'

IF KEYWORD_SET(type1) THEN BEGIN
  num1 = SPTYPE2NUM(type1)
  IF num1 LT 5 THEN MESSAGE, 'Spectral type cannot be earlier than O5'
  IF num1 GT 66 THEN MESSAGE, 'Spectral type cannot be later than M6'
  num_med = (num0+num1)/2D0
  num_max = MAX([num0,num1])
  num_min = MIN([num0,num1])
  num_range = (num_max - num_min) 
ENDIF ELSE BEGIN
  num_med = num0
ENDELSE

; Read the pseudomag table for the appopriate band
ptable = PMAGTABLE(band)

; Extract the appropriate line from the pseudomag table
pseudomag = INTERPOL(ptable[1,*], ptable[0,*], num_med)
IF KEYWORD_SET(num_range) THEN sig_pseudomag = SQRT(VARIANCE(ptable[1,num_min-ptable[0,0]:num_max-ptable[0,0]]) + ptable[2,num_med-ptable[0,0]]^2)$
  ELSE sig_pseudomag = ptable[2,num_med-ptable[0,0]]
 
diam = 10D0^(-0.2*(cv*mir - cir*mv)/(cv-cir)+pseudomag)

diam_err = diam * ALOG(10D0)* SQRT(0.04*(cv^2*e_mir^2 + cir^2*e_mv^2)/(cir-cv)^2 + sig_pseudomag^2) 

RETURN, [diam, diam_err]

END
