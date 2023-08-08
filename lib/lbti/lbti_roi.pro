;+
; NAME: LBTI_ROI
; 
; PURPOSE:
;   Returns the region of interest (ROI) for a given subframe
;
; KEYWORDS:
;   INSTRUM    : Specify the instrument (e.g., 'nomic', 'lmircam')
;
; OUTPUTS:
;   roi: 4-element vector with the X,Y coordinates of the ROI (i.e., x_min, y_min, x_max, y_max)
; 
; NOTES:
;   This will have to include the CLKTABLE at some point
; 
; MODIFICATION HISTORY:
;   Version 1.0,  12-APR-2017, by Denis Defrère, Steward Observatory (denis@lbti.org)
;   Version 1.1,  11-AUG-2018, DD: added new ROIs

FUNCTION LBTI_ROI, subsectm, INSTRUM=instrum

; Default ROI
roi = [0,0,0,0]

IF STRLOWCASE(instrum) EQ 'nomic' THEN BEGIN
  CASE subsectm OF 
    'Ch6-9_256x256': roi = [0,0,127,255]
    'Ch4-7_256x256': roi = [0,0,127,255]
    'full_frame': roi = [100,100,900,900] ; Avoid vignetted part
    ELSE: roi = [0,0,0,0]
  ENDCASE
ENDIF

RETURN, roi
END