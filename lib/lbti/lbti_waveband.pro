;+
; NAME: LBTI_WAVEBAND
;
; PURPOSE:
;   Main procedure to compute the wavelength and bandwidth based on the input filters
;
; KEYWORDS:
;   INSTRUME   :  Set to the instrument name
;   LMIR_FW1   :  LMIRCAM filter wheel 1 (superseed keyword in the input files)
;   LMIR_FW2   :  LMIRCAM filter wheel 2 (superseed keyword in the input files)
;   LMIR_FW3   :  LMIRCAM filter wheel 3 (superseed keyword in the input files)
;   LMIR_FW4   :  LMIRCAM filter wheel 4 (superseed keyword in the input files)
;   NOM_FW1    :  NOMIC filter wheel 1 (superseed keyword in the input files)
;   NOM_FW2    :  NOMIC filter wheel 2 (superseed keyword in the input files)
;   
; OUTPUTS
;   lam_cen    :  central wavelength 
;   bandwidth  :  bandwidth 
;
; KEYWORDS:
;
; MODIFICATION HISTORY:
;   Version 1.0, 14-JAN-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 22-DEC-2015, now continue if unknown instrument name.

PRO LBTI_WAVEBAND, lam_cen, bandwidth, INSTRUM=instrum, LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2

; Compute wavelength and bandwidth
CASE instrum OF 
  'NOMIC'  : BEGIN
                ; Rough LBTI/NOMIC bandpass (without filter)
                lam_min = 7D-6
                lam_max = 2D-5
                ; Scan through the filter wheels
                CASE STRTRIM(nom_fw1,2) OF
                  'W08699-9'  : BEGIN
                    lam_min1  = 8.13D-6
                    lam_max1  = 9.35D-6
                  END
                  'W08699-9_122' : BEGIN
                    lam_min1  = 8.13D-6
                    lam_max1  = 9.35D-6
                  END
                  'Nprime'  : BEGIN
                    lam_min1  = 9.81D-6
                    lam_max1  = 12.41D-6
                  END
                  'N/A': BEGIN
                    lam_min1  = 0D
                    lam_max1  = 1D
                  END
                  ELSE : BEGIN
                    lam_min1  = 0D
                    lam_max1  = 1D
                  END
                ENDCASE
                lam_min2  = lam_min
                lam_max2  = lam_max
                lam_min3  = lam_min
                lam_max3  = lam_max
                lam_min4  = lam_min
                lam_max4  = lam_max
             END
  'LMIRCAM': BEGIN
                ; Rough LBTI/LMIRCAM bandpass (without filter)
                lam_min = 1D-6
                lam_max = 6D-6
                ; Scan through the filter wheels
                lam_min1  = lam_min
                lam_max1  = lam_max
                CASE STRTRIM(lmir_fw2,2) OF
                  'L-cont2'  : BEGIN
                    lam_min2  = 3.39D-6
                    lam_max2  = 3.54D-6
                  END
                  'L-cont4'  : BEGIN
                    lam_min2  = 3.68D-6
                    lam_max2  = 3.88D-6
                  END
                  'SX-Half-moon': BEGIN
                    lam_min2  = 1D-6
                    lam_max2  = 6D-6
                  END
                  'DX-Half-moon': BEGIN
                    lam_min2  = 1D-6
                    lam_max2  = 6D-6
                  END
                  'Open': BEGIN
                    lam_min2  = 1D-6
                    lam_max2  = 6D-6
                  END
                  'ND2.0-T1': BEGIN
                    lam_min2  = 1D-6
                    lam_max2  = 6D-6
                  END
                  'N/A': BEGIN
                    lam_min2  = 0D
                    lam_max2  = 1D
                  END
                  'ND1.0-T10': BEGIN
                    lam_min2  = 0D
                    lam_max2  = 1D
                  END
                  ELSE : BEGIN
                    lam_min2  = 0D
                    lam_max2  = 1D
                  END
                ENDCASE
                CASE STRTRIM(lmir_fw3,2) OF
                  'Lp'  : BEGIN
                    lam_min3  = 2.81D-6
                    lam_max3  = 4.00D-6
                  END
                  'Std-M': BEGIN
                    lam_min3  = 4.60D-6
                    lam_max3  = 4.97D-6
                  END
                  'L-short': BEGIN
                    lam_min3  = 3.10D-6
                    lam_max3  = 3.51D-6
                  END
                  'Blank': BEGIN
                    lam_min3  = 0D
                    lam_max3  = 1D
                  END
                  'Open': BEGIN
                    lam_min3  = 1D-6
                    lam_max3  = 6D-6
                  END
                  ELSE : BEGIN
                    lam_min3  = 0D
                    lam_max3  = 1D
                  END
                ENDCASE
                CASE STRTRIM(lmir_fw4,2) OF
                  'H2O-Ice2'  : BEGIN
                    lam_min4  = 3.01D-6
                    lam_max4  = 3.15D-6
                  END
                  'Std-L'  : BEGIN
                    lam_min4  = 3.41D-6
                    lam_max4  = 3.99D-6
                  END
                  'Br-Alpha-On': BEGIN
                    lam_min4  = 4.02D-6
                    lam_max4  = 4.09D-6
                  END
                  'Br-Alpha-Off': BEGIN
                    lam_min4  = 3.972D-6
                    lam_max4  = 4.04D-6
                  END
                  'Open': BEGIN
                    lam_min4  = 1D-6
                    lam_max4  = 6D-6
                  END
                  ELSE: BEGIN
                    lam_min4  = 0D
                    lam_max4  = 1D
                  END
                ENDCASE
             END
  ELSE : BEGIN
      MESSAGE, 'Undefined instrument name', /CONTINUE
      lam_min1 = 0 & lam_min2 = 0 & lam_min3 = 0 & lam_min4 = 0
      lam_max1 = 0 & lam_max2 = 0 & lam_max3 = 0 & lam_max4 = 0
     END
ENDCASE

; Define central wavelength and bandwidth
lam_min   = lam_min1 > lam_min2 > lam_min3 > lam_min4 
lam_max   = lam_max1 < lam_max2 < lam_max3 < lam_max4
bandwidth = lam_max-lam_min
lam_cen   = lam_min + 0.5*bandwidth

; Consistency checj
IF lam_cen GT 5D-5 THEN MESSAGE, 'There is a problem with the wavelength computation'
END