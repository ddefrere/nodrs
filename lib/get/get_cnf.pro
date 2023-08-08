;+
; NAME: GET_CNF
; 
; PURPOSE:
;   Initializes instrumental and detector parameters
;
; KEYWORDS:
;   BANDWIDTH  : User-defined bandwidth (in m, sueprseed the hard-coded values) 
;   LAMBDA     : User-defined wavelength (in m, sueprseed the hard-coded values) 
;   INSTRUM    : Specify the instrument (e.g., 'nomic', 'lmircam')
;
; OUTPUTS:
;   cnf: structure containing all relevant parameters.
;
; MODIFICATION HISTORY:
;   Version 1.0,  22-MAR-2013, by Denis Defrère, Steward Observatory (ddefrere@email.arizona.edu)
;   Version 1.1,  28-JUN-2013, DD, added keywords 'LAMBDA' and 'BANDWIDTH'
;   Version 1.2,  03-APR-2015, DD, added quantum efficiency and read-noise
;   Version 1.3,  10-APR-2015, DD, added adu2e
;   Version 1.4,  15-NOV-2015, DD, replaced NOMIC and LMIRCAM keywords by INSTRUM
;   Version 1.5,  12-FEB-2016, DD, added X_CHAN and Y_CHAN to output structure
;   Version 1.6,  26-MAR-2016, DD, added MIN_ADU and MAX_ADU to output structure
;   Version 1.7,  19-OCT-2016, DD, added adu2e for the new preamp
;   Version 1.8,  19-JUL-2018, DD, reversed the order of det2adu since the preamp is now used by default
;      
; TO-DO LIST
;   Modify to read instead an instrumental configuration file

PRO GET_CNF, cnf, LAMBDA=lambda, BANDWIDTH=bandwidth, INSTRUM=instrum

; Keyword check
IF NOT KEYWORD_SET(INSTRUM) THEN MESSAGE, 'You must specifiy one detector!"

; Define everything here (temporary implementation)
tel_diam  = 8.25           ; m  (we doesn't see the full primary which is 8.4-m wide)
cen_obs   = 0.1*tel_diam   ; central obscuration 
base      = 14.4           ; m
h         = 6.6262D-34     ; [J s]

; Wavelength definition
IF instrum NE 'lmircam' THEN BEGIN
  IF NOT KEYWORD_SET(lambda)    THEN lambda    = 1.11D-5         ; m
  IF NOT KEYWORD_SET(bandwidth) THEN bandwidth = 2.6D-6         ; m
  n_chan = N_ELEMENTS(lambda)   
  jy2ph  = 1D-26/h*(bandwidth/lambda)  ; [ph/s/m2] converts Jy to photons per second per m2    
  
  ; Detector parameters (Raytheon)
  nomic    = 1
  lmircam  = 0
  pix_size = 17.85          ; mas
  qe       = 0.40           ; quantum efficiency
  ron      = 450            ; read noise in e-
  well_det = 1.1D7          ; e-       
  min_adu  = [3290.3500]    ; Approximate minimum number of ADU in linear regime (assuming 1 co-add). For level shift and preamp respectively.
  max_adu  = [10500,14600]   ; Approximate maximum number of ADU in linear regime (assuming 1 co-add). For level shift and preamp respectively. (10500 instead of 8800 for old data!)
  adu2e    = [67.,670.,145.,1500.]   ; ADU to electron conversion factor with the the preAmp (high and low gain) and the level shift (high and low gains) 
  adu2ph   = FIX(adu2e/qe)
  x_chan   = 512            ; width of each channel in the X direction
  y_chan   = 128            ; width of each channel in the Y direction
ENDIF ELSE BEGIN
  IF NOT KEYWORD_SET(lambda)    THEN lambda    = 3.5D-6         ; m
  IF NOT KEYWORD_SET(bandwidth) THEN bandwidth = 2.5D-6         ; m
  n_chan = N_ELEMENTS(lambda)
  jy2ph  = 1D-26/h*(bandwidth/lambda)  ; [ph/s/m2] converts Jy to photons per second per m2
  
  ; Detector parameters (Raytheon)
  nomic    = 0
  lmircam  = 1
  pix_size = 10.68          ; mas
  qe       = 0              ; ?, quantum efficiency
  ron      = 0              ; ?, read noise
  well_det = 1.1D7          ; e-
  min_adu  = [0,0]          ; Approximate minimum number of ADU in linear regime (assuming 1 co-add). For level shift and preamp respectively.
  max_adu  = [10000,48000]  ; Approximate maximum number of ADU in linear regime (assuming 1 co-add). For level shift and preamp respectively.
  adu2e    = [154.,1500.]   ; need to be revised for LMIRCam (but not critical)
  adu2ph   = FIX(adu2e/qe)
  x_chan   = 64             ; width of each channel in the X direction
  y_chan   = 1024           ; width of each channel in the Y direction
ENDELSE

; Compute size of the PSF
psf_rad  = 1.22*lambda/tel_diam              ; [rad], psf radius
psf_mas  = psf_rad*206264806.246D0           ; [mas]
psf_pix  = psf_mas/pix_size                  ; [pixels]
psf_fwhm = 1.028/1.22*psf_pix

; Output structure
cnf={ADU2E:adu2e, ADU2PH:adu2ph, BANDWIDTH:bandwidth, BASE:base, CEN_OBS:cen_obs, JY2PH:jy2ph, LAMBDA:lambda, LMIRCAM:lmircam, MAX_ADU: max_adu, MIN_ADU: min_adu, NOMIC:nomic, $
     PIX_SIZE:pix_size, PSF_RAD:psf_rad, PSF_MAS:psf_mas, PSF_PIX:psf_pix, PSF_FWHM:psf_fwhm, QE: qe, RON: ron, TEL_DIAM:tel_diam, X_CHAN: x_chan, Y_CHAN: y_chan, WELL_DET:well_det}    
END
