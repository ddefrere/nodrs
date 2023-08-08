;+
; NAME: LBTI_IMG2VIS
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
;   PLOT          :  Set this keyword to plot the data to eps files
;
; MODIFICATION HISTORY:
;   Version 1.0, 18-MAY-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (adapted from former routine 'remove_bck.pro')

PRO LBTI_IMG2VIS, img_in, hdr_in, LOG_FILE=log_file, INFO=info, NO_SAVE=no_save, PLOT=plot, OB_IN=ob_in

; Define operational parameters
COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs, log

; Keyword sanity check
IF NOT KEYWORD_SET(LOG_FILE) THEN lun = -1 ELSE lun = log_file
IF NOT KEYWORD_SET(OB_IN)    THEN ob_in =0

; Number of frames and pixels
n_img  = N_ELEMENTS(img_in[0,0,*])
n_ypix = N_ELEMENTS(img_in[0,*,0])
n_xpix = N_ELEMENTS(img_in[*,0,0])

; Running paramaters
obstype = hdr_in.header.obstype
n_spec  = 1024;256                  ; Size of the spectrum in pixels
n_wav   = n_spec
y_off   = 75                   ; Y offset between the measured spectrum center and the actual center
n_bin   = drs.n_bin > 1
n_clip  = drs.n_clip < n_xpix

; Number of beams in the image
IF obstype EQ 0 THEN n_beam = 2 ELSE n_beam = 1

; Find spectrum position based on median image
img_med = MEDIAN(img_in, DIMENSION=3)
pro_row = TOTAL(img_med, 1)
pro_col = TOTAL(img_med, 2)
x_spe   = TOTAL(INDGEN(n_xpix)*pro_col^2)/TOTAL(pro_col^2)
y_spe   = 512;TOTAL(INDGEN(n_ypix)*pro_row^2)/TOTAL(pro_row^2)+y_off

; Resample the time on an evenly spaced grid and compute the frequency coordinates (see FFT manual)
x      = (FINDGEN((n_clip - 1)/2) + 1)
period = 1  ; [pix]
IF n_clip MOD 2 THEN freq = [-(n_clip/2 + 1) + x, 0.0, x]/(n_clip*period) $
                ELSE freq = [-n_clip/2 + x,0.0, x, n_clip/2]/(n_clip*period)
                
; Define wav vector (apprximated)
pix2m   = hdr_in.header.bandwidth/n_spec
wav_pix = (n_spec-1)*DINDGEN(n_wav)/(n_wav-1)
wav_m   = hdr_in.header.lam_cen - 0.5*hdr_in.header.bandwidth + wav_pix*pix2m
                
; Prepare output data
img = FLTARR(n_clip, n_spec, n_img, n_beam, /NOZERO)
vis = FLTARR(n_img, n_wav, n_beam)
phi = FLTARR(n_img, n_wav, n_beam)

; Loop over the frames and compute phase
; idx0 = [INDGEN(0.5*n_clip-4),0.5*n_clip+4+INDGEN(0.5*n_clip-4)]
FOR i_img = 0, n_img-1 DO BEGIN
  FOR i_b = 0, n_beam-1 DO BEGIN
    ; Extract sub image
    img_crp        = EXTRAC(REFORM(img_in[*,*,i_img]), x_spe-0.5*n_clip, y_spe-0.5*n_spec, n_clip, n_spec)
    ;img_crp[idx0,*] = 0
    img[0,0,i_img,i_b] = img_crp    
    ; Compute 1-d FFT
    img_tmp = SHIFT(img_crp, 0.5*n_clip, 0)
    img_fft = FFT(img_tmp, DIMENSION=1, /DOUBLE)
    img_fft = SHIFT(img_fft, 0.5*n_clip, 0)
    ; Loop over the spectral channels
    FOR i_wav=0, n_wav-1 DO BEGIN
      fft_wav = img_fft[*,wav_pix[i_wav]]
      fft_abs = ABS(fft_wav)
      fft_pha = ATAN(Imaginary(fft_wav),Real_part(fft_wav))
      idx_ok  = WHERE(freq GT 0.1)
      idx_max = WHERE(fft_abs[idx_ok] EQ MAX(fft_abs[idx_ok]))
      vis[i_img,i_wav,i_b] = fft_abs[idx_ok[idx_max[0]]]
      phi[i_img,i_wav,i_b] = fft_pha[idx_ok[idx_max[0]]]
    ENDFOR
    ; Shift spectrum
  ENDFOR
ENDFOR

; Save results
IF NOT KEYWORD_SET(NO_SAVE) THEN BEGIN
  ; Derive the number of OBs
  ob_uniq = hdr_in.data(UNIQ(hdr_in.data.ob_id, SORT(hdr_in.data.ob_id))).ob_id
  n_ob    = N_ELEMENTS(ob_uniq) ; remove background and photometry frames

  ; Loop over the OBs
  FOR i_ob=0, n_ob-1 DO BEGIN
    ; Derive ob files
    idx_ob = WHERE(hdr_in.data.ob_id EQ ob_uniq[i_ob], n_frob)
    IF n_frob LT 1 THEN GOTO, SKIP_OB

    ; Extract OB data
    hdr_ob  = hdr_in.header
    data_ob = hdr_in.data[idx_ob]

    ; Photometric OBs
    IF hdr_ob.obstype EQ 1 THEN BEGIN
      ; Keep only closed frames (AO and PHASE)
      idx_keep = WHERE(data_ob.dloopon EQ 1 AND data_ob.sloopon EQ 1, n_keep)
      IF n_keep LE 0 THEN GOTO, SKIP_VIS
      ; Loop over the 2 beams
      FOR i_b=0, n_beam-1 DO BEGIN
        ; Store beam ID
        hdr_ob.beam_id = i_b
        ; Extract data t save
        vis_ob = vis[idx_keep,*,i_b]
        phi_ob = phi[idx_keep,*,i_b]
        ; Data to save
        img_vis  = REFORM(img[*,*,idx_keep,i_b])
        data_vis = data_ob[idx_keep]
        vis_meas = {WAV: wav_m, VIS: vis_ob, PHASE: phi_ob}
        ; Save data
        LBTI_SAVEL1IMG, img_vis, hdr_ob, data_vis, FILE_ID=ob_uniq[i_ob]
        LBTI_SAVEL1VIS, vis_meas, hdr_ob, data_vis, FILE_ID=ob_uniq[i_ob]
      ENDFOR
    ENDIF
    ; Jump point if no phot frame
    SKIP_VIS:
    SKIP_OB:
  ENDFOR
ENDIF

END
