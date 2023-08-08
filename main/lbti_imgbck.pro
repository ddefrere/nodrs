;+
; NAME: LBTI_IMGBCK
; 
; PURPOSE:
;   This procedure performs the chop/nod background subtraction and save the data.
;
; INPUT:
;   img_in        :  An image data cube
;   hdr_in        :  Corresponding header information
;
; KEYWORDS
;   BCKG_PATH     :
;   LOG_FILE      :  Set this keyword the path of a file where the log will be printed
;   NO_SAVE       :  Set this keyword to not save the images
;   INFO          :  Define the level of information printed to screen:
;                       - 0: completely silent execution;
;                       - 1: minimum level of information;
;                       - 2: nominal level of information;
;                       - 3: debugging use.
;   PLOT          :  Set this keyword to plot the data to eps files
;
; OUTPUT
;   Data cube with the background-subtracted images
;
; MODIFICATION HISTORY:
;   Version 1.0,  23-FEB-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (based on former routine "lbti_imgcen.pro")
;   Version 1.1,  23-MAY-2014, DD: removed keyword BCKG_MODE now passed through the drs structure
;   Version 1.2,  18-AUG-2014, DD: unused frames are now properly removed from the output structure
;   Version 1.3,  19-SEP-2014, DD: removed part that finds the star
;   Version 1.4,  01-NOV-2014, DD: added keyword NO_MEDIAN + improved speed
;   Version 1.5,  05-NOV-2014, DD: now find beam for each instrumental setup seperately + added keyword NO_SAVE
;   Version 1.6,  07-NOV-2014, DD: added slope to output header
;   Version 1.7,  11-NOV-2014, DD: added PHASECam setpoint in the definition of instrumental configurations
;   Version 1.8,  22-NOV-2014, DD: removed keyword NO_MEDIAN
;   Version 1.9,  04-APR-2015, DD: Added nodding period
;   Version 2.0,  04-OCT-2015, DD: Added pre-crop function
;   Version 2.1,  10-OCT-2015, DD: Added pactcdly to config definition
;   Version 2.2,  29-OCT-2015, DD: Removed definition of config ID (now assign right from the start and pass through header information)
;   Version 2.3,  21-DEC-2015, DD: Minor bug corrected + added quick mode
;   Version 2.4,  24-DEC-2015, DD: Updated to decrease memory usage + removed obsolete keyword HDR_OUT
;   Version 2.5,  21-JAN-2016, DD: Added best-fit FWHM to output data
;   Version 2.6,  28-MAR-2016, DD: Removed cropping, now done earlier!
;   Version 2.7,  25-APR-2016, DD: Now BCKG_MODE=0 means no background subtraction!
;   Version 2.8,  06-JUL-2016, DD: input arrays are not modified anymore (more memory consuming but speed up LBTI_DRS.pro)
;   Version 2.9,  15-JUL-2016, DD: added cursor
;   Version 3.0,  11-FEB-2017, DD: now save both raw and background-subtracted cubes
;   Version 3.1,  11-APR-2017, DD: implemented ROIs

PRO LBTI_IMGBCK, img_in, hdr_in, BCKG_PATH=bckg_path,  LOG_FILE=log_file, NO_SAVE=no_save, INFO=info, PLOT=plot  
                                                                                                                                     ;     
; Define operational parameters
COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs
ON_ERROR, 0

; At this point, data might still be LONG or INT if no flat fielding
; Convert to FLOAT or DOUBLE. Make a copy because we don't want to modify img_in
IF drs.precision THEN img_cur = DOUBLE(img_in) $
                 ELSE img_cur = FLOAT(img_in)

; Extract background data
idx_bck = WHERE(hdr_in.bck_nod NE 0, COMPLEMENT=idx_nod, n_bck)
IF n_bck LE 0 THEN BEGIN
  MESSAGE, 'No background frames found!', /CONTINUE
  skip_bck = 1
ENDIF ELSE IF drs.bckg_mode EQ 0 THEN skip_bck = 1 ELSE skip_bck = 0
IF N_ELEMENTS(idx_nod) LE 0 THEN MESSAGE, 'No science frames found!!!'
 
; Extract nod and background HDR data
hdr_nod = hdr_in[idx_nod]
hdr_bck = hdr_in[idx_bck]

; Compute nod "parity" 
nod_id  = hdr_nod[0].nod_id
idx_tmp = WHERE(hdr_bck.nod_id LT nod_id, nl)
idx_tmp = WHERE(hdr_bck.nod_id GT nod_id, nu)
hdr_nod.bck_nod = (-nl+nu)/FLOAT(nl+nu)

; Compute mean/median if negative background mode
IF KEYWORD_SET(bckg_path) THEN BEGIN
   img_bck = LBTI_READDATA(bckg_path, HDR_DATA=hdr_bck, INFO=info)
   med_bck = LBTI_IMGMED(img_bck, HDR_DATA=hdr_bck, MEAN=drs.img_mode, PRECISION=drs.precision, INFO=0)
ENDIF ELSE BEGIN
  IF drs.bckg_mode LE 0 THEN med_bck = LBTI_IMGMED(img_cur[*,*,idx_bck], HDR_DATA=hdr_bck, MEAN=drs.img_mode, PRECISION=drs.precision, INFO=0) $
                        ELSE med_bck = img_cur[*,*,idx_bck]
ENDELSE

; If bckg_mode is 0, don't do background subtraction
IF drs.bckg_mode EQ 0 THEN hdr_nod.n_frbck = 0 ELSE hdr_nod.n_frbck = n_bck
   
; Derive the number of instrumental configurations of the current nod
cfg_id   = hdr_nod.cfg_id
cfg_uniq = cfg_id[UNIQ(cfg_id,  SORT(cfg_id))]
n_cfg    = N_ELEMENTS(cfg_uniq)

; Loop over the instrumental configurations
FOR i_cfg = 0, n_cfg-1 DO BEGIN
  
  ; Extract data of the current configuration
  idx_cfg = WHERE(cfg_id EQ cfg_uniq[i_cfg], nfr, /NULL)
  hdr_cfg = hdr_nod[idx_cfg]     ; HDR for image of the current nod and config
  idx_img = idx_nod[idx_cfg]     ; IDX for image of the current nod and config
    
  ; Save RAW image cube
  IF NOT KEYWORD_SET(NO_SAVE) THEN LBTI_SAVEL0RED, img_cur[*,*,idx_img], hdr_cfg, SUB_DIR='raw', /SAVEMEDIAN ;(don't save median to speed up code)
  
  ; Now perform background subtraction and find beam position
  IF NOT skip_bck THEN BEGIN
    ; Index of background frames of this configuration
    idx_tmp = WHERE(hdr_cfg[0].cfg_id EQ hdr_bck.cfg_id, na)
    
    ; If background frames found, process
    IF na GT 0 THEN BEGIN  
      ; Loop over the frames, subtract nod image, and clean frame
      ; i_fr MOD na to account for both background modes
      FOR ifr = 0, nfr-1 DO img_cur[0,0,idx_img[ifr]] = LBTI_CLEANIMG((img_cur[*,*,idx_img[ifr]]-med_bck[*,*,idx_tmp[ifr MOD na]]), NOMIC=cnf.nomic, LMIRCAM=cnf.lmircam) 
                                         
      ; Compute mean nodding frequency (in seconds)...TO BE CORRECTED 
      hdr_cfg.nod_frq = ABS(MEAN(hdr_cfg.mjd_obs)-MEAN(hdr_bck[idx_tmp].mjd_obs))*24*60*60
      
      ; Find beam for this configuration and store results (if not a background frame)  
      ; Good centroid must be computed on close-loop frames!!! (but we save all images for diagnostic)
      ; THIS SHOULD BE IMPROVED (SEE LBTI_MASTERLOG AND HOW TO FLAG OPEN LOOP FRAMES)
      obstype = hdr_cfg[0].obstype
      CASE obstype  OF
        ; 0: PHOTOMETRY. We request that both AO loops are closed
        0: idx_keep = WHERE(hdr_cfg.dloopon EQ 1 AND hdr_cfg.sloopon EQ 1, nfr)
        ; 1: FIZEAU. We request that both AO loops are closed
        1: idx_keep = WHERE(hdr_cfg.dloopon EQ 1 AND hdr_cfg.sloopon EQ 1, nfr)
        ; 2: NULLING. We request that both AO loops and the phase loop are closed
        2: idx_keep = WHERE(hdr_cfg.dloopon EQ 1 AND hdr_cfg.sloopon EQ 1 AND hdr_cfg.pcclosed EQ 1, nfr)
        ; 3: BACKGROUND. No conditions
        3: idx_keep = LINDGEN(nfr)
        ; 4: TEST. No conditions
        4: idx_keep = LINDGEN(nfr)
        ELSE: MESSAGE, 'Unknown OBSTYPE'
      ENDCASE 
      
      ; If enough closed-loop frames, find centroid
      IF nfr LE 0 THEN BEGIN
        MESSAGE, 'Not enough closed-loop frames!! This nod will be discarded from now on.', /CONTINUE
        hdr_cfg[*].xcen[*] = 0
        hdr_cfg[*].ycen[*] = 0
      ENDIF ELSE BEGIN
        ; Define ROIs
        IF MAX(hdr_cfg[0].roi) EQ 0 THEN roi = [0,0,N_ELEMENTS(img_cur[*,0,0])-1,N_ELEMENTS(img_cur[0,*,0])-1] ELSE roi = hdr_cfg[0].roi
        ; Compute MEDIAN
        IF nfr GT 1 THEN img_med = MEDIAN(img_cur[roi[0]:roi[2],roi[1]:roi[3],idx_img[idx_keep]], DIMENSION=3) ELSE img_med = img_cur[roi[0]:roi[2],roi[1]:roi[3],idx_img[idx_keep]]
        ; Compute number of beams
        IF obstype EQ 0 AND drs.overlap EQ 0 THEN n_beam = 2 ELSE n_beam = 1
        ; run beam finding routine
        IF drs.xcen[0] EQ 0 AND drs.ycen[0] EQ 0 AND drs.no_find NE 1 AND obstype NE 3 THEN BEGIN
          data = IMG_FINDBEAM(img_med, cnf.psf_pix, AUTO_BEAM=drs.auto_beam, CURSOR=drs.cursor, EDGE=[cnf.x_chan,cnf.y_chan], FIT_MODE=3, $
                              N_BEAM=n_beam, PIXSIZE=cnf.pix_size, INFO=info, PLOT=plot)
          FOR i_beam=0, n_beam-1 DO BEGIN
            hdr_cfg[*].xcen[i_beam]   = data.beam_pos[i_beam,0] + roi[0]
            hdr_cfg[*].ycen[i_beam]   = data.beam_pos[i_beam,1] + roi[1]
            hdr_cfg[*].fwhm_x[i_beam] = data.fwhm_fit[i_beam,0]
            hdr_cfg[*].fwhm_y[i_beam] = data.fwhm_fit[i_beam,1]
            hdr_cfg[*].slope[i_beam]  = data.slope[i_beam]
          ENDFOR
        ENDIF ELSE BEGIN
          hdr_cfg[*].xcen = drs.xcen
          hdr_cfg[*].ycen = drs.ycen
        ENDELSE
      ENDELSE ; End on the number of closed-loop frames
    ENDIF ELSE MESSAGE, "No corresponding background frames for file " + hdr_cfg[0].filename, /CONTINUE
    
    ; Save background-subtracted image cube
    IF NOT KEYWORD_SET(NO_SAVE) THEN LBTI_SAVEL0RED, img_cur[*,*,idx_img], hdr_cfg, SUB_DIR='bckg', /SAVEMEDIAN ;(don't save median to speed up code)
  ENDIF 

  ; Update datalog file
  sav_path = pth.l0fits_path + drs.date_obs + pth.sep
  IF NOT FILE_TEST(sav_path) THEN FILE_MKDIR, sav_path
  LBTI_DATALOG, sav_path + 'datalog.sav', hdr_cfg
  WRITE_DATALOG, drs.date_obs, INSTRUM=drs.instrum
ENDFOR ; end of the loop on the instrumental setups of this NOD 
END