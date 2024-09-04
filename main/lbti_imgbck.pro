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
;   Version 3.2,  04-SEP-2024, DD: create different folder for different IMG_MODE
;
pro LBTI_IMGBCK, img_in, hdr_in, bckg_path = bckg_path, log_file = log_file, no_save = no_save, info = info, plot = plot
  compile_opt idl2
  ;
  ; Define operational parameters
  common GLOBAL, prm, cnf, wav, tgt, pth, drs
  on_error, 0

  ; At this point, data might still be LONG or INT if no flat fielding
  ; Convert to FLOAT or DOUBLE. Make a copy because we don't want to modify img_in
  if drs.precision then img_cur = double(img_in) $
  else img_cur = float(img_in)

  ; Extract background data
  idx_bck = where(hdr_in.bck_nod ne 0, complement = idx_nod, n_bck)
  if n_bck le 0 then begin
    message, 'No background frames found!', /continue
    skip_bck = 1
  endif else if drs.bckg_mode eq 0 then skip_bck = 1 else skip_bck = 0
  if n_elements(idx_nod) le 0 then message, 'No science frames found!!!'

  ; Extract nod and background HDR data
  hdr_nod = hdr_in[idx_nod]
  hdr_bck = hdr_in[idx_bck]

  ; Compute nod "parity"
  nod_id = hdr_nod[0].nod_id
  idx_tmp = where(hdr_bck.nod_id lt nod_id, nl)
  idx_tmp = where(hdr_bck.nod_id gt nod_id, nu)
  hdr_nod.bck_nod = (-nl + nu) / float(nl + nu)

  ; Compute mean/median if negative background mode
  if keyword_set(bckg_path) then begin
    img_bck = LBTI_READDATA(bckg_path, hdr_data = hdr_bck, info = info)
    med_bck = LBTI_IMGMED(img_bck, hdr_data = hdr_bck, mean = drs.img_mode, precision = drs.precision, info = 0)
  endif else begin
    if drs.bckg_mode le 0 then med_bck = LBTI_IMGMED(img_cur[*, *, idx_bck], hdr_data = hdr_bck, mean = drs.img_mode, precision = drs.precision, info = 0) $
    else med_bck = img_cur[*, *, idx_bck]
  endelse

  ; If bckg_mode is 0, don't do background subtraction
  if drs.bckg_mode eq 0 then hdr_nod.n_frbck = 0 else hdr_nod.n_frbck = n_bck

  ; Derive the number of instrumental configurations of the current nod
  cfg_id = hdr_nod.cfg_id
  cfg_uniq = cfg_id[uniq(cfg_id, sort(cfg_id))]
  n_cfg = n_elements(cfg_uniq)

  ; Loop over the instrumental configurations
  for i_cfg = 0, n_cfg - 1 do begin
    ; Extract data of the current configuration
    idx_cfg = where(cfg_id eq cfg_uniq[i_cfg], nfr, /null)
    hdr_cfg = hdr_nod[idx_cfg] ; HDR for image of the current nod and config
    idx_img = idx_nod[idx_cfg] ; IDX for image of the current nod and config

    ; Save RAW image cube
    if not keyword_set(no_save) then LBTI_SAVEL0RED, img_cur[*, *, idx_img], hdr_cfg, sub_dir = 'raw', /savemedian ; (don't save median to speed up code)

    ; Now perform background subtraction and find beam position
    if not skip_bck then begin
      ; Index of background frames of this configuration
      idx_tmp = where(hdr_cfg[0].cfg_id eq hdr_bck.cfg_id, na)

      ; If background frames found, process
      if na gt 0 then begin
        ; Loop over the frames, subtract nod image, and clean frame
        ; i_fr MOD na to account for both background modes
        for ifr = 0, nfr - 1 do img_cur[0, 0, idx_img[ifr]] = LBTI_CLEANIMG((img_cur[*, *, idx_img[ifr]] - med_bck[*, *, idx_tmp[ifr mod na]]), nomic = cnf.nomic, lmircam = cnf.lmircam)

        ; Compute mean nodding frequency (in seconds)...TO BE CORRECTED
        hdr_cfg.nod_frq = abs(mean(hdr_cfg.mjd_obs) - mean(hdr_bck[idx_tmp].mjd_obs)) * 24 * 60 * 60

        ; Find beam for this configuration and store results (if not a background frame)
        ; Good centroid must be computed on close-loop frames!!! (but we save all images for diagnostic)
        ; THIS SHOULD BE IMPROVED (SEE LBTI_MASTERLOG AND HOW TO FLAG OPEN LOOP FRAMES)
        obstype = hdr_cfg[0].obstype
        case obstype of
          ; 0: PHOTOMETRY. We request that both AO loops are closed
          0: idx_keep = where(hdr_cfg.dloopon eq 1 and hdr_cfg.sloopon eq 1, nfr)
          ; 1: FIZEAU. We request that both AO loops are closed
          1: idx_keep = where(hdr_cfg.dloopon eq 1 and hdr_cfg.sloopon eq 1, nfr)
          ; 2: NULLING. We request that both AO loops and the phase loop are closed
          2: idx_keep = where(hdr_cfg.dloopon eq 1 and hdr_cfg.sloopon eq 1 and hdr_cfg.pcclosed eq 1, nfr)
          ; 3: BACKGROUND. No conditions
          3: idx_keep = lindgen(nfr)
          ; 4: TEST. No conditions
          4: idx_keep = lindgen(nfr)
          else: message, 'Unknown OBSTYPE'
        endcase

        ; If enough closed-loop frames, find centroid
        if nfr le 0 then begin
          message, 'Not enough closed-loop frames!! This nod will be discarded from now on.', /continue
          hdr_cfg[*].xcen[*] = 0
          hdr_cfg[*].ycen[*] = 0
        endif else begin
          ; Define ROIs
          if max(hdr_cfg[0].ROI) eq 0 then roi = [0, 0, n_elements(img_cur[*, 0, 0]) - 1, n_elements(img_cur[0, *, 0]) - 1] else roi = hdr_cfg[0].ROI
          ; Compute MEDIAN
          if nfr gt 1 then img_med = median(img_cur[roi[0] : roi[2], roi[1] : roi[3], idx_img[idx_keep]], dimension = 3) else img_med = img_cur[roi[0] : roi[2], roi[1] : roi[3], idx_img[idx_keep]]
          ; Compute number of beams
          if obstype eq 0 and drs.overlap eq 0 then n_beam = 2 else n_beam = 1
          ; run beam finding routine
          if drs.xcen[0] eq 0 and drs.ycen[0] eq 0 and drs.no_find ne 1 and obstype ne 3 then begin
            data = IMG_FINDBEAM(img_med, cnf.psf_pix, auto_beam = drs.auto_beam, cursor = drs.cursor, edge = [cnf.x_chan, cnf.y_chan], fit_mode = 3, $
              n_beam = n_beam, pixsize = cnf.pix_size, info = info, plot = plot)
            for i_beam = 0, n_beam - 1 do begin
              hdr_cfg[*].xcen[i_beam] = data.beam_pos[i_beam, 0] + roi[0]
              hdr_cfg[*].ycen[i_beam] = data.beam_pos[i_beam, 1] + roi[1]
              hdr_cfg[*].fwhm_x[i_beam] = data.fwhm_fit[i_beam, 0]
              hdr_cfg[*].fwhm_y[i_beam] = data.fwhm_fit[i_beam, 1]
              hdr_cfg[*].slope[i_beam] = data.slope[i_beam]
            endfor
          endif else begin
            hdr_cfg[*].xcen = drs.xcen
            hdr_cfg[*].ycen = drs.ycen
          endelse
        endelse ; End on the number of closed-loop frames
      endif else message, 'No corresponding background frames for file ' + hdr_cfg[0].filename, /continue

      ; Save background-subtracted image cube
      case drs.img_mode of  ;Image combination mode (0: median, 1: mean, 2: resistant mean)
        0: sub_dir = 'med'
        1: sub_dir = 'mean'
        2: sub_dir = 'resm'
        else: message, 'Unknown IMG_MODE'
      endcase
      if not keyword_set(no_save) then LBTI_SAVEL0RED, img_cur[*, *, idx_img], hdr_cfg, sub_dir = sub_dir, /savemedian ; (don't save median to speed up code)
    endif

    ; Update datalog file
    sav_path = pth.l0Fits_path + drs.date_obs + pth.sep
    if not file_test(sav_path) then file_mkdir, sav_path
    LBTI_DATALOG, sav_path + 'datalog.sav', hdr_cfg
    WRITE_DATALOG, drs.date_obs, instrum = drs.instrum
  endfor ; end of the loop on the instrumental setups of this NOD
end