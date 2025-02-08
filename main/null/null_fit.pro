;+
; NAME: NULL_FIT
;
; PURPOSE:
;   This function fit null data given as input.
;
; INPUTS:
;   data_path     :  Path to the level 2 fits files (calibrated sequences of null).
;
; KEYWORDS
;    L1TYPE        : set this keyword to the type of L1 file to fit (OB, PT, or BIAS). Pt is default.
;    PLOT          : set this keyword to have plot the results
;    PS            : set this keyword to save the plot in an eps figure
;
; MODIFICATION HISTORY:
;   Version 1.0,  16-JAN-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  14-APR-2015, DD: corrected bug on statistical error computation
;   Version 1.2,  14-APR-2015, DD: modified for ORR version + added LOG_FILE
;   Version 1.3,  14-APR-2015, DD: now read HA from L! file rather than computing it
;   Version 1.4,  04-NOV-2016, DD: added keyword L1 type
;   Version 1.5,  04-APR-2017, DD: adapted for new formlaism of GET_TGT
;   Version 1.6,  03-AUG-2017, DD: improved output formalism
;   Version 1.7,  05-AUG-2017, DD: now use AVGERR instead of AVGSDV
;   Version 1.8,  29-SEP-2017, DD: updated touse the new keywords in the L2 files (i.e., NULL_EST and SIG_SCA)
;   Version 1.9,  29-JUL-2018, DD: added comments to error computation
;   Version 2.0,  31-JUL-2018, DD: improved consitency with NULL_CALIB.pro
;   Version 2.1,  24-MAY-2024, DD: Update file permission

pro NULL_FIT, date, cfg_file, info = info, log_file = log_file, l1type = l1type, plot = plot, ps = ps
  compile_opt idl2

  ; Start actual code
  on_error, 0
  if not keyword_set(info) then info = 1
  drs_version = 2.1
  drs_date = '24-MAY-2024'

  ; DEFINE GLOBAL VARIABLES
  ; ***********************

  ; Astronomical and physical constants
  GET_PRM, prm

  ; Obtain the definition of the configuration
  GET_CNF, cnf, instrum = 'nomic'

  ; Read config file with reduction parameters
  GET_DRS, drs, 'nodrs/cfg/' + cfg_file

  ; Recover the IDL running path
  DECLARE_PATH, pth, instrum = 'nomic'

  ; COMPUTE FILE PATH
  ; *****************

  ; Get long version of the date
  date_lng = '20' + strmid(date, 0, 2) + '-' + strmid(date, 2, 2) + '-' + strmid(date, 4, 2)

  ; Define operational parameters
  common LBTI_PARAM, tel_diam, base, lambda, pix_size, psf_pix

  ; L2 directory
  if max(drs.null_rad) ne 0 and not strmatch(drs.dir_label, '*_APR') then drs.dir_label = drs.dir_label + '_APR'
  l2_dir = pth.l2Fits_path + pth.sep + date_lng + drs.dir_label + pth.sep ;+ '_APR' + pth.sep
  if not file_test(l2_dir) then file_mkdir, l2_dir

  ; Init log file
  if keyword_set(log_file) then begin
    log_file = l2_dir + date_lng + '.txt'
    openw, lun, log_file, /get_lun, width = 800, /append
    printf, lun, ' '
    printf, lun, '******************************************************************************************************* '
    printf, lun, ' '
    printf, lun, 'NULL_FIT.pro version ' + string(drs_version, format = '(F3.1)') + ' -- ' + drs_date + ' -- Denis Defrère - Steward Observatory (denis@lbti.org)'
    printf, lun, 'Data fitted on ' + systime()
    printf, lun, ' '
    printf, lun, '******************************************************************************************************* '
    printf, lun, ' '
  endif else lun = -1

  ; START CALIBRATION
  ; *****************

  ; Search for L1 file and read it
  if drs.cal_mode eq 1 then file_type = 'OB' else file_type = 'PT'
  if keyword_set(l1type) then file_type = file_type + '-' + l1type
  l2file = l2_dir + 'UT' + date_lng + '_calib_' + file_type + '.fits'
  if not file_test(l2file) then RETURN ; NOT FOUND
  l2data = MRDFITS(l2file, 1, header, /silent)
  n_data = n_elements(l2data.mjd_obs)
  tmp = MRDFITS(l2file, 0, header, /silent)

  ; Read header information
  null_est = FXPAR(header, 'NULL_EST', /nocontinue)
  sig_sca = FXPAR(header, 'SIG_SCA', /nocontinue)
  cal_mode = FXPAR(header, 'CAL_MODE', /nocontinue)
  file_ver = FXPAR(header, 'FILE_VER', /nocontinue)
  ; bck_mode = FXPAR(header, 'BCK_MODE', /NOCONTINUE)

  ; Print some info to screen
  aper_rad = l2data.aper_rad
  bck_irad = l2data.bck_irad
  bck_orad = l2data.bck_orad
  flag = l2data.quality_flg
  ; IF info GT 0 THEN PRINT, 'Aperture radius : ' + STRING(aper_rad, FORMAT='(I0)')

  ; Read data
  wav = mean(l2data.wav_eff)
  bdw = mean(l2data.bandwidth)
  tgt = l2data.objname

  ; Wavelength range
  n_lam = 20
  lambda_chnl = wav - 0.5 * bdw + bdw * (1. + 1. / n_lam) * dindgen(n_lam) / (n_lam - 1)

  ; Extract LBTI transmission
  trans = 1. ; Temporary implementation

  ; Find unique targets
  tgt_uniq = tgt[uniq(tgt, sort(tgt))]
  n_tgt = n_elements(tgt_uniq)

  ; Loop over the targets
  for i_tgt = 0, n_tgt - 1 do begin
    ; extract target
    idx_tgt = where(tgt eq tgt_uniq[i_tgt], n_sci)
    sci_time = l2data[idx_tgt].mjd_obs
    sci_para = l2data[idx_tgt].lbt_para
    sci_ha = l2data[idx_tgt].hr_angle
    idx_neg = where(sci_para gt -180 and sci_para lt -90, nn)
    ; IF nn GT 0 THEN sci_para[idx_neg] += 360
    null_sci = l2data[idx_tgt].null_cal ; - 0.005
    null_esyst = l2data[idx_tgt].null_cal_sys
    null_estat = l2data[idx_tgt].null_cal_sta
    null_etot = l2data[idx_tgt].null_cal_err
    null_eraw = l2data[idx_tgt].null_meas_err
    null_flo_err = l2data[idx_tgt].null_flo_err

    ; Extract target information
    GET_TGT, tgt_uniq[i_tgt], star, database = pth.input_path + 'hosts_db.dat'
    th_mean = star.ldm
    th_min = star.ldm - star.lde
    th_max = star.ldm + star.lde
    uld = star.uld

    ; Convert the time to UT hours
    time_mean = mean(sci_time)
    DAYCNV, time_mean + 2400000.5d0, yr, mn, day, hr
    ; title_day = STRING(yr,FORMAT='(I0)') + '/' + STRING(mn,FORMAT='(I0)') + '/' + STRING(day,FORMAT='(I0)')
    JDCNV, yr, mn, day, 0., t0day
    t0day = t0day - 2400000.5d0
    scitime = sci_time - t0day

    ; Generate stellar spectrum for broadband null computation
    star_flux = BLACKBODY(star.temp, lambda_chnl, standard = 3) ; flux in photons/s/m²/Hz/sr

    ; Compute the expected null for this star
    norm = total((trans * star_flux) ^ 2)
    null_lam = GEOMETRIC_NULL(th_mean, cnf.base, lambda_chnl, uld = uld)
    null_wb = total((trans * star_flux) ^ 2 * null_lam) / norm
    ; Compute minimum null (with minimum diameter)
    null_min = GEOMETRIC_NULL(th_min, cnf.base, lambda_chnl, uld = uld)
    null_wb_min = total((trans * star_flux) ^ 2 * null_min) / norm
    ; Compute maximum null (with maximum diameter)
    null_max = GEOMETRIC_NULL(th_max, cnf.base, lambda_chnl, uld = uld)
    null_wb_max = total((trans * star_flux) ^ 2 * null_max) / norm
    ; Derive the error bar on the estimated null
    null_err = (null_wb_max - null_wb_min) / 2d0

    ; Compute statistical error bar on SCI OBs
    rms_stat = 1. / sqrt(total(1. / null_eraw ^ 2)) ; statistical error on SCI measurements

    ; Compute systematic and total error bar on TF
    rms_tf_tot = total(null_flo_err) / n_sci
    rms_syst = sqrt(total(null_esyst / n_sci) ^ 2 + null_err ^ 2)

    ; Compute weighted-average null and error bars
    if n_sci le 1 then begin
      null_avg = null_sci
      disp_m_null = 0
      disp_m_null_we = 0
      rms_exc = 0.
    endif else begin
      AVGERR, null_sci, null_etot, null_avg_unw, null_avg_wei, rms_stat_null, disp_m_null, disp_m_null_we, rms_exc, kappa = sig_sca, idx_out = idx_out, n_out = n_out
      if null_est eq 0 then null_avg = null_avg_unw else null_avg = null_avg_wei
    endelse

    ; Definition of systematic errors
    case (abs(drs.err_mode) mod 3) of
      0: rms_sys = disp_m_null
      1: rms_sys = disp_m_null_we
      2: rms_sys = rms_exc
      else: message, 'Undefined error mode'
    endcase

    ; Add errors (if drs.err_mode GT 3, use only the greater of the two)
    if drs.err_mode ge 3 then rms_tot = sqrt(null_err ^ 2 + rms_tf_tot ^ 2 + rms_stat ^ 2) > sqrt(null_err ^ 2 + rms_tf_tot ^ 2 + rms_sys ^ 2) $
    else rms_tot = sqrt(null_err ^ 2 + rms_tf_tot ^ 2 + rms_stat ^ 2 + rms_sys ^ 2) ; null err is the error on the stellar diameter, rms_tf_tot is the total error on the TF at that position, rms_stat is the stat error on the SCI OBs, and rms_sys is the SYS error

    ; Compute stat error on TF
    rms_tf_stat = mean(sqrt(null_estat ^ 2 - null_eraw ^ 2))

    ; Print results
    if info gt 0 then begin
      print, tgt_uniq[i_tgt] + ' [%]:'
      print, 'CALIBRATION PARAMETERS'
      print, '   - Photometric aperture  :', aper_rad[i_tgt]
      print, '   - Calibration mode      :', cal_mode
      print, '   - Null estimator        :', null_est
      print, '   - Sigma threshold limit :', sig_sca
      print, 'CALIBRATED NULL'
      print, '   Estimated null for photopshere         : ', string((null_wb) * 1d2, format = '(F8.5)') + ' +/- ' + string(null_err * 1d2, format = '(F7.5)')
      print, '   Calibrated null (photosphere corected) : ', string((null_avg - null_wb) * 1d2, format = '(F8.5)') + ' +/- ' + string(rms_tot * 1d2, format = '(F7.5)')
      print, '   - Statistical error on SCI OBs      : ', string(rms_stat * 1d2, format = '(F7.5)')
      print, '   - Scatter error on SCI OBs          : ', string(disp_m_null * 1d2, format = '(F7.5)')
      print, '   - Error due to diameter uncertainty : ', string(null_err * 1d2, format = '(F7.5)')
      print, '   - Error on TF                       : ', string(rms_tf_tot * 1d2, format = '(F7.5)')
      if lun gt 0 then begin
        printf, lun, ' '
        printf, lun, tgt_uniq[i_tgt] + ' [%]:'
        printf, lun, 'CALIBRATION PARAMETERS'
        printf, lun, '   - Photometric aperture  :', aper_rad[i_tgt]
        printf, lun, '   - Calibration mode      :', cal_mode
        printf, lun, '   - Null estimator        :', null_est
        printf, lun, '   - Sigma threshold limit :', sig_sca
        printf, lun, 'CALIBRATED NULL'
        printf, lun, '   Estimated null for photopshere         : ', string((null_wb) * 1d2, format = '(F8.5)') + ' +/- ' + string(null_err * 1d2, format = '(F7.5)')
        printf, lun, '   Calibrated null (photosphere corected) : ', string((null_avg - null_wb) * 1d2, format = '(F8.5)') + ' +/- ' + string(rms_tot * 1d2, format = '(F7.5)')
        printf, lun, '   - Statistical error on SCI OBs      : ', string(rms_stat * 1d2, format = '(F7.5)')
        printf, lun, '   - Scatter error on SCI OBs          : ', string(disp_m_null * 1d2, format = '(F7.5)')
        printf, lun, '   - Error due to diameter uncertainty : ', string(null_err * 1d2, format = '(F7.5)')
        printf, lun, '   - Error on TF                       : ', string(rms_tf_tot * 1d2, format = '(F7.5)')
      endif

      ; Write summary file
      if keyword_set(log_file) then begin
        ; Create log file
        ; log_file2 =  pth.l2fits_path + pth.sep + 'results_hosts.txt'
        log_file2 = pth.l2Fits_path + pth.sep + 'results_hosts2018.txt'
        if not file_test(log_file2) then begin
          openw, lun2, log_file2, /get_lun, width = 800, /append
          printf, lun2, 'Column signification'
          printf, lun2, '1: Date'
          printf, lun2, '2: Target'
          printf, lun2, '3: Aperture radius'
          printf, lun2, '4: Background outer radius'
          printf, lun2, '5: Background inner radius'
          printf, lun2, '6: Null estimator mode'
          printf, lun2, '7: Number of sigmas for scatter rejection of OBs'
          printf, lun2, '8: Calibration mode (0 per pointing; 1 per OB)'
          printf, lun2, '9: Calibrated null [%]'
          printf, lun2, '10: Error on calibrated null [%]'
          printf, lun2, '11: Quality flag'
          printf, lun2, '12: L1 file version'
          printf, lun2, '13: Execution date'
          printf, lun2, ' 1 ', '  ', '   2  ', '  ', '                  3     ', '  ', '   4   ', '  ', '   5   ', '  ', '   6   ', '  ', '   7   ', '  ', '   8   ', '  ', '   9   ', '  ', '   10   ', '   11   ', '  ', '   12   ', '   13   '
          printf, lun2, '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
        endif else openw, lun2, log_file2, /get_lun, width = 800, /append
        printf, lun2, date, '  ', tgt_uniq[i_tgt], aper_rad[idx_tgt[0]], '  ', bck_irad[idx_tgt[0]], '  ', bck_orad[idx_tgt[0]], '  ', null_est, '  ', sig_sca, '  ', cal_mode, '  ', $
          string((null_avg - null_wb) * 1d2, format = '(F8.5)'), '  ', string(rms_tot * 1d2, format = '(F7.5)'), '  ', flag[i_tgt], '  ', string(file_ver, format = '(I0)'), '  ', systime()
        close, lun2
        free_lun, lun2
      endif
    endif

    ; Fit a binary model (very simple approach!!)
    n_tet = 1.
    tet_max = 500. * prm.m2R
    n_phi = 1.
    phi_max = !dpi
    n_ctr = 100.
    ctr_max = 0.20
    tet = 252.21 * prm.m2R ; DINDGEN(n_tet)/(n_tet-1.) * tet_max  ; Angular separation
    phi = 244. * !dtor ; DINDGEN(n_phi)/(n_phi-1.) * phi_max  ; Position angle
    ctr = dindgen(n_ctr) / (n_ctr - 1.) * ctr_max
    ; n_tet = 50. & tet_max = 500. * prm.m2r
    ; n_phi = 50. & phi_max = !Dpi
    ; n_ctr = 50. & ctr_max = 0.20
    ; tet = DINDGEN(n_tet)/(n_tet-1.) * tet_max  ; Angular separation
    ; phi = DINDGEN(n_phi)/(n_phi-1.) * phi_max  ; Position angle
    ; ctr = DINDGEN(n_ctr)/(n_ctr-1.) * ctr_max

    ; Target visibilities
    v1 = (1 - abs(null_wb)) / (1 + abs(null_wb))
    v2 = 1.00

    ; Measured baseline position angles
    phi_base = atan(l2data[idx_tgt].u_coord, l2data[idx_tgt].v_coord)

    ; Predict baseline position angle based on hour angle (consitancy check)
    ra = star.ra
    dec = star.dec
    ra *= 15
    precess, ra, dec, 2000.0, 2014.0
    ra /= 15
    par_ang = parangle(sci_ha, dec, TEN(32, 42, 05.87))
    ; PRINT, par_ang-sci_para
    ; PRINT, sci_para

    ; Loop over the parameters
    r_chi2 = dblarr(n_tet, n_phi, n_ctr)
    for i_tet = 0, n_tet - 1 do begin
      for i_phi = 0, n_phi - 1 do begin
        for i_ctr = 0, n_ctr - 1 do begin
          ; tm                        = SIN(!Dpi*cnf.base*tet[i_tet]/wav*COS(phi_meas-phi[i_phi]))^2*ctr[i_ctr] ; old approach
          vis = (v1 + ctr[i_ctr] * v2 * exp(-dcomplex(0, 2 * !dpi * cnf.base * tet[i_tet] / wav * cos(phi[i_phi] - phi_base)))) / (1 + ctr[i_ctr])
          null = (1. - abs(vis)) / (1. + abs(vis))
          r_chi2[i_tet, i_phi, i_ctr] = total((null_sci - null) ^ 2 / null_etot ^ 2) / (n_elements(null_sci) - 2)
        endfor
      endfor
    endfor
    idx_min = where(r_chi2 eq min(r_chi2), n_min)
    idx_chi = array_indices(r_chi2, idx_min)
    ; PRINT, 'Best chi2: ', MIN(r_chi2)
    ; PRINT, 'Best ang. sep.:', tet[idx_chi[0]] / prm.m2r
    ; PRINT, 'Best PA      :', phi[idx_chi[1]] / !DTOR
    ; PRINT, 'Best contrast:', ctr[idx_chi[2]]
    ctr_opt = ctr[idx_chi[2]]

    ; Compute error bar
    r_chi2 = dblarr(n_tet, n_phi, n_ctr)
    null_sci_up = null_sci + null_etot
    for i_tet = 0, n_tet - 1 do begin
      for i_phi = 0, n_phi - 1 do begin
        for i_ctr = 0, n_ctr - 1 do begin
          ; tm                        = SIN(!Dpi*cnf.base*tet[i_tet]/wav*COS(phi_meas-phi[i_phi]))^2*ctr[i_ctr] ; old approach
          vis = (v1 + ctr[i_ctr] * v2 * exp(-dcomplex(0, 2 * !dpi * cnf.base * tet[i_tet] / wav * cos(phi[i_phi] - phi_base)))) / (1 + ctr[i_ctr])
          null = (1. - abs(vis)) / (1. + abs(vis))
          r_chi2[i_tet, i_phi, i_ctr] = total((null_sci_up - null) ^ 2 / null_etot ^ 2) / (n_elements(null_sci) - 2)
        endfor
      endfor
    endfor
    idx_min = where(r_chi2 eq min(r_chi2), n_min)
    idx_chi = array_indices(r_chi2, idx_min)
    ctr_up = ctr[idx_chi[2]]

    r_chi2 = dblarr(n_tet, n_phi, n_ctr)
    null_sci_do = null_sci - null_etot
    for i_tet = 0, n_tet - 1 do begin
      for i_phi = 0, n_phi - 1 do begin
        for i_ctr = 0, n_ctr - 1 do begin
          ; tm                        = SIN(!Dpi*cnf.base*tet[i_tet]/wav*COS(phi_meas-phi[i_phi]))^2*ctr[i_ctr] ; old approach
          vis = (v1 + ctr[i_ctr] * v2 * exp(-dcomplex(0, 2 * !dpi * cnf.base * tet[i_tet] / wav * cos(phi[i_phi] - phi_base)))) / (1 + ctr[i_ctr])
          null = (1. - abs(vis)) / (1. + abs(vis))
          r_chi2[i_tet, i_phi, i_ctr] = total((null_sci_do - null) ^ 2 / null_etot ^ 2) / (n_elements(null_sci) - 2)
        endfor
      endfor
    endfor
    idx_min = where(r_chi2 eq min(r_chi2), n_min)
    idx_chi = array_indices(r_chi2, idx_min)
    ctr_do = ctr[idx_chi[2]]

    ; Error bar
    print, 'Error on contrast', 0.5 * abs(ctr_do - ctr_up)

    ; Use theoretical values
    ctr_opt = 0.0355
    ctr_do = ctr_opt - 0.0022
    ctr_up = ctr_opt + 0.0022

    ; Plot the calibrated nulls
    if keyword_set(plot) then begin
      result_path = pth.result_path + pth.sep + 'fit' + pth.sep
      if not file_test(result_path) then file_mkdir, result_path
      plotname = result_path + 'UT' + date + '_' + strtrim(tgt_uniq[i_tgt], 2) + '_calib.eps'
      if keyword_set(ps) then begin
        PREP_PS, /bold
        device, filen = plotname, /encaps, /color, xsize = 17.8, ysize = 13.2, /times
      endif $
        else WINDOW, !d.window + 1, title = 'Calibrated null'
      loadct, 12, /silent
      pa_range = abs(max(sci_para) - min(sci_para))
      xrange = [min(sci_para) - 0.1 * pa_range, max(sci_para) + 0.1 * pa_range]
      yrange = [-0.5, 4.]
      PLOT, [0.], [0.], xtitle = 'Parallactic angle [deg]', ytitle = 'Calibrated null [%]', title = ' ', xstyle = 1, ystyle = 1, xrange = xrange, yrange = yrange, color = 0
      ; Overplot the estimated null values for calibrators measurements
      loadct, 13, /silent
      usersym, [1, 1, -1, -1, 1] * 0.7, [1, -1, -1, 1, 1] * 0.7, thick = 3.0, /fill
      oplot, [sci_para], 100. * [null_sci], psym = 8, color = 255
      errplot, [sci_para], 100. * ([null_sci] - [null_etot]), 100. * ([null_sci] + [null_etot]), color = 255
      loadct, 12, /silent
      ; Overplot the stellar photosphere
      n_plot = 1000
      n_surf = 1000
      plot_gleak = dblarr(n_plot, n_surf) ; Prepare surfaces for 3-sigma filled contour plot
      ysurf = dindgen(n_surf) / (n_surf - 1d0) * (yrange[1] - yrange[0]) + yrange[0]
      plot_pa = dindgen(n_plot) / (n_plot - 1d0) * (xrange[1] - xrange[0]) + xrange[0]
      for iplot = 0, n_plot - 1d0 do begin
        plot_gleak[iplot, *] = (abs(ysurf - 100. * null_wb) lt 3. * 100. * null_err)
      endfor
      oplot, plot_pa, replicate(100. * null_wb, n_plot), color = 0., linestyle = 2
      CONTOUR, plot_gleak, plot_pa, ysurf, /cell_fill, /overplot, c_colors = [100], levels = [1]
      ; Overplot the best model
      phi_plot = 0.5 * !dpi + plot_pa * !dtor
      idx_neg = where(plot_pa gt 0, n_neg)
      if n_neg gt 0 then phi_plot[idx_neg] = -0.5 * !dpi + plot_pa[idx_neg] * !dtor
      vis = (v1 + ctr_opt * v2 * exp(-dcomplex(0, 2 * !dpi * cnf.base * tet[idx_chi[0]] / wav * cos(phi[idx_chi[1]] - phi_plot)))) / (1 + ctr_opt)
      bin_fit = (1. - abs(vis)) / (1. + abs(vis))
      oplot, plot_pa, 100. * bin_fit, linestyle = 0, color = 100
      vis = (v1 + ctr_up * v2 * exp(-dcomplex(0, 2 * !dpi * cnf.base * tet[idx_chi[0]] / wav * cos(phi[idx_chi[1]] - phi_plot)))) / (1 + ctr_up)
      bin_fit = (1. - abs(vis)) / (1. + abs(vis))
      oplot, plot_pa, 100. * bin_fit, linestyle = 1, color = 100
      vis = (v1 + ctr_do * v2 * exp(-dcomplex(0, 2 * !dpi * cnf.base * tet[idx_chi[0]] / wav * cos(phi[idx_chi[1]] - phi_plot)))) / (1 + ctr_do)
      bin_fit = (1. - abs(vis)) / (1. + abs(vis))
      oplot, plot_pa, 100. * bin_fit, linestyle = 1, color = 100
      ; Print info to plot
      ; xpos0 = 0.04*(xrange[1]-xrange[0]) + xrange[0]
      ; ypos0 = 0.95*(yrange[1]-yrange[0]) + yrange[0]
      ; XYOUTS, xpos0, ypos0*0.98, 'Angular separation = ' + STRING(tet[idx_chi[1]]/m2r,FORMAT='(I0)') + ' mas'
      ; XYOUTS, xpos0, ypos0*0.92, 'Contrast = ' + STRING(ctr[idx_chi[2]],FORMAT='(F5.2)') + '%'
      ; XYOUTS, xpos0, ypos0*0.86, '!4v!3!dr!u2!n = ' + STRING(MIN(r_chi2),FORMAT='(F6.2)')
      loadct, 0, /silent
      if keyword_set(ps) then begin
        device, /close
        END_PS
      endif
    endif
  endfor

  if keyword_set(log_file) then begin
    printf, lun, ' '
    printf, lun, ' '
    close, lun
    free_lun, lun
  endif
end
