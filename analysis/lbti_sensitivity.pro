; +
; NAME: LBTI_SENSITIVITY
;   
; PURPOSE:
;   Main sensitivity analysis routine. Sensitivity is computed as the standard deviation of a series of flux measurements,
;   assuming nulling mode and nulling data reduction.
;
;   Contain sub-routine TEST_SENSITIVITY2 used for the nulling analysis paper. 
;
; KEYWORDS:
;   DATA_OUT : Output structure with various infos on the files
;   FILES    : Files to process
;   FLUX     : Flux of the star in Jy (superseed automatic computation)
;   INFO     : Print info to screen
;   LOG      : Log info to file
;
; MODIFICATION HISTORY:
;   Version 1.0,  04-APR-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  04-APR-2014, DD: added factor 2 for background error
;   Version 1.2,  20-AUG-2015, DD: added keyword DATA_OUT
;   Version 1.3,  25-OCT-2016, DD: added keyword LOG

FUNCTION LBTI_SENSITIVITY, DATA_OUT=data_out, FILES=files, FLUX=flux, INFO=info, LOG=log

; Keyword check
IF NOT KEYWORD_SET(files) THEN BEGIN
  path  = '/Volumes/hosts_results/nomic/l1_fits/'
  files = [path + '2014-02-12/UT2014-02-12_ID14_SCI_eta_crv_DIT-85ms_11um_PHOT2.fits', path + '2014-02-15/UT2014-02-15_ID61_SCI_alf_Lyr_DIT-42ms_11um_PHOT2.fits', + $
           path + '2014-03-17/UT2014-03-17_ID05_SCI_alf_Lyr_DIT-27ms_11um_PHOT2.fits', path + '2014-11-10/UT2014-11-10_ID14_SCI_eps_eri_DIT-15ms_11um_PHOT2.fits', + $
           path + '2015-02-08/UT2015-02-08_ID13_SCI_bet_Leo_DIT-60ms_11um_PHOT1.fits', path + '2015-02-08/UT2015-02-08_ID13_SCI_bet_Leo_DIT-60ms_11um_PHOT2.fits']
ENDIF

; Get astronomical parameters
GET_PRM, prm

; Prepare scaling of stellar flux to the aperture radius used
n_pt  = 1D+3

; Number of files
n_files = N_ELEMENTS(files)

; Choose aperture radius 
aper_rad0 = 8

; Loop over the files
FOR i_f = 0, n_files-1 DO BEGIN
  ; Read file
  data   = MRDFITS(files[i_f], 1, header, /SILENT)
  header = HEADFITS(files[i_f], /SILENT)
  n_data = N_ELEMENTS(data.mjd_obs)
       
  ; Look for corresponding background file
  file_bckg = files[i_f]     
       
  ; Find closest aperture radius
  n_apr    = N_ELEMENTS(data[0].flx_tot)
  apr_all  = FXPAR(header, 'APERRAD' + STRING(0, FORMAT='(I0)'))
  FOR i_r  = 1, n_apr-1 DO apr_all = [apr_all, FXPAR(header, 'APERRAD' + STRING(i_r, FORMAT='(I0)'))]
  tmp      = MIN(ABS(apr_all-aper_rad0), i_aper)
  aper_rad = apr_all[i_aper]
  
  ; Normalization factor. Background and sensitivity computation are done assuming both shutter opened and nulling mode
  flx_norm = 1
  bck_norm = 1
  
  ; Extract data  
  time     = data.mjd_obs
  flx_tot  = data.flx_tot[i_aper,*]*flx_norm
  flx_err  = data.flx_err[i_aper,*]*SQRT(flx_norm)
  bck_flx  = data.bck_tot[i_aper,*]*flx_norm
  bck_err  = data.bck_err[i_aper,*]*SQRT(bck_norm)  
  
  ; Header information    
  int_time = FXPAR(header, 'EXPTIME')
  acq_time = FXPAR(header, 'ACQTIME') > int_time
  tgt_name = FXPAR(header, 'OBJNAME')
  date     = FXPAR(header, 'DATE_OBS')
  flag     = FXPAR(header, 'FLAG')
  lambda   = FXPAR(header, 'WAV_EFF')
  bandwidth= FXPAR(header, 'BANDWIDT')
  detmode  = FXPAR(header, 'DETMODE')
  bckmode  = FXPAR(header, 'BCK_MODE')
  no_dark  = FXPAR(header, 'NO_DARK')
  n_xpix   = FXPAR(header, 'N_XPIX')
  xcen     = FXPAR(header, 'XCEN')
  bck_irad = FXPAR(header, 'BCK_IRAD')
  bck_orad = FXPAR(header, 'BCK_ORAD')
  xcen     = FXPAR(header, 'XCEN')
  xpos     = xcen + 0.5*(1024-n_xpix)
    
  ; Get config parameters
  GET_CNF, cnf, LAMBDA=lambda, BANDWIDTH=bandwidth, INSTRUM='nomic'
  adu2e_lhft = [145,1500]
  ;det_adu2e = FXPAR(header, 'EPERADU')  ; DOUBLE CHECK!! Must be 67 with the preamp
  date_shrt =  LONG(STRMID(date, 2, 2) + STRMID(date, 5, 2) +  STRMID(date, 8, 2))
  IF date_shrt LE 161016 THEN det_adu2e = adu2e_lhft[detmode] ELSE det_adu2e = cnf.adu2e[detmode]; temporary hack. This assumes that the level-shift was used before 161016 (must be correct)

  ; Stellar flux
  IF NOT KEYWORD_SET(FLUX) THEN BEGIN
    star_flux = HOSTS_FLX(tgt_name)
    IF star_flux EQ -1 THEN BEGIN
       star      = GET_TGT(tgt_name)
       star_flux = BLACKBODY(star.temp, lambda, STANDARD=0) * !Dpi*(0.5*star.ldm*prm.m2r)^2        ; flux in Jy
       PRINT, 'Computed stellar flux : ', star_flux
    ENDIF
  ENDIF ELSE star_flux = flux
  star_ph   = star_flux*cnf.jy2ph*int_time*2.*!Dpi*((0.5*cnf.tel_diam)^2-(0.5*cnf.cen_obs)^2) ; flux in ph/read before instrument (throughput = 1)
  star_e    = star_ph                                                                         ; flux in e-/read before instrument (throughput = 1 and QE=1)
  star_adu  = star_e/det_adu2e                                                                ; flux in ADU/read before instrument (throughput = 1 and QE=1)
  
  ; Save results
  data_out = {ACQ_TIME: acq_time, DIT: int_time, FLAG: STRCOMPRESS(flag, /REMOVE_ALL), FLX: star_flux, LAMBDA: lambda}
      
  ; Scale flux to the aperture radius
  psf_pix     = aper_rad^2*!Dpi                                  ; Pixel in the photometric aperture
  theta       = aper_rad/cnf.psf_pix*cnf.psf_rad                ; radial distance from center (rad)
  r           = 2.*!Dpi/lambda*0.5*cnf.tel_diam*theta
  t           = (1+DINDGEN(n_pt))/n_pt*r
  ratio       = cnf.cen_obs/cnf.tel_diam                     ; secondary to primary ratio
  energy_frac = 1./(1.-ratio^2)*(1.-BESELJ(r,0)^2-BESELJ(r,1)^2+ratio^2*(1-BESELJ(ratio*r,0)^2-BESELJ(ratio*r,1)^2)-$
                4.*ratio*INT_TABULATED(t,BESELJ(t,1)*BESELJ(ratio*t,1)/t))
  
  ; Convert measured flux to Jy (sensitivity is defined at constructive peak)
  AVGSDV, flx_tot, flx_avg, flx_rms_tim, KAPPA=5, WEI=(1./flx_err^2)
  flx_con = 2.*2.*flx_avg       ; Star flux computed for both aperture at peak 
  jy2adu  = flx_con/star_flux   ; Flux at peak of a 1Jy star
  
  ; Compute spatial noise
  flx_rms_spa =  MEAN(bck_err)*SQRT(psf_pix)
  
  ; Account for background subtraction
  ; Negative background mode means that we have subtracted the mean of the background frames.
  ; Positive background mode means that we have subtracted the background frame by frame
  ; Hence, we multiply by SQRT(2) for negative background mode. We want here the sensitivity in nulling mode!
  IF bckmode LT 0 THEN BEGIN
    flx_rms_spa *= SQRT(2)
    flx_rms_tim *= SQRT(2)
  ENDIF
  
  ; Compute transmission
  trans = flx_con*det_adu2e/(star_ph*energy_frac)
  
  ; Compute sensitivity in mJy (per frame as stellar flux in Jy divided by SNR)  
  sen_spa = 1D+3*star_flux/(flx_con/flx_rms_spa)
  sen_tim = 1D+3*star_flux/(flx_con/flx_rms_tim) 
  sen_ron = 1D+3*cnf.ron/det_adu2e/jy2adu*SQRT(psf_pix)*SQRT(2)
    
  ; Compute mean background and date for this day
  bck_jy     = bck_flx/jy2adu
  idx_ok     = WHERE(bck_jy NE 0., n_ok)
  bck_date   = MEAN(bck_jy[idx_ok])
    
  ; Compute sensitivity in mJy/10min, including overhead
  t_tot      = 10.*60.
  eff        = 1.
  fac        = (acq_time/(eff*t_tot))^0.5  
  
  ; Print info to screen
  IF KEYWORD_SET(INFO) THEN BEGIN
    PRINT, ' '
    PRINT, ' '
    PRINT, date
    PRINT, 'Target name: ', tgt_name
    PRINT, 'Observation parameters' 
    PRINT, ' - Integration time [s]               : ', STRING(int_time, FORMAT='(F6.4)')
    PRINT, ' - Acquisitionn time [s]              : ', STRING(acq_time, FORMAT='(F6.4)')
    PRINT, ' - Wavelength [um]                    : ', STRING(lambda*1D+6, FORMAT='(F6.3)')
    PRINT, ' - Bandwidth [um]                     : ', STRING(bandwidth*1D+6, FORMAT='(F6.4)')
    PRINT, 'Detector parameters'
    PRINT, ' - Detector mode                      : ', STRING(detmode, FORMAT='(I0)')
    PRINT, ' - Electrons per ADU                  : ', STRING(det_adu2e, FORMAT='(I0)')
    PRINT, ' - Quantum efficiency                 : ', STRING(cnf.qe, FORMAT='(F4.2)')
    PRINT, ' - Estimated read noise               : ', STRING(cnf.ron, FORMAT='(I0)')
    PRINT, ' - Frame size [pix]                   : ', STRING(n_xpix, FORMAT='(I0)')
    PRINT, ' - Beam position [pix]                : ', STRING(xpos, FORMAT='(I0)')
    ;PRINT, ' - Detector gain                     : ', STRING(detmode, FORMAT='(I0)'), '       (0: high gain, 1: low gain)'
    PRINT, 'Reduction parameters                '
    PRINT, ' - Background mode                    : ', STRING(bckmode, FORMAT='(I0)')
    PRINT, ' - Aperture radius [pix]              : ', STRING(aper_rad, FORMAT='(I0)')
    PRINT, ' - Background inner radius [pix]      : ', STRING(bck_irad, FORMAT='(I0)')
    PRINT, ' - Background outer radius [pix]      : ', STRING(bck_orad, FORMAT='(I0)')
    PRINT, ' - Fraction energy in phot. aperture  : ', STRING(energy_frac, FORMAT='(F6.4)')
    PRINT, 'Estimated target flux (constructive peak, optical throughput=100%, and QE=1)'
    PRINT, ' - Jansky                             : ', STRING(star_flux, FORMAT='(F6.3)')
    PRINT, ' - Photons/read                       : ', STRING(star_ph, FORMAT='(G9.4)')
    PRINT, ' - Electrons/read                     : ', STRING(star_e, FORMAT='(G9.4)')
    PRINT, ' - ADU/read                           : ', STRING(star_adu, FORMAT='(I0)')
    PRINT, 'Measured flux (constructive peak) '
    PRINT, ' - In photometric aperture [ADU/read] : ', STRING(flx_con, FORMAT='(I0)')
    PRINT, ' - Total [ADU/read]                   : ', STRING(flx_con/energy_frac, FORMAT='(I0)')
    PRINT, 'Measured noise (constructive peak) '
    PRINT, ' - In photometric aperture [ADU/read] : ', STRING(flx_rms_tim, FORMAT='(I0)') + ', ' + STRING(flx_rms_spa, FORMAT='(I0)'), ' (temporal and spatial)'
    PRINT, ' - SNR                                : ', STRING(flx_con/flx_rms_tim, FORMAT='(I0)') + ', ' + STRING(flx_con/flx_rms_spa, FORMAT='(I0)'), ' (temporal and spatial)'
    PRINT, 'Flux at peak for 1Jy star [ADU/read]  : ', STRING(jy2adu, FORMAT='(I0)')
    PRINT, 'Flux at peak for 1Jy star [ADU/s]     : ', STRING(jy2adu/int_time, FORMAT='(I0)')
    PRINT, 'Estimated throughput                 '
    PRINT, ' - Optical throughput [%]             : ', STRING(100.*trans/cnf.qe, FORMAT='(F6.3)')
    PRINT, ' - Overall throughput [%]             : ', STRING(100.*trans, FORMAT='(F6.3)')
    ;IF no_dark EQ 0 THEN PRINT, 'Background [Jy]                      : ', STRING(bck_date, FORMAT='(F6.4)')  , '       (only valid if data reduced with darks)'
    PRINT, 'Expected sensitivity (RON only)       : '
    PRINT, ' - spatial                            : ', STRING(sen_ron, FORMAT='(F7.3)') + ', ' +  STRING(fac*sen_ron, FORMAT='(F7.3)'), ' (mJy/read and mJy/10min)'
    ;PRINT, ' - temporal                           : ', STRING(1D+3*sen_ron, FORMAT='(F7.3)') + ', ' +  STRING(1D+3*sen_ron10, FORMAT='(F7.3)'), ' (mJy/read and mJy/10min)'    
    PRINT, 'Measured sensitivity on target        : '
    PRINT, ' - spatial                            : ', STRING(sen_spa, FORMAT='(F7.3)') + ', ' +  STRING(fac*sen_spa, FORMAT='(F7.3)'), ' (mJy/read and mJy/10min)'
    PRINT, ' - temporal                           : ', STRING(sen_tim, FORMAT='(F7.3)') + ', ' +  STRING(fac*sen_tim, FORMAT='(F7.3)'), ' (mJy/read and mJy/10min)'
  ENDIF
  
  ; Print info to log
  IF KEYWORD_SET(LOG) THEN BEGIN
    ; Log directory will be the same as the one with the L1 files
    DECLARE_PATH,  pth, INSTRUM='NOMIC', PATH_FILE='path.cfg'
    log_dir  = pth.result_path + pth.sep + 'diagnostic' + pth.sep + 'sensitivity' + pth.sep + date + pth.sep
    IF NOT FILE_TEST(log_dir) THEN FILE_MKDIR, log_dir
    ; Create log file
    log_file =  log_dir + date + '_' + tgt_name + '.txt'
    OPENW, lun, log_file, /GET_LUN, WIDTH=800, /APPEND
    PRINTF, lun, ' '
    PRINTF, lun, ' '
    PRINTF, lun, date
    PRINTF, lun, 'Target name: ', tgt_name
    PRINTF, lun, 'Observation parameters'
    PRINTF, lun, ' - Integration time [s]               : ', STRING(int_time, FORMAT='(F6.4)')
    PRINTF, lun, ' - Acquisitionn time [s]              : ', STRING(acq_time, FORMAT='(F6.4)')
    PRINTF, lun, ' - Wavelength [um]                    : ', STRING(lambda*1D+6, FORMAT='(F6.3)')
    PRINTF, lun, ' - Bandwidth [um]                     : ', STRING(bandwidth*1D+6, FORMAT='(F6.4)')
    PRINTF, lun, 'Detector parameters'
    PRINTF, lun, ' - Electrons per ADU                  : ', STRING(det_adu2e, FORMAT='(I0)')
    PRINTF, lun, ' - Quantum efficiency                 : ', STRING(cnf.qe, FORMAT='(F4.2)')
    PRINTF, lun, ' - Estimated read noise               : ', STRING(cnf.ron, FORMAT='(I0)')
    PRINTF, lun, ' - Frame size [pix]                   : ', STRING(n_xpix, FORMAT='(I0)')
    PRINTF, lun, ' - Beam position [pix]                : ', STRING(xpos, FORMAT='(I0)')
    ;PRINTF, lun, ' - Detector gain                     : ', STRING(detmode, FORMAT='(I0)'), '       (0: high gain, 1: low gain)'
    PRINTF, lun, 'Reduction parameters                '
    PRINTF, lun, ' - Background mode                    : ', STRING(bckmode, FORMAT='(I0)')
    PRINTF, lun, ' - Aperture radius [pix]              : ', STRING(aper_rad, FORMAT='(I0)')
    PRINTF, lun, ' - Background inner radius [pix]      : ', STRING(bck_irad, FORMAT='(I0)')
    PRINTF, lun, ' - Background outer radius [pix]      : ', STRING(bck_orad, FORMAT='(I0)')
    PRINTF, lun, ' - Fraction energy in phot. aperture  : ', STRING(energy_frac, FORMAT='(F6.4)')
    PRINTF, lun, 'Estimated target flux (constructive peak, optical throughput=100%, and QE=1)'
    PRINTF, lun, ' - Jansky                             : ', STRING(star_flux, FORMAT='(F6.3)')
    PRINTF, lun, ' - Photons/read                       : ', STRING(star_ph, FORMAT='(G9.4)')
    PRINTF, lun, ' - Electrons/read                     : ', STRING(star_e, FORMAT='(G9.4)')
    PRINTF, lun, ' - ADU/read                           : ', STRING(star_adu, FORMAT='(I0)')
    PRINTF, lun, 'Measured flux (constructive peak) '
    PRINTF, lun, ' - In photometric aperture [ADU/read] : ', STRING(flx_con, FORMAT='(I0)')
    PRINTF, lun, ' - Total [ADU/read]                   : ', STRING(flx_con/energy_frac, FORMAT='(I0)')
    PRINTF, lun, 'Measured noise (constructive peak) '
    PRINTF, lun, ' - In photometric aperture [ADU/read] : ', STRING(flx_rms_tim, FORMAT='(I0)') + ', ' + STRING(flx_rms_spa, FORMAT='(I0)'), ' (temporal and spatial)'
    PRINTF, lun, ' - SNR                                : ', STRING(flx_con/flx_rms_tim, FORMAT='(I0)') + ', ' + STRING(flx_con/flx_rms_spa, FORMAT='(I0)'), ' (temporal and spatial)'
    PRINTF, lun, 'Flux at peak for 1Jy star [ADU/read]  : ', STRING(jy2adu, FORMAT='(I0)')
    PRINTF, lun, 'Flux at peak for 1Jy star [ADU/s]     : ', STRING(jy2adu/int_time, FORMAT='(I0)')
    PRINTF, lun, 'Estimated throughput                 '
    PRINTF, lun, ' - Optical throughput [%]             : ', STRING(100.*trans/cnf.qe, FORMAT='(F6.3)')
    PRINTF, lun, ' - Overall throughput [%]             : ', STRING(100.*trans, FORMAT='(F6.3)')
    ;IF no_dark EQ 0 THEN PRINTF, lun, 'Background [Jy]                      : ', STRING(bck_date, FORMAT='(F6.4)')  , '       (only valid if data reduced with darks)'
    PRINTF, lun, 'Expected sensitivity (RON only)       : '
    PRINTF, lun, ' - spatial                            : ', STRING(sen_ron, FORMAT='(F7.3)') + ', ' +  STRING(fac*sen_ron, FORMAT='(F7.3)'), ' (mJy/read and mJy/10min)'
    ;PRINTF, lun, ' - temporal                           : ', STRING(1D+3*sen_ron, FORMAT='(F7.3)') + ', ' +  STRING(1D+3*sen_ron10, FORMAT='(F7.3)'), ' (mJy/read and mJy/10min)'
    PRINTF, lun, 'Measured sensitivity on target        : '
    PRINTF, lun, ' - spatial                            : ', STRING(sen_spa, FORMAT='(F7.3)') + ', ' +  STRING(fac*sen_spa, FORMAT='(F7.3)'), ' (mJy/read and mJy/10min)'
    PRINTF, lun, ' - temporal                           : ', STRING(sen_tim, FORMAT='(F7.3)') + ', ' +  STRING(fac*sen_tim, FORMAT='(F7.3)'), ' (mJy/read and mJy/10min)'
    CLOSE, lun
    FREE_LUN, lun
  ENDIF
ENDFOR

RETURN, [trans/cnf.qe, trans, fac*sen_tim, fac*sen_ron, jy2adu]
END
; ----------- TEST ROUTINE ------------
PRO TEST_SENSITIVITY
  ; Keyword check
  path = '/Volumes/hosts_results/nomic/l1_fits/'
  ;path = '/Volumes/LaCie/nodrs/results/nomic/l1_fits/'
  info = 1
  
;  ; Feb 2014 Run
;    PRINT, '====================================================================== '
;  files = path + ['2014-02-12/UT2014-02-12_ID20_SCI_eta_crv_DIT-85ms_11um_PHOT1.fits', '2014-02-12/UT2014-02-12_ID20_SCI_eta_crv_DIT-85ms_11um_PHOT2.fits']
;  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=1.76*1.2, INFO=info)
;  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=1.76*1.2, INFO=info)
;  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F6.4)')
;  PRINT, '====================================================================== '
;  
;  ; Mar 2014 Run
;    PRINT, '====================================================================== '
;  files = path + ['2014-03-17/UT2014-03-17_ID05_SCI_alf_Lyr_DIT-27ms_11um_PHOT1.fits', '2014-03-17/UT2014-03-17_ID05_SCI_alf_Lyr_DIT-27ms_11um_PHOT2.fits']
;  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=38.55, INFO=info)
;  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=38.55, INFO=info)
;  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F6.4)')
;  PRINT, '====================================================================== '
;  
;  ; Nov 2014 Run 1
;    PRINT, '====================================================================== '
;  files = path + ['2014-11-09/UT2014-11-09_ID15_SCI_eps_eri_DIT-15ms_11um_PHOT1.fits', '2014-11-09/UT2014-11-09_ID15_SCI_eps_eri_DIT-15ms_11um_PHOT2.fits']
;  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=7.39, INFO=info)
;  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=7.39, INFO=info)
;  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F6.4)')
;  PRINT, '====================================================================== '
;  
;  ; Nov 2014 Run 2
;    PRINT, '====================================================================== '
;  files = path + ['2014-11-10/UT2014-11-10_ID14_SCI_eps_eri_DIT-15ms_11um_PHOT1.fits', '2014-11-10/UT2014-11-10_ID14_SCI_eps_eri_DIT-15ms_11um_PHOT2.fits']
;  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=7.39, INFO=info)
;  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=7.39, INFO=info)
;  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F6.4)')
;  PRINT, '====================================================================== '
;  
;  ; Jan 2015 Run (30ms)
;    PRINT, '====================================================================== '
;  files = path + ['2015-01-10/UT2015-01-10_ID02_SCI_betaLeo_DIT-30ms_11um_PHOT1.fits', '2015-01-10/UT2015-01-10_ID02_SCI_betaLeo_DIT-30ms_11um_PHOT2.fits']
;  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=6.85, INFO=info)
;  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=6.85, INFO=info)
;  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F6.4)')
;  PRINT, '====================================================================== '
;  
;  ; Jan 2015 Run (60ms)
;    PRINT, '====================================================================== '
;  files = path + ['2015-01-10/UT2015-01-10_ID03_SCI_betaLeo_DIT-60ms_11um_PHOT1.fits', '2015-01-10/UT2015-01-10_ID03_SCI_betaLeo_DIT-60ms_11um_PHOT2.fits']
;  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=6.85, INFO=info)
;  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=6.85, INFO=info)
;  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F6.4)')
;  PRINT, '====================================================================== '
  
;  ; Feb 2015 Run
  PRINT, '====================================================================== '
  files = path + ['2015-02-08_APR/UT2015-02-08_ID013_SCI_bet_Leo_DIT-60ms_11um_PHOT1.fits', '2015-02-08_APR/UT2015-02-08_ID013_SCI_bet_Leo_DIT-60ms_11um_PHOT2.fits']
  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=6.85, INFO=info)
  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=6.85, INFO=info)
  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F6.4)')
  PRINT, '====================================================================== '
;  
;  ; Feb 2016 Run
;  ; eta Boo is w3=1.447. This is 8.35 Jy in W3. For NOMIC, using the same scaling as Gamma Ser (other G0iV), this gives 9.59
;  PRINT, '====================================================================== '
;  files = path + ['2016-02-17/UT2016-02-17_ID014_SCI_eta_boo_DIT-35ms_11um_PHOT1.fits', '2016-02-17/UT2016-02-17_ID014_SCI_eta_boo_DIT-35ms_11um_PHOT2.fits']
;  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=9.59, INFO=info)
;  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=9.59, INFO=info)
;  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F6.4)')
;  PRINT, '====================================================================== '
;  
;  ; April 2016 Run
;  ; Vega
;  PRINT, '====================================================================== '
;  files = path + ['2016-04-18_APR/UT2016-04-18_ID007_SCI_alf_Lyr_DIT-43ms_11um_PHOT1.fits', '2016-04-18_APR/UT2016-04-18_ID007_SCI_alf_Lyr_DIT-43ms_11um_PHOT2.fits']
;  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=38.55, INFO=info)
;  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=38.55, INFO=info)
;  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F8.4)')
;  PRINT, '====================================================================== '
  
  ; Nay 2016 Run
  ; Vega
;  PRINT, '====================================================================== '
;  files = path + ['2016-06-19_APR/UT2016-06-19_ID000_SCI_Vega_DIT-1004ms_11um_PHOT1.fits', '2016-06-19_APR/UT2016-06-19_ID001_SCI_Vega_DIT-1004ms_11um_PHOT2.fits']
;  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=38.55, INFO=info)
;  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=38.55, INFO=info)
;  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
;  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F8.4)')
;  PRINT, '====================================================================== '
  
  ; Oct 2016 Run
  ; Alf Cep
  PRINT, '====================================================================== '
  files = path + ['2016-10-16_APR/UT2016-10-16_ID009_SCI_alf_Cep_DIT-53ms_11um_PHOT1.fits', '2016-10-16_APR/UT2016-10-16_ID009_SCI_alf_Cep_DIT-53ms_11um_PHOT2.fits']
  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=7.04, INFO=info)
  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=7.04, INFO=info)
  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F8.4)')
  PRINT, '====================================================================== '
  
  ; May 2018 Run
  ; Tet Boo
  PRINT, '====================================================================== '
  files = path + ['2018-05-23_APR/UT2018-05-23_ID015_SCI_tet_Boo_DIT-45ms_11um_PHOT1.fits', '2018-05-23_APR/UT2018-05-23_ID015_SCI_tet_Boo_DIT-45ms_11um_PHOT2.fits']
  data_p1 = LBTI_SENSITIVITY(FILES=files[0], FLUX=3.15, INFO=info)
  data_p2 = LBTI_SENSITIVITY(FILES=files[1], FLUX=3.15, INFO=info)
  PRINT, 'Mean optical throughput [%]   : ', STRING(100*(data_p1[0]+data_p2[0])*0.5, FORMAT='(F6.3)')
  PRINT, 'Mean overall throughput [%]   : ', STRING(100*(data_p1[1]+data_p2[1])*0.5, FORMAT='(F6.3)')
  PRINT, 'Mean sensitivity [mJy/10min]  : ', STRING((data_p1[2]+data_p2[2])*0.5, FORMAT='(F8.4)')
  PRINT, '====================================================================== '
  
END

; ----------- PIPELINE PAPER ROUTINE (NOW OBSOLOTE BUT KEEP AS IS) ------------
PRO TEST_SENSITIVITY2

; Date to process
date = ['2013-12-31','2014-02-12','2014-02-15','2014-03-17','2014-11-09','2014-11-10','2015-01-10','2015-02-03','2015-02-04','2015-02-08','2016-10-16']
;date = ['2015-02-08']

; Aperture radius and others
GET_CNF, cnf, INSTRUM='NOMIC', LAMBDA=lambda, BANDWIDTH=bandwidth
aper   = 8
i_aper = 2
thru   = 0.033
ron    = cnf.ron
lambda = 1.1D-5

; Compute energy fraction
n_pt        = 1D4
theta       = aper/cnf.psf_pix*cnf.psf_rad                ; radial distance from center (rad)
r           = 2.*!Dpi/lambda*0.5*cnf.tel_diam*theta
t           = (1+DINDGEN(n_pt))/n_pt*r
ratio       = cnf.cen_obs/cnf.tel_diam                     ; secondary to primary ratio
energy_frac = 1./(1.-ratio^2)*(1.-BESELJ(r,0)^2-BESELJ(r,1)^2+ratio^2*(1-BESELJ(ratio*r,0)^2-BESELJ(ratio*r,1)^2)-$
              4.*ratio*INT_TABULATED(t,BESELJ(t,1)*BESELJ(ratio*t,1)/t))

; Data path
path = '/Volumes/nodrs/nodrs/results/nomic/l1_fits/'
info = 0

; Init plot
; Plot results
charthick = 4.0
charsize  = 1.3
result_path = '/Volumes/nodrs/nodrs/results/nomic/diagnostic/sensitivity/'
IF NOT FILE_TEST(result_path) THEN FILE_MKDIR, result_path

plotname =  result_path + 'sensitivity_all.eps'
PREP_PS, /BOLD
DEVICE, FILEN=plotname, /ENCAPS, /COLOR,  XSIZE=17.8, YSIZE=13.2, /TIMES
xrange     = [5.,100]
yrange     = [0.1,10.]
PLOT, [0], [0], XTITLE='Detector integration time [ms]', YTITLE='Photometric sensitivity [mJy/10min]', TITLE=title_day, XSTYLE=1, YSTYLE=1, /YLOG, /XLOG, XRANGE=xrange, YRANGE=yrange;, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize
LOADCT, 12, /SILENT
; Loop over the dates
n_date = N_ELEMENTS(date)
n      = 0 ; number of iterations
thru   = 0 ; will contain the mean throughput
FOR id = 0, n_date-1 DO BEGIN
  ; Search PHOT files
  phot1_files = FILE_SEARCH(path + date[id] + '/', '*PHOT1.fits', COUNT=n_phot1)
  ; Loop over the files
  FOR ip=0, n_phot1-1 DO BEGIN
    ; Find PHOT2 files
    phot1_file = phot1_files[ip]
    phot2_file = STRMID(phot1_file, 0, STRPOS(phot1_file, 'PHOT1')) + 'PHOT2.fits'
    ; Process file
    data_p1 = LBTI_SENSITIVITY(DATA_OUT=data_out1, FILES=phot1_file, FLUX=ff, INFO=info)
    data_p2 = LBTI_SENSITIVITY(DATA_OUT=data_out2, FILES=phot2_file, FLUX=ff, INFO=info)
    ; Background file
;    bckg_file = STRMID(phot1_file, 0, STRPOS(phot1_file, 'PHOT1')) + 'BCKG.fits'
;    data_bckg = MRDFITS(bckg_file, 1, /SILENT)
;    bck_tot_phot = REFORM(data_bckg.flx_tot[i_aper]) & bck_err_phot = REFORM(data_bckg.flx_err[i_aper])     ; used in the statistical reduction technique
;    ; Remove obvious outliers
;    AVGSDV, bck_tot_phot, bck_avg, bck_rms, bck_rms_mean, KAPPA=5 & idx_ok = WHERE(ABS(bck_tot_phot-bck_avg) LE 5*bck_rms)
;    bck_tot_phot = bck_tot_phot[idx_ok] & bck_err_phot = bck_err_phot[idx_ok]
;    ; Compute mean and rms values
;    AVGSDV, bck_tot_phot, bck_avg, bck_rms, bck_rms_mean, KAPPA=5, WEI=(1/bck_err_phot^2)
;    bck_rms1 = bck_rms*cnf.adu2e[0]/thru/energy_frac/data_out1.dit/(cnf.jy2ph*2.*!Dpi*((0.5*cnf.tel_diam)^2-(0.5*cnf.cen_obs)^2))/SQRT(60.*10/data_out1.acq_time)      ; Convert background noise to ph/pix
;    USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
    ;OPLOT, [1D3*data_out1.dit], 1D3*[bck_rms1], PSYM=8, COLOR=0
    data_bckg = MRDFITS(phot1_file, 1, /SILENT)
    bck_tot_phot = REFORM(data_bckg.flx_tot[i_aper]) & bck_err_phot = REFORM(data_bckg.flx_err[i_aper])     ; used in the statistical reduction technique
    ; Remove obvious outliers
    AVGSDV, bck_tot_phot, bck_avg, bck_rms, bck_rms_mean, KAPPA=5 & idx_ok = WHERE(ABS(bck_tot_phot-bck_avg) LE 5*bck_rms)
    bck_tot_phot = bck_tot_phot[idx_ok] & bck_err_phot = bck_err_phot[idx_ok]
    ; Compute mean and rms values
    AVGSDV, bck_tot_phot, bck_avg, bck_rms, bck_rms_mean, KAPPA=5, WEI=(1/bck_err_phot^2)
    bck_rms1 = bck_rms*cnf.adu2e[0]/thru/energy_frac/data_out1.dit/(cnf.jy2ph*2.*!Dpi*((0.5*cnf.tel_diam)^2-(0.5*cnf.cen_obs)^2))/SQRT(60.*10/data_out1.acq_time)      ; Convert background noise to ph/pix
    USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
    data_bckg = MRDFITS(phot2_file, 1, /SILENT)
    bck_tot_phot = REFORM(data_bckg.flx_tot[i_aper]) & bck_err_phot = REFORM(data_bckg.flx_err[i_aper])     ; used in the statistical reduction technique
    ; Remove obvious outliers
    AVGSDV, bck_tot_phot, bck_avg, bck_rms, bck_rms_mean, KAPPA=5 & idx_ok = WHERE(ABS(bck_tot_phot-bck_avg) LE 5*bck_rms)
    bck_tot_phot = bck_tot_phot[idx_ok] & bck_err_phot = bck_err_phot[idx_ok]
    ; Compute mean and rms values
    AVGSDV, bck_tot_phot, bck_avg, bck_rms, bck_rms_mean, KAPPA=5, WEI=(1/bck_err_phot^2)
    bck_rms2 = bck_rms*cnf.adu2e[0]/thru/energy_frac/data_out1.dit/(cnf.jy2ph*2.*!Dpi*((0.5*cnf.tel_diam)^2-(0.5*cnf.cen_obs)^2))/SQRT(60.*10/data_out1.acq_time)      ; Convert background noise to ph/pix
    USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
    ;OPLOT, [1D3*data_out1.dit], 1D3*0.5*[bck_rms1+bck_rms2], PSYM=8, COLOR=0
    ; Keep only non super-monster stars
    sen10  = 0.5*1D3*(data_p1[2]+data_p2[2])
    ron10  = 0.5*1D3*(data_p1[3]+data_p2[3])
    idx_ok = WHERE(data_out2.flx LT 200 AND sen10 LT 6.3 AND sen10 GT 0.2 AND data_out2.lambda GT 1D-5, n_ok)
    IF n_ok GT 0 THEN BEGIN
      ; Plot results
      IF data_out1.flag EQ 'CAL' THEN BEGIN
        USERSYM, [1,1,-1,-1,1]*0.7, [1,-1,-1,1,1]*0.8, THICK=1.5, /FILL
        OPLOT, [1D3*data_out1.dit], [sen10], PSYM=8, COLOR=100
      ENDIF ELSE BEGIN
        USERSYM, [1,0,-1,0,1], [0,1,0,-1,0], THICK=1.5, /FILL
        OPLOT, [1D3*data_out1.dit], [sen10], PSYM=8, COLOR=200
      ENDELSE
      ;USERSYM, [1,0,-1,0,1], [0,1,0,-1,0], THICK=1.5, /FILL
      ;OPLOT, [1D3*data_out1.dit], [ron10], PSYM=8, COLOR=0
    ENDIF
    ;ERRPLOT, caltime_pt*fac, (null_avg_pt[idx_cal_pt]+null_err_pt[idx_cal_pt])*1D2, (null_avg_pt[idx_cal_pt]-null_err_pt[idx_cal_pt])*1D2, COLOR=100
    n =+ 1
    thru =+ 0.5*(data_p1[1]+data_p2[1])
  ENDFOR
ENDFOR
LOADCT, 0, /SILENT
thru = thru/n
PRINT, MEAN(thru)
  
; Overplot expected sensitivity for readout noise only;
; Sensitivity in ph/s over DIT is SQRT(ron^2+bckg*T)/T where T is the DIT
int_time = 1D-3*(MIN(xrange) + (MAX(xrange)-MIN(xrange))*DINDGEN(100)/99.)
bckg_m   = 0                                        ; Measured background in adu/s/pix 
bckg     = bckg_m*cnf.adu2e[0]/cnf.qe*int_time      ; Convert background noise to ph/pix
bckg_ph  = SQRT(bckg)                               ; Background noise in ph/pix
ron_ph   = ron/cnf.qe                          ; RON in ph/pix
sen_ph   = SQRT(ron_ph^2+bckg_ph^2)/int_time        ; Sensitivity in ph/s
sen_jy   = SQRT(4)*SQRT(!dpi*aper^2)*sen_ph*cnf.qe/thru/energy_frac/(cnf.jy2ph*2.*!Dpi*((0.5*cnf.tel_diam)^2-(0.5*cnf.cen_obs)^2)) ;SQRT(4) to account for double background subtraction 
sen_jy10 = sen_jy/SQRT(60.*10/int_time)
OPLOT, 1D3*int_time, 1D3*sen_jy10, LINESTYLE=0, COLOR=0

; Overplot expected sensitivity for readout noise and background noise;
bckg_m   = 91733.3                                  ; Measured background in adu/s/pix 
PRINT, 'Estimated background flux per pixel', bckg_m*cnf.adu2e[0]/cnf.qe/thru/(cnf.jy2ph*2.*!Dpi*((0.5*cnf.tel_diam)^2-(0.5*cnf.cen_obs)^2))
bckg     = bckg_m*cnf.adu2e[0]/cnf.qe*int_time      ; Convert background noise to ph/pix
bckg_ph  = SQRT(bckg)                               ; Background noise in ph/pix
ron_ph   = ron/cnf.qe                          ; RON in ph/pix
sen_ph   = SQRT(ron_ph^2+bckg_ph^2)/int_time        ; Sensitivity in ph/s
sen_jy   = SQRT(4)*SQRT(!dpi*aper^2)*sen_ph*cnf.qe/thru/energy_frac/(cnf.jy2ph*2.*!Dpi*((0.5*cnf.tel_diam)^2-(0.5*cnf.cen_obs)^2)) ;SQRT(4) to account for double background subtraction 
sen_jy10_max = sen_jy/SQRT(60.*10/int_time)

bckg_m   = 91733.3/2                                  ; Measured background in adu/s/pix
PRINT, 'Estimated background flux per pixel', bckg_m*cnf.adu2e[0]/cnf.qe/thru/(cnf.jy2ph*2.*!Dpi*((0.5*cnf.tel_diam)^2-(0.5*cnf.cen_obs)^2))
bckg     = bckg_m*cnf.adu2e[0]/cnf.qe*int_time      ; Convert background noise to ph/pix
bckg_ph  = SQRT(bckg)                               ; Background noise in ph/pix
ron_ph   = ron/cnf.qe                          ; RON in ph/pix
sen_ph   = SQRT(ron_ph^2+bckg_ph^2)/int_time        ; Sensitivity in ph/s
sen_jy   = SQRT(4)*SQRT(!dpi*aper^2)*sen_ph*cnf.qe/thru/energy_frac/(cnf.jy2ph*2.*!Dpi*((0.5*cnf.tel_diam)^2-(0.5*cnf.cen_obs)^2)) ;SQRT(4) to account for double background subtraction
sen_jy10_min = sen_jy/SQRT(60.*10/int_time)
POLYFILL, 1D3*[int_time,REVERSE(int_time)],1D3*[sen_jy10_min,REVERSE(sen_jy10_max)], COLOR=180

bckg_m   = (91733.3/2 + 91733.3)/2                  ; Measured background in adu/s/pix
PRINT, 'Estimated background flux per pixel', bckg_m*cnf.adu2e[0]/cnf.qe/thru/(cnf.jy2ph*2.*!Dpi*((0.5*cnf.tel_diam)^2-(0.5*cnf.cen_obs)^2))
bckg     = bckg_m*cnf.adu2e[0]/cnf.qe*int_time      ; Convert background noise to ph/pix
bckg_ph  = SQRT(bckg)                               ; Background noise in ph/pix
ron_ph   = ron/cnf.qe                          ; RON in ph/pix
sen_ph   = SQRT(ron_ph^2+bckg_ph^2)/int_time        ; Sensitivity in ph/s
sen_jy   = SQRT(4)*SQRT(!dpi*aper^2)*sen_ph*cnf.qe/thru/energy_frac/(cnf.jy2ph*2.*!Dpi*((0.5*cnf.tel_diam)^2-(0.5*cnf.cen_obs)^2)) ;SQRT(4) to account for double background subtraction
sen_jy10_mean = sen_jy/SQRT(60.*10/int_time)

; Find DIT where background becomes larger than RON
junk = MIN(ABS(SQRT(sen_jy10_mean^2-sen_jy10^2)-sen_jy10), idx)
OPLOT,  1D3*[int_time[idx],int_time[idx]], yrange, LINESTYLE=1, THICK=3
PRINT, 1D3*[int_time[idx],int_time[idx]]

; Close plot
DEVICE, /CLOSE & END_PS
END

; --- New routine to reduce several nights and plot the results ---
PRO TEST_SENSITIVITY3

  ; Date to process
  date = ['2013-12-31','2014-02-12','2014-02-15','2014-03-17','2014-11-09','2014-11-10','2015-01-10','2015-02-03','2015-02-04','2015-02-08','2016-10-16']
  ;date = ['2015-02-08']
  n_date = N_ELEMENTS(date)
  
  ; Data path
  root        = '/Volumes/nodrs/nodrs/results/nomic/'
  data_path   = root + 'l1_fits/'
  result_path = root + 'diagnostic/sensitivity/'
    
  ; Init plot parameters
  charthick = 4.0
  charsize  = 1.3
  ps        = 4    ; PSYM

  IF NOT FILE_TEST(result_path) THEN FILE_MKDIR, result_path
  
  ; Prepare arrays
  jy2adu  = FLTARR(n_date) & jy2adu_err  = jy2adu
  thruput = jy2adu         & thruput_err = jy2adu
  sen     = jy2adu         & sen_err     = jy2adu
   
  ; Loop over the date
  FOR id = 0, n_date-1 DO BEGIN
    
    ; Search PHOT1 files
    phot1_files = FILE_SEARCH(data_path + date[id] + '/', '*PHOT1.fits', COUNT=n_phot1) 
    
    ; Find number of independent ones (one per pointing)
    ob_id   = STRMID(phot1_files, STRPOS(phot1_files[0], '_ID')+3, 3)
    IF STRMATCH(ob_id[0], '*_') EQ 1 THEN ob_id = STRMID(phot1_files, STRPOS(phot1_files[0], '_ID')+3, 2) ; backward compatibility
    objname = STRMID(phot1_files, STRPOS(phot1_files[0], '_ID')+10, (STRPOS(phot1_files[0], '_DIT')-(STRPOS(phot1_files[0], '_ID')+10)))
    objname = objname[SORT(ob_id)]
    idx_ind = WHERE(objname NE SHIFT(objname, 1), n_ind) 
    IF n_ind LE 0 THEN idx_ind = [0] & n_ind = n_ind > 1
    
    ; Prepare arrays and loop
    star_flx  = FLTARR(n_ind)
    jy2adu_d  = FLTARR(n_ind,2)
    thruput_d = jy2adu_d
    sen_d     = jy2adu_d
    FOR i = 0, n_ind-1 DO BEGIN
      ; Photometry files
      phot1_file = phot1_files[idx_ind[i]]
      phot2_file = STRMID(phot1_file, 0, STRPOS(phot1_file, 'PHOT1')) + 'PHOT2.fits'
      ; Process file
      data_p1 = LBTI_SENSITIVITY(DATA_OUT=data_out1, FILES=phot1_file, FLUX=ff, INFO=info, /LOG)
      data_p2 = LBTI_SENSITIVITY(DATA_OUT=data_out2, FILES=phot2_file, FLUX=ff, INFO=info, /LOG)
      ; Store results
      star_flx[i]    = data_out1.flx
      jy2adu_d[i,0]  = data_p1[4] & jy2adu_d[i,1]  = data_p2[4]
      thruput_d[i,0] = data_p1[1] & thruput_d[i,1] = data_p2[1]
      sen_d[i,0]     = data_p1[2] & sen_d[i,1]     = data_p2[2]
    ENDFOR
    
    ; Compute mean values
    idx = WHERE(star_flx LT 50 AND sen_d LT 4)
    AVGSDV, 0.5*(jy2adu_d[idx,0]+jy2adu_d[idx,1]), avg, rms, KAPPA=3
    jy2adu[id] = avg & jy2adu_err[id] = rms
    AVGSDV, 0.5*(thruput_d[idx,0]+thruput_d[idx,1]), avg, rms, KAPPA=3
    thruput[id] = avg & thruput_err[id] = rms
    AVGSDV, 0.5*(sen_d[idx,0]+sen_d[idx,1]), avg, rms, KAPPA=3
    sen[id] = avg & sen_err[id] = rms
    
    ; Plot results
    PREP_PS, /BOLD
    plot_path =  result_path + date[id] + '/'
    IF NOT FILE_TEST(plot_path) THEN FILE_MKDIR, plot_path
    DEVICE, FILEN=plot_path + 'jy2adu.eps', /ENCAPS, /COLOR,  XSIZE=17.8, YSIZE=13.2, /TIMES
    xrange     = [5.,100]
    yrange     = [0.1,10.]
    PLOT, [0], [0], XTITLE='Stellar flux [Jy]', YTITLE='Jy to ADU', TITLE=title_day, XSTYLE=1, YSTYLE=1, XRANGE=[0,1.5*MAX(star_flx)], YRANGE=[0,1.5*MAX(jy2adu_d)];, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize
    LOADCT, 13, /SILENT
    OPLOT, star_flx, jy2adu_d[*,0], PSYM=ps, COLOR=90 
    OPLOT, star_flx, jy2adu_d[*,1], PSYM=ps, COLOR=255
    ; Close plot
    LOADCT, 0, /SILENT
    DEVICE, /CLOSE 
    DEVICE, FILEN=plot_path + 'thruput.eps', /ENCAPS, /COLOR,  XSIZE=17.8, YSIZE=13.2, /TIMES
    xrange     = [5.,100]
    yrange     = [0.1,10.]
    PLOT, [0], [0], XTITLE='Stellar flux [Jy]', YTITLE='Throughput [%]', TITLE=title_day, XSTYLE=1, YSTYLE=1, XRANGE=[0,1.5*MAX(star_flx)], YRANGE=[0,1.5*MAX(thruput_d)];, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize
    LOADCT, 13, /SILENT
    OPLOT, star_flx, thruput_d[*,0], PSYM=ps, COLOR=90
    OPLOT, star_flx, thruput_d[*,1], PSYM=ps, COLOR=255
    ; Close plot
    LOADCT, 0, /SILENT
    DEVICE, /CLOSE 
    DEVICE, FILEN=plot_path + 'sensitivity.eps', /ENCAPS, /COLOR,  XSIZE=17.8, YSIZE=13.2, /TIMES
    xrange     = [5.,100]
    yrange     = [0.1,10.]
    PLOT, [0], [0], XTITLE='Stellar flux [Jy]', YTITLE='Photometric sensitivity [mJy/10min]', TITLE=title_day, XSTYLE=1, YSTYLE=1, XRANGE=[0,1.5*MAX(star_flx)], YRANGE=[0,1.5*MAX(sen_d)];, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize
    LOADCT, 13, /SILENT
    OPLOT, star_flx, sen_d[*,0], PSYM=ps, COLOR=90
    OPLOT, star_flx, sen_d[*,1], PSYM=ps, COLOR=255
    ; Close plot
    LOADCT, 0, /SILENT
    DEVICE, /CLOSE & END_PS
  ENDFOR
  
  ; Plot results
  x = 1+INDGEN(n_date)
  thruput *= 1D2
  thruput_err *= 1D2
  PREP_PS, /BOLD
  plot_path =  result_path
  DEVICE, FILEN=plot_path + 'jy2adu.eps', /ENCAPS, /COLOR,  XSIZE=17.8, YSIZE=13.2, /TIMES
  xrange     = [5.,100]
  yrange     = [0.1,10.]
  PLOT, [0], [0], XTITLE='Observing run day', YTITLE='Jy to ADU', TITLE=title_day, XSTYLE=1, YSTYLE=1, XRANGE=[0,n_date+1], YRANGE=[0,1.5*MAX(jy2adu)];, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize
  LOADCT, 13, /SILENT
  OPLOT, x, jy2adu, PSYM=ps, COLOR=90
  ERRPLOT, x, jy2adu-jy2adu_err,jy2adu+jy2adu_err, COLOR=90
  FOR i=0, n_date-1 DO XYOUTS, x[i], 1.2*MAX(jy2adu), date[i], ORIENTATION=90, CHARSIZE=0.8
  ; Close plot
  LOADCT, 0, /SILENT
  DEVICE, /CLOSE
  DEVICE, FILEN=plot_path + 'thruput.eps', /ENCAPS, /COLOR,  XSIZE=17.8, YSIZE=13.2, /TIMES
  xrange     = [5.,100]
  yrange     = [0.1,10.]
  PLOT, [0], [0], XTITLE='Observing run day', YTITLE='Overall throughput [%]', TITLE=title_day, XSTYLE=1, YSTYLE=1, XRANGE=[0,n_date+1], YRANGE=[0,1.5*MAX(thruput)];, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize
  LOADCT, 13, /SILENT
  OPLOT, x, thruput, PSYM=ps, COLOR=90
  ERRPLOT, x, thruput-thruput_err,thruput+thruput_err, COLOR=90
  FOR i=0, n_date-1 DO XYOUTS, x[i], 1.2*MAX(thruput), date[i], ORIENTATION=90, CHARSIZE=0.8
  ; Close plot
  LOADCT, 0, /SILENT
  DEVICE, /CLOSE
  DEVICE, FILEN=plot_path + 'sensitivity.eps', /ENCAPS, /COLOR,  XSIZE=17.8, YSIZE=13.2, /TIMES
  xrange     = [5.,100]
  yrange     = [0.1,10.]
  PLOT, [0], [0], XTITLE='Observing run day', YTITLE='Photometric sensitivity [mJy/10min]', TITLE=title_day, XSTYLE=1, YSTYLE=1, XRANGE=[0,n_date+1], YRANGE=[0,1.5*MAX(sen)];, XTHICK=4, YTHICK=4, CHARTHICK=charthick, CHARSIZE=charsize
  LOADCT, 13, /SILENT
  OPLOT, x, sen, PSYM=ps, COLOR=90
  ERRPLOT, x, sen-sen_err,sen+sen_err, COLOR=90
  FOR i=0, n_date-1 DO XYOUTS, x[i], 1.2*MAX(sen), date[i], ORIENTATION=90, CHARSIZE=0.8
  ; Close plot
  LOADCT, 0, /SILENT
  DEVICE, /CLOSE & END_PS
END
