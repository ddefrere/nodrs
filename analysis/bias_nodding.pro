PRO BIAS_NODDING
  ; Basic routine for background bias analysis
  date_obs  = '2015-02-08' ; '2018-03-30' ;'2017-02-09' ;'2015-02-08'  
  which     = 'bckg'        ; use 'bckg/' for background subtracted images
  
  ; Read image cube
  data_path = '/Volumes/hosts_results/nomic/l0_fits/' + date_obs + '/'
  datalog   = data_path + 'datalog.sav'
  data_path = data_path + which + '/'     
    
  ; File list
  files = FILE_SEARCH(data_path, '*IMG.fits', COUNT=n_files) 

  ; Image and PSF parameters (pick center of bottom right quadrant)
  xcen  = 256-64;+n_pix/2;+5
  ycen  = 64;+n_pix/2-0.4
  fwhm  = 16
  rad   = fwhm/2
  bck_irad1 = 31
  bck_irad2 = 31
  bck_orad1 = 40
  bck_orad2 = 40
  clipsig   = 5
  sky_col1  = 1
  sky_col2  = 0
  
  flx_avg = DBLARR(n_files, 2)
  bck_avg = flx_avg
  flx_rms = flx_avg
  bck_rms = flx_avg
  flx_all = 0
  bck_all = 0
  el_all  = 0
  ; Loop over files
  PRINT, 'Number of image cubes to process', n_files
  FOR i_f = 0, n_files-1 DO BEGIN
    img   = READFITS(files[i_f], hdr, /SILENT)
    n_img = N_ELEMENTS(img[0,0,*])
    tmp_avg = DBLARR(n_img, 2)
    tmp_rms = flx_avg
    tmp_avg2 = DBLARR(n_img, 2)
    tmp_rms2 = flx_avg
    FOR i_img = 0, n_img-1 DO BEGIN
      shft_x = 0;i_img
      shft = 0;i_img; RANDOMN(seed)
      k = img[*,*,i_img];+ sig_noise*RANDOMN(seed, n_pix, n_pix);+flx_star*PSF_GAUSSIAN(NPIXEL=[n_pix,n_pix], FWHM=fwhm, /double, /NORMALIZE)
      APER_WEIGHT, k, xcen+shft_x, ycen+shft, flx_tot3, flx_err0, bck_flx3, bck_err0, 154/0.4, [rad], [bck_irad1,bck_orad1], [-10000,20000], SKYLIM=skylim, /FLUX, /EXACT, /SILENT, /MEANBACK, CLIPSIG=clipsig, SKY_COL=sky_col1
      tmp_avg[i_img, 0] = flx_tot3   
      tmp_avg2[i_img, 0] = bck_flx3
      APER_WEIGHT, k, xcen+shft_x, ycen+shft, flx_tot4, flx_err0, bck_flx4, bck_err0, 154/0.4, [rad], [bck_irad2,bck_orad2], [-10000,20000], SKYLIM=skylim, /FLUX, /EXACT, /SILENT, /MEANBACK, CLIPSIG=clipsig, SKY_COL=sky_col2;, /sky_weight
      tmp_avg[i_img, 1] = flx_tot4
      tmp_avg2[i_img, 1] = bck_flx4
      flx_all = [flx_all, flx_tot4]
      bck_all = [bck_all, bck_flx4]   
    ENDFOR   
    data = MRDFITS(files[i_f], 1, hdr, /SILENT)
    el_all  = [el_all, data.lbt_alt]
    AVGSDV, tmp_avg[*,1], avg, rms, KAPPA=5
    flx_avg[i_f,1] = avg
    flx_rms[i_f,1] = rms/SQRT(n_img)
    AVGSDV, tmp_avg2[*,1], avg, rms, KAPPA=5
    bck_avg[i_f,1] = avg
    bck_rms[i_f,1] = rms/SQRT(n_img)
    AVGSDV, tmp_avg[*,0], avg, rms, KAPPA=5
    flx_avg[i_f,0] = avg
    flx_rms[i_f,0] = rms/SQRT(n_img)
    AVGSDV, tmp_avg2[*,0], avg, rms, KAPPA=5
    bck_avg[i_f,0] = avg
    bck_rms[i_f,0] = rms/SQRT(n_img)
    PRINT, i_f 
  ENDFOR
  
  flx_all = flx_all[1:*]
  bck_all = bck_all[1:*]
  
  coeff = POLY_FIT(bck_all, flx_all, 2)
  
  ; Now correct!!
  ; Restore datalog
  RESTORE, datalog
  pt_id = data_r.pt_id
  pt_uniq = pt_id[UNIQ(pt_id, SORT(pt_id))]
  n_pt = N_ELEMENTS(pt_uniq)
  rms_pt = DBLARR(n_pt,2)
  rms_pt2 = DBLARR(n_pt,2)
  flx_cor = flx_avg
  f = LINDGEN(n_files)
  FOR i_pt = 0, n_pt-1 DO BEGIN
    idx = f[WHERE(pt_id EQ i_pt, nt)]
    AVGSDV, flx_avg[idx,0], avg, rms, KAPPA=5
    rms_pt[i_pt, 0] = rms
    AVGSDV, flx_avg[idx,1], avg, rms, KAPPA=5
    rms_pt[i_pt, 1] = rms
    ; Correct!
    FOR ii = 0, nt-1 DO BEGIN
      flx_cor[idx[ii],0] = flx_avg[idx[ii],0]-(coeff[2]*bck_avg[idx[ii],0]^2+coeff[1]*bck_avg[idx[ii],0]+coeff[0])
      flx_cor[idx[ii],1] = flx_avg[idx[ii],1]-(coeff[2]*bck_avg[idx[ii],1]^2+coeff[1]*bck_avg[idx[ii],1]+coeff[0])
    ENDFOR
    AVGSDV, flx_cor[idx,0], avg, rms, KAPPA=5
    rms_pt2[i_pt, 0] = rms
    AVGSDV, flx_cor[idx,1], avg, rms, KAPPA=5
    rms_pt2[i_pt, 1] = rms
  ENDFOR
  
  ; Now compute MEAN RMS over all pointings
  AVGSDV, rms_pt[*,0], avg1, rms, KAPPA=5
  AVGSDV, rms_pt[*,1], avg2, rms, KAPPA=5
  PRINT, 'POINTING-BASED APPROACH'
  PRINT, 'Uncorrected RMS of background estimates per pointing (SKY COLUMNS, ANNULUS) :', avg1, avg2
  AVGSDV, rms_pt2[*,0], avg1, rms
  AVGSDV, rms_pt2[*,1], avg2, rms
  PRINT, 'Corrected RMS of background estimates per pointing (SKY COLUMNS, ANNULUS) :', avg1, avg2
  
  WINDOW, 1
  PLOT, INDGEN(n_files), flx_avg[*,0], PSYM=4, XTITLE='Nod ID', YTITLE='Flux measurements per NOD';, YRANGE=[-700,100]
  ERRPLOT, INDGEN(n_files), flx_avg[*,0]-flx_rms[*,0], flx_avg[*,0]+flx_rms[*,0]
  LOADCT, 13, /SILENT
  ;OPLOT, INDGEN(n_files), flx_avg[*,1], PSYM=5, COLOR=240
  ;ERRPLOT, INDGEN(n_files), flx_avg[*,1]-flx_rms[*,1], flx_avg[*,1]+flx_rms[*,1], COLOR=240
  AVGSDV, flx_avg[*,0], avg0, rms0
  AVGSDV, flx_avg[*,1], avg1, rms1
  PRINT, 'OB-BASED APPROACH'
  PRINT, 'RMS of background estimates (SKY COLUMNS, ANNULUS) :', rms0, rms1
  AVGSDV, flx_rms[*,0], avg0, rms0
  AVGSDV, flx_rms[*,1], avg1, rms1
  PRINT, 'Expected RMS of background estimates (SKY COLUMNS, ANNULUS) :', avg0, avg1  

  n = N_ELEMENTS(flx_all)
  WINDOW, 2
  PLOT, LINDGEN(n), flx_all, XTITLE='File ID', YTITLE='Flux measurements'
  
  WINDOW, 3
  PLOT, LINDGEN(n), bck_all, YRANGE=[7000,8000], YSTYLE=1 , XTITLE='File ID', YTITLE='Mean background in annulus'
  
  WINDOW, 4
  PLOT, bck_all, flx_all, PSYM=4, XTITLE='Mean background in annulus', YTITLE='FLux measurements'
  
  WINDOW, 5
  PLOT, LINDGEN(n), el_all, XTITLE='File ID', YTITLE='LBT Elevation'
  
  WINDOW, 6
  PLOT, el_all, flx_all, XTITLE='LBT Elevation', YTITLE='FLux measurements'
  
  WINDOW, 7
  PLOT, bck_avg[*,0], flx_avg[*,0], PSYM=4, XTITLE='Mean background level in annulus', YTITLE='Background bias'
  ERRPLOT, bck_avg[*,0], flx_avg[*,0]-flx_rms[*,0], flx_avg[*,0]+flx_rms[*,0]
END