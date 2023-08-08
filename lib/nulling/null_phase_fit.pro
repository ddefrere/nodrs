@nodrs/lib/plot/colorbar.pro
;+
; NAME: NULL_PHASE_FIT
; 
; PURPOSE:
;   Main routine to call the statistical reduction technique.
;
; INPUTS:
;   data_in : a structure containing the four following fields
;        .i1        : the photometric signal for the first beam
;        .i2        : the photometric signal for the second beam
;        .flx       : the interferometric signal
;        .spdthpos  : the injected phase signal
;        .pcphmcos_t: mean COS(phase) over integration time T (superseed phi_var_t if not null)
;        .pcphmsin_t: mean SIN(phase) over integration time T (superseed phi_var_t if not null)
;
; KEYWORDS
;
; OUTPUT
;
; MODIFICATION HISTORY:
;  Version 1.0, 09-OCT-2017, by Denis Defrère, University of Arizona, denis@lbti.org

FUNCTION NULL_PHASE_FIT, data_in

  ; Keyword sanity check
  
  ; Extract useful data
  n_pt = N_ELEMENTS(data_in.flx)
  I1   = data_in.i1
  I2   = data_in.i2
  phot = I1+I2+2*SQRT(I1*I2)
  
  ; 1. Assign cycle ID number. 
  ; *************************
   
  ; Find unique values
  ph_pos      = data_in.spdthpos
  ph_pos_uniq = ph_pos[UNIQ(ph_pos,  SORT(ph_pos))]  
  IF N_ELEMENTS(ph_pos_uniq) LE 1 THEN BEGIN
    MESSAGE, 'No phase pattern detected', /CONTINUE
    RETURN, 0
  ENDIF

  ; Find where there is a jump to MAX phase
  ph_max  = MAX(ph_pos_uniq)
  idx_jmp = WHERE(ph_pos NE (SHIFT(ph_pos, -1)), n_jmp) + 1  ; Add 1 to have the starting frame
    
  ; Remove last frame if marked as jump
  IF idx_jmp[n_jmp-1] GE n_pt-1 THEN idx_jmp = idx_jmp[0:(n_jmp-2)]
  
  ; Phase value at changes
  ph_jmp  = ph_pos[idx_jmp]
  idx_max = WHERE(ph_jmp EQ ph_max)
  
  ; Add first frame if the sequence starts with ph_max (uncaptured by above algorithm)
  IF ph_pos[0] EQ ph_max AND idx_max[0] NE 0 THEN BEGIN  
    idx_max = [0, idx_max] 
    idx_jmp = [0, idx_jmp]
  ENDIF
   
  ; Now assign cycle ID
  ; Mark the files that don't follow this scheme as -1.
  n_cyc  = N_ELEMENTS(idx_max)
  cyc_id = INTARR(n_pt)-1
  FOR i_c = 0, n_cyc-2 DO cyc_id[idx_jmp[idx_max[i_c]]:idx_jmp[idx_max[i_c+1]]] = i_c
  cyc_id[idx_jmp[idx_max[n_cyc-1]]:*] = n_cyc-1
  
  ; 2. Now compute null estimate for each cycle
  ; *******************************************
  
  ; Prepare arrays
  null  = DBLARR(n_cyc)
  bias  = null
  setp  = null
  n_ph  = 1D5
  phase = -!DPI + 2*!DPI*DINDGEN(n_ph)/(n_ph-1)
  
  ; Loop over the cycle
  FOR ic = 0, n_cyc-1 DO BEGIN
    ; Extract data for this cycle
    idx_cyc = WHERE(cyc_id EQ ic)
    
    ; Identify center phase pattern
    ph_pos_cyc  = ph_pos[idx_cyc]
    ph_pos_uniq = ph_pos_cyc[UNIQ(ph_pos_cyc,  SORT(ph_pos_cyc))]
    n_pos       = N_ELEMENTS(ph_pos_uniq)
    IF n_pos NE 3 THEN BEGIN
      MESSAGE, 'Number of phase pattern should be 3!', /CONTINUE
      cyc_id[idx_cyc] = -1
    ENDIF ELSE BEGIN
      ; To follow Denis' math notations
      ph_pos0 = ph_pos_uniq[1] ;& PRINT, ph_pos0
      ph_pos1 = ph_pos_uniq[0]
      ph_pos2 = ph_pos_uniq[2]
      ; Compute phase dither
      ph_dith1 = ph_pos0-ph_pos1
      ph_dith2 = ph_pos2-ph_pos0
      IF ph_dith1 EQ ph_dith2 THEN BEGIN
        ; Index for each sequence of this cycle
        idx_seq0 = idx_cyc[WHERE(ph_pos_cyc EQ ph_pos0)] ;& PRINT, ph_pos[idx_seq0]
        idx_seq1 = idx_cyc[WHERE(ph_pos_cyc EQ ph_pos1)] ;& PRINT, ph_pos[idx_seq1]
        idx_seq2 = idx_cyc[WHERE(ph_pos_cyc EQ ph_pos2)] ;& PRINT, ph_pos[idx_seq2]
        ; Now extract sequence
        f0 = data_in.flx[idx_seq0] & AVGSDV, f0, m0, s0, sm0, KAPPA=5 ;& m0 = MEDIAN(f0)
        f1 = data_in.flx[idx_seq1] & AVGSDV, f1, m1, s1, sm1, KAPPA=5 ;& m1 = MEDIAN(f1)
        f2 = data_in.flx[idx_seq2] & AVGSDV, f2, m2, s2, sm2, KAPPA=5 ;& m2 = MEDIAN(f2)
        c0 = data_in.pcphmcos[idx_seq0] & AVGSDV, c0, mc0, sc0, sc0_m, KAPPA=5
        c1 = data_in.pcphmcos[idx_seq1] & AVGSDV, c1, mc1, sc1, sc1_m, KAPPA=5
        c2 = data_in.pcphmcos[idx_seq2] & AVGSDV, c2, mc2, sc2, sc2_m, KAPPA=5
        s0 = data_in.pcphmsin[idx_seq0] & AVGSDV, s0, ms0, ss0, ss0_m, KAPPA=5
        s1 = data_in.pcphmsin[idx_seq1] & AVGSDV, s1, ms1, ss1, ss1_m, KAPPA=5
        s2 = data_in.pcphmsin[idx_seq2] & AVGSDV, s2, ms2, ss2, ss2_m, KAPPA=5
        
;       ; Now compute values (ask Denis for math description)
;       ; 1. Background bias, which has an analytical solutio
;       ; Math approach (obsolote!!)
;       ; **************************
        
        bias[ic] = 0.5*(m1+m2-2*COS(ph_dith1)*m0)/(1-COS(ph_dith1))-(I1+I2)
        ; 2. Find best setpoint
        vec1 = (COS(phase)*ms0+SIN(phase)*mc0)/(COS(phase)*mc0+SIN(phase)*ms0)
        vec2 = 0.5*(m1-m2)/(m0-I1-I2-bias[ic])/SIN(ph_dith1)
        jnk  = MIN(ABS(vec1-vec2), idx_p)
        setp[ic] = phase[idx_p]
        ; 3. Visibility
        vis = 0.25*(m2-m1)/SQRT(I1*I2)/SIN(ph_dith1)/(COS(setp[ic])*ms0+SIN(setp[ic])*mc0)
        null[ic] = (1.-ABS(vis))/(1+ABS(vis))
        ; Plot diagnostic tool
        x = [ph_pos0, ph_pos1, ph_pos2]
        y = [f0, f1, f2]/phot*100
        s = [sm0, sm1, sm2]/phot*100
        xrange = 2*[MIN(x),MAX(x)]
        yrange = [0,5];[-2*ABS(MIN(y)),2*MAX(y)]
        d = I1+I2+2*ABS(vis)*SQRT(I1*I2)*(COS(phase-setp[ic])*mc0-SIN(phase-setp[ic])*ms0) + bias[ic]
        PLOT, x, y, XRANGE=xrange, YRANGE=yrange, PSYM=4, XSTYLE=1, YTITLE='Null depth [%]', XTITLE='Phase dither'
        OPLOT, phase, d, LINESTYLE=0      
        OPLOT, xrange, [bias[ic],bias[ic]], LINESTYLE=1
        ERRPLOT, x, y-s, y+s
        WAIT, 2
        ; Print info
        PRINT, 'Null [%], back bias [%], setpoint :', 1D2*null[ic], 1D2*bias[ic]/phot, setp[ic]/!Dpi*5.55

         ; Fit approach!
         ;contr_disk[i0,i1,i2,i3] = MPFIT('ERRNULL', guess_disc, BESTNORM=chi2, DOF=dof, FUNCTARGS={X:REFORM(data_base),Y:REFORM(data_vsqr),ERR:REFORM(data_evsq)}, PARINFO=parinfo, PERROR=perror, /QUIET)

                 
      ENDIF ELSE BEGIN
        MESSAGE, 'Phase dither should be uniform! Here ' + STRING(ph_dith1, FORMAT='(F5.2)') + STRING(ph_dith2, FORMAT='(F5.2)'), /CONTINUE 
        cyc_id[idx_cyc] = -1
      ENDELSE
    ENDELSE
  ENDFOR
  
  PRINT, bias/(I1+I2+2*SQRT(I1*I2))
  
  data_out = {NULL: null, BIAS: bias, SETP: setp}
  RETURN, data_out
END

PRO TEST_NULL
  i_aper     = 4
  ob_id      = 108 ; this is Vega for 2017-04-06_APR
  data_path  = '/Volumes/hosts_results/nomic/l1_fits/2017-04-06_APR/'
  ;ob_id      = 80  ; this is Altair for 2017-05-12_APR
  ;data_path  = '/Volumes/hosts_results/nomic/l1_fits/2017-05-12_APR/'
  data_bckg  = LBTI_READL1DATA(data_path, ob_id, 'BCKG', APER=apr_rad, IDX_APER=i_aper, /FILTER)
  data_phot1 = LBTI_READL1DATA(data_path, ob_id, 'PHOT1', APER=apr_rad, IDX_APER=i_aper, /FILTER)
  data_phot2 = LBTI_READL1DATA(data_path, ob_id, 'PHOT2', APER=apr_rad, IDX_APER=i_aper, /FILTER)
  data_null  = LBTI_READL1DATA(data_path, ob_id, 'NULL', APER=apr_rad, IDX_APER=i_aper, /FILTER)
  
  ; Prepare array
  AVGSDV, data_bckg.flx_tot[i_aper], mb, sb, KAPPA=5
  AVGSDV, data_phot1.flx_tot[i_aper], i1, s1, KAPPA=5
  AVGSDV, data_phot2.flx_tot[i_aper], i2, s2, KAPPA=5
  flx      = data_null.flx_tot[i_aper]
  pcphmean = data_null.pcphmean*!dtor*2.2/11.1
  pcphmcos = data_null.pcphmcos*COS(pcphmean)+data_null.pcphmsin*SIN(pcphmean)  ; mean cos(phase-<phase>) over DIT + backward compatibility
  pcphmsin = data_null.pcphmsin*COS(pcphmean)-data_null.pcphmcos*SIN(pcphmean)
  spdthpos = data_null.spdthpos*!dtor*2.2/11.1
  
  idx_ok = WHERE(pcphmcos GE -1 AND pcphmcos LE 1 AND pcphmsin GE -1 AND pcphmsin LE 1)  
  data_in = {i1: i1, i2: i2, flx: flx[idx_ok], spdthpos: spdthpos[idx_ok], pcphmcos: pcphmcos[idx_ok], pcphmsin: pcphmsin[idx_ok]}
  

  plot,spdthpos
  d = NULL_PHASE_FIT(data_in)
  PLOT, d.null
  STOP
  
  ; Compute dither pattern for 5 sigma
  flx_star = 38.55 ; Vega
  sig_null = sb/(i1+i2+2*SQRT(i1*i2))
  
  n_pt    = 5
  flx_all = DINDGEN(40)
  rms_all = sig_null*flx_star/flx_all/SQRT(n_pt)
  PLOT, flx_all, SQRT(4*5*rms_all), XTITLE='Target flux [Jy]', YTITLE='N-band phase [rad]', TITLE='Minimum dither phase step'
  OPLOT, flx_all, SQRT(4*3*rms_all), LINESTYLE=1
  
  WINDOW, 2
  d = data_null.pcphmean-data_null.spdthpos
  diff = d-SHIFT(d,1)
  plot, diff, YTITLE='K-band phase difference [deg]', XTITLE='Frama ID'
  AVGSDV, diff, m, s
  PRINT, 'Difference between succesive frames [deg]', s
  
END

  



