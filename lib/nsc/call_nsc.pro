@nodrs/lib/plot/colorbar.pro
;+
; NAME: CALL_NSC
; 
; PURPOSE:
;   Main routine to call the statistical reduction technique.
;
; INPUTS:
;   arraydata : a structure containing at least the four following fields
;        .dark      : the background measurements
;        .a         : the photometric signal for the first beam
;        .b         : the photometric signal for the second beam
;        .ab        : the interferometric signal
;        .phi_avg_t : mean phase over integration time T
;        .phi_var_t : phase variance over integration time T
;        .pcphmcos_t: mean COS(phase) over integration time T (superseed phi_var_t if not null)
;        .pcphmsin_t: mean SIN(phase) over integration time T (superseed phi_var_t if not null)
;
; KEYWORDS
;   BANDWIDTH    :  Bandwidth [m]
;   CUBE_SIZE    :  3-element vector containing the number of elements in the chi2 cube along the null, mean phase and rms phase directions respectively
;                   If not set, the default values of STAT_RED.pro will be used.
;   FILE_ID      :  File ID number for plots
;   N_BOOTSTRAP  :  Number of bootstrap samples. If 0 or 1, the error bar is computed by a Bayesian approach (DEFAULT approach now)
;   NBINS_FACTOR :  Nominaly, if the number of recorded null points is Np, you use sqrt(Np) bins in the initial histogram,
;                   which corresponds to nbins_factor=1. The number of histogram bins can in general be set to nbins_factor*sqrt(Np)
;   NO_MULTI     :  Set this keyword to turn off multi-threading
;   NULL_RANGE   :  Source null range (to scan with NSC)    
;   OCC_MIN      :  Minimum number of occurences per bin for the fit
;   SEED_OPD     :  Defines the seed of the random OPD generator used in the simulated null sequence and histogram.
;                   Typically, seed2 is undefined so that it uses the computer's system clock
;   VARIABLE     :  Bin size for NSC reduction (0: constant, 1: variable)
;   WAV_EFF      :  Effective wavelength [m]
;   INFO         :  Set this keyword to print info to screen
;   EPS          :  Save files as EPS files rather than plot them on screen
;   PLOT         :  Turn on/off plots
;   FILE_PATH    :  Set this keyword to the path where EPS figures will be saved (workspace directory by default)
;
; OUTPUT
;   An 4-element structure with the following sub-structure:
;       .null_opt:  Optimum bets-fit values from likelihood estimator
;       .null_bay:  Optimum values after marginalization (Bayesian approach)
;       .null_avg:  Average of bootstrap samples
;       .null_mod:  Mode of bootstrap samples
;  Each of these sub-structures is a structure containing the following elements:
;       .nas          : astrophysical null
;       .nas_err_low  : inferior error
;       .nas_err_sup  : superior error
;       .mu
;       .mu_err_low
;       .mu_err_sup
;       .sig
;       .sig_err_low
;       .sir_err_sup
;       .kdrk         : background RMS multiplication factor
;       .kdrk_err_low  
;       .kdrk_err_sup
;       .chi2
;
; MODIFICATION HISTORY:
;  Version 1.0, 28-MAY-2014, by Denis Defrère, University of Arizona, denis@lbti.org
;  Version 1.1, 28-NOV-2014, DD: added keyword FILE_ID and EPS. Now return best-fit phase and chi2
;  Version 1.2, 29-NOV-2014, DD: added keyword CUBE_SIZE and FILE_PATH + add distribution plots
;  Version 1.3, 13-DEC-2014, DD: added standard errors to output structure
;  Version 1.4, 14-JAN-2015, DD: added keyword SEED_OPD
;  Version 1.5, 13-MAR-2015, DD: added multi-threading functionality
;  Version 1.6, 14-MAR-2015, DD: implemtented Bayesian approach
;  Version 1.7, 01-APR-2015, DD: added PDF plot
;  Version 1.8, 31-JUL-2015, DD: modified comments, cleaned implementation, and modified multi-threading approach (now on the cube rather than
;                                the bootstrap sample => work with Bayesian approach)
;  Version 1.9, 01-AUG-2015, DD: added coarse approach to find approximate initial values for each parameter
;  Version 2.0, 04-AUG-2015, DD: added keyword NBINS_FACTOR and OCC_MIN which proves critical!
;  Version 2.1, 08-AUG-2015, DD: now make sure that the inputs in arraydata have the right size!
;  Version 2.2, 24-AUG-2015, DD: adapted for new format of arraydata
;  Version 2.3, 15-NOV-2015, DD: added plot of the coarse chi2 maps
;  Version 2.3, 03-APR-2016, DD: now works over SSH + work for n_btsp of 0 or 1
;  Version 2.4, 06-JUL-2016, DD: now display results properly when n_btsp is 0 or 1
;  Version 2.5, 08-JUL-2016, DD: added keyword NULL_RANGE
;  Version 2.6, 08-NOV-2016, DD: now returns two null estimates: the mode and the mean of the null distributions
;  Version 2.6, 15-NOV-2016, DD: added PCPHMCOS_T and PCPHMSIN_T to input structure
;  Version 2.7, 29-NOV-2016, DD: improved error bar computation
;  Version 2.8, 31-JAN-2017, DD: Added dither pattern and optimized code speed
;  Version 2.9, 13-FEB-2017, DD: Updated for new format of STAT_RED
;  Version 3.0, 22-FEB-2017, DD: Added best-fit background multiplication factor to output structure

FUNCTION CALL_NSC, arraydata, BANDWIDTH=bandwidth, CUBE_SIZE=cube_size, FILE_ID=file_id, N_BOOTSTRAP=n_bootstrap, NBINS_FACTOR=nbins_factor, NO_MULTI=no_multi, NULL_RANGE=null_range, OCC_MIN=occ_min, $
                   SEED_OPD=seed_opd, VARIABLE=variable, WAV_EFF=wav_eff,  $                                                   ; Input keywords
                   EPS=eps, PLOT=plot, FILE_PATH=file_path, INFO=info                                                          ; Plot and print keywords
  
; Sanity checks
IF NOT KEYWORD_SET(FILE_ID)       THEN file_id      = 0
IF NOT KEYWORD_SET(N_BOOTSTRAP)   THEN n_btsp       = 1 ELSE n_btsp = n_bootstrap
IF NOT KEYWORD_SET(NBINS_FACTOR)  THEN nbins_factor = 1
IF NOT KEYWORD_SET(NO_MULTI)      THEN no_multi     = 0
IF NOT KEYWORD_SET(NULL_RANGE)    THEN null_range   = [-0.01,0.04]  ; make sure it works for every star (-2% seems to be necessary for old data where the high=frequency phase noise is overestimated)
IF N_ELEMENTS(null_range) NE 2    THEN MESSAGE, 'NULL_RANGE must be a two elements vector. Check your config file'
IF NOT KEYWORD_SET(OCC_MIN)       THEN occ_min      = 0             ; Set to remove bins with less than truncate occurences from the fit
IF NOT KEYWORD_SET(VARIABLE)      THEN variable     = 0  
IF NOT KEYWORD_SET(INFO)          THEN info         = 0
IF NOT KEYWORD_SET(PLOT)          THEN plot         = 0
IF KEYWORD_SET(EPS)               THEN ps           = 1 

; Consistancy checks
n_null = N_ELEMENTS(arraydata.ab)
IF N_ELEMENTS(arraydata.phi_avg_t)  NE n_null THEN MESSAGE, 'The mean phase vector should have the same size as the null vector'
IF N_ELEMENTS(arraydata.phi_var_t)  NE n_null THEN MESSAGE, 'The variance phase vector should have the same size as the null vector'
IF N_ELEMENTS(arraydata.pcphmcos_t) NE n_null THEN MESSAGE, 'The COS phase vector should have the same size as the null vector'
IF N_ELEMENTS(arraydata.pcphmsin_t) NE n_null THEN MESSAGE, 'The SIN phase vector should have the same size as the null vector'
IF N_ELEMENTS(arraydata.spdthpos)   NE n_null THEN MESSAGE, 'The dither pattern vector should have the same size as the null vector'

; Running paramaters
neg      = 1    ; Allow negative values for astro null  
               
; Open NSC log file
OPENW, lun, file_path + 'nsc_log.txt', /GET_LUN, WIDTH=800
PRINTF,lun, ' '
PRINTF,lun, 'NSC log file'
PRINTF,lun, ' '

; Ensure that the photometries and the background have at least the same number of points as the null.
; If not, increase their size. Also check that the photometries and the background have the same size
n_phot1 = N_ELEMENTS(arraydata.a)
n_phot2 = N_ELEMENTS(arraydata.b)
n_bckg  = N_ELEMENTS(arraydata.dark2)
;n_aux   = n_bckg < n_phot1 < n_phot2
IF n_phot1 LT n_null THEN a = [arraydata.a, arraydata.a[FIX(RANDOMU(!NULL,n_null-n_phot1)*(n_phot1-1))]]       ELSE a = arraydata.a[0:(n_null-1)]
IF n_phot2 LT n_null THEN b = [arraydata.b, arraydata.b[FIX(RANDOMU(!NULL,n_null-n_phot2)*(n_phot2-1))]]       ELSE b = arraydata.b[0:(n_null-1)]
IF n_bckg  LT n_null THEN c = [arraydata.dark2, arraydata.dark2[FIX(RANDOMU(!NULL,n_null-n_bckg)*(n_bckg-1))]] ELSE c = arraydata.dark2[0:(n_null-1)]

; Store in array
nas_min   = arraydata.nas_min 
phi0_min  = 1D-5
phie_min  = 1D-5;arraydata.phi_rms_min > 0.01 DON'T USE PHI_MIN RIGHT NOW BECAUSE UNRELIABLE WITH FRINGE JUMPS
spdthpos  = arraydata.spdthpos
tmpdata   = {DARK2: c, A: a, B: b, AB: arraydata.ab, Npts: n_null, PHI_AVG_T: arraydata.phi_avg_t, PHI_VAR_T: arraydata.phi_var_t, PCPHMCOS_T: arraydata.pcphmcos_t, PCPHMSIN_T: arraydata.pcphmsin_t, SPDTHPOS: spdthpos}
arraydata = tmpdata

; Compute maximum null precision based on background noise and number of NULL points
; This will be used to compute the steps in the NSC cube
phot_tot = MEAN(a+b+2*SQRT(ABS(a*b))) 
AVGSDV, c, avg_bck, rms_bck, KAPPA=3
null_pre = rms_bck/phot_tot/SQRT(n_null)

; 1. Run COARSE NSC to get closer to the optimium values
; ******************************************************

; Init time cpunter
t0   = SYSTIME(1)

; Minimum phase jitter (measured by frimge tracker)
null_range[0] = null_range[0] > nas_min

; Ranges for the coarse runs (here, we don't scan over the OPD seed and the DARK rms)
cube_coarse = [100,100,100,1,1]   ; Hard-coded size (takes ~10 seconds)
null_lim_c  = null_range          ; From config file. Should work for all stars (but be careful for eta Crv and zet Lep)
phi0_lim_c  = [phi0_min,1.5]      ; This is an offset of ~lambda/2 off the null. Should be safe (if fit doesn't work with that range, the point is likely bad anyway. Check chi2!) 
phie_lim_c  = [phie_min,1.0]      ; This is an OPD jitter of ~3.5 microns. Should be safe in all cases. 

; Run coarse computation (no need to bootstrap, nor to multi-thread. We are only interested in approximate starting values for the three parameters).
; If derived null is fishy, inspect Nas vs chi2 plot to see if it's not stuck in a local minimum. 
STAT_RED, arraydata, param_coarse, $
          AUTO_FWHM=1, BANDWIDTH=bandwidth, BIN_FACT=nbins_factor, BOOTSTRAP=0, CUBE_SIZE=cube_coarse, FILE_ID=file_id, NEG=neg, NFWHM_LEFT=nfwhm_left, NFWHM_RIGHT=nfwhm_right, $
          NULL_LIM=null_lim_c, PHI0_LIM=phi0_lim_c, PHIE_LIM=phie_lim_c, SEED_BOOT=seed, SEED_OPD=seed_opd, SKIPEND=skipend, TRSH=occ_min, VARIABLE=variable, WAV_EFF=wav_eff,$   ; Input keywords                                                                                                          ; Output keywords
          INFO=0, PLOT=plot, FILE_PATH= file_path, FLAG='COARSE', /PS

; Define ranges for FINE run below. The rule is to keep only the paramater space where chi2 < MIN(chi2) + 5
nas_coarse = param_coarse.nas.lklh.bay
mu_coarse  = param_coarse.mu.lklh.bay
sig_coarse = param_coarse.sig.lklh.bay
null_lim   = nas_coarse + 10*[-param_coarse.nas.lklh.err_low, param_coarse.nas.lklh.err_sup] & null_lim[0] = null_lim[0] > nas_min  ;& null_lim[1] = null_lim[1] < null_lim_c[1]
phi0_lim   = mu_coarse  + 10*[-param_coarse.mu.lklh.err_low, param_coarse.mu.lklh.err_sup]   & phi0_lim[0] = phi0_lim[0] > phi0_min ;& phi0_lim[1] = phi0_lim[1] < null_lim_c[1]    
phie_lim   = sig_coarse + 10*[-param_coarse.sig.lklh.err_low, param_coarse.sig.lklh.err_sup] & phie_lim[0] = phie_lim[0] > phie_min ;& phie_lim[1] = phie_lim[1] < null_lim_c[1]

; Now compute the cube size if automatic
; Also make sure there is enough memory by imposing a limit on the cube size 
; Also impose a lower limit on the number of elements
IF cube_size[0] EQ 0 THEN n1 = (ROUND((null_lim[1]-null_lim[0])/null_pre) < 250) > 50                ELSE n1 = cube_size[0]
IF cube_size[1] EQ 0 THEN n2 = (ROUND((phi0_lim[1]-phi0_lim[0])/(2*null_pre/mu_coarse)) < 250) > 50  ELSE n2 = cube_size[1] ; N ~ mu^2/4 => dN = mu*dmu/2 => dmu = 2*dN/mu
IF cube_size[2] EQ 0 THEN n3 = (ROUND((phie_lim[1]-phie_lim[0])/(2*null_pre/sig_coarse)) < 200) > 25 ELSE n3 = cube_size[2] ; N ~ sig^2/4 => dN = sig*dsig/2 => dsig = 2*dN/sig
IF cube_size[3] EQ 0 THEN n4 = 50 ELSE n4 = cube_size[3]
IF cube_size[4] EQ 0 THEN n5 = 1  ELSE n5 = cube_size[4]
cube_fine = [n1,n2,n3,n4,n5] ;& PRINT, cube_fine

; Print info to screen
IF info GT 0 THEN PRINT, 'Time for coarse NSC reduction of this OB [s]', SYSTIME(1)-t0


; 2. Run FINE NSC and compute error bar
; *************************************

; Initiate time counter
t1 = SYSTIME(1)

; Run coarse computation (no need to bootstrap, nor to multi-thread. We are only interested in approximate starting values for the three parameters).
; If derived null is fishy, inspect Nas vs chi2 plot to see if it's not stuck in a local minimum.
STAT_RED, arraydata, param_fine, $
          AUTO_FWHM=1, BANDWIDTH=bandwidth, BIN_FACT=nbins_factor, BOOTSTRAP=0, CUBE_SIZE=cube_fine, FILE_ID=file_id, NEG=neg, NFWHM_LEFT=nfwhm_left, NFWHM_RIGHT=nfwhm_right, $
          NULL_LIM=null_lim, PHI0_LIM=phi0_lim, PHIE_LIM=phie_lim, SEED_BOOT=seed, SEED_OPD=seed_opd, SKIPEND=skipend, TRSH=occ_min, VARIABLE=variable, WAV_EFF=wav_eff,$   ; Input keywords                                                                                                          ; Output keywords
          INFO=0, PLOT=plot, FILE_PATH=file_path, FLAG='FINE', /PS

; Print info to screen
IF info GT 0 THEN PRINT, 'Time for fine NSC reduction of this OB [s]', SYSTIME(1)-t1

; 3. Compute error bar by bootstrapping if requested
; **************************************************

; Now compute error bar (or not if n_btstp is 0 or 1)
IF n_btsp GT 1 THEN BEGIN
  MESSAGE, 'COMPUTING ERROR BARS BY BOOTSTRAPPING IS DEPRECIATED. ASK DENIS IF NEEDED. EXPERIENCE SHOWS THAT BAYESIAN APPROACH LEADS TO SIMILAR ERROR BARS AS BOOTSTRASPPING.'
  ; Print info to screen
  IF info GT 0 THEN PRINT, 'Running multi-thread NSC reduction. This may take a while...'
  
  ; Run // reductions
  ; First, find best OPD vector for each position on the cube
  n_cpu = 2^(FLOOR(ALOG(!CPU.HW_NCPU-1)/ALOG(2))) < n_btsp  ; By experience, use 2^n cores where 2^n is smaller than the number of available cores
  n_btsp -= (n_btsp MOD n_cpu)
  btsp_lst = INDGEN(n_btsp)
  SPLIT_FOR, 0, n_btsp-1, commands=[$
            'STAT_RED, arraydata, nas, mu, sig, chi2, lklh, AUTO_FWHM=1, BANDWIDTH=bandwidth, BIN_FACT=nbins_factor, BOOTSTRAP=0, CUBE_SIZE=cube_fine, FILE_ID=file_id, FILEPHI=filephi, NEG=neg,'+ $
            'NFWHM_LEFT=nfwhm_left, NFWHM_RIGHT=nfwhm_right, NULL_LIM=null_lim, PHI0_LIM=phi0_lim, PHIE_LIM=phie_lim, DARK_RMS=0, RNDM_OPD=rndm_opd, SEED_BOOT=seed,' +$
            'SEED_OPD=btsp_lst[i], SKIPEND=skipend, TRSH=occ_min, VARIABLE=variable, WAV_EFF=wav_eff,' +$      ; Input keywords
            'HIST_MODEL=hist_model_opt, CHI2CUBE=cube, LKLH_CUBE=lcube, DARK_SCALE=drk, NULL_EXPL=nullexpl, MU_EXPL=muexpl, SIG_EXPL=sigexpl, NBINS_KEPT=nbins_kept, HISTOBS=histobs, NULL_LOC=null_loc,'+$                                               ; Output keywords
            'INFO=info, PLOT=plot, FILE_PATH=file_path, PS=ps' + $
            '& IF N_ELEMENTS(nasopt)  EQ 0 THEN nasopt   = nas   ELSE nasopt   = [nasopt,nas]' + $
            '& IF N_ELEMENTS(muopt)   EQ 0 THEN muopt    = mu    ELSE muopt    = [muopt,mu]' + $
            '& IF N_ELEMENTS(sigopt)  EQ 0 THEN sigopt   = sig   ELSE sigopt   = [sigopt,sig]' + $
            '& IF N_ELEMENTS(chi2opt) EQ 0 THEN chi2opt  = chi2  ELSE chi2opt  = [chi2opt,chi2]'+ $
            '& IF N_ELEMENTS(lklhopt) EQ 0 THEN lklhopt  = lklh  ELSE lklhopt  = [lklhopt,lklh]'+ $
            '& IF N_ELEMENTS(drkopt)  EQ 0 THEN drkopt   = drk   ELSE drkopt   = [drkopt,drk]'+ $
            '& IF N_ELEMENTS(lcubeopt) EQ 0 THEN lcubeopt = lcube ELSE lcubeopt = [lcubeopt,lcube]'+ $
            '& IF N_ELEMENTS(cubeopt) EQ 0 THEN cubeopt  = cube  ELSE cubeopt  = [cubeopt,cube]'], $
            struct2pass1=arraydata,$
            varnames=['nbins_factor','file_path','file_id','neg','cube_fine','occ_min','variable','wav_eff','rndm_opd','null_lim','phi0_lim','phie_lim','btsp_lst'],$ 
            outvar=['nasopt','muopt','sigopt','chi2opt','lklhopt','lcubeopt','cubeopt','nullexpl','muexpl','sigexpl','drkopt'], NSPLIT=n_cpu, /SILENT
  ; Print status to screenprint,
  IF info GT 0 THEN PRINT, 'Multi-threading done'
  ; Parse results
  bestnull  = [nasopt0]
  mu_phase  = [muopt0]
  sig_phase = [sigopt0]
  bestlklh  = [lklhopt0]
  bestchi2  = [chi2opt0]
  drk_scale = [drkopt0]
  chi2cube  = [cubeopt0]
  lklhcube  = [lcubeopt0]
  null_expl = [nullexpl0]
  mu_expl = [muexpl0]
  sig_expl = [sigexpl0]
  i         = 1
  WHILE EXECUTE('tmp = nasopt'  + STRING(i, FORMAT='(I0)')) DO BEGIN
    IF EXECUTE('tmp = nasopt'   + STRING(i, FORMAT='(I0)')) THEN bestnull  = [bestnull,tmp]
    IF EXECUTE('tmp = muopt'    + STRING(i, FORMAT='(I0)')) THEN mu_phase  = [mu_phase,tmp]
    IF EXECUTE('tmp = sigopt'   + STRING(i, FORMAT='(I0)')) THEN sig_phase = [sig_phase,tmp]
    IF EXECUTE('tmp = drkopt'   + STRING(i, FORMAT='(I0)')) THEN drk_scale = [drk_scale,tmp]
    IF EXECUTE('tmp = chi2opt'  + STRING(i, FORMAT='(I0)')) THEN bestchi2  = [bestchi2,tmp]
    IF EXECUTE('tmp = lklhopt'  + STRING(i, FORMAT='(I0)')) THEN bestlklh  = [bestlklh,tmp]
    IF EXECUTE('tmp = lcubeopt' + STRING(i, FORMAT='(I0)')) THEN lklhcube  = [lklhcube,tmp]
    IF EXECUTE('tmp = cubeopt'  + STRING(i++, FORMAT='(I0)')) THEN chi2cube  = [chi2cube,tmp]
  ENDWHILE
  
  ; Find minimum chi2cube along bootstrap position (different for each cube position)
  chi2out = FLTARR(cube_fine[0],cube_fine[1],cube_fine[2],/NOZERO)
  chi2idx = INTARR(cube_fine[0],cube_fine[1],cube_fine[2],/NOZERO)
  n_btsp  = N_ELEMENTS(chi2cube[*,0,0])/cube_fine[0]
  FOR i=0L, cube_fine[0]-1 DO BEGIN
    chi2slice      = chi2cube[i+cube_fine[0]*LINDGEN(n_btsp),*,*]
    chi2out[i,*,*] = MIN(chi2slice, idxmin, DIMENSION=1)
    chi2idx[i,*,*] = (ARRAY_INDICES(chi2slice, idxmin))[0,*,*]
  ENDFOR
  chi2cube = chi2out
  
  ; Find best paramaters
  chi2_fine = MIN(bestchi2, idx_best)
  nas_chi2  = bestnull[idx_best]
  mu_chi2   = mu_phase[idx_best]
  sig_chi2  = sig_phase[idx_best]
  
  ; Do the same for the likelyhood cube
  lklhout = FLTARR(cube_fine[0],cube_fine[1],cube_fine[2],/NOZERO)
  lklhidx = INTARR(cube_fine[0],cube_fine[1],cube_fine[2],/NOZERO)
  FOR i=0L, cube_fine[0]-1 DO BEGIN
    lklhslice      = lklhcube[i+cube_fine[0]*LINDGEN(n_btsp),*,*]
    lklhout[i,*,*] = MIN(lklhslice, idxmin, DIMENSION=1)
    lklhidx[i,*,*] = (ARRAY_INDICES(lklhslice, idxmin))[0,*,*]
  ENDFOR
  lklhcube = lklhout
  
  ; Marginalize null
  CHI2STAT, chi2cube, null_expl, mu_expl, sig_expl, nas_chi2_fine, mu_chi2_fine, sig_chi2_fine, DOF=(nbins_kept-4), PLOT_PATH=file_path, TAG='CHI2_FINE'
  LKLHSTAT, EXP(-0.5*lklhcube), null_expl, mu_expl, sig_expl, nas_lklh_fine, mu_lklh_fine, sig_lklh_fine, PLOT_PATH=file_path, TAG='LKLH_FINE
  
  ; Extract value
  null_bay = nas_lklh_fine.nas_bay & null_bay_err_low = nas_lklh_fine.nas_err_low & null_bay_err_sup = nas_lklh_fine.nas_err_sup
  mu_bay   = mu_lklh_fine.mu_bay   & mu_bay_err_low   = mu_lklh_fine.mu_err_low   & mu_bay_err_sup   = mu_lklh_fine.mu_err_sup
  sig_bay  = sig_lklh_fine.sig_bay & sig_bay_err_low  = sig_lklh_fine.sig_err_low & sig_bay_err_sup  = sig_lklh_fine.sig_err_sup
    
  ; Chi2 at max L
  idx_tmp = WHERE(lklhcube EQ MIN(lklhcube))
  lklhchi2_1 = chi2cube[idx_tmp[0]]
  chi2_bay = lklhchi2_1
  
  ; Find best paramaters
  lklh_1   = MIN(bestlklh, idx_best)
  nas_lklh = bestnull[idx_best]
  mu_lklh  = mu_phase[idx_best]
  sig_lklh = sig_phase[idx_best]
    
  ; Print info to screen
  IF info GT 0 THEN PRINT, 'Time for OPD optimization of this OB [s]', SYSTIME(1)-t1
  
  ; Find best bootstrap sample for each position on the grid
  tmp = MIN(lklhcube,index)
  index3D = ARRAY_INDICES(lklhcube, index[0]) ; translate the index of the minimum into the 3D coordinates within the chi2 cube
  iopt = index3D[0] ; index of nas that minimizes chi2 -- min Na
  jopt = index3D[1] ; for the min(na) ,value of mu that minimizes chi2 -- min mean phase
  kopt = index3D[2] ; value of sigma that minimizes chi2 -- min  phase_rms
  
  ; Now run data bootstrap to find error bars. I used to pass a cube as best_opd_seed but this is ~9x more time consuming.
  ; Now I just pass the best OPD seed as a single value. MARGINLAZE because it usually provides a better estimator
  ;best_opd_seed = lklhidx ;chi2idx ; old approach more time consuming
  best_opd_seed = lklhidx[iopt,jopt,kopt]
  SPLIT_FOR, 0, n_btsp-1, commands=[$
            'STAT_RED, arraydata, nas, mu, sig, chi2, lklh, AUTO_FWHM=1, BANDWIDTH=bandwidth, BIN_FACT=nbins_factor, BOOTSTRAP=bootstrap, CUBE_SIZE=cube_fine, FILE_ID=file_id, FILEPHI=filephi, NEG=neg,'+ $
            'MARGINALIZE=1, NFWHM_LEFT=nfwhm_left, NFWHM_RIGHT=nfwhm_right, NULL_LIM=null_lim, PHI0_LIM=phi0_lim, PHIE_LIM=phie_lim, DARK_RMS=dark_rms, RNDM_OPD=rndm_opd, SEED_BOOT=seed,' +$
            'SEED_OPD=best_opd_seed, SKIPEND=skipend, TRSH=occ_min, VARIABLE=variable, WAV_EFF=wav_eff,' +$      ; Input keywords
            'HIST_MODEL=hist_model_opt, CHI2CUBE=cube, LKLH_CUBE=lcube, NULL_EXPL=nullexpl, NBINS_KEPT=nbins_kept, HISTOBS=histobs, NULL_LOC=null_loc,'+$                                                                  ; Output keywords
            'INFO=info, PLOT=plot, FILE_PATH=file_path, PS=ps' + $
            '& IF N_ELEMENTS(nasopt)  EQ 0 THEN nasopt  = nas  ELSE nasopt  = [nasopt,nas]' + $
            '& IF N_ELEMENTS(muopt)   EQ 0 THEN muopt   = mu   ELSE muopt   = [muopt,mu]' + $
            '& IF N_ELEMENTS(sigopt)  EQ 0 THEN sigopt  = sig  ELSE sigopt  = [sigopt,sig]' + $
            '& IF N_ELEMENTS(lklhopt) EQ 0 THEN lklhopt = lklh ELSE lklhopt = [lklhopt,lklh]'+ $
            '& IF N_ELEMENTS(chi2opt) EQ 0 THEN chi2opt = chi2 ELSE chi2opt = [chi2opt,chi2]'], $
            struct2pass1=arraydata,$
            varnames=['nbins_factor','file_path','file_id','bootstrap','neg','cube_fine','occ_min','variable','wav_eff','rndm_opd','null_lim','phi0_lim','phie_lim','best_opd_seed'],$
            outvar=['nasopt','muopt','sigopt','chi2opt','lklhopt'], NSPLIT=n_cpu, /SILENT
  ; Print status to screenprint,
  IF info GT 0 THEN PRINT, 'Multi-threading done'
  ; Parse results
  nas_boot  = [nasopt0]
  mu_boot   = [muopt0]
  sig_boot  = [sigopt0]
  chi2_boot = [chi2opt0]
  lklh_boot = [lklhopt0]
  i         = 1
  WHILE EXECUTE('tmp = nasopt' + STRING(i, FORMAT='(I0)')) DO BEGIN
    IF EXECUTE('tmp = nasopt'  + STRING(i, FORMAT='(I0)')) THEN nas_boot = [nas_boot,tmp]
    IF EXECUTE('tmp = muopt'   + STRING(i, FORMAT='(I0)')) THEN mu_boot = [mu_boot,tmp]
    IF EXECUTE('tmp = sigopt'  + STRING(i, FORMAT='(I0)')) THEN sig_boot = [sig_boot,tmp]
    IF EXECUTE('tmp = lklhopt' + STRING(i, FORMAT='(I0)')) THEN lklh_boot = [lklh_boot,tmp]
    IF EXECUTE('tmp = chi2opt' + STRING(i++, FORMAT='(I0)')) THEN chi2_boot = [chi2_boot,tmp]
  ENDWHILE
  
  ; Extract best-fit mean values  
  AVGSDV, nas_boot, null_mean, null_mean_err, tmp, KAPPA=3
  AVGSDV, mu_boot, mu_mean, mu_mean_err, tmp, KAPPA=3
  AVGSDV, sig_boot, sig_mean, sig_mean_err, tmp, KAPPA=3
  AVGSDV, chi2_boot, chi2_mean, chi2_mean_err, tmp, KAPPA=3
  AVGSDV, lklh_boot, lklh_mean, lklh_mean_err, tmp, KAPPA=3
      
  ; Extract best-fit mode values
  MODSDV, nas_boot, null_mode, null_mode_err_low, null_mode_err_sup  
  MODSDV, mu_boot, mu_mode, mu_mode_err_low, mu_mode_err_sup         
  MODSDV, sig_boot, sig_mode, sig_mode_err_low, sig_mode_err_sup     
  MODSDV, chi2_boot, chi2_mode, chi2_mode_err_low, chi2_mode_err_sup
ENDIF

; Print info to screen
IF info GT 0 THEN PRINT, 'Time for NSC reduction of this OB [s]', SYSTIME(1)-t0
       
; Print info to screen
IF info GT 0 THEN BEGIN
  PRINT, ' ' 
  PRINT, ' COARSE RUN '
  PRINT, ' '
  PRINT, 'Parameters'
  PRINT, ' - null range [%]             :', STRING(1D2*null_lim_c[0], FORMAT='(F6.2)'), STRING(1D2*null_lim_c[1], FORMAT='(F6.2)')
  PRINT, ' - mu range [rad]             :', STRING(phi0_lim_c[0], FORMAT='(F6.2)'), STRING(phi0_lim_c[1], FORMAT='(F6.2)')
  PRINT, ' - phi RMS range [rad]        :', STRING(phie_lim_c[0], FORMAT='(F6.2)'), STRING(phie_lim_c[1], FORMAT='(F6.2)')
  PRINT, ' - size of coarse cube        :', cube_coarse
  PRINT, ' '
  PRINT, 'Results chi2 (minimizatiom)'
  PRINT, ' - best-fit Nas [%]           :', STRING(1D2*param_coarse.nas.chi2.opt, FORMAT='(F6.3)') 
  PRINT, ' - best-fit phase [rad]       :', STRING(param_coarse.mu.chi2.opt, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase error [rad] :', STRING(param_coarse.sig.chi2.opt, FORMAT='(F6.3)')
  PRINT, ' - best-fit background factor :', STRING(param_coarse.kdark.chi2.opt, FORMAT='(F6.3)')
  PRINT, ' - reduced chi2               :', STRING(param_coarse.chi2, FORMAT='(F6.3)')
  PRINT, ' '
  PRINT, 'Results likelihood (maximization)'
  PRINT, ' - best-fit Nas [%]           :', STRING(1D2*param_coarse.nas.lklh.opt, FORMAT='(F6.3)') 
  PRINT, ' - best-fit phase [rad]       :', STRING(param_coarse.mu.lklh.opt, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase error [rad] :', STRING(param_coarse.sig.lklh.opt, FORMAT='(F6.3)')
  PRINT, ' - best-fit background factor :', STRING(param_coarse.kdark.lklh.opt, FORMAT='(F6.3)')
  PRINT, ' '
  PRINT, 'Results chi2 (Bayesian approach)'
  PRINT, ' - best-fit Nas [%]           :', STRING(1D2*param_coarse.nas.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(1D2*param_coarse.nas.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(1D2*param_coarse.nas.chi2.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase [rad]       :', STRING(param_coarse.mu.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.mu.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.mu.chi2.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase error [rad] :', STRING(param_coarse.sig.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.sig.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.sig.chi2.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit background factor :', STRING(param_coarse.kdark.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.kdark.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.kdark.chi2.err_low, FORMAT='(F6.3)')
  PRINT, ' '
  PRINT, 'Results likelihood (Bayesian approach)'
  PRINT, ' - best-fit Nas [%]           :', STRING(1D2*param_coarse.nas.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(1D2*param_coarse.nas.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(1D2*param_coarse.nas.lklh.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase [rad]       :', STRING(param_coarse.mu.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.mu.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.mu.lklh.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase error [rad] :', STRING(param_coarse.sig.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.sig.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.sig.lklh.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit background factor :', STRING(param_coarse.kdark.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.kdark.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.kdark.lklh.err_low, FORMAT='(F6.3)')
  PRINT, ' '
  PRINT, ' FINE RUN '
  PRINT, ' '
  PRINT, 'Parameters'
  PRINT, ' - maximum null precision [%] :', STRING(1D2*null_pre, FORMAT='(F6.3)')
  PRINT, ' - null range [%]             :', STRING(1D2*null_lim[0], FORMAT='(F6.2)'), STRING(1D2*null_lim[1], FORMAT='(F6.2)')
  PRINT, ' - mu range [rad]             :', STRING(phi0_lim[0], FORMAT='(F6.2)'), STRING(phi0_lim[1], FORMAT='(F6.2)')
  PRINT, ' - phi RMS range [rad]        :', STRING(phie_lim[0], FORMAT='(F6.2)'), STRING(phie_lim[1], FORMAT='(F6.2)')
  PRINT, ' - size of fine cube          :', cube_fine
  PRINT, ' '
  PRINT, 'Results chi2 (minimizatiom)'
  PRINT, ' - best-fit Nas [%]           :', STRING(1D2*param_fine.nas.chi2.opt, FORMAT='(F6.3)') 
  PRINT, ' - best-fit phase [rad]       :', STRING(param_fine.mu.chi2.opt, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase error [rad] :', STRING(param_fine.sig.chi2.opt, FORMAT='(F6.3)')
  PRINT, ' - best-fit background factor :', STRING(param_fine.kdark.chi2.opt, FORMAT='(F6.3)')
  PRINT, ' - reduced chi2               :', STRING(param_fine.chi2, FORMAT='(F6.3)')
  PRINT, ' '
  PRINT, 'Results likelihood (maximization)'
  PRINT, ' - best-fit Nas [%]           :', STRING(1D2*param_fine.nas.lklh.opt, FORMAT='(F6.3)') 
  PRINT, ' - best-fit phase [rad]       :', STRING(param_fine.mu.lklh.opt, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase error [rad] :', STRING(param_fine.sig.lklh.opt, FORMAT='(F6.3)')
  PRINT, ' - best-fit background factor :', STRING(param_fine.kdark.lklh.opt, FORMAT='(F6.3)')
  PRINT, ' '
  PRINT, 'Results chi2 (Bayesian approach)'
  PRINT, ' - best-fit Nas [%]           :', STRING(1D2*param_fine.nas.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(1D2*param_fine.nas.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(1D2*param_fine.nas.chi2.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase [rad]       :', STRING(param_fine.mu.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.mu.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.mu.chi2.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase error [rad] :', STRING(param_fine.sig.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.sig.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.sig.chi2.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit background factor :', STRING(param_fine.kdark.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.kdark.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.kdark.chi2.err_low, FORMAT='(F6.3)')
  PRINT, ' '
  PRINT, 'Results likelihood (Bayesian approach)'
  PRINT, ' - best-fit Nas [%]           :', STRING(1D2*param_fine.nas.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(1D2*param_fine.nas.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(1D2*param_fine.nas.lklh.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase [rad]       :', STRING(param_fine.mu.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.mu.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.mu.lklh.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit phase error [rad] :', STRING(param_fine.sig.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.sig.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.sig.lklh.err_low, FORMAT='(F6.3)')
  PRINT, ' - best-fit background factor :', STRING(param_fine.kdark.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.kdark.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.kdark.lklh.err_low, FORMAT='(F6.3)')
  PRINT, ' '
  IF n_btsp GT 1 THEN BEGIN
    PRINT, 'Results from bootstrap (mean)'
    PRINT, ' - best-fit Nas               :', STRING(1D2*null_mean, FORMAT='(F6.3)'),' +/- ',STRING(1D2*null_mean_err, FORMAT='(F6.3)')
    PRINT, ' - best-fit phase [rad]       :', STRING(mu_mean, FORMAT='(F6.3)'),' +/- ',STRING(mu_mean_err, FORMAT='(F6.3)')
    PRINT, ' - best-fit phase error [rad] :', STRING(sig_mean, FORMAT='(F6.3)'),' +/- ',STRING(sig_mean_err, FORMAT='(F6.3)')
    PRINT, ' - reduced chi2               :', STRING(chi2_mean, FORMAT='(F6.3)'),' +/- ',STRING(chi2_mean_err, FORMAT='(F6.3)')
    PRINT, 'Results from bootstrap (mode)'
    PRINT, ' - best-fit Nas               :', STRING(1D2*null_mode, FORMAT='(F6.3)'),' + ',STRING(1D2*null_mode_err_sup, FORMAT='(F6.3)') ,' - ',STRING(1D2*null_mode_err_low, FORMAT='(F6.3)')
    PRINT, ' - best-fit phase [rad]       :', STRING(mu_mode, FORMAT='(F6.3)'),' + ',STRING(mu_mode_err_sup, FORMAT='(F6.3)'), ' - ',STRING(mu_mode_err_low, FORMAT='(F6.3)')
    PRINT, ' - best-fit phase error [rad] :', STRING(sig_mode, FORMAT='(F6.3)'),' + ',STRING(sig_mode_err_sup, FORMAT='(F6.3)'),' - ',STRING(sig_mode_err_low, FORMAT='(F6.3)')
    PRINT, ' - reduced chi2               :', STRING(chi2_mode, FORMAT='(F6.3)'),' + ',STRING(chi2_mode_err_sup, FORMAT='(F6.3)'),' - ',STRING(chi2_mode_err_low, FORMAT='(F6.3)')
  ENDIF
  ; Now print to LOG file
  PRINTF, lun, ' ' 
  PRINTF, lun, ' COARSE RUN '
  PRINTF, lun, ' '
  PRINTF, lun, 'Parameters'
  PRINTF, lun, ' - null range [%]             :', STRING(1D2*null_lim_c[0], FORMAT='(F6.2)'), STRING(1D2*null_lim_c[1], FORMAT='(F6.2)')
  PRINTF, lun, ' - mu range [rad]             :', STRING(phi0_lim_c[0], FORMAT='(F6.2)'), STRING(phi0_lim_c[1], FORMAT='(F6.2)')
  PRINTF, lun, ' - phi RMS range [rad]        :', STRING(phie_lim_c[0], FORMAT='(F6.2)'), STRING(phie_lim_c[1], FORMAT='(F6.2)')
  PRINTF, lun, ' - size of coarse cube        :', cube_coarse
  PRINTF, lun, ' '
  PRINTF, lun,  'Results chi2 (minimizatiom)'
  PRINTF, lun,  ' - best-fit Nas [%]           :', STRING(1D2*param_coarse.nas.chi2.opt, FORMAT='(F6.3)') 
  PRINTF, lun,  ' - best-fit phase [rad]       :', STRING(param_coarse.mu.chi2.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase error [rad] :', STRING(param_coarse.sig.chi2.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit background factor :', STRING(param_coarse.kdark.chi2.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' - reduced chi2               :', STRING(param_coarse.chi2, FORMAT='(F6.3)')
  PRINTF, lun,  ' '
  PRINTF, lun,  'Results likelihood (maximization)'
  PRINTF, lun,  ' - best-fit Nas [%]           :', STRING(1D2*param_coarse.nas.lklh.opt, FORMAT='(F6.3)') 
  PRINTF, lun,  ' - best-fit phase [rad]       :', STRING(param_coarse.mu.lklh.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase error [rad] :', STRING(param_coarse.sig.lklh.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit background factor :', STRING(param_coarse.kdark.lklh.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' '
  PRINTF, lun,  'Results chi2 (Bayesian approach)'
  PRINTF, lun,  ' - best-fit Nas [%]           :', STRING(1D2*param_coarse.nas.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(1D2*param_coarse.nas.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(1D2*param_coarse.nas.chi2.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase [rad]       :', STRING(param_coarse.mu.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.mu.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.mu.chi2.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase error [rad] :', STRING(param_coarse.sig.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.sig.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.sig.chi2.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit background factor :', STRING(param_coarse.kdark.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.kdark.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.kdark.chi2.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' '
  PRINTF, lun,  'Results likelihood (Bayesian approach)'
  PRINTF, lun,  ' - best-fit Nas [%]           :', STRING(1D2*param_coarse.nas.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(1D2*param_coarse.nas.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(1D2*param_coarse.nas.lklh.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase [rad]       :', STRING(param_coarse.mu.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.mu.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.mu.lklh.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase error [rad] :', STRING(param_coarse.sig.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.sig.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.sig.lklh.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit background factor :', STRING(param_coarse.kdark.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_coarse.kdark.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_coarse.kdark.lklh.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' '
  PRINTF, lun,  ' FINE RUN '
  PRINTF, lun,  ' '
  PRINTF, lun,  'Parameters'
  PRINTF, lun,  ' - maximum null precision [%] :', STRING(1D2*null_pre, FORMAT='(F6.3)')
  PRINTF, lun,  ' - null range [%]             :', STRING(1D2*null_lim[0], FORMAT='(F6.2)'), STRING(1D2*null_lim[1], FORMAT='(F6.2)')
  PRINTF, lun,  ' - mu range [rad]             :', STRING(phi0_lim[0], FORMAT='(F6.2)'), STRING(phi0_lim[1], FORMAT='(F6.2)')
  PRINTF, lun,  ' - phi RMS range [rad]        :', STRING(phie_lim[0], FORMAT='(F6.2)'), STRING(phie_lim[1], FORMAT='(F6.2)')
  PRINTF, lun,  ' - size of fine cube          :', cube_fine
  PRINTF, lun,  ' '
  PRINTF, lun,  'Results chi2 (minimizatiom)'
  PRINTF, lun,  ' - best-fit Nas [%]           :', STRING(1D2*param_fine.nas.chi2.opt, FORMAT='(F6.3)') 
  PRINTF, lun,  ' - best-fit phase [rad]       :', STRING(param_fine.mu.chi2.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase error [rad] :', STRING(param_fine.sig.chi2.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit background factor :', STRING(param_fine.kdark.chi2.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' - reduced chi2               :', STRING(param_fine.chi2, FORMAT='(F6.3)')
  PRINTF, lun,  ' '
  PRINTF, lun,  'Results likelihood (maximization)'
  PRINTF, lun,  ' - best-fit Nas [%]           :', STRING(1D2*param_fine.nas.lklh.opt, FORMAT='(F6.3)') 
  PRINTF, lun,  ' - best-fit phase [rad]       :', STRING(param_fine.mu.lklh.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase error [rad] :', STRING(param_fine.sig.lklh.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit background factor :', STRING(param_fine.kdark.lklh.opt, FORMAT='(F6.3)')
  PRINTF, lun,  ' '
  PRINTF, lun,  'Results chi2 (Bayesian approach)'
  PRINTF, lun,  ' - best-fit Nas [%]           :', STRING(1D2*param_fine.nas.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(1D2*param_fine.nas.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(1D2*param_fine.nas.chi2.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase [rad]       :', STRING(param_fine.mu.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.mu.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.mu.chi2.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase error [rad] :', STRING(param_fine.sig.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.sig.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.sig.chi2.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit background factor :', STRING(param_fine.kdark.chi2.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.kdark.chi2.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.kdark.chi2.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' '
  PRINTF, lun,  'Results likelihood (Bayesian approach)'
  PRINTF, lun,  ' - best-fit Nas [%]           :', STRING(1D2*param_fine.nas.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(1D2*param_fine.nas.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(1D2*param_fine.nas.lklh.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase [rad]       :', STRING(param_fine.mu.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.mu.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.mu.lklh.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit phase error [rad] :', STRING(param_fine.sig.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.sig.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.sig.lklh.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' - best-fit background factor :', STRING(param_fine.kdark.lklh.bay, FORMAT='(F6.3)'),' + ',STRING(param_fine.kdark.lklh.err_sup, FORMAT='(F6.3)') ,' - ',STRING(param_fine.kdark.lklh.err_low, FORMAT='(F6.3)')
  PRINTF, lun,  ' '
  IF n_btsp GT 1 THEN BEGIN
    PRINTF, lun, 'Results from bootstrap (mean)'
    PRINTF, lun, ' - best-fit Nas               :', STRING(1D2*null_mean, FORMAT='(F6.3)'),' +/- ',STRING(1D2*null_mean_err, FORMAT='(F6.3)')
    PRINTF, lun, ' - best-fit phase [rad]       :', STRING(mu_mean, FORMAT='(F6.3)'),' +/- ',STRING(mu_mean_err, FORMAT='(F6.3)')
    PRINTF, lun, ' - best-fit phase error [rad] :', STRING(sig_mean, FORMAT='(F6.3)'),' +/- ',STRING(sig_mean_err, FORMAT='(F6.3)')
    PRINTF, lun, ' - reduced chi2               :', STRING(chi2_mean, FORMAT='(F6.3)'),' +/- ',STRING(chi2_mean_err, FORMAT='(F6.3)')
    PRINTF, lun, 'Results from bootstrap (mode)'
    PRINTF, lun, ' - best-fit Nas               :', STRING(1D2*null_mode, FORMAT='(F6.3)'),' + ',STRING(1D2*null_mode_err_sup, FORMAT='(F6.3)') ,' - ',STRING(1D2*null_mode_err_low, FORMAT='(F6.3)')
    PRINTF, lun, ' - best-fit phase [rad]       :', STRING(mu_mode, FORMAT='(F6.3)'),' + ',STRING(mu_mode_err_sup, FORMAT='(F6.3)'), ' - ',STRING(mu_mode_err_low, FORMAT='(F6.3)')
    PRINTF, lun, ' - best-fit phase error [rad] :', STRING(sig_mode, FORMAT='(F6.3)'),' + ',STRING(sig_mode_err_sup, FORMAT='(F6.3)'),' - ',STRING(sig_mode_err_low, FORMAT='(F6.3)')
    PRINTF, lun, ' - reduced chi2               :', STRING(chi2_mode, FORMAT='(F6.3)'),' + ',STRING(chi2_mode_err_sup, FORMAT='(F6.3)'),' - ',STRING(chi2_mode_err_low, FORMAT='(F6.3)')
  ENDIF
ENDIF

; Close log file
CLOSE, lun
FREE_LUN, lun

; Arrange output data structure (LBTI_FLX2NULL format is different)
null_opt = {NAS: DOUBLE(param_fine.nas.lklh.opt), NAS_ERR_LOW: 0D, NAS_ERR_SUP: 0D, MU: FLOAT(param_fine.mu.lklh.opt), MU_ERR_LOW: 0., MU_ERR_SUP: 0., SIG: FLOAT(param_fine.sig.lklh.opt), $
            SIG_ERR_LOW: 0., SIG_ERR_SUP: 0., KDRK: FLOAT(param_coarse.kdark.lklh.opt), KDRK_ERR_LOW: 0., KDRK_ERR_SUP: 0., CHI2: 0., CHI2_ERR_LOW: 0., CHI2_ERR_SUP: 0.}

null_bay = {NAS: DOUBLE(param_fine.nas.lklh.bay), NAS_ERR_LOW: DOUBLE(param_fine.nas.lklh.err_low), NAS_ERR_SUP: DOUBLE(param_fine.nas.lklh.err_sup), MU: FLOAT(param_fine.mu.lklh.bay), MU_ERR_LOW: FLOAT(param_fine.mu.lklh.err_low), MU_ERR_SUP: FLOAT(param_fine.mu.lklh.err_sup), $
            SIG: FLOAT(param_fine.sig.lklh.bay), SIG_ERR_LOW: FLOAT(param_fine.sig.lklh.err_low), SIG_ERR_SUP: FLOAT(param_fine.sig.lklh.err_sup), KDRK: FLOAT(param_fine.kdark.lklh.bay), KDRK_ERR_LOW: FLOAT(param_fine.kdark.lklh.err_low),$
            KDRK_ERR_SUP: FLOAT(param_fine.kdark.lklh.err_sup), CHI2: FLOAT(param_fine.chi2), CHI2_ERR_LOW: 0., CHI2_ERR_SUP: 0.}

IF n_btsp GT 1 THEN BEGIN         
  null_avg = {NAS: DOUBLE(null_mean), NAS_ERR_LOW: DOUBLE(null_mean_err), NAS_ERR_SUP: DOUBLE(null_mean_err), MU: FLOAT(mu_mean), MU_ERR_LOW: FLOAT(mu_mean_err), MU_ERR_SUP: FLOAT(mu_mean_err), SIG: FLOAT(sig_mean), $
              SIG_ERR_LOW: FLOAT(sig_mean_err), SIG_ERR_SUP: FLOAT(sig_mean_err), KDRK: 0., KDRK_ERR_LOW: 0., KDRK_ERR_SUP: 0., CHI2: FLOAT(chi2_mean), CHI2_ERR_LOW: FLOAT(chi2_mean_err), CHI2_ERR_SUP: FLOAT(chi2_mean_err)}
    
  null_mod = {NAS: DOUBLE(null_mode), NAS_ERR_LOW: DOUBLE(null_mode_err_low), NAS_ERR_SUP: DOUBLE(null_mode_err_sup), MU: FLOAT(mu_mode), MU_ERR_LOW: FLOAT(mu_mode_err_low), MU_ERR_SUP: FLOAT(mu_mode_err_sup), SIG: FLOAT(sig_mode), $
              SIG_ERR_LOW: FLOAT(sig_mode_err_low), SIG_ERR_SUP: FLOAT(sig_mode_err_sup), KDRK: 0., KDRK_ERR_LOW: 0., KDRK_ERR_SUP: 0., CHI2: FLOAT(chi2_mode), CHI2_ERR_LOW: FLOAT(chi2_mode_err_low), CHI2_ERR_SUP: FLOAT(chi2_mode_err_sup)}
ENDIF ELSE BEGIN
  null_avg = null_opt
  null_mod = null_opt
ENDELSE

; Prepare output data structure
data_out = {NULL_OPT: null_opt, NULL_BAY: null_bay, NULL_AVG: null_avg, NULL_MOD: null_mod}

; Plot distributions
IF KEYWORD_SET(PLOT) THEN BEGIN
  ; Prepare plots
  plotname = file_path + 'bootstrap_ID' + STRING(file_ID, FORMAT='(I0)') + '_NBIN-FAC' + STRING(nbins_factor, FORMAT='(F3.1)')
  fit = 20./1720.
  thick  = 5.0
  xthick = 5.0
  ythick = xthick
  cthick = 3.5
  csize  = 1.3
  
  ; Plot OPD seed samples
  ; *********************

  IF n_btsp GT 1 THEN BEGIN
    ; Prepare plots
    x  = INDGEN(N_ELEMENTS(bestnull)) ; Bootstrap sample ID
    ; 1. Plot distribution of astrophysical nulls
    PLOTALL, x, 1D2*bestnull, 0, NAME=plotname, TAG='OPD_NULL', XTITLE='OPD seed samples', YTITLE='Astrophysical null [%]', TITLE='', EPS=eps, /NO_FFT, /KERNEL
    ; 2. Plot distribution of best-fit chi2
    PLOTALL, x, bestchi2, 0, NAME=plotname, TAG='OPD_CHI2', XTITLE='OPD seed samples', YTITLE='Best-fit chi2', TITLE='', EPS=eps, /NO_FFT, /KERNEL
    ; 3. Plot distribution of best-fit phase
    PLOTALL, x, mu_phase, sig_phase, NAME=plotname, TAG='OPD_MUPHI', XTITLE='OPD seed samples', YTITLE='Best-fit mean phase [rad]', TITLE='', EPS=eps, /NO_FFT, /KERNEL
    ; 3. Plot distribution of best-fit phase
    PLOTALL, x, sig_phase, 0, NAME=plotname, TAG='OPD_SIGPHI', XTITLE='OPD seed samples', YTITLE='Best-fit phase jitter [rad]', TITLE='', EPS=eps, /NO_FFT, /KERNEL
  ENDIF   

  ; Plot bootstrap samples
  ; **********************
  
  IF n_btsp GT 1 THEN BEGIN
    ; Prepare plots
    x  = INDGEN(N_ELEMENTS(nas_boot)) ; Bootstrap sample ID
    ; 1. Plot distribution of astrophysical nulls
    PLOTALL, x, 1D2*nas_boot, 0, NAME=plotname, TAG='BOOT_NULL', XTITLE='Bootstrap sample', YTITLE='Astrophysical null [%]', TITLE='', EPS=eps, /NO_FFT, /KERNEL
    ; 2. Plot distribution of best-fit chi2
    PLOTALL, x, chi2_boot, 0, NAME=plotname, TAG='BOOT_CHI2', XTITLE='Bootstrap sample', YTITLE='Best-fit chi2', TITLE='', EPS=eps, /NO_FFT, /KERNEL
    ; 3. Plot distribution of best-fit phase
    PLOTALL, x, mu_boot, sig_phase, NAME=plotname, TAG='BOOT_MUPHI', XTITLE='Bootstrap sample', YTITLE='Best-fit mean phase [rad]', TITLE='', EPS=eps, /NO_FFT, /KERNEL
    ; 3. Plot distribution of best-fit phase
    PLOTALL, x, sig_boot, 0, NAME=plotname, TAG='BOOT_SIGPHI', XTITLE='Bootstrap sample', YTITLE='Best-fit phase jitter [rad]', TITLE='', EPS=eps, /NO_FFT, /KERNEL
  ENDIF ;ELSE PLOTALL, 1D2*null_expl, napdf, 0, NAME=plotname, TAG='PDF', XTITLE='Astrophysical null [%]', YTITLE='Proba Density', TITLE='Marginalized PDF', EPS=eps, /NO_FFT 
ENDIF

RETURN, data_out
END

; --------------
PRO TEST_CALLNSC1
  ; This procedures is used to check the behavior of CALL_NSC.pro regarding error bars
  ; Paramaters
  i_aper   = 2
  n_btstrp = 100
  n_run    = 5        ; Number of reduction (each time with a different OPD seed)
  lam_cen  = 1.1D-5
  bdwdth   = 2.6D-6
  cube_size = [80,80,80]
  
  ; File with a good chi2
  file_path  = '/Volumes/nodrs/nodrs/results/nomic/l1_fits/2015-02-08/'
  ob_id      = 0;27
  null_file  = FILE_SEARCH(file_path, '*ID' + STRING(ob_id, FORMAT='(I03)') + '*NULL.fits')
  phot1_file = FILE_SEARCH(file_path, '*ID' + STRING(ob_id, FORMAT='(I03)') + '*PHOT1.fits')
  phot2_file = FILE_SEARCH(file_path, '*ID' + STRING(ob_id, FORMAT='(I03)') + '*PHOT2.fits')
  bckg_file  = FILE_SEARCH(file_path, '*ID' + STRING(ob_id, FORMAT='(I03)') + '*BCKG.fits')
  
  ; Read files
  data_null  = MRDFITS(null_file, 1, /SILENT)
  data_phot1 = MRDFITS(phot1_file, 1, /SILENT) 
  data_phot2 = MRDFITS(phot2_file, 1, /SILENT)
  data_bckg  = MRDFITS(bckg_file, 1, /SILENT)
  
  ; Compute photometry
  phot_tot_phot1 = REFORM(data_phot1.flx_tot[i_aper]) & phot_err_phot1 = REFORM(data_phot1.flx_err[i_aper])
  phot_tot_phot2 = REFORM(data_phot2.flx_tot[i_aper]) & phot_err_phot2 = REFORM(data_phot2.flx_err[i_aper])
  phot_tot = 2*(MEAN(phot_tot_phot1)+MEAN(phot_tot_phot2)) 
  
  ; Compute mean fringe SNR and keep only the best one (low SNR usualy means bad phase setpoint)
  sig_out = 5
  min_fr  = 500
  null_lim = [0.005,0.20]
  MEANCLIP, data_null.pcmsnr, avg_snr, rms_snr, CLIPSIG=5
  idx_null = WHERE(data_null.pcmsnr GE (avg_snr-3*rms_snr), n_null)
  ; Compute mean phase noise and reject high noise
  MEANCLIP, data_null.pcphstd, avg_phstd, rms_phstd, CLIPSIG=5
  idx_null = idx_null[WHERE(data_null[idx_null].pcphstd LE (avg_phstd+3*rms_phstd), n_null)]
  ; Only keep the frame in null_lim range
  IF N_ELEMENTS(null_lim) EQ 2 THEN BEGIN
    IF null_lim[0] NE 0 OR null_lim[1] NE 0 THEN BEGIN
      flx_tmp  = data_null[idx_null].flx_tot[i_aper]
      idx_null = idx_null[WHERE(flx_tmp/phot_tot GE null_lim[0] AND flx_tmp/phot_tot LE null_lim[1], n_null, /NULL)]
    ENDIF
  ENDIF
  ; Remove frames more than 5 sigma away from the median
  flx_tmp = data_null[idx_null].flx_tot[i_aper]
  AVGSDV, flx_tmp, avg_tmp, rms_tmp, rms_tmp_m
  idx_null = idx_null[WHERE(flx_tmp LT MEDIAN(flx_tmp) + 5.*rms_tmp, n_null, /NULL)]
  ; Use background in surrounding regions to remove outliers
  flx_tmp = data_null[idx_null].bck_tot[i_aper]
  AVGSDV, flx_tmp, avg_tmp, rms_tmp, rms_tmp_m
  idx_null = idx_null[WHERE(ABS(flx_tmp-avg_tmp) LE sig_out*rms_tmp)]
  ; Extract data
  data_null = data_null[idx_null]
    
  ; Extract data
  null_tot_phot = REFORM(data_null.flx_tot[i_aper])  & null_err_phot  = REFORM(data_null.flx_err[i_aper])   
  bck_tot_phot  = REFORM(data_bckg.flx_tot[i_aper])  & bck_err_phot   = REFORM(data_bckg.flx_err[i_aper])  
  
  ; Remove mean background
  AVGSDV, bck_tot_phot, avg_bck, rms_bck, KAPPA=5
  bck_tot_phot -= avg_bck
  
  phi_rms_dit = REFORM(data_null.pcphstd)*!Dpi/180.*2.2D-6/lam_cen;/1.6 
    
  ; Prepare and run NSC reduction
  chi2      = FLTARR(n_run)
  nas       = chi2
  err       = chi2
  arraydata = {DARK2: bck_tot_phot, A: phot_tot_phot1, B: phot_tot_phot2, AB: null_tot_phot, Npts: N_ELEMENTS(null_tot_phot), PHIT: phi_rms_dit^2}
  FOR i=0, n_run-1 DO BEGIN
    null_star = CALL_NSC(arraydata, CUBE_SIZE=cube_size, N_BOOTSTRAP=n_btstrp, FILE_ID=ob_id, WAV_EFF=lam_cen, BANDWIDTH=bdwdth, SEED_OPD=i+1)
    chi2[i]   = null_star.chi2
    nas[i]    = null_star.nas
    err[i]    = null_star.nas_err
    PRINT, i, nas[i], err[i], chi2[i]
  ENDFOR
  
  ; Comoute results dispersion
  AVGSDV, nas, avg_nas, rms_nas, rms_mean
  AVGSDV, err, avg_err, rms_err, rms_mean
  PRINT, 'Nas and error on nas stability: ', rms_nas, rms_err 
    
  ; Plot results
  err *= 1D+2
  nas *= 1D+2
  LOADCT, 0, /SILENT
  chi2_range = MAX(chi2)-MIN(chi2)
  err_range  = MAX(err)-MIN(err)
  xrange     = [MIN(chi2)-0.1*chi2_range,MAX(chi2)+0.1*chi2_range]
  yrange     = [MIN(err)-0.1*err_range,MAX(err)+0.1*err_range]
  PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*20./1720., FILENAME = 'NSC_err_vs_chi2_SEED-OPD.eps' 
  PLOTXY, chi2, err, YLOG=ylog, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE=title, XTITLE='Chi2', YTITLE='Error on Nas [%]', GRID=0, $
          XSTYLE=1, YSTYLE=1, XTHICK=4, YTHICK=4, THICK=4, CHARTHICK=1.2, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*20./1720.
  LOADCT, 13, /SILENT
  PLOTXY, chi2, err, /ADD, SYMBOL=2, COLOR=90, THICK=4
  LOADCT, 0, /SILENT
  PLOTXY, xrange, [MEAN(err),MEAN(err)], /ADD, LINESTYLE=1, COLOR=0, THICK=4
  PLOTXY, /FIN
  LOADCT, 0, /SILENT

  nas_range  = MAX(nas)-MIN(nas)
  yrange     = [MIN(nas)-0.1*nas_range,MAX(nas)+0.1*nas_range]
  PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*20./1720., FILENAME = 'NSC_null_vs_chi2_SEED-OPD.eps'
  PLOTXY, chi2, nas, YLOG=ylog, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE=title, XTITLE='Chi2', YTITLE='Nas [%]', GRID=0, $
    XSTYLE=1, YSTYLE=1, XTHICK=4, YTHICK=4, THICK=4, CHARTHICK=1.2, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*20./1720.
  LOADCT, 13, /SILENT
  PLOTXY, chi2, nas, /ADD, SYMBOL=2, COLOR=90, THICK=4
  LOADCT, 0, /SILENT
  PLOTXY, xrange, [MEAN(nas),MEAN(nas)], /ADD, LINESTYLE=1, COLOR=0, THICK=4
  PLOTXY, /FIN
END

; --------------
PRO TEST_CALLNSC2
  ; This procedures is used to check the behavior of CALL_NSC.pro regarding NBINS_FAC
  ; Paramaters
  i_aper   = 2
  n_btstrp = 1
  n_run    = 1        ; Number of reduction (each time with a different OPD seed)
  lam_cen  = 1.1D-5
  bdwdth   = 2.6D-6
  ;cube_size = [400,200,2100]
  nbins_fac = 1.0 + 0.5*FINDGEN(n_run)

  ; File with a good chi2
  file_path  = '/Volumes/nodrs/nodrs/results/nomic/l1_fits/2015-02-08/'
  ob_id      = 48
  null_file  = FILE_SEARCH(file_path, '*ID' + STRING(ob_id, FORMAT='(I0)') + '*NULL.fits')
  phot1_file = FILE_SEARCH(file_path, '*ID' + STRING(ob_id, FORMAT='(I0)') + '*PHOT1.fits')
  phot2_file = FILE_SEARCH(file_path, '*ID' + STRING(ob_id, FORMAT='(I0)') + '*PHOT2.fits')
  bckg_file  = FILE_SEARCH(file_path, '*ID' + STRING(ob_id, FORMAT='(I0)') + '*BCKG.fits')

  ; Read files
  data_null  = MRDFITS(null_file, 1, /SILENT) 
  data_phot1 = MRDFITS(phot1_file, 1, /SILENT)
  data_phot2 = MRDFITS(phot2_file, 1, /SILENT)
  data_bckg  = MRDFITS(bckg_file, 1, /SILENT)
  
  ; Extract data
  null_tot_phot  = data_null.flx_tot[i_aper,*]  & null_err_phot  = data_null.flx_err[i_aper,*]
  null_tot_bias  = data_null.bck_bias[i_aper,*] & null_err_bias  = data_null.bck_ebias[i_aper,*]
  phot_tot_phot1 = data_phot1.flx_tot[i_aper,*] & phot_err_phot1 = data_phot1.flx_err[i_aper,*]
  phot_tot_phot2 = data_phot2.flx_tot[i_aper,*] & phot_err_phot2 = data_phot2.flx_err[i_aper,*]
  bck_tot_phot   = data_bckg.flx_tot[i_aper,*]  & bck_err_phot   = data_bckg.flx_err[i_aper,*]

  ; Prepare and run NSC reduction
  chi2      = FLTARR(n_run)
  nas       = chi2
  err       = chi2
  arraydata = {DARK2: TRANSPOSE(null_tot_bias), A: phot_tot_phot1, B: phot_tot_phot2, AB: TRANSPOSE(null_tot_phot), Npts: N_ELEMENTS(null_tot_phot), PHIT: 0}
  FOR i=0, n_run-1 DO BEGIN
    null_star = CALL_NSC(arraydata, CUBE_SIZE=cube_size, N_BOOTSTRAP=n_btstrp, FILE_ID=ob_id, WAV_EFF=lam_cen, BANDWIDTH=bdwdth, NBINS_FACTOR=nbins_fac[i], OCC_MIN=occ_min, SEED_OPD=seed_opd, VARIABLE=variable, /PLOT, /EPS)
    chi2[i]   = null_star.chi2
    nas[i]    = null_star.nas
    err[i]    = null_star.nas_err
  ENDFOR

  ; Print results
  PRINT, nas*1D2 
  PRINT, err*1D2

  ; Plot results
  IF n_run GT 1 THEN BEGIN
    LOADCT, 0, /SILENT
    nbins_fac_range = MAX(nbins_fac)-MIN(nbins_fac)
    nas_range  = MAX(nas)-MIN(nas)
    xrange     = [MIN(nbins_fac)-0.1*nbins_fac_range,MAX(nbins_fac)+0.1*nbins_fac_range]
    yrange     = [MIN(nas)-0.1*nas_range,MAX(nas)+0.1*nas_range]
    PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*20./1720., FILENAME = 'NSC_err_vs_chi2_NBINS_FAC.eps'
    PLOTXY, nbins_fac, nas, YLOG=ylog, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE=title, XTITLE='nbins_fac', YTITLE='Nas', GRID=0, $
            XSTYLE=1, YSTYLE=1, XTHICK=4, YTHICK=4, THICK=4, CHARTHICK=1.2, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*20./1720.
    LOADCT, 13, /SILENT
    PLOTXY, nbins_fac, nas, /ADD, SYMBOL=2, COLOR=90, THICK=4
    ERRPLOT, nbins_fac, nas-err, nas+err, COLOR=90, THICK=4
    LOADCT, 0, /SILENT
    ;PLOTXY, xrange, [MEAN(err),MEAN(err)], /ADD, LINESTYLE=1, COLOR=0, THICK=4
    PLOTXY, /FIN
  ENDIF
END

; ---------------
PRO TEST_CALLNSC3 ; this procedure is a blind test for CALL NSC!
  ; Define parameters
  n_data   = 2D3
  phi0     = !pi
  dl       = 2.5d-6
  lam      = 1.1D-5
  dsigma   = dl/(lam^2)
  phot1    = 5500.
  phot2    = 6000.
  nas      = 0.01D-2
  opd_off  = 0.1d-6
 
  rms_bckg = 150.
  rms_phot = 30. 
  pha_trm1 = 1.
  pha_trm2 = 0.
  nas_min  = -1D-2
  phi_rms_min = 0
  cube_size   = [0,0,0,20,1]
  
  ns = 1
  rms_opd = 0.6d-6
  ;rms_opd = INDGEN(ns)*rms_opd         
  
  ; Generate background sequence
  err_nas = DBLARR(ns)
  err_mu  = err_nas
  err_sig = err_nas
  FOR i=0, ns-1 DO BEGIN
    opd            = opd_off + rms_opd[i]*RANDOMN(58, n_data)  ; seed=0 is probed
    dark           = rms_bckg*RANDOMN(seed1, n_data)
    phi_avg_dit    = DBLARR(n_data)
    phi_rms_dit    = DBLARR(n_data)
    pcphmcos       = 1 + DBLARR(n_data)
    pcphmsin       = DBLARR(n_data)
    spdthpos       = DBLARR(n_data)
    phot_tot_phot1 = phot1 + rms_phot*RANDOMN(seed2, n_data)
    phot_tot_phot2 = phot2 + rms_phot*RANDOMN(seed3, n_data)
    null_tot_phot  = dark + phot1 + phot2 + 2*SQRT(ABS(phot1*phot2)) * (1 - nas)/(1 + nas) * SIN(!dpi*dsigma*opd)/(!dpi*dsigma*opd) * (COS(phi0+2*!pi*opd/lam) * pha_trm1 - SIN(phi0+2*!pi*opd/lam) * pha_trm2)
  
    arraydata = {DARK2: dark, A: phot_tot_phot1, B: phot_tot_phot2, AB: null_tot_phot, NAS_MIN: nas_min, PHI_RMS_MIN: phi_rms_min, PHI_AVG_T: phi_avg_dit, PHI_VAR_T: phi_rms_dit^2, PCPHMCOS_T: pcphmcos, PCPHMSIN_T: pcphmsin, SPDTHPOS: spdthpos}
    
    path = '/Users/denis/IDLWorkspace84/nodrs/analysis/results/'
    data = CALL_NSC(arraydata, CUBE_SIZE=cube_size, N_BOOTSTRAP=0, NBINS_FACTOR=1, NO_MULTI=no_multi, NULL_RANGE=[-1d-2,2d-2], FILE_ID=0, WAV_EFF=lam, BANDWIDTH=dl,$
                    OCC_MIN=0, VARIABLE=0, FILE_PATH=path, /EPS, /PLOT, /INFO)
    mu  = opd_off/lam*2*!dpi
    sig = rms_opd[i]/lam*2*!dpi
    print, 'Good values'
    print, nas
    print, mu
    print, sig
    err_nas[i] = ABS(nas-data.NULL_BAY.nas)
    err_mu[i]  = ABS(mu-data.NULL_BAY.mu)
    err_sig[i] = ABS(sig-data.NULL_BAY.sig)
  ENDFOR
  PLOT, rms_opd, err_mu
  
END

