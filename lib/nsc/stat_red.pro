;+
; NAME: STAT_RED
; 
; DESCRIPTION
;   Computes the best astrophysical null (nas), mean phase (muopt) and rms phase (sigopt) and gives the corresponding chi2 and likelyhood.
;
; INPUT
;   arraydata  : a structure containing at least the four following fields
;        .dark      : the background measurements
;        .a         : the photometric signal for the first beam
;        .b         : the photometric signal for the second beam
;        .ab        : the interferometric signal
;        .phi_avg_t : mean phase over integration time T
;        .phi_var_t : phase variance over integration time T
;        .pcphmcos_t: mean COS(phase) over integration time T (superseed phi_var_t if not null)
;        .pcphmsin_t: mean SIN(phase) over integration time T (superseed phi_var_t if not null)
;        .dit_pat   : OPD dither pattern
;
; KEYWORDS
;   AUTO_FWHM  : Use the automatic cutoffs with the FWHM of the histogram
;   BANDWIDTH  : The bandwidth of the observations [m]
;   BIN_FACT   : Nominaly, if the number of recorded null points is Np, you use sqrt(Np) bins in the initial histogram,
;                which corresponds to nbins_factor=1. The number of histogram bins can in general be set to nbins_factor*sqrt(Np)
;   BOOTSTRAP  : Set this keyword if you want the analysis to based on a "boostrapped" version of the input data (random sampling with replacement, not used by default)
;   CUBE_SIZE  : 5-element vector containing the number of elements in the chi2 cube along the null, mean phase, rms phase directions, OPD see, and multiplicative factor on the abckground RMS respectively
;                If not set, the default values (200, 100, 100, 1, 1) will be used.
;   FILE_ID    : in output, names the graphics according to the star used
;   FLAG       : flag to add to plot names
;   NEG        : Allows negative values for the astrophysical null
;   NFWHM_LEFT : Defines the range of the histogram to be kept on the short (left) end, expressed in terms of FHWM (default = -1: include all data at short end)
;   NFWHM_RIGHT: Defines the range of the histogram to be kept on the long (right) end, expressed in terms of FHWM (default = -1: include all data at long end)
;   NULL_LIM   : 2-element vector containing the min and max null to be explored in the chi2 cube. If not set, will be "best guessed" by the routine.
;   PHI0_LIM   : On input, the range of mean phases to explore (in radian). Derived from the null if not set.
;   PHIE_LIM   : On input, the range of phase errors to explore (in radian). Derived from the null if not set.
;   SEED_BOOT  : Defines the seed of the random generator for the bootstrapping to be used on the observed data
;   SEED_OPD   : Defines the seed of the random OPD generator used in the simulated null sequence and histogram.
;                Typically, seed2 is undefined so that it uses the computer's system clock.
;   SKIPEND    : Set this keyword to the number of bins to skip at the end of the histogram
;   TRSH       : Minimum number of occurences per bin used in the fit.
;   VARIABLE   : Set this keyword if you want to use histrograms with variable bin size (constant bin size used by default).
;   WAV_EFF    : The effective wavelength [m]
;   INFO       : Set this keyword to print info to screen
;   PLOT       : Set this keyword to plot various results
;   FILE_PATH  : Set this keyword to the path where EPS figures will be saved (workspace directory by default)
;   PS         : Set this keyword to sent to graphics to a postscript output instead of the screen
;
; OUTPUTS
;   A structure with the following sub-structures:
;     - nas      : Best-fit and bayesian astrophysical null
;     - mu       : Best-fit and bayesian phase
;     - sig      : Best-fit and bayesian phase jitter
;     - seed_opd : Best-fit and bayesian OPD seed
;     - kdark    : Best-fit and bayesian scale facot on backgroun RMS
;     
;   Each sub-structure contains the following information:
;     - chi2     : chi2-based results
;     - lklh     : likelyhood-based results
;     
;   Each of this strucure contains the following information
;     - opt      : best-fit value 
;     - bay      : bayesian result
;     - err_low  : low error on bay 
;     - err_sup  : sup error on bay
;     
; MODIFICATION HISTORY
;   Version 1.0, 26-OCT-2012, by Lindsay Marion, ULg, marion@astro.ulg.ac.be (based on REDCHOPFN, initally developped by B. Mennesson)
;   Version 1.1, 30-OCT-2012, OA: introduced keywords NULL_LIM and CUBE_SIZE + cleaned up the code
;   Version 1.2, 06-NOV-2012, LM/OA: introduced keyword NBINS in call to HISTO_VARIABLE, added 2D plots of the chi2 cube, cleaned up code
;   Version 1.3, 21-NOV-2012, LM/OA: introduced keyword SKIPEND, correction of bootstrap
;   Version 1.4, 27-NOV-2012, LM: automatic way to select most relevant bins in histogram (constant and variable bin size), removed "threshold"
;   Version 1.5, 06-DEC-2012, LM: cleaned up the code + improved selection of most relevant bins based on FWHM of histogram
;   Version 1.6, 18-DEC-2012, OA: further cleaning of the code, re-introduced skipend to manually skip bins if necessary, merged thresholds in constant and variable cases
;   Version 1.7, 15-NOV-2013, LM: correcting bug in the histogram plot
;   Version 1.8, 22-NOV-2013, LM: add a bootstrap for the histogram and chi2cube : allows to plot a mean histogram and to have a mean chi2cube
;   Version 1.9, 20-JAN-2014, LM: add clear comments + cleaning up the code
;   Version 1.10, 21-JAN-2014,OA: final cleaning before commit to SVN
;   Version 1.11, 17-FEB-2014,LM: clean up the code and add keywords add keywords "neg" to allow searching for negative nulls and "dark_rms" to use the dark_rms as an additionnal free parameter if needed
;   Version 1.12, 3-MAR-2014, LM: adding a new use of the histogram : first the according to FWHM, then recalculate the hole histogram and fit this new resized histogram instead of keeping
;                                     the original one as in first version
;   Version 1.13, 14-MAR-2014,LM/DD: simplifying the choose of data to reduce : add keyword betauma, betaleo and alfser : just place all observations for 1 star in
;                                    one file, choose the star you want and let the computer work!
;   Version 1.14, 18-MAR-2014, LM: clean up the code : remove useless graphics
;   Version 1.15, 01-MAY-2014, LM: add keyword "phifile" : can select the phase from the file
;   Version 1.16, 06-MAY-2014, DD: wavelength and bandwidth are now passed through the WAV_EFF and BANDWIDTH keywords
;   Version 1.17, 06-MAY-2014, DD: added keywords PLOT and INFO
;   Version 1.18, 28-NOV-2014, DD: added keyword FILE_PATH and updated figure file names
;   Version 1.19, 08-JAN-2015, DD: improved keyword description
;   Version 1.20, 1&-MAR-2015, DD: shifted bin location to center of the bins for constant bin size (for display only)
;   Version 1.21, 12-MAR-2015, DD: now works for negative nulls at peak
;   Version 1.22, 13-MAR-2015, DD: added phase noise within DIT to null model (following B. Mennesson approach in reduclbti.pro) + prepared tip/tilt in model but 0 for now
;   Version 1.23, 13-APR-2015, DD: added MODE information
;   Version 1.24, 31-JUL-2015, DD: modified for Bayesian approach (added keyword RDNM_OPD) + added keywords PHI0_LIM and PHIE_LIM
;   Version 1.25, 04-AUG-2015, DD: corrected bug with nbins_factor and constant bin size (nbins_factor was omitted after the histogram resize)
;   Version 1.26, 06-AUG-2015, DD: the TRSH keyword is now used to pass an user-defined minimum of occurences per bin
;   Version 1.27, 08-AUG-2015, DD: now allow different sizes between the auxiliary data sequences (photometries and background) and the null sequence (i.e., rescale the modeled histogram).
;   Version 1.28, 10-AUG-2015, DD: now use a quadratic grid for muphi and sigphi (null ~ mu^2 and sig^2)
;   Version 1.29, 24-AUG-2015, DD: removed quadratic grid and added mean phase over DIT to null model
;   Version 1.30, 30-AUG-2015, DD: now works if cube_size contains 1(s) + allow different OPD seeds over the grid (i.e., seed_opd can be an array)
;   Version 1.31, 11-SEP-2015, DD: changed the range used when fitting the background RMS and add keyword DRK_SCALE
;   Version 1.32, 15-NOV-2016, DD: added PCPHMCOS_T and PCPHMSIN_T to input structure
;   Version 2.00, 24-NOV-2016, DD: added maximum likelyhood keyword (i.e., LKLHCUBE). Maximum likelyhood is now used by default to compute nasopt, muopt, and sigopt.
;   Version 2.01, 05-DEC-2016, DD: added keywords MU_EXPL and SIG_EXPL
;   Version 2.02, 13-DEC-2016, DD: added keyword MARGINALIZE
;   Version 2.03, 30-JAN-2017, DD: removed keywords FILE_PHI, RNDM_OPD, DARK_RMS, and DARK_SCALE + added FLAG keyword + CUBE_SIZE is now a 5 element vector + added OPD dither pattern.
;   Version 2.04, 13-FEB-2017, DD: now marginalize over OPD_SEED and background multiplifcation factor too (+ removed keyword MARGINALIZED)

PRO STAT_RED, arraydata, param_out, $
              AUTO_FWHM=auto_fwhm, BANDWIDTH=bandwidth, BIN_FACT=nbins_factor, BOOTSTRAP=bootstrap, CUBE_SIZE=cube_size, FILE_ID=file_id, NEG=neg, $
              NFWHM_LEFT=nfwhm_left, NFWHM_RIGHT=nfwhm_right, NULL_LIM=null_lim, PHI0_LIM=phi0_lim, PHIE_LIM=phie_lim, SEED_BOOT=seed_boot, $
              SEED_OPD=seed_opd, SKIPEND=skipend, TRSH=trsh, VARIABLE=variable, WAV_EFF=wav_eff, $  ; Input keywords                                                                                                                                                                                              ; Output keywords
              INFO=info, PLOT=plot, FILE_PATH=file_path, FLAG=flag, PS=ps                                                                                                                                                                                             ; Print/plot keywords

; Keyword sanity check
IF NOT KEYWORD_SET(SKIPEND)      THEN skipend      = 0
IF NOT KEYWORD_SET(NBINS_FACTOR) THEN nbins_factor = 1
IF NOT KEYWORD_SET(NFWHM_LEFT)   THEN nfwhm_left   = -1
IF NOT KEYWORD_SET(NFWHM_RIGHT)  THEN nfwhm_right  = -1
IF NOT KEYWORD_SET(FILE_ID)      THEN file_id      = 0
IF NOT KEYWORD_SET(FLAG)         THEN flag         = ''   ELSE flag = '_' + flag
IF NOT KEYWORD_SET(WAV_EFF)      THEN lam          = 11.1 ELSE lam  = wav_eff*1D6 
IF NOT KEYWORD_SET(BANDWIDTH)    THEN dl           = 2.6  ELSE dl   = bandwidth*1D6 
IF NOT KEYWORD_SET(PLOT)         THEN plot         = 0 
IF NOT KEYWORD_SET(FILE_PATH)    THEN file_path    = '' 
IF NOT KEYWORD_SET(INFO)         THEN info         = 0 
IF NOT KEYWORD_SET(TIP)          THEN tip          = 0 
IF NOT KEYWORD_SET(TILT)         THEN tilt         = 0 

; Rearrange the input data
dark       = arraydata.dark2 ;background
a          = arraydata.a     ;photometric signal "a"
b          = arraydata.b     ;photometric signal "b"
ab         = arraydata.ab    ;interferometric signal a+b
Npts       = N_ELEMENTS(ab)
spdthpos   = arraydata.spdthpos
IF MAX(arraydata.pcphmcos_T) NE 0 THEN BEGIN 
  pha_trm1 = arraydata.pcphmcos_T
  pha_trm2 = arraydata.pcphmsin_T
ENDIF ELSE BEGIN
  pha_trm1 = 1.-0.5*arraydata.phi_var_t+0.125*arraydata.phi_var_t^2
  pha_trm2 = arraydata.phi_avg_t
ENDELSE

;********************************
; Histogram with normalised null  ---  just for information, the best fit will be done with ab and not ab/peak
;********************************

;We choose a number of bins = sqrt(Npts) and multiply it by nbins_factor (more bins => best sampling)
phot_tot       = MEAN(a+b+2*SQRT(ABS(a*b))) 
null_normalise = ab/phot_tot
moy_null       = (MEAN(ab))/(MEAN(a + b + 2*SQRT(ABS(a*b))))
nbinss         = ROUND(nbins_factor * SQRT(FLOAT(Npts)))
di             = 0.5*(b-a)/(a+b) ;differential in intensity
nph            = di^2 ; estimate the contribution of intensity variations to the mesured null
histobs_full   = HISTOGRAM(null_normalise, NBINS=nbinss, LOCATIONS=null_localise, MIN=MIN(null_normalise), MAX=MAX(null_normalise))
nbins_real     = nbinss-skipend ; skip the last bins if specified
histobs0       = histobs_full[0:nbins_real-1]

; Print info to sreen
IF info GT 0 THEN BEGIN
  PRINT, ''
  PRINT,'Minimum null value observed (in counts) =',MIN(ab)
  PRINT,'Minimum dark value observed (in counts) =',MIN(dark)
  PRINT,'Quick and dirty null estimate =',(MIN(ab) - MIN(dark))/(mom_a[0] + mom_b[0] + 2*SQRT(mom_a[0]*mom_b[0]))
  mod_null = MODE_WEIGHT(ab)/(MEAN(a + b + 2*SQRT(ABS(a*b))))
  PRINT,'Mean null level=',moy_null, mod_null ,  '(mode)'
  PRINT,'Mean null level expected from photometric imbalance=',MEAN(nph)
ENDIF

; Plot histogram
IF KEYWORD_SET(PLOT) THEN BEGIN
  IF KEYWORD_SET(ps) THEN BEGIN
    PREP_PS
    DEVICE, FILENAME= file_path + 'histo_nullnorm_ID' + STRING(file_ID, FORMAT='(I0)') +'.eps', /ENCAPSULATE;, /COL, BITS_PER_PIXEL=8
  ENDIF ELSE WINDOW,2
  null_loc_shft = null_localise + 0.5*(null_localise[1]-null_localise[0]) ; Shift to bin center
  PLOT, null_loc_shft, histobs_full, title='Observed null signal histogram', xtitle='Null signal', ytitle='Nb of occurrences', psym=10;, xrange=[0.0,0.2]
  IF KEYWORD_SET(PS) THEN BEGIN
    DEVICE, /CLOSE
    END_PS
  ENDIF
ENDIF

;*******************
; Histogram with ab 
;*******************

null_measured = ab ; signal close to null (NOT normalized by peak)

; WARNING, RANDOMN function update the seed!!!
; Store seeds in separate variables
IF KEYWORD_SET(SEED_BOOT) THEN seed = seed_boot 

; Re-arrange the data vector if boostrapping is required
IF KEYWORD_SET(bootstrap) THEN BEGIN
  ind_resamp = MYRESAMPLE(Npts,seed)
  null_measured = null_measured[ind_resamp]
  spdthpos      = spdthpos[ind_resamp]
ENDIF 

; Resample auxiliary data.
;ind_resamp = MYRESAMPLE(Npts,seed2)
;;a = a[ind_resamp]
;b = b[ind_resamp]
;dark       = dark[ind_resamp]
;phi_avg_T  = phi_avg_T[ind_resamp]
;phi_var_T  = phi_var_T[ind_resamp]

; Plot the constant-bin-size histogram for the fringe signal 
histobs_full = HISTOGRAM(null_measured, NBINS=nbinss, LOCATIONS=null_localise, MIN=MIN(null_measured), MAX=MAX(null_measured)) 
histobs      = histobs_full[0:nbins_real-1]

; Plot histogram
IF KEYWORD_SET(PLOT) THEN BEGIN
  IF KEYWORD_SET(ps) THEN BEGIN
    PREP_PS
    DEVICE, FILENAME = file_path + 'histo_signal_ID'+ STRING(file_ID, FORMAT='(I0)') +'.eps', /ENCAPSULATE;, /COL, BITS_PER_PIXEL=8
  ENDIF ELSE WINDOW,4
  null_loc_shft = null_localise + 0.5*(null_localise[1]-null_localise[0]) ; Shift to bin center
  PLOT, null_loc_shft, histobs_full, title='Observed fringe signal histogram', xtitle='Interferometric signal [counts]', ytitle='Nb of occurrences', psym=10;, xrange = [0,500]
  IF KEYWORD_SET(skipend) THEN BEGIN
    LOADCT, 1, /SILENT
    OPLOT, null_loc_shft[0:nbins_real-1], histobs, psym=-4, COLOR=150 ; overplot the selected point in blue
    LOADCT, 0
  ENDIF
  IF KEYWORD_SET(PS) THEN BEGIN
    DEVICE, /CLOSE
    END_PS
  ENDIF
ENDIF

IF KEYWORD_SET(auto_fwhm) THEN BEGIN

  ;*************************************
  ; Prepare variable bin size histogram 
  ;*************************************
  
  ; after several tests, variable bin size histogram is useless
  IF KEYWORD_SET(variable) THEN BEGIN 
  
    ; SET_VARSIZE_BIN used for each bootstrap realization of the observed vector of null values.
    ; Loses a bit of computing time, but this is more exact. 
    binpos = SET_VARSIZE_BINS(null_measured, N_BIN=nbinss)
    histobs_var = HISTO_VARIABLE(null_measured, binpos, NBINS=nbinss)
    nbins_real = nbinss-skipend ; skip the last bins if specified
    histobs = histobs_var[0:nbins_real-1]
    
    ; Define the range of nulls on which the histogram will need to be fitted
    IF (nfwhm_left NE -1) OR (nfwhm_right NE -1) THEN BEGIN
    
      ; Find the histogram peak, i.e., the smallest bin
      ; Since, with a variable bin size, all the bins will be equally populated, we use the smallest bin, i.e. the one that was, the most populated before resizing to
      ;compute an artificial full width at half maximum that we will copy left and right of the smallest bin.  
      binsize = (binpos - SHIFT(binpos,1))[1:*]
      ; It may happen that the max value of the histogram is repeated many time --> last values of minibin1 are 0 
      ;--> remove to find true minibin (happens with vegap36deg2 for instance in case of PFN)
      tmp = WHERE(binsize EQ 0, COMPL=index_nonzero)
      binsize_nonzero = binsize[index_nonzero]
      binsize_smooth = SMOOTH(binsize_nonzero, 5)
      minbinsize = MIN(binsize_smooth, indtmp)
      index_peak = index_nonzero[indtmp]
      peak = (binpos[index_peak] + binpos[index_peak+1]) / 2 ; use the center of the bin with minimum size as a reference for the peak position
      nullmod = (peak)/(MEAN(a+b+2*SQRT(ABS(a*b)))) ; null normalisation
      
      ; Defining an equivalent FWHM, using the bins where size(bin) = 2*minbinsize
      tmp = MIN(ABS(index_nonzero - index_peak), index_peak_nonzero) ; define the index of the peak in the index_nonzero vector
      bin1 = MIN(ABS(binsize_smooth[0:index_peak]-2*minbinsize),index_fwhm_short) ; find the bin that is twice the minimum bin size to the left (short end)
      bin2 = MIN(ABS(binsize_smooth[index_peak:*]-2*minbinsize),index_tmp) ; find the bin that is twice the minimum bin size to the right (long end)
      index_fwhm_long = index_peak + index_tmp
      FWHME = binpos[index_fwhm_long] - binpos[index_fwhm_short]
  
      ; Print info to screen
      IF info GT 0 THEN BEGIN
        PRINT, 'Size of the smallest bin in the variable bin-size histogram: ', minbinsize
        PRINT, 'Interferometric signal (in counts) at peak of distribution: ', peak
        PRINT, 'Null level at peak of distribution: ', nullmod
        PRINT, 'Equivalent FWHM of the histogram: ', FWHME
      ENDIF
      
      ; Values of null that will be the cutoff min and max
      IF nfwhm_left NE -1 THEN BEGIN ; define the min cutoff
        IF peak - nfwhm_left*FWHME LE 0 THEN null_lim_left = binpos[0] $
        ELSE null_lim_left = peak - nfwhm_left*FWHME
        tmp1 = MIN(ABS(binpos - null_lim_left),index_min)
      ENDIF ELSE index_min = 0
      IF nfwhm_right NE -1 THEN BEGIN ; define the max cutoff
        IF peak + nfwhm_right*FWHME GE binpos[nbins_real-1] THEN null_lim_right = binpos[nbins_real-1] $
        ELSE null_lim_right = peak + nfwhm_right*FWHME
        tmp2 = MIN(ABS(binpos - null_lim_right),index_max)
      ENDIF ELSE index_max = nbins_real-1
        
    ENDIF ELSE BEGIN
      index_min = 0
      index_max = nbins_real-1
    ENDELSE
    
    ; Print info to screen
    cutoff_min = binpos[index_min] ; left cutoff in term of null
    cutoff_max = binpos[index_max] ; right cutoff in term of null
    nbins_kept = index_max-index_min
    IF info GT 0 THEN BEGIN
      PRINT, 'index_min',index_min
      PRINT, 'index_max', index_max
      PRINT, 'cutoff min', cutoff_min
      PRINT, 'cutoff max', cutoff_max
    ENDIF
  
    ; Resize histogram according to cutoffs
    histobs = histobs[index_min:index_max]
    
    ; Plot histogram
    IF KEYWORD_SET(PLOT) THEN BEGIN
      IF KEYWORD_SET(ps) THEN BEGIN
        PREP_PS
        DEVICE, FILENAME = file_path + 'histo_signal_var_'+ STRING(file_ID, FORMAT='(I0)') +'.eps', /ENCAPSULATE;, /COL, BITS_PER_PIXEL=8
      ENDIF ELSE WINDOW, 3
      PLOT, binpos, histobs_var, title='Variable bin size histogram (selected: in blue)', xtitle='Interferometric signal [counts]', ytitle='Nb of occurrences', psym=10
      LOADCT, 1, /SILENT
      OPLOT, binpos[index_min:index_max], histobs, PSYM=10, COLOR=150
      LOADCT, 0, /SILENT
      IF KEYWORD_SET(PS) THEN BEGIN
        DEVICE, /CLOSE
        END_PS
      ENDIF
    ENDIF   
    
  ;*************************************
  ; Prepare constant bin size histogram 
  ;*************************************
  
  ENDIF ELSE BEGIN
    ; Smooth the histogram => avoid abrupt fluctuations
    histosmooth = SMOOTH(histobs,3)
  
    ; Find the maximum of the histogram
    max_histo = MAX(histosmooth,indmax)
    halfmax   = max_histo/2                         ; define half of the peak
    peak      = null_localise[indmax[0]]            ; maximum position
    nullmod   = (peak)/(MEAN(a+b+2*SQRT(ABS(a*b)))) ; normalisation du null
    
    ; Print info to screen
    IF info GT 0 THEN BEGIN
      PRINT,'Interferometric signal (in counts) at peak of distribution: ',peak
      PRINT,'Null level at peak of distribution:',nullmod
    ENDIF
    
    ; Define the range of nulls on which the histogram will need to be fitted
    IF (nfwhm_left NE -1) OR (nfwhm_right NE -1) THEN BEGIN
    
      ; Find index corresponding to halfmax before peak, and then derive the corresponding null value (using linear interpolation for better accuracy)
      tmp = MIN(ABS(histosmooth[0:indmax]-halfmax), index_left1)
      IF (histosmooth[index_left1] LE halfmax) THEN BEGIN
        index_left2 = index_left1 + 1 
      ENDIF ELSE BEGIN
        IF index_left1 EQ 0 THEN BEGIN 
          testindex = 1
        ENDIF ELSE BEGIN
          index_left2 = index_left1
          index_left1 = index_left1 - 1
          ENDELSE
      ENDELSE
      ;IF index_left1 EQ -1 THEN index_left1 = 0
      IF KEYWORD_SET(testindex) THEN BEGIN
        nullhalf_left = null_localise[0]
      ENDIF ELSE BEGIN
        bound_lo = histosmooth[index_left1]
        bound_up = histosmooth[index_left2]
        IF histosmooth[index_left1] EQ histosmooth[index_left2] THEN bound_up = histosmooth[index_left2+1]
        nullhalf_left = (FLOAT(halfmax-bound_lo)) / (FLOAT(bound_up-bound_lo)) * (null_localise[index_left2]-null_localise[index_left1]) + null_localise[index_left1]
      ENDELSE
      ; Idem for halfmax after peak
      tmp = MIN(ABS(histosmooth[indmax:*]-halfmax),index_tmp)
      index_right1 = indmax + index_tmp
      IF (histosmooth[index_right1] LE halfmax) THEN index_right2 = index_right1 - 1 ELSE BEGIN
        index_right2 = index_right1
        index_right1 = index_right1 - 1
      ENDELSE
      bound_lo = histosmooth[index_right1]
      bound_up = histosmooth[index_right2]
      nullhalf_right = (FLOAT(halfmax-bound_lo)) / (FLOAT(bound_up-bound_lo)) * (null_localise[index_right2]-null_localise[index_right1]) + null_localise[index_right1]
    
      ; Derive FWHM in terms of null and related indices
      FWHM = FLOAT(nullhalf_right) - FLOAT(nullhalf_left)
      IF info GT 0 THEN PRINT, 'FWHM of the histogram: ', FWHM
      IF nfwhm_left NE -1 THEN BEGIN
        IF peak-nfwhm_left*FWHM LE 0 THEN null_lim_left = null_localise[0] ELSE null_lim_left = peak-nfwhm_left*FWHM
        tmp1 = MIN(ABS(null_localise-null_lim_left),index_min)
        nindmin = N_ELEMENTS(index_min)
        IF nindmin GT 1 THEN index_min = index_min[nindmin/2]
      ENDIF ELSE BEGIN
      index_min =0
      ENDELSE
      IF nfwhm_right NE -1 THEN BEGIN
        IF peak+nfwhm_right*FWHM GE null_localise[nbins_real-1] THEN null_lim_rig = null_localise[nbins_real-1] ELSE null_lim_rig = peak+nfwhm_right*FWHM
        tmp2 = MIN(ABS(null_localise-null_lim_rig),index_max)
        nindmax = N_ELEMENTS(index_max)
        IF nindmax GT 1 THEN index_max = index_max[nindmax/2]
      ENDIF ELSE BEGIN
        index_max = nbins_real-1
      ENDELSE
          
    ENDIF ELSE BEGIN
      index_min = 0
      index_max = nbins_real-1
    ENDELSE
    
    ; Print info to screen
    cutoff_min = null_localise[index_min] ; left cutoff in term of null
    cutoff_max = null_localise[index_max] ; right cutoff in term of null
    IF info GT 0 THEN BEGIN
      PRINT, 'index_min',index_min
      PRINT, 'index_max', index_max
      PRINT, 'cutoff min', cutoff_min
      PRINT, 'cutoff max', cutoff_max
    ENDIF
  
    ; Resize histogram according to cutoffs : keep the number of occurrences between the cutoffs and recompute an histogram between the cutoffs
    histo_resize = histobs[index_min:index_max]
    N_occ = TOTAL(histo_resize)
    nbins_kept = ROUND(nbins_factor*SQRT(FLOAT(N_occ)))
    null_measured_sort = null_measured[SORT(null_measured)]
    MIN_null = MIN(ABS(null_measured_sort - cutoff_min), index_mini)
    Max_null = MIN(ABS(null_measured_sort - cutoff_max), index_maxi)
  
    null_measured_kept = null_measured_sort[index_mini:index_maxi]
    histobs = HISTOGRAM(null_measured_kept,nbins=nbins_kept,MIN=MIN(null_measured_kept), MAX=MAX(null_measured_kept),location=null_loc)
    
    binsize_init = (Max(null_measured)-MIN(null_measured))/nbinss-1
    binsize_new = (MAX(null_measured_kept)-MIN(null_measured_kept))/nbins_kept-1
    
    ; Plot histogram
    IF KEYWORD_SET(PLOT) THEN BEGIN
      IF KEYWORD_SET(ps) THEN BEGIN
        PREP_PS
        DEVICE, FILENAME = file_path + 'histo_signal_con_ID'+ STRING(file_ID, FORMAT='(I0)') +'.eps', /ENCAPSULATE;, /COL, BITS_PER_PIXEL=8
      ENDIF ELSE WINDOW,5
      null_loca_shft = null_localise + 0.5*(null_normalise[1]-null_normalise[0])
      PLOT, null_loc_shft, histobs_full, title='Constant bin size histogram (selected: in blue)', xtitle='Interferometric signal [counts]', ytitle='Nb of occurrences', psym=10;,xrange = [-500,1500]
      LOADCT, 1, /SILENT
      OPLOT, null_loc, histobs*(binsize_init/binsize_new), PSYM=10, COLOR=150
      LOADCT, 0, /SILENT
      IF KEYWORD_SET(PS) THEN BEGIN
        DEVICE, /CLOSE
        END_PS
      ENDIF
    ENDIF  
  ENDELSE
ENDIF  ; End IF on AUTO_FWHM

; Remove bins for which the occurence is below the threshold 
IF KEYWORD_SET(trsh) THEN BEGIN
 
  treshold  = trsh
  pixel_max = MAX(histobs_full,indmax);index représente l'abscisse du pixel le plus peuplé
  peak      = null_localise[indmax[0]] ; position du maximum
  nullmod   = (peak)/(MEAN(a+b+2*SQRT(ABS(a*b))));normalisation du null
  
  ; Print info to screen
  IF info GT 0 THEN BEGIN
    PRINT,'interferometric signal (in counts) at peak of distribution=',peak
    PRINT,'null level at peak of distribution=',nullmod
  ENDIF
  
  index = WHERE(histobs_full GT treshold) 
  ; imax will contain the maximum index such that histobs0 > treshold in a "contiguous" way
  ; index>indmax - SHIFT(index>indmax,1) gives the separation of two consecutive indices (valid only beyond the peak -- 0 before the peak)
  ; We look where this expression is greater than 1, which represents a discontinuity
  ; Then, we take the first element where this condition is fullfilled, and subtract 1 to stay in the contiguous part
  IF KEYWORD_SET(trsh) THEN BEGIN
    cont = WHERE((index>indmax - SHIFT(index>indmax,1)) GT 1)
    cont = cont[0]
    IF cont LT 0 THEN imax=MAX(index) ELSE imax = MAX(index[cont-1])
  ENDIF ELSE BEGIN
    imax = index[(WHERE(index>indmax - SHIFT(index>indmax,1)) GT 1)[0]-1]
    imax = max(index)
  ENDELSE
  
  cutoff_min = null_localise[index[0]]
  cutoff_max = null_localise[imax]
  nbins_kept = imax - index[0] + 1
  
  histobs = HISTOGRAM(null_measured, nbins=nbins_kept, locations=null_localise, min=cutoff_min, max=cutoff_max)
  
  ; Plot histogram
  IF KEYWORD_SET(PLOT) THEN BEGIN
    IF KEYWORD_SET(ps) THEN BEGIN
      PREP_PS
      DEVICE, FILENAME = file_path + 'histo_signal_threshold_'+STRING(file_ID, FORMAT='(I0)')+'.eps', /ENCAPSULATE;, /COL, BITS_PER_PIXEL=8
    ENDIF ELSE WINDOW,5
    PLOT, null_localise, histobs, title='Observed signal histogram (where occur > threshold)', xtitle='Signal',ytitle='Nb of occurrences', psym=10, yrange=[0,MAX(histobs)]
    IF KEYWORD_SET(PS) THEN BEGIN
      DEVICE, /CLOSE
      END_PS
    ENDIF
  ENDIF
ENDIF

; Final number of elements in observed histogram
n_obs = TOTAL(histobs)

;***********************
; Monte-Carlo Algorithm
;***********************

;3 parameters are free : astrophysical null, mean phase and rms phase. The others are supposed to be known (a,b and dark)
a_model    = a   ; OK for PFN data. (for KI data am2, bm2 are different from a and b because of the need for interpolation etc)
b_model    = b
dark_model = dark

IF KEYWORD_SET(cube_size) THEN BEGIN
  Nf   = cube_size[0]
  Nmp  = cube_size[1]
  Nrms = cube_size[2]
  Nseed= cube_size[3]
  Ndrms= cube_size[4]
ENDIF ELSE BEGIN
  Nf    = 200 ; number of sampled values for Nas 
  Nmp   = 100 ; number of sampled values for mean phase (or  opd)
  Nrms  = 100 ; number of sampled values for phase rms over the sequence
  Nseed = 1   ; number of OPD seed
  Ndrms = 1   ; number of dark RMS multiplication factor
ENDELSE

; DEFINE PRIORS
; The contribution to the null of the phase is given by : N= delta phi^2/4 ==> delta phi =2sqrt(N) (assuming no other problems)
IF NOT KEYWORD_SET(phi0_lim) THEN BEGIN
  offset    = MEAN(arraydata.phi_var_T)/4D           ; estimated null from the high-frequency phase noise
  phi_mean0 = 6*SQRT(moy_null-offset)                ; maximum offset from pi given the observed mean null
  phi_mean  = phi_mean0*(FINDGEN(Nmp)+1)/(Nmp)       ; ONLY NEEDS TO BE BETWEEN 0 AND PI, SINCE PHI AND 2PI-PHI GIVE SAME COSINUS (avoid 0)
  ; in fact this is the mean offset from constructive, so it'd better be closer to pi.
ENDIF ELSE phi_mean = MIN(phi0_lim) + (MAX(phi0_lim)-MIN(phi0_lim))*FINDGEN(Nmp)/((Nmp-1)>1);ALOG((FINDGEN(Nmp)+1))/ALOG(Nmp)   ; ALOG because null propor. to mu^2 => dN ~ mudmu

IF NOT KEYWORD_SET(phie_lim) THEN BEGIN 
  phase_rms0 = 6*SQRT(moy_null-offset) ;sqrt(2*moy_null) is the max value for the phase rms given the observed mean null value moy_null; 
  phi_rms    = phase_rms0*(FINDGEN(Nrms)+1)/(Nrms) ;+/- 20% from zabcd measured value (avoid 0)
ENDIF ELSE phi_rms = MIN(phie_lim) + (MAX(phie_lim)-MIN(phie_lim))*FINDGEN(Nrms)/((Nrms-1)>1);ALOG((FINDGEN(Nrms)+1))/ALOG(Nrms)  ; ALOG because null propor. to sig^2 => dN ~ sigdsig

;conversion radian-microns
opd_mean = phi_mean*lam/(2*!pi) ;can be allowed to exceed lambda, just choose phimean > 2pi
opd_rms  = phi_rms*lam/(2*!pi)
spdthpos = spdthpos*lam/(2*!pi)

; This will explore Nas values betweeen -1% and +1% from the null value at the peak of the distribution.
IF KEYWORD_SET(neg)          THEN null_min = -0.02 ELSE null_min = 0
IF NOT KEYWORD_SET(null_lim) THEN null_ok  = [0.00>null_min, 0.03] ELSE null_ok = null_lim       ; don't replace null_lim!!!
null_expl = MIN(null_ok) + (MAX(null_ok)-MIN(null_ok))*FINDGEN(Nf)/(FLOAT(Nf-1) > 1)
IF info GT 0 THEN PRINT,'NAs step size = ',null_expl[1]-null_expl[0]

;Visibility of the fringes
vis = (1-null_expl)/(1+null_expl)

; Dsigma is bandwidth
dsigma = dl/(lam^2) ; delta_sigma of WL channel in micron^-1 ;use to be 0.28=dlambda, replaced by 0.06 for spec channels
phi0   = !pi        ;assumes no dispersion at phase = phi0 (only valid for PFN ; on KI at 2 microns, I assumed phiacro=0)

; Prior on dark_rms (hard coded for now)
IF Ndrms GT 1 THEN BEGIN   
  kmin  = 0.5
  kmax  = 2.00 
  kdark = kmin+(kmax-kmin)*float(findgen(Ndrms))/(float(Ndrms-1))
ENDIF ELSE kdark = [1]

; Generate random OPD vector
n_gen = FLOAT(N_ELEMENTS(dark)<N_ELEMENTS(a)<N_ELEMENTS(b))
IF KEYWORD_SET(seed_opd) THEN seed2 = REPLICATE(seed_opd, Nseed) ELSE seed2 = LINDGEN(Nseed)*1D5

;chi2cube = cube of chi2, dimensions: Nf x Nmp x Nrms x Nseed z Ndrms
;Where Nf is the number of sampling for the null, Nmp is the number of sampling for the mean phase
;and Nrms is the number of sampling of the rms phase
chi2cube  = FLTARR(Nf,Nmp,Nrms,Nseed,Ndrms,/NOZERO)
lklh_cube = chi2cube

; *********************
; Compute the chi2 cube 
; *********************

; Histogram scaling factor
scale = FLOAT(n_obs/n_gen)

; Constant variables accross the loops
tmp1 = 2*SQRT((a_model*b_model))

; Select between variable step and constant step (with or without background fitting)
IF KEYWORD_SET(variable) THEN BEGIN 
  FOR m = 0, Ndrms-1 DO BEGIN
    tmp0 = kdark[m]*dark_model+(a_model+b_model)
    ; Loop over OPD seed
    FOR l = 0, Nseed-1 DO BEGIN
      ; Generate RNDM  OPD vector
      rndm = RANDOMN(seed2[l], n_gen)
      IF KEYWORD_SET(bootstrap) THEN rndm = rndm[ind_resamp]
      ; Loop over phase jitter
      FOR k = 0, Nrms-1 DO BEGIN
        ; OPD randomn sequence
        opd_seq = opd_rms[k]*rndm + spdthpos
        ; Loop over mean phase
        FOR j = 0, Nmp-1 DO BEGIN
          opd  = opd_mean[j] + opd_seq 
          tmp2 = tmp1 * (SIN(!dpi*dsigma*opd)/(!dpi*dsigma*opd)) * (COS(phi0+2*!pi*opd/lam) * pha_trm1 - SIN(phi0+2*!pi*opd/lam) * pha_trm2)  ;* SIN(0.000001+!Dpi*tilt)/(0.000001+!Dpi*tilt)   ; remove call to SINC function to save computing time (make sure opd doesn't go to zero)
          ; Loop over astro nulls
          FOR i = 0, Nf-1 DO BEGIN
            n_model    = tmp0 + tmp2 * vis(i)
            hist_model = HISTO_VARIABLE(n_model,binpos,NBINS=nbinss)
            hist_model = scale*hist_model[index_min:index_max]    ; Not sure about the rescaling here...
            chi2cube[i,j,k,l,m] = TOTAL(((FLOAT(hist_model)-FLOAT(histobs))^2) / (FLOAT(hist_model)+1)) ;/ FLOAT(nbins_kept-4)  ; now do normalization later to save a little time
            lklh_cube[i,j,k,l,m] = -2*TOTAL(FLOAT(histobs)*ALOG((1.0+FLOAT(hist_model))/FLOAT(n_obs)))    ; This is a max likelihood based estimator (likelihood of data given a set of model parameters) as given by Pelat p 286
          ENDFOR
        ENDFOR
      ENDFOR
    ENDFOR
  ENDFOR
ENDIF ELSE BEGIN
   FOR m = 0, Ndrms-1 DO BEGIN
     tmp0 = kdark[m]*dark_model+(a_model+b_model)
     ; Loop over OPD seed
     FOR l = 0, Nseed-1 DO BEGIN 
      ; Generate RNDM  OPD vector
      rndm = RANDOMN(seed2[l], n_gen) 
      IF KEYWORD_SET(bootstrap) THEN rndm = rndm[ind_resamp]
      ; Loop over phase jitter
      FOR k = 0, Nrms-1 DO BEGIN
        ; OPD randomn sequence
        opd_seq = opd_rms[k]*rndm + spdthpos
        ; Loop over mean phase
        FOR j = 0, Nmp-1 DO BEGIN
          opd  = opd_mean[j] + opd_seq 
          tmp2 = tmp1 * (SIN(!dpi*dsigma*opd)/(!dpi*dsigma*opd)) * (COS(phi0+2*!pi*opd/lam) * pha_trm1 - SIN(phi0+2*!pi*opd/lam) * pha_trm2)  ;* SIN(0.000001+!Dpi*tilt)/(0.000001+!Dpi*tilt)   ; remove call to SINC function to save computing time (make sure opd doesn't go to zero)
          ; Loop over astro nulls
          FOR i = 0, Nf-1 DO BEGIN
            n_model    = tmp0 + tmp2 * vis(i)
            hist_model = scale*HISTOGRAM(n_model, nbins=nbins_kept, min=cutoff_min, max=cutoff_max) 
            chi2cube[i,j,k,l,m] = TOTAL(((FLOAT(hist_model)-FLOAT(histobs))^2) / (FLOAT(hist_model)+1)) ;/ FLOAT(nbins_kept-4)  ; now do normalization later 
            lklh_cube[i,j,k,l,m] = -2*TOTAL(FLOAT(histobs)*ALOG((1.0+FLOAT(hist_model))/FLOAT(n_obs)))  ; This is a max likelihood based estimator as given by Pelat p 286 (-2 ln(Likelihood) + K to be mopre specific)
          ENDFOR
        ENDFOR
      ENDFOR
     ENDFOR
   ENDFOR  
ENDELSE

; Normalze chi2 and likelyhood cube
lklh_cube -= MIN(lklh_cube)    ; Remove constant K from likelyhood
lklhcube = EXP(-0.5*lklh_cube) ; compute likelyhood for Bayesian stats

; Now Bayesian analysis for both chi2 and lklh approach
dof   = FLOAT(nbins_kept-4)
chi2  = MIN(chi2cube, idx1)/dof
lklh  = MAX(lklhcube, idx2)
idx   = [idx1, idx2]

; Prepare plots
fit = 20./1720.
inv     = 0
!P.FONT = 0
thick   = 4.0
xthick  = 3.5
ythick  = xthick
cthick  = 3
label1   = ['CHI2', 'LKLH']
label2   = ['NAS', 'MU', 'SIG', 'SEED_OPD', 'KDARK']
xtitle   = ['Tested null [%]', 'Tested mean phase [rad]', 'Tested phase jitter [rad]', 'Tested seed OPD', 'Test background scaling factor']
plotname = file_path

; Define output array (5 dim even if not used later)
n_dim = 5
opt   = DBLARR(2, n_dim)
bay   = opt
elow  = opt
esup  = opt

; DIMENSION of 1 for the last elemetns are not created in CUBE by for e.g. DBLARR 
FOR i = 1, 5 DO BEGIN
  IF cube_size[N_ELEMENTS(cube_size)-i] EQ 1 THEN BEGIN
    n_dim -= 1 
    ; Select data
    CASE i OF
      5: p = null_expl*1D2 ; convert null to %
      4: p = phi_mean
      3: p = phi_rms
      2: p = INDGEN(Nseed)
      1: p = kdark
    ENDCASE
    opt[*,N_ELEMENTS(cube_size)-i] = p
    bay[*,N_ELEMENTS(cube_size)-i] = p
  ENDIF ELSE i = 6
ENDFOR
dim   = INDGEN(n_dim) 

; Now compute optimum and Bayesian parameters
FOR i = 0, 1 DO BEGIN 
  ; Choose between chi2 and lklh
  IF i NE 0 THEN BEGIN
    cube    = lklh_cube
    pro_den = lklhcube
  ENDIF ELSE BEGIN
    cube = chi2cube/dof
    q = !quiet & !quiet = 1
    pro_den = CHI2_PDF(DOUBLE(chi2cube),DF=dof) ; chi2_pdf Computes proba density for chi-sq distribution
    !quiet = q
  ENDELSE
  
  ; Normalization factor
  norm = TOTAL(pro_den)   
  
  ; Index of optimum value
  idx5D = ARRAY_INDICES(pro_den, idx[i])
  
  ; Loop over dimensions of the data hypercube
  FOR j = 0, n_dim-1 DO BEGIN  
    ; Select data
    CASE j OF
      0: p = null_expl*1D2 ; convert null to %
      1: p = phi_mean
      2: p = phi_rms
      3: p = INDGEN(Nseed)
      4: p = kdark
    ENDCASE
      
    ; Only if dimension higher than 1
    IF cube_size[j] GT 1 THEN BEGIN  
           
      ; Optimum value
      opt[i,j] = DOUBLE(p[idx5D[j]])  
      
      ; Compute Bayesian value
      dim0 = dim[WHERE(dim NE j, n_dim0)]
      pdf  = pro_den
      FOR id = 1, n_dim0 DO pdf = TOTAL(pdf, dim0[n_dim0-id]+1)  ; TOTAL is 1-based!!
      pdf  /= norm
      PDFSDV, pdf, p, bay0, err_low0, err_sup0, BAR_H=bar0
      bay[i,j]  = bay0
      elow[i,j] = err_low0
      esup[i,j] = err_sup0
      
      ; Plot results if requested
      IF KEYWORD_SET(PLOT) THEN BEGIN
        ; 1. Plot PDF
        x_range = [MIN(p),MAX(p)]
        y_range = [MIN(pdf),MAX(pdf)]
        PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plotname + 'PDF_' + label1[i] + '_' + label2[j] + flag + '.eps'
        LOADCT, 0, /SILENT
        PLOTXY, p, pdf, YLOG=ylog, /NEW, XRANGE=x_range, YRANGE=y_range, TITLE=title, XTITLE=xtitle[j], YTITLE='Marginalized PDF', GRID=0, $
                XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
        LOADCT, 13, /SILENT
        PLOTXY, p, pdf, /ADD, SYMBOL=symbol, COLOR=90, THICK=thick
        PLOTXY, p, bar0+INTARR(N_ELEMENTS(p)), /ADD, SYMBOL=symbol, COLOR=0, THICK=thick, LINESTYLE=1   ; Overplot BAR
        PLOTXY, [opt[i,j],opt[i,j]], y_range, /ADD, SYMBOL=symbol, COLOR=255, THICK=thick, LINESTYLE=1   ; Overplot optimum value
        XYOUTS, x_range[0]+0.1*(x_range[1]-x_range[0]), y_range[1]-0.1*(y_range[1]-y_range[0]), 'MODE   : ' + STRING(bay0, FORMAT='(F5.2)'), CHARSIZE=1.0, CHARTHICK=2.5
        XYOUTS, x_range[0]+0.1*(x_range[1]-x_range[0]), y_range[1]-0.2*(y_range[1]-y_range[0]), 'ERR_LOW: ' + STRING(err_low0, FORMAT='(F5.2)'), CHARSIZE=1.0, CHARTHICK=2.5
        XYOUTS, x_range[0]+0.1*(x_range[1]-x_range[0]), y_range[1]-0.3*(y_range[1]-y_range[0]), 'ERR_SUP: ' + STRING(err_sup0, FORMAT='(F5.2)'), CHARSIZE=1.0, CHARTHICK=2.5
        LOADCT, 0, /SILENT
        PLOTXY, /FIN
        
        ; 2. Plot maps
        FOR k = j+1, n_dim-1 DO BEGIN
          ; Compute map for this pair
          dim1 = dim[WHERE(dim NE j AND dim NE k, n_dim1)]
          map = cube
          FOR id = 1, n_dim1 DO map = MIN(map, DIMENSION=dim1[n_dim1-id]+1)  ; MIN is 1 based
          
          ; Select data
          CASE k OF
            0: p = null_expl*1D2 ; convert null to %
            1: p = phi_mean
            2: p = phi_rms
            3: p = INDGEN(Nseed)
            4: p = kdark
          ENDCASE  
          
          ; Impose a maximum value for plotting clarity
          IF i EQ 0 THEN BEGIN
            idx10 = WHERE(map GT 10, n10)
            IF n10 GT 0 THEN map[idx10] = 10
            lk_thre = MIN(lklh_cube[idx10])        ; corresponding threshold on likelyhood
          ENDIF ELSE BEGIN
            idx10 = WHERE(map GT lk_thre, n10)
            IF n10 GT 0 THEN map[idx10] = lk_thre
          ENDELSE
  
          ; 1. First, the fine maps
          tmp     = MIN(map, idxmin)
          idx_min = ARRAY_INDICES(map, idxmin)
          y_range = [MIN(p),MAX(p)]
          PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1000, 800]*20./1720., FILENAME = plotname + 'MAP' + '_' + label1[i] + '_' + label2[j] + '-vs-' + label2[k] + flag + '.eps'
          LOADCT, 13, /SILENT
          COLORBAR, XTITLE='', YTITLE='', RANGE=[MIN(map),MAX(map)+0.0001], FORMAT='(F6.2)', CHARSIZE=charsize, CHARTHICK=1.2, /VERTICAL, /RIGHT, POSITION=[0.8,0.12,0.86,0.87]
          PLOTXY, map, YLOG=ylog, /NEW, XRANGE=x_range, YRANGE=y_range, TITLE=title, XTITLE=xtitle[j], YTITLE=xtitle[k], GRID=0, $
            XSTYLE=1, YSTYLE=1, XTHICK=4, YTHICK=4, THICK=4, CHARTHICK=1.2, /NOERASE, WINDOW=[150,100,750,700]*20./1720.
          LOADCT, 0, /SILENT
          PLOTXY, [x_range[0]+(x_range[1]-x_range[0])*idx_min[0]/(N_ELEMENTS(map[*,0])-1)], [y_range[0]+(y_range[1]-y_range[0])*idx_min[1]/(N_ELEMENTS(map[0,*])-1)], COLOR=255, SYMBOL=2, THICK=8, /ADD
          PLOTXY, /FIN
        ENDFOR
      ENDIF ELSE BEGIN
        opt[i,j] = p
        bay[i,j] = p
      ENDELSE
    ENDIF
  ENDFOR          
ENDFOR

; Define output structure
stru1     = {OPT: 0D, BAY: 0D, ERR_LOW: 0D, ERR_SUP: 0D}
stru2     = {CHI2: stru1, LKLH: stru1}
param_out = {NAS: stru2, MU: stru2, SIG: stru2, SEED_OPD: stru2, KDARK: stru2, CHI2: 0.}

; Convert NULL back to "non-%"
opt[*,0] /= 1D2
bay[*,0] /= 1D2
elow[*,0] /= 1D2
esup[*,0] /= 1D2
 
; Parse results to output structure
; 1. chi2 results
i = 0
j = 0 & param_out.nas.chi2.opt = opt[i,j]      & param_out.nas.chi2.bay = bay[i,j]      & param_out.nas.chi2.err_low = elow[i,j]      & param_out.nas.chi2.err_sup = esup[i,j]
j = 1 & param_out.mu.chi2.opt = opt[i,j]       & param_out.mu.chi2.bay = bay[i,j]       & param_out.mu.chi2.err_low = elow[i,j]       & param_out.mu.chi2.err_sup = esup[i,j]
j = 2 & param_out.sig.chi2.opt = opt[i,j]      & param_out.sig.chi2.bay = bay[i,j]      & param_out.sig.chi2.err_low = elow[i,j]      & param_out.sig.chi2.err_sup = esup[i,j]
j = 3 & param_out.seed_opd.chi2.opt = opt[i,j] & param_out.seed_opd.chi2.bay = bay[i,j] & param_out.seed_opd.chi2.err_low = elow[i,j] & param_out.seed_opd.chi2.err_sup = esup[i,j]
j = 4 & param_out.kdark.chi2.opt = opt[i,j]    & param_out.kdark.chi2.bay = bay[i,j]    & param_out.kdark.chi2.err_low = elow[i,j]    & param_out.kdark.chi2.err_sup = esup[i,j]
param_out.chi2 = chi2

; 2. lklh results
i = 1
j = 0 & param_out.nas.lklh.opt = opt[i,j]      & param_out.nas.lklh.bay = bay[i,j]      & param_out.nas.lklh.err_low = elow[i,j]      & param_out.nas.lklh.err_sup = esup[i,j]
j = 1 & param_out.mu.lklh.opt = opt[i,j]       & param_out.mu.lklh.bay = bay[i,j]       & param_out.mu.lklh.err_low = elow[i,j]       & param_out.mu.lklh.err_sup = esup[i,j]
j = 2 & param_out.sig.lklh.opt = opt[i,j]      & param_out.sig.lklh.bay = bay[i,j]      & param_out.sig.lklh.err_low = elow[i,j]      & param_out.sig.lklh.err_sup = esup[i,j]
j = 3 & param_out.seed_opd.lklh.opt = opt[i,j] & param_out.seed_opd.lklh.bay = bay[i,j] & param_out.seed_opd.lklh.err_low = elow[i,j] & param_out.seed_opd.lklh.err_sup = esup[i,j]
j = 4 & param_out.kdark.lklh.opt = opt[i,j]    & param_out.kdark.lklh.bay = bay[i,j]    & param_out.kdark.lklh.err_low = elow[i,j]    & param_out.kdark.lklh.err_sup = esup[i,j]

; Plot observed histogram and best-fit histogram
; **********************************************

; 1. OPTIMUM PARAMETERS
; Optimum parameters (use bayesian and likelihood)
nas_opt = param_out.nas.chi2.opt
mu_opt  = param_out.mu.chi2.opt
sig_opt = param_out.sig.chi2.opt
seed_opt= param_out.seed_opd.chi2.opt
drk_opt = param_out.kdark.chi2.opt ; not used anymore

; Compute model null
dsigma  = dl/(lam^2)
phi0    = !pi
opd_opt = (mu_opt + sig_opt*RANDOMN(seed_opt, n_gen))*lam/(2*!pi) + spdthpos
n_model = (drk_opt*dark_model + ABS(a+b) + 2*SQRT(ABS(a*b)) * (1 - nas_opt)/(1 + nas_opt) * SIN(!dpi*dsigma*opd_opt)/(!dpi*dsigma*opd_opt) * (COS(phi0+2*!pi*opd_opt/lam) * pha_trm1 - SIN(phi0+2*!pi*opd_opt/lam) * pha_trm2)) / phot_tot
n_obs   = arraydata.ab / phot_tot 

; Compute data range
data_min = MIN(n_obs) < MIN(n_model)
data_max = MAX(n_obs) > MAX(n_model)
data_min -= 0.5*0.05*(data_max-data_min)
data_max += 0.5*0.05*(data_max-data_min)
xrange   = 100.*[data_min,data_max]

; Compute histogram
n_bin      = ROUND(SQRT(n_gen))
model_hist = HISTOGRAM(n_model, NBINS=n_bin, MAX=data_max, MIN=data_min, LOCATIONS=bins)
obs_hist   = HISTOGRAM(n_obs, NBINS=n_bin, MAX=data_max, MIN=data_min, LOCATIONS=bins)
scale      = TOTAL(model_hist)/TOTAL(obs_hist)
model_hist = model_hist/scale
bins       = bins + 0.5*(bins[1]-bins[0]) ; Shift to bin center
yrange     = [0.,MAX(obs_hist)*1.2]

PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plotname + 'histrograms_OPT' + flag + '.eps'
PLOTXY, 100.*bins, model_hist, PSYM=10, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE=title, XTITLE='Instantaneous null depth [%]', YTITLE= 'Number of occurence', GRID=0, $
  XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
LOADCT, 13, /SILENT
PLOTXY, 100.*bins, model_hist, /ADD, PSYM=10, COLOR=90, THICK=7.0, LINESTYLE=1
PLOTXY, 100.*bins, obs_hist, /ADD, PSYM=10, COLOR=0, THICK=thick
LOADCT, 0, /SILENT
PLOTXY, /FIN

; 2. CHI2 paramaters
; Optimum parameters (use bayesian and likelihood)
nas_opt = param_out.nas.chi2.bay
mu_opt  = param_out.mu.chi2.bay
sig_opt = param_out.sig.chi2.bay
seed_opt= param_out.seed_opd.chi2.bay
drk_opt = param_out.kdark.chi2.bay ; not used anymore

; Compute model null
opd_opt = (mu_opt + sig_opt*RANDOMN(seed_opt, n_gen))*lam/(2*!pi) + spdthpos
n_model = (drk_opt*dark_model + ABS(a+b) + 2*SQRT(ABS(a*b)) * (1 - nas_opt)/(1 + nas_opt) * SIN(!dpi*dsigma*opd_opt)/(!dpi*dsigma*opd_opt) * (COS(phi0+2*!pi*opd_opt/lam) * pha_trm1 - SIN(phi0+2*!pi*opd_opt/lam) * pha_trm2)) / phot_tot

; Compute histogram
n_bin      = ROUND(SQRT(n_gen))
model_hist = HISTOGRAM(n_model, NBINS=n_bin, MAX=data_max, MIN=data_min, LOCATIONS=bins)
model_hist = model_hist/scale
bins       = bins + 0.5*(bins[1]-bins[0]) ; Shift to bin center

PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plotname + 'histrograms_CHI2' + flag + '.eps'
PLOTXY, 100.*bins, model_hist, PSYM=10, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE=title, XTITLE='Instantaneous null depth [%]', YTITLE= 'Number of occurence', GRID=0, $
  XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
LOADCT, 13, /SILENT
PLOTXY, 100.*bins, model_hist, /ADD, PSYM=10, COLOR=90, THICK=7.0, LINESTYLE=1
PLOTXY, 100.*bins, obs_hist, /ADD, PSYM=10, COLOR=0, THICK=thick
LOADCT, 0, /SILENT
PLOTXY, /FIN

; 3. Likelihood parameters
; Optimum parameters (use bayesian and likelihood)
nas_opt = param_out.nas.lklh.bay
mu_opt  = param_out.mu.lklh.bay
sig_opt = param_out.sig.lklh.bay
seed_opt= param_out.seed_opd.lklh.bay
drk_opt = param_out.kdark.lklh.bay ; not used anymore

; Compute model null
opd_opt = (mu_opt + sig_opt*RANDOMN(seed_opt, n_gen))*lam/(2*!pi) + spdthpos
n_model = (drk_opt*dark_model + ABS(a+b) + 2*SQRT(ABS(a*b)) * (1 - nas_opt)/(1 + nas_opt) * SIN(!dpi*dsigma*opd_opt)/(!dpi*dsigma*opd_opt) * (COS(phi0+2*!pi*opd_opt/lam) * pha_trm1 - SIN(phi0+2*!pi*opd_opt/lam) * pha_trm2)) / phot_tot

; Compute histogram
n_bin      = ROUND(SQRT(n_gen))
model_hist = HISTOGRAM(n_model, NBINS=n_bin, MAX=data_max, MIN=data_min, LOCATIONS=bins)
scale      = TOTAL(model_hist)/TOTAL(obs_hist)
model_hist = model_hist/scale
bins       = bins + 0.5*(bins[1]-bins[0]) ; Shift to bin center

PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plotname + 'histrograms_LKLH' + flag + '.eps'
PLOTXY, 100.*bins, model_hist, PSYM=10, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE=title, XTITLE='Instantaneous null depth [%]', YTITLE= 'Number of occurence', GRID=0, $
  XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
LOADCT, 13, /SILENT
PLOTXY, 100.*bins, model_hist, /ADD, PSYM=10, COLOR=90, THICK=7.0, LINESTYLE=1
PLOTXY, 100.*bins, obs_hist, /ADD, PSYM=10, COLOR=0, THICK=thick
LOADCT, 0, /SILENT
PLOTXY, /FIN

; 2. Model with negative best-fit Nas
diagnostic = 1
IF diagnostic EQ 1 THEN BEGIN
  cube = REFORM(lklhcube[0,*,*,*,*])  
  n_dim   = (SIZE(cube))[0]
  nas_opt = MIN(null_expl)
  tmp     = MAX(cube, idx3D)
  idx     = ARRAY_INDICES(cube, idx3D)
  mu_opt  = phi_mean[idx[0]]
  sig_opt = phi_rms[idx[1]]
  IF N_ELEMENTS(idx) GT 2 THEN seed_opt= idx[2] ELSE seed_opt = 0
  IF N_ELEMENTS(idx) GT 3 THEN drk_opt = kdark[idx[3]] ELSE drk_opt = 1

  ; Compute model null
  opd_opt = (mu_opt + sig_opt*RANDOMN(seed_opt, n_gen))*lam/(2*!pi) + spdthpos
  n_model = (drk_opt*dark_model + ABS(a+b) + 2*SQRT(ABS(a*b)) * (1 - nas_opt)/(1 + nas_opt) * SIN(!dpi*dsigma*opd_opt)/(!dpi*dsigma*opd_opt) * (COS(phi0+2*!pi*opd_opt/lam) * pha_trm1 - SIN(phi0+2*!pi*opd_opt/lam) * pha_trm2)) / phot_tot

  n_bin      = ROUND(SQRT(n_gen))
  model_hist = HISTOGRAM(n_model, NBINS=n_bin, MAX=data_max, MIN=data_min, LOCATIONS=bins)
  scale      = TOTAL(model_hist)/TOTAL(obs_hist)
  model_hist = model_hist/scale
  bins       = bins + 0.5*(bins[1]-bins[0]) ; Shift to bin center

  PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plotname + 'histrograms_MIN(NAS)' + flag + '.eps'
  PLOTXY, 100.*bins, model_hist, PSYM=10, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE=title, XTITLE='Instantaneous null depth [%]', YTITLE= 'Number of occurence', GRID=0, $
    XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
  LOADCT, 13, /SILENT
  PLOTXY, 100.*bins, model_hist, /ADD, PSYM=10, COLOR=90, THICK=7.0, LINESTYLE=1
  PLOTXY, 100.*bins, obs_hist, /ADD, PSYM=10, COLOR=0, THICK=thick
  LOADCT, 0, /SILENT
  PLOTXY, /FIN
ENDIF

; Plot observed and best-fit null sequence (based on likelyhood)
; *************************************************************
   
x_obs   = INDGEN(N_ELEMENTS(n_obs))
x_model = INDGEN(N_ELEMENTS(n_model))
PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = plotname + 'null_sequence_LKLH' + flag + '.eps'
PLOTXY, x_obs, 1D2*n_obs, PSYM=10, /NEW, XRANGE=[MIN(x_obs),MAX(x_obs)], YRANGE=[-1,2*1D2*MAX(n_obs)], TITLE=title, YTITLE='Instantaneous null depth [%]', XTITLE= 'File #', GRID=0, $
        XSTYLE=1, YSTYLE=1, XTHICK=xthick, YTHICK=ythick, THICK=thick, CHARTHICK=cthick, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
LOADCT, 13, /SILENT
PLOTXY, x_obs, 1D2*n_obs, /ADD, PSYM=10, COLOR=0, THICK=thick;, LINESTYLE=1
PLOTXY, x_model, 1D2*n_model, /ADD, PSYM=10, COLOR=90, THICK=1
LOADCT, 0, /SILENT
PLOTXY, /FIN     
END