;+
; NAME: SET_VARSIZE_BINS
; 
; DESCRIPTION
;   Determines the bin positions in the variable-bin-size histogram
;
; INPUT
;   vector: the input data
;
; KEWORD
;   N_BIN : the number of bins desired in the histogram (if not set, the SQRT of the number of elements in vector will be used)
;
; OUTPUT
;   Positions of the bins and number of bins after resampling
;
; MODIFICATION HISTORY
;   Version 1.0, 27-OCT-2012, by Olivier Absil, ULg, absil@astro.ulg.ac.be
;   Version 1.1, 30-OCT-2012, OAb: corrected the definitions of occmean and lims
;   Version 1.2, 27-NOV-2015, OAb: cleaned up comments

FUNCTION SET_VARSIZE_BINS, vector, N_BIN=n_bin

n_vec = N_ELEMENTS(vector)
IF NOT KEYWORD_SET(n_bin) THEN n_bin = ROUND(SQRT(n_vec)) ; number of bins

;The CEIL function returns the closest integer greater than or equal to its argument. 
occmean = FLOAT(n_vec)/FLOAT(n_bin) ; mean number of occurence per bin for all bins except the last one

vecsort = vector[SORT(vector)] ; sort in increasing order

lims = ROUND(occmean*INDGEN(n_bin)) ; define the indices of the bin 
;boundaries in vecsort --> [0, occmean, 2*occmean, ..., n_bin*occmean]

binpos = (vecsort[(lims-1)>0]+vecsort[lims])/2 ; compute the bin 
;positions as the mean of the 2 surrounding vecsort values -- 
;(limits-1)>0 make sure that the first element is properly computed
; Note: the last element of binpos is equal to the starting point of the 
;last bin (which extends to +infinity)
;nbins_resampled = N_ELEMENTS(binpos)
;PRINT, nbins_resampled

RETURN, binpos

END