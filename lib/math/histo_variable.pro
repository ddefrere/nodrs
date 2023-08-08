;+
; NAME: HISTO_VARIABLE
; 
; DESCRIPTION
;   Compute the histogram with more or less the same number of occurence per bin
;
; INPUT
;   vector : input data
;   binpos : vector containing the null values at the boundary of the bins in the variable-bin-size histogram
;
; KEYWORD
;   nbins : set this keyword if you want the histogram to have at least a number "nbins" of bins (zeros are padded if histogram too short)
;
; OUTPUT
;   The variable-bin-size histogram
;
; MODIFICATION HISTORY
;   Version 1.0, 27-OCT-2012, by Olivier Absil, ULg, absil@astro.ulg.ac.be
;   Version 1.1, 06-NOV-2012, LM/OA: introduced keyword NBINS and implemented zero-padding at the end of the histogram
;   Version 1.2, 06-NOV-2015, DD: cleaned up comments

FUNCTION HISTO_VARIABLE, vector, binpos, NBINS=nbins

;The VALUE_LOCATE function finds the intervals within a given monotonic vector 
;that brackets a given set of one or more search values. This function is useful 
;for interpolation and table-lookup, and is an adaptation of the locate() routine 
;in Numerical Recipes. VALUE_LOCATE uses the bisection method to locate the interval. 
;Syntax 
;Result = VALUE_LOCATE ( Vector, Value [, /L64 ] ) 
;Each return value, Result [i], is an index, j, into Vector, corresponding to the interval into which the given Value [i] falls. 

; FIRST VERSION
;vecsort = vector[SORT(vector)] 
;index = VALUE_LOCATE(vecsort,binpos)
;hist = index - SHIFT(index,1);step :nb of elts/bin = cst
;hist = [hist[1:*],N_ELEMENTS(vecsort)-MAX(index)]

; SECOND VERSION (somewhat faster)
mappedData = VALUE_LOCATE(binpos, vector)
hist = HISTOGRAM(mappedData, MIN=0, OMAX=maxbin)
IF KEYWORD_SET(nbins) THEN BEGIN
  IF maxbin LT nbins-1 THEN hist = [TEMPORARY(hist), DBLARR(nbins-maxbin)]
ENDIF

RETURN, hist

END