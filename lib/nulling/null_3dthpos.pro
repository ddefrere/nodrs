;+
; NAME: NULL_3DTHPOS
; 
; PURPOSE:
;   Compute null from 3-position dither pattern
;
; INPUTS:
;   null          :  Input null measurements
;   null_err      :  Corresponding errors
;
; KEYWORDS
;   RATIO         :  Ratio of the input data to bin
;   MODE          :  If set, return the mode of the input data
;   INFO
;   PLOT
;
; OUTPUT
;   
;
; MODIFICATION HISTORY:
;   Version 1.0,  01-AUG-2018, by Denis Defrère, University of Liege, denis@lbti.org

FUNCTION NULL_3DTHPOS, i1, i2, null, dthpos, pcphmean, PLOT=plot

; Size of the sequence
n_data = N_ELEMENTS(null)

; Find where the ditther pattern changes
idx_seq = WHERE(dthpos NE SHIFT(dthpos, 1), n_seq)
IF n_seq LE 2 THEN RETURN, 0
IF dthpos[0] EQ dthpos[-1] THEN idx_seq =[0, idx_seq]

; Number of elements
n_seq   = N_ELEMENTS(idx_seq)
idx_seq = [idx_seq, n_data]

;Check whether it can be divided by three

; Loop through the sequence
is = 0
bias_out = 0
phi_out  = 0
null_out = 0
WHILE is LE (n_seq-3) DO BEGIN
  ; Extract dither angles
  dthpos1 = dthpos[idx_seq[is]]
  dthpos2 = dthpos[idx_seq[is+1]]
  dthpos3 = dthpos[idx_seq[is+2]]
  ; Extract dither angles
  pcphmean1 = pcphmean[idx_seq[is]:(idx_seq[is+1]-1)]
  pcphmean2 = pcphmean[idx_seq[is+1]:(idx_seq[is+2]-1)]
  pcphmean3 = pcphmean[idx_seq[is+2]:(idx_seq[is+3]-1)]
  ; If all different, continue
  IF dthpos1 NE dthpos2 AND dthpos1 NE dthpos3 AND dthpos2 NE dthpos3 THEN BEGIN
    ; Check whether at least one is 0 and two other are opposite
    IF (dthpos1 EQ 0 OR dthpos2 EQ 0 OR dthpos3 EQ 0) AND (ABS(dthpos1) EQ ABS(dthpos2) OR ABS(dthpos1) EQ ABS(dthpos3) OR ABS(dthpos2) EQ ABS(dthpos3)) THEN BEGIN
      ; Extract null values
      IF dthpos1 EQ 0 THEN BEGIN
        dither    = ABS(dthpos2)
        null_pos1 = MEAN(null[idx_seq[is]:(idx_seq[is+1]-1)])   ; this is the 0 position
        null_pos2 = MEAN(null[idx_seq[is+1]:(idx_seq[is+2]-1)]) 
        null_pos3 = MEAN(null[idx_seq[is+2]:(idx_seq[is+3]-1)])
      ENDIF
      IF dthpos2 EQ 0 THEN BEGIN
        dither    = ABS(dthpos3)
        null_pos2 = MEAN(null[idx_seq[is]:(idx_seq[is+1]-1)])   
        null_pos1 = MEAN(null[idx_seq[is+1]:(idx_seq[is+2]-1)]) ; this is the 0 position
        null_pos3 = MEAN(null[idx_seq[is+2]:(idx_seq[is+3]-1)])
      ENDIF
      IF dthpos3 EQ 0 THEN BEGIN
        dither    = ABS(dthpos1)
        null_pos2 = MEAN(null[idx_seq[is]:(idx_seq[is+1]-1)])
        null_pos3 = MEAN(null[idx_seq[is+1]:(idx_seq[is+2]-1)]) 
        null_pos1 = MEAN(null[idx_seq[is+2]:(idx_seq[is+3]-1)]) ; this is the 0 position
      ENDIF    ; Extract null sequence
      ; Compute visibility, bias, and phase
      bias = 0.5*(null_pos2+null_pos3-2.*COS(dither)*null_pos1)/(1.-COS(dither))-(i1+i2)
      phi  = ATAN((null_pos2-null_pos3),2.*(null_pos1-i1-i2-bias)*SIN(dither))
      phi  = (null_pos2-null_pos3)/(2.*(null_pos1-i1-i2-bias)*SIN(dither))
      vis  = 0.5*(null_pos1-i1-i2-bias)/(COS(dither)*SQRT(i1*i2))
      nas  = (1.-ABS(vis))/(1.+ABS(vis))
      ; Store value
      bias_out = [bias_out, bias]
      phi_out  = [phi_out, bias]
      ;null_out = [null_out, nas]
      ; Diagnostic plot to make sure that the best model goes throught the 3 measurements
      ;IF KEYWORD_SET(plot) THEN BEGIN
      ;  n_pt  = 1000
      ;  alpha = -5*dither+10*dither*DINDGEN(n_pt)/(n_pt-1)
      ;  model = i1+i2+2*vis*SQRT(i1*i2)*COS(alpha+phi)+bias
      ;  PLOT, [dthpos1, dthpos2, dthpos3], [MEAN(null[idx_seq[is]:(idx_seq[is+1]-1)]), MEAN(null[idx_seq[is+1]:(idx_seq[is+2]-1)]), MEAN(null[idx_seq[is+2]:(idx_seq[is+3]-1)])]
      ;  OPLOT, alpha, model
      ;ENDIF
      ; POLY_FIT approach
      coeff = POLY_FIT([dthpos1+pcphmean1, dthpos2+pcphmean2, dthpos3+pcphmean3], [null[idx_seq[is]:(idx_seq[is+1]-1)], null[idx_seq[is+1]:(idx_seq[is+2]-1)], null[idx_seq[is+2]:(idx_seq[is+3]-1)]], 2, /DOUBLE)
      IF coeff[2] GT 0 THEN BEGIN
        n_pt  = 1000
        alpha = -5*dither+10*dither*DINDGEN(n_pt)/(n_pt-1)
        model = coeff[2]*alpha^2+coeff[1]*alpha+coeff[0]
        nas   = MIN(model)/(i1+i2+2*SQRT(I1*i2))
        IF KEYWORD_SET(plot) THEN BEGIN          
          ;PLOT, [dthpos1+pcphmean1, dthpos2+pcphmean2, dthpos3+pcphmean3], [null[idx_seq[is]:(idx_seq[is+1]-1)], null[idx_seq[is+1]:(idx_seq[is+2]-1)], null[idx_seq[is+2]:(idx_seq[is+3]-1)]]
          PLOT, [dthpos1+MEAN(pcphmean1), dthpos2+MEAN(pcphmean2), dthpos3+MEAN(pcphmean3)], [null_pos1, null_pos2, null_pos3], XTITLE="Dither angle", YTITLE="Neasured null", PSYM=4
          PLOT, [dthpos1, dthpos2, dthpos3], [null_pos1, null_pos2, null_pos3], XTITLE="Dither angle", YTITLE="Neasured null", PSYM=4
          OPLOT, alpha, model
          PRINT, 100*nas
        ENDIF        
        null_out = [null_out, nas]
        ; Increment to the next three sequences
        is += 3
      ENDIF ELSE is += 1
    ENDIF ELSE is += 1 ;only increment one in this case
  ENDIF ELSE is += 1 ;only increment one in this case
ENDWHILE

; Remove firt entry
bias_out = bias_out[1:*]
phi_out  = phi_out[1:*]
null_out = null_out[1:*]

; Compute final value
AVGSDV, null_out, null_avg, null_rms, null_err, KAPPA=5

; Output structure
null = {NAS: null_avg, NAS_ERR: null_err}

RETURN, null
END