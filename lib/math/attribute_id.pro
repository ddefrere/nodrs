;+
; NAME: ATTRIBUTE_ID
;   
; PURPOSE:
;   Attribute an ID to each line of the input array such as each identical line receives the same ID.
;   The result is stored in a vector that contains as many elements as the number of lines in 'data'.
;
; INPUTS:
;   data    : 1-D or 2-D input array
;
; KEYWORD:
;   PERMUTE : Set this keyword to the index of columns within which the element order does not matter (0-based);
;
; MODIFICATION HISTORY:
;   Version 1.0, 03-APR-2012, by Denis Defrere, ddefrere@mpifr.de
;   Version 1.1, 03-APR-2013, DD, added keyword PERMUTE and simplified implementation

FUNCTION ATTRIBUTE_ID, data, PERMUTE= permute

; First, find the number of columns and lines in the input data array
n_col    = N_ELEMENTS(data[*,0]) 
n_lin    = N_ELEMENTS(data[0,*]) 
data_idx = DBLARR(n_lin) & data_idx2 = data_idx
    
; If PERMUTE is set, sort first the corresponding columns 
IF KEYWORD_SET(PERMUTE) THEN BEGIN
  n_per = N_ELEMENTS(permute)
  IF n_per GT 1 THEN BEGIN
    sub_arr = data[permute,*]
    FOR i_lin = 0, n_lin-1 DO BEGIN
      idx_srt             = SORT(sub_arr[*,i_lin]) 
      data[permute,i_lin] = sub_arr[idx_srt,i_lin]
    ENDFOR
  ENDIF
ENDIF

; Loop over each other column
n_idx = 1                                  ; initiate the number of different index
FOR i_col = 0, n_col-1 DO BEGIN
  ; Extract the current 'i_col' column
  column  = data[i_col,*]                  
  idx_cnt = 0                              ; initiate index counter 
  ; Loop over the different IDs of previous column
  FOR i_idx = 0, n_idx-1 DO BEGIN
   ; Extract data of each indix
   idx     = WHERE(data_idx EQ i_idx, ni)                         
   sub_col = column[idx]   
   ; Find unique element in the sub column
   col_dis = sub_col(UNIQ(sub_col, SORT(sub_col)))      ; distinct input elements of current sub-column    
   n_idx2  = N_ELEMENTS(col_dis)                        ; number of distinct input elements of current column, i.e. the number of different IDs
   ; Now loop over the distinct values in this sub-sub-column
   FOR j_idx = 0, n_idx2-1 DO BEGIN
     idx2            = idx[WHERE(sub_col EQ col_dis[j_idx], nj)]   ; Extract data of the same value
     data_idx2[idx2] = j_idx + idx_cnt                             ; Attribute the indix
   ENDFOR
   idx_cnt = idx_cnt + n_idx2
  ENDFOR 
  data_idx = data_idx2                                            ; Store the new data_idx table 
  n_idx    = MAX(data_idx) + 1                                    ; Re-initiate the number of different index
ENDFOR

RETURN, data_idx
END