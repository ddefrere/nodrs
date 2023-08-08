FUNCTION FIND_FREQIDX, data, THRESHOLD=threshold, MIN_SEQ = min_seq

; PRINCIPLE
; 1. Find the array index where data changes
; 2. Make the difference of successive positions where it changes 
; 3. Look at the difference of the difference. If it changes by more than the threshold over at least "MIN_SEQ", create a new sequence!

IF NOT KEYWORD_SET(threshold) THEN threshold = ABS(MAX(data)-MIN(data))/2.

; Find where an element is not equal to its next neighbor
idx = WHERE(data NE SHIFT(data, 1), n_diff)

; If the first and last element are not the same, remove first element
IF data[0] NE data[-1] THEN idx = idx[1:*]

; If no difference, return
IF n_diff LE 0 THEN RETURN, 0

; Compute the difference between the index of changing frames
diff = idx-SHIFT(idx,1)
diff = diff[1:*]  ; Remove the first element which is always an artefact
 
; Compute where this difference is above the threshold
idx2 = WHERE(ABS(diff-SHIFT(diff, 1)) GE threshold, n_diff)

; If no difference, return
IF n_diff LE 0 THEN RETURN, 0

; If the first and last element are not the same, remove first element
IF diff[0] NE diff[-1] AND n_diff GT 1 THEN idx2 = idx2[1:*]

; Keep only the sequences with at least min_seq consecutive repetition
idx_good = WHERE(ABS(idx2-SHIFT(idx2, 1)) GE min_seq, n_good)
IF n_good LE 0 THEN RETURN, 0

; Remove duplicate (just to be sure)
idx_seq = [0, idx[idx2[idx_good]]]
idx_seq = idx_seq[UNIQ(idx_seq, SORT(idx_seq))]

RETURN, idx_seq
END