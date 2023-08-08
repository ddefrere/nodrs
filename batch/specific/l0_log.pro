
; REDUCTION NOTES
; ***************

; Feb. 12nd 2014 (Weird transcient behavior at the beginning of each sequence) => remove the first 15 (=n_bad) frames of each new sequence
; FLAT_IDX=[64541,64550], DARK_IDX=[64131,64530]
'140212':
n_trans = 20

; March 17 2014 (data acquistion started too soon for each nod, AO loop was not yet closed)
'140317':
n_trans = 20

'141109':
LBTI_FITSUTIL, '141109', [35434,35454], LBT_LXOS=0.5, LBT_LYOS=-1, LBT_RXOS=0.5, LBT_RYOS=-1
bad_idx = 35434 + INDGEN(11)
remove_ob = [20,21] ; 20,21: bad background

'131231':
;data_idx=[10667,16666]
;nod_idx=[[10667,14666],[14667,16666]]
;

'150203'
LBTI_FITSUTIL, '150203', OBJNAME='HD112769'

'150204'
LBTI_FITSUTIL, '150204', [0,800], OBJNAME='mu_gem'
LBTI_FITSUTIL, '150204', [17000,27246], OBJNAME='HD48433
LBTI_FITSUTIL, '150204', [9418,11350], OBSTYPE=3, /LMIRCAM

'150208'
LBTI_FITSUTIL, '150208', [40807,48806], FLAG='CAL'
LBTI_FITSUTIL, '150208', [60807,63346], OBJNAME='DOME', DATATYPE=1 ; 'Dark frames'
