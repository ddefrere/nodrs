; +
; NAME: LBTI_DATALOG
; 
; PURPOSE:
;   This procedure saves calibrated L0 images and relevant keywords in FITS files.
;
; MANDATRY INPUTS:
;   log_file  :  string with the name of the IDL log file (must be '.sav')
;   hdr_in    :  structure with the corresponding header information
;
; KEYWORDS:
;
; MODIFICATION HISTORY:
;   Version 1.0, 18-SEP-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 04-APR-2015, DD: Added nodding period
;   Version 1.2, 14-OCT-2015, DD: Modified configuration definition
;   Version 1.3, 29-OCT-2015, DD: Removed definition of config ID (now assign right from the start and pass through header information)
;   Version 1.4, 02-DEC-2015, DD: Added pointing ID to output structure
;   Version 1.5, 21-JAN-2016, DD: Corrected implementation of config ID

PRO LBTI_DATALOG, log_file, hdr_in

; Extract main parameters 
data = {CFG_ID: hdr_in[0].cfg_id[0], NOD_ID: hdr_in[0].nod_id[0], PT_ID:hdr_in[0].pt_id[0], UT_TIME: hdr_in[0].lbt_utc[0], OBJNAME: hdr_in[0].objname[0], FLAG: hdr_in[0].flag[0], DATATYPE: hdr_in[0].datatype[0], $
        OBSTYPE: hdr_in[0].obstype[0], INT_TIME: hdr_in[0].int_time[0], N_COADD: hdr_in[0].n_coadd[0], SMPLMODE: hdr_in[0].smplmode[0], PAGAIN: hdr_in[0].pagain[0], $
        PABANDW: hdr_in[0].pabandw[0], DETMODE: hdr_in[0].detmode[0], PACTCDLY: hdr_in[0].pactcdly[0], LAM_CEN: hdr_in[0].lam_cen[0], BANDWDTH: hdr_in[0].bandwidt[0],$
        XCEN_SX: MEAN(hdr_in.xcen[0]), YCEN_SX: MEAN(hdr_in.ycen[0]), XCEN_DX: MEAN(hdr_in.xcen[1]), YCEN_DX: MEAN(hdr_in.ycen[1]), MIN_FID: MIN(hdr_in.file_id),$
        MAX_FID: MAX(hdr_in.file_id), PCPLSP: hdr_in[0].pcplsp[0], PCTIPSP: hdr_in[0].pctipsp[0], PCTLTSP: hdr_in[0].pctltsp[0], PID: hdr_in[0].pid[0], DATE_RED: systime()}
       
; Update processed L0 log file 
IF FILE_TEST(log_file) THEN BEGIN
  ; Restore file
  RESTORE, log_file
  ; Check whether this nod and config ID have already been saved
  idx_nod = WHERE(data.nod_id EQ data_r.nod_id AND data.cfg_id EQ data_r.cfg_id, n_match)
  ; If already exists, replace the entry. If not, add a new line for this nod/cfg_id pair
  IF n_match GT 0 THEN data_r[idx_nod] = data ELSE data_r = [data_r, data]
ENDIF ELSE data_r = data

; Sort data according to MIN_FID (chronological order)
data_r = data_r[SORT(data_r.MIN_FID)]
SAVE, data_r, FILENAME=log_file         
END