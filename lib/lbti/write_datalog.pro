;+
; NAME: WRITE_DATALOG
; 
; PURPOSE:
;   Create a datalog file that goes with the intermediate L0 files and used for the L1 processing.
;
; MANDATORY INPUT:
;   date  : date of the processed data (long format: 20YY-MM-DD)
;
; KEYWORDS
;   INSTRUM :
;   LOGFILE :
;
; MODIFICATION HISTORY:
;   Version 1.0, 11-FEB-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 06-JUL-2016, DD: date is now long format
;   Version 1.2, 06-OCT-2016, DD: added second output to datalog (one line per pointing instead of nod)

PRO WRITE_DATALOG, date_obs, INSTRUM=instrum, LOGFILE=logfile
   
  ; Find path
  DECLARE_PATH, pth, INSTRUM=instrum
  
  ; Define log file and restore it
  IF NOT KEYWORD_SET(LOGFILE) THEN logfile = pth.l0fits_path + date_obs + pth.sep + 'datalog.sav'
  RESTORE, logfile
  
  ; 1. Create new log per nod
  OPENW,  log,  pth.l0fits_path + date_obs + pth.sep + 'datalog.txt', /GET_LUN, WIDTH=500
  PRINTF, log, 'Data log created on ' + SYSTIME()
  PRINTF, log, ''
  PRINTF, log, 'DATALOG PER NOD'
  PRINTF, log, 'NOD_ID;CFG_ID;UT_TIME;OBJNAME;FLAG;DATATYPE;OBSTYPE;INT_TIME;LAM_CEN;BANDWDTH;XCEN_SX;YCEN_SX;XCEN_DX;YCEN_DX;MIN_FID;MAX_FID;PLSTP;TLTSP;PID;DATE_RED'
  n_line = N_ELEMENTS(data_r.nod_id)
  ; Print first part of the log
  FOR i=0, n_line-1 DO PRINTF, log, STRING(data_r[i].nod_id, FORMAT='(I03)'), ';', STRING(data_r[i].cfg_id, FORMAT='(I0)'), ';', STRCOMPRESS(STRING(data_r[i].ut_time), /REMOVE_ALL), ';', $
                       STRCOMPRESS(STRING(data_r[i].objname), /REMOVE_ALL), ';', STRCOMPRESS(STRING(data_r[i].flag), /REMOVE_ALL), ';', STRING(data_r[i].datatype, FORMAT='(I0)'), ';', $
                       STRING(data_r[i].obstype, FORMAT='(I0)'), ';', STRING(data_r[i].int_time, FORMAT='(E8.2)'), ';', STRING(data_r[i].lam_cen, FORMAT='(E8.2)'), ';', $
                       STRING(data_r[i].bandwdth, FORMAT='(E8.2)'), ';', STRING(ROUND(data_r[i].xcen_sx), FORMAT='(I0)'), ';', STRING(ROUND(data_r[i].ycen_sx), FORMAT='(I0)'), ';', $
                       STRING(ROUND(data_r[i].xcen_dx), FORMAT='(I0)'), ';', STRING(ROUND(data_r[i].ycen_dx), FORMAT='(I0)'), ';', STRING(data_r[i].min_fid, FORMAT='(I0)'), ';',$
                       STRING(data_r[i].max_fid, FORMAT='(I0)'), ';', STRING(data_r[i].pcplsp, FORMAT='(I0)'), ';', STRING(data_r[i].pctltsp, FORMAT='(I0)'), ';', $
                       STRING(data_r[i].pid, FORMAT='(I0)'), ';', data_r[i].date_red                       
  ; Now second part of the log
  ; Derive number of pointing
  pt_id   = data_r.pt_id
  pt_uniq = pt_id[UNIQ(pt_id,  SORT(pt_id))]
  n_pt    = N_ELEMENTS(pt_uniq)
  PRINTF, log, ''
  PRINTF, log, 'DATALOG PER POINTING'
  PRINTF, log, 'PT_ID;UT_TIME;OBJNAME;FLAG;DATATYPE;N_NOD;N_CFG;MIN_FID;MAX_FID;PID
  FOR i_pt = 0, n_pt-1 DO BEGIN
    idx_pt   = WHERE(pt_id EQ pt_uniq[i_pt])
    ut_time  = data_r[idx_pt[0]].ut_time
    objname  = data_r[idx_pt[0]].objname
    flag     = data_r[idx_pt[0]].flag
    datatype = data_r[idx_pt[0]].datatype
    nod_id   = data_r[idx_pt].nod_id
    nod_uniq = nod_id[UNIQ(nod_id,  SORT(nod_id))] & n_nod = N_ELEMENTS(nod_uniq)
    cfg_id   = data_r[idx_pt].cfg_id
    cfg_uniq = cfg_id[UNIQ(cfg_id,  SORT(cfg_id))] & n_cfg = N_ELEMENTS(cfg_uniq)
    min_fid  = MIN(data_r[idx_pt].min_fid)
    max_fid  = MAX(data_r[idx_pt].max_fid)
    pid      = data_r[idx_pt[0]].pid
    PRINTF, log, STRING(pt_uniq[i_pt], FORMAT='(I03)'),  ';', STRCOMPRESS(STRING(ut_time), /REMOVE_ALL), ';', STRCOMPRESS(STRING(objname), /REMOVE_ALL), ';', STRCOMPRESS(STRING(flag), /REMOVE_ALL), ';', $
                 STRING(datatype, FORMAT='(I0)'), ';', STRING(n_nod, FORMAT='(I0)'), ';', STRING(n_cfg, FORMAT='(I0)'), ';', STRING(min_fid, FORMAT='(I0)'), ';', STRING(max_fid, FORMAT='(I0)'), ';', $
                 STRING(pid, FORMAT='(I0)')
  ENDFOR  
  CLOSE, log
  FREE_LUN, log
    
END