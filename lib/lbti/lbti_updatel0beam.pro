PRO LBTI_UPDATEL0BEAM, date, nod, pos_sx, pos_dx
  
  IF NOT KEYWORD_SET(pos_dx) THEN pos_dx = [0,0]
  
  ; Declare path
  DECLARE_PATH, pth, INSTRUM='NOMIC'
  
  ; --- Long-format date
  date_obs = '20' + STRMID(date, 0, 2) + '-' + STRMID(date, 2, 2) + '-' + STRMID(date, 4, 2)
  
  ; --- Restore LO log file and read useful data
  datalog = pth.l0fits_path + date_obs + pth.sep + 'datalog.sav'
  IF NOT FILE_TEST(datalog) THEN BEGIN
    MESSAGE, 'No data log file found!', /CONTINUE
    RETURN
  ENDIF ELSE RESTORE, datalog

  ; -- Update position
  idx_nod = WHERE(data_r.nod_id EQ nod, n_nod)
  IF n_nod LE 0 THEN MESSAGE, 'No corresponding nod!'
  
  ; --- If ok, make a copy
  datalog_bu = pth.l0fits_path + date_obs + pth.sep + 'datalog_bu.sav'
  IF NOT FILE_TEST(datalog_bu) THEN FILE_COPY, datalog, datalog_bu 
  utc = data_r[idx_nod].ut_time
  data_r[idx_nod].xcen_sx = pos_sx[0]
  data_r[idx_nod].ycen_sx = pos_sx[1]
  data_r[idx_nod].xcen_dx = pos_dx[0]
  data_r[idx_nod].ycen_dx = pos_dx[1]
  
  ; --- Save file
  SAVE, data_r, FILENAME=datalog
  
  ; --- Write DATALOG
  WRITE_DATALOG, date_obs, INSTRUM='NOMIC'
  
  ; --- Now update fits file
  files = FILE_SEARCH(pth.l0fits_path + date_obs + pth.sep + 'bckg' + pth.sep, '*N' + STRING(nod, FORMAT='(I03)') + '*IMG.fits', COUNT=n_files)
  FOR i=0, n_files-1 DO BEGIN
   hdr = HEADFITS(files[i])
   SXADDPAR,hdr,'XCEN_SX',pos_sx[0]
   SXADDPAR,hdr,'YCEN_SX',pos_sx[1]
   SXADDPAR,hdr,'XCEN_DX',pos_dx[0]
   SXADDPAR,hdr,'YCEN_DX',pos_dx[1]
   MODFITS,files[i],0,hdr
   data = MRDFITS(files[i], 1, /SILENT)
   data.xcen_sx[*] = pos_sx[0]
   data.ycen_sx[*] = pos_sx[1]
   data.xcen_dx[*] = pos_dx[0]
   data.ycen_dx[*] = pos_dx[1]
   MODFITS,files[i],data, /exten
  ENDFOR
  
  ; --- Now print log
  logfile =  pth.l0fits_path + date_obs + pth.sep + 'changelog.txt'
  IF NOT FILE_TEST(logfile) THEN OPENW, log,  logfile, /GET_LUN, WIDTH=1000 $
                            ELSE OPENW, log, logfile, /GET_LUN, /APPEND
  PRINTF, log, SYSTIME() + ' ; NOD ' + STRING(nod, FORMAT='(I0)') + ' (UTC ', utc, ') -- new position : ' + STRING(pos_sx[0], FORMAT='(I0)') + ',' + STRING(pos_sx[1], FORMAT='(I0)') + '--' + STRING(pos_dx[0], FORMAT='(I0)') + ',' + STRING(pos_dx[1], FORMAT='(I0)')
  CLOSE, log
  FREE_LUN, log

END