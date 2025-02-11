pro LBTI_UPDATEL0BEAM, date, nod, pos_sx, pos_dx
  ;
  ; PURPOSE:
  ; Update/force the beam positon both in the data log and in the FITS files.
  ;
  ; INPUTS:
  ; date   :  Date of the data
  ; nod    :  Nod ID to be updated
  ; pos_sx :  Position of the left beam (x, y)
  ; pos_dx :  Position of the left beam (x, y)
  ;
  ; MODIFICATION HISTORY:
  ; Version 1.0, 12-JAN-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
  ; Version 1.1, 11-FEB-2025, DD: added date of backup to the backup datalog file
  compile_opt idl2

  if not keyword_set(pos_dx) then pos_dx = [0, 0]

  ; Declare path
  DECLARE_PATH, pth, instrum = 'NOMIC'

  ; --- Long-format date
  date_obs = '20' + strmid(date, 0, 2) + '-' + strmid(date, 2, 2) + '-' + strmid(date, 4, 2)

  ; --- Restore LO log file and read useful data
  datalog = pth.l0Fits_path + date_obs + pth.sep + 'datalog.sav'
  if not file_test(datalog) then begin
    message, 'No data log file found!', /continue
    RETURN
  endif else restore, datalog

  ; -- Update position
  idx_nod = where(data_r.nod_id eq nod, n_nod)
  if n_nod le 0 then message, 'No corresponding nod!'

  ; --- If ok, make a copy
  caldat, julday(), m, d, y, hh, mm, ss
  yymmdd = strcompress(string(y) + string(m, format = '(I2.2)') + string(d, format = '(I2.2)'), /remove_all) ; Concatenate into "YYYYMMDD"
  datalog_bu = pth.l0Fits_path + date_obs + pth.sep + 'backup' + pth.sep + '/datalog_bu-' + yymmdd + '.sav'
  if not file_test(datalog_bu) then file_copy, datalog, datalog_bu
  utc = data_r[idx_nod].ut_time
  xcen_sx = data_r[idx_nod].xcen_sx
  ycen_sx = data_r[idx_nod].ycen_sx
  xcen_dx = data_r[idx_nod].xcen_dx
  ycen_dx = data_r[idx_nod].ycen_dx
  data_r[idx_nod].xcen_sx = pos_sx[0]
  data_r[idx_nod].ycen_sx = pos_sx[1]
  data_r[idx_nod].xcen_dx = pos_dx[0]
  data_r[idx_nod].ycen_dx = pos_dx[1]

  ; --- Save file
  save, data_r, filename = datalog

  ; --- Backup and write datalog txt file
  datalog = pth.l0Fits_path + date_obs + pth.sep + 'datalog.txt'
  datalog_bu = pth.l0Fits_path + date_obs + pth.sep + 'backup' + pth.sep + 'datalog_bu-' + yymmdd + '.txt'
  if not file_test(datalog_bu) then file_copy, datalog, datalog_bu
  WRITE_DATALOG, date_obs, instrum = 'NOMIC'

  ; --- Now update fits file
  files = file_search(pth.l0Fits_path + date_obs + pth.sep + 'bckg' + pth.sep, '*N' + string(nod, format = '(I03)') + '*IMG.fits', count = n_files)
  for i = 0, n_files - 1 do begin
    hdr = HEADFITS(files[i])
    SXADDPAR, hdr, 'XCEN_SX', pos_sx[0]
    SXADDPAR, hdr, 'YCEN_SX', pos_sx[1]
    SXADDPAR, hdr, 'XCEN_DX', pos_dx[0]
    SXADDPAR, hdr, 'YCEN_DX', pos_dx[1]
    MODFITS, files[i], 0, hdr
    data = MRDFITS(files[i], 1, /silent)
    data.xcen_sx[*] = pos_sx[0]
    data.ycen_sx[*] = pos_sx[1]
    data.xcen_dx[*] = pos_dx[0]
    data.ycen_dx[*] = pos_dx[1]
    MODFITS, files[i], data, /exten
  endfor

  ; --- Now print log
  logfile = pth.l0Fits_path + date_obs + pth.sep + 'changelog.txt'
  if not file_test(logfile) then openw, log, logfile, /get_lun, width = 1000 $
  else openw, log, logfile, /get_lun, /append
  printf, log, systime() + ' ; NOD ' + string(nod, format = '(I0)') + ' (UTC ', utc, ') -- old position : ' + string(xcen_sx, format = '(I0)') + ',' + string(ycen_sx, format = '(I0)') + '-' + string(xcen_dx, format = '(I0)') + ',' + string(ycen_dx, format = '(I0)') + ' -- new position : ' + string(pos_sx[0], format = '(I0)') + ',' + string(pos_sx[1], format = '(I0)') + '-' + string(pos_dx[0], format = '(I0)') + ',' + string(pos_dx[1], format = '(I0)')
  close, log
  free_lun, log
end
