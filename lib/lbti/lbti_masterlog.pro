; +
; NAME: LBTI_MASTERLOG
; 
; PURPOSE:
;   Read header of all FITS files located in data_path and create a summary file containing a line for each of these files.
;   This "masterlog" file is saved in data_path as well as a config ID file, which summarizes the different instrumental configurations. 
;   WARNING: running this routine on a large amount of files is quite time consuming. 
;
; INPUTS:
;   data_path  :  String vector with the path to the fits files with data.
;
; KEYWORDS
;   APPEND     :  If set, the code will append new lines to the masterlog file if it already exists
;   IDL        :  If set, force the use of IDL to compute the masterlog
;   
; CAVEATS
;   The nod_id computation doesn't work for NOMIC if RA/DEC offsets are used instead of DETXY.
;   This is because, for nulling, we take data between nods and the offset and AO keywords are changed at different times.
;
; MODIFICATION HISTORY:
;   Version 1.0,  14-JAN-2014, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  16-JAN-2014, DD: Added Time stamp to the masterlofg file
;   Version 1.2,  08-FEB-2014, DD: Improved speed
;   Version 1.3,  10-FEB-2014, DD: Added keyword append
;   Version 1.4,  11-SEP-2014, DD: Cleaned code and added /NOCONTINUE to SXPAR calls
;   Version 1.5,  16-SEP-2014, DD: Added call to C routine which is the default now (IDL slower!)
;   Version 1.6,  23-SEP-2014, DD: Added keyword IDL
;   Version 1.7,  23-OCT-2014, DD: Minor bug corrected
;   Version 1.8,  11-NOV-2014, DD: Improved robustness to dead telescope server
;   Version 1.9,  29-OCT-2015, DD: Added config ID number and file
;   Version 2.0,  01-DEC-2015, DD: Added pointing ID and elevation to masterlog file (also implemented a nod threshold)
;   Version 2.1,  02-DEC-2015, DD: Replaced FXPAR by SXPAR which is ~15% faster
;   Version 2.2,  20-DEC-2015, DD: Minor bug corrected
;   Version 2.3,  22-DEC-2015, DD: Now handle cases where the header is corrupted
;   Version 2.4,  23-DEC-2015, DD: Updated to use new AO offset keywords
;   Version 2.5,  20-JAN-2016, DD: Implemented a more clever way to determine whether the telescope has been nodded
;   Version 2.6,  20-FEB-2016, DD: Now doesn't consider AO keywords if NOMIC frames (only DET XY are considered to compute the nod ID). 
;   Version 2.7,  23-FEB-2016, DD: Updated to deal with AO loop status of 2 (which means AOPause)
;   Version 2.8,  22-MAR-2016, DD: Fixed bug with config_id file when /APPEND is set 
;   Version 2.9,  27-MAR-2016, DD: Now force DATATYPE to 1 when NOMIC's filter wheel 2 is 'Blank+tape'
;   Version 3.0,  11-MAY-2016, DD: Now increase nod ID if OBSTYPE is modified
;   Version 3.1,  12-MAY-2016, RMG: Added keyword ALT_OUT_DIR
;   Version 3.2,  26-MAY-2016, DD: Added a column with mean background level
;   Version 3.3,  05-OCT-2016, DD: Replaced STRTRIM by STRCOMPRESS to remove middle white spaces too
;   Version 3.4,  08-OCT-2016, DD: Nod ID and pointing ID are now starting at 1 (instead of 0)
;   Version 3.5,  16-OCT-2016, DD: Added nod_id log file + use AO status keywords 
;   Version 3.6,  23-MAR-2017, DD: Now flag non-linear frames as BAD 
;   Version 3.7,  16-MAY-2017, DD: Implemented cfid_file.sav to be comptaible with real-time running
;   Version 3.7,  23-MAY-2017, DD: Now compute central value over ROI
;   Version 3.8,  23-NOV-2017, DD: Corrected minor bug for large number of coadds
;   Version 3.9,  16-APR-2018, DD: Improved robustness against bad header inputs

PRO LBTI_MASTERLOG, data_path, APPEND=append, IDL=idl, ALT_OUT_DIR=alt_out_dir
  
  ; Hard-coded paramaters
  nod_thre = 0.5    ; In arcsec, this is the threshold to assign a new nod number
  width    = 138    ; Number of columns in the output log      
        
  ; C version is currently obsolete
  IF NOT KEYWORD_SET(IDL) THEN idl = 1
  
  ; Retrieve FITS files in the input directory
  data_files = FILE_SEARCH(data_path,'*.fits', COUNT=n_files) 
  IF n_files LT 1 THEN MESSAGE, 'Input data path empty, skipping it.' 
  
  ; Masterlog, config ID, and nod ID files
  IF NOT KEYWORD_SET(ALT_OUT_DIR) THEN log_path = data_path ELSE log_path = alt_out_dir
  mlog_file = log_path + 'masterlog.dat'
  cfid_file = log_path + 'config_id.dat'
  cfid_sav  = log_path + 'config_id.sav'
  nid_file  = log_path + 'nod_id.dat' 	
  
  ; If the C library is present, go for it! Otherwise, use much slower IDL.
  c_path = 'nodrs/c/'
  IF FILE_TEST(c_path + 'lbti_masterlog.so') AND NOT KEYWORD_SET(IDL) THEN BEGIN
      
      ; Make sure the file array has the right format
      IF (SIZE(data_files,/TYPE) NE 7) THEN MESSAGE, 'ERROR: data_files must be an array of strings.'
      ; Prepare output arrays for C
      iflen   = STRLEN(data_files)
      n_files = SIZE(data_files,/n_elements)
      ifn     = n_files[0]
      strloc  = 0
      ; Arrange the file array for c
      inputfile4c =BYTARR(TOTAL(iflen) + ifn)
      FOR i=0,ifn-1 DO BEGIN
        inputfile4c[strloc:strloc+iflen[i]-1] = BYTE(data_files[i])
        inputfile4c[strloc+iflen[i]] = BYTE(0)
        strloc += (iflen[i] + 1)
     ENDFOR
     ; Prepare masterlog file (must have a null string at the end to be properly interpreted by C)
     mlog_byte = [BYTE(mlog_file),BYTE(0)]
     ; Call c routine
     tmp  = CALL_EXTERNAL(c_path + 'lbti_masterlog.so', 'lbti_masterlog', mlog_byte, inputfile4c, n_files, 4L)
     
  ENDIF ELSE BEGIN
  
    ; If append is set, read the existing file to find the latest entry
    ; Al three log files must exist. Otherwise, start from 0 (ensure backward compatibility).
    IF KEYWORD_SET(APPEND) AND FILE_TEST(mlog_file) AND FILE_TEST(cfid_file) AND FILE_TEST(nid_file) THEN BEGIN
      ; 1. Read masterlog file
      READ_TABLE, mlog_file, file_id, time_obs, pt_id, nod_id, FIRST=2, SKIP=[1,3,4,5,7,9,10,11,12,13,14,15], STRING_ARRAY=[1], SEP=';'      
      id_min   = MAX(file_id)+1
      IF id_min GE N_ELEMENTS(data_files) THEN RETURN
      start_id = MIN(file_id[WHERE(nod_id EQ MAX(nod_id))])
      nod_id   = MAX(nod_id)
      pt_id    = MAX(pt_id)
      OPENW, lun, mlog_file, /GET_LUN, /APPEND, WIDTH=width
      header = HEADFITS(data_files[id_min-1], /SILENT)
      ; 2. Read config_ID file
      IF FILE_TEST(cfid_sav) THEN RESTORE, cfid_sav ELSE READ_TABLE, cfid_file, cfg_id, wav_id, bdw_id, dit_id, ncoadd_id, smplmode_id, pagain_id, detmode_id, pabandw_id, pactcdly_id, FIRST=2, STRING_ARRAY=[3], SEP=';'
      OPENW, lun2, cfid_file, /GET_LUN, /APPEND, WIDTH=width
      ; 3. Open NOD ID file
      OPENW, lun3, nid_file, /GET_LUN, /APPEND, WIDTH=width  ; Open file for append new confid IDs
    ENDIF ELSE BEGIN
      ; 1. Create masterlog file
      OPENW, lun, mlog_file, /GET_LUN, WIDTH=width
      PRINTF,lun, 'Masterlog created on ' + SYSTIME() + ' by LBTI_DRS pipeline.'
      PRINTF,lun, 'ID;name;time_obs;elevation;lam_eff[m];int_time[s];PT_ID;CFG_ID;NOD_ID;CHP_ID;OBSTYPE;DATATYPE;CV;N_XPIX;N_YPIX;N_FRAME;STATUS '
      id_min   = 0  ; First file to read is 0
      nod_id   = 1  ; First nod will be given number 1
      pt_id    = 1  ; First pointing will be given number 1
      start_id = 0  ; First file
      header   = HEADFITS(data_files[id_min], /SILENT)
      ;2. Create config ID file
      OPENW, lun2, cfid_file, /GET_LUN, WIDTH=width
      PRINTF,lun2, 'Config ID file created on ' + SYSTIME()
      PRINTF,lun2, 'CFG_ID;lam_eff[m];bndw[m];int_time[s];ncoadd;smplmode;pagaini;detmode;pabandwi;pactcdly'
      cfg_id = 0 & wav_id=0 & bdw_id = 0 & dit_id = 0 & ncoadd_id = 0 & smplmode_id = 0 & pagain_id = 0 & detmode_id = 0 & pabandw_id = 0 & pactcdly_id = 0
      ; 3. Create NOD ID file
      OPENW, lun3, nid_file, /GET_LUN, WIDTH=width
      PRINTF,lun3, 'Nod ID file created on ' + SYSTIME()
      PRINTF,lun3, 'NOD_ID;PT_ID;name;time_obs;elevation;lam_eff[m];int_time[s];OBSTYPE;N_XPIX;N_YPIX;START_ID;END_ID
    ENDELSE
             
    ; Init nodding position variables (read first file)
    ; We have to check both the LBT_**OS keywords, which is updated only for DETXY offsets, and *OFFSET* keywords, which are updated for both DETXY and DETRADEC
    ; but can remain below the threshold if the nod failed.
    chp_id    = 0
    nod_id0   = nod_id
    pt_id0    = pt_id
    lbt_lxao0 = FLOAT(SXPAR(header, 'LOFFSETX', /NOCONTINUE)) 
    lbt_lyao0 = FLOAT(SXPAR(header, 'LOFFSETY', /NOCONTINUE)) 
    lbt_rxao0 = FLOAT(SXPAR(header, 'ROFFSETX', /NOCONTINUE))  
    lbt_ryao0 = FLOAT(SXPAR(header, 'ROFFSETY', /NOCONTINUE)) 
    lbt_lxos0 = FLOAT(SXPAR(header, 'LBT_LXOS', /NOCONTINUE))
    lbt_lyos0 = FLOAT(SXPAR(header, 'LBT_LYOS', /NOCONTINUE))
    lbt_rxos0 = FLOAT(SXPAR(header, 'LBT_RXOS', /NOCONTINUE))
    lbt_ryos0 = FLOAT(SXPAR(header, 'LBT_RYOS', /NOCONTINUE))
    n_x0      = FIX(SXPAR(header, 'NAXIS1', /NOCONTINUE))
    n_y0      = FIX(SXPAR(header, 'NAXIS2', /NOCONTINUE))
    time_obs0 = SXPAR(header, 'TIME-OBS', /NOCONTINUE)
    lbt_alt0  = DOUBLE(SXPAR(header, 'LBT_ALT', /NOCONTINUE) < 90)
    int_time0 = DOUBLE(SXPAR(header, 'EXPTIME', /NOCONTINUE))
    objname0  = STRCOMPRESS(STRING(SXPAR(header, 'OBJNAME', /NOCONTINUE)), /REMOVE_ALL)
    obstype0  = SXPAR(header, 'OBSTYPE', /NOCONTINUE) & IF STRTRIM(STRING(obstype0)) EQ 'object' THEN obstype = 0 ELSE obstype0 = FIX(obstype0) ; backwards compatibility
    
    ; Loop over the files and append to masterlog file
    FOR i_f=id_min, n_files-1 DO BEGIN
      
      ; Read header and check integrity
      img = READFITS(data_files[i_f], header, /SILENT)
      IF (SIZE(img))[0] EQ 0 THEN GOTO, SKIP_FRAME
      
      ; Compute central value over ROI
      instrum  = STRUPCASE(STRTRIM(SXPAR(header, 'INSTRUME'), 2))
      subsectm = FXPAR(header, 'SUBSECTM')
      IF !err NE -1 THEN roi = LBTI_ROI(subsectm, INSTRUM=instrum) ELSE roi = [0,0,N_ELEMENTS(img[*,0,0])-1,N_ELEMENTS(img[0,*,0])-1]
      IF MAX(roi) EQ 0 THEN roi =  [0,0,N_ELEMENTS(img[*,0,0])-1,N_ELEMENTS(img[0,*,0])-1]
      cv = MEDIAN(img[roi[0]:roi[2],roi[1]:roi[3]])
      
      ; Read header information
      file_id  = LONG(STREGEX(STREGEX(data_files[i_f],'_[0-9]+.fits',/EXTRACT),'[0-9]+',/EXTRACT)) 
      objname  = STRCOMPRESS(STRING(SXPAR(header, 'OBJNAME', /NOCONTINUE)), /REMOVE_ALL)
      time_obs = SXPAR(header, 'TIME-OBS', /NOCONTINUE)
      n_coadds = FIX(SXPAR(header, 'NCOADDS', /NOCONTINUE))
      smplmode = FIX(SXPAR(header, 'SMPLMODE', /NOCONTINUE))
      pagain   = FIX(SXPAR(header, 'PAGAINI', /NOCONTINUE))
      detmode  = FIX(SXPAR(header, 'DETMODE', /NOCONTINUE))
      pabandw  = FIX(SXPAR(header, 'PABANDWI', /NOCONTINUE))
      pactcdly = FIX(SXPAR(header, 'PACTCDLY', /NOCONTINUE))
      lbt_alt  = DOUBLE(SXPAR(header, 'LBT_ALT', /NOCONTINUE) < 90)    ; Failsafe in case 999
      lbt_lxao = FLOAT(SXPAR(header, 'LOFFSETX', /NOCONTINUE)) 
      lbt_lyao = FLOAT(SXPAR(header, 'LOFFSETY', /NOCONTINUE)) 
      lbt_rxao = FLOAT(SXPAR(header, 'ROFFSETX', /NOCONTINUE))  
      lbt_ryao = FLOAT(SXPAR(header, 'ROFFSETY', /NOCONTINUE)) 
      lbt_lxos = FLOAT(SXPAR(header, 'LBT_LXOS', /NOCONTINUE))
      lbt_lyos = FLOAT(SXPAR(header, 'LBT_LYOS', /NOCONTINUE))
      lbt_rxos = FLOAT(SXPAR(header, 'LBT_RXOS', /NOCONTINUE))
      lbt_ryos = FLOAT(SXPAR(header, 'LBT_RYOS', /NOCONTINUE))
      flag     = STRTRIM(STRING(SXPAR(header, 'FLAG', /NOCONTINUE)), 2)
      IF !err NE -1 THEN BEGIN
        CASE flag OF
          'SCI': datatype = 0
          'DRK': datatype = 1
          'FLT': datatype = 2
          'CAL': datatype = 3
          'BAD': datatype = -1
          ELSE : datatype = 0
        ENDCASE
      ENDIF ELSE datatype = FIX(SXPAR(header, 'DATATYPE'))
      int_time = DOUBLE(SXPAR(header, 'EXPTIME', /NOCONTINUE))
      obstype  = SXPAR(header, 'OBSTYPE', /NOCONTINUE)          & IF STRTRIM(STRING(obstype)) EQ 'object' THEN obstype = 0 ELSE obstype = FIX(obstype) ; will be populated later
      n_x      = FIX(SXPAR(header, 'NAXIS1', /NOCONTINUE))
      n_y      = FIX(SXPAR(header, 'NAXIS2', /NOCONTINUE))
      n_fr     = FIX(SXPAR(header, 'NAXIS3', /NOCONTINUE))      & IF !err EQ -1 THEN n_fr = 1
        
      ; Derive wavelength
      CASE instrum OF
        'NOMIC'  : BEGIN
            nom_fw1  = STRTRIM(SXPAR(header, 'NOMICFW1', /NOCONTINUE), 2) & IF !err EQ -1 THEN nom_fw1  = STRTRIM(SXPAR(header, 'NOM_FW1'), 2)
            nom_fw2  = STRTRIM(SXPAR(header, 'NOMICFW2', /NOCONTINUE), 2) & IF !err EQ -1 THEN nom_fw2  = STRTRIM(SXPAR(header, 'NOM_FW2'), 2)
            ; Force DATATYPE to 1 if  nom_fw2 is 'Blank+tape' (in some old data, FLAG is not correctly set)
            IF nom_fw2 EQ 'Blank+tape' THEN BEGIN
              datatype = 1
              flag     = 'DRK'
            ENDIF
          END
        'LMIRCAM': BEGIN
            lmir_fw1 = STRTRIM(SXPAR(header, 'LMIR_FW1', /NOCONTINUE), 2) & IF !err EQ -1 THEN lmir_fw1 = STRTRIM(SXPAR(header, 'LMIR_FW1'), 2)
            lmir_fw2 = STRTRIM(SXPAR(header, 'LMIR_FW2', /NOCONTINUE), 2) & IF !err EQ -1 THEN lmir_fw2 = STRTRIM(SXPAR(header, 'LMIR_FW2'), 2)
            lmir_fw3 = STRTRIM(SXPAR(header, 'LMIR_FW3', /NOCONTINUE), 2) & IF !err EQ -1 THEN lmir_fw3 = STRTRIM(SXPAR(header, 'LMIR_FW3'), 2)
            lmir_fw4 = STRTRIM(SXPAR(header, 'LMIR_FW4', /NOCONTINUE), 2) & IF !err EQ -1 THEN lmir_fw4 = STRTRIM(SXPAR(header, 'LMIR_FW4'), 2)
           END
        ELSE: BEGIN
            lam_cen = 0 & bdw = 0 & flag = 'B'
            GOTO, SKIP_FRAME
          END
      ENDCASE
      LBTI_WAVEBAND, lam_cen, bdw, INSTRUM=instrum, LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1 = nom_fw1, NOM_FW2 = nom_fw2 ; derive wavelength and bandwidth based on the header
             
      ; Flag non-linear frames
      GET_CNF, cnf, INSTRUM=instrum 
      interface = STRLOWCASE(STRTRIM(SXPAR(header, 'INTERFAC'), 2))
      IF interface EQ 'preamp' THEN sat_lim = cnf.max_adu[1] ELSE sat_lim = cnf.max_adu[0]
      IF cv GE LONG(n_coadds)*sat_lim THEN flag = 'BAD'
         
      ; Compute pointing ID and nod ID numbers (increase nod ID if new pointing, i.e. new OBJNAME)
      ; If same OBJNAME, a new nod number is assigned if one of the telescopes has moved by more than nod_thre
      ; If flag is BAD, don't do anything here
      IF flag NE 'BAD' THEN BEGIN
        IF NOT STRMATCH(objname, objname0, /FOLD_CASE) THEN BEGIN
          pt_id  += 1
          nod_id += 1
        ENDIF ELSE BEGIN
          ; Increase nod ID if OBSTYPE has changed (but ignore an obstype of 4)
          ; If not, look at the AO nodding keywords
          IF obstype EQ obstype0 THEN BEGIN
            ; Ignore the file if offset fields report 999 (i.e., dead telescope server)
            IF lbt_lxos NE 999 AND lbt_lyos NE 999 AND lbt_rxos NE 999 AND lbt_ryos NE 999 THEN BEGIN
              ; Here we check whether an XY offset has been sent to the telescope (default for nulling).
              ; If not, we check whether there was an AO offset (Bayside stages). If yes, this is considered a new nod only if the XY telescope offset keywords haven't changed and no nulling data for which we only use DET XY.
              ; If they have changed, this means that it wasn't an RA-DEC telescope offset (for which we don't have keyword) but probably a failed nod.
              IF ABS(lbt_lxos-lbt_lxos0) GT nod_thre OR ABS(lbt_lyos-lbt_lyos0) GT nod_thre OR ABS(lbt_rxos-lbt_rxos0) GT nod_thre OR ABS(lbt_ryos-lbt_ryos0) GT nod_thre THEN nod_id += 1 $
              ELSE IF ABS(lbt_lxao-lbt_lxao0) GT nod_thre OR ABS(lbt_lyao-lbt_lyao0) GT nod_thre OR ABS(lbt_rxao-lbt_rxao0) GT nod_thre OR ABS(lbt_ryao-lbt_ryao0) GT nod_thre THEN $
                IF lbt_lxos EQ lbt_lxos0 AND lbt_lyos EQ lbt_lyos0 AND lbt_rxos EQ lbt_rxos0 AND lbt_ryos EQ lbt_ryos0 AND instrum NE 'NOMIC' THEN nod_id += 1
            ENDIF
          ENDIF ELSE nod_id += 1  
        ENDELSE
      ENDIF    
      
      ; Derive chop ID (not yet implemented)
      ; chp_id = 0
      
      ; Set loop status flag (for obstype 0,1,2 only)
      IF flag NE 'BAD' THEN BEGIN
        IF obstype NE 3 THEN BEGIN
          flag     = 'C' ; Assume frame closed by default
          ; LLoopOn and RLoopOn cannot be trusted, use also status
          IF FIX(SXPAR(header, 'RLOOPON', /NOCONTINUE)) NE 0 OR STRCOMPRESS(SXPAR(header, 'RSTATUS', /NOCONTINUE), /REMOVE_ALL) EQ 'AORunning' THEN dloopon = 1 ELSE dloopon = 0 & IF !err EQ -1 THEN dloopon = 1  ; Assume loop is closed if not present
          IF FIX(SXPAR(header, 'LLOOPON', /NOCONTINUE)) NE 0 OR STRCOMPRESS(SXPAR(header, 'LSTATUS', /NOCONTINUE), /REMOVE_ALL) EQ 'AORunning' THEN sloopon = 1 ELSE sloopon = 0 & IF !err EQ -1 THEN sloopon = 1  ; Assume loop is closed if not present
          pcclosed = FIX(SXPAR(header, 'PCCLOSED', /NOCONTINUE)) & IF !err EQ -1 THEN plc_status = 1 
          IF obstype EQ 0 THEN BEGIN
            IF sloopon NE 1 AND dloopon NE 1 THEN flag = 'O'  ; Open-loop frame
            IF sloopon NE 1 AND dloopon EQ 1 THEN flag = 'R'  ; Right side only
            IF sloopon EQ 1 AND dloopon NE 1 THEN flag = 'L'  ; Left side only
          ENDIF ELSE IF sloopon NE 1 OR dloopon NE 1 OR pcclosed EQ 0 THEN flag = 'O' ; Open-loop frame
        ENDIF ELSE flag = 'O' ; Frame open
      ENDIF ELSE flag = 'B'
            
      ; Compute config ID and append to config file if new
      idx_id   = WHERE(lam_cen EQ wav_id AND bdw EQ bdw_id AND int_time EQ dit_id AND n_coadds EQ ncoadd_id AND smplmode EQ smplmode_id AND pagain EQ pagain_id AND detmode EQ detmode_id AND pabandw EQ pabandw_id AND pactcdly EQ pactcdly_id, n_id)
      IF n_id LT 1 THEN BEGIN
        ; Create new entry to the config id list
        new_id = MAX(cfg_id + 1)
        cfg_id = [cfg_id,new_id] & wav_id = [wav_id,lam_cen] & bdw_id = [bdw_id,bdw] & dit_id = [dit_id,int_time] & ncoadd_id = [ncoadd_id,n_coadds] & smplmode_id = [smplmode_id,smplmode] 
        pagain_id = [pagain_id,pagain] & detmode_id = [detmode_id,detmode] & pabandw_id = [pabandw_id,pabandw] & pactcdly_id = [pactcdly_id,pactcdly]
        ; Save
        SAVE, cfg_id, wav_id, bdw_id, dit_id, ncoadd_id, smplmode_id, pagain_id, detmode_id, pabandw_id, pactcdly_id, FILENAME=cfid_sav    
        ; Append new config to config file
        PRINTF, lun2, STRING(new_id, FORMAT='(I0)'), ';', STRING(lam_cen, FORMAT='(E8.2)'),';', STRING(bdw, FORMAT='(E8.2)'), ';', STRING(int_time, FORMAT='(E8.2)'), ';', STRING(n_coadds, FORMAT='(I0)'), ';', $
                STRING(smplmode, FORMAT='(I0)'), ';', STRING(pagain, FORMAT='(I0)'), ';', STRING(detmode, FORMAT='(I0)'), ';', STRING(pabandw, FORMAT='(I0)'), ';', STRING(pactcdly, FORMAT='(I0)')
      ENDIF ELSE new_id = cfg_id[idx_id]
      
      ; If new NOD, write nod_id file and reset nod variables
      ; Also print to display to keep user entertained
      IF nod_id NE nod_id0 THEN BEGIN
        PRINT, 'Nod ID ' + STRING(nod_id0, FORMAT='(I0)') + ' : files ' + STRING(start_id, FORMAT='(I0)') + ' to ' + STRING(file_id-1, FORMAT='(I0)') + ' (' + STRING(100*FLOAT(file_id)/n_files, FORMAT='(F4.1)') + '%)' 
        PRINTF, lun3, STRING(nod_id0, FORMAT='(I0)'), ';', STRING(pt_id0, FORMAT='(I0)'), ';', STRCOMPRESS(STRING(objname0), /REMOVE_ALL), ';', STRCOMPRESS(STRING(time_obs0), /REMOVE_ALL), ';', STRING(lbt_alt0, FORMAT='(F6.3)'), ';', $
                      STRING(int_time0, FORMAT='(E8.2)'), ';', STRING(obstype0, FORMAT='(I0)'), ';', STRING(n_x0, FORMAT='(I0)'), ';', STRING(n_y0, FORMAT='(I0)'), ';', STRING(start_id, FORMAT='(I0)'), ';', STRING(file_id-1, FORMAT='(I0)')
        start_id  = file_id  
        nod_id0   = nod_id & pt_id0 = pt_id & objname0 = objname & time_obs0 = time_obs & lbt_alt0 = lbt_alt & int_time0 = int_time & obstype0 = obstype & n_x0 = n_x & n_y0 = n_y
        lbt_lxos0 = lbt_lxos & lbt_lyos0 = lbt_lyos & lbt_rxos0 = lbt_rxos & lbt_ryos0 = lbt_ryos
        lbt_lxao0 = lbt_lxao & lbt_lyao0 = lbt_lyao & lbt_rxao0 = lbt_rxao & lbt_ryao0 = lbt_ryao   
      ENDIF
      
      ; Append to masterlog file
      PRINTF, lun, STRING(file_id, FORMAT='(I0)'), ';', STRCOMPRESS(STRING(objname), /REMOVE_ALL), ';', STRCOMPRESS(STRING(time_obs), /REMOVE_ALL), ';', STRING(lbt_alt, FORMAT='(F6.3)'), ';', STRING(lam_cen, FORMAT='(E8.2)'), ';', $
                   STRING(int_time, FORMAT='(E8.2)'), ';', STRING(pt_id, FORMAT='(I0)'), ';', STRING(new_id, FORMAT='(I0)'), ';', STRING(nod_id, FORMAT='(I0)'), ';', STRING(chp_id, FORMAT='(I0)'), ';', STRING(obstype, FORMAT='(I0)'), ';',$
                   STRING(datatype, FORMAT='(I0)'), ';', STRING(cv, FORMAT='(I0)'), ';', STRING(n_x, FORMAT='(I0)'), ';', STRING(n_y, FORMAT='(I0)'), ';', STRING(n_fr, FORMAT='(I0)'), ';', flag
                    
      ; Jump point
      SKIP_FRAME:
    ENDFOR
    
    ; When done, write the last line of nod_id file
    IF N_ELEMENTS(file_id) EQ 1 THEN BEGIN
      PRINT, 'Nod ID ' + STRING(nod_id0, FORMAT='(I0)') + ' : files ' + STRING(start_id, FORMAT='(I0)') + ' to ' + STRING(file_id-1, FORMAT='(I0)') + ' (' + STRING(100*FLOAT(file_id)/n_files, FORMAT='(F4.1)') + '%)' 
      PRINTF, lun3, STRING(nod_id0, FORMAT='(I0)'), ';', STRING(pt_id0, FORMAT='(I0)'), ';', STRCOMPRESS(STRING(objname0), /REMOVE_ALL), ';', STRCOMPRESS(STRING(time_obs0), /REMOVE_ALL), ';', STRING(lbt_alt0, FORMAT='(F6.3)'), ';', $
                    STRING(int_time0, FORMAT='(E8.2)'), ';', STRING(obstype0, FORMAT='(I0)'), ';', STRING(n_x0, FORMAT='(I0)'), ';', STRING(n_y0, FORMAT='(I0)'), ';', STRING(start_id, FORMAT='(I0)'), ';', STRING(file_id-1, FORMAT='(I0)')
    ENDIF
    
    ; Close masterlog and config ID files
    CLOSE, lun  & FREE_LUN, lun
    CLOSE, lun2 & FREE_LUN, lun2
    CLOSE, lun3 & FREE_LUN, lun3
  ENDELSE ; Endelse on IDL or C
END
