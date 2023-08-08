;+
; NAME: LBTI_FITSUTIL
;
; PURPOSE:
;   This procedure edits the header of the raw L0 files of a given date.
;   If set, the new L0 file header will be modified with the value(s) of the input keywords.
;
; MADATORY INPUT
;   date          : The date of the files (6-digit format)
;   data_idx      : The file number range to process
;
; KEYWORDS
;   LMIRCAM       : Set for LMIRCam
;   NOMIC         : Set for NOMIC
;
; MODIFICATION HISTORY:
;   Version 1.0,  08-JUL-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  08-OCT-2013, DD: added DATATYP, OBSTYPE, and OBJNAME
;   Version 1.2,  13-FEB-2016, DD: updated for new filter wheel name
;   Version 1.3,  19-APR-2016, DD: added keywords RLOOPON and LLOOPON
;   Version 1.4,  06-OCT-2016, DD: added keyword PID

PRO LBTI_FITSUTIL, date, data_idx, OBJNAME=objname, LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, $
                   NOM_APW=nom_apw, PHA_FW1=pha_fw1, PHA_FW2=pha_fw2, PHA_IMG=pha_img, NIL_DIC=nil_dic,$
                   LBT_LXOS=lbt_lxos, LBT_LYOS=lbt_lyos, LBT_RXOS=lbt_rxos, LBT_RYOS=lbt_ryos, DATATYPE=datatype, OBSTYPE=obstype, FLAG=flag, LMIRCAM=lmircam, CROP=crop, $
                   LLOOPON=lloopon, RLOOPON=rloopon, PID=pid
                   

; Declare path
IF KEYWORD_SET(LMIRCAM) THEN instrum = 'lmircam' ELSE instrum = 'nomic'
DECLARE_PATH, pth, INSTRUM=instrum

; Retrieve FITS files in the input directory
path = pth.root_data
data_files = FILE_SEARCH(path + date,'*.fits') & n_files = N_ELEMENTS(data_files)
IF n_files LT 1 THEN MESSAGE, 'Input data path empty'

; Only reads the file in the 'data_idx' range. For some reason, it's not always working without the LONG!!!!
; Warning, this only works assuming the LBTI file definition (i.e., 'n_130424_000123.fits')
data_nbr = LONG(STREGEX(STREGEX(data_files,'_[0-9]+.fits',/EXTRACT),'[0-9]+',/EXTRACT))
IF KEYWORD_SET(DATA_IDX) THEN idx_out = WHERE(data_nbr GE data_idx[0] AND data_nbr LE data_idx[1], n_files) ELSE idx_out = INDGEN(n_files)

; Loop over the files
FOR i_f=0, n_files-1 DO BEGIN
  ; Read file
  img_file = READFITS(data_files[idx_out[i_f]], hdr, /SILENT)
 
  IF KEYWORD_SET(CROP) THEN BEGIN
    ; Derive image size
    n_img  = N_ELEMENTS(img_file[0,0,*])
    n_xpix = N_ELEMENTS(img_file[*,0,0])
    n_ypix = N_ELEMENTS(img_file[0,*,0])
    ; Bias correct the image if it's a full frame
    IF KEYWORD_SET(LMIRCAM) AND n_xpix EQ 1024 AND n_ypix EQ 1024 THEN img_file = LBTI_BIASSUBTRACT(TEMPORARY(img_file), /LMIRCAM)
    ; Crop the frame
    IF n_img GT 1 THEN img_file = EXTRAC(TEMPORARY(img_file), crop[0], crop[1], 0, crop[2], crop[3], n_img) ELSE img_file = EXTRAC(TEMPORARY(img_file), crop[0], crop[1], crop[2], crop[3])
    ; Parse header with the new image size
    SXADDPAR, hdr, 'NAXIS1', crop[2]
    SXADDPAR, hdr, 'NAXIS2', crop[3]
  ENDIF
  
  ; Parse new header information
  IF KEYWORD_SET(OBJNAME)  THEN SXADDPAR, hdr, 'OBJNAME',  STRING(objname, FORMAT='(A)')   ; Format is required in case the name starts by a number (e.g., 1_Ori)
  IF KEYWORD_SET(LMIR_FW1) THEN SXADDPAR, hdr, 'LMIR_FW1', lmir_fw1
  IF KEYWORD_SET(LMIR_FW2) THEN SXADDPAR, hdr, 'LMIR_FW2', lmir_fw2
  IF KEYWORD_SET(LMIR_FW3) THEN SXADDPAR, hdr, 'LMIR_FW3', lmir_fw3
  IF KEYWORD_SET(LMIR_FW4) THEN SXADDPAR, hdr, 'LMIR_FW4', lmir_fw4
  IF KEYWORD_SET(NOM_FW1)  THEN SXADDPAR, hdr, 'NOMICFW1', nom_fw1
  IF KEYWORD_SET(NOM_FW2)  THEN SXADDPAR, hdr, 'NOMICFW2', nom_fw2
  IF KEYWORD_SET(NOM_APW)  THEN SXADDPAR, hdr, 'NOM_APW' , nom_apw
  IF KEYWORD_SET(PHA_FW1)  THEN SXADDPAR, hdr, 'PHA_FW1' , pha_fw1
  IF KEYWORD_SET(PHA_FW2)  THEN SXADDPAR, hdr, 'PHA_FW2' , pha_fw2
  IF KEYWORD_SET(LBT_LXOS) THEN SXADDPAR, hdr, 'LBT_LXOS', lbt_lxos
  IF KEYWORD_SET(LBT_LYOS) THEN SXADDPAR, hdr, 'LBT_LYOS', lbt_lyos
  IF KEYWORD_SET(LBT_RXOS) THEN SXADDPAR, hdr, 'LBT_RXOS', lbt_rxos
  IF KEYWORD_SET(LBT_RYOS) THEN SXADDPAR, hdr, 'LBT_RYOS', lbt_ryos
  IF KEYWORD_SET(RLOOPON)  THEN SXADDPAR, hdr, 'RLOOPON',  rloopon
  IF KEYWORD_SET(LLOOPON)  THEN SXADDPAR, hdr, 'LLOOPON',  lloopon
  IF KEYWORD_SET(DATATYPE) THEN SXADDPAR, hdr, 'DATATYPE', datatype > 0   ; DATATYPE of 0 are passed as -1 
  IF KEYWORD_SET(OBSTYPE)  THEN SXADDPAR, hdr, 'OBSTYPE',  obstype > 0    ; OBSTYPE of 0 are passed as -1
  IF KEYWORD_SET(FLAG)     THEN SXADDPAR, hdr, 'FLAG',     flag
  IF KEYWORD_SET(PID)      THEN SXADDPAR, hdr, 'PID',      pid > 0 
  
  ; Save fits
  WRITEFITS, data_files[idx_out[i_f]], img_file, hdr 
END         
END

