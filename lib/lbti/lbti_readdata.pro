;+
; NAME: LBTI_READDATA
; 
; PURPOSE:
;   Read FITS files obtained with LBTI and return a data cube with the images.
;   FITS files are searched in the 'data_path' directory and keyword information is returned in HDR_DATA.
;
; INPUTS:
;   data_path  :  String vector with the path to the fits files with data.
;
; KEYWORDS
;   DATA_IDX   :  Vector with the file numbers of the data files (e.g., [20,119], superseed header information). If two-element vector, the first element
;              :  corresponds to the lower limit number and the second element to the upper limit number (all files in between are considered) .
;   TGT_NAME   :  The name of the star (superseed the name given in the input files)
;   FLAG       :  Set this keyword to the frame type (e.g., SCI, DARK, FLAT,...superseed the value given in the input files)
;   OBSTYPE    :  Set this keyword according to the type of observations:
;                      - 0: regular imaging;
;                      - 1: coherent imaging;
;                      - 2: nulling interferometry;
;                      - 3: background.
;   LAMBDA_CEN :  [m], central wavelength (superseed keyword in the input files)
;   BANDWIDTH  :  [m], bandwidth (superseed keyword in the input files)
;   CROP       :  If set, crop the frames (4-element vector: x_min, y_min, x_max, y_max)
;   LMIR_FW1   :  LMIRCAM filter wheel 1 (superseed keyword in the input files)
;   LMIR_FW2   :  LMIRCAM filter wheel 2 (superseed keyword in the input files)
;   LMIR_FW3   :  LMIRCAM filter wheel 3 (superseed keyword in the input files)
;   LMIR_FW4   :  LMIRCAM filter wheel 4 (superseed keyword in the input files)
;   MEAN       :  If set, return a cube and header data that have been "resistant" averaged (much slower than /MEDIAN).
;              :  If MLOG_DATA is set, average over the config IDs.
;   MEDIAN     :  If set, return a cube and header data that have been median-combined over the config ID.
;              :  If MLOG_DATA is set, median-combined over the config IDs.
;   MLOG_DATA  :  If set, parse the masterlog data to HDR_DATA (see LBTI_READMASTERLOG)
;   NOM_FW1    :  NOMIC filter wheel 1 (superseed keyword in the input files)
;   NOM_FW2    :  NOMIC filter wheel 2 (superseed keyword in the input files)
;   SCALE      :  If set, scale the output MEDIAN or MEAN image to one coadd
;   HDR_DATA   :  On output, a structure containing the information
;   INFO       :  Set this keyword to print info to screen
;   PLOT       :  Set this keyword to plot the results
;
; OUTPUT
;   A data cube with the images contained in the input FITS files. The header info is stored in "DATA_INFO".
;
; MODIFICATION HISTORY:
;   Version 1.0,  17-FEB-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  02-MAY-2013, DD: Added keyword GPU
;   Version 1.2,  06-MAY-2013, DD: Added keyword GOOD_PIX and BAD_PIX
;   Version 1.3,  11-JUN-2013, DD: Added all filter keywords that can superseed the header information
;   Version 1.4,  01-AUG-2013, DD: Modified use of parameter FLAG (unrelated to datatype from now on)
;   Version 1.5,  02-AUG-2013, DD: Speed-up code by a factor ~15 by avoiding stupid concatenations of large data cubes
;   Version 1.6,  08-AUG-2013, DD: Updated according to new formalism of 'declare_path.pro'
;   Version 1.7,  13-AUG-2013, DD: Implemented automatic nod ID attribution
;   Version 1.8,  28-OCT-2013, DD: Improved memory usage and clean up code
;   Version 1.9,  04-NOV-2013, DD: Added keyword 'BAD_IDX'
;   Version 2.0,  13-JAN-2014, DD: Implemented masterlog file
;   Version 2.1,  14-JAN-2014, DD: Added flexibility to 'DATA_IDX'
;   Version 2.2,  15-JAN-2014, DD: Modified loop index number to long to avoid stupid IDL mistake!
;   Version 2.3,  08-FEB-2014, DD: Improved speed, updated for the use with new LBTI_READL0, and removed obsolete filter keywords
;   Version 2.4,  13-MAR-2014, DD: Replaced NO_OVERLAP by OVERLAP
;   Version 2.5,  23-MAY-2014, DD: Added /NO_ZERO to speed-up memory allocation of main data cubes
;   Version 2.6,  22-SEP-2014, DD: Removed XCEN, YCEN, and OVERLAP keywords
;   Version 2.7,  21-OCT-2014, DD: Replaced call to function MATCH by built-in IDL routine "VALUE_LOCATE"
;   Version 2.8,  22-OCT-2014, DD: Minor bug corrected
;   Version 2.9,  03-FEB-2015, DD: Corrected stupid bug with DATA_IDX!
;   Version 2.9,  16-FEB-2015, DD: Removed keywords BAD_IDX, BCKG_IDX, DARK_IDX, FLAT_IDX, BAD_PIX and GOOD_PIX
;   Version 3.0,  03-DEC-2015, DD: Added keyword mapping call to LBTI_READL0.pro
;   Version 3.1,  23-DEC-2015, DD: Now allocate memory depending on array type + added keyword MLOG_DATA + added keyword MEAN
;   Version 3.2,  26-MAR-2016, DD: added keyword SCALE
;   Version 3.3,  26-MAR-2016, DD: added keyword CROP
;   Version 3.4,  03-JUL-2016, DD: more robust access to mlog_data structure + added central value to HDR_DATA
;   Version 3.5,  06-JUL-2016, DD: replace VALUE_LOCATE by MATCH!

FUNCTION LBTI_READDATA, data_path, DATA_IDX=data_idx, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, $
                        LAMBDA_CEN=lambda_cen, BANDWIDTH=bandwidth, CROP=crop, UBC_DXSP=ubc_dxsp, UBC_SXSP=ubc_sxsp, NIC_FSTP=nic_fstp, NIC_BEAM=nic_beam,$
                        LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, MEAN=mean, MEDIAN=median, MLOG_DATA=mlog_data, $
                        NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, NOM_APW=nom_apw, PHA_FW1=pha_fw1, PHA_FW2=pha_fw2, PHA_IMG=pha_img, NIL_DIC=nil_dic, SCALE=scale,$  ; Optional input keywords superseeding keywords in the header
                        HDR_DATA=hdr_data, $                                                                                                                  ; Outputs
                        INFO=info, PLOT=plot   

; Retrieve FITS files in the input directory
data_files = FILE_SEARCH(data_path,'*.fits') & n_files = N_ELEMENTS(data_files)
IF n_files LT 1 THEN MESSAGE, 'Input data path empty'

; Only reads the file in the 'data_idx' range. 
; Warning, this only works assuming the LBTI file definition (i.e., 'n_130424_000123.fits')
data_nbr = LONG(STREGEX(STREGEX(data_files,'_[0-9]+.fits',/EXTRACT),'[0-9]+',/EXTRACT))
IF KEYWORD_SET(DATA_IDX) THEN BEGIN
  IF N_ELEMENTS(DATA_IDX) EQ 2 THEN idx_out = WHERE(data_nbr GE data_idx[0] AND data_nbr LE data_idx[1], n_files) $
                               ELSE MATCH, data_nbr, data_idx, idx_out, idx_tmp 
ENDIF ELSE idx_out = LINDGEN(n_files) ; use LINDGEN because n_files is likely to be too large for INDGEN

; Extract the files
n_files = N_ELEMENTS(idx_out)
IF n_files GT 0 THEN data_files = data_files[idx_out] ELSE MESSAGE, 'File number not found in the given range'

; Don't read the filter wheels if the wavelength and bandwidth are already known
IF KEYWORD_SET(LAMBDA_CEN) AND KEYWORD_SET(BANDWIDTH) THEN no_filter=1

; Read the first file and initiate arrays. LBTI files always come as long integers.
; Also map the keyword position for faster access in the loop below (KEY_MAP is returned here and used in the next call to LBTI_READL0).
img0  = LBTI_READL0(data_files[0], HDR_DATA=hdr0, KEY_MAP=key_map, LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, NO_FILTER=no_filter)     ; Read FITS file
n_x0  = hdr0.n_xpix
n_y0  = hdr0.n_ypix
n_fr0 = hdr0.n_fr  

; Crop the frame if CROP is set
IF KEYWORD_SET(CROP) THEN img0 = TEMPORARY(img0[crop[0]:crop[2],crop[1]:crop[3],*])

; Determine array value type (don't trust the value given in the fits header, which are always LONG)
; This assumes that the images returned by the camera are always in positive integer numbers
; Don't use unsigned types because of background subtraction!
IF MAX(img0) LT 32767 THEN type = 2 ELSE type = 3
img_in   = MAKE_ARRAY((SIZE(img0))[1], (SIZE(img0))[2], n_fr0*n_files, TYPE=type, /NOZERO) ; This assume that all the frames have the same format                      
hdr_data = REPLICATE(hdr0[0], n_fr0*n_files)                                               ; This assume that all the frames have the same format  

; Parse first file in the output arrays and init variables
img_in[0,0,0] = img0
hdr_data[0]   = hdr0
idx_img       = LONG(n_fr0)
n_fr          = 0 

; Loop over the remaining files and append data 
FOR i_f=1L, n_files-1 DO BEGIN
  ; Read image
  img0 = LBTI_READL0(data_files[i_f], HDR_DATA=hdr0, KEY_MAP=key_map, LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, NO_FILTER=no_filter)   ; Read FITS file and save header
  ; Skip frame if bad 
  IF (SIZE(img0))[0] NE 0 THEN BEGIN
    ; Check if original size is the same as that of the first file of this nod
    IF hdr0.n_xpix EQ n_x0 AND hdr0.n_ypix EQ n_y0 THEN BEGIN
      ; Crop the frame if requested
      IF KEYWORD_SET(CROP) THEN img0 = TEMPORARY(img0[crop[0]:crop[2],crop[1]:crop[3],*])
      ; Parse to output arrays (use IDL trick to avoid memory-consuming [*,*] on the left side of the equation)
      img_in[0,0,idx_img] = img0
      hdr_data[idx_img]   = hdr0
      n_fr                = N_ELEMENTS(img0[0,0,*]) 
      ; Increase element index
      idx_img += n_fr
    ENDIF ELSE PRINT, "  SKIPPED: this file doesn't have the same size has the first file of this nod : " + data_files[i_f] 
  ENDIF ELSE PRINT, '  SKIPPED: corrupted file : ' + data_files[i_f]
ENDFOR

; Truncate array if too many frames were initially allocated
IF idx_img NE n_fr0*n_files THEN BEGIN
  idx_keep = LINDGEN(idx_img)
  img_in   = TEMPORARY(img_in[*,*,idx_keep])
  hdr_data = hdr_data[idx_keep]
ENDIF

; Total number of frames
n_tot = N_ELEMENTS(hdr_data.objname)
                                          
; Parse input keywords to header data if set
IF KEYWORD_SET(TGT_NAME)   THEN hdr_data.objname   = tgt_name
IF KEYWORD_SET(OBSTYPE)    THEN hdr_data.obstype   = obstype
IF KEYWORD_SET(LAMBDA_CEN) THEN hdr_data.lam_cen   = lambda_cen
IF KEYWORD_SET(BANDWIDTH)  THEN hdr_data.bandwidth = bandwidth
IF KEYWORD_SET(FLAG)       THEN BEGIN
  hdr_data.flag = flag
  CASE flag OF
    'SCI':  hdr_data.datatype = 0
    'DRK':  hdr_data.datatype = 1
    'FLT':  hdr_data.datatype = 2
    'FLAT': hdr_data.datatype = 2
    'DARK': hdr_data.datatype = 1
    'CAL':  hdr_data.datatype = 3
    'LIN':  hdr_data.datatype = 0
    ELSE: MESSAGE, 'Undefined flag: ' + flag
  ENDCASE
ENDIF

; Parse masterlog data
IF KEYWORD_SET(MLOG_DATA) THEN BEGIN
  MATCH, mlog_data.file_id, hdr_data.file_id, idx_fid, idx_tmp
  hdr_data.nod_id = (mlog_data.nod_id)[idx_fid]
  hdr_data.cfg_id = (mlog_data.cfg_id)[idx_fid]
  hdr_data.pt_id  = (mlog_data.pt_id)[idx_fid]
  hdr_data.cv     = (mlog_data.cv)[idx_fid]
ENDIF

; Compute median
IF KEYWORD_SET(MEDIAN) OR KEYWORD_SET(MEAN) THEN img_in = LBTI_IMGMED(TEMPORARY(img_in), HDR_DATA=hdr_data, MEAN=mean, SCALE=scale)

RETURN, img_in
END

; Test routine to read only a sub frame of a fits file
; Results show that REAFITS is faster (even if it reads the whole image)
PRO TEST_SUBFITS
  ; Define subframe and read once to put array in shared memory
  lim = [0,0,127,255]
  img = READFITS('/Volumes/LaCie/data/LBTI/nomic/160326/n_160326_002255.fits', /SILENT)
  ; First approach. Read everything with REAFITS and then crop
  t0  = SYSTIME(1)
  img = READFITS('/Volumes/LaCie/data/LBTI/nomic/160326/n_160326_002255.fits', /SILENT)
  img = TEMPORARY(img[lim[0]:lim[2],lim[1]:lim[3]])
  t1  = SYSTIME(1)
  ; Approach 2, read only a subframe with MRDFITS
  ;row  = lim[0] + INDGEN(lim[2]-lim[0])
  ;col  = lim[1] + INDGEN(lim[3]-lim[1])
  ;img2 = MRDFITS('/Volumes/LaCie/data/LBTI/nomic/160326/n_160326_002255.fits', COLUMNS=col, ROWS=row, /SILENT)
  img2 = MRDFITS('/Volumes/LaCie/data/LBTI/nomic/160326/n_160326_002255.fits', /SILENT)
  t2   = SYSTIME(1)
  PLOTXY, img
  PLOTXY, img2
  PRINT, 'execution time :', t1-t0, t2-t1
END