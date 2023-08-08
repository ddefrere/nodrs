;+
; NAME: DECLARE_PATH
;
; PURPOSE:
;   This procedure reads the "path.cfg" file and parse the result in the output structure.
;
; MANDATORY INPUT KEYWORD:
;   INSTRUM : Set to the instrument name
; 
; OUTPUT:
;   pth     :  2-dimension array with the image to process
;   
; MODIFICATION HISTORY:
;   Version 1.0,  17-APR-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  08-AUG-2013, DD: output is now a structure
;   Version 1.2,  07-OCT-2013, DD: updated 'result_path' for each instrument
;   Version 1.3,  24-OCT-2013, DD: updated to use of the external hard drive if present
;   Version 1.4,  24-MAY-2014, DD: added dark_path and flat_path to output structure
;   Version 1.5,  11-SEP-2014, DD: now use PATH_SEP to retrieve separator for the current system
;   Version 1.6,  13-JAN-2015, DD: added bpm_path
;   Version 1.7,  16-FEB-2015, DD: now read path from a config file
;   Version 1.8,  15-SEP-2015, DD: added Pickles path
;   Version 1.9,  15-NOV-2015, DD: replaced NOMIC and LMIRCAM keywords by INSTRUM
;   Version 2.0,  29-NOV-2015, DD: updated for new Pickles path
;   Version 2.1,  15-JAN-2016, DD: updated to handle the full paths in path.cfg. 
;   Version 2.1,  20-JUL-2016, DD: added keyword PATH_FILE
;   Version 2.2,  09-APR-2017, DD: added INPUT_PATH to output structure

PRO DECLARE_PATH, pth, INSTRUM=instrum, PATH_FILE=path_file

; Keyword check
IF NOT KEYWORD_SET(INSTRUM)   THEN MESSAGE, 'You must specifiy one detector!"
IF NOT KEYWORD_SET(PATH_FILE) THEN path_file = 'path.cfg'

; Retrieve OS version
sep = PATH_SEP()

; Read config file
pth_cfg = READ_CONFIG('nodrs/cfg/' + path_file)

; Read data path
; Find position of the %instrum string and replace if found.
disk_path = pth_cfg.data_path1
IF STRMATCH(disk_path, '*%instrum*') THEN disk_path = STRJOIN(STRSPLIT(disk_path, '%instrum', /extract, /regex, /preserve_null), instrum)

; Check if directory exist. If not, look at secondary directory.
IF NOT FILE_TEST(disk_path) THEN BEGIN
  disk_path = pth_cfg.data_path2
  IF STRMATCH(disk_path, '*%instrum*') THEN disk_path = STRJOIN(STRSPLIT(disk_path, '%instrum', /extract, /regex, /preserve_null), instrum)
ENDIF

; If the second directory is not defined either, use the code directory
IF NOT FILE_TEST(disk_path) THEN disk_path = GET_PATH('lbti_drs.pro', N_DIR_UP=2)

; Add instrument name (now obsolete)
root_data = disk_path ;+ instrum + sep

; Path to Pickles library
input_path   = 'nodrs/input/'
pickles_path = 'nodrs/input/pickles/'

; Read result path
result_path = pth_cfg.result_path1
IF NOT FILE_TEST(result_path) THEN result_path = pth_cfg.result_path2
result_path = result_path +  instrum + sep              ; Path to result folder
l0fits_path = result_path + 'l0_fits' + sep             ; Path to the processed L0 FITS file
l1fits_path = result_path + 'l1_fits' + sep             ; Path to the fits level 1 files
l2fits_path = result_path + 'l2_fits' + sep             ; Path to the fits level 2 files
bpm_path    = result_path + 'bpm'     + sep             ; Path to bpm folder
dark_path   = result_path + 'dark'    + sep             ; Path to dark folder
flat_path   = result_path + 'flat'    + sep             ; Path to flat folder

; Create directory if does not exist (don't create the directory here anymore)
;IF NOT FILE_TEST(result_path) THEN FILE_MKDIR, result_path
;IF NOT FILE_TEST(bpm_path)    THEN FILE_MKDIR, bpm_path
;IF NOT FILE_TEST(dark_path)   THEN FILE_MKDIR, dark_path
;IF NOT FILE_TEST(flat_path)   THEN FILE_MKDIR, flat_path
;IF NOT FILE_TEST(l0fits_path) THEN FILE_MKDIR, l0fits_path
;IF NOT FILE_TEST(l1fits_path) THEN FILE_MKDIR, l1fits_path
;IF NOT FILE_TEST(l2fits_path) THEN FILE_MKDIR, l2fits_path

; Output structure
pth = {SEP:sep, ROOT_DATA:root_data, DATA_PATH: '', BPM_PATH: bpm_path, DARK_PATH: dark_path, FLAT_PATH: flat_path, INPUT_PATH: input_path, $
      L0FITS_PATH:l0fits_path, L1FITS_PATH:l1fits_path, L2FITS_PATH:l2fits_path, PICKLES_PATH: pickles_path, RESULT_PATH:result_path}
END