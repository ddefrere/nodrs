;+
; NAME: NULL_CALIBIMG
; 
; PURPOSE:
;   This function calibrates null data given as input.
;
; INPUTS:
;   date           :  The date to be calibrated
;   cfg_file       :  String with the name of the config file with the reduction parameters
;
; KEYWORDS
;   LOG_FILE       :
;   NO_INSET       : set this keyword to not display the bottom inset with the background null in the TF plot
;   PLOT           : set this keyword to have plot the results
;
; MODIFICATION HISTORY:
;   Version 1.0,  04-APR-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 24-MAY-2024, DD: Update file permission

PRO NULL_CALIBIMG, date, cfg_file, LOG_FILE=log_file, NO_INSET=no_inset, INFO=info, PLOT=plot

drs_version = 1.0
drs_date    = '04-APR-2015'

COMMON GLOBAL, prm, cnf, wav, tgt, pth, drs

; Start actual code
ON_ERROR, 0

; Keyword sanity check
IF NOT KEYWORD_SET(info)       THEN info       = 1
IF NOT KEYWORD_SET(no_inset)   THEN no_inset   = 0

; Define running and plotting paramaters
n_lam     = 20                     ;Number of wavelength bins within the bandwidth
n_time    = 1000                   ;Number of points in the interpolated TF curve
charsize  = 1.2
charthick = 3.0

; DEFINE GLOBAL VARIABLES
; ***********************

; Astronomical and physical constants
GET_PRM, prm

; Read config file with reduction parameters
GET_DRS, drs, 'nodrs/cfg/' + cfg_file
IF drs.null_weight GT 0 THEN wei_exp = 2 ELSE wei_exp = 0                      ; O to turn off weights (2 to turn on)

; Obtain the definition of the configuration
GET_CNF, cnf, INSTRUM=drs.instrum

; Recover the IDL running path
DECLARE_PATH, pth, INSTRUM=drs.instrum

n_crop = 64
n_bin  = 1

; COMPUTE FILE PATH AND READ DATA
; *******************************

; Derive long version of the date
date_lng    = '20' + STRMID(date, 0, 2) + '-' + STRMID(date, 2, 2) + '-' + STRMID(date, 4, 2)
l1fits_path = pth.l1fits_path + date_lng + pth.sep
IF NOT FILE_TEST(l1fits_path) THEN MESSAGE, 'No L1 data for ' + date_lng

; Output directory
l2_dir = pth.l2fits_path + pth.sep + date_lng + pth.sep
IF NOT FILE_TEST(l2_dir) THEN FILE_MKDIR, l2_dir

; Init log file
IF KEYWORD_SET(log_file) THEN BEGIN
  log_file =  l2_dir + date_lng + '_IMG.txt'
  OPENW, lun, log_file, /GET_LUN, WIDTH=800, /APPEND
  PRINTF,lun, ' '
  PRINTF,lun, 'NULL_CALIBIMG.pro version ' + STRING(drs_version, FORMAT='(F3.1)') +  ' -- ' + drs_date + ' -- Denis Defrère - Steward Observatory (ddefrere@email.arizona.edu)'
  PRINTF,lun, 'Proceeding at your own risk now, on '+ SYSTIME()
  PRINTF,lun, ' '
ENDIF ELSE lun = -1

; Merge CAL files
filename = pth.l1fits_path + date_lng + pth.sep + date_lng + '_CAL'
LBTI_MERGEL1FILES, FILE_SEARCH(pth.l1fits_path + date_lng + pth.sep, '*CAL*NULL_IMG.fits', COUNT=n0), FILENAME=filename
LBTI_IMGSEL, filename + '_IMG_ALL.fits', filename + '_DATA_ALL.fits', INFO=info, PLOT=plot;,/MEDIAN

; Now used that image as PSF
psf_file = filename + '_IMG_SEL-MED.fits'

; Merge SCI files
filename = pth.l1fits_path + date_lng + pth.sep + date_lng + '_SCI'
LBTI_MERGEL1FILES, FILE_SEARCH(pth.l1fits_path + date_lng + pth.sep, '*SCI*NULL_IMG.fits', COUNT=n0), FILENAME=filename
LBTI_IMGSEL, filename + '_IMG_ALL.fits', filename + '_DATA_ALL.fits', INFO=info, PLOT=plot;,/MEDIAN
sci_file = pth.l1fits_path + date_lng + pth.sep + date_lng + '_SCI_IMG_SEL.fits'

; Derotate image
LBTI_IMGDEROT, sci_file, 0, PLOT=plot, VERBOSE=info, PSF_FILE=psf_file

END
