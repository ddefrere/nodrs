;+
; NAME: LBTI_READL0
; 
; PURPOSE:
;   Read a single FITS file produced by LBTI and return the image data cube.
;   Keyword information is returned by the keyword 'HDR_DATA'.
;
; INPUTS:
;   lbti_file  :  String vector with the path to the fits file.
;
; KEYWORDS
;   HDR_DATA    :  On output, a structure containing the information
;   KEY_MAP     :  
;   NO_FILTER   :  If set, the filter wheels are not read (WARNING: the wavelength and the bandwidth are hence not computed).
;   INFO        :  Set this keyword to print info to screen
;
; OUTPUT
;   A data cube with the image contained in the input FITS file.
;
; MODIFICATION HISTORY:
;   Version 1.0,  17-APR-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1,  27-JUN-2013, DD: now account for old nomenclature of filter wheels
;   Version 1.2,  28-JUN-2013, DD: added 'lam_cen' and 'bandwidth' to output structure
;   Version 1.3,  08-JUL-2013, DD: added nod and chop positions to output structure
;   Version 1.4,  28-JUL-2013, DD: now header data output as an array of structure rather than a structre of arrays
;   Version 1.4,  28-JUL-2013, DD: added parameter flag
;   Version 1.5,  28-JUL-2013, DD: added file ID number to output structure
;   Version 1.6,  13-AUG-2013, DD: added LBT offsets to output structure
;   Version 1.7,  19-SEP-2013, DD: added the OB number to output structure
;   Version 1.8,  10-OCT-2013, DD: added AO and PHASEcam parameters to the output structure
;   Version 1.9,  17-OCT-2013, DD: added LMIRCam aperture wheel keyword
;   Version 2.0,  02-NOV-2013, DD: adapted to read L1 files
;   Version 2.1,  12-DEC-2013, DD: added bias correction
;   Version 2.2,  08-FEB-2014, DD: improved speed
;   Version 2.2,  09-MAR-2014, DD: added new PHASEcam outputs
;   Version 2.3,  13-MAR-2014, DD: Replaced NO_OVERLAP by OVERLAP
;   Version 2.4,  31-MAR-2014, DD: Added MJD computation
;   Version 2.5,  27-MAY-2014, DD: Added keyword NO_FILTER
;   Version 2.6,  09-SEP-2014, DD: Added option /NO_CONTINUE to header reading to speed-up routine
;   Version 2.7,  22-SEP-2014, DD: Added bck_nod to output structure
;   Version 2.8,  01-OCT-2014, DD: Updated AO keywords with new names
;   Version 2.9,  28-OCT-2014, DD: Now subtract reference pixels only if non-CDS data + added PID to output structure
;   Version 3.0,  07-NOV-2014, DD: Added slope to output structure
;   Version 3.1,  11-NOV-2014, DD: Added new PHASECam fields to output structure
;   Version 3.2,  30-NOV-2014, DD: Force PID to be integer
;   Version 3.3,  08-JAN-2015, DD: SAOMODE and DAOMODE no longer in the L0 fits header
;   Version 3.4,  08-FEB-2015, DD: Minor bug corrected in PHASECam header
;   Version 3.5,  16-FEB-2015, DD: Added new PHASECam fieds (required by NSC)
;   Version 3.6,  18-FEB-2015, DD: Moved bias subtraction to external routine
;   Version 3.7,  20-FEB-2015, DD: Removed 'slope' from the output structure
;   Version 3.8,  04-APR-2015, DD: Added nodding period and number of frames in background
;   Version 3.9,  10-OCT-2015, DD: Added PACTCDLY to output structure
;   Version 4.0,  19-OCT-2015, DD: Added new names for UBC image stops (RMGSTP and LMGSTP)
;   Version 4.1,  23-NOV-2015, DD: Added pointing ID to output structure
;   Version 4.2,  02-DEC-2015, DD: Added keyword KEY_MAP and added call to DEFINE_HDRDATA
;   Version 4.3,  22-DEC-2015, DD: Corrected minor bug
;   Version 4.4,  16-FEB-2016, DD: Added data consistancy checks
;   Version 4.5,  13-MAR-2016, DD: Corrected bug with backward compatibility
;   Version 4.6,  06-JUL-2016, DD: Now reset keyword map if it has changed
;   Version 4.7,  05-OCT-2016, DD: Replaced STRTRIM by STRCOMPRESS to remove middle white spaces too
;   Version 4.8,  19-OCT-2016, DD: Now use AO status keywords to check whether the AO loop is open or closed
;   Version 4.9,  20-OCT-2016, DD: Now read keyword INTERFAC
;   Version 5.0,  26-OCT-2016, DD: Added keyword RSLOPRMS and LSLOPRMS
;   Version 5.1,  30-JAN-2017, DD: Added keyword SPDTHPOS (OPD dither position)
;   Version 5.2,  12-APR-2017, DD: added region of interest!
;   Version 5.3,  12-MAY-2017, DD: updated PCPLSP to PCPLSP1
;   Version 5.4,  03-MAR-2019, DD: Now handle LMIRCAM_DETECTOR as instrument name

FUNCTION LBTI_READL0, lbti_file, HDR_DATA=hdata, KEY_MAP=key_map, LMIR_FW1=lmir_fw1, LMIR_FW2=lmir_fw2, LMIR_FW3=lmir_fw3, LMIR_FW4=lmir_fw4, NOM_FW1=nom_fw1, NOM_FW2=nom_fw2, NO_DETECTOR=no_detector, NO_FILTER=no_filter, INFO=info

; 0. Read detector image 
IF NOT FILE_TEST(lbti_file) THEN MESSAGE, 'Input file missing :' + lbti_file
img_in   = READFITS(lbti_file, header, /SILENT)  ; READ_FITS faster than MRDFITS and FITS_READ
size_img = SIZE(img_in)
n_key    = N_ELEMENTS(header)

; 1. Keyword extraction paramater
; If the keyword map is not given, compute it in this pass
pre_check  = 0. ; look keyword one position before the given start position
post_check = 0. ; look keyword one position after the given start position
IF KEYWORD_SET(KEY_MAP) THEN BEGIN
  ; Check than the header length has not changed (saved as p99 in previous pass)
  IF n_key EQ key_map.p99 THEN BEGIN
    IF KEYWORD_SET(KEY_MAP) THEN DEFINE_KEYMAP, key_map, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29,$
                                                p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58,$
                                                p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87,$
                                                p88, p89, p90, p91, p92, p93, p94, p95, p96, p97, p98, p99 
  ENDIF ELSE key_map = 0 ; reset keyword map since it has changed  
ENDIF



; 2. Define header structure
DEFINE_HDRDATA, hdata

; 3. Check data integrity
hdata.instrum  = STRUPCASE(STRTRIM(FXPAR(header, 'INSTRUME', START=p1, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2))
IF hdata.instrum EQ 'LMIRCAM_DETECTOR' THEN  hdata.instrum = 'LMIRCAM'
hdata.acq_time = FXPAR(header, 'ACQTIME', START=p2, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) 
IF hdata.instrum NE 'NOMIC' AND hdata.instrum NE 'LMIRCAM' AND hdata.instrum NE 'MIRAC' OR NOT FINITE(hdata.acq_time) THEN RETURN, 0

; 4. Derive file ID number
hdata.file_id  = STREGEX(STREGEX(lbti_file,'_[0-9]+.fits',/EXTRACT),'[0-9]+',/EXTRACT)

; 5. Read FITS header. 
; The double or float is to avoid conflicting structures (sometimes these parameters are float, and sometimes double!!!!)
; Image and timing parameters
hdata.filename = FXPAR(header, 'FILENAME', START=p3, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) 
hdata.date_obs = FXPAR(header, 'DATE-OBS', START=p4, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE);& IF !err EQ -1 THEN date_obs = STRTRIM(FXPAR(header, 'DATE_OBS'), 2)
hdata.time_end = FXPAR(header, 'TIME-END', START=p5, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.pid      = FXPAR(header, 'PID', START=p6, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) & IF !err EQ -1 THEN hdata.pid = -1 ; -1 indicates that no PID has been found           

; Target parameters
hdata.objname = STRCOMPRESS(STRING(FXPAR(header, 'OBJNAME', START=p7, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)), /REMOVE_ALL)
hdata.flag    = STRTRIM(STRING(FXPAR(header, 'FLAG', START=p8, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)), 2)
IF !err NE -1 THEN BEGIN
  CASE hdata.flag OF
    'SCI': hdata.datatype = 0
    'DRK': hdata.datatype = 1
    'FLT': hdata.datatype = 2
    'CAL': hdata.datatype = 3
    ELSE : hdata.datatype = 0
  ENDCASE
ENDIF ELSE BEGIN  ; Means old file format
  hdata.datatype = FIX(FXPAR(header, 'DATATYPE', /NOCONTINUE))
  CASE hdata.datatype OF
    0: hdata.flag='SCI'
    1: hdata.flag='DRK'
    2: hdata.flag='FLT'
    3: hdata.flag='CAL'
    ELSE : hdata.flag='UND'
  ENDCASE
ENDELSE
obstype = FXPAR(header, 'OBSTYPE', START=p9, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) & IF STRTRIM(STRING(obstype)) EQ 'object' THEN obstype = 0 ELSE obstype = FIX(obstype) 
hdata.obstype  = obstype
hdata.lbt_utc  = FXPAR(header, 'LBT_UTC', START=p10, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.lbt_lst  = FXPAR(header, 'LBT_LST', START=p11, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.lbt_ra   = FXPAR(header, 'LBT_RA', START=p12, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.lbt_dec  = FXPAR(header, 'LBT_DEC', START=p13, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.lbt_alt  = FXPAR(header, 'LBT_ALT', START=p14, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.lbt_az   = FXPAR(header, 'LBT_AZ', START=p15, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.lbt_para = FXPAR(header, 'LBT_PARA', START=p16, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)

; Detector and integration
hdata.smplmode = FXPAR(header, 'SMPLMODE', START=p17, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.detmode  = FXPAR(header, 'DETMODE', START=p18, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.int_time = FXPAR(header, 'EXPTIME', START=p19, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.n_coadd  = FXPAR(header, 'NCOADDS', START=p20, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.pagain   = FXPAR(header, 'PAGAIN', START=p21, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.pabandw  = FXPAR(header, 'PABANDW', START=p22, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.pactcdly = FXPAR(header, 'PACTCDLY', START=p23, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.detbias  = FXPAR(header, 'DETBIAS', START=p24, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.eperadu  = FXPAR(header, 'EPERADU', START=p25, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.interfac = FXPAR(header, 'INTERFAC', START=p26, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)

; Compute MJD
hr  = DOUBLE(STRMID(hdata.time_end, 0, 2)) + DOUBLE(STRMID(hdata.time_end, 3, 2))/60D + DOUBLE(STRMID(hdata.time_end, 6, 2))/3.6D+3 + DOUBLE(STRMID(hdata.time_end, 9, 3))/3.6D+6
JDCNV, FIX(STRMID(hdata.date_obs, 0, 4)), FIX(STRMID(hdata.date_obs, 5, 2)), FIX(STRMID(hdata.date_obs, 8, 2)), hr, jd
hdata.mjd_obs = jd - 2400000.5D  ; mjd

; Reduction parameters (most will be populated later)
hdata.chp_id  = -99             ; Chop ID number
hdata.nod_id  = -99             ; Nod ID number, -99 for unused frames (-1 is for background frames and 0,1,2,.. forn regular nod frames)
hdata.n_xpix  = size_img[1]     ; Initial number of X pixels
hdata.n_ypix  = size_img[2]     ; Initial number of Y pixels
IF size_img[0] EQ 3 THEN hdata.n_fr = size_img[3] ELSE hdata.n_fr = 1

; Read filter wheels if N_FILTER is not set (the IFs account for old nomenclature)
IF NOT KEYWORD_SET(NO_FILTER) THEN BEGIN
  IF NOT KEYWORD_SET(LMIR_FW1) THEN hdata.lmir_fw1 = STRTRIM(FXPAR(header, 'LMIR_FW1', START=p27, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2) & IF !err EQ -1 THEN hdata.lmir_fw1 = STRTRIM(FXPAR(header, 'LMIR_FW1', START=p27, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2)  ; !err check for backward compatibility
  IF NOT KEYWORD_SET(LMIR_FW2) THEN hdata.lmir_fw2 = STRTRIM(FXPAR(header, 'LMIR_FW2', START=p28, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2) & IF !err EQ -1 THEN hdata.lmir_fw2 = STRTRIM(FXPAR(header, 'LMIR_FW2', START=p28, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2)
  IF NOT KEYWORD_SET(LMIR_FW3) THEN hdata.lmir_fw3 = STRTRIM(FXPAR(header, 'LMIR_FW3', START=p29, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2) & IF !err EQ -1 THEN hdata.lmir_fw3 = STRTRIM(FXPAR(header, 'LMIR_FW3', START=p29, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2)
  IF NOT KEYWORD_SET(LMIR_FW4) THEN hdata.lmir_fw4 = STRTRIM(FXPAR(header, 'LMIR_FW4', START=p30, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2) & IF !err EQ -1 THEN hdata.lmir_fw4 = STRTRIM(FXPAR(header, 'LMIR_FW4', START=p30, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2)
  IF NOT KEYWORD_SET(NOM_FW1)  THEN hdata.nom_fw1  = STRTRIM(FXPAR(header, 'NOMICFW1', START=p31, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2) & IF !err EQ -1 THEN hdata.nom_fw1  = STRTRIM(FXPAR(header, 'NOM_FW1', START=p31, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2)
  IF NOT KEYWORD_SET(NOM_FW2)  THEN hdata.nom_fw2  = STRTRIM(FXPAR(header, 'NOMICFW2', START=p32, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2) & IF !err EQ -1 THEN hdata.nom_fw2  = STRTRIM(FXPAR(header, 'NOM_FW2', START=p32, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE), 2)
  hdata.ubc_dxsp = STRTRIM(STRING(FXPAR(header, 'RMGSTP', START=p33, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)), 2) & IF !err EQ -1 THEN hdata.nom_fw2  = STRTRIM(STRING(FXPAR(header, 'UBC_DXSP', START=p33, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)), 2)
  hdata.ubc_sxsp = STRTRIM(STRING(FXPAR(header, 'LMGSTP', START=p34, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)), 2) & IF !err EQ -1 THEN hdata.nom_fw2  = STRTRIM(STRING(FXPAR(header, 'UBC_SXSP', START=p34, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)), 2)
  hdata.nic_fstp = STRING(FXPAR(header, 'NIC_FSTP', START=p35, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE))
  hdata.nic_beam = STRING(FXPAR(header, 'NIC_BEAM', START=p36, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE))
  hdata.lmir_awl = STRING(FXPAR(header, 'LMIR_AWL', START=p37, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE))
  hdata.nom_apw  = STRING(FXPAR(header, 'NOM_APW',  START=p38, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE))
  hdata.pha_fw1  = STRING(FXPAR(header, 'PHA_FW1',  START=p39, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE))
  hdata.pha_fw2  = STRING(FXPAR(header, 'PHA_FW2',  START=p40, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE))
  hdata.pha_img  = STRING(FXPAR(header, 'PHA_IMG',  START=p41, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE))
  hdata.nil_dic  = ' ';STRING(STRTRIM(FXPAR(header, 'NIL_DIC'), 2))  & IF !err EQ -1 THEN nil_dic = 'N/A'  ; not used for now
  LBTI_WAVEBAND, lam_cen, bandwidth, INSTRUM=hdata.instrum, LMIR_FW1=hdata.lmir_fw1, LMIR_FW2=hdata.lmir_fw2, LMIR_FW3=hdata.lmir_fw3, LMIR_FW4=hdata.lmir_fw4, NOM_FW1=hdata.nom_fw1, NOM_FW2=hdata.nom_fw2   ; derive wavelength and bandwidth based on the header
  hdata.lam_cen   = lam_cen
  hdata.bandwidt  = bandwidth
ENDIF

; AO parameters
; As of Sep. 2016, loop ON keywords are not always reliable! They can be 0 while the loop is closed (ex. nod 19, 2016-10-16).
; To avoid rejecting good frames, I check the AOstatus keywords from now on (must be 'AORunning' if the loop is closed).
; Frame selection is done later based on the SLOPE RMS and the fitted FWHM (when applicable)
hdata.daomode  = 'ACE-AO';FXPAR(header, 'Rmode') ; no longer in the header 
hdata.daostrhl = -1 ; FXPAR(header, 'RStrehl', START=p, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) ; no longer in header
hdata.dcmodes  = FXPAR(header, 'RCMODES', START=p42, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.dloopgn  = FXPAR(header, 'RGAIN', START=p43, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
IF FIX(FXPAR(header, 'RLOOPON', START=p44, /NOCONTINUE)) NE 0 OR STRCOMPRESS(FXPAR(header, 'RSTATUS', START=p45, /NOCONTINUE), /REMOVE_ALL) EQ 'AORunning' THEN hdata.dloopon = 1 ELSE hdata.dloopon = 0 & IF !err EQ -1 THEN hdata.dloopon = 1  ; Assume loop is closed if not present
hdata.dwfscfrq = FXPAR(header, 'RWccdFrq', START=p46, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.dwfscbin = FXPAR(header, 'RWccdBin', START=p47, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.dsloprms = FXPAR(header, 'RSLOPRMS', START=p48, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.dslopmax = FXPAR(header, 'RSLOPMAX', START=p49, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.dsafeskp = FXPAR(header, 'RSAFESKP', START=p50, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.saomode  = 'ACE-AO' ;FXPAR(header, 'Lmode') ; no longer in the header 
hdata.saostrhl = -1 ; FXPAR(header, 'LStrehl', START=p, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) ;  no longer in header
hdata.scmodes  = FXPAR(header, 'LCMODES', START=p51, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.sloopgn  = FXPAR(header, 'LGAIN', START=p52, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) 
IF FIX(FXPAR(header, 'LLOOPON', START=p53, /NOCONTINUE)) NE 0 OR STRCOMPRESS(FXPAR(header, 'LSTATUS', START=p54, /NOCONTINUE), /REMOVE_ALL) EQ 'AORunning' THEN hdata.sloopon = 1 ELSE hdata.sloopon = 0 & IF !err EQ -1 THEN hdata.sloopon = 1  ; Assume loop is closed if not present
hdata.swfscfrq = FXPAR(header, 'LWccdFrq', START=p55, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.swfscbin = FXPAR(header, 'LWccdBin', START=p56, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.ssloprms = FXPAR(header, 'LSLOPRMS', START=p57, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.sslopmax = FXPAR(header, 'LSLOPMAX', START=p58, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.ssafeskp = FXPAR(header, 'LSAFESKP', START=p59, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)

; Phasecam parameters
; Lots of backward compatibility here
hdata.pcclosed   = FXPAR(header, 'PCCLOSED', START=p60, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) & IF !err EQ -1 THEN hdata.pcclosed = 1  ; Assume loop is closed if not present
hdata.spc_pist   = FXPAR(header, 'SPCPIST',  START=p61, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.spc_az     = 0. & p62 = 0 ;FXPAR(header, 'SPCAZ',    START=p62, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) ; not yet included
hdata.spc_el     = 0. & p63 = 0 ;FXPAR(header, 'SPCEL',    START=p63, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) ; not yet included
hdata.fpc_pistm  = FXPAR(header, 'FPCPISTM', START=p64, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.fpc_pists  = FXPAR(header, 'FPCPISTS', START=p65, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.fpc_azm    = FXPAR(header, 'FPCAZM',   START=p66, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.fpc_azs    = FXPAR(header, 'FPCAZS',   START=p67, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.fpc_elm    = FXPAR(header, 'FPCELM',   START=p68, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.fpc_els    = FXPAR(header, 'FPCELS',   START=p69, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.pcplsp     = FXPAR(header, 'PCPLSP1',   START=p70, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) & IF !err EQ -1 THEN hdata.pcplsp     = FXPAR(header, 'PCPLSP',   START=p70, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.pctipsp    = FXPAR(header, 'PCTIPSP',  START=p71, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.pctltsp    = FXPAR(header, 'PCTLTSP',  START=p72, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.pcfjmps    = FXPAR(header, 'PCFJMPS',  START=p73, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.spdthpos   = FXPAR(header, 'SPDTHPOS', START=p74, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.pcphmean   = FXPAR(header, 'PCPHMEN1', START=p75, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)  
IF !err NE -1 THEN BEGIN 
  hdata.pcphstd   = FXPAR(header, 'PCPHSTD1', START=p76, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
  hdata.pcphmcos  = FXPAR(header, 'PCPHMCS1', START=p77, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
  hdata.pcphmsin  = FXPAR(header, 'PCPHMSN1', START=p78, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
  hdata.pcmsnr    = FXPAR(header, 'PCMSNR1',  START=p79, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) 
  hdata.pcphmean2 = FXPAR(header, 'PCPHMEN2', START=p80, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
  hdata.pcphstd2  = FXPAR(header, 'PCPHSTD2', START=p81, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
  hdata.pcphmcos2 = FXPAR(header, 'PCPHMCS2', START=p82, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
  hdata.pcphmsin2 = FXPAR(header, 'PCPHMSN2', START=p83, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
  hdata.pcmsnr2   = FXPAR(header, 'PCMSNR2',  START=p84, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
ENDIF ELSE BEGIN 
  hdata.pcphmean  = FXPAR(header, 'PCPHMEAN', START=p75, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
  hdata.pcphstd   = FXPAR(header, 'PCPHSTD',  START=p76, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
  hdata.pcphmcos  = FXPAR(header, 'PCPHMCOS', START=p77, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
  hdata.pcphmsin  = FXPAR(header, 'PCPHMSIN', START=p78, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
  hdata.pcmsnr    = FXPAR(header, 'PCMSNR',   START=p79, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) & IF !err EQ -1 THEN hdata.pcmsnr = FXPAR(header, 'PCMSWR', START=p79, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE) ; Old name
  hdata.pcphmean2 = 0. & p80 = 0
  hdata.pcphstd2  = 0. & p81 = 0
  hdata.pcphmcos2 = 0. & p82 = 0
  hdata.pcphmsin2 = 0. & p83 = 0
  hdata.pcmsnr2   = 0. & p84 = 0
ENDELSE

; Weather conditions
hdata.seeing   = FXPAR(header, 'SEEING',  START=p85, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.smttau   = FXPAR(header, 'SMTTAU',  START=p86, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.lbttemp  = FXPAR(header, 'LBTTEMP', START=p87, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.winddir  = FXPAR(header, 'WINDDIR', START=p88, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
hdata.windspd  = FXPAR(header, 'WINDSPD', START=p89, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)

; Extract region of interest
subsectm  = FXPAR(header, 'SUBSECTM',  START=p90, PRECHECK=pre_check, POSTCHECK=post_check, /NOCONTINUE)
IF !err NE -1 THEN hdata.roi = LBTI_ROI(subsectm, INSTRUM=hdata.instrum) 

; 5. Save key_map if undefimne
p91 = 0 & p92 = 0 & p93 = 0 & p94 = 0 & p95 = 0 & p96 = 0 & p97 = 0 & p98 = 0 & p99 = 0  ; unused keys
IF NOT KEYWORD_SET(KEY_MAP) THEN DEFINE_KEYMAP, key_map, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, $
                                                p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, $
                                                p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, $
                                                p88, p89, p90, p91, p92, p93, p94, p95, p96, p97, p98, n_key, /SAVE

; 6. Additional manipulations if it's a CDS L0 file
; If CDS, only keep the difference frame and adjust integration time
; OR keep both frames if the integration time is equal to the frame time
; If not CDS, subtract the reference pixels with routine LBTI_BIASSUBTRACT
IF hdata.smplmode EQ 6 OR hdata.n_fr GT 1 THEN BEGIN          ; Second condition is a temporary hack. The new data does not seem to have a Sample mode. How do we recognize CDS?
   frmtime = DOUBLE(FXPAR(header, 'FRMTIME', /NOCONTINUE))    ; Backward compatibility (for data before Summer 2018, before electronics upgrade)
   IF !err NE -1 THEN BEGIN  
     dit_eff = hdata.int_time - frmtime
     IF dit_eff NE 0 THEN BEGIN
        img_in     = img_in[*,*,2]-img_in[*,*,1]
        hdata.n_fr = 1
        hdata[*].int_time = dit_eff
     ENDIF ELSE BEGIN
        img_in     = TEMPORARY(img_in[*,*,[1,2]])
        hdata.n_fr = 2
     ENDELSE
   ENDIF ELSE BEGIN
      frmtime = DOUBLE(FXPAR(header, 'FRAME', /NOCONTINUE)) 
      dit_eff = hdata.int_time - frmtime
      IF dit_eff NE 0 THEN BEGIN
        img_in     = img_in[*,*,1]-img_in[*,*,0]
        hdata.n_fr = 1
        hdata[*].int_time = dit_eff
      ENDIF ELSE BEGIN
        img_in     = TEMPORARY(img_in[*,*,[0,1]])
        hdata.n_fr = 2
      ENDELSE 
   ENDELSE
ENDIF

; 7. Replicate header
IF hdata.n_fr GT 1 THEN hdata = REPLICATE(hdata, hdata.n_fr)

; 8. Attribute frame ID
hdata[*].fr_id = INDGEN(hdata.n_fr)

RETURN, img_in
END
