;+
; NAME: LBTI_READL0RED
; 
; PURPOSE:
;   Read a single pre-processed L0 FITS file and return the image data cube.
;   Keyword information is returned by the keyword 'HDR_DATA'.
;
; INPUTS:
;   lbti_file  :  String vector with the path to the fits file.
;
; KEYWORDS
;   HDR_DATA   :  On output, a structure containing the information
;   INFO       :  Set this keyword to print info to screen
;
; OUTPUT
;   A data cube with the image contained in the input FITS file.
;
; MODIFICATION HISTORY:
;   Version 1.0,  08-FEB-2014, Denis Defrère, based on a more general routine LBTI_READFILE
;   Version 1.1,  09-MAR-2014, DD: added new PHASEcam outputs
;   Version 1.2,  24-MAY-2014, DD: removed fixed reduction paramaters from the output structure
;   Version 1.3,  25-MAY-2014, DD: improved speed
;   Version 1.4,  25-OCT-2014, DD: added PID to output structure
;   Version 1.5,  15-JAN-2015, DD: added NO_FLAT, NO_DARK, and NO_BPM to output structure
;   Version 1.6,  20-DEB-2015, DD: added slopes to output structure
;   Version 1.7,  04-APR-2015, DD: Added nodding period and frames per background
;   Version 1.8,  04-MAY-2015, DD: Added background OB id
;   Version 1.9,  24-MAY-2015, DD: Added PHASECam information in the header
;   Version 2.0,  06-OCT-2015, DD: Added PCWAV and PCSPEC to output header
;   Version 2.1,  06-NOV-2015, DD: Added config ID, pointing ID, and BCKG_SEL
;   Version 2.2,  23-DEC-2015, DD: Added precision and image mode to output header
;   Version 2.3,  16-OCT-2016, DD: Added FWHM fields to output structure

FUNCTION LBTI_READL0RED, lbti_file, HDR_DATA=hdr_data, INFO=info

; 1. Read detector image
IF NOT FILE_TEST(lbti_file) THEN MESSAGE, 'Input file missing :' + lbti_file
img_in   = READFITS(lbti_file, header, /SILENT) 
size_img = SIZE(img_in)
n_xpix   = size_img[1]     ; Initial number of X pixels
n_ypix   = size_img[2]     ; Initial number of Y pixels
IF size_img[0] EQ 3 THEN n_fr = size_img[3] ELSE n_fr = 1L     

; 2. Check data integrity
;IF !err EQ -1 THEN RETURN, 0

; 3. Read FITS header
; Read object information
;telescop = STRTRIM(SXPAR(header, 'TELESCOP'), 2)  
instrume = STRTRIM(SXPAR(header, 'INSTRUME'), 2)  
date_obs = STRTRIM(SXPAR(header, 'DATE_OBS'), 2)  
objname  = STRCOMPRESS(STRTRIM(STRING(SXPAR(header, 'OBJNAME')), 2), /REMOVE_ALL)
pid      = 4; temporary hack!!FIX(SXPAR(header, 'PID'))
datatype = FIX(SXPAR(header, 'DATATYPE'))
obstype  = FIX(SXPAR(header, 'OBSTYPE'))
flag     = STRTRIM(STRUPCASE((SXPAR(header, 'FLAG'))))
lbt_ra   = SXPAR(header, 'LBT_RA')
lbt_dec  = SXPAR(header, 'LBT_DEC')
cfg_id   = FIX(SXPAR(header, 'CFG_ID'))
nod_id   = FIX(SXPAR(header, 'NOD_ID'))
pt_id    = FIX(SXPAR(header, 'PT_ID'))
nod_frq  = FIX(SXPAR(header, 'NOD_FRQ'))
n_frbck  = FIX(SXPAR(header, 'N_FRBCK'))

; Detector and integration parameters
pixscale = FLOAT(SXPAR(header, 'PIXSCALE'))
detmode  = FLOAT(SXPAR(header, 'DETMODE'))
smplmode = FLOAT(SXPAR(header, 'SMPLMODE'))
int_time = DOUBLE(SXPAR(header, 'EXPTIME'))
acq_time = DOUBLE(SXPAR(header, 'ACQTIME'))
n_coadd  = FIX(SXPAR(header, 'NCOADDS'))
frmtime  = DOUBLE(SXPAR(header, 'FRMTIME'))
pagain   = LONG(SXPAR(header, 'PAGAIN'))
pabandw  = LONG(SXPAR(header, 'PABANDW'))
detbias  = FLOAT(SXPAR(header, 'DETBIAS'))
eperadu  = FLOAT(SXPAR(header, 'EPERADU'))  
     
; AO information
daomode  = SXPAR(header, 'DAOMODE', /NOCONTINUE)
dcmodes  = SXPAR(header, 'DCMODES', /NOCONTINUE)
dloopgn  = SXPAR(header, 'DLOOPGN', /NOCONTINUE)
dwfscfrq = SXPAR(header, 'DWFSCFRQ', /NOCONTINUE)
dwfscbin = SXPAR(header, 'DWFSCBIN', /NOCONTINUE)
saomode  = SXPAR(header, 'SAOMODE', /NOCONTINUE)
scmodes  = SXPAR(header, 'SCMODES', /NOCONTINUE)
sloopgn  = SXPAR(header, 'SLOOPGN', /NOCONTINUE)
swfscfrq = SXPAR(header, 'SWFSCFRQ', /NOCONTINUE)
swfscbin = SXPAR(header, 'SWFSCBIN', /NOCONTINUE)
    
; Filter information
lam_cen   = FLOAT(SXPAR(header, 'WAVELENG'))
bandwidth = FLOAT(SXPAR(header, 'BANDWIDT'))
ubc_dxsp  = STRTRIM(STRING(SXPAR(header, 'UBC_DXSP')), 2)
ubc_sxsp  = STRTRIM(STRING(SXPAR(header, 'UBC_SXSP')), 2)
nic_fstp  = STRTRIM(STRING(SXPAR(header, 'NIC_FSTP')), 2)
nic_beam  = STRTRIM(STRING(SXPAR(header, 'NIC_BEAM')), 2)
lmir_fw1  = STRTRIM(STRING(SXPAR(header, 'LMIR_FW1')), 2)
lmir_fw2  = STRTRIM(STRING(SXPAR(header, 'LMIR_FW2')), 2)
lmir_fw3  = STRTRIM(STRING(SXPAR(header, 'LMIR_FW3')), 2)
lmir_fw4  = STRTRIM(STRING(SXPAR(header, 'LMIR_FW4')), 2)
nom_fw1   = STRTRIM(STRING(SXPAR(header, 'NOM_FW1')), 2)
nom_fw2   = STRTRIM(STRING(SXPAR(header, 'NOM_FW2')), 2)
nom_apw   = STRTRIM(STRING(SXPAR(header, 'NOM_APW')), 2)
lmir_awl  = STRTRIM(STRING(SXPAR(header, 'LMIR_AWL')), 2)
nom_apw   = STRTRIM(STRING(SXPAR(header, 'NOM_APW')), 2)
pha_fw1   = STRTRIM(STRING(SXPAR(header, 'PHA_FW1')), 2)
pha_fw2   = STRTRIM(STRING(SXPAR(header, 'PHA_FW2')), 2)
pha_img   = STRTRIM(STRING(SXPAR(header, 'PHA_IMG')), 2)

; Weather information
seeing    = FLOAT(SXPAR(header, 'SEEING'))
smttau    = FLOAT(SXPAR(header, 'SMTTAU'))
lbttemp   = FLOAT(SXPAR(header, 'LBTTEMP'))
winddir   = FLOAT(SXPAR(header, 'WINDDIR'))
windspd   = FLOAT(SXPAR(header, 'WINDSPD'))

; Read image reduction parameters
drs_vers  = SXPAR(header, 'DRS_VERS')
drs_date  = SXPAR(header, 'DRS_DATE')
date_red  = SXPAR(header, 'DATE_RED')
bck_mode  = FIX(SXPAR(header, 'BCK_MODE'))
bck_sel   = FIX(SXPAR(header, 'BCK_SEL'))
bck_nod   = FIX(SXPAR(header, 'BCK_NOD'))
img_mode  = FIX(SXPAR(header, 'IMG_MODE'))
no_bpm    = FIX(SXPAR(header, 'NO_BPM'))
no_dark   = FIX(SXPAR(header, 'NO_DARK'))
no_flat   = FIX(SXPAR(header, 'NO_FLAT'))
precision = FIX(SXPAR(header, 'PRECISIO'))

; Will be populated later
bck_ob    = 0
beam_id   = 0 
n_rej     = 0

; Parse information into structure
header = {INSTRUM: instrume, DATE_OBS:date_obs, OBJNAME:objname, PID: pid, DATATYPE:datatype, OBSTYPE:obstype, FLAG:flag, LBT_RA:lbt_ra, LBT_DEC:lbt_dec, CFG_ID: cfg_id, NOD_ID: nod_id, $
          NOD_FRQ: nod_frq, N_FRBCK: n_frbck, PT_ID: pt_id, $                                                                                                                                         ; Target and instrument info
          PIXSCALE: pixscale, SMPLMODE:smplmode, DETMODE:detmode, ACQ_TIME:acq_time, INT_TIME:int_time, N_COADD:n_coadd, PAGAIN:pagain, PABANDW:pabandw, DETBIAS:detbias, EPERADU:eperadu, $          ; Detector information 
          DAOMODE: daomode, DCMODES: dcmodes, DLOOPGN: dloopgn, DWFSCFRQ: dwfscfrq, DWFSCBIN: dwfscbin, SAOMODE: saomode, SCMODES: scmodes, SLOOPGN: sloopgn, SWFSCFRQ: swfscfrq, SWFSCBIN: swfscbin, $   ; AO information                                                                                                                                                              ;
          LAM_CEN:lam_cen, BANDWIDTH:bandwidth, UBC_DXSP:ubc_dxsp, UBC_SXSP:ubc_sxsp, NIC_FSTP:nic_fstp, NIC_BEAM:nic_beam, LMIR_FW1:lmir_fw1, LMIR_FW2:lmir_fw2, LMIR_FW3:lmir_fw3, $
          LMIR_FW4:lmir_fw4, LMIR_AWL:lmir_awl, NOM_FW1:nom_fw1, NOM_FW2:nom_fw2, NOM_APW:nom_apw, PHA_FW1:pha_fw1, PHA_FW2:pha_fw2, PHA_IMG:pha_img,  $                                              ; Filter information
          SEEING: seeing, SMTTAU: smttau, LBTTEMP: lbttemp, WINDDIR: winddir, WINDSPD: windspd, $                                                                                                     ; Weather information
          DRS_VERS: drs_vers, DRS_DATE: drs_date, DATE_RED: date_red, BCK_MODE: bck_mode, BCK_SEL: bck_sel, N_XPIX: n_xpix, N_YPIX: n_ypix, BCK_NOD: bck_nod, IMG_MODE: img_mode, NO_BPM: no_bpm,$
          NO_DARK: no_dark, NO_FLAT: no_flat, PRECISION: precision, BCK_OB: bck_ob, BEAM_ID: beam_id, N_REJ: n_rej}                                                                                   ; Reduction information
                                                                                                                                                                   
                                                                                                                   
; 4. Read data that changes frame by frame
data = MRDFITS(lbti_file, 1, junk, /SILENT) 

; Add new fields for later
struct_add_field, data, 'ob_id',    INTARR(N_ELEMENTS(data.file_id))
struct_add_field, data, 'fwhmx_sx', FLTARR(N_ELEMENTS(data.file_id))
struct_add_field, data, 'fwhmy_sx', FLTARR(N_ELEMENTS(data.file_id))
struct_add_field, data, 'fwhmx_dx', FLTARR(N_ELEMENTS(data.file_id))
struct_add_field, data, 'fwhmy_dx', FLTARR(N_ELEMENTS(data.file_id))
struct_add_field, data, 'slope_sx', FLTARR(N_ELEMENTS(data.file_id))
struct_add_field, data, 'slope_dx', FLTARR(N_ELEMENTS(data.file_id))
IF NOT TAG_EXIST(data, 'dsloprms') THEN struct_add_field, data, 'dsloprms', FLTARR(N_ELEMENTS(data.file_id))  ; backwards compatibility
IF NOT TAG_EXIST(data, 'ssloprms') THEN struct_add_field, data, 'ssloprms', FLTARR(N_ELEMENTS(data.file_id))  ; backwards compatibility
IF NOT TAG_EXIST(data, 'spdthpos') THEN struct_add_field, data, 'spdthpos', FLTARR(N_ELEMENTS(data.file_id))  ; backwards compatibility


; Parse FITS header and DATA in the same structure
hdr_data = {HEADER: header, DATA: data}

RETURN, img_in
END
