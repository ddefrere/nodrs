;+
; NAME: DEFINE_HDRDATA
;
; PURPOSE:
;   Simple procedure to define the array structure that is used to read the L0 header
;
; OUTPUT
;   data : the data structure
;
; MODIFICATION HISTORY:
;   Version 1.0, 02-DEC-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.1, 21-JAN-2016, DD: added FWHM_X, FWHM_Y, and SLOPE to output structure
;   Version 1.2, 21-JAN-2016, DD: added PHASECam information for beam 2
;   Version 1.3, 03-JUL-2016, DD: added central value
;   Version 1.4, 20-OCT-2016, DD: added INTERFAC
;   Version 1.5, 26-OCT-2016, DD: added DSLOPRMS and SLOPRMS
;   Version 1.6, 08-NOV-2016, DD: added D/SSAFESKP and S/DLOPMAX 
;   Version 1.7, 30-JAN-2017, DD: added SPDTHPOS
;   Version 1.8, 29-JUL-2018, DD: BCK_AVG and BCK_RMS are now 64-element vector (64 is the maximum number of possible channels)

PRO DEFINE_HDRDATA, data

data = { FILENAME:' ', DATE_OBS:' ', LBT_UTC:' ', LBT_LST:' ', LBT_RA:' ', LBT_DEC:' ', LBT_ALT:0., LBT_AZ:0., LBT_PARA:0.,  MJD_OBS: 0D, OBJNAME:' ', TIME_END:' ', $; File and target info
		     FLAG:' ', DATATYPE: 0, OBSTYPE:0, $
         BCK_NOD: 0., BCK_AVG: DBLARR(64), BCK_RMS: DBLARR(64), CFG_ID: 0, CHP_ID: 0, FILE_ID: 0L, FR_ID: 0, FWHM_X: [0.,0.], FWHM_Y: [0.,0.], N_FRBCK: 0, N_XPIX: 0, N_YPIX: 0, N_FR: 0, ROI: [0,0,0,0], $
         NOD_FRQ: 0., NOD_ID: 0, PT_ID: 0, OVERLAP: 0, PID: 0, SLOPE: [0.,0.],  XCEN: [0.,0.], YCEN: [0.,0.], $                                                                     ; Reduction information
         INSTRUM: ' ', CV: 0., SMPLMODE: 0., DETMODE: 0., ACQ_TIME: 0D, INT_TIME: 0D, N_COADD: 0., PAGAIN: 0L, PABANDW: 0L, INTERFAC:' ', PACTCDLY: 0L, DETBIAS: 0., EPERADU:0., $  ; Detector and integration information
         LAM_CEN: 0., BANDWIDT: 0., UBC_DXSP:' ', UBC_SXSP:' ', NIC_FSTP:' ', NIC_BEAM:' ', LMIR_FW1:' ', LMIR_FW2: ' ', LMIR_FW3:' ', LMIR_FW4:' ', LMIR_AWL:' ', $
         NOM_FW1:' ', NOM_FW2:' ', NOM_APW:' ', PHA_FW1:' ', PHA_FW2:' ', PHA_IMG:' ', NIL_DIC:' ',  $                                                                ; Filter information
         DAOMODE:' ', DAOSTRHL: 0., DCMODES: 0, DLOOPGN:' ', DLOOPON:0, DWFSCFRQ: 0., DWFSCBIN: 0, DSLOPRMS: 0., DSLOPMAX: 0., DSAFESKP: 0., SAOMODE:' ', SAOSTRHL: 0., SCMODES: 0, SLOOPGN:' ',$
         SLOOPON:0, SWFSCFRQ: 0., SWFSCBIN: 0, SSLOPRMS: 0., SSLOPMAX: 0., SSAFESKP: 0., $  ; AO information
         PCCLOSED: 0, SPC_PIST: 0., SPC_EL: 0., SPC_AZ: 0., FPC_PISTM: 0., FPC_PISTS: 0., FPC_AZM: 0., FPC_AZS: 0., FPC_ELM: 0., FPC_ELS: 0.,$
         PCPLSP: 0., PCTIPSP: 0., PCTLTSP: 0., SPDTHPOS: 0., PCPHMEAN: 0., PCPHSTD: 0., PCPHMCOS: 0., PCPHMSIN: 0., PCFJMPS: 0., PCMSNR: 0., PCPHMEAN2: 0., PCPHSTD2: 0., PCPHMCOS2: 0., PCPHMSIN2: 0., PCMSNR2: 0., $  ; Phasecam information
         SEEING: 0., SMTTAU: 0., LBTTEMP: 0., WINDDIR: 0., WINDSPD: 0. $                                                                                              ; Weather information
        }                                                                                                                      
END
