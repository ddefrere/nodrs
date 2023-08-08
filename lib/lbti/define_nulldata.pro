;+
; NAME: DEFINE_NULLDATA
;
; PURPOSE:
;   Simple procedure to define the array structure that will be saved in the L1 sumamry file
;
; OUTPUT
;   data :  the data structure
;
; MODIFICATION HISTORY:
;   Version 1.0, 02-DEC-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (adapted from former routine 'nomic_null.pro')
;   Version 1.1, 04-DEC-2015, DD: added output PLC_SPEC and PIXSCALE
;   Version 1.2, 23-DEC-2015, DD: Added precision and image mode to output data + corrected FIT_MODE
;   Version 1.3, 03-FEB-2015, DD: hanged a couple of values from FLOAT to DOUBLE (+ added int_err)
;   Version 1.4, 05-FEB-2015, DD: Added pointing ID
;   Version 1.5, 05-NOV-2016, DD: Added null offset
;   Version 1.6, 12-JAN-2017, DD: New null structure
;   Version 1.7, 27-MAR-2017, DD: Added frame mode
;   Version 1.8, 27-JUL-2017, DD: Background floor mode
;   Version 1.8, 31-JUL-2018, DD: Added DITH_PER to output structure

PRO DEFINE_NULLDATA, data

; Structure with null results
null_str = {NAS: 0D, NAS_ERR_LOW: 0D, NAS_ERR_SUP: 0D, MU: 0., MU_ERR_LOW: 0., MU_ERR_SUP: 0., SIG: 0., $
            SIG_ERR_LOW: 0., SIG_ERR_SUP: 0., KDRK: 0., KDRK_ERR_LOW: 0., KDRK_ERR_SUP: 0., CHI2: 0., CHI2_ERR_LOW: 0., CHI2_ERR_SUP: 0.}
             
null_star = {NULL_OPT: null_str, NULL_BAY: null_str, NULL_AVG: null_str, NULL_MOD: null_str}

data = {OBJNAME:'', PID: 0L, FLAG:'', LBT_UTC:'', LBT_LST:'', LBT_RA:'', LBT_DEC:'', LBT_ALT:0D, LBT_AZ:0D, LBT_PARA:0D,$                         ; Target info
        DATE_OBS:' ', MJD_OBS:0D, FILE_ID:0, BCK_MODE: 0, BFL_MODE: 0, FIT_MODE: 0, FLX_MODE: 0, FRA_MODE: 0, IMG_MODE: 0, OB_MODE: 0, PRECISION: 0, NOD_ID: 0, CHOP_ID: 0, OB_ID: 0, PT_ID: 0, $  ; File, mode, and timing info
        INSTRUM:'', SMPLMODE:0, DETMODE:0, ACQ_TIME:0., INT_TIME:0., N_COADD:0., PAGAIN:0., PABANDW:0., DETBIAS:0., DITH_PER: 0, $
        EPERADU: 0., N_XPIX: 0., N_YPIX: 0., PIXSCALE: 0D, $                                                                                      ; Detector and integration information
        XCEN: 0., YCEN: 0., APER_RAD: 0., BCK_IRAD: 0., BCK_ORAD: 0., $                                                                           ; Reduction information
        LAM_CEN: 0D, BANDWIDTH: 0D, UBC_DXSP:'', UBC_SXSP:'', NIC_FSTP:'', NIC_BEAM:'', LMIR_FW1:'', LMIR_FW2:'', LMIR_FW3:'', $
        LMIR_FW4:'', NOM_FW1:'', NOM_FW2:'', NOM_APW:'', PHA_FW1:'', PHA_FW2:'', PHA_IMG:'', NIL_DIC:'',  $                                       ; Filter information
        DAOMODE: '', DAOSTRHL:0., DCMODES: 0, DLOOPON:0., DLOOPGN:0., DWFSCFRQ:0., DWFSCBIN:0., SAOMODE: '', SAOSTRHL: 0., SCMODES: 0, SLOOPON: 0., $
        SLOOPGN: 0., SWFSCFRQ: 0., SWFSCBIN: 0., $                                                                                               ; AO information
        PCCLOSED: 0., PLC_WAV: 0., PLC_SPEC: ' ', PCPLSP: 0., PCTIPSP: 0., PCTLTSP: 0., SPC_PIST: 0., SPC_EL: 0., SPC_AZ: 0., FPC_PISTM: 0., FPC_PISTS: 0., $
        FPC_AZM: 0., FPC_AZS: 0., FPC_ELM: 0., FPC_ELS:0., PCPHMEAN: 0., PCPHMEAN_ERR: 0., PCPHSTD: 0., PCPHSTD_ERR: 0., PCMSNR: 0., $            ; Phasecam information
        SEEING:0., SMTTAU:0., LBTTEMP: 0., WINDDIR: 0., WINDSPD: 0., BCKG_AVG: 0., BCKG_RMS: 0., BCKG_ERR: 0., BIAS_UNC: 0., BIAS_COR: 0., BIAS_ERR: 0., $
        NULL_STAR: null_star, NULL_MEAS: 0D, NULL_ERR: 0D, NULL_OFFSET: 0D, NULL_SNR: 0., NULL_AVG: 0., NULL_RMS: 0., NSC_PHAVG: 0., NSC_PHRMS: 0., NSC_KDRK: 0., NSC_CHI2: 0., NSC_PHAVG_ERR: 0., NSC_PHRMS_ERR: 0., NSC_KDRK_ERR: 0., NSC_CHI2_ERR: 0., $
        BCKG_MEAS: 0., BCKG_EMEAS: 0., PHOTDX_AVG: 0., PHOTDX_RMS: 0., PHOTDX_SNR: 0., PHOTSX_AVG: 0., PHOTSX_RMS: 0., PHOTSX_SNR: 0., RMS_OPD: 0D, RMS_TT: [0.,0.], PHOT_AVG:0., PHOT_ERR:0., INT_ERR: 0.,$
        NFR_IN: 0, NFR_OB: 0, NFR_REJ: 0, N_FRBCK: 0, NOD_FRQ: 0., QUA_FLG: 0}                                                                    ; Results of the reduction

END