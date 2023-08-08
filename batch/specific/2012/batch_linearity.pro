;
; +
; PURPOSE:
;   Batch routine to plot linearity data.
;   Set reduce to re-reduce the data (data restored from disk otherwise).
;
; LAST MODIFICATION:
;   Version 1.0,  04-MAY-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu

PRO BATCH_LINEARITY, date, DATA_IDX=data_idx, DARK_IDX=dark_idx

; Running parameter
plot      = 1
info      = 3
aper_rad  = 0.
obstype   = 0.     ; 2 is for photometric data (1 for interferometric data)
no_gpu    = 1.     ; Damned, gpuFFT does not work with the free version!!!
no_save   = 0
no_center = 1
skip_cal  = 0
masterlog = 0

; Date to reduce
tgt_name = 'DOME'
flag     = 'LIN'
 
; Assign nod position (for noise calculation)
;FOR i=0, 39 DO LBTI_FITSUTIL, date, [data_idx[0],data_idx[0]]+i*10, LBT_LXOS=1, /RAW_DATA
;FOR i=0, 99 DO LBTI_FITSUTIL, date, [003010,003010]+i*10, LBT_LXOS=1, /RAW_DATA
;LBTI_FITSUTIL, date, [126795,127594], OBSTYPE=0, /RAW_DATA

xc       = [80]
yc       = [190]
aper_rad = 1.
bck_irad = 2.
bck_orad = 12.
LBTI_DRS, date, DATA_IDX=data_idx, BCKG_IDX=bckg_idx, DARK_IDX=dark_idx, FLAT_IDX=flat_idx, NOD_IDX=nod_idx,$                             ; Date to be reduced
          ;BCKG_PATH=bckg_folder, $;DATA_PATH=data_folder, DARK_PATH=dark_folder, FLAT_PATH=flat_folder, $                                ; Path to data directories
          XCEN=xc, YCEN=yc, OVERLAP=overlap, TGT_NAME=tgt_name, FLAG=flag, OBSTYPE=obstype, $
          LAMBDA_CEN=lambda, BANDWIDTH=bandwidth, NOM_FW2=nom_fw2, $                                                                      ; Optional input keywords superseeding keywords in the header
          APER_RAD=aper_rad, BCK_IRAD=bck_irad, BCK_ORAD=bck_orad, N_BIN=n_bin, N_CLIP=n_clip, NO_GPU=no_gpu, NO_CENTER=no_center,  $     ; Optional input keywords critical for data reduction
          LOG_FILE=log_file, INFO=info, MASTERLOG=masterlog, MOVIE=movie, PLOT=plot, NO_SAVE=no_save, SKIP_CAL=skip_cal                   ; Optional input keywords for output information
END
