PRO BATCH_MERAK
  date = '131225'
  ; Reduce PSF
  ;nodrs, DATE='131225', DATA_IDX=[2222,2621], TGT_NAME='PSF'
  ; Reduce target
  nodrs, DATE=date, DATA_IDX=[197,2221], TGT_NAME='merak', /ADI
END