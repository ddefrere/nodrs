PRO PREP_PS, BOLD=bold, RUNUTS=runuts,VERYBOLD=verybold

IF KEYWORD_SET(bold) THEN BEGIN
  !P.CHARSIZE=1.3
  !P.CHARTHICK=3.5
  !X.THICK=5.0
  !Y.THICK=5.0
  !P.THICK=5.0
ENDIF ELSE IF KEYWORD_SET(verybold) THEN BEGIN
  !P.CHARSIZE=1.5
  !P.CHARTHICK=6.0
  !X.THICK=5.
  !Y.THICK=5.
  !P.THICK=5.
ENDIF ELSE IF KEYWORD_SET(runuts) THEN BEGIN
  !P.CHARSIZE=1.5
  !P.CHARTHICK=6.0
  !X.THICK=10.
  !Y.THICK=10.
  !P.THICK=10.
ENDIF ELSE BEGIN
  !P.CHARSIZE=1.2
  !P.CHARTHICK=1.5
  !X.THICK=2
  !Y.THICK=2
  !P.THICK=2
ENDELSE

!P.FONT=0    ; Must be set to 0 to accept /TIMES
SET_PLOT,'ps'
END
