;
; NAME: READ_TABLE
;
; PURPOSE:
;   Input of the data from an ASCII file with up to ten columns from a file on
;   disk.
;
; CATEGORY:
;   Input / output routines.
;
; CALLING SEQUENCE:
;   READ_TABLE, filename, A0 [, A1, A2, A3, A4, A5, A6, A7, A8, A9, AA, AB, AC, AD, AE, AF]
;
; INPUTS:
;   filename:
;     String variable or string containing the name of the file with the
;     experimental data. The format is an ASCII file with columns of floats.
;
; KEYWORDS:
;   FIRST, LAST
;     Number of first, resp. last line to be read. If not specified, the last
;     line is computed from the input file size. The line counter starts with 0,
;     not with 1!
;
;   MAXSIZE:
;     Initial size of arrays to be allocated. If not specified, the value is set
;     to 1E+6.
;
;   N_COLUMNS:
;     When set to the number of columns in the list, collect all columns in A0.
;     Works with SKIP. Does not work with a mixture of STRING and FLOAT columns
;
;   SEPARATOR:
;     Specifies the separator between columns in the table. Default is SPACE.
;     Set to '&', ',', '|' or whatever original idea your collegue has come up with.
;     Non-standard characters: 'TAB'
;
;   SKIP:
;     Array with the numbers of the columns to be skipped.
;
;   STRING_ARRAY:
;     Array with the numbers of the output arrays that are to be read in STRING format.
;     Remember to keep track of input columns to be skipped
;
; OUTPUTS:
;   A0,...,AF:
;     FLOAT arrays containing the respective columns.
;
; COMMON BLOCKS:
;   None.
;
; SIDE EFFECTS:
;   None.
;
; RESTRICTIONS:
;   None.
;
; SUBROUTINES:
;   None.
;
; PROCEDURE:
;
;   IDL> READ_TABLE,'myfile.dat', A0,A1,A2,A3,A4,A5, $
;                   SKIP=[0,2,4,6], FIRST=501, LAST=10500, MAXSIZE=10000
;   will read entries 501 through 10500 of columns 1, 3, 5, 7, 8, and 9 from
;   the ASCII file 'myfile.dat', and write them into separate FLOAT arrays
;   A0,...,A5.
;
; MODIFICATION HISTORY:
;   Created  06-DEC-1996, by Roland den Hartog, Aurora, ESTEC SA division, NL
;   Modified 16-FEB-1998, RdH: TAB detection and filesize detection implemented
;   Modified 06-DEC-1999, RdH: Separator specification and collection into one array
;   Modified 25-FEB-2004, RdH: String columns can be identified and read
;   Modified 07-APR-2004, RdH: In line 153, only leading blanks are removed, trailing blanks sometimes have a meaning
;   Modified 28-APR-2005, RdH: Solved problem with redundant declared columns
;   Modified 18-OCT-2005, RdH: Extended number of columns that can be read to 16. Removed duplication of inputs due to empty lines in input file.
;   Modified 05-DEC-2005, RdH: Repaired small bug involving definition of parameter 'str'

; CALL:
PRO READ_TABLE, filename, A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, AA, AB, AC, AD, AE, AF, $
  DOUBLE=double, ERROR=err, FIRST=first, LAST=last, MAXSIZE=maxsize, N_COLUMNS=n_columns, $
  SEPARATOR=separator, SILENT=silent, SKIP=skip, STRING_ARRAY=string_array, TEST=test

; ERROR HANDLING:
ON_ERROR,2                  ; Return to caller on error detection

; CODE:

; Open input file
OPENR, unit, filename, /GET_LUN, ERROR=err
IF err NE 0 THEN BEGIN
  IF NOT KEYWORD_SET(silent) THEN PRINT, !ERROR_STATE.MSG
  RETURN
ENDIF

; Set conversion from position in array to time of measurement
IF N_ELEMENTS(maxsize) LE 0 THEN maxsize=1E+6
IF N_ELEMENTS(first) LE 0 THEN first=0
IF N_ELEMENTS(last) LE 0 THEN BEGIN
  f = FSTAT(unit)           ; f is a structure defined by FSTAT
  last = f.size
ENDIF

IF N_ELEMENTS(separator) LE 0 THEN separator=' '
IF STRPOS(separator,'TAB',0) GE 0 OR STRPOS(separator,'tab',0) GE 0 THEN separator=STRING(9B)

Narr=N_PARAMS(0)-1
IF Narr LT 1 THEN MESSAGE,'No output arrays specified.'

; Allocate memory space
Nskp=N_ELEMENTS(skip)
Nlen=(last-first+1) < maxsize
str=INTARR(16)
IF N_ELEMENTS(n_columns) GT 0 THEN Nalloc=N_columns $
ELSE BEGIN
  IF N_ELEMENTS(string_array) GT 0 THEN BEGIN
    s=STRARR(Nlen)
    str[string_array]=1
  ENDIF
  Nalloc=Narr+Nskp
  IF KEYWORD_SET(double) THEN BEGIN
    IF Nlen LE 1 THEN a=0D0 ELSE a=DBLARR(Nlen)
  ENDIF ELSE BEGIN
    IF Nlen LE 1 THEN a=0. ELSE a=FLTARR(Nlen)
  ENDELSE
  IF str[0] THEN A0=s ELSE A0=a
  IF Narr GT  1 THEN IF str[ 1] THEN A1=S ELSE A1=A
  IF Narr GT  2 THEN IF str[ 2] THEN A2=S ELSE A2=A
  IF Narr GT  3 THEN IF str[ 3] THEN A3=S ELSE A3=A
  IF Narr GT  4 THEN IF str[ 4] THEN A4=S ELSE A4=A
  IF Narr GT  5 THEN IF str[ 5] THEN A5=S ELSE A5=A
  IF Narr GT  6 THEN IF str[ 6] THEN A6=S ELSE A6=A
  IF Narr GT  7 THEN IF str[ 7] THEN A7=S ELSE A7=A
  IF Narr GT  8 THEN IF str[ 8] THEN A8=S ELSE A8=A
  IF Narr GT  9 THEN IF str[ 9] THEN A9=S ELSE A9=A
  IF Narr GT 10 THEN IF str[10] THEN AA=S ELSE AA=A
  IF Narr GT 11 THEN IF str[11] THEN AB=S ELSE AB=A
  IF Narr GT 12 THEN IF str[12] THEN AC=S ELSE AC=A
  IF Narr GT 13 THEN IF str[13] THEN AD=S ELSE AD=A
  IF Narr GT 14 THEN IF str[14] THEN AE=S ELSE AE=A
  IF Narr GT 15 THEN IF str[15] THEN AF=S ELSE AF=A
  a=0 & s=0
ENDELSE

; Determine which columns to get
get=INTARR(Nalloc)+1
IF Nskp GT 0 THEN BEGIN
  get(skip)=0
ENDIF

; Input loop
item=0L
iline=0L
line=' '

temp=STRARR(512)
WHILE NOT EOF(unit) DO BEGIN

; Get a line from the input file in STRING format
  READF,unit,line,FORMAT='(1024A)'
  blank = STRTRIM(line,2) EQ ''

  IF iline GE first AND iline LE last AND NOT blank THEN BEGIN

;   Obtain the columns from the input line
    cnt=0L
    eol=0B
    line=separator+STRTRIM(line,1)   ; Remove leading blanks, except first
    len=STRLEN(line)
    IF KEYWORD_SET(test) THEN PRINT,'new: ',line
    FOR col=0L,Nalloc-1L DO BEGIN
      i=0 & j=0L
      next:
      i=STRPOS(line,separator,i)+1  ; Search for the separator
      IF i GT 1024 THEN GOTO,jump
      IF (STRMID(line,i,1) EQ separator) THEN GOTO,next ELSE BEGIN
        a=STRMID(line,i,1024)
        IF get[col] GT 0 AND STRLEN(a) GT 0 THEN BEGIN
          ;PRINT,col,cnt,i,STRLEN(a),' ',a
          j=STRPOS(a,separator)
          IF j GT 0 THEN temp[cnt]=STRMID(a,0,j) ELSE IF eol THEN temp[cnt]='' ELSE BEGIN & temp[cnt]=a & eol=1B & ENDELSE
          cnt=cnt+1L
        ENDIF
      ENDELSE
      line=a
      IF KEYWORD_SET(test) THEN PRINT,i,j,cnt,': ',a,':',STRLEN(a),STRPOS(a,separator),':',temp[cnt]
    ENDFOR
    jump:

    IF N_ELEMENTS(n_columns) GT 0 THEN BEGIN
      IF cnt GT 0 THEN IF item LE 0 THEN A0=temp[0:cnt-1] ELSE A0=[A0,temp[0:cnt-1]]
    ENDIF ELSE BEGIN
      IF Nlen GT 1 THEN BEGIN
        IF str[0] THEN A0[item]=temp[0] ELSE BEGIN
          READS, temp[0], value
          A0[item]=value
        ENDELSE
        IF Narr GT 1 THEN IF str[1] THEN A1[item]=temp[1] ELSE BEGIN
          READS, temp[1], value
          A1[item]=value
        ENDELSE
        IF Narr GT 2 THEN IF str[2] THEN A2[item]=temp[2] ELSE BEGIN
          READS, temp[2], value
          A2[item]=value
        ENDELSE
        IF Narr GT 3 THEN IF str[3] THEN A3[item]=temp[3] ELSE BEGIN
          READS, temp[3], value
          A3[item]=value
        ENDELSE
        IF Narr GT 4 THEN IF str[4] THEN A4[item]=temp[4] ELSE BEGIN
          READS, temp[4], value
          A4[item]=value
        ENDELSE
        IF Narr GT 5 THEN IF str[5] THEN A5[item]=temp[5] ELSE BEGIN
          READS, temp[5], value
          A5[item]=value
        ENDELSE
        IF Narr GT 6 THEN IF str[6] THEN A6[item]=temp[6] ELSE BEGIN
          READS, temp[6], value
          A6[item]=value
        ENDELSE
        IF Narr GT 7 THEN IF str[7] THEN A7[item]=temp[7] ELSE BEGIN
          READS, temp[7], value
          A7[item]=value
        ENDELSE
        IF Narr GT 8 THEN IF str[8] THEN A8[item]=temp[8] ELSE BEGIN
          READS, temp[8], value
          A8[item]=value
        ENDELSE
        IF Narr GT 9 THEN IF str[9] THEN A9[item]=temp[9] ELSE BEGIN
          READS, temp[9], value
          A9[item]=value
        ENDELSE
        IF Narr GT 10 THEN IF str[10] THEN AA[item]=temp[10] ELSE BEGIN
          READS, temp[10], value
          AA[item]=value
        ENDELSE
        IF Narr GT 11 THEN IF str[11] THEN AB[item]=temp[11] ELSE BEGIN
          READS, temp[11], value
          AB[item]=value
        ENDELSE
        IF Narr GT 12 THEN IF str[12] THEN AC[item]=temp[12] ELSE BEGIN
          READS, temp[12], value
          AC[item]=value
        ENDELSE
        IF Narr GT 13 THEN IF str[13] THEN AD[item]=temp[13] ELSE BEGIN
          READS, temp[13], value
          AD[item]=value
        ENDELSE
        IF Narr GT 14 THEN IF str[14] THEN AE[item]=temp[14] ELSE BEGIN
          READS, temp[14], value
          AE[item]=value
        ENDELSE
        IF Narr GT 15 THEN IF str[15] THEN AF[item]=temp[15] ELSE BEGIN
          READS, temp[15], value
          AF[item]=value
        ENDELSE
      ENDIF ELSE BEGIN
        IF str[0] THEN A0=temp[0] ELSE READS,temp[0],A0
        IF Narr GT  1 THEN IF str[ 1] THEN A1=temp[ 1] ELSE READS,temp[ 1],A1
        IF Narr GT  2 THEN IF str[ 2] THEN A2=temp[ 2] ELSE READS,temp[ 2],A2
        IF Narr GT  3 THEN IF str[ 3] THEN A3=temp[ 3] ELSE READS,temp[ 3],A3
        IF Narr GT  4 THEN IF str[ 4] THEN A4=temp[ 4] ELSE READS,temp[ 4],A4
        IF Narr GT  5 THEN IF str[ 5] THEN A5=temp[ 5] ELSE READS,temp[ 5],A5
        IF Narr GT  6 THEN IF str[ 6] THEN A6=temp[ 6] ELSE READS,temp[ 6],A6
        IF Narr GT  7 THEN IF str[ 7] THEN A7=temp[ 7] ELSE READS,temp[ 7],A7
        IF Narr GT  8 THEN IF str[ 8] THEN A8=temp[ 8] ELSE READS,temp[ 8],A8
        IF Narr GT  9 THEN IF str[ 9] THEN A9=temp[ 9] ELSE READS,temp[ 9],A9
        IF Narr GT 10 THEN IF str[10] THEN AA=temp[10] ELSE READS,temp[10],AA
        IF Narr GT 11 THEN IF str[11] THEN AB=temp[11] ELSE READS,temp[11],AB
        IF Narr GT 12 THEN IF str[12] THEN AC=temp[12] ELSE READS,temp[12],AC
        IF Narr GT 13 THEN IF str[13] THEN AD=temp[13] ELSE READS,temp[13],AD
        IF Narr GT 14 THEN IF str[14] THEN AE=temp[14] ELSE READS,temp[14],AE
        IF Narr GT 15 THEN IF str[15] THEN AF=temp[15] ELSE READS,temp[15],AF
      ENDELSE
    ENDELSE
    item=item+1
  ENDIF
  iline=iline+1L
ENDWHILE
ntot=item-1

; Truncate the output arrays
IF N_ELEMENTS(n_columns) LE 0 AND Nlen GT 1 THEN BEGIN
  A0=A0[0:ntot]
  IF Narr GT  1 THEN A1=A1[0:ntot]
  IF Narr GT  2 THEN A2=A2[0:ntot]
  IF Narr GT  3 THEN A3=A3[0:ntot]
  IF Narr GT  4 THEN A4=A4[0:ntot]
  IF Narr GT  5 THEN A5=A5[0:ntot]
  IF Narr GT  6 THEN A6=A6[0:ntot]
  IF Narr GT  7 THEN A7=A7[0:ntot]
  IF Narr GT  8 THEN A8=A8[0:ntot]
  IF Narr GT  9 THEN A9=A9[0:ntot]
  IF Narr GT 10 THEN AA=AA[0:ntot]
  IF Narr GT 11 THEN AB=AB[0:ntot]
  IF Narr GT 12 THEN AC=AC[0:ntot]
  IF Narr GT 13 THEN AD=AD[0:ntot]
  IF Narr GT 14 THEN AE=AE[0:ntot]
  IF Narr GT 15 THEN AF=AF[0:ntot]
ENDIF

; Close the input file
CLOSE,unit
FREE_LUN,unit
END
