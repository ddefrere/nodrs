; NAME: PLOTXY
;
; PURPOSE:
;   Multiple purpose routine to plot data in a diagram with an x and y axis.
;   Adaptable for the following circumstances:
;   - x and y points as well as images
;   - in a screen window or in a postscript file
;   - single or multiple plots on a page
;   - histograms, curves or scatterplots, with a large variety of symbols
;   - full range of titles, subtitles and insets
;
; CATEGORY:
;   Array processing and Miscellaneous Math routines
;
; CALLING SEQUENCE:
;   PlotXY, x [, y]
;
; INPUTS:
;   x, y:
;     In the case of data pairs: one-dimensional arrays containing the data to
;     be displayed. If only one array is specified, it is treated as y-data.
;     The x-data is then generated using INDGEN.
;     In the case of an image, x must be a two-dimensional array.
;
; KEYWORDS:
;   FILENAME:
;     Name of the output file.
;
;   PS:
;     Send output to a postscript file. If a valid filename has been given to
;     FILENAME the output file will have this name, otherwise the name is
;     'plotxy.ps'
;
;   EPS:
;     If set, the format of the PS file is Encapsulated Post Script.
;
;   HP:
;     Send output yo a file in HP-GL format, like PS a pen-driven format, but
;     unlike PS, recognized by MS Word.
;
;   TIFF:
;     Store output to a file in TIFF format. If a valid filename has been given
;     to FILENAME the output file will have this name, otherwise the name is
;     'plotxy.tif'. Cannot be used simultaneously with PSFILE because the
;     creation of the TIFF file is TV driven.
;
;   BMP:
;     Store output to a file in TIFF format. If a valid filename has been given
;     to FILENAME the output file will have this name, otherwise the name is
;     'plotxy.tif'. Cannot be used simultaneously with PSFILE because the
;     creation of the TIFF file is TV driven.
;
;   USD:
;     If set, write bitmap up-side down. Used for e.g. MS applications.
;
;   INV:
;     If set, invert bitmap, i.e. black background, white plotlines.
;
;   MIR:
;     If set, mirror bitmap.
;
;   TITLE:
;     Title to appear above plot, e.g. identification of the image.
;
;   SUBTITLE:
;     Subtitle to appear below the X-axis.
;
;   XTITLE:
;     Title to appear below the x-axis.
;
;   YTITLE:
;     Title to appear left of y-axis.
;
;   INSET_LL, INSET_UL, INSET_LR, INSET_UR:
;     Generate insets in resp. the lower left, upper left, lower right and
;     upper right corner of the diagram. Consecutive calls to PLOTXY with the
;     ADD keyword on, generate separated lines of insets.
;     The special character !P has been defined to plot either a line or a
;     symbol specified as outlined under SYMBOL,
;     e.g. PLOTXY, ..., INSET_LL='!P-3 : a dash-dotted line'
;
;   XRANGE, YRANGE:
;     Restrict the range of the x-axis resp. y-axis to the specified values.
;     Both keywords can be set either to a 2-element array or a 4-element array.
;     XRANGE=[xmin,xmax] contains the minimum and maximum x-value to appear on
;     the axes.
;     XRANGE=[xmin0,xmax0,xmin1,xmax1] forces a different scaling on the lower
;     and upper axes, given resp. by the former and latter pair of values.
;
;   SYMBOL:
;     Defines the symbol used to represent the data.
;
;     SYMBOL = 0 : solid line (default)
;     SYMBOL < 0 : lines
;     SYMBOL > 0 : markers
;
;     Value:       Result:
;
;                  LINE-TYPES
;     -5           Long Dashes
;     -4           Dash Dot Dot Dot
;     -3           Dash Dot
;     -2           Dashed Line
;     -1           Dotted Line
;      0           Solid Line
;
;                  VERTEX-CONNECTED SYMBOLS:
;      1           .
;      2           +
;      3           x
;      4           *
;      5           arrow, up
;      6           arrow, down
;      7           arrow, left
;      8           arrow, right
;
;                  OPEN SYMBOLS:
;     10           triangle, up
;     11           triangle, down
;     12           triangle, left
;     13           triangle, right
;     14           square
;     15           diamond
;     16           vertex 5-point star
;     17           convex 5-point star
;     18           circle, low resolution
;     19           circle, high resolution
;
;     +10          FILLED SYMBOLS (SOLID BLACK)
;     +20          FILLED SYMBOLS (SHADED)
;     >40          Freely defined symbols, via previous call to USERSYM
;
;    XLOG, YLOG:
;      If set, the respective axis will be drawn in logarithmic scaling.
;      Not compatible with image data.
;
;    SCATTER:
;      If set, produce a scatter plot (same as setting SYMBOL to 1).
;      Takes priority over SYMBOL.
;
;    HISTO:
;      If set, produce an histogram. Takes priority over SCATTER.
;
;    GRID:
;      Add a grid to the plot. By default the grid has solid bars, but the
;      linestyle can be changed according to the same conventions as SYMBOL.
;
;    NODATA:
;      If set, suppress the plotting of data points.
;
;    NOAXES:
;      If set, suppress the axes plus annotations
;
;    NOTICKS, NOXTICKS, NOYTICKS:
;      If set, suppress the tick marks for resp. both, the X- or the Y-axis.
;
;    NOTICKV, NOXTICKV, NOYTICKV:
;      If set, suppress the annotations for resp. both, the X- or the Y-axis.
;
;    FALSE / BW:
;      Plot in false colors resp. black & white.
;
;    GREY:
;      If set to a value > 0 and < 256, the colour of the characters, lines
;      etc. is shifted progressively towards white. GREY=255 corresponds to
;      white.
;
;    CONTOUR:
;      Make a contour plot. Works only with images.
;
;    NLEVELS:
;      Number of contour levels. Default is 4.
;
;    LEVELS:
;      Actual values of the levels, in image greyscale values.
;
;    WINDOW:
;      4-element array that contains the coordinates [X0,Y0,X1,Y1] of resp. the
;      lower left (X0,Y0) and upper right (X1,Y1) corner of the plot window.
;      Units are in cm. for PS format, and in pixels for on screen display.
;
;    ERASE:
;      Indicates that previous plot must be erased. By default the plotting is
;      in NOERASE mode.
;
;    INIT:
;      If set, initialize new plot page, without doing any plotting. The size of
;      the plotpage is determined by the WINDOW keyword. Allows to set !P.MULTI
;      first.
;
;    MULTI:
;      Used in combination with INIT. Depending on its dimensions MULTI gives
;      the number of rows, (# rows, # columns) or the full definition of the
;      !P.MULTI graphic keyword. Overruled by subsequent use of the WINDOW
;      keyword. DOES NOT WORK AT THE MOMENT!
;
;    NEW:
;      Skip the page initialization and plot a single diagram.
;      The position of the diagram on the page can now be controlled by the
;      WINDOW keyword.
;
;    HOLD:
;      When generating a PS file, refrain from closing the device. To be used
;      together with the ADD keyword.
;
;    ADD:
;      Add a set of datapoints by overplotting, or add an additional inset.
;      Assumes that a previous plot is present in the output medium.
;
;    FINAL:
;      Close the output device, i.e. finish the PS file.
;
; OUTPUTS:
;   Window on screen or plot in PS format, containing the diagram with the data.
;
; CALLS:
;   PLOTSYM.PRO, external routine, modified by RdH, that uses USERSYM to
;   generate a set of the most common predefined symbols, that are not standard
;   in IDL.
;
; COMMON BLOCKS:
;   PLOTXY: passes on information between consecutive calls to PLOTXY
;
; SIDE EFFECTS:
;   None.
;
; RESTRICTIONS:
;   None.
;
; EXAMPLE:
;   Inverting colors on the screen
;   IDL> LOADCT, 0, /SILENT
;   IDL> x=FINDGEN(1000) & y=SIN(x/100.)
;   IDL> PLOTXY, /INIT, /MSWIN, WINDOW=[0,0,500,400], /INV
;   IDL> PLOTXY, x, y, /NEW, WINDOW=[70,40,490,390], /NOERASE, COLOR=0.1
;
; MODIFICATION HISTORY:
;   Created 09-APR-1996, by Roland den Hartog, Aurora, ESTEC SA division, NL
;   Version 23-AUG-1996, by RdH: removed problems with positioning of window
;   Version 14-NOV-1996, by RdH: implemented XTITLE and YTITLE keywords
;   Version 14-APR-1997, by RdH: impl. SCATTER, HISTO, WINDOW keywords
;   Version 15-APR-1997, by RdH: impl. INIT, HOLD, NEW, ADD and FINAL kwds
;   Version 27-MAY-1997, by RdH: impl. different axes on either side of plot
;   Version 28-MAY-1997, by RdH: impl. _EXTRA and TIFF kwd
;   Version 25-JUN-1997, by RdH: impl. EPS and other bitmap kwds
;   Version 18-MAY-1998, by RdH: impl. SHOW image-display functionality
;   Version 26-MAY-1998, by RdH: improved calculation of INSET positions
;   Version 09-JUL-1998, by RdH: impl. HP format
;   Version 19-AUG-1998, by RdH: various improvements
;   Version 16-JAN-1999, by RdH: MS Windows default
;   Version 02-MAY-1999, by RdH: PLOTSYM procedure included to avoid problems
;                                with calls to this routine
;   Version 02-JUN-1999, by RdH: Changed size of plot window on screen
;   Version 01-JUL-1999, by RdH: Solved the colour problem: !P.COLOR=255 --> red!
;   Version 14-NOV-1999, by RdH: Refined tick suppression
;   Version 01-MAR-1999, by RdH: Problem with colours solved
;   Version 21-FEB-2001, by RdH: Problem with CONVERT_COORD solved
;   Version 21-SEP-2001, by RdH: JPEG output implemented, but it doesn't work properly yet
;   Version 12-OCT-2001, by RdH: Put an end to Didier's ordeal by making a full reset of all system parameters
;   Version 20-FEB-2002, by RdH: Fixed problems with grey shading
;   Version 25-APR-2002, by RdH: Solved contrast problem with tiff file generation
;   Version 15-APR-2003, by RdH: removed PLOT_INSET from internal routine list, as it is defined as a separate routine.
;   Version 08-DEC-2003, by RdH: cracked the !P.MULTI problem!
;   Version 19-JAN-2005, by RdH: cracked the inverse color problem!
;   Version 29-MAR-2011, by DDe: restored default plot title to 0 
;   Version 30-APR-2011, by DDe: added keywords XTHICK and YTHICK
;   Version 13-JUN-2013, by DDe: if NWINDOW is set together with INIT, overplot in the window #nwindow rather than creating a new one
;   Version 28-MAR-2016, by DDe: Added /TIMES to DEVICE's call + improved speed for ssh use
;
;-------------------------------------------------------------------------------

; SUBROUTINE
;   Define useful plotting symbols not in the standard !PSYM definitions.
;   After symbol has been defined with PLOTSYM, a plotting command should
;   follow with either PSYM = 8 or !P.PSYM = 8 (see USERSYM)

PRO PLOTSYM, psym, psize, FILL=fill, COLOR=color, THICK=thick
ON_ERROR,2

IF N_ELEMENTS(psym) LT 1 THEN BEGIN
     PRINT,'Syntax:  PLOTSYM, psym, [ size, /FILL ]'
     PRINT,'PSYM values:'
     PRINT,' 0 - circle, low accuracy'
     PRINT,' 1 - circle, high accuracy'
     PRINT,' 2 - up triangle'
     PRINT,' 3 - down triangle'
     PRINT,' 4 - left-pointing triangle'
     PRINT,' 5 - right-pointing triangle'
     PRINT,' 6 - square'
     PRINT,' 7 - diamond'
     PRINT,' 8 - up arrow'
     PRINT,' 9 - down arrow'
     PRINT,'10 - left arrow'
     PRINT,'11 - right arrow'
     PRINT,'12 - vertex-connected 5-pointed star'
     PRINT,'13 - convex-connected 5-pointed star'
  RETURN
ENDIF

IF N_ELEMENTS(psize) LT 1 THEN psize = 1 ELSE psize = psize > 0.1
IF KEYWORD_SET(fill) THEN fill=1 ELSE fill=0

CASE psym OF
 0: BEGIN                                     ; Circle, low accuracy
    ang = 2*!PI*FINDGEN(19)/18.               ; Get position every 20 deg
    xarr = COS(ang) & yarr = SIN(ang)
    END

 1: BEGIN                                     ; Circle, high accuracy
    ang = 2*!PI*FINDGEN(49)/48.               ; Get position every 5 deg
    xarr = COS(ang) & yarr = SIN(ang)
    END

 2: BEGIN                                     ; Up-pointing triangle
    xarr = [-1.,1.,0.,-1.] & yarr = [-1.,-1.,1.,-1.]
    END

 3: BEGIN                                     ; Down-pointing triangle
    xarr = [-1.,1.,0.,-1.] & yarr = [ 1.,1.,-1.,1.]
    END

 4: BEGIN                                     ; Left-pointing triangle
    xarr = [1.,1.,-1.,1.] & yarr = [-1.,1.,0.,-1.]
    END

 5: BEGIN                                     ; Right-pointing triangle
    xarr = [-1.,-1.,1.,-1.] & yarr = [-1.,1.,0.,-1.]
    END

 6: BEGIN                     ; Square
    xarr = [-1., 1., 1., -1., -1.]
    yarr = [-1., -1., 1., 1., -1.]
    END

 7: BEGIN                     ; Diamond
    xarr = [-1., 0., 1., 0., -1.]*SQRT(2.0)
    yarr = [ 0., 1., 0., -1., 0.]*SQRT(2.0)
    END

 8: BEGIN                                     ; Up arrow
    xarr = [0,0,.5,0,-.5]
    yarr = [0,2,1.4,2,1.4]
    END

 9: BEGIN                                     ; Down arrow
    xarr = [0,0,.5,0,-.5]
    yarr = [0,-2,-1.4,-2,-1.4]
    END

10: BEGIN                                     ; Left pointing arrow
    yarr = [0, 0, 0.5, 0, -.5]
    xarr = [0,-2,-1.4,-2,-1.4]
    END

11: BEGIN                                     ; Left pointing arrow
    yarr = [ 0, 0, 0.5, 0, -.5]
    xarr = [ 0, 2, 1.4, 2, 1.4]
    END

12: BEGIN                                     ; Vertex-connected star
    ang = (FINDGEN(6)*144. + 90.) / !RADEG    ; Define star angles every 144 deg
    xarr = 2.0*COS(ang) & yarr = 2.0*SIN(ang)
    END

13: BEGIN                                     ; Convex-connected star
    ang = (FINDGEN(6)*72.+90.) / !RADEG       ; Define star angles every 72 deg
    xarr = 2.0*COS(ang) & yarr = 2.0*SIN(ang)
    END

ELSE: MESSAGE,'Unknown plotting symbol value: psym'
ENDCASE

;IF fill NE 0 THEN thick=2 ELSE thick=1
USERSYM, xarr*psize, yarr*psize, FILL=fill, THICK=thick, COLOR=color

RETURN
END

;------------------------------------------------------------------------------
; SUBROUTINE to define plot symbols or linestyle

PRO SETSYM, symbol, THICK=thick
ON_ERROR, 2
IF N_ELEMENTS(thick) GT 0 THEN !P.THICK=thick
IF symbol GT 0 THEN BEGIN
  CASE symbol OF
    1: !P.PSYM=3
    2: !P.PSYM=1
    3: !P.PSYM=7
    4: !P.PSYM=2
    ELSE: !P.PSYM=8
  ENDCASE
  IF symbol GE 30 THEN color=0.5*!D.N_COLORS
  IF symbol GT 4 THEN BEGIN
    p=!P.SYMSIZE
    CASE symbol OF
      5:  PLOTSYM, 8, p, THICK=thick
      6:  PLOTSYM, 9, p, THICK=thick
      7:  PLOTSYM,10, p, THICK=thick
      8:  PLOTSYM,11, p, THICK=thick
      10: PLOTSYM, 2, p, THICK=thick
      11: PLOTSYM, 3, p, THICK=thick
      12: PLOTSYM, 4, p, THICK=thick
      13: PLOTSYM, 5, p, THICK=thick
      14: PLOTSYM, 6, p, THICK=thick
      15: PLOTSYM, 7, p, THICK=thick
      16: PLOTSYM,12, p, THICK=thick
      17: PLOTSYM,13, p, THICK=thick
      18: PLOTSYM, 0, p, THICK=thick
      19: PLOTSYM, 1, p, THICK=thick
      20: PLOTSYM, 2, p, /FILL, THICK=thick
      21: PLOTSYM, 3, p, /FILL, THICK=thick
      22: PLOTSYM, 4, p, /FILL, THICK=thick
      23: PLOTSYM, 5, p, /FILL, THICK=thick
      24: PLOTSYM, 6, p, /FILL, THICK=thick
      25: PLOTSYM, 7, p, /FILL, THICK=thick
      26: PLOTSYM,12, p, /FILL, THICK=thick
      27: PLOTSYM,13, p, /FILL, THICK=thick
      28: PLOTSYM, 0, p, /FILL, THICK=thick
      29: PLOTSYM, 1, p, /FILL, THICK=thick
      30: PLOTSYM, 2, p, /FILL, COLOR=color, THICK=thick
      31: PLOTSYM, 3, p, /FILL, COLOR=color, THICK=thick
      32: PLOTSYM, 4, p, /FILL, COLOR=color, THICK=thick
      33: PLOTSYM, 5, p, /FILL, COLOR=color, THICK=thick
      34: PLOTSYM, 6, p, /FILL, COLOR=color, THICK=thick
      35: PLOTSYM, 7, p, /FILL, COLOR=color, THICK=thick
      36: PLOTSYM,12, p, /FILL, COLOR=color, THICK=thick
      37: PLOTSYM,13, p, /FILL, COLOR=color, THICK=thick
      38: PLOTSYM, 0, p, /FILL, COLOR=color, THICK=thick
      39: PLOTSYM, 1, p, /FILL, COLOR=color, THICK=thick
    ELSE: RETURN
    ENDCASE
  ENDIF
ENDIF ELSE !P.LINESTYLE=ABS(symbol)

RETURN
END

;------------------------------------------------------------------------------

PRO PLOT_INSET, inset, xynorm, sign, prev, COLOR=color
ON_ERROR,2

xmin=xynorm(0,0) < xynorm(0,1) & xmax=xynorm(0,0) > xynorm(0,1)
ymin=xynorm(1,0) < xynorm(1,1) & ymax=xynorm(1,0) > xynorm(1,1)

charsize=!P.CHARSIZE
thick=!P.THICK
dx=!D.X_CH_SIZE*1.2*charsize & dy=!D.Y_CH_SIZE*1.2*charsize

i=1.5-sign[0]-0.5*sign[1]
IF sign[0] GT 0 THEN xp=xmin+dx ELSE xp=xmax-dx
IF sign[1] GT 0 THEN yp=ymin+(prev(i)+1.5)*dy $
                ELSE yp=ymax-(prev(i)+1.)*dy
prev(i)=prev(i)+1

!P.LINESTYLE=0
symsize=!P.SYMSIZE
i=STRPOS(inset,'!P')
IF i GE 0 THEN BEGIN
  IF i EQ 0 THEN line=STRMID(inset,i+4,256) $
            ELSE line=STRMID(inset,0,i-1)+STRMID(inset,i+4,256)
  s=FIX(STRMID(inset,i+2,2))
  SETSYM, s, THICK=thick
  IF s LE 0 THEN BEGIN
    !P.PSYM=0
    PLOTS, [xp,xp+sign[0]*3.8*dx], [yp,yp], /DEVICE, SYMSIZE=symsize;, COLOR=color
    xp=xp+sign[0]*0.1*dx
    !P.LINESTYLE=0
  ENDIF ELSE BEGIN
    psym=!P.PSYM
    PLOTS, xp+sign[0]*dx, yp+0.2*dy, /DEVICE, PSYM=psym, SYMSIZE=symsize;, COLOR=color
    line='  '+line
    !P.PSYM=0
  ENDELSE
ENDIF ELSE BEGIN
  line=inset
ENDELSE

charthick=!P.CHARTHICK
XYOUTS, xp, yp, line, ALIGN=0.5-sign[0]*0.5, /DEVICE, $
  CHARSIZE=charsize, CHARTHICK=charthick

RETURN
END

;===============================================================================
; MAIN:

PRO PlotXY, xin, yin, ADD=add, BMP=bmp, BW=bw, CHARSIZE=charsize, CHARTHICK=charthick, $
  COLOR=color, CONTOUR=contour, EPS=eps, FALSE=false, FILENAME=filename, $
  FINAL=final, GREY=grey, GRID=grid, HISTO=histo, HOLD=hold, HP=hp, INIT=init, $
  INSET_LL=inset_ll, INSET_LR=inset_lr, INSET_UL=inset_ul, INSET_UR=inset_ur, $
  INV=inv, JPEG=jpeg, LEVELS=levels, MIR=mir, MSWIN=mswin, MULTI=multi, NEW=new, NLEVELS=nlevels, $
  NOAXES=noaxes, NODATA=nodata, NOERASE=noerase,  NO_TICK=no_tick, NOTICKS=noticks, NOTICKV=notickv, $
  NO_XTICK=no_xtick, NOXTICKS=noxticks, NOXTICKV=noxtickv, NO_YTICK=no_ytick, $
  NOYTICKS=noyticks, NOYTICKV=noytickv, NWINDOW=nwindow, PS=ps, PROTECT=protect, SCATTER=scatter, $
  SUBTITLE=subtitle, SYMBOL=symbol, SYMSIZE=symsize, THICK=thick, XTHICK=xthick, YTHICK=ythick, TICKSIZE=ticksize, $
  TIFF=tiff, TITLE=title, USD=usd, WINDOW=window, XTITLE=xtitles, YTITLE=ytitles, $
  XLOG=xlog, XRANGE=xrange, YLOG=ylog, YRANGE=yrange, _EXTRA=extra

ON_ERROR,2                            ;Return to caller if an error is detected
; Use a common block to pass on information from one PLOTXY session to the next
COMMON PLOTXY, first, pkeep, xkeep, ykeep, previous, olddev, bitmap, mapname, xynorm, many, defcol

; For PCs, this is necessary to have the right colour allocation:
IF KEYWORD_SET(final) THEN GOTO, final
IF N_ELEMENTS(first) LE 0 THEN first=0
IF first EQ 0 THEN BEGIN
  pkeep=!P & xkeep=!X & ykeep=!Y
  olddev=!D.NAME
  IF KEYWORD_SET(ps)    THEN BEGIN & SET_PLOT,'PS'  & bitmap=0 & ENDIF
  IF KEYWORD_SET(eps)   THEN BEGIN & SET_PLOT,'PS'  & bitmap=0 & ENDIF
  IF KEYWORD_SET(hp)    THEN BEGIN & SET_PLOT,'HP'  & bitmap=0 & ENDIF
  IF KEYWORD_SET(mswin) THEN BEGIN & SET_PLOT,'WIN' & bitmap=-1 & ENDIF
  IF KEYWORD_SET(xwin)  THEN BEGIN & SET_PLOT,'X'   & bitmap=-1 & ENDIF
  DEVICE, DECOMPOSE=0
ENDIF
first=first+1
IF KEYWORD_SET(false) THEN LOADCT, 5, /SILENT
IF KEYWORD_SET(bw)    THEN LOADCT, 0, /SILENT

; Extract basic plot window information
IF NOT KEYWORD_SET(window) THEN BEGIN
  ; Specify a default window for single PS plots
  IF KEYWORD_SET(init) THEN window=[0., 0., 13., 10.]
  IF KEYWORD_SET(new) THEN window=[2.5, 1.5, 12.75, 9.25]
ENDIF
IF KEYWORD_SET(window) THEN BEGIN
  xpos=window(0) & ypos=window(1)
  xsize=ABS(window(2)-window(0)) & ysize=ABS(window(3)-window(1))
ENDIF

IF KEYWORD_SET(init) THEN GOTO, init

;-------------------------------------------------------------------------------
; Initialization of the data plotting

; Determine the size of the data arrays, and copy input arrays into workspace
sx=SIZE(xin)

ndim=sx(0)
i=FIX(sx(0))+1
type=sx(i)
npixx=FLOAT(sx(i+1))

CASE ndim OF
1: BEGIN                             ; One-dimensional data
     IF N_PARAMS(0) LT 2 THEN BEGIN
       x=INDGEN(npixx)
       y=xin
     ENDIF ELSE BEGIN
       sy=SIZE(yin)
       i=FIX(sy(0))+1
       type=sy(i)
       npixy=FLOAT(sy(i+1))

       IF npixx NE npixy THEN $
         MESSAGE,'WARNING, x and y are not equally large: ' + $
         'Nx = '+STRTRIM(LONG(npixx),2)+', Ny = '+STRTRIM(LONG(npixy),2)+ $
         '.  Using the smallest of both';,/INFORM
       npix=npixx < npixy
       x=xin
       y=yin
     ENDELSE
   END
2: BEGIN
     IF KEYWORD_SET(inv) THEN im=MAX(xin)+MIN(xin)-xin ELSE im=xin
     IF KEYWORD_SET(mir) THEN im=ROTATE(im,5)
     IF KEYWORD_SET(usd) THEN im=ROTATE(im,2)
     xdim=sx(1) & ydim=sx(2)
   END
ELSE: IF NOT KEYWORD_SET(nodata) THEN $
      MESSAGE,'Input data has improper dimensions'
ENDCASE

IF KEYWORD_SET(add) THEN BEGIN
  GOTO, add
ENDIF
IF KEYWORD_SET(new) THEN BEGIN
  previous=INTARR(4)
  GOTO, new
ENDIF

;-------------------------------------------------------------------------------
; Initialization of the plot page


; Branch point for initialization of plotpage
init:

; Set new device type
bitmap=-2 &  mapname=''

IF KEYWORD_SET(ps)    THEN BEGIN & SET_PLOT,'PS'  & bitmap=0 & ENDIF
IF KEYWORD_SET(eps)   THEN BEGIN & SET_PLOT,'PS'  & bitmap=0 & ENDIF
IF KEYWORD_SET(hp)    THEN BEGIN & SET_PLOT,'HP'  & bitmap=0 & ENDIF
IF KEYWORD_SET(mswin) THEN BEGIN & SET_PLOT,'WIN' & bitmap=-1 & ENDIF
IF KEYWORD_SET(xwin)  THEN BEGIN & SET_PLOT,'X'   & bitmap=-1 & ENDIF
; Default is on screen
IF bitmap LT -1 THEN BEGIN
    IF !VERSION.OS_FAMILY EQ 'unix' THEN SET_PLOT,'X' ELSE SET_PLOT,'WIN'
    bitmap=-1
ENDIF

IF KEYWORD_SET(tiff)  THEN bitmap=1
IF KEYWORD_SET(bmp)   THEN bitmap=2
IF KEYWORD_SET(jpeg)  THEN BEGIN
  bitmap=3
  IF jpeg GT 1 THEN quality=jpeg ELSE quality=75
ENDIF

; Set default color
IF bitmap NE 0 THEN defcol=255 ELSE defcol=0

; Define default plotpage
IF NOT KEYWORD_SET(window) THEN BEGIN
  IF bitmap EQ 0 THEN BEGIN
    xpos=0 & ypos=0 & xsize=21 & ysize=30
  ENDIF ELSE BEGIN
    xpos=0 & ypos=0
    IF ndim GT 1 THEN BEGIN
      xsize=xdim & ysize=ydim & noaxes=1
    ENDIF ELSE BEGIN
      xsize=550 & ysize=500
    ENDELSE
  ENDELSE
ENDIF

; Define default filename
IF N_ELEMENTS(filename) LE 0 THEN filename=''

; Set default title
IF N_ELEMENTS(subtitle) LE 0 THEN subtitle=''

IF bitmap EQ 0 THEN BEGIN
  IF KEYWORD_SET(ps) AND STRPOS(filename,'.ps') LT 0 THEN filename='plotxy.ps'
  IF KEYWORD_SET(hp) AND STRPOS(filename,'.hp') LT 0 THEN filename='plotxy.hp'
  IF KEYWORD_SET(eps) AND STRPOS(filename,'.eps') LT 0 THEN filename='plotxy.eps'

  DEVICE, FILENAME=filename
  IF KEYWORD_SET(eps) THEN DEVICE, /ENCAPSULATED, /COLOR, bits_per_pixel=8, /TIMES
  DEVICE, XSIZE=xsize, YSIZE=ysize, XOFFSET=xpos, YOFFSET=ypos
  ; Set the plotpage for positioning of the title
  IF N_ELEMENTS(title) LE 0 OR N_ELEMENTS(init) LE 0 THEN htitle=' ' ELSE BEGIN
    htitle=title
    IF NOT KEYWORD_SET(charsize) THEN charsize=2.0
    IF NOT KEYWORD_SET(charthick) THEN charthick=1.5
    !P.POSITION=[0.,0.,1.,1.25]
    PLOT, [0.,1.], [0.,1.], /NODATA, /NOERASE, TITLE=htitle, SUBTITLE=subtitle,$
          XSTYLE=4, YSTYLE=4, CHARSIZE=charsize, CHARTHICK=charthick
  ENDELSE
ENDIF ELSE BEGIN

  ; Plot diagram in new window with landscape A4 scale. 
  hold_window = 0
  IF NOT KEYWORD_SET(nwindow) THEN nwindow=!D.WINDOW+1 ELSE IF nwindow EQ !D.WINDOW THEN hold_window = 1
  IF N_ELEMENTS(title) LE 0 THEN title='PLOTXY '+STRTRIM(nwindow,2)
  
  ; If window number 'nwindow' is already open, overplot in the same window
  IF hold_window NE 1 THEN BEGIN
    ; Use RETAIN=2 to avoid problems with TVRD. See IDL user manual.
    IF xpos GE 0 AND ypos GE 0 THEN BEGIN
      WINDOW, nwindow, XPOS=xpos, YPOS=ypos, XSIZE=xsize, YSIZE=ysize, TITLE=title, RETAIN=2
    ENDIF ELSE BEGIN
      WINDOW, nwindow, XSIZE=xsize, YSIZE=ysize, TITLE=title, RETAIN=2
    ENDELSE
  ENDIF

  ; For inverted image (black-on-white) use trick to create window with white background
  IF KEYWORD_SET(inv) THEN BEGIN
    d=BYTARR(xsize,ysize)+255
    TV, d
    d=0
  ENDIF

  ; Bitmaps are created using TVRD, so the filename must be kept to the last.
  IF KEYWORD_SET(tiff) AND STRPOS(filename,'.tif') LT 0 THEN filename='plotxy.tif'
  IF KEYWORD_SET(jpeg) AND STRPOS(filename,'.jpg') LT 0 THEN filename='plotxy.jpg'
  IF KEYWORD_SET(bmp) AND STRPOS(filename,'.bmp') LT 0 THEN filename='plotxy.bmp'
  mapname=filename
ENDELSE

; Reset values for page initialization
previous=INTARR(4)
!P.POSITION=[0,0,0,0]

IF NOT KEYWORD_SET(many) THEN BEGIN
  CASE N_ELEMENTS(multi) OF
  0:    BEGIN & !P.MULTI=0 &     many=0 & END
  1:    BEGIN & !P.MULTI=0 &     many=0 & END
  ELSE: BEGIN & !P.MULTI=multi & many=1 & END
  ENDCASE
ENDIF

IF KEYWORD_SET(init) THEN RETURN

;-------------------------------------------------------------------------------
; Plot the data

; Branch point for a new plot on the same page
new:

; Take care of this stupid factor 1000 necessary in ps plot positioning.
IF KEYWORD_SET(new) AND KEYWORD_SET(window) AND NOT KEYWORD_SET(many) THEN $
  IF bitmap EQ 0 THEN !P.POSITION=window*1000. ELSE !P.POSITION=window

IF N_ELEMENTS(title) LE 0 THEN title=' '       & !P.TITLE=title
IF N_ELEMENTS(subtitle) LE 0 THEN subtitle=' ' & !P.SUBTITLE=subtitle
axtitle=' ' & aytitle=' '
CASE N_ELEMENTS(xtitles) OF
  0: !X.TITLE='X'
  1: !X.TITLE=xtitles
  ELSE: BEGIN & !X.TITLE=xtitles(0) & axtitle=xtitles(1) & END
ENDCASE
CASE N_ELEMENTS(ytitles) OF
  0: !Y.TITLE='Y'
  1: !Y.TITLE=ytitles
  ELSE: BEGIN & !Y.TITLE=ytitles(0) & aytitle=ytitles(1) & END
ENDCASE

axes=0
IF KEYWORD_SET(xlog) THEN BEGIN
  axes=axes+1
  IF KEYWORD_SET(xrange) THEN xmin=xrange(0) ELSE BEGIN
    w0=WHERE(x GT 0)
    xmin=MIN(x(w0))
  ENDELSE
  x = x > 0.1*xmin
ENDIF
IF KEYWORD_SET(ylog) THEN BEGIN
  axes=axes+2
  IF KEYWORD_SET(yrange) THEN ymin=yrange(0) ELSE BEGIN
    w0=WHERE(y GT 0)
    ymin=MIN(y(w0))
  ENDELSE
  y = y > 0.1*ymin
ENDIF

; Determine the begin and end points of the axes

CASE N_ELEMENTS(xrange) OF
4: BEGIN
     !X.STYLE=9
     xmin=xrange(0) & xmax=xrange(1) & dx=ABS(xmax-xmin)
     axrange=[xrange(2),xrange(3)]
   END
2: BEGIN
     !X.STYLE=1
     xmin=xrange(0) & xmax=xrange(1) & dx=ABS(xmax-xmin)
   END
ELSE: BEGIN
     !X.STYLE=1
     IF ndim EQ 1 THEN BEGIN
       xmin=MIN(x) & xmax=MAX(x)
     ENDIF ELSE BEGIN
       xmin=0 & xmax=xdim-1
     ENDELSE
     dx=ABS(xmax-xmin)
     IF KEYWORD_SET(xlog) THEN xrange=[xmin, xmin+1.05*dx] $
     ELSE IF ndim EQ 1 THEN xrange=[xmax-1.05*dx,xmin+1.05*dx] $
                       ELSE xrange=[xmin, xmax]
   END
ENDCASE
!X.RANGE=xrange(0:1)

CASE N_ELEMENTS(yrange) OF
4: BEGIN
     !Y.STYLE=9
     ymin=yrange(0) & ymax=yrange(1) & dy=ABS(ymax-ymin)
     ayrange=[yrange(2),yrange(3)]
   END
2: BEGIN
     !Y.STYLE=1
     ymin=yrange(0) & ymax=yrange(1) & dy=ABS(ymax-ymin)
   END
ELSE: BEGIN
     !Y.STYLE=1
     IF ndim EQ 1 THEN BEGIN
       ymin=MIN(y) & ymax=MAX(y)
     ENDIF ELSE BEGIN
       ymin=0 & ymax=ydim-1
     ENDELSE
     dy=ABS(ymax-ymin)
     IF KEYWORD_SET(ylog) THEN yrange=[ymin,ymin+1.05*dy] $
     ELSE IF ndim EQ 1 THEN yrange=[ymax-1.05*dy,ymin+1.05*dy] $
                       ELSE yrange=[ymin,ymax]
   END
ENDCASE
!Y.RANGE=yrange(0:1)

; Branch point for additional overplotting
add:

; Reset values for plot initialization
IF KEYWORD_SET(noerase) THEN !P.NOERASE=1

; Set the symbol or linestyle!
!P.PSYM=0
!P.LINESTYLE=0
IF KEYWORD_SET(symsize) THEN !P.SYMSIZE=symsize ELSE !P.SYMSIZE=1
IF KEYWORD_SET(thick) THEN !P.THICK=thick ELSE !P.THICK=1
IF KEYWORD_SET(xthick) THEN !X.THICK=thick ELSE !X.THICK=1
IF KEYWORD_SET(ythick) THEN !Y.THICK=thick ELSE !Y.THICK=1
IF KEYWORD_SET(symbol) THEN SETSYM, symbol, THICK=thick
IF KEYWORD_SET(scatter) THEN !P.PSYM=3
IF KEYWORD_SET(histo) THEN !P.PSYM=10
IF KEYWORD_SET(charsize) THEN !P.CHARSIZE=charsize ELSE !P.CHARSIZE=1
IF KEYWORD_SET(charthick) THEN !P.CHARTHICK=charthick ELSE !P.CHARTHICK=1
IF KEYWORD_SET(color) THEN !P.COLOR=color ELSE !P.COLOR=defcol
IF KEYWORD_SET(noaxes) THEN BEGIN & !X.STYLE=12 & !Y.STYLE=12 & ENDIF
IF KEYWORD_SET(no_tick) THEN BEGIN & noticks=1 & notickv=1 & ENDIF
IF KEYWORD_SET(no_xtick) THEN BEGIN & noxticks=1 & noxtickv=1 & ENDIF
IF KEYWORD_SET(no_ytick) THEN BEGIN & noyticks=1 & noytickv=1 & ENDIF
IF KEYWORD_SET(notickv) THEN BEGIN & noxtickv=1 & noytickv=1 & ENDIF
IF KEYWORD_SET(noxtickv) THEN !X.TICKNAME=REPLICATE(' ',30)
IF KEYWORD_SET(noytickv) THEN !Y.TICKNAME=REPLICATE(' ',30)
IF KEYWORD_SET(noticks) THEN BEGIN & noxticks=1 & noyticks=1 & ENDIF
IF KEYWORD_SET(noxticks) THEN !X.TICKS=1
IF KEYWORD_SET(noyticks) THEN !Y.TICKS=1

axticklen=0 & ayticklen=0
IF N_ELEMENTS(ticksize) GT 0 THEN BEGIN
  !X.TICKLEN=ticksize(0)
  axticklen =ticksize(0)
ENDIF
IF N_ELEMENTS(ticksize) GT 1 THEN BEGIN
  !Y.TICKLEN=ticksize(1)
  ayticklen =ticksize(1)
ENDIF
IF N_ELEMENTS(ticksize) GT 2 THEN axticklen =ticksize(2)
IF N_ELEMENTS(ticksize) GT 3 THEN ayticklen =ticksize(3)

IF KEYWORD_SET(grid) THEN BEGIN
  ticklen=!P.TICKLEN
  !P.TICKLEN=1
  IF grid NE 1 THEN BEGIN
    !X.GRIDSTYLE=ABS(grid) & !Y.GRIDSTYLE=ABS(grid)
  ENDIF
ENDIF

; Do the actual plotting
IF ndim EQ 1 THEN BEGIN
  IF KEYWORD_SET(add) THEN BEGIN
    IF NOT KEYWORD_SET(nodata) THEN OPLOT, x, y, _EXTRA=extra
  ENDIF ELSE BEGIN
    CASE axes OF
      0: PLOT, x, y, /DEVICE, _EXTRA=extra, NODATA=nodata
      1: PLOT, x, y, /XLOG, /DEVICE, _EXTRA=extra, NODATA=nodata
      2: PLOT, x, y, /YLOG, /DEVICE, _EXTRA=extra, NODATA=nodata
      3: PLOT, x, y, /XLOG, /YLOG, /DEVICE, _EXTRA=extra, NODATA=nodata
    ENDCASE
  ENDELSE
ENDIF ELSE BEGIN
  IF KEYWORD_SET(contour) THEN BEGIN
    nc=!D.N_COLORS-1                       ; brightest color
    c_color=[0,0,0,0,0,0,0,0,0]            ; color vector
    maxim=MAX(im,MIN=minim) & irange=maxim-minim
    IF NOT KEYWORD_SET(nlevels) THEN nlevels=4 & nlevels=nlevels+1
    IF NOT KEYWORD_SET(levels) THEN levels=minim+FINDGEN(nlevels)/nlevels*irange
    x=xrange(0)+FINDGEN(xdim)*ABS(xrange(1)-xrange(0))/FLOAT(xdim-1)
    y=yrange(0)+FINDGEN(ydim)*ABS(yrange(1)-yrange(0))/FLOAT(ydim-1)
    CONTOUR,im,x,y,LEVELS=levels, /FOLLOW, /DEVICE, /NOERASE, NODATA=nodata, $
     _EXTRA=extra
  ENDIF ELSE BEGIN
    IF bitmap EQ 0 THEN BEGIN
      TVSCL,im,xpos,ypos,XSIZE=xsize,YSIZE=ysize,/CENTIMETERS
      IF many EQ 0 THEN !P.POSITION=[xpos,ypos,xpos+xsize,ypos+ysize]*1000.
    ENDIF ELSE BEGIN
      IF xdim NE xsize OR ydim NE ysize THEN $
        temp=CONGRID(im,xsize,ysize,/CUBIC) ELSE temp=im
      IF NOT KEYWORD_SET(nodata) THEN TVSCL,temp,xpos,ypos
      IF many EQ 0 THEN !P.POSITION=[xpos,ypos,xpos+xsize,ypos+ysize]
    ENDELSE
    IF NOT KEYWORD_SET(noaxes) THEN $
      PLOT, [0.,1.], [0.,1.], /NODATA, /DEVICE, /NOERASE, _EXTRA=extra
  ENDELSE
ENDELSE


; Add the additional axes if two axes were specified
!X.STYLE=1 & !Y.STYLE=1
!X.TICKLEN=axticklen & !Y.TICKLEN=ayticklen
!X.TICKS=0 & !Y.TICKS=0
!X.TICKNAME='' & !Y.TICKNAME=''
IF N_ELEMENTS(axrange) GT 0 THEN $
  CASE axes OF
    0: AXIS, XAXIS=1, XTITLE=axtitle, XRANGE=axrange
    1: AXIS, XAXIS=1, XTITLE=axtitle, XRANGE=axrange, /XLOG
    2: AXIS, XAXIS=1, XTITLE=axtitle, XRANGE=axrange
    3: AXIS, XAXIS=1, XTITLE=axtitle, XRANGE=axrange, /XLOG
  ENDCASE

IF N_ELEMENTS(ayrange) GT 0 THEN $
  CASE axes OF
    0: AXIS, YAXIS=1, YTITLE=aytitle, YRANGE=ayrange
    1: AXIS, YAXIS=1, YTITLE=aytitle, YRANGE=ayrange
    2: AXIS, YAXIS=1, YTITLE=aytitle, YRANGE=ayrange, /YLOG
    3: AXIS, YAXIS=1, YTITLE=aytitle, YRANGE=ayrange, /YLOG
  ENDCASE

; Add the insets. Note that stacked sets of insets can be made through added
; calls to PLOTXY
; Establish the normal coordinate system for later use
IF NOT KEYWORD_SET(add) AND NOT KEYWORD_SET(protect) THEN xynorm=CONVERT_COORD(xrange,yrange,/DATA,/TO_DEVICE)
IF KEYWORD_SET(inset_ll) THEN PLOT_INSET, inset_ll, xynorm, [+1,+1], previous
IF KEYWORD_SET(inset_ul) THEN PLOT_INSET, inset_ul, xynorm, [+1,-1], previous
IF KEYWORD_SET(inset_lr) THEN PLOT_INSET, inset_lr, xynorm, [-1,+1], previous
IF KEYWORD_SET(inset_ur) THEN PLOT_INSET, inset_ur, xynorm, [-1,-1], previous

; Reset some of the plot parameters
!P.PSYM=0
!P.LINESTYLE=0
!P.THICK=1
!X.THICK=1
!Y.THICK=1
!P.TICKLEN=0.02
!P.CHARSIZE=1
!P.CHARTHICK=1
!P.SYMSIZE=1
!X.TICKLEN=0 & !Y.TICKLEN=0
!X.GRIDSTYLE=0 & !Y.GRIDSTYLE=0
!P.COLOR=defcol
!P.NOERASE=0
!P.TITLE=''

; Leave the routine
IF KEYWORD_SET(new) OR KEYWORD_SET(hold) OR KEYWORD_SET(add) THEN RETURN

;------------------------------------------------------------------------------
; Return the system in its original state

; Branch point for closing the device, and finishing the ps file
final:

IF bitmap GT 0 THEN BEGIN
  IF bitmap LT 2 THEN map=ROTATE(255-TVRD(),7) ELSE map=255-TVRD()
  IF KEYWORD_SET(usd) THEN map=REVERSE(map,2)
  IF KEYWORD_SET(inv) THEN map=255-map
ENDIF

CASE bitmap OF
  0: DEVICE, /CLOSE
  1: WRITE_TIFF, mapname, map
  2: WRITE_BMP,  mapname, map
  3: WRITE_JPEG, mapname, map, QUALITY=quality
  ELSE:
ENDCASE

; Restore internal plot parameters as they were before calling PLOTXY
SET_PLOT,olddev
!P=pkeep
!X=xkeep
!Y=ykeep

; Free memory
first=0 & pkeep=0 & xkeep=0 & ykeep=0 & previous=0 & keep=0 & bitmap='' & mapname='' & xynorm=0 & many=0 & defcol=0
x=0 & y=0 & im=0

END

;------------------------------------------------------------------------------------------------------------------

PRO PLOT_TEST

;PLOTXY, RANDOMU(seed,1000), RANDOMU(seed,1000), /SCATTER

;PLOTXY, /INIT, /TIFF, WIN=[0,0,500,500], FILE='test.tif'
;PLOTXY, /NEW, RANDOMU(seed,1000), RANDOMU(seed,1000), /SCATTER
;PLOTXY, /FIN

PLOTXY, /INIT, WIN=[0,0,15.,15.], MULTI=[0, 2, 2], /PS, FILE='test.ps'
PRINT, !P.MULTI

PLOTXY, /NEW, RANDOMU(seed,1000), RANDOMU(seed,1000), /SCATTER
PLOTXY, /NEW, 0.001*DINDGEN(1000), 0.001*DINDGEN(1000)
PLOTXY, /NEW, 0.001*DINDGEN(1000), 1D0-0.001*DINDGEN(1000)
PLOTXY, /NEW, RANDOMU(seed,1000), RANDOMU(seed,1000), /SCATTER

PLOTXY, /FIN
PRINT, !P.MULTI

END