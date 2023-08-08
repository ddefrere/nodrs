PRO NOMIC_CHANNEL_STAT, img_in

; Recover the IDL running path
IF !VERSION.OS_FAMILY EQ 'unix' THEN sep='/' ELSE sep='\'
path        = GET_PATH('nomic_null.pro', N_DIR_UP=1)           ; Main path of LBTI software
dark_path   = path + 'results' + sep + 'detector' + sep        ; Path to plot the results
IF NOT FILE_TEST(dark_path) THEN FILE_MKDIR, bckg_path         ; Create directory if it does not exist

; Nomic infos
col = [0, 512, 1024]                               & n_col = N_ELEMENTS(col) - 1
row = [0, 128, 256, 384, 512, 640, 768, 896, 1024] & n_row = N_ELEMENTS(row) - 1

; Inititate plot      
LOADCT, 0, /SILENT  
fit      = 20./1720.
inv      = 0.
xrange   = [MIN(img_in),MAX(img_in)]
PLOTXY, /INIT, INV=inv, /COLOR, /EPS, WINDOW=[0, 0, 1200, 800]*fit, FILENAME = dark_path + 'nomic_behavior.eps'
PLOTXY, [0,0], [0,0], /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='NOMIC detector display', GRID=0, $
        XSTYLE=1, YSTYLE=1, /NODATA, /NOERASE, WINDOW=[150,100,1150,700]*fit;, INSET_UR='a.'
FOR i_col = 0, n_col-1 DO BEGIN
  FOR i_row = 0, n_row-1 DO BEGIN
    col_span = col[i_col+1] - col[i_col]
    row_span = row[i_row+1] - row[i_row]
    img_tmp = EXTRAC(img_in, col[i_col], row[i_row], col_span, row_span) 
    AUTOHIST, REFORM(img_tmp, N_ELEMENTS(img_tmp)), x_lin, y_lin , x_cen , y_cen, /NOPLOT
    ; Now plot the histogram
    xpos0=0.1+0.8*(i_col/n_col) & xpos1=0.1+0.8*((1+i_col)/n_col)
    ypos0=0.1+0.8*(i_row/n_row) & ypos1=0.1+0.8*((1+i_row)/n_row) 
    LOADCT, 0, /SILENT
    PLOTXY,  x_cen, y_cen, /ADD, /HISTO, /NOERASE, COLOR=240, YRANGE=[0.,1.2*MAX(y_cen)], XRANGE=xrange, XTICKS=5, XTICKNAME=REPLICATE(' ', 6), WINDOW=[xpos0,ypos0,xpos1,ypos1]*fit*1150
    LOADCT, 27, /SILENT
    OPLOT, x_cen, y_cen, PSYM=10, COLOR=240
  ENDFOR
ENDFOR 
PLOTXY, /FIN
LOADCT, 0, /SILENT 
 
END