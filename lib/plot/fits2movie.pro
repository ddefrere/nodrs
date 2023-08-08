PRO FITS2MOVIE, FILES=files, DIRECTORY=directory, XC=xc, YC=yc, N_CLIP=n_clip, SCALE_X=scale_x, SCALE_Y=scale_y, N_FR=n_fr

; PURPOSE:
;   This routine produces an MPEG movie from fits files. 
;
; MODIFICATION HISTORY:
;   Version 1.0, 01-NOV-2012, by Denis Defr�re, University of Arizona, ddefrere@email.arizona.edu

; Read fits files
IF KEYWORD_SET(directory) AND NOT KEYWORD_SET(files) THEN files = FILE_SEARCH(directory,'*.fits')
;files = files[200*20+indgen(600)]
n_files = N_ELEMENTS(files)

; Recover the IDL runnig path
IF !VERSION.OS_FAMILY EQ 'unix' THEN sep='/' ELSE sep='\'
path         = GET_PATH('fits2movie.pro', N_DIR_UP=0)                 ; Main path of LBTI software
movie_path   = path + 'movie' + sep                           ; Path to the result folder
IF NOT FILE_TEST(movie_path) THEN FILE_MKDIR, movie_path       ; Create directory if it does not exist

; Compute the time ratio with mpeg movies
fr_rate  = 60                           ; [fps]
mpg_rate = 60                           ; [fps], standard mpeg fps
ratio    = FLOOR(mpg_rate/fr_rate)
          
; Read the files
LOADCT, 3, /SILENT
temp_path = movie_path + sep + 'temp' + sep      ; temporary folder to store the images for the movie
FILE_MKDIR, temp_path                            ; create directory
i_c = 0                                          ; initialize file counter
FOR i_f=0, n_files-1 DO BEGIN
    ; Extract fits image
    img0 = READFITS(files[i_f], header, /SILENT)
    
    ; Find the number of frames
    n_fr = N_ELEMENTS(img0[0,0,*])
    
    FOR i_fr = 0, n_fr-1 DO BEGIN
      img_tmp = img0[*,*,i_fr]
      IF KEYWORD_SET(N_CLIP) THEN img_tmp = EXTRAC(img_tmp, xc-0.5* scale_x*n_clip, yc-0.5*scale_y*n_clip, scale_x*n_clip, scale_y*n_clip)
      
      ; Remove bad pixels
      ;box_width = 20 
      ;Nsigma    = 10
      ;img_tmp   = SIGMA_FILTER(img_tmp, box_width, N_SIGMA=Nsigma, MONITOR=monitor, RADIUS=radius, $
      ;                         N_CHANGE=nchange, VARIANCE_IMAGE=imvar, DEVIATION_IMAGE=imdev, /KEEP)
       
      ; Congrid image for display purposes
      IF KEYWORD_SET(SCALE_X) AND KEYWORD_SET(SCALE_Y) THEN img_tmp = CONGRID(img_tmp, scale_x*256, scale_y*256) 
                                                                
      ; Plot image and save it          
      PLOTXY, img_tmp;, WINDOW=[0, 0, scale_x*256, scale_y*256]
      FOR i_d=0, ratio-1 DO BEGIN 
        img_name = 'temp_' + STRING(i_c, FORMAT='(I6.6)')+ '.jpg'         
        SAVEIMAGE, temp_path + img_name, /JPEG, QUALITY=QUALITY, DITHER=DITHER, CUBE=CUBE, /QUIET    
        i_c = i_c + 1
      ENDFOR
      WDELETE   
    ENDFOR  
ENDFOR
PRINT, i_c

; Define the number of frames per movie
IF NOT KEYWORD_SET(n_fr) THEN n_fr = n_files
n_seq = 1;FLOOR(n_files/n_fr)

; Now produce the movie
n_fr = n_fr*ratio
FOR i_s = 0, n_seq-1 DO make_mpg, prefix = temp_path + 'temp_', suffix='jpg', n_start=i_s*n_fr, n_end=(i_s+1)*(n_fr-1), frame_rate=5, digits=6, mpeg_file = movie_path + sep + 'fits_movies_seq-' + STRING(i_s, FORMAT='(I0)') + '.mpg'
FILE_DELETE, temp_path, /RECURSIVE                  ; delete temporary directory
END