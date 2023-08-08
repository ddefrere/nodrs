;+
; NAME: PLOT_FFT
; 
; PURPOSE:
;   Standard routine to plot an image, its FFT and its PSF
;
; INPUTS:
;   img_in        : The input image
;
; KEYWORDS
;   EPS           : Saved image as an EPS file
;
; MODIFICATION HISTORY:
;   Version 1.0,  01-MAR-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu

PRO PLOT_FFT, img_in, EPS=eps

; Recover the data path
IF !VERSION.OS_FAMILY EQ 'unix' THEN sep='/' ELSE sep='\'
path          = GET_PATH('nomic_null.pro', N_DIR_UP=1)           ; Main path of LBTI software
result_path   = path + 'results' + sep + 'sensitivity'           ; Path to the result folder
IF NOT FILE_TEST(result_path) THEN FILE_MKDIR, result_path       ; Create directory if it does not exist

; Transform the image into the frequency domain and shift the zero-frequency location from (0,0) to
; the center of the data.
img_fft = FFT(img_in, /CENTER)
 
; Compute the power spectrum of the transform and apply a log scale.
psf_img  = ABS(img_fft)^2
psf_log  = ALOG10(psf_img)
 
; Scale the power spectrum to make its maximum value equal to 0.
psf_sca = psf_log - MAX(psf_log)
  
; Apply a mask to remove values less than -7, just below the peak of the power spectrum. The data type
; of the array returned by the FFT function is complex, which contains real and imaginary parts. In image 
; processing, we are more concerned with the amplitude, which is the real part. The amplitude is the only part
; represented in the surface and displays the results of the transformation.
mask = REAL_PART(psf_sca) GT -7
 
; Apply the mask to the transform to exclude the noise.
psf_new = psf_sca*mask
 
; Perform an inverse FFT to the masked transform, to transform it back to the spatial domain.
img_new = REAL_PART(FFT(psf_new, /INVERSE))
 
; Plot the input image, the FFT, and the cleaned image
LOADCT, 0, /SILENT
IF KEYWORD_SET(EPS) THEN fit = 20./1720. ELSE fit = 1.
IF KEYWORD_SET(EPS) THEN PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1900, 600]*fit, FILENAME=plot_name $
   ELSE PLOTXY, /INIT, WINDOW=[0, 0, 1900, 600]*fit
PLOTXY, img_in, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='Detector frame', XTITLE='', YTITLE='', GRID=0, $
        CHARSIZE=1.0, CHARTHICK=1.5, THICK=3.0, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=[100,50,600,550]*fit;, INSET_UR='a.'
PLOTXY, psf_log, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='Power spectrum', XTITLE='', YTITLE='', GRID=0, $
        CHARSIZE=1.0, CHARTHICK=1.5, THICK=3.0, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=[715,50,1215,550]*fit;, INSET_UR='a.'
PLOTXY, img_new, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='Cleaned image', XTITLE='', YTITLE='', GRID=0, $
        CHARSIZE=1.0, CHARTHICK=1.5, THICK=3.0, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=[1330,50,1830,550]*fit;, INSET_UR='a.'        
; Close plot
LOADCT, 0, /SILENT
PLOTXY, /FIN
END