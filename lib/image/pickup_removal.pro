;+
; NAME: PICKUP_REMOVAL
; 
; PURPOSE:
;   This function removes high frequency noise from the input image and returns the corresponding filtered image.
;   If an image cube is provided, the FFT is performed along the time direction rather than spatially.
;   This is particularly useful to estimate and remove pickup noise. Standard plot of the image, its FFT and its PSF are also displayed.
;
; INPUTS:
;  The input image or image cube.
;
; KEYWORDS
;   FREQ_CUT      : The frequency cut (either in Hz or in 1/pix)
;   PIXSIZE       : The pixel size (in mas)
;   PERIOD        : The time between each image in the data cube (in sec)
;   NO_GPU        : Set this keyword to prevent the use of GPU-assisted computation (only available for NVIDIA GPUs).
;   PLOT          : Set this keyword to plot the input and output detector frames
;
; MODIFICATION HISTORY:
;   Version 1.0,  01-MAR-2013, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu
;   Version 1.0,  13-JUN-2013, DD: implemented FFT computing by GPU
;   Version 1.1,  13-JUN-2013, DD: added keywords INIT and ADD

FUNCTION PICKUP_REMOVAL, img_in, FREQ_CUT=freq_cut, PIXSIZE=pixsize, PERIOD=period, NO_GPU=no_gpu, PLOT=plot

; Check whether GPU-assisted computation can be enabled on this machine and initialize it
IF NOT KEYWORD_SET(NO_GPU) THEN gpu = gpu_detect() ELSE gpu = 0 ; return 1 if yes

; Compute the number of frames in the input data
n_fr = N_ELEMENTS(img_in[0,0,*])

; Do the spatial or temporal FFT
IF n_fr EQ 1 THEN BEGIN 
    ; Compute array center
    npix = N_ELEMENTS(img_in[*,0])
    xcen = 0.5*npix-0.5 & ycen = xcen
    
    ; Compute frequency mask. According to FFT help, frequency is sampled up to the critical Nyquist
    ; frequency (i.e., 1/(2*T)). Here, T is one pixel.
    DIST_CIRCLE, dist_map, npix, xcen, ycen, /DOUBLE
    mask    = dist_map LE freq_cut
    idx_out = WHERE(dist_map GT freq_cut, n_out)
    
    ; Start FFT computation (with or without GPU)
    IF gpu NE 1 THEN BEGIN
      ; Transform the image into the frequency domain and shift the zero-frequency location from (0,0) 
      ; the center of the data.
      img_fft = FFT(img_in, /CENTER)
      
      ; Compute the power spectrum of the transform and apply a log scale.
      psf_img  = ABS(img_fft)^2
      psf_log  = ALOG10(psf_img)
      
      ; Scale the power spectrum to make its maximum value equal to 0.
      psf_sca = psf_log - MAX(psf_log)
    
      ; Apply the mask to the transform to exclude the noise.
      img_new = img_fft*mask
      ;IF n_out GT 0 THEN img_fft[idx_out] = MEAN(img_fft[idx_out])
      ;img_new = img_fft
      
      ; Perform an inverse FFT to the masked transform, to transform it back to the spatial domain.
      img_out = REAL_PART(FFT(img_new, /INVERSE, /CENTER))
    ENDIF ELSE BEGIN
      ; Transfer data to GPU
      gpuPutArr, img_in, img_gpu
      gpuPutArr, mask, mask_gpu
      
      ; Transform the image into the frequency domain and shift the zero-frequency location from (0,0) to
      ; the center of the data.
      img_tmp = gpufft(img_gpu, -1)
      img_fft = gpushift(img_tmp, -npix/2., -npix/2.)
      
      ; Apply the mask to the transform to exclude the noise.
      img_new = gpuMatrix_Multiply(img_fft,mask)
      
      ; Perform an inverse FFT to the masked transform, to transform it back to the spatial domain.
      img_tmp  = gpufft(img_new, 1)
      img_real = gpuabs(img_tmp)
      img_shft = gpushift(img_real, -npix_gpu/2., -npix_gpu/2.)
     
      ; transfer results back to CPU
      img_out = gpugetarr(img_tmp)
      psf_log = ALOG10(gpugetarr(img_fft))

      ; Free memory
      gpufree, [img_gpu, mask_gpu, img_fft, img_tmp, img_real, img_shft]
    ENDELSE
ENDIF ELSE BEGIN
    ; Actual FFT for each pixel
    fft_time = FFT(img_in, DIMENSION=3, /CENTER)
    ; Median FFT per pixel
    fft_time = MEDIAN(img_in, DIMENSION=3)
ENDELSE

; Plot the input image, the FFT, and the cleaned image
IF KEYWORD_SET(PLOT) THEN BEGIN         
  ; Display time or spatial FFT
  IF n_fr GT 1 THEN BEGIN
    ; Compute one-side PSD
    IF KEYWORD_SET(PERIOD) THEN x_coord = 1./period*INDGEN(0.5*n_fr)/(0.5*n_fr-1) ELSE x_coord = INDGEN(0.5*n_fr)
    psd     = REVERSE(ABS(fft_time[x_coord])^2) 
    ; Now plot
    LOADCT, 0, /SILENT
    xrange  = [MIN(x_coord[1:*]),MAX(x_coord)]
    yrange  = [0.5*MIN(psd),1.5*MAX(psd)]
    IF KEYWORD_SET(EPS) THEN fit = 20./1720. ELSE fit = 1.
    IF KEYWORD_SET(EPS) THEN PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 600, 600]*fit, FILENAME=plot_name $
       ELSE PLOTXY, /INIT, WINDOW=[0, 0, 800, 600]*fit
    ; Plot input image
    PLOTXY, x_coord, psd, /YLOG, /XLOG, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='Temporal PSD', XTITLE='Frequency [Hz]', YTITLE='Amplitude', GRID=0, $
            CHARSIZE=1.0, CHARTHICK=1.5, THICK=3.0, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=[100,50,750,550]*fit;, INSET_UR='a.'
    LOADCT, 0, /SILENT
    PLOTXY, /FIN
  ENDIF ELSE BEGIN
      ; Display spatial FFT
      LOADCT, 0, /SILENT
      fit = 0.7
      PLOTXY, /INIT, NWINDOW=1, WINDOW=[0, 0, 1900, 600]*fit
      ; Plot input image
      IF KEYWORD_SET(PIXSIZE) THEN fac = pixsize ELSE fac = 1.
      xrange = [-1.,1.]*0.5*npix*fac
      PLOTXY, img_in, /NEW, NWINDOW=1, XRANGE=xrange, YRANGE=xrange, TITLE='Detector frame', XTITLE='Angular separation [mas]', YTITLE='Angular separation [mas]', GRID=0, $
        CHARSIZE=1.2, CHARTHICK=1.5, THICK=3.0, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=[100,50,600,550]*fit;, INSET_UR='a.'
      ; Overplot circles for PSF and background computation
      LOADCT, 13, /SILENT
      LOADCT, 0, /SILENT
      ; Plot the PSD
      PLOTXY, psf_log, /NEW, NWINDOW=1, XRANGE=[-0.5,0.5]/fac, YRANGE=[-0.5,0.5]/fac, TITLE='Spatial PSD (log scale)', XTITLE='Spatial frequency [1/mas]', YTITLE='Spatial frequency [1/mas]', GRID=0, $
        CHARSIZE=1.2, CHARTHICK=1.5, THICK=3.0, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=[715,50,1215,550]*fit;, INSET_UR='a.'
      TVCIRCLE, freq_cut/(npix*fac), 0., 0., THICK=2., LINESTYLE=0, /DATA
      ; Plot the output image
      PLOTXY, img_out, /NEW, NWINDOW=1, XRANGE=xrange, YRANGE=xrange, TITLE='Cleaned frame', XTITLE='Angular separation [mas]', YTITLE='Angular separation [mas]', GRID=0, $
        CHARSIZE=1.2, CHARTHICK=1.5, THICK=3.0, XSTYLE=1, YSTYLE=1, /NOERASE, WINDOW=[1330,50,1830,550]*fit;, INSET_UR='a.'
      ; Overplot circles for PSF and background computation
      LOADCT, 13, /SILENT
      ; Close plot
      LOADCT, 0, /SILENT
      PLOTXY, /FIN    
  ENDELSE
ENDIF

RETURN, img_out
END

;------------------------------------------
PRO TEST_PICKUP_REMOVAL

; Recover the IDL running path
IF !VERSION.OS_FAMILY EQ 'unix' THEN sep='/' ELSE sep='\'
path          = GET_PATH('test_pickup_removal.pro', N_DIR_UP=1)           ; Main path of LBTI software
l1fits_path   = path + 'results' + sep + 'l1_fits' + sep         ; Path to the fits level 1 files

; Data file
data_file = l1fits_path + '2013-07-08_DOME_psd0_DIT-13ms_IMG.fits'

; Extract detector image
img_tmp = READFITS(data_file, header, /SILENT)

; Take only one frame for testing
img_in  = img_tmp

; Remove pickup noise
img_out = PICKUP_REMOVAL(img_in[*,*,1], FREQ_CUT=5., PIXSIZE=18., PERIOD=0.055, NO_GPU=1, /PLOT)

END