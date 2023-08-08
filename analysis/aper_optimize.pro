; +
; NAME: APER_OPTIMIZE
;   
; PURPOSE
;   Test routine to compute the optimum radius for aperture photometry (+ background annulus radii).
;
; LAST MODIFICATION
;   17-APR-2015, Denis Defrère, University of Arizona, ddefrere@email.arizona.edu

PRO APER_OPTIMIZE

; Running parameters
n_pt     = 500
n_aper   = 1000
m2r      = 4.848136811D-9        ; mas to rad conversion factor

; Telescope and detector parameters
tel_diam = 8.4
cen_obs  = 0.1
lambda   = 1D-5
ratio    = cen_obs/tel_diam
pix_size = 18.                ; mas
rms_pix  = 1.

; Target info
star_flux = 1.

; Define the aperture radius array
aper_rad = (DINDGEN(n_aper)+1)/n_aper*20.*lambda/tel_diam      ; rad
aper_stp = aper_rad[1]-aper_rad[0]
nsky     = 100.

; Define limits
aper_lim  = 8*pix_size*m2r/aper_stp
ring_ilim = 3800*m2r/aper_stp
ring_olim = 4200*m2r/aper_stp

; Loop over the aperture radii
snr1 = DBLARR(n_aper) & snr2 = DBLARR(n_aper)
FOR i_aper = 0, n_aper-1 DO BEGIN
  ; Compute encyrcled energy
  r           = 2.*!Dpi/lambda*0.5*tel_diam*aper_rad[i_aper]   
  t           = (1+DINDGEN(n_pt))/n_pt*r
  energy_frac = 1./(1.-ratio^2)*(1.-BESELJ(r,0)^2-BESELJ(r,1)^2+ratio^2*(1-BESELJ(ratio*r,0)^2-BESELJ(ratio*r,1)^2)-$
                4.*ratio*INT_TABULATED(t,BESELJ(t,1)*BESELJ(ratio*t,1)/t))
  ;PRINT, aper_rad, energy_frac
  signal      = star_flux*energy_frac
  
  ; Store values
  IF i_aper LT aper_lim  THEN flx_apr  = energy_frac
  IF i_aper LT ring_ilim THEN flx_rin  = energy_frac
  IF i_aper LT ring_olim THEN flx_rout = energy_frac
  
  ; Compute number of pixel in this region
  n_pix       = !Dpi*(aper_rad[i_aper]/(pix_size*m2r))^2
  noise1      = rms_pix*SQRT(n_pix)
  noise2      = rms_pix/nsky*n_pix
  
  ; Compute SNR
  snr1[i_aper] =  signal/noise1 
  snr2[i_aper] =  signal/noise2                
ENDFOR

; PRINT info
PRINT,aper_lim, ring_ilim, ring_olim
PRINT, 'Energy enclosed in photometric aperture :', flx_apr
PRINT, 'Fraction in outer ring [%] :', (flx_rout-flx_rin)/flx_apr*100
r_surf = aper_lim^2/(ring_olim^2-ring_ilim^2)
PRINT, 'Fraction in a photometric aperture located in outer ring [%] :', r_surf*(flx_rout-flx_rin)/flx_apr*100
PRINT, 'surface ratio :', r_surf

; Find maximum
tmp     = MAX(snr1, idx_max)
sep_max = aper_rad[idx_max]/(lambda/tel_diam)

; Plot results
LOADCT, 0, /SILENT
fit     = 20./1720.
symsize = 0.9
yrange  = [0., 1.]
xrange  = [0.,5.]
plot_name = 'aperphot_optimize.eps'
PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME=plot_name
PLOTXY, aper_rad/(lambda/tel_diam), snr1/MAX(snr1), /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='', XTITLE='Photometric aperture radius [lam/D]', YTITLE='Normalized SNR', GRID=0, $
  CHARSIZE=1.0, CHARTHICK=2.5, THICK=4.5, XTHICK=4.5, YTHICK=4.5, YSTYLE=1, /NOERASE, /NODATA, WINDOW=[250,100,1350,900]*fit;, INSET_UR='a.'$
LOADCT, 13, /SILENT
PLOTXY, aper_rad/(lambda/tel_diam), snr1/MAX(snr1), THICK=4.5, COLOR=90, /ADD
PLOTXY, aper_rad/(lambda/tel_diam), snr2/MAX(snr2), THICK=4.5, COLOR=250, /ADD
XYOUTS, 0.5*xrange[1],0.78*yrange[1], 'Best aperture radius: ' + STRING(sep_max, FORMAT='(F4.2)') + 'lam/D', CHARSIZE=1.0, CHARTHICK=2.5
LEGEND, ['Background RMS','Mean background RMS'], LINESTYLE=[0,0], COLOR=[90,250], THICK=4.5, PTHICK=4.5, CHARSIZE=1.0, CHARTHICK=2.5, /TOP_LEGEND, /RIGHT_LEGEND
; Close plot
LOADCT, 0, /SILENT
PLOTXY, /FIN

; Compute optimal sky radius
n_sky = 100
nsky  = (DINDGEN(n_sky)+1)/(n_sky)*10.

; Number of pixels in the aperture
n_pix1 = !Dpi*(aper_rad[idx_max]/(pix_size*m2r))^2
n_pix2 = !Dpi*(0.2*aper_rad[idx_max]/(pix_size*m2r))^2
  
; Loop over the sky values
noise_ratio1 = DBLARR(n_sky)
noise_ratio2 = DBLARR(n_sky)
FOR i_sky = 0, n_sky-1 DO BEGIN
  noise_ratio1[i_sky] = (rms_pix/nsky[i_sky])/(rms_pix*SQRT(n_pix1))
  noise_ratio2[i_sky] = (rms_pix/nsky[i_sky])/(rms_pix*SQRT(n_pix2))
ENDFOR

; Plot results
LOADCT, 0, /SILENT
fit     = 20./1720.
symsize = 0.9
yrange  = [0.01, MAX(noise_ratio2)]
xrange  = [0.,10.]
plot_name = 'aperphot_optimize2.eps'
PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME=plot_name
PLOTXY, nsky, noise_ratio1, /YLOG, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='', XTITLE='Number of sky pixels/number of aperture pixels', YTITLE='Mean background RMS/background RMS', GRID=0, $
  CHARSIZE=1.0, CHARTHICK=2.5, THICK=4.5, XTHICK=4.5, YTHICK=4.5, YSTYLE=1, /NOERASE, /NODATA, WINDOW=[250,100,1350,900]*fit;, INSET_UR='a.'$
LOADCT, 13, /SILENT
PLOTXY, nsky, noise_ratio1, THICK=4.5, COLOR=90, /ADD
PLOTXY, nsky, noise_ratio2, THICK=4.5, COLOR=250, /ADD
LEGEND, ['Aperture radius: ' + STRING(sep_max, FORMAT='(F4.2)') + ' lam/D','Aperture radius: ' + STRING(0.2*sep_max, FORMAT='(F4.2)') + 'lam/D'], LINESTYLE=[0,0], COLOR=[90,250], THICK=4.5, PTHICK=4.5, CHARSIZE=1.0, CHARTHICK=2.5, /TOP_LEGEND, /RIGHT_LEGEND
; Close plot
LOADCT, 0, /SILENT
PLOTXY, /FIN

; Define background limits
n_bck = 100.
rin   = 0.5+(DINDGEN(n_bck)+1)/n_bck*3.5      ; in lambda/D
rout  = 2.00+(DINDGEN(n_bck)+1)/n_bck*6.0

; Loop over the size
bck_avg = DBLARR(n_bck,n_bck)
FOR i_rin = 0, n_bck-1 DO BEGIN
  FOR i_rout = i_rin+1, n_bck-1 DO BEGIN
    theta   = (rin[i_rin] + (DINDGEN(n_aper))/(n_aper-1)*(rout[i_rout]-rin[i_rin]))*lambda/tel_diam
    r       = 2.*!Dpi/lambda*0.5*tel_diam*theta
    bck_avg[i_rin,i_rout] = 2.*!Dpi*INT_TABULATED(theta,theta*(2*BESELJ(r,1)/r)^2)/(!Dpi*(MAX(theta)^2-MIN(theta)^2))   ; *theta because this is a polar intagration!
  ENDFOR
ENDFOR

 ; Plot results
LOADCT, 13, /SILENT
fit     = 20./1720.
symsize = 0.9
xrange  = [MIN(rin), MAX(rin)]
yrange  = [MIN(rout), MAX(rout)]
plot_name = 'aperphot_optimize3.eps'
PLOTXY, /INIT, /EPS, /INV, WINDOW=[0, 0, 1400, 1000]*fit, FILENAME=plot_name
PLOTXY, bck_avg, /NEW, XRANGE=xrange, YRANGE=yrange, TITLE='', XTITLE='Inner radius (!4k!3/D)', YTITLE='Outer radius (!4k!3/D)', GRID=0, $
        CHARSIZE=1.2, CHARTHICK=2.5, THICK=4.5, XTHICK=4.5, YTHICK=4.5, YSTYLE=1, /NOERASE, WINDOW=[250,120,1050,920]*fit;, INSET_UR='a.'$
!X.TITLE=''
!Y.TITLE=''
COLORBAR, CHARSIZE=1.2, COLOR=color, DIVISIONS=divisions,FORMAT='(F5.2)', POSITION=[0.78, 0.12, 0.83, 0.92], MAXRANGE=100.*MAX(bck_avg), MINRANGE=100.*MIN(bck_avg), TITLE='Relative error [%]', /RIGHT, /VERTICAL
LOADCT, 0, /SILENT
CONTOUR, bck_avg, rin, rout, LEVELS=[0.1, 1.0, 3.0]*1D-2, C_ANNOTATION=['0.1%','1.0%','3.0%'], COLOR=255, /OVERPLOT
PLOTXY, /FIN
END