;+
; NAME:
;   SPCA_ADI_V13
;   
; PURPOSE:
;   Smart principal component analysis for an ADI image cube
;   
; CALLING SEQUENCE:
;   Result = SPCA_ADI_V13(cx_init, cy_init, dim, dim_init, fwhm, plsc, truncate_pca, rin_init, step_init, delta, right_handed, $
;                         NORMAL=normal, ONEZONE=onezone, DISPLAY=display, VERBOSE=verbose, WCS=wcs)
;                        
; DESCRIPTION:
;   Based on Soummer, Pueyo & Larkin, 2012, ApJ.
;   Implementation of the smart PCA for ADI cubes: 
;     - the FOV is divided into annuli
;     - for each annulus, each frame is considered separately to define a PSF library
;     - the PSF library is built frame by frame by choosing only the frames where the planet is far enough from its current position
;     - PCA is applied frame by frame in each annulus --> warning: CPU intensive!
;  
; INPUT VARIABLES:
;   cx_init      : x position of the star in the input image cube [pixels]
;   cy_init      : y position of the star in the input image cube [pixels]
;   dim          : dimension of output image, on which PCA will be carried out
;   dim_init     : dimension of input image
;   fwhm         : FWHM of the telescope PSF in pixels
;   plsc         : plate scale [arcsec per pixel]
;   truncate_pca : number of principal components retained (< nobj)
;   rin_init     : inner radius where PCA starts, in resolution elements (fwhm)
;   step_init    : number of FWHM in the definition of each of the annuli
;   delta        : requested separation of the planet in the library compared to current position [in FWHM]
;   right_handed : for right_handed WCS (rare, e.g. NICI)
;   
; INPUT KEYWORDS:
;   normal  : normalisation factor. If not set, no normalisation used. If set to 0, looks for 'vec_'+tg_name+'_photometry.fits'.
;   onezone : set this keyword to use only one radial zone for the PCA. The value given to "onezone" will set the location where the sPCA is to be optimised [in fwhm].
;   display : set this keyword to display intermediate and final results
;   verbose : display various outputs/status in the console
;   wcs     : set this keyword for WCS data
;  
; COMMON BLOCK:
;   common_name : contains target name => implies that the data is available in datadir, 
;                 in the form 'img_'+tg_name+'_dc.fits'
;                 the parallactic angle and photometry vectors (optional) should also be available in datadir 
;                 in the form 'vec_'+tg_name+'_paral.fits' and 'vec_'+tg_name+'_photometry.fits'
; 
; OUTPUTS:
;   Returns the final, derotated, median image after sPCA processing.
;   In addition, two fits outputs: 
;     - the same image is saved in a fits file ("*_median.fits") in the "procdir" directory
;     - the full processed cube before derotation is also saved in a fits file ("*_dc.fits") in the "procdir" directory
;   
;  DEPENDENCIES:
;    AstroLib: READFITS, WRITEFITS 
;
; MODIFICATION HISTORY:
;   Version 1.0, 21-MAY-2013, by Dimitri Mawet, ESO - JPL
;   Version 1.1, 24-MAY-2013, OA: added header, cleaned up code
;   Version 1.2, 06-JUN-2013, OA: added keyword ONEZONE, removed bug in computation of parallactic angle exclusion
;   Version 1.3, 30-JUL-2013, OA: removed rounding of FWHM
;-

FUNCTION SPCA_ADI_V13, cx_init, cy_init, dim, dim_init, fwhm, plsc, truncate_pca, rin_init, step_init, delta, right_handed, $
                       NORMAL=normal, ONEZONE=onezone, DISPLAY=display, VERBOSE=verbose, WCS=wcs

COMMON common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir


;Rounding and formatting PCA input params
;----------------------------------------
truncate_pca=float(round(truncate_pca))
dim=float(dim)
fwhm_rounded=ceil(fwhm)

if keyword_set(verbose) then begin
  print, '***************'
  print, 'Performing sPCA'
  print, '***************'
  print, 'Double-check sPCA params (truncate_pca, delta, rin_init, step_init): ', truncate_pca, delta, rin_init, step_init
endif


;Reading object input file
;-------------------------
obj=readfits(procdir+'img_'+tg_name_dyn+'_dc.fits',header)
nobj=(size(obj))[3]


;Parallactic angle loading
;-------------------------
paral=readfits(datadir+'vec_'+tg_name_bas+'_paral.fits'); reading parallactic angle
if keyword_set(display) then begin
  wset, 0
  plot, paral, xtitle='Frame #', ytitle='Parallactic Angle'
endif


; Load photometry
; ---------------
IF NOT KEYWORD_SET(normal) THEN normal = 1
IF normal EQ 0 THEN normal_v = READFITS(procdir+'vec_'+tg_name_bas+'_photometry.fits') $
ELSE normal_v = FLTARR(nobj) + normal  


;Mask center (for coronagraphy or saturated images)
;--------------------------------------------------
if rin_init ne 0 then begin
  mask_t=shift(dist(dim),dim/2,dim/2)
  mask= mask_t GE fwhm*rin_init
  for i=0,nobj-1 do obj[*,*,i]=obj[*,*,i]*mask
endif


;Initializing PCA geometrical parameters
;---------------------------------------
objtmp=fltarr(dim*dim,nobj)
mean_obj=fltarr(dim,dim,nobj)
objt=fltarr(dim,dim,nobj)
t=systime(1)
rho=shift(dist(dim),dim/2,dim/2) ; radial coordinate in pixels 

rin=rin_init*fwhm ;in pixel, starting at rin_init
rout=dim/2-1 ;finishing at edge of input image
n_annuli=floor((rout-rin)/(fwhm*step_init))
IF KEYWORD_SET(onezone) THEN n_annuli = 1 ; redefine one single optimisation zone 


;Building zones and performing the PCA
;-------------------------------------
for iann = 0,n_annuli-1 do begin
	;Define annulus #i
	in = rin + iann*step_init*fwhm;-dr/2+1 ; inside limit of the zone s in pixel
	IF KEYWORD_SET(onezone) THEN BEGIN
	  out = rout
	  step = out - in
	ENDIF ELSE BEGIN
  	step = (rin + (iann+1)*step_init*fwhm) - in ; radial thickness of the zone s in pixel (can be varying)
  	out = in + step
  ENDELSE
	print, in, out, step
  index_annulus = WHERE( (rho lt out) AND (rho ge in), npix_annulus )
 	obj_annulus = FLTARR(npix_annulus,nobj) ; create a 2D matrix that will contain the info of the considered annulus for each frame

	;Delta_deg for reference selection (requested separation in terms of azimuth)
	IF KEYWORD_SET(onezone) THEN delta_deg = delta / onezone / !dpi*180 $
	ELSE delta_deg = delta / (in/fwhm) / !dpi*180 ; computed at the inner edge of annulus for safety
	
	if keyword_set(verbose) then print, 'Cutting and saving annulus of '+strcompress(step)+'pix ('+strcompress(step/fwhm)+$
	                                    ' fwhm) at a radius of'+strcompress(in)+' pix ('+strcompress(in*plsc)+' arcsec).' 
	
	;Fill the 2D matrix with the pixels intensities from the considered annulus
	for fr=0, nobj-1 do obj_annulus[*,fr] = (obj[*,*,fr])[index_annulus]

  ;Conditioning data (removing mean, necessary for PCA)
  if keyword_set(verbose) then print, 'Doing sPCA on annulus '+strcompress(iann+1)+' of '+strcompress(n_annuli)
  data=0 & covMatrix=0 & eigenval=0 & eigenvect=0 
	obj_annulus_mean=fltarr(nobj)
	for fr=0, nobj-1 do begin
	  obj_annulus_mean[fr]=mean(obj_annulus[*,fr])
		obj_annulus[*,fr]=obj_annulus[*,fr]-obj_annulus_mean[fr]
	endfor
	
	; Doing the sPCA analysis frame by frame
 	for fr=0,nobj-1 do begin 
  	index_ref = where((paral lt paral[fr]-delta_deg) or (paral gt paral[fr]+delta_deg)) ; pick out the frames that are not in the parallactic exclusion zone 
  	nref = n_elements(index_ref)
  	print, strcompress(iann+1)+'/'+strcompress(n_annuli), fr, nref
  	data = obj_annulus[*,index_ref]
  	
  	;PCA	
  	covMatrix = matrix_multiply(data, data, /ATRANSPOSE)
  	eigenval = LA_EIGENQL(covMatrix, EIGENVECTORS=eigenvect, range=[nref-truncate_pca,nref-1], /DOUBLE)
    eigenval = reverse(eigenval)
  	eigenvect = reverse(eigenvect,2)
  	;if keyword_set(verbose) then print, 'PCA Eigenvalues: ', eigenval
  	
    ;Principal components 
   	pc = matrix_multiply(eigenvect, data, /ATRANSPOSE, /BTRANSPOSE)
   	supref=fltarr(npix_annulus)	
   	for k=0,truncate_pca-1 do begin
  		pc[k,*]=pc[k,*] / sqrt(eigenval[k])
  		coeff=total(obj_annulus[*,fr]*reform(pc[k,*]))
  		;print, coeff
  		supref=supref+coeff*reform(pc[k,*])
  	endfor
  	
   	objtmp[index_annulus,fr]=obj_annulus[*,fr]-supref;+obj_annulus_mean[fr]
	endfor
		
endfor
print, 'Timing for annulus sPCA: ', systime(1)-t

objt=reform(objtmp,dim,dim,nobj)


;Reconstructing and derotating
;-----------------------------
for i=0, nobj-1 do begin
if keyword_set(verbose) then print, 'Derotate by: ', paral[i]
	if right_handed eq 1 then paral[i]=-paral[i]
		objt[*,*,i]=rot(objt[*,*,i],-paral[i],1.0,dim/2,dim/2,cubic=-0.5,/pivot)
endfor

tg_name_dyn=tg_name_dyn+'_pca'


;Final images and outputs
;------------------------
if keyword_set(wcs) then begin
  resm_fin=fltarr(dim_init,dim_init)
  resm_fin[dim_init/2-dim/2:dim_init/2+dim/2-1,dim_init/2-dim/2:dim_init/2+dim/2-1]=median(objt,dim=3)/median(normal_v)
  resm_fin=shift(resm_fin,cx_init-dim_init/2,cy_init-dim_init/2)
  writefits, procdir+'img_'+tg_name_dyn+'_median.fits',resm_fin, header
endif else begin
  resm_fin=median(objt,dim=3)/median(normal_v)
  writefits, procdir+'img_'+tg_name_dyn+'_median.fits',resm_fin, header
endelse

writefits, procdir+tg_name_dyn+'_dc.fits', objt

if keyword_set(display) then begin
  wset,0
  tvscl, congrid(median(objt,dim=3),512,512) 
endif

return, resm_fin

END
