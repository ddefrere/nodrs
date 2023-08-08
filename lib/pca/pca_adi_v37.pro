;+
; NAME:
;   PCA_ADI_V37
;   
; PURPOSE:
;   Principal component analysis for an ADI image cube
;   
; CALLING SEQUENCE:
;   Result = PCA_ADI_V37(cx_init, cy_init, dim, dim_init, fwhm, truncate_pca, Na, rin_init, right_handed, $
;                        normal=normal, display=display, verbose=verbose, wcs=wcs)
;                        
; DESCRIPTION:
;   Based on Soummer, Pueyo & Larkin, 2012, ApJ
;  
; INPUT VARIABLES:
;   cx_init      : x position of the star in the input image cube [pixels]
;   cy_init      : y position of the star in the input image cube [pixels]
;   dim          : dimension of output image, on which PCA will be carried out
;   dim_init     : dimension of input image
;   fwhm         : FWHM of the telescope PSF in pixels
;   truncate_pca : number of principal components retained (< nobj)
;   Na           : number of pixels inside the PCA tile (note: an odd number of boxes should better be used if star centred in the images)
;   rin_init     : inner radius where PCA starts, in resolution elements (fwhm) 
;   right_handed : for right_handed WCS (rare, e.g. NICI)
;   
; INPUT KEYWORDS:
;   normal  : if set to 1, looks for 'vec_'+tg_name+'_photometry.fits'
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
;   Returns the final, derotated, median image after PCA processing.
;   In addition, two fits outputs: 
;     - the same image is saved in a fits file ("*_median.fits") in the "procdir" directory
;     - the full processed cube before derotation is also saved in a fits file ("*_dc.fits") in the "procdir" directory
;   
;  DEPENDENCIES:
;    AstroLib: READFITS, WRITEFITS 
;
; MODIFICATION HISTORY:
;   Version 1.0, 01-JUL-2012, by Dimitri Mawet, ESO - JPL
;   (...)
;   Version 3.6, 16-MAY-2013, OA: added header, cleaned up code
;   Version 3.7, 06-JUN-2013, OA: now returns the final median image as an output
;-


function pca_adi_v37, cx_init, cy_init, dim, dim_init, fwhm, truncate_pca, Na, rin_init, right_handed, $
                      normal=normal, display=display,  verbose=verbose, wcs=wcs

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

;Rounding and formatting PCA input params
;----------------------------------------
truncate_pca=float(round(truncate_pca))
Na=float(Na)
dim=float(dim)

if keyword_set(verbose) then begin
print, '**************'
print, 'Performing PCA'
print, '**************'
print, 'Double-check PCA params (# refs imgs, Na): ', truncate_pca, Na
endif

;Reading object input file
;-------------------------
obj=readfits(procdir+'img_'+tg_name_dyn+'_dc.fits',header)
nobj=(size(obj))[3]

;Parallactic angle loading
;-------------------------
paral=readfits(procdir+'vec_'+tg_name_bas+'_paral.fits'); reading parallactic angle
if keyword_set(display) then begin
  wset, 0
  plot, paral, xtitle='Frame #', ytitle='Parallactic Angle'
endif

;Loading photometry
;------------------
if normal eq 0 then begin
  normal_v=readfits(procdir+'vec_'+tg_name_bas+'_photometry.fits') 
endif else begin
  normal_v=fltarr(nobj)
  normal_v[*]=normal
endelse

;Mask center (for coronagraphy or saturated images)
if rin_init ne 0 then begin
  mask_t=shift(dist(dim),dim/2,dim/2)
  mask= mask_t ge (fwhm*rin_init)
  mask(where(mask eq 0))=0
  for i=0,nobj-1 do obj[*,*,i]=obj[*,*,i]*mask
endif

;Initializing PCA geometrical parameters
objt=fltarr(dim,dim,nobj)
mean_obj=fltarr(dim,dim,nobj)
bs=floor(Na)
objt=fltarr(dim,dim,nobj)
t=systime(1)

for x=0,floor(dim/Na)-1 do begin
	for y=0, floor(dim/Na)-1 do begin
		data=0 & covMatrix=0 & eigenval=0 &	eigenvect=0 
		;Conditioning data (removing mean, necessary for PCA)
		for i=0, nobj-1 do begin
		dum=0
		dum=reform(obj[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i])
		dum=dum-mean(dum)
		obj[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i]=dum
		endfor
	;PCA	
	data=transpose(reform(obj[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,*],bs*bs,nobj))
	covMatrix = matrix_multiply(data, data, /BTRANSPOSE)
	eigenval = LA_EIGENQL(covMatrix, EIGENVECTORS=eigenvect, range=[nobj-1-truncate_pca,nobj-1], /DOUBLE)
	eigenvect = reverse(eigenvect,2)
	;eigenval = EIGENQL(covMatrix, EIGENVECTORS=eigenvect, /DOUBLE)
	if keyword_set(verbose) then print, 'PCA Eigenvalues: ', eigenval
    ;Principal component 
 	pc=fltarr(nobj,bs*bs)
 	pc = matrix_multiply(eigenvect,data,/ATRANSPOSE)
 	supref = reform(matrix_multiply(eigenvect[*,0:truncate_pca-1], pc[0:truncate_pca-1,*]),nobj, bs, bs)
 	for i=0,nobj-1 do begin 
 	objt[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i]=reform(obj[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i])-reform(supref[i,*,*])
	endfor
	endfor
endfor
print, 'Timing for optimized tiled PCA: ', systime(1)-t

;Derotating
for i=0, nobj-1 do begin
if keyword_set(verbose) then print, 'Derotate by: ', paral[i]
	if right_handed eq 1 then paral[i]=-paral[i]
		objt[*,*,i]=rot(objt[*,*,i],-paral[i],1.0,dim/2,dim/2,cubic=-0.5,/pivot)
endfor

tg_name_dyn=tg_name_dyn+'_pca'
;tg_name_ori=tg_name_ori+'_pca'

;Outputs
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
