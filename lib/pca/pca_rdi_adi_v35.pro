;-------------------------------------------------------------------
;Principal Component Analysis, following Soummer, Pueyo, Larkin 2012
;-------------------------------------------------------------------
;v3, D. Mawet, (ESO - JPL), July 2012
;
;Input parameters
;----------------
;dim: dimension of input array (assumes star is centered at dim/2, dim/2)
;fwhm: in pixels
;truncate_cpa: number of principal components retained (< nobj)
;Na: box size (typically 8 pixels - whole image)
;/right_hand: for right_handed WCS (rare, e.g. NICI)
;/display: display on
;normalisation factor (photometry): if set to 1, looks for 'vec_'+tg_name+'_photometry.fits'
;/verbose: display various outputs/status
;
;Common handling: contains target name => implies that the data is available in datadir, 
;in the form 'img_'+tg_name+'_dc.fits'
;the parallactic angle and photometry vectors (optional) should also be available in datadir 
;in the form 'vec_'+tg_name+'_paral.fits' and 'vec_'+tg_name+'_photometry.fits'
;
;Output
;------
;1: writefits, procdir+'proj_'+tg_name_dyn+'_dc.fits', proj_obj_dc
;2: writefits, procdir+tg_name_dyn+'_final.fits', filter_image(median(objt,dim=3),fwhm=fwhm)/normal, header
;3: writefits, procdir+tg_name_dyn+'_dc.fits', objt
;
;1: principal components datacube
;2: result => derotated median-combined reduced datacube (smoothed by FWHM and normalized)
;3: derotated, reduced datacube
;
;Example
;-------
;dim=200
;fwhm=round(1.65e-6/8*206265/plsc)
;truncate_pca=8
;Na=64
;normal=1

;status=pca_adi_v35(cx_init, cy_init, dim, dim_init, fwhm,$
;	 	truncate_pca, Na, rin_init, right_handed, /speed, normal=1, /verbose)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function pca_rdi_adi_v35, ref_name, cx_init, cy_init, dim, dim_init, fwhm, truncate_pca_init, Na, rin_init, right_handed, normal=normal, display=display,  verbose=verbose, wcs=wcs

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

;Rounding and formatting PCA input params
;----------------------------------------
truncate_pca_init=float(round(truncate_pca_init))
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

ref=readfits(datadir+ref_name)
nref=(size(ref))[3]

;Parallactic angle loading
;-------------------------
paral=readfits(datadir+'vec_'+tg_name_bas+'_paral.fits'); reading parallactic angle
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
for i=0,nobj-1 do begin
obj[*,*,i]=obj[*,*,i]*mask
endfor
for i=0,nref-1 do begin
ref[*,*,i]=ref[*,*,i]*mask
endfor
endif

;Initializing PCA geometrical parameters
objt=fltarr(dim,dim,nobj)
mean_obj=fltarr(dim,dim,nobj)
bs=floor(Na)
objt=fltarr(dim,dim,nobj)
t=systime(1)

for x=0,floor(dim/Na)-1 do begin
for y=0, floor(dim/Na)-1 do begin

if x eq floor(dim/Na)/2 and y eq floor(dim/Na)/2 then truncate_pca=truncate_pca_init*6 else truncate_pca=truncate_pca_init


		data=0 & covMatrix=0 & eigenval=0 &	eigenvect=0 
		;Conditioning data (removing mean, necessary for PCA)
		for i=0, nobj-1 do begin
		dum=0
		dum=reform(obj[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i])
		dum=dum-mean(dum)
		obj[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i]=dum
		endfor
		for i=0, nref-1 do begin
		dum=0
		dum=reform(ref[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i])
		dum=dum-mean(dum)
		ref[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i]=dum
		endfor

	;PCA	
	data_ref=transpose(reform(ref[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,*],bs*bs,nref))
;	data_obj=transpose(reform(obj[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,*],bs*bs,nobj))
	covMatrix = matrix_multiply(data_ref, data_ref, /BTRANSPOSE)
	eigenval = LA_EIGENQL(covMatrix, EIGENVECTORS=eigenvect, range=[nref-truncate_pca,nref-1], /DOUBLE)
	
	eigenval = reverse(eigenval)
	eigenvect = reverse(eigenvect,2)
	;eigenval = EIGENQL(covMatrix, EIGENVECTORS=eigenvect, /DOUBLE)
	if keyword_set(verbose) then print, 'PCA Eigenvalues: ', eigenval
    ;Principal component 
 	pc=fltarr(truncate_pca,bs*bs)
 	pc = matrix_multiply(eigenvect,data_ref,/ATRANSPOSE)
 	pc = reform(transpose(pc), bs, bs, truncate_pca)

 	for k=0,truncate_pca-1 do begin
	pc[*,*,k]=pc[*,*,k] / sqrt(eigenval[k])
	endfor 	
  	;stop
  	for i=0,nobj-1 do begin 
  	    supref=fltarr(bs,bs)	
 		for k=0,truncate_pca-1 do begin
		coeff=total(reform(obj[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i])*reform(pc[*,*,k]))
		supref=supref+coeff*reform(pc[*,*,k])
	;	stop
		endfor
 	objt[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i]=reform(obj[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i])-reform(supref)
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

return, 1

END
