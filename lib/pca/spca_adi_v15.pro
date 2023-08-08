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
function spca_adi_v15, cx_init, cy_init, dim, dim_init, fwhm, plsc, truncate_pca, rin_init, step_init, delta, right_handed, normal=normal, display=display,  verbose=verbose, wcs=wcs

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

;Rounding and formatting PCA input params
;----------------------------------------
truncate_pca=float(round(truncate_pca))
dim=float(dim)
fwhm_rounded=ceil(fwhm)

if keyword_set(verbose) then begin
print, '**************'
print, 'Performing PCA'
print, '**************'
print, 'Double-check PCA params (truncate_pca, delta, rin_init, step_init): ', truncate_pca, delta, rin_init, step_init
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
mask= mask_t ge (fwhm_rounded*rin_init)
mask(where(mask eq 0))=0
for i=0,nobj-1 do begin
obj[*,*,i]=obj[*,*,i]*mask
endfor
endif

;Initializing PCA geometrical parameters
objtmp=fltarr(dim*dim,nobj)
mean_obj=fltarr(dim,dim,nobj)
objt=fltarr(dim,dim,nobj)
t=systime(1)
rho=shift(dist(dim),dim/2,dim/2)

rin=rin_init*fwhm_rounded;in pixel, starting at rin_init
rout=dim/2-1;finishing at edge of input image
n_annuli=floor((rout-rin)/(fwhm_rounded*step_init))
delta_deg=fltarr(n_annuli)

;Making temporary directory
if keyword_set(verbose) then begin
print, '********************************************************************************'
print, 'Del (for safety, should be already del if executed successfully), mk tmpb subdir'
print, '********************************************************************************'
endif
spawn, 'pwd', path
spawn, 'rm -d -r tmpb'
spawn, 'mkdir ' + path + '/tmpb/'


;Building zones
;--------------
for i=0,n_annuli-1 do begin
	;Annulus #i
	in = rin + floor((i*step_init*fwhm_rounded));-dr/2+1 ; inside limit of the zone s
	step = (rin + floor(((i+1)*step_init*fwhm_rounded))) - in ; radial thickness of the zone s (varying)
	out = in + step
	print, in, out, step
    index_annulus=where( (rho lt out) AND (rho ge in) )
    npix_annulus_i = n_elements(index_annulus)
   	obj_annulus_i= fltarr(npix_annulus_i,nobj)
     
	;Delta_deg for reference selection	
	delta_deg[i]= delta / (in/fwhm_rounded) / !dpi*180
	
	if keyword_set(verbose) then begin
	print, 'Cutting and saving annuli of size (pix, fwhm): ', strcompress(step),$
	strcompress(step/fwhm_rounded), ', at radius (pix, arcsec): ', strcompress(in), ', ', strcompress(in*plsc)
	endif
	
	for fr=0, nobj-1 do begin
		o=obj[*,*,fr]
		obj_annulus_i[*,fr]=o[index_annulus]
	endfor
	filename=strcompress(path) + '/tmpb/' + 'obj_'+tg_name_dyn+'_annulus_'+strcompress(i,/remove_all)+'.fits'	
	MWRFITS, obj_annulus_i, filename
	MWRFITS, index_annulus, filename		
	obj_annulus_i=0
endfor


;PCA core
for i=0,n_annuli-1 do begin
	
	obj_annulus_i = readfits(strcompress(path) + '/tmpb/' + 'obj_'+tg_name_dyn+'_annulus_'+$
							strcompress(i,/remove_all)+'.fits',exten=0)
    index_annulus = readfits(strcompress(path) + '/tmpb/' + 'obj_'+tg_name_dyn+'_annulus_'+$
							strcompress(i,/remove_all)+'.fits',exten=1)
    npix_annulus_i = n_elements(index_annulus)
    
	if keyword_set(verbose) then begin
	print, 'Doing sPCA on annuli, : ', i
	endif
	
	;Conditioning data (removing mean, necessary for PCA)
	
	;Remove temporal mean
	obj_annulus_i_mean=median(obj_annulus_i,dim=2)
	for npix=0, npix_annulus_i-1 do begin
	obj_annulus_i[npix,*]=obj_annulus_i[npix,*]-obj_annulus_i_mean[npix]
	endfor
    obj_annulus_i_mean=fltarr(nobj)
	;Remove spatial mean
	for fr=0, nobj-1 do begin
	    obj_annulus_i_mean[fr]=median(obj_annulus_i[*,fr])
		obj_annulus_i[*,fr]=obj_annulus_i[*,fr]-obj_annulus_i_mean[fr]
		;if keyword_set(verbose) then begin
	    ;print, 'Control of mean for annulus #, object # : ', obj_annulus_i_mean[fr], i , fr
	    ;endif
	endfor
	
 	for fr=0,nobj-1 do begin 
 	index_ref=0 & ref_library=0 & nref=0 & data=0 & covMatrix=0 & eigenval=0 & eigenvect=0 & pc=0
	;index_ref = where( (abs(paral) lt (abs(paral[fr])-delta_deg[i])) or (abs(paral) gt (abs(paral[fr])+delta_deg[i])))
	index_ref = where( ((paral) lt ((paral[fr])-delta_deg[i])) or ((paral) gt ((paral[fr])+delta_deg[i])))
	nref = n_elements(index_ref)
	if keyword_set(verbose) then begin
	print, 'Annulus #, Object #, Number of reference images in the library, : ', i, fr, nref
	wset, 0
	!psym=1
	plot, index_ref, index_ref, yrange=[0,nref-1], xrange=[0,nref-1] 
	;stop
	endif
	ref_library = obj_annulus_i[*,index_ref]	
	;PCA	
	data=transpose(ref_library)
	covMatrix = matrix_multiply(data, data, /BTRANSPOSE)
	;eigenval = LA_EIGENQL(covMatrix, EIGENVECTORS=eigenvect, range=[nref-truncate_pca,nref-1], /DOUBLE)
    eigenval = LA_EIGENQL(covMatrix, EIGENVECTORS=eigenvect, range=[nref-truncate_pca,nref-1], /DOUBLE)
    eigenval = reverse(eigenval)
    if truncate_pca eq 1 then eigenvect = transpose(eigenvect) else	eigenvect = reverse(eigenvect,2)
	if keyword_set(verbose) then print, 'PCA Eigenvalues: ', eigenval[0]
    ;Principal component 
 	pc = matrix_multiply(eigenvect,data,/ATRANSPOSE)
 	;pc = transpose(pc)
 	;stop
 	supref=fltarr(npix_annulus_i)	
 	
 		for k=0,truncate_pca-1 do begin
 		coeff=0
		pc[k,*]=pc[k,*] / sqrt(eigenval[k])
		coeff=total(reform(obj_annulus_i[*,fr]*reform(pc[k,*])))
		;coeff=median(reform(obj_annulus_i[*,fr]*reform(pc[k,*])))*nref
		supref=supref+coeff*reform(pc[k,*])
		endfor
	
 	obj_annulus_i[*,fr]=obj_annulus_i[*,fr]-supref
 	test_mean_supref=mean(supref)
 	test_mean_ann=mean(obj_annulus_i[*,fr])
 		if keyword_set(verbose) then begin
 		print, 'Test means supref, obj_annulus_i', test_mean_supref, test_mean_ann  
 		endif
 	objtmp[index_annulus,fr]=obj_annulus_i[*,fr];+obj_annulus_i_mean[fr]
	endfor
		
endfor
print, 'Timing for annulus sPCA: ', systime(1)-t

objt=reform(objtmp,dim,dim,nobj)

;Reconstructing and Derotating

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
