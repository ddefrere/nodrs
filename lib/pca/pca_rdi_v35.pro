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
;Na: box size (typically 10-20 pixels)
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
;4: writefits, procdir+'supref_dc.fits', supreft
;5: writefits, procdir+'obj_minus_mean.fits', obj
;
;1: principal components datacube
;2: result => derotated median-combined reduced datacube (smoothed by FWHM and normalized)
;3: derotated, reduced datacube
;4: super reference (obtained by projection on the principal components) data cube
;5: normalized object
;
;Example
;-------
;dim=200
;fwhm=round(1.65e-6/8*206265/plsc)
;truncate_pca=8
;Na=16
;normal=1

;status=pca_adi_v3(dim, plsc, fwhm,$
;	 	truncate_pca, Na, normal=normal, /verbose)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function pca_rdi_v35, cx, cy, cxref, cyref, dim, fwhm, plsc, truncate_pca_low, truncate_pca_high, Na, rin_init, fake_companion, filter, display=display, recenter=recenter, normal=normal, remove_mean=remove_mean, verbose=verbose

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

;Rounding and formatting LOCI input params
;-----------------------------------------
truncate_pca=(truncate_pca_high-truncate_pca_low)+1
Na=float(Na)
dim=float(dim)

if keyword_set(verbose) then begin
print, '**************'
print, 'Performing PCA'
print, '**************'
print, 'Double-check PCA params (# refs imgs, Na): ', truncate_pca, Na
endif

;Merging object and ref
;----------------------
if file_test('tg_'+tg_name_dyn+'_dc.fits') eq 0 then begin
obj_ttmp=readfits('img_'+tg_name_dyn+'_obj'+'_dc.fits',header)
nobj=(size(obj_ttmp))[3]
obj_ttmp=obj_ttmp[cx-dim/2:cx+dim/2-1,cy-dim/2:cy+dim/2-1,0:nobj-1]

ref_ttmp=readfits('img_'+tg_name_dyn+'_ref'+'_dc.fits',header2)
nref=(size(ref_ttmp))[3]
obj_tmp=fltarr((size(obj_ttmp))[1],(size(obj_ttmp))[2],nobj+nref)
obj_tmp[*,*,0:nobj-1]=obj_ttmp
for i=0,nref-1 do begin
obj_tmp[*,*,nobj+i]=ref_ttmp[cxref-dim/2:cxref+dim/2-1,cyref-dim/2:cyref+dim/2-1,i];*mask_defect
endfor

;Recenter
if keyword_set(recenter) then  begin
tmp_ref=obj_tmp[*,*,0]
for i=0, nobj+nref-1 do begin
tmp=obj_tmp[*,*,i]
recenter_frames_v5, tmp, ref, dim/2,dim/2, fwhm, 6
obj_tmp[*,*,i]=tmp
endfor
SXADDPAR, header,'NREF', nref
writefits, 'img_'+tg_name_dyn+'_dc.fits', obj_tmp, header
endif
endif

if filter ne 0 then begin
obj=obj_tmp[*,*,0:nobj-1]
for i=0,nobj-1 do begin
obj[*,*,i]=median(obj[*,*,i]-median(obj[*,*,i],filter*fwhm),3)
endfor
ref=obj_tmp[*,*,nobj:nref+nobj-1]
for i=0, nref-1 do begin
ref[*,*,i]=median(ref[*,*,i]-median(ref[*,*,i],filter*fwhm),3)
endfor
endif else begin
obj=obj_tmp[*,*,0:nobj-1]
ref=obj_tmp[*,*,nobj:nref+nobj-1]
endelse

;Mask
if rin_init ne 0 then begin
mask_t=shift(dist(dim),dim/2,dim/2)
mask= mask_t ge (fwhm*rin_init)
mask(where(mask eq 0))=0
for i=0,nobj-1 do begin
obj[*,*,i]=obj[*,*,i]*mask
endfor
for i=0, nref-1 do begin
ref[*,*,i]=ref[*,*,i]*mask
endfor
endif


;Fake companion
;--------------
if total(fake_companion) ne 0 then begin
	print, '**********************************************************************************'
	print, 'Inj. fake comp. @ 0".09, 0".12, 0".15, 0".2, 0".25, wt flux: ', strcompress(fake_companion)
	print, '**********************************************************************************'
	fcp_t=readfits('fake_companion_naco_L_l27.fits')
	sz=(size(fcp_t))[2]
	fcp=fltarr(dim,dim)
	fcp[dim/2-sz/2:dim/2+sz/2-1,dim/2-sz/2:dim/2+sz/2-1]=fcp_t
	fc=fake_companion
	ks=(size(fc))[2]
	test=fltarr(dim,dim,nobj)
	for i=0, nobj-1 do begin
		for k=0, ks-1 do begin
		test[*,*,i]=test[*,*,i]+fc[0,k]*shifti(fcp,fc[1,k]/plsc,0)*normal
		obj[*,*,i]=obj[*,*,i]+fc[0,k]*shifti(fcp,fc[1,k]/plsc,0)*normal
		print, fc[0,k]*normal, fc[1,k]/plsc
		endfor
	endfor
	;stop
endif

;Loading photometry
;------------------
if normal eq 1 then begin
normal_v=readfits(datadir+'vec_'+tg_name+'_photometry.fits') 
endif else begin
normal_v=fltarr(nobj)
normal_v[*]=normal
endelse

;Initializing PCA geometrical parameters
pc_tile=fltarr(dim,dim,truncate_pca)
objt=fltarr(dim,dim,nobj)
supref_save=fltarr(dim,dim,nobj)
mean_obj=fltarr(dim,dim,nobj)
test_th=fltarr(dim,dim,nobj)

bs=floor(Na)
objt=fltarr(dim,dim,nobj)

for x=0,floor(dim/Na)-1 do begin
	for y=0, floor(dim/Na)-1 do begin
		data=0 & covMatrix=0 & eigenval=0 &	eigenvect=0 & index=0
		
		if keyword_set(remove_mean) then begin
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
		endif
		
	;PCA of REF
	data=transpose(reform(ref[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,*],bs*bs,nref))
	covMatrix = matrix_multiply(data, data, /BTRANSPOSE)
	;eigenval = LA_EIGENQL(covMatrix, EIGENVECTORS=eigenvect, EigenIndex=index, /DOUBLE)
	eigenval = LA_EIGENQL(covMatrix, EIGENVECTORS=eigenvect, range=[nref-1-truncate_pca,nref-1], /DOUBLE)
	eigenvect = reverse(eigenvect,2)
	eigenval = reverse(eigenval)

	if keyword_set(verbose) then print, 'PCA Eigenvalues: ', eigenval
    ;Principal component 
 	pc=fltarr(bs*bs,nref)
 	
 	for k=0,nobj-1 do begin
	;res[*,k]=data ## eigenvect[*,k] / sqrt(eigenval[k])
	pc[*,k]=matrix_multiply(eigenvect[*,k],data) / sqrt(eigenval[k])
	endfor
	
   	pc=Reform((pc), bs, bs, nref)
   	pc=pc(*,*,truncate_pca_low:truncate_pca_high)	  
	pc_tile[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,*]=pc[*,*,*]
	
	;Projection of target image on the KL basis 
	for i=0, nobj-1 do begin
	supref=fltarr(bs,bs)
	supref_test=fltarr(bs,bs)

		for k=0,truncate_pca-1 do begin
		proj_coeff=total(reform(obj[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i])*$
					reform(pc_tile[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,k]))
		supref=supref+proj_coeff*reform(pc_tile[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,k])
		endfor
	supref_save[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i]=supref
	;Subtraction of the super reference from the target image
	objt[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i]=reform(obj[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i])-supref
	
	if total(fake_companion) ne 0 then begin
	test_th[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i]=reform(test[x*bs:x*bs+bs-1,y*bs:y*bs+bs-1,i])-supref_test
	endif
	endfor
	endfor
	if keyword_set(display) then begin
	wset,0
	tvscl, congrid((objt[*,*,0]),512,512)>0<10
	endif
endfor

tg_name_dyn=tg_name_dyn+'_pca'
tg_name_bas=tg_name_ori+'_pca'

;Outputs
writefits, procdir+'pc_'+tg_name_dyn+'_dc.fits', pc_tile
writefits, procdir+'img_'+tg_name_dyn+'_median.fits', median(objt,dim=3)/normal, header
if total(fake_companion) ne 0 then begin
writefits, procdir+tg_name_dyn+'_test_th.fits', median(test_th,dim=3)/normal, header
endif
writefits, procdir+'img_'+tg_name_dyn+'_dc.fits', objt
writefits, procdir+'supref_dc.fits', supref_save
writefits, procdir+'obj_minus_mean.fits', obj

if keyword_set(display) then begin
wset,0
tvscl, congrid(median(objt,dim=3),512,512) 
endif

return, 1

END
