PRO merging_window_v5b, cx, cy, cxref, cyref, dim, resize_only=resize_only

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

;Merging object and ref
;----------------------
if keyword_set(resize_only) then begin
obj_ttmp=readfits(procdir+'img_'+tg_name_bas+'_dc.fits',header)
	if (size(obj_ttmp))[1] ne dim then begin
	tg_name_dyn=tg_name_dyn+'_rsz'
	nobj=(size(obj_ttmp))[3]
	obj_ttmp=obj_ttmp[cx-dim/2:cx+dim/2-1,cy-dim/2:cy+dim/2-1,0:nobj-1]
	writefits, procdir+'img_'+tg_name_dyn+'_dc.fits', obj_ttmp, header
	endif else begin
	tg_name_dyn=tg_name_dyn+'_rsz'
	writefits, procdir+'img_'+tg_name_dyn+'_dc.fits', obj_ttmp, header
	endelse
	
endif else begin
	;This part has to be revised, might be obsolete.
	if file_test(procdir+'img_'+tg_name_dyn+'_dc.fits') eq 0 then begin
	obj_ttmp=readfits(procdir+'tg_'+tg_name_ori+'_obj'+'_dc.fits',header)
	nobj=(size(obj_ttmp))[3]
	obj_ttmp=obj_ttmp[cx-dim/2:cx+dim/2-1,cy-dim/2:cy+dim/2-1,0:nobj-1]
	ref_ttmp=readfits(procdir+'tg_'+tg_name_ori+'_ref'+'_dc.fits',header2)
	nref=(size(ref_ttmp))[3]
	ref_ttmp=ref_ttmp[cxref-dim/2:cxref+dim/2-1,cyref-dim/2:cyref+dim/2-1,0:nref-1]
	obj_tmp=fltarr((size(obj_ttmp))[1],(size(obj_ttmp))[2],nobj+nref)
	obj_tmp[*,*,0:nobj-1]=obj_ttmp
	obj_tmp[*,*,nobj:*]=ref_ttmp
	SXADDPAR, header,'NREF', nref
	writefits, procdir+'img_'+tg_name_dyn+'_dc.fits', obj_tmp, header
	endif else begin
	obj_ttmp=readfits(procdir+'tg_'+tg_name_dyn+'_dc.fits',header)
		if (size(obj_ttmp))[3] ne dim then begin
		print, 'Warning: Image dimension has changed, remerging'
		obj_ttmp=readfits(procdir+'tg_'+tg_name_ori+'_obj'+'_dc.fits',header)
		nobj=(size(obj_ttmp))[3]
		obj_ttmp=obj_ttmp[cx-dim/2:cx+dim/2-1,cy-dim/2:cy+dim/2-1,0:nobj-1]		
		ref_ttmp=readfits(procdir+'tg_'+tg_name_ori+'_ref'+'_dc.fits',header2)
		nref=(size(ref_ttmp))[3]
		ref_ttmp=ref_ttmp[cxref-dim/2:cxref+dim/2-1,cyref-dim/2:cyref+dim/2-1,0:nref-1]		
		obj_tmp=fltarr((size(obj_ttmp))[1],(size(obj_ttmp))[2],nobj+nref)
		obj_tmp[*,*,0:nobj-1]=obj_ttmp
		obj_tmp[*,*,nobj:*]=ref_ttmp
		SXADDPAR, header,'NREF', nref
		writefits, procdir+'img_'+tg_name_dyn+'_dc.fits', obj_tmp, header
		endif
	endelse
endelse
END