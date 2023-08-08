function adinirc2_v5, cx, cy, dim

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

print, '***************************************'
print, 'Performing basic ADI median subtraction'
print, '***************************************'
img_dc=readfits(procdir+'img_'+tg_name_dyn+'_dc.fits', Headers_init)
paral_dc_parang=readfits(procdir+'vec_'+tg_name_bas+'_paral.fits')
normal=readfits(procdir+'vec_'+tg_name_bas+'_photometry.fits')
tg_name_dyn=tg_name_dyn + '_adi'
;tg_name_ori=tg_name_ori + '_adi'

print, '---------------------------------------------'
print, 'Computing median pupil-pinned speckle pattern'
print, '---------------------------------------------'
nobj=(size(img_dc))[3]; Number of frames in the input object datacube

for i=0,nobj-1 do begin
img_dc[*,*,i]=img_dc[*,*,i]/normal[i]
endfor

img_med=median(img_dc,dim=3)

res=fltarr(dim,dim,nobj)
res_tmp=fltarr(dim,dim)
res_fin=fltarr(dim,dim)
for i=0,nobj-1 do begin
obj_tmp=img_dc[*,*,i]
print, 'Subtracting median from Frame # ', strcompress(i), ', derotating and adding'
res[*,*,i]=rot(obj_tmp-img_med,-paral_dc_parang[i],1.0,cx,cy,cubic=-0.5,/pivot)
res_tmp=res_tmp+res[*,*,i]/nobj
wset,0
tvscl, congrid(res_tmp,500,500)
endfor
res_fin=median(res,dim=3);/median(normal)

SXADDPAR, Headers_init,'PARANG', 0
writefits, procdir+'img_'+tg_name_dyn+'.fits',res_fin, Headers_init
writefits, procdir+'img_'+tg_name_dyn+'_median.fits',img_med
writefits, procdir+'img_'+tg_name_dyn+'_dc_minus_median.fits',res, Headers_init

return, res_fin

END