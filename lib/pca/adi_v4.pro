PRO adi_v4, cx, cy, dim, fwhm,normal=normal, right_handed=right_handed

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

print, '***************************************'
print, 'Performing basic ADI median subtraction'
print, '***************************************'
img_dc=readfits(procdir+'img_'+tg_name_dyn+'_dc.fits', Headers_init)
paral=readfits(procdir+'vec_'+tg_name_ori+'_paral.fits')
;normal=readfits(procdir+'vec_'+tg_name_ori+'_photometry.fits')
tg_name_dyn=tg_name_dyn + '_adi'
tg_name_bas=tg_name_ori + '_adi'

print, '---------------------------------------------'
print, 'Computing median pupil-pinned speckle pattern'
print, '---------------------------------------------'
nobj=(size(img_dc))[3]; Number of frames in the input object datacube

for i=0,nobj-1 do begin
img_dc[*,*,i]=img_dc[*,*,i]/normal
endfor

img_med=median(img_dc,dim=3)

res=fltarr(dim,dim,nobj)
res_tmp=fltarr(dim,dim)
res_fin=fltarr(dim,dim)
for i=0,nobj-1 do begin
tmp=img_dc[*,*,i]
;destrip
mask_strip_th=2*abs(stddev(tmp))
mask_strip=tmp ge (-mask_strip_th) and tmp le (mask_strip_th)
tmp_strip=replicate(1d, dim)#(median(tmp*mask_strip,dim=2))
obj_tmp=tmp-transpose(tmp_strip)


print, 'Subtracting median from Frame # ', strcompress(i), ', derotating and adding'
if keyword_set(right_handed) then paral[i]=-paral[i]
res[*,*,i]=rot(obj_tmp-img_med,-paral[i],1.0,cx,cy,cubic=-0.5,/pivot)
;res[*,*,i]=rot(median(obj_tmp,fwhm),-paral[i],1.0,cx,cy,cubic=-0.5,/pivot)

res_tmp=res_tmp+res[*,*,i]/nobj
wset,0
tvscl, congrid(res_tmp,500,500)
endfor
res_fin=median(res,dim=3);/median(normal)

SXADDPAR, Headers_init,'PARANG', 0
writefits, procdir+'img_'+tg_name_dyn+'.fits',res_fin, Headers_init
writefits, procdir+'img_'+tg_name_dyn+'_median.fits',img_med
writefits, procdir+'img_'+tg_name_dyn+'_dc_minus_median.fits',res, Headers_init

delvar, res_fin, res_tmp, res, img_dc, img_med, obj_tmp

END