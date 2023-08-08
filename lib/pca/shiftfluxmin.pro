FUNCTION shiftfluxmin, arg

COMMON frames, frame_init_com, obj_com, mask_com

shx_init=arg[0]
shy_init=arg[1]
;flux=arg[2]
flux=1
sz=(size(frame_init_com))[1]
res=mask_com*(frame_init_com-flux*shifti(obj_com,shx_init,shy_init))

buffer=2;*ceil(max(abs([shx_init,shy_init])))
;if buffer lt 1 then buffer=1
res_buffer=res[buffer:sz-1-buffer, buffer:sz-1-buffer]
sz_buffer=(size(res_buffer))[1]

wset,2 
tvscl, congrid((res_buffer),500,500)>(-1000)<1000

;figureofmerit=total(abs(res_buffer)^2)/(sz_buffer)^2
figureofmerit=stddev(res_buffer)

;print, 'x,y, flux: ', shx_init, shy_init, flux 
print, 'Figure of Merit: ', figureofmerit
return, figureofmerit

END