;Basic cosmetic data reduction for Keck
;D. Mawet, 2010-2012
;V1-V4, 2010-2011, DM, Adapted from Palomar WCS basic data reduction
;V5, OCT-2012, DM: cleanup, retro-compatibility to IDL 7.0-8.0
;V6, 14-SEP-2013, DM: REFINE peak photometry, updated PARANG COMPUTATION ACCORDING TO NIRC2 WEBPAGE
;    http://www2.keck.hawaii.edu/inst/nirc2/nirc2_ao.html

PRO make_basic_nirc2_v6, cx, cy, st_ob, nd_ob, pfmt, dim, fwhm, nirc2_filter, t_init, n_init, Header_start, flat, bpix

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

;Basic treatment
;---------------
;Specific input parameters for basic treatment
nfr=nd_ob-st_ob+1; calculated number of frames
THSZINIT = 12	; subwindow size for recentering

;Init array (datacube and parallactic angle vector)
img_dc 	   	= fltarr(dim,dim,nfr)
paral_dc_parang	= fltarr(nfr)	;? tracking parallactic ang?

;Subraster mode
;--------------
if dim eq 512 then begin
  dim  = 1024
  flat = flat[dim/2-dim/4:dim/2+dim/4-1,dim/2-dim/4:dim/2+dim/4-1]
  bpix = bpix[dim/2-dim/4:dim/2+dim/4-1,dim/2-dim/4:dim/2+dim/4-1]
  dim  = 512
endif

if dim eq 256 then begin
  dim  = 1024
  flat = flat[dim/2-dim/8:dim/2+dim/8-1,dim/2-dim/8:dim/2+dim/8-1]
  bpix = bpix[dim/2-dim/8:dim/2+dim/8-1,dim/2-dim/8:dim/2+dim/8-1]
  dim  = 256
endif

;Main loop
;---------
;Sky frame; sky frame obsolete for LOCI

;Reference frame
Headers_init = headfits(datadir+STRING(st_ob,f=pfmt))
paral_init   = SXPAR(Headers_init,'PARANG')

;Target/object frame
flno = FIX([st_ob, nd_ob])
frno = flno[0] + LINDGEN(flno[1]-flno[0]+1)
normal=fltarr(nfr)
index_dc=0

for i=0, nfr-1 do begin
   print, 'Doing frame #: ', frno[i]
   ;File reading - Headers reading for parallactic angle calculation
   objtmp=readfits(datadir+STRING(frno[i],f=pfmt))
   Headers=headfits(datadir+STRING(frno[i],f=pfmt))
   t = SXPAR(Headers,'ITIME')
   n = SXPAR(Headers,'COADDS')
   tg_name_check=strcompress(SXPAR(Header_start,'OBJECT'),/remove_all)
   nirc2_filter_check=strcompress(SXPAR(Header_start,'FWINAME'),/remove_all)

	;Fool-proof Checks
	;-----------------
	if tg_name_check ne tg_name_bas then begin
	print, '!!! Warning !!! Current frame not nominal target !!! Stopping !!!'
	stop
	endif
	if nirc2_filter_check ne nirc2_filter then begin
	print, '!!! Warning !!! Current filter not nominal !!! Stopping !!!'
	stop
	endif
	if t ne t_init then begin
	print, 'Warning: integration time of current frame not nominal: ', strcompress(t)
	read, dummy, PROMPT='Do you want to keep the current frame ? Press "1" if yes, "2" if no, then enter '
	if dummy eq 2 then begin
	index_dc=index_dc
	GOTO, SKIP
	endif
	endif
	if n ne n_init then begin
	print, 'Warning: # of coadds of current frame not nominal: ', strcompress(n)
	read, dummy, PROMPT='Do you want to keep the current frame ? Press "1" if yes, "2" if no, then enter '
	if dummy eq 2 then begin
	index_dc=index_dc
	GOTO, SKIP
	endif
	endif

   ;Flat-fielding, sky subtraction
   objtmp=(objtmp)/flat

   ;bad pixel - cosmic ray removal
   objtmp = BPIXFIX(objtmp, bpix)
   objtmp=sigma_filter(objtmp,5,N_SIGMA=3,/monitor)

   ;Recentering!!!!!!!!!!!!!!!!!!!!
   if i eq 0 then begin
   frame_initt = objtmp
   endif
   recenter_frames_v6, objtmp, frame_initt, cx, cy, fwhm, 0, /use_cursor   
   ;Photometric vector
   if nirc2_filter eq 'H' then begin
   normal[index_dc]=max(filter_image(objtmp[cx-1.5*fwhm:cx+1.5*fwhm,cy-1.5*fwhm:cy+1.5*fwhm],fwhm=fwhm))/0.0012
   print, 'Peak: ',normal[index_dc]
   endif
   if nirc2_filter eq 'Kp' then begin
   normal[index_dc]=max(filter_image(objtmp[cx-1.5*fwhm:cx+1.5*fwhm,cy-1.5*fwhm:cy+1.5*fwhm],fwhm=fwhm))/0.0022
   print, 'Peak: ',normal[index_dc]
   endif
   if nirc2_filter eq 'Ks' then begin
   normal[index_dc]=max(filter_image(objtmp[cx-1.5*fwhm:cx+1.5*fwhm,cy-1.5*fwhm:cy+1.5*fwhm],fwhm=fwhm))/0.0022
   print, 'Peak: ',normal[index_dc]
   endif
   if nirc2_filter eq 'K' then begin
   normal[index_dc]=max(filter_image(objtmp[cx-1.5*fwhm:cx+1.5*fwhm,cy-1.5*fwhm:cy+1.5*fwhm],fwhm=fwhm))/0.0022
   print, 'Peak: ',normal[index_dc]
   endif
   wset, 0
   tvscl, congrid(objtmp,500,500)>(-10000)<10000
   img_dc[*,*,index_dc]=objtmp
   ;nirc2PA = PARANG + ROTPOSN - INSTANGL; ONLY IN VERTANG MODE
   paral_dc_parang(index_dc)=SXPAR(Headers,'PARANG')+SXPAR(Headers,'ROTPOSN')-SXPAR(Headers,'INSTANGL')
   index_dc=index_dc+1
   SKIP:
endfor

img_dc=img_dc[*,*,0:index_dc-1]
paral_dc_parang=paral_dc_parang[0:index_dc-1]
normal=normal[0:index_dc-1]

print, 'Total Parallactic angle variation', (max(paral_dc_parang)-min(paral_dc_parang))
print, 'Total integration time', n*t*index_dc
SXADDPAR, Headers_init,'PARALTOT', abs(max(paral_dc_parang)-min(paral_dc_parang))
SXADDPAR, Headers_init,'TINTTOT', n*t*index_dc
 
;Write output file
;-----------------
writefits, procdir+'img_'+tg_name_bas+'_dc.fits',img_dc, Headers_init
writefits, procdir+'vec_'+tg_name_bas+'_paral.fits', paral_dc_parang;-paral_init
writefits, procdir+'vec_'+tg_name_bas+'_photometry.fits', normal;-paral_init
img_dc=0
objtmp=0

END
