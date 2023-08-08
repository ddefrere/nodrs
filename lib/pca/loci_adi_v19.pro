;------------------------------------------------------------------
;Locally Optimized Combination of Images (LOCI) - General Procedure
;------------------------------------------------------------------
;
;Specific LOCI input parameters
;------------------------------
;delta   = 1.0	; Delta in resolution element for ADI reference image exclusion
		; FM: Frames within these +/- deg of rotation will not be used
		; in correlation. excluding frames close in time
		; s.t. the source/planet signal is not there to be exluded.
		; Increase if field rotation is slow.  Typical 0.5 to 1 ok, used 1.5.
		; Tradeoff: Speckles are better subtracted from close-in-time frames
		; Do < 1.5 if source in noisy area.
;truncate_init=32; truncation of the reference image pool sorted by correlation
		; FM: # of frames used for correlation, selection from ref. frames
		; correlated zone by zone defines best inversion of...
;dc_sz   = 350	; Output image size, 300, 500~1.5", 1000~3"? needs to be smaller than (dim-cx)*1 or (dim-cy)*2
;Na      = 150	; number of resolution element inside the optimization zone
		; too small can kill companion signal
;g       = 1.0   ; geometry of the optimization zone (see Lafreniere et al. 2007)
		; 0.75* = elongated along theta, squashed radially
		; 1     = sqr
		; 2     = elongated radially, 2xlong
;rin_init=6    ; inner radius where LOCI starts, in resolution elements (fwhm) 
;step_init=0.5  	; radial step size for the section s in resel (fwhm)
;compression_factor = 1 ; rebin of image by this factor (speed), choose between 1 (no compression) or 2 (4 times as fast)
;fake_companion=6e-5 ; level of injected fake companions
;
;
;Keywords
;--------
;	 
;	 	;/filter: apply high-pass filter to the data, with a cutoff at 4*fwhm
;		;/nnls = non-negative least square inversion, which forces the coefficients
;		; 	to be all positive so it always subtracts.
;		;	==> SVD is the default 
;		;/mask_zone_s = masks the s zone in the coefficient computation, to prevent companion flux losses
;		;/display, plot coefficients and display zone o
;		;/contrast_curve, generate contrast curve plot at the end 
;		;fake_companion=X, inject fake companions at level X
;
;
;Written output files
;--------------------
;written by loci_adi_v18 as follows  
;writefits, 'img_'+tg_name_dyn+'_loci_dc.fits', resb, header
;writefits, 'img_'+tg_name_dyn+'_loci_median.fits', resm, header
;writefits, 'img_'+tg_name_dyn+'_loci_mean.fits', rest, header
;
;Example
;-------
;
;status = loci_adi_v16_e(tg_name, cx, cy, dc_sz, plsc, fwhm, rin_init,$
;	 	step_init, truncate_init, Na, g, delta, compression_factor, /H, fake_companion=fake_companion, 		     
;	 	contrast_curve=contrast_curve, /filter)
;
;V1.8: March 2012, D. Mawet (ESO-JPL)
;---------------------------------------
function loci_adi_v19, cx_init, cy_init, dim, dim_init, plsc, fwhm, rin_init, step_init, truncate_inittt, Natt, gtt, deltatt, right_handed, mask_zone_s=mask_zone_s, display=display, NNLS=NNLS, bvls=bvls, svd=svd, damped=damped, normal=normal, verbose=verbose

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir

;Rounding and formatting LOCI input params
;-----------------------------------------
truncate_init=float(round(truncate_inittt))
Na=float(Natt)
g=float(gtt)
delta=float(deltatt)

if keyword_set(verbose) then begin
print, '**************************************************'
print, 'Performing Locally Optimized combination of Images'
print, '**************************************************'
print, 'Double-check LOCI params (# refs imgs, Na, g, detla): ', truncate_init, Na, g, delta
endif

;Reading input file
obj=readfits(procdir+'img_'+tg_name_dyn+'_dc.fits',header)
nobj=(size(obj))[3]

;Parallactic angle loading and unwrapping
;----------------------------------------
paral=readfits(procdir+'vec_'+tg_name_ori+'_paral.fits'); reading parallactic angle

wset, 0
plot, paral, xtitle='Frame #', ytitle='Parallactic Angle'
;read, dummy, PROMPT='Double check parallactic angle variation, hit "1" if OK, "2" if not, then enter: '
;if dummy eq 2 then begin
;print, 'Unwrapping PARANG !'
;paral=(UNWRAP_PHASE(paral/180*!dpi)/!dpi*180); unwrapping PA (!!! needs double checking TBD)
;wset, 0
;plot, paral, xtitle='Frame #', ytitle='Parallactic Angle'
;print, 'Check Again !!!'
;endif

;Loading photometry
;------------------
if normal eq 1 then begin
normal_v=readfits(procdir+'vec_'+tg_name_ori+'_photometry.fits') 
endif else begin
normal_v=fltarr(nobj)
normal_v[*]=normal
endelse

;Initializing LOCI geometrical parameters
;----------------------------------------
m=1
gain=1d
dr=floor (sqrt(!dpi*g*Na*fwhm^2/4)); radial "thickness" of the optimization region
rin=rin_init*fwhm ; inner radius of the LOCI
rout=dim/2-1 - dr ; outer radius max (limited by box O size)
step_law=1.3 ; exponent of the step size increase law 
nstep=floor((rout-rin)^(1d/step_law)/(step_init*fwhm)) ; total number of radial steps
dphimax=(g/2.+(2d*rout/fwhm)*sqrt(float(g)/(!dpi*Na)))^(-1d)
dphimax=(dphimax/!dpi*180)
fc=360d
nqmax=floor(fc/dphimax)
if (nqmax mod 2) EQ 0 then begin
	nqmax=nqmax+1
endif
Na_fwhm=(Na+Na)*fwhm^2 ;+ 100

;Initializing arrays
;-------------------
y = replicate(1,dim) # findgen(dim) - (dim/2) 
y = y / (dim/4)
y=y*2.*!dpi
x=transpose(y)
thetamag=atan(temporary(y),temporary(x))/!dpi*180 + 180
rho=shift(dist(dim),dim/2,dim/2)
nq=findgen(dim/2)
resb=fltarr(dim,dim,nobj)
supref=fltarr(dim,dim)
objr=fltarr(dim,dim)
in=fltarr(nstep)
step=fltarr(nstep)
delta_deg=fltarr(nstep)
ind_an_qt_s=fltarr(Na_fwhm,nstep,nqmax)
if keyword_set(mask_zone_s) then ind_an_qt_o_s=fltarr(Na_fwhm,nstep,nqmax)$
else ind_an_qt_o=fltarr(Na_fwhm,nstep,nqmax)

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
for i=0,nstep-1 do begin
	in[i] = rin + floor((i*step_init*fwhm)^step_law);-dr/2+1 ; inside limit of the zone s
	step[i] = (rin + floor(((i+1)*step_init*fwhm)^step_law)) - in[i] ; radial thickness of the zone s (varying)
	;Delta_deg for reference selection	
	delta_deg[i]= delta / (in[i]/fwhm) / !dpi*180
	
	if keyword_set(verbose) then begin
	print, 'Building boxes for Zones S and O, cutting and saving annuli of size (pix, fwhm): ', strcompress(step[i]),$
	 strcompress(step[i]/fwhm), ', at radius (pix, arcsec): ', strcompress(in[i]), ', ', strcompress(in[i]*plsc)
	endif
	
	;Camemberts
	dphi=(g/2+(2*in[i]/fwhm)*sqrt(g/(!dpi*Na)))^(-1d)
	dphi=(dphi/!dpi*180)
	nq[i]=floor(fc/dphi)
		if (nq[i] mod 2) EQ 0 then begin
		nq[i]=nq[i]+1
		endif

	obj_s_tmp=fltarr(Na_fwhm,nq[i],nobj)
	if keyword_set(mask_zone_s) then obj_o_s_tmp=fltarr(Na_fwhm,nq[i],nobj)$
	else obj_o_tmp=fltarr(Na_fwhm,nq[i],nobj)

	for q=0,nq[i]-1 do begin
	;Zone S
	out=in[i]+step[i]
	;stop
	ind_an_qt_s[0,i,q]=where( (thetamag GE (q*fc/nq[i])) AND (thetamag LE ((q+1)*fc/nq[i])) AND $
	(rho lt out) AND (rho ge in[i]))
	if keyword_set(mask_zone_s) then begin
	;Zone O - ZONE S
	out=in[i]+dr	
	ind_an_qt_o_s[0,i,q]=where( (thetamag GE (q*fc/nq[i])) AND (thetamag LE ((q+1)*fc/nq[i])) AND $
	(rho lt out) AND (rho ge in[i]+step[i]))
	endif else begin
	;Zone O
	out=in[i]+dr	
	ind_an_qt_o[0,i,q]=where( (thetamag GE (q*fc/nq[i])) AND (thetamag LE ((q+1)*fc/nq[i])) AND $
	(rho lt out) AND (rho ge in[i]))
	endelse
		for im=0, nobj-1 do begin
		o=obj[*,*,im]
		obj_s_tmp[*,q,im]=o[ind_an_qt_s[*,i,q]]
		if keyword_set(mask_zone_s) then obj_o_s_tmp[*,q,im]=o[ind_an_qt_o_s[*,i,q]]$
		else obj_o_tmp[*,q,im]=o[ind_an_qt_o[*,i,q]]
		endfor
	endfor
	writefits,strcompress(path) + '/tmpb/' + 'obj_'+tg_name_dyn+'_tr_s_'+strcompress(i,/remove_all)+'.fits', obj_s_tmp
	if keyword_set(mask_zone_s) then begin
	writefits,strcompress(path) + '/tmpb/' + 'obj_'+tg_name_dyn+'_tr_o_s_'+strcompress(i,/remove_all)+'.fits', obj_o_s_tmp
	endif else begin
    writefits,strcompress(path) + '/tmpb/' + 'obj_'+tg_name_dyn+'_tr_o_'+strcompress(i,/remove_all)+'.fits', obj_o_tmp
    endelse
	obj_o_tmp=0
	obj_s_tmp=0
	obj_o_s_tmp=0
endfor


;Main LOCI 
;---------
;Loop on "target images"
;***********************
for im=0,nobj-1 do begin

if keyword_set(verbose) then begin
print, '----------------------------------------------------------'
print, 'Doing image #: ', strcompress(im), ', at parallactic angle: ', strcompress(paral[im])
print, '----------------------------------------------------------'
endif

;Reinit variables
objr=fltarr(dim,dim)
supref=fltarr(dim,dim)

;Loop on annulii
;***************
for i=0,nstep-1 do begin
obj_s_tmp=(readfits(strcompress(path) + '/tmpb/' + 'obj_'+tg_name_dyn+'_tr_s_'+strcompress(i,/remove_all)+'.fits',/silent))[*,*,*]	
if keyword_set(mask_zone_s) then $
 obj_o_s_tmp=(readfits(strcompress(path) + '/tmpb/' + 'obj_'+tg_name_dyn+'_tr_o_s_'+strcompress(i,/remove_all)+'.fits',/silent))[*,*,*]$
 else obj_o_tmp=(readfits(strcompress(path) + '/tmpb/' + 'obj_'+tg_name_dyn+'_tr_o_'+strcompress(i,/remove_all)+'.fits',/silent))[*,*,*]
	
index=where( (abs(paral) lt (abs(paral[im])-delta_deg[i])) or (abs(paral) gt (abs(paral[im])+delta_deg[i])) )
nref=n_elements(index)

if nref eq 1 OR index[0] eq -1 then begin
	delta_deg_tmp=delta_deg
	if keyword_set(verbose) then print, '!!!Exclusion law too severe:reducing cone opening automatically!!!'
	while nref eq 1 OR index[0] eq -1 do begin
	delta_deg_tmp=delta_deg_tmp/2
	
	index=where( (abs(paral) lt (abs(paral[im])-delta_deg_tmp[i])) or (abs(paral) gt (abs(paral[im])+delta_deg_tmp[i])) )
	nref=n_elements(index)
	if keyword_set(verbose) then print, 'Exclusion cone opening=', strcompress(2*delta_deg[i]), ', # reference images left in the pool:', strcompress(nref)
	endwhile
endif

if keyword_set(verbose) then begin
print, 'Radius (pix, arcsec):', strcompress(in[i]), ', ',strcompress(in[i]*plsc), ', exclusion cone opening=', strcompress(2*delta_deg[i]), ', # reference images left in the pool:', strcompress(nref)
endif
 
;Reinit variables
if truncate_init gt nref then begin
	truncate_svd = nref
	A=fltarr(truncate_svd,truncate_svd)
	b=fltarr(truncate_svd)
endif else begin
	truncate_svd = truncate_init
	A=fltarr(truncate_svd,truncate_svd)
	b=fltarr(truncate_svd)
endelse
correlation_b_r=fltarr(nref)
b_o_t=fltarr(nref)
b_s_t=fltarr(nref)


	;Loop on camemberts
	;******************
	for q=0,nq[i]-1 do begin
			
	;Build Vector b
	;--------------
	if keyword_set(mask_zone_s) then begin 
	b_obj_o=obj_o_s_tmp[*,q,im] 
	b_obj_s=obj_s_tmp[*,q,im] 
	endif else begin 
	b_obj_o=obj_o_tmp[*,q,im]
	b_obj_s=obj_s_tmp[*,q,im] 
	endelse
	
	for j=0,nref-1 do begin	
	if keyword_set(mask_zone_s) then $
	b_ref_o=obj_o_s_tmp[*,q,index[j]] else $ 
	b_ref_o=obj_o_tmp[*,q,index[j]]
	b_ref_s=obj_s_tmp[*,q,index[j]] 

	;Selection based on Correlation
	b_o_t[j]=total(b_obj_o*b_ref_o)
	b_s_t[j]=total(b_obj_s*b_ref_s)
	endfor
	
	index_svd=(reverse(sort(b_o_t)))[0:truncate_svd-1]
	b_o=b_o_t[index_svd]
	b_s=b_s_t[index_svd]
	b_obj_2=total(b_obj_o*b_obj_o)

	;Build Matrix A
	;--------------
	for j=0,truncate_svd-1 do begin
		if keyword_set(mask_zone_s) then $ 
		rj=obj_o_s_tmp[*,q,index(index_svd[j])] else $ 
		rj=obj_o_tmp[*,q,index(index_svd[j])]
			for k=0,truncate_svd-1 do begin
				if keyword_set(mask_zone_s) then $ 
				rk=obj_o_s_tmp[*,q,index(index_svd[k])] else $ 
				rk=obj_o_tmp[*,q,index(index_svd[k])]
			A[k,j]=total(rj*rk) 
			endfor
	endfor
	
	;Coefficients - Least square minimization
	;----------------------------------------
	if keyword_set(NNLS) OR keyword_set(damped) then begin
		;NNLS: non-negative least squares
		;-------------------------------
		m=truncate_svd
		n=m
		indx=findgen(n)
		indx[*]=0d
		w=findgen(n)
		w[*]=1d
		sol=fltarr(truncate_svd)
		if keyword_set(damped) then begin
		;Still need to nest a non-linear 1D min of J1 with damped parameter
		;Currently assuming uniform lagrange (damped) multiplier
		b=b_o-damped/2*b_s;introducing the damping here 
		A_s=A;saving the matrix A because nnls reorders it
		b_s=b;saving the vector because nnls reorders it
		nnls, A, m, n, b, sol, rnorm, w, indx, mode	
		J1=transpose(sol)#A_s#sol-2*b_o##transpose(sol)+b_obj_2+damped*b_s##transpose(sol)
		if keyword_set(verbose) then begin
		print, 'Cost function', J1, 'Lagrange multiplier', damped
		endif
		endif else begin
		b=b_o
		A_s=A
		b_s=b
		nnls, A, m, n, b, sol, rnorm, w, indx, mode	
		endelse
	endif 
	
	if keyword_set(SVD) then begin
		;SVD: singular value decomposition
		;---------------------------------
		if truncate_svd eq 1 then sol=1 else begin
		LA_SVD, A, W, U, V
		sol=SVSOL(U,W,V,b_o)
		;reg=1e6
		;Ap=[A,reg*identity((size(A))[2])]
		;bp=[b_o,fltarr((size(A))[2])]
		;LA_SVD, Ap, W, U, V
		;SVDC, Ap, W, U, V
		;sol=SVSOL(U,W,V,b_o,/double)
		;plot, sol
		;stop
		endelse
	endif
	
	if keyword_set(bvls) then begin
		;BVLS
		BND=fltarr(2,(size(b_o))[1])
		BND[0,*]=-1d
		BND[1,*]=1d
	   	BVLS, A, b_o, BND, sol, $
        EPS=eps, /FASTNORM, IERR=ierr, INDEX=indexbvls, ITER=iter, $
        ITMAX=itmax, NSETP=nsetp, RNORM=rnorm, W=w
	endif
	
	;Super PSF construction, and subtraction from target frame
	;---------------------------------------------------------
		;t=systime(1)
		if keyword_set(display) then begin
		if truncate_svd ne 1 then begin
		wset, 1
		plot, sol
		endif
		endif
	
	;Master PSF	
	for s=0,truncate_svd-1 do begin 
		r=fltarr(dim,dim)
		r[ind_an_qt_s[*,i,q]]=obj_s_tmp[*,q,index(index_svd[s])]
		supref[ind_an_qt_s[*,i,q]] = supref[ind_an_qt_s[*,i,q]]+ sol[s] * r[ind_an_qt_s[*,i,q]]
	endfor
	
	;Master PSF subtraction from current object frame
	o=fltarr(dim,dim)
	o[ind_an_qt_s[*,i,q]]=obj_s_tmp[*,q,im]
	objr[ind_an_qt_s[*,i,q]]=objr[ind_an_qt_s[*,i,q]]+(o[ind_an_qt_s[*,i,q]]-gain*supref[ind_an_qt_s[*,i,q]])
			if keyword_set(display) then begin
			o_d=fltarr(dim,dim)
			objr_d=fltarr(dim,dim)
			o_d[ind_an_qt_o[*,i,q]]=obj_o_tmp[*,q,im]
			objr_d[ind_an_qt_o[*,i,q]]=objr_d[ind_an_qt_o[*,i,q]]+(o_d[ind_an_qt_o[*,i,q]]-gain*supref[ind_an_qt_o[*,i,q]])
			wset, 3
			tvscl, congrid(objr+objr_d,500,500)>(-1000)<10000
			endif		
	endfor		
wset, 0
tvscl, congrid(objr,500,500)>(-1000)<10000
endfor
	o=obj[*,*,im]
	noise_init=stddev(o(where(abs(o) le 2d^16)))/normal_v[im]
	if keyword_set(verbose) then begin
	print, '*************************************'
	print, 'Initial image noise: ', noise_init
	endif
	noise_final=stddev(objr(where(abs(objr) le 2d^16)))/normal_v[im]
	if keyword_set(verbose) then begin
	print, 'Final image noise: ', noise_final
	print, 'Noise attenuation factor: ', noise_init/noise_final
	print, '*************************************'
	print, ''
	print, 'Derotating and stacking image #: ', strcompress(im)
	print, ''
	endif
	if keyword_set(right_handed) then paral[im]=-paral[im]
	resb[*,*,im]=rot(objr,-paral[im],1.0,dim/2,dim/2,cubic=-0.5,/pivot)
if nobj ne 1 then begin
wset, 2
tvscl, congrid(smooth(total(resb,3)/(im+1),3),500,500)>(-1000)<10000
endif
endfor

;Updating parang keyword in header
SXADDPAR, header,'PARANG', 0

;Writing output files
;--------------------
if keyword_set(verbose) then begin
print, '**************************************************************'
print, 'Saving output images under ', tg_name_dyn, '_loci_median/mean.fits'
print, '**************************************************************'
endif
resm=median(resb,dim=3)/median(normal_v)

tg_name_dyn=tg_name_dyn+'_loci'
tg_name_bas=tg_name_ori+'_loci'


;Outputs
if keyword_set(wcs) then begin
resm_fin=fltarr(dim_init,dim_init)
resm_fin[dim_init/2-dim/2:dim_init/2+dim/2-1,dim_init/2-dim/2:dim_init/2+dim/2-1]=resm
resm_fin=shift(resm_fin,cx_init-dim_init/2,cy_init-dim_init/2)
writefits, procdir+'img_'+tg_name_dyn+'_median.fits',resm_fin, header
endif else begin
writefits, procdir+'img_'+tg_name_dyn+'_median.fits',resm, header
endelse


;writefits, procdir+'proj_'+tg_name_dyn+'_dc.fits', proj_obj_dc
;writefits, procdir+tg_name_dyn+'_dc.fits', objt
;writefits, procdir+'supref_dc.fits', supreft
;writefits, procdir+'obj_minus_mean.fits', obj


if keyword_set(verbose) then begin
print, '********************************************'
print, 'Deleting temporary directory, LOCI completed'
print, '********************************************'
endif

spawn, 'rm -d -r tmpb'

return, 1

END
