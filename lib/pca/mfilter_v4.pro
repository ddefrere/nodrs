;Image filter routine
;D. Mawet, February 2012
;The cutoff inputs must be in pixels

PRO mfilter_v4, highpass=highpass, lowpass=lowpass, bandpass_low=bandpass_low, bandpass_high=bandpass_high, median_f=median_f,$ 
mean_f=mean_f, gauss_f=gauss_f, fourier_ideal=fourier_ideal, fourier_gauss=fourier_gauss, fourier_butterworth=fourier_butterworth

common common_name, tg_name_dyn, tg_name_ori, tg_name_bas, datadir, procdir


print, '*********************************'
print, 'Performing high-pass filtering'
print, '*********************************'

objtmp=readfits(procdir+'img_'+tg_name_dyn+'_dc.fits',header)
tg_name_dyn=tg_name_dyn+'_filt'

for i=0,(size(objtmp))[3]-1 do begin
im=objtmp[*,*,i]

if keyword_set(median_f) then begin
print, 'Doing median filtering for frame #', i
	if keyword_set(highpass) then im=im-filter_image( im, MEDIAN = highpass)
	if keyword_set(lowpass) then im=filter_image( im, MEDIAN = lowpass)
endif

if keyword_set(mean_f) then begin
print, 'Doing smooth filtering for frame #', i
	if keyword_set(highpass) then im=im-filter_image( im, SMOOTH = highpass)
	if keyword_set(lowpass) then im=filter_image( im, SMOOTH = lowpass)
endif

if keyword_set(gauss_f) then begin
print, 'Doing Gaussian filtering for frame #', i
	if keyword_set(highpass) then im=im-filter_image( im, FWHM = highpass)
	if keyword_set(lowpass) then im=filter_image( im, FWHM = lowpass)
endif

if keyword_set(fourier_ideal) then begin
print, 'Doing pure Fourier filtering for frame #', i
	if keyword_set(highpass) then im=bandpass_filter(im, 1d/highpass, 1d, /ideal)
	if keyword_set(lowpass) then im=bandpass_filter(im, 0d, 1d/lowpass, /ideal)
	if keyword_set(bandpass_low) AND keyword_set(bandpass_high) then $
	im=bandpass_filter(im, bandpass_low=bandpass_low, bandpass_high=bandpass_high, /ideal)	
endif

if keyword_set(fourier_gauss) then begin
print, 'Doing Gaussian Fourier filtering for frame #', i
	if keyword_set(highpass) then im=bandpass_filter(im, 1d/highpass, 1d, /gauss)
	if keyword_set(lowpass) then im=bandpass_filter(im, 0d, 1d/lowpass, /gauss)
	if keyword_set(bandpass_low) AND keyword_set(bandpass_high) then $
	im=bandpass_filter(im, bandpass_low=bandpass_low, bandpass_high=bandpass_high, /gauss)	
endif

if keyword_set(fourier_butterworth) then begin
print, 'Doing Butterworth (n=1) Fourier filtering for frame #', i
	if keyword_set(highpass) then im=bandpass_filter(im, 1d/highpass, 1d, butterworth=1)
	if keyword_set(lowpass) then im=bandpass_filter(im, 0d, 1d/lowpass, butterworth=1)
	if keyword_set(bandpass_low) AND keyword_set(bandpass_high) then $
	im=bandpass_filter(im, bandpass_low=bandpass_low, bandpass_high=bandpass_high, butterworth=1)	
endif

objtmp[*,*,i]=im

wset, 1
tvscl, congrid(im,512,512)
endfor

writefits, procdir+'img_'+tg_name_dyn+'_dc.fits', objtmp, header

END

