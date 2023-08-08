;+
; NAME: HOUR_ANGLE
; 
; DESCRIPTION
;   Computes the hour angle for the observation of a given star at a given time
;   Two different techniques are used and checked wrt each other
;
; INPUTS
;   mjd:     : the Modified Julian Date (JD - 2,400,000.5) -- can be scalar of array
;   ra:      : the right ascension of the object in hours (float, e.g. for Vega: ra = 18+36/60D0+56.3364/3600D0)
;   dec:     : the declination of the object in degrees (float, e.g. for Vega: dec = 38+47/60D0+01.291/3600D0)
;
; KEYWORDS
;   lat      : the latitude of the observing site in degrees (default is Mount Graham, AZ)
;   lon      : the longitude of the observing site in degrees (default is Mount Graham, AZ)
;   altitude : the altitude of the observing site in meters (default is Mount Graham, AZ)
;   star_alt : on output, gives the stellar altitude at the time of observation [degrees]
;   star_az  : on output, gives the stellar azimuth (measured East from North) at the time of observation [degrees]
;   info     : set this parameter to print the two different estimations on the screen
;
; OUTPUT
;   Returns the hour angle in degrees (scalar or array)
;
; MODIFICATION HISTORY
;   Version 1.0, 08-AUG-2006, by Olivier Absil, University of Liège,  absil@astro.ulg.ac.be
;   Version 1.1, 09-APR-2008, OA: stellar altitude and azimuth can now be retrieved through keywords
;   Version 1.2, 10-NOV-2009, DD: added output information
;   Version 1.3, 12-DEC-2012, DD: modified default site to Mount Graham

FUNCTION HOUR_ANGLE, mjd, ra, dec, LAT=lat, LON=lon, ALTITUDE=altitude, STAR_ALT=alt, STAR_AZ=az, INFO=info

; Check input parameters
IF NOT KEYWORD_SET(lat) THEN lat = TEN(32,42,04.69)
IF NOT KEYWORD_SET(lon) THEN lon = -TEN(109,53,21.09)
IF NOT KEYWORD_SET(altitude) THEN altitude = 3267D0
IF NOT KEYWORD_SET(info) THEN info = 0

IF (SIZE(mjd))[0] GT 0 THEN mjd = REFORM(mjd)
n_mjd = N_ELEMENTS(mjd)
IF n_mjd GT 1 THEN BEGIN
  ra2 = CONGRID([ra], n_mjd)
  dec2 = CONGRID([dec], n_mjd)
ENDIF ELSE BEGIN
  ra2 = ra
  dec2 = dec
ENDELSE

; First method with Astrolib
EQ2HOR, ra2*360D0/24D0, dec2, mjd + 2400000.5D0, alt, az, ha, LAT=lat, LON=lon, ALTITUDE=altitude
ha = ha - (ha GT 180)*360
IF info THEN PRINT, 'Hour angle: ', ha

; Second method with the sidereal time
CT2LST, lst, lon, 0, mjd + 2400000.5D0
;PRINT, 'Local sidereal time: ', lst
ha_2 = (lst - ra)*360D0/24D0 ; in degrees
ha_2 = ha_2 - (ha_2 GT 180)*360 + (ha_2 LT -180)*360
IF info THEN PRINT, 'Hour angle (check): ', ha_2

IF TOTAL(ha - ha_2) GT 1.0 THEN MESSAGE, 'Problem with hour angle computation.'

; Antoine's method
;DAYCNV, 2400000.5D0 + MJD, yr, mn, day, UT
;PRINT, yr, mn, day, UT
;lst = 100.46 + 0.985647 * (MJD-51544.5D0) + lon + 15 * UT
;PRINT, 'Local sidereal time: ', lst
;ha_3 = lst - ra
;PRINT, 'Hour angle: ', (ha_3*360D0/24D0) MOD 360

RETURN, ha
END