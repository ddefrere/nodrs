;+
; NAME: QUERY_TGT
; 
; DESCRIPTION:
;   This function gives on output a structure with the following information on the output target
;             - name     : name of the target
;             - sp_type  : spectral type
;             - spectrum : spectrum file of the target (Pickles library)
;             - ra       : target right ascension
;             - dec      : target declination
;             - dist     : target distance in pc
;             - vmag     : V-band magnitude
;             - e_vmag   : error on V-band magnitude
;             - hmag     : H-band magnitude
;             - e_hmag   : error on H-band magnitude
;             - kmag     : K-band magnitude
;             - e_kmag   : error on K-band magnitude
;             - lmag     : L-banb magnitude
;             - e_lmag   : error on L_band magnitude
;             - mass     : target mass in solar mass
;             - temp     : target temperature in K
;             - vsini    : target projected rotationnal velocity in km/s (not supported in current version)
;             - loggg    : log g
;             - ldm      : target limb-darkened angular diameter (in mas)
;             - lde      : error on the limb-darkened angular diameter (in mas)
;             - uldh     : H-band limb darkening linear coefficient
;             - uldk     : K-band limb darkening linear coefficient
;             - pa       : position angle of the stellar photosphere (not supported in current version)
;
; INPUTS:
;   input_name:    Name of the target to retrieve
;
; KEYWORD
;   CFA:           Force to start with CFA server
;   
; MODIFICATION HISTORY:
;   Version 1.0,  03-APR-2017, Denis Defrère, ddefrere@ulg.ac.be

FUNCTION QUERY_TGT, input_name, CFA=cfa

   ; 1. QUERY SIMBAD
   QUERYSIMBAD, input_name, ra, de, id, parallax=p, CFA=cfa, /silent, ERRMSG=errmsg
   IF KEYWORD_SET(errmsg) THEN BEGIN
    IF KEYWORD_SET(CFA) THEN cfa = 0 ELSE cfa = 1
    QUERYSIMBAD, input_name, ra, de, id, parallax=p, /silent, CFA=cfa, ERRMSG=errmsg2
    IF KEYWORD_SET(errmsg2) THEN BEGIN
      MESSAGE, 'Object ' + input_name + ' not found. SKIPPED.', /CONTINUE
      RETURN, {NAME: input_name, RA: 99., DEC: 99., SP_TYPE: 'unknown', SPECTRUM: 'notvalid', DIST: 0., VMAG: 0., E_VMAG: 0., VMAG_REF: 'None00', HMAG: 0., E_HMAG: 0., HMAG_REF: 'None00', $
               KMAG: 0., E_KMAG: 0., KMAG_REF: 'None00', LMAG: 0., E_LMAG: 0., LMAG_REF: 'None00', TEFF: 0, TEFF_REF: 'None00', MASS: 0., MASS_REF: 'None00', LDM: 0., $
               LDE: 0., LDM_REF: 'None00', LOGG: 0., LOGG_REF: 'None00', VSINI: 0., VSINI_REF: 'None00', ULDH: 0., ULDH_REF: 'None00', ULDK: 0., ULDK_REF: 'None00', PA: 0., PA_REF: 'None00'}
    ENDIF
   ENDIF
   ra /= 15. ; deg to hour
   IF NOT KEYWORD_SET(p) THEN READ, p, PROMPT='Enter distance in pc for ' + input_name + ' : '
      
   ; 2. QUERY KAHRCHENKO
   data = QUERYVIZIER('I/280B', input_name, /allcolumns, /silent, CFA=cfa)
   ;IF data EQ -1 THEN MESSAGE, 'Target not found!'
   IF N_ELEMENTS(data._r) GT 1 THEN BEGIN
     tmp      = MIN(data._r, idx_ok) 
     data     = data[idx_ok]
   ENDIF
   sp_type  = data.sptype
   vmag     = data.vmag
   e_vmag   = data.e_vmag
   vmag_ref = "Ka2009"
   hmag     = data.hmag
   e_hmag   = data.e_hmag
   hmag_ref = "Ka2009"
   kmag     = data.kmag
   e_kmag   = data.e_kmag
   kmag_ref = "Ka2009"
   lmag     = 0.
   e_lmag   = 0.
   lmag_ref = "Not found"
   
   ; 3. QUERY DUCATI 2002 for H, K and L band magnitudes
   data = QUERYVIZIER('II/237', input_name, CFA=cfa, /allcolumns, /silent) ; Ducati (report K-V!!!)
   IF (SIZE(data))[0] NE 0 THEN BEGIN
     IF N_ELEMENTS(data._r) EQ 1 THEN BEGIN ; If more than one star, one cannot tell which one is the right one!!
       IF FINITE(data.k_v) THEN BEGIN
         IF data.ksig NE 0 THEN BEGIN
           kmag     = data.k_v + vmag
           e_kmag   = data.ksig*0.01  ; Ducati returns 100*error
           kmag_ref = "Du2002"
           k_ok     = 1
         ENDIF
       ENDIF
       IF FINITE(data.h_v) THEN BEGIN
         IF data.hsig NE 0 THEN BEGIN
           hmag     = data.h_v + vmag
           e_hmag   = data.hsig*0.01  ; Ducati returns 100*error
           hmag_ref = "Du2002"
           h_ok     = 1
         ENDIF
       ENDIF
       IF FINITE(data.l_v) THEN BEGIN
         IF data.lsig NE 0 THEN BEGIN
           lmag     = data.l_v + vmag
           e_lmag   = data.lsig*0.01  ; Ducati returns 100*error
           lmag_ref = "Du2002"
           l_ok     = 1
         ENDIF
       ENDIF
     ENDIF
   ENDIF
   
   ; 4. If not found (or error null) look into GEZARI 1999
   IF NOT KEYWORD_SET(k_ok) THEN BEGIN   
     data = QUERYVIZIER('II/225/catalog', input_name, CFA=cfa, /allcolumns, /silent)
     IF (SIZE(data))[0] NE 0 THEN BEGIN
       idx = WHERE(data.lambda EQ 2.2 AND data.x_f_ir EQ 'M', n)
       IF n GT 1 THEN BEGIN
         km_data  = MOMENT(data[idx].f_ir)
         kmag     = km_data[0]
         e_kmag   = SQRT(km_data[1]) > 0.01
         kmag_ref = "Ge1999"
       ENDIF
       idx = WHERE(data.lambda EQ 1.65 AND data.x_f_ir EQ 'M', n)
       IF n GT 1 THEN BEGIN
         hm_data  = MOMENT(data[idx].f_ir)
         hmag     = hm_data[0]
         e_hmag   = SQRT(hm_data[1]) > 0.01
         hmag_ref = "Ge1999"
       ENDIF
     ENDIF 
   ENDIF
    
   ; 5. Compute diameter with surface brightness relations
   diam     = SBR_CHELLI(vmag, e_vmag, kmag, e_kmag, sp_type)
   diam_ref = "SBRChelli"
   
   ; 5. QUERY VALENTI for teff and mass (mass is not important for exozodi, only for binary fitting)
   teff_ref = 'None00'
   logg = 0 & logg_ref = 'None00'
   mass = 0 & mass_ref = 'None00'
   data = QUERYVIZIER('J/ApJS/159/141', input_name, CFA=cfa, /allcolumns, /silent)
   IF (SIZE(data))[0] NE 0 THEN BEGIN
     teff = data.teff  & teff_ref = "Va2005"
     mass = data.mass  & mass_ref = "Va2005"
     logg = data.log_g & logg_ref = "Va2005"
   ENDIF ELSE BEGIN
     ; If not found, query Gray 2006
     data = QUERYVIZIER('J/AJ/132/161', input_name, CFA=cfa, /allcolumns, /silent)
     IF (SIZE(data))[0] NE 0 THEN BEGIN
       teff = data.teff  & teff_ref = "Gr2006"
       logg = data.log_g & logg_ref = "Gr2006"
       mass = 0
     ENDIF ELSE BEGIN
       ; If not found, query Espramer
       data = QUERYVIZIER('J/A+A/398/1121', input_name, CFA=cfa, /allcolumns, /silent)
       IF (SIZE(data))[0] NE 0 THEN BEGIN
         teff = data.teff  & teff_ref = "Er2003"
         logg = data.log_g & logg_ref = "Er2003"
         mass = 0
       ENDIF ELSE BEGIN
         ; If not found, query Soubiran 2011
         data = QUERYVIZIER('J/A+A/530/A138', input_name, CFA=cfa, /allcolumns, /silent)
         IF (SIZE(data))[0] NE 0 THEN BEGIN
           IF N_ELEMENTS(data._r) GT 1 THEN BEGIN
             tmp  = MIN(data._r, idx_ok)
             data = data[idx_ok]
           ENDIF
           IF FINITE(data.teff) THEN BEGIN
             teff = data.teff 
             teff_ref = "So2011"
           ENDIF ELSE teff = 0
           IF FINITE(data.logg) THEN BEGIN
             logg = data.logg
             logg_ref = "So2011"
           ENDIF ELSE teff = 0               
           mass = 0
         ENDIF
       ENDELSE
     ENDELSE
   ENDELSE  
   
   ; 6. Associate correct spectrum file (iteratively look for a match!)
   sp_type = STRLOWCASE(sp_type)
   s       = STRMID(sp_type, 0, 1)   
   vec1    = ['0','1','2','3','4','5','6','7','8','9'] ; Possible values for the spectral number
   vec2    = ['i','ii','iii','iv','v']                 ; Possible values for the luminosity class
   i1      = FIX(STRMID(sp_type, 1, 2))
   i2      = 4
   IF STREGEX(sp_type, 'i')   GE 0 THEN i2 = 0
   IF STREGEX(sp_type, 'ii')  GE 0 THEN i2 = 1
   IF STREGEX(sp_type, 'iii') GE 0 THEN i2 = 2
   IF STREGEX(sp_type, 'v')   GE 0 THEN i2 = 4
   IF STREGEX(sp_type, 'iv')  GE 0 THEN i2 = 3
   spectrum = 'uk' + s + vec1[i1] + vec2[i2] + '.dat'
   
   pickles_files = ['ukb2ii.dat'   ,'ukf2iii.dat'  ,'ukg5iv.dat'   ,'ukk4i.dat'    ,'ukm5iii.dat'  ,'ukrk4iii.dat', 'ukb2iv.dat'   ,'ukf2v.dat'    ,'ukg5v.dat'    ,'ukk4iii.dat'  ,'ukm5v.dat'    ,'ukrk5iii.dat', $
                     'uka0i.dat'    ,'ukb3i.dat'    ,'ukf5i.dat'    ,'ukg8i.dat'    ,'ukk4v.dat'    ,'ukm6iii.dat'  ,'ukwf5v.dat'   ,'uka0iii.dat'  ,'ukb3iii.dat'  ,'ukf5iii.dat'  ,'ukg8iii.dat'  ,'ukk5iii.dat' , $
                     'ukm6v.dat'    ,'ukwf8v.dat'   ,'uka0iv.dat'   ,'ukb3v.dat'    ,'ukf5iv.dat'   ,'ukg8iv.dat'   ,'ukk5v.dat'    ,'ukm7iii.dat'  ,'ukwg0v.dat'   ,'uka0v.dat'    ,'ukb57v.dat'   ,'ukf5v.dat',    $
                     'ukg8v.dat'    ,'ukk7v.dat'    ,'ukm8iii.dat'  ,'ukwg5iii.dat' ,'uka2i.dat'    ,'ukb5i.dat'    ,'ukf6v.dat'    ,'ukk01ii.dat'  ,'ukm0iii.dat'  ,'ukm9iii.dat'  ,'ukwg5v.dat'   ,'uka2v.dat',    $
                     'ukb5ii.dat'   ,'ukf8i.dat'    ,'ukk0iii.dat'  ,'ukm0v.dat'    ,'uko5v.dat'    ,'ukwg8iii.dat' ,'uka3iii.dat'  ,'ukb5iii.dat'  ,'ukf8iv.dat'   ,'ukk0iv.dat'   ,'ukm10iii.dat' ,'uko8iii.dat',  $
                     'ukwk0iii.dat' ,'uka3v.dat'    ,'ukb6iv.dat'   ,'ukf8v.dat'    ,'ukk0v.dat'    ,'ukm1iii.dat'  ,'uko9v.dat'    ,'ukwk1iii.dat' ,'uka47iv.dat'  ,'ukb8i.dat'    ,'ukg0i.dat'    ,'ukk1iii.dat',  $
                     'ukm1v.dat'    ,'ukrf6v.dat'   ,'ukwk2iii.dat' ,'uka5iii.dat'  ,'ukb8v.dat'    ,'ukg0iii.dat'  ,'ukk1iv.dat'   ,'ukm2.5v.dat'  ,'ukrf8v.dat'   ,'ukwk3iii.dat' ,'uka5v.dat'    ,'ukb9iii.dat',  $
                     'ukg0iv.dat'   ,'ukk2i.dat'    ,'ukm2i.dat'    ,'ukrg0v.dat'   ,'ukwk4iii.dat' ,'uka7iii.dat'  ,'ukb9v.dat'    ,'ukg0v.dat'    ,'ukk2iii.dat'  ,'ukm2iii.dat'  ,'ukrg5iii.dat' ,'uka7v.dat',    $
                     'ukf02iv.dat'  ,'ukg2i.dat'    ,'ukk2v.dat'    ,'ukm2v.dat'    ,'ukrg5v.dat'   ,'ukb0i.dat'    ,'ukf0i.dat'    ,'ukg2iv.dat'   ,'ukk34ii.dat'  ,'ukm3ii.dat'   ,'ukrk0iii.dat' ,'ukb0v.dat',    $
                     'ukf0ii.dat'   ,'ukg2v.dat'    ,'ukk3i.dat'    ,'ukm3iii.dat'  ,'ukrk0v.dat'   ,'ukb12iii.dat' ,'ukf0iii.dat'  ,'ukg5i.dat'    ,'ukk3iii.dat'  ,'ukm3v.dat'    ,'ukrk1iii.dat' ,'ukb1i.dat',    $
                     'ukf0v.dat'    ,'ukg5ii.dat'   ,'ukk3iv.dat'   ,'ukm4iii.dat'  ,'ukrk2iii.dat' ,'ukb1v.dat'    ,'ukf2ii.dat'   ,'ukg5iii.dat'  ,'ukk3v.dat'    ,'ukm4v.dat'    ,'ukrk3iii.dat']
  
   c1      = 0 ; Counter 1
   c2      = 0 ; Counter 2
   spec_ok = 0 
   WHILE spec_ok EQ 0 DO BEGIN
     idx_ok = WHERE(spectrum EQ pickles_files, n_ok)
     IF n_ok EQ 0 THEN BEGIN
       c2 += 1
       spectrum = 'uk' + s + vec1[(i1+c1) MOD 10] + vec2[(i2+c2) MOD 5] + '.dat'
       IF c2 EQ 4 THEN BEGIN
         c2 = 0
         c1 += 1
       ENDIF
     ENDIF ELSE spec_ok = 1
   ENDWHILE
   
   ; 5. Get TEFF if not found
   IF NOT KEYWORD_SET(teff) THEN BEGIN
     IF i2 GE 3 THEN BEGIN
       sp_class = ["O5","O6","O7","O8","O9","B0","B1","B2","B3","B5","B6","B7","B8","B9","A0","A1","A2","A3","A4","A5","A7","F0","F2","F3","F5","F6","F7","F8","G0","G1","G2","G5","G8","K0","K1","K2","K3","K4","K5","K7","K8","M0","M1","M2","M3","M4","M5","M6","M7","M8","L0","L3","L8","T2","T6","T8"]
       t_class  = [54000,45000,43300,40600,37800,29200,23000,21000,17600,15200,14300,13500,12300,11400,9600,9330,9040,8750,8480,8310,7920,7350,7050,6850,6700,6550,6400,6300,6050,5930,5800,5660,5440,5240,5110,4960,4800,4600,4400,4000,3850,3750,3700,3600,3500,3400,3200,3100,2900,2700,2600,2200,1500,1400,1000,800]
       ;m_class  = [32.7,29.5,28.3,25.2,22.6,17.8,11.7,10,7.32,5.46,4.79,4.1,3.41,2.91,2.48,2.35,2.25,2.13,2.03,1.98,1.86,1.59,1.46,1.39,1.33,1.22,1.12,1.1,1.05,1.03,1,0.91,0.82,0.76,0.72,0.7,0.67,0.61,0.58,0.53,0.49,0.44,0.38,0.36,0.34,0.28,0.22,0.19,0.17,0.098,0.078,0.052,0.048,0.032,0.028]
     ENDIF
     IF i2 LE 3 THEN BEGIN
       sp_class = ["G5","G8","K0","K1","K2","K3","K4","K5","M0","M1","M2","M3","M4","M5","M6"]
       t_class  = [5010,4870,4720,4580,4460,4210,4010,3780,3660,3600,3500,3300,3100,2950,2800]
     ENDIF
     idx  = WHERE(STRMID(sp_type, 0, 2) EQ STRLOWCASE(sp_class))
     teff = MEAN(t_class[idx])
     teff_ref = "Allen"
   ENDIF ELSE BEGIN
     IF N_ELEMENTS(teff) GT 1 THEN teff = teff[WHERE(teff NE -1)] & IF N_ELEMENTS(teff) GT 1 THEN teff = teff[0]
     IF N_ELEMENTS(mass) GT 1 THEN mass = mass[WHERE(mass NE -1)] & IF N_ELEMENTS(mass) GT 1 THEN mass = mass[0]
     IF N_ELEMENTS(logg) GT 1 THEN logg = logg[WHERE(logg NE -1)] & IF N_ELEMENTS(logg) GT 1 THEN logg = logg[0]
   ENDELSE
   
   ; 6. Get uld from Claret
   IF logg NE 0 AND teff NE 0 THEN BEGIN   
     logg_tab = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
     teff_tab = [3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000, 8250, 8500, 8750, 9000, 9250, 9500, 9750, 10000, 10500, 11000, 11500, 12000, 12500, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000, 32000, 33000, 34000, 35000, 37500, 40000, 42500, 45000, 47500, 50000]
  
     grid_h = [[0.482, 0.465, 0.442, 0.417, 0.392, 0.368, 0.343, 0.321, 0.334, 0.326, 0.322, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.478, 0.463, 0.441, 0.417, 0.392, 0.368, 0.343, 0.320, 0.297, 0.317, 0.307, 0.302, 0.300, 0.301, 0.305, 0.311, 0.321, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.475, 0.462, 0.441, 0.417, 0.392, 0.368, 0.345, 0.322, 0.298, 0.279, 0.303, 0.293, 0.285, 0.280, 0.278, 0.278, 0.280, 0.283, 0.289, 0.301, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.472, 0.460, 0.441, 0.418, 0.393, 0.370, 0.347, 0.325, 0.301, 0.280, 0.262, 0.290, 0.280, 0.271, 0.264, 0.260, 0.258, 0.256, 0.255, 0.257, 0.261, 0.266, 0.272, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.470, 0.458, 0.440, 0.418, 0.395, 0.372, 0.349, 0.328, 0.305, 0.284, 0.266, 0.251, 0.280, 0.269, 0.259, 0.251, 0.246, 0.241, 0.237, 0.234, 0.233, 0.234, 0.236, 0.238, 0.239, 0.240, 0.240, 0.238, 0.300, 0.300, 0.300, 0.300, 0.300, 0.242, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.465, 0.456, 0.440, 0.419, 0.397, 0.374, 0.352, 0.331, 0.309, 0.289, 0.271, 0.254, 0.242, 0.270, 0.259, 0.249, 0.240, 0.233, 0.227, 0.222, 0.218, 0.216, 0.215, 0.215, 0.215, 0.215, 0.214, 0.211, 0.206, 0.201, 0.198, 0.195, 0.193, 0.192, 0.192, 0.195, 0.201, 0.210, 0.222, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.454, 0.452, 0.439, 0.420, 0.398, 0.376, 0.355, 0.334, 0.314, 0.294, 0.276, 0.260, 0.245, 0.236, 0.261, 0.250, 0.240, 0.230, 0.222, 0.215, 0.209, 0.205, 0.202, 0.201, 0.200, 0.199, 0.198, 0.195, 0.190, 0.184, 0.179, 0.175, 0.171, 0.166, 0.163, 0.162, 0.162, 0.164, 0.167, 0.172, 0.177, 0.182, 0.187, 0.193, 0.201, 0.211, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.396, 0.440, 0.436, 0.420, 0.400, 0.379, 0.358, 0.338, 0.319, 0.300, 0.282, 0.266, 0.251, 0.238, 0.230, 0.253, 0.242, 0.231, 0.221, 0.212, 0.205, 0.199, 0.195, 0.192, 0.190, 0.189, 0.187, 0.184, 0.180, 0.175, 0.169, 0.164, 0.159, 0.153, 0.148, 0.145, 0.143, 0.142, 0.142, 0.143, 0.145, 0.149, 0.151, 0.153, 0.155, 0.157, 0.160, 0.163, 0.166, 0.168, 0.170, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.300, 0.391, 0.425, 0.418, 0.400, 0.381, 0.361, 0.342, 0.323, 0.305, 0.288, 0.272, 0.257, 0.245, 0.234, 0.225, 0.246, 0.234, 0.224, 0.214, 0.204, 0.196, 0.191, 0.187, 0.184, 0.182, 0.180, 0.177, 0.173, 0.169, 0.163, 0.158, 0.153, 0.145, 0.139, 0.135, 0.132, 0.130, 0.129, 0.128, 0.129, 0.131, 0.133, 0.135, 0.136, 0.137, 0.137, 0.138, 0.139, 0.139, 0.139, 0.138, 0.137, 0.137, 0.135, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.248, 0.308, 0.384, 0.408, 0.400, 0.382, 0.364, 0.345, 0.327, 0.310, 0.293, 0.277, 0.263, 0.250, 0.239, 0.229, 0.222, 0.240, 0.229, 0.219, 0.208, 0.198, 0.190, 0.185, 0.181, 0.178, 0.176, 0.172, 0.168, 0.164, 0.160, 0.155, 0.150, 0.142, 0.135, 0.130, 0.127, 0.124, 0.122, 0.121, 0.120, 0.121, 0.123, 0.124, 0.126, 0.127, 0.127, 0.127, 0.127, 0.127, 0.126, 0.124, 0.121, 0.119, 0.117, 0.111, 0.107, 0.300, 0.300, 0.300, 0.300],$
               [0.229, 0.253, 0.315, 0.376, 0.392, 0.382, 0.365, 0.347, 0.329, 0.313, 0.296, 0.282, 0.268, 0.256, 0.246, 0.236, 0.227, 0.221, 0.217, 0.225, 0.215, 0.204, 0.192, 0.186, 0.180, 0.177, 0.174, 0.169, 0.165, 0.162, 0.158, 0.153, 0.149, 0.141, 0.134, 0.129, 0.124, 0.121, 0.119, 0.117, 0.116, 0.116, 0.117, 0.119, 0.121, 0.122, 0.123, 0.122, 0.122, 0.121, 0.119, 0.117, 0.114, 0.111, 0.108, 0.101, 0.095, 0.093, 0.092, 0.092, 0.092]]
  
     grid_k = [[0.399, 0.387, 0.372, 0.354, 0.335, 0.316, 0.296, 0.278, 0.285, 0.275, 0.272, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.395, 0.386, 0.370, 0.353, 0.333, 0.316, 0.296, 0.277, 0.255, 0.268, 0.259, 0.254, 0.251, 0.253, 0.256, 0.260, 0.267, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.393, 0.384, 0.370, 0.352, 0.333, 0.315, 0.298, 0.279, 0.257, 0.239, 0.256, 0.247, 0.240, 0.235, 0.234, 0.233, 0.234, 0.235, 0.239, 0.248, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.391, 0.383, 0.370, 0.353, 0.334, 0.317, 0.299, 0.282, 0.261, 0.241, 0.226, 0.247, 0.237, 0.229, 0.223, 0.219, 0.217, 0.215, 0.213, 0.213, 0.215, 0.218, 0.224, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.389, 0.382, 0.370, 0.354, 0.336, 0.318, 0.301, 0.284, 0.265, 0.246, 0.230, 0.218, 0.239, 0.229, 0.220, 0.214, 0.209, 0.204, 0.200, 0.197, 0.194, 0.194, 0.195, 0.196, 0.197, 0.198, 0.198, 0.196, 0.300, 0.300, 0.300, 0.300, 0.300, 0.200, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.386, 0.381, 0.370, 0.356, 0.339, 0.321, 0.304, 0.287, 0.269, 0.252, 0.236, 0.222, 0.213, 0.231, 0.221, 0.213, 0.206, 0.200, 0.194, 0.189, 0.184, 0.181, 0.180, 0.179, 0.178, 0.178, 0.177, 0.175, 0.171, 0.167, 0.163, 0.161, 0.159, 0.158, 0.159, 0.161, 0.166, 0.174, 0.185, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.376, 0.378, 0.371, 0.357, 0.340, 0.323, 0.306, 0.290, 0.273, 0.257, 0.241, 0.228, 0.217, 0.210, 0.224, 0.215, 0.206, 0.199, 0.192, 0.185, 0.179, 0.175, 0.172, 0.170, 0.168, 0.167, 0.165, 0.162, 0.158, 0.154, 0.149, 0.145, 0.142, 0.138, 0.135, 0.134, 0.135, 0.136, 0.139, 0.144, 0.149, 0.155, 0.159, 0.164, 0.170, 0.178, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.324, 0.368, 0.369, 0.358, 0.342, 0.326, 0.309, 0.294, 0.278, 0.263, 0.247, 0.234, 0.222, 0.212, 0.206, 0.218, 0.209, 0.200, 0.192, 0.185, 0.178, 0.172, 0.167, 0.164, 0.162, 0.160, 0.158, 0.155, 0.151, 0.147, 0.142, 0.137, 0.133, 0.128, 0.124, 0.121, 0.119, 0.118, 0.119, 0.120, 0.123, 0.126, 0.129, 0.132, 0.133, 0.135, 0.137, 0.139, 0.141, 0.142, 0.144, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.247, 0.322, 0.358, 0.357, 0.343, 0.328, 0.312, 0.297, 0.282, 0.267, 0.252, 0.239, 0.226, 0.217, 0.209, 0.202, 0.212, 0.203, 0.195, 0.186, 0.178, 0.171, 0.166, 0.162, 0.159, 0.156, 0.154, 0.150, 0.147, 0.142, 0.138, 0.133, 0.129, 0.123, 0.118, 0.114, 0.111, 0.109, 0.108, 0.108, 0.109, 0.111, 0.114, 0.116, 0.118, 0.119, 0.119, 0.120, 0.120, 0.119, 0.118, 0.117, 0.116, 0.116, 0.116, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],$
               [0.207, 0.254, 0.318, 0.347, 0.343, 0.329, 0.315, 0.300, 0.284, 0.271, 0.256, 0.243, 0.231, 0.221, 0.212, 0.205, 0.199, 0.207, 0.199, 0.190, 0.181, 0.172, 0.166, 0.161, 0.157, 0.154, 0.152, 0.148, 0.144, 0.140, 0.136, 0.131, 0.127, 0.120, 0.115, 0.111, 0.107, 0.105, 0.103, 0.102, 0.102, 0.103, 0.105, 0.107, 0.109, 0.111, 0.111, 0.111, 0.110, 0.109, 0.108, 0.106, 0.103, 0.101, 0.100, 0.096, 0.094, 0.300, 0.300, 0.300, 0.300],$
               [0.192, 0.212, 0.261, 0.314, 0.335, 0.329, 0.315, 0.301, 0.286, 0.273, 0.259, 0.246, 0.234, 0.224, 0.216, 0.209, 0.202, 0.198, 0.194, 0.194, 0.186, 0.177, 0.167, 0.161, 0.156, 0.153, 0.150, 0.146, 0.142, 0.138, 0.134, 0.130, 0.126, 0.120, 0.114, 0.109, 0.106, 0.103, 0.101, 0.100, 0.099, 0.099, 0.101, 0.103, 0.105, 0.106, 0.107, 0.106, 0.106, 0.104, 0.103, 0.100, 0.097, 0.094, 0.092, 0.087, 0.083, 0.082, 0.081, 0.081, 0.081]]
     tmp  = MIN(ABS(logg_tab-logg), i1)
     tmp  = MIN(ABS(teff_tab-teff), i2)
     uldh = grid_h[i2,i1] & uldh_ref  = 'Clar95'
     uldk = grid_k[i2,i1] & uldk_ref  = 'Clar95'
   ENDIF ELSE BEGIN
     uldh  = 0. & uldh_ref  = 'none00'
     uldk  = 0. & uldk_ref  = 'none00'
   ENDELSE
   
   ; 7. Get PA and vsini
   vsini = 0. & vsini_ref = 'none00'
   pa    = 0. & pa_ref    = 'none00' 
 
   RETURN, {NAME: input_name, RA: ra, DEC: de, SP_TYPE: sp_type, SPECTRUM: spectrum, DIST: FLOAT(1D3/p), VMAG: FLOAT(vmag), E_VMAG: FLOAT(e_vmag), VMAG_REF: vmag_ref, HMAG: FLOAT(Hmag), E_HMAG: FLOAT(e_hmag), HMAG_REF: hmag_ref, $
            KMAG: FLOAT(kmag), E_KMAG: FLOAT(e_kmag), KMAG_REF: kmag_ref, LMAG: FLOAT(lmag), E_LMAG: FLOAT(e_lmag), LMAG_REF: lmag_ref, TEFF: FIX(teff), TEFF_REF: teff_ref, MASS: FLOAT(mass), MASS_REF: mass_ref, LDM: FLOAT(diam[0]), $ 
            LDE: FLOAT(diam[1]), LDM_REF: diam_ref, LOGG: FLOAT(logg), LOGG_REF: logg_ref, VSINI: FLOAT(vsini), VSINI_REF: vsini_ref, ULDH: FLOAT(uldh), ULDH_REF: uldh_ref, ULDK: FLOAT(uldk), ULDK_REF: uldk_ref, PA: FLOAT(pa), PA_REF: pa_ref}
END