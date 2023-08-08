;+
; NAME: GET_TGT
; 
; DESCRIPTION:
;   This procedure reads the star database and returns the information in an array of structure (one for each target).
;   If the target is not found, it will query SIMBAD and VIZIER to find the information.
;
; INPUTS:
;   names:  Name of the target to retrieve
;   
; KEYWORD
;   DATABASE: Set this keyword to the file name with your star database (will be created if non existant)
;   HBAND   : Set this keyword to return the H-band limb-darkining coefficient (K-band by default)
;   
; OUTPUT
;   tgt:    An array of structures with each the following 12 elements:
;             - name     : name of the target
;             - spectrum : spectrum file of the target
;             - ra       : target right ascension
;             - dec      : target declination
;             - dist     : target distance in pc
;             - mass     : target mass in solar mass
;             - temp     : target temperature in K
;             - vsini    : target projected rotationnal velocity in km/s
;             - ldm      : target limb-darkened angular diameter in mas
;             - lde      : error on the limb-darkened angular diameter
;             - uld      : limb darkeining linear coefficient (WARNING: only valid at one wavelength!)
;             - pa       : position angle of the stellar photosphere
;
; MODIFICATION HISTORY:
;   Version 1.0, 04-MAY-2012, Denis Defrère, ddefrere@ulg.ac.be
;   Version 1.1, 04-DEC-2015, DD: added calfor to output structure
;   Version 1.2, 03-APR-2017, DD: now read database and query online catalogs if not found
;   Version 1.3, 06-APR-2017, DD: added keyword DATABASE

PRO GET_TGT, names, tgt, DATABASE=database, HBAND=hband

  ; Read star database
  IF NOT KEYWORD_SET(DATABASE) THEN database = 'star_database.dat'
  IF FILE_TEST(database) THEN READ_TABLE,database,db_names,ra,dec,spectrum,dist,temp,mass,ldm,lde,vsini,uldh,uldk,pa,FIRST=2,SKIP=[3,6+INDGEN(12),19,21,24,26,27,28,30,32,34,35],STRING_ARRAY=[0,3], SEP=';' $
                         ELSE db_names = ['none']
  db_names = strlowcase(db_names)
                           
  ; Number of target to find
  n_tgt = N_ELEMENTS(names)
  
  ; Loop over the target names
  FOR i_tgt = 0, n_tgt-1 DO BEGIN
    ; Look in the database
    idx_ok = WHERE(STRCOMPRESS(STRING(STRLOWCASE(names[i_tgt])), /REMOVE_ALL) EQ db_names, n_ok)
    ; If not found, QUERY SIMBAD/VIZIER and save results
    IF n_ok LE 0 THEN BEGIN
      ; Print info
      PRINT, "QUERYING SIMBAD AND VIZIER FOR TARGET INFO (" + names[i_tgt] + ")"
      d = QUERY_TGT(names[i_tgt])
      ; H or K band uld?
      IF KEYWORD_SET(HBAND) THEN uld = d.uldh ELSE uld = d.uldk
      ; Create output structre
      tgt_tmp = {name: STRCOMPRESS(STRING(STRLOWCASE(d.name)), /REMOVE_ALL), spectrum:d.spectrum, ra:FLOAT(d.ra), dec:FLOAT(d.dec), dist:FLOAT(d.dist), mass:FLOAT(d.mass), temp:FIX(d.teff), vsini:FLOAT(d.vsini), ldm:FLOAT(d.ldm), lde:FLOAT(d.lde), uld:FLOAT(uld), pa:FLOAT(d.pa)}
      ; Append results to database (create it if inexistant)
      IF NOT FILE_TEST(database) THEN BEGIN
        OPENW,  dbu, database, /GET_LUN, WIDTH=1000
        PRINTF, dbu, 'STAR DATABASE'
        PRINTF, dbu, '0-NAME;1-RA;2-DEC;3-SPT;4-SPEC;5-DST;6-VMAG;7-EVMAG;8-VREF;9-HMAG;10-EHMAG;11-HREF;12-KMAG;13-EKMAG;14-KREF;15-LMAG;16-ELMAG;17-LREF;18-TEFF;19-TREF;20-MASS;21-MREF;22-LDM;23-LDE;24-LDM_REF;25-VSINI;26-VSINI_REF;27-LOGG;28-LOGG-REF;29-ULDH;30-ULDH_REF;31-ULDK;32-ULDK_REF;33-PA;34-PA_REF;35-DATE'        
      ENDIF ELSE OPENW, dbu, database, /GET_LUN, /APPEND, WIDTH=1000
      PRINTF, dbu, STRCOMPRESS(STRING(d.name), /REMOVE_ALL), ';', STRING(d.ra, FORMAT='(F6.2)'), ';', STRING(d.dec, FORMAT='(F5.1)'), ';', STRCOMPRESS(STRING(d.sp_type), /REMOVE_ALL), ';', STRCOMPRESS(STRING(d.spectrum), /REMOVE_ALL), ';', $
                   STRING(d.dist, FORMAT='(F6.1)'), ';', STRING(d.vmag, FORMAT='(F5.2)'), ';', STRING(d.e_vmag, FORMAT='(F5.3)'), ';', STRCOMPRESS(STRING(d.vmag_ref), /REMOVE_ALL), ';', STRING(d.hmag, FORMAT='(F5.2)'), ';', STRING(d.e_hmag, FORMAT='(F5.3)'), ';', $
                   STRCOMPRESS(STRING(d.hmag_ref), /REMOVE_ALL), ';', STRING(d.kmag, FORMAT='(F5.2)'), ';', STRING(d.e_kmag, FORMAT='(F5.3)'), ';', STRCOMPRESS(STRING(d.kmag_ref), /REMOVE_ALL), ';', STRING(d.lmag, FORMAT='(F5.2)'), ';', $
                   STRING(d.e_lmag, FORMAT='(F5.3)'), ';', STRCOMPRESS(STRING(d.lmag_ref), /REMOVE_ALL),';', STRING(d.teff, FORMAT='(I05)'), ';', STRCOMPRESS(STRING(d.teff_ref), /REMOVE_ALL), ';', STRING(d.mass, FORMAT='(F6.2)'), ';', STRCOMPRESS(STRING(d.mass_ref), /REMOVE_ALL), ';', $
                   STRING(d.ldm, FORMAT='(F5.2)'), ';', STRING(d.lde, FORMAT='(F4.2)'), ';', STRCOMPRESS(STRING(d.ldm_ref), /REMOVE_ALL), ';', STRING(d.vsini, FORMAT='(F5.2)'), ';', STRCOMPRESS(STRING(d.vsini_ref), /REMOVE_ALL), ';', STRING(d.logg, FORMAT='(F5.3)'), ';', STRCOMPRESS(STRING(d.logg_ref), /REMOVE_ALL), ';', $
                   STRING(d.uldh, FORMAT='(F4.2)'), ';', STRCOMPRESS(STRING(d.uldh_ref), /REMOVE_ALL), ';', STRING(d.uldk, FORMAT='(F4.2)'), ';', STRCOMPRESS(STRING(d.uldk_ref), /REMOVE_ALL), ';', STRING(d.pa, FORMAT='(F5.2)'), ';', STRCOMPRESS(STRING(d.pa_ref), /REMOVE_ALL),  ';', $
                   SYSTIME()
      CLOSE, dbu
      FREE_LUN, dbu            
    ENDIF ELSE BEGIN     
      IF KEYWORD_SET(HBAND) THEN uld = uldh ELSE uld = uldk ; H or K band uld?
      tgt_tmp = {name: (db_names[idx_ok])[0], spectrum: (spectrum[idx_ok])[0], ra: (FLOAT(ra[idx_ok]))[0], dec:(FLOAT(dec[idx_ok]))[0], dist:(FLOAT(dist[idx_ok]))[0], mass:(FLOAT(mass[idx_ok]))[0], temp:(FIX(temp[idx_ok]))[0], vsini:(FLOAT(vsini[idx_ok]))[0], ldm:(FLOAT(ldm[idx_ok]))[0], lde:(FLOAT(lde[idx_ok]))[0], $
                 uld:(FLOAT(uld[idx_ok]))[0], pa:(FLOAT(pa[idx_ok]))[0]}
    ENDELSE
    ; Now append to output data
    IF i_tgt EQ 0 THEN tgt = tgt_tmp ELSE tgt = [tgt, tgt_tmp]
  ENDFOR
END

PRO BUILD_HOSTS
  stars = ['6_Cet','HD4628','ups_and','107_Psc','tau_cet','GJ_75','GJ_105','13_Per','tau01_Eri','iot_per','eps_eri','LHS1569','tau06_eri','omi02_eri','1_Ori','HD32147','lam_aur','V1119_Tau','gam_lep','ksi_gem',$
          'sig02_uma', 'GJ_364', 'NSV_4765', '40_Leo', '36_uma', '47_Uma', '61_Uma', 'bet_Vir', '61_Vir', 'tau_boo', 'tet_boo', 'KX_Lib', 'lam_ser', 'LHS_3127', 'gam_ser', 'V2133_Oph', 'V2215_Oph', 'w_Her',$
          '58_Oph','110_Her','sig_Dra','GJ_785','psi_cap','61_CygA', '61_CygB','ksi_Peg','HD219134', 'iot_Psc','bet_Eri', 'zet_lep','eta_lep','h_UMa','beta_UMa','del_Leo','bet_leo','gam_UMa','alf_Crv','del_UMa',$
          'del_Crv','eta_Crv', 'sig_Boo','107_Vir','zet_Ser','Vega','alf_Lyr','alp_Lyr','Altair','alf_Aql','alp_Aql','Alderamin','Alf_Cep','tet_Peg','Fomalhaut', 'alf_Psa']
  GET_TGT, stars, tgt, DATABASE='hosts_db.dat'
END

PRO BUILD_PIONIER
  stars = ['hd203','hd2834','hd3126','hd4113','hd9672','hd10269','hd10939','hd15427','hd17848','hd23484','hd24649','hd28287','hd29137','hd36187','hd37306','hd37484','hd38949','hd41278','hd44524',$
           'hd60491','hd61005','hd71722','hd76143','hd80883','hd89886','hd90781','hd90874','hd92945','hd105850','hd105912','hd109573','hd109704','hd112603','hd117716','hd118972','hd136544',$
           'hd141943','hd161612','hd174474','hd179520','hd181327','hd185615','hd191089','hd192758','hd196141','hd205674','hd220476','hd224228',$
           'hd142','hd1581','hd2262','hd3302','hd3823','hd4150','hd7570','hd7788','hd10647','hd11171','hd14412','hd15008','hd15798','hd16555','hd17051','hd17925','hd19107','hd20766','hd20794',$
           'hd20807','hd22001','hd23249','hd25457','hd28355','hd29388','hd30495','hd31295','hd31925','hd33111','hd33262','hd34721','hd38858','hd39060','hd40307','hd43162','hd45184','hd53705',$
           'hd56537','hd69830','hd71155','hd72673','hd76151','hd76932','hd82434','hd88955','hd90132','hd91324','hd99211','hd102365','hd104731','hd108767','hd109787','hd115617','hd120136','hd128898',$
           'hd129502','hd130109','hd134083','hd135379','hd136202','hd139664','hd141891','hd149661','hd152391','hd160032','hd160915','hd164259','hd165777','hd172555','hd178253','hd182572','hd188228',$
           'hd192425','hd195627','hd197157','hd197692','hd202730','hd203608','hd206860','hd207129','hd210049','hd210277','hd210302','hd210418','hd213845','hd214953','hd215648','hd215789','hd216435',$
           'hd219482','hd219571','hd224392']
  GET_TGT, stars, tgt, DATABASE='pionierH_database.dat', /HBAND
END

; Old verions of the routine
PRO GET_TGT_OLD,  input_name, CALIB=calib
; Definition of targets
; *********************
IF STRTRIM(input_name) EQ 'Vega' THEN input_name = 'alf_Lyr'

alf_crb = {name:'alf_Crb', spectrum:'uka0v.dat', ra:TEN(15,34,41), dec:TEN(26,42,52), dist:22.9, mass:2.58, temp:9790., vsini:138., ldm:1.202, lde:0.024, uld:0.173, pa:0.0, calfor:'none'}
gam_oph = {name:'gam_Oph', spectrum:'uka0v.dat', ra:TEN(17,47,54), dec:TEN(02,42,26), dist:29.1, mass:2.18, temp:9100., vsini:210., ldm:0.616, lde:0.012, uld:0.183, pa:0.0, calfor:'none'}

; To be checked (temperature, mass, ld)
zet_aql = {name:'zet_Aql', spectrum:'uka0v.dat', ra:TEN(19,05,25), dec:TEN(13,51,49), dist:25.5, mass:2.37, temp:9400., vsini:317., ldm:0.880, lde:0.018, uld:0.162, pa:0.0, calfor:'none'}


;alf_lyr = {name:'alf_Lyr', spectrum:'uka0v.dat', ra:TEN(18,36,56), dec:TEN(38,47,01), dist:7.76, mass:2.3, temp:10000., vsini:21.9, ldm:3.202, lde:0.005, uld:0.0, pa:40.0, calfor:'none'} ; from Absil 2006
;alf_lyr = {name:'alf_Lyr', spectrum:'uka0v.dat', ra:TEN(18,36,56), dec:TEN(38,47,01), dist:7.76, mass:2.3, temp:10000., vsini:21.9, ldm:3.198, lde:0.003, uld:0.0, pa:40.0, calfor:'none'} ; UD model from Aufdenberg 2006
;alf_lyr = {name:'alf_Lyr', spectrum:'uka0v.dat', ra:TEN(18,36,56), dec:TEN(38,47,01), dist:7.76, mass:2.3, temp:10000., vsini:21.9, ldm:3.329, lde:0.006, uld:0.47, pa:40.0, calfor:'none'} ; linear LD model derived from Aufdenberg 2006
;alf_lyr = {name:'alf_Lyr', spectrum:'uka0v.dat', ra:TEN(18,36,56), dec:TEN(38,47,01), dist:7.76, mass:2.3, temp:10000., vsini:21.9, ldm:3.313, lde:0.006, uld:0.47, pa:40.0, calfor:'none'} ; linear LD model derived from Aufdenberg 2006

;alf_psa = {name:'alf_PsA', spectrum:'uka4v.dat', ra:TEN(22,57,39), dec:-TEN(29,37,20), dist:7.7, mass:2.0, temp:8760., vsini:93., ldm:2.178, lde:0.111, uld:0.189, pa:156.0, calfor:'none'}
;alf_psa = {name:'alf_PsA', spectrum:'uka4v.dat', ra:TEN(22,57,39), dec:-TEN(29,37,20), dist:7.7, mass:2.0, temp:8760., vsini:93., ldm:2.197, lde:0.023, uld:0.0, pa:156.0, calfor:'none'} ; mean UD from Di Folco 04
;alf_psa = {name:'alf_PsA', spectrum:'uka4v.dat', ra:TEN(22,57,39), dec:-TEN(29,37,20), dist:7.7, mass:2.0, temp:8760., vsini:93., ldm:2.228, lde:0.023, uld:0.189, pa:156.0, calfor:'none'} ; mean LD from Di Folco 04
; According to Di Folco, UD diameter is 2.191 mas (LD=2.222mas) along the minor axis (30° azimuth of siderostat baseline --> parallel to photospheric minor axis)
; Assuming an oblateness of 1.021, this gives a geometric mean UD diameter of 2.214 mas (LD=2.245mas)
;alf_psa = {name:'alf_PsA', spectrum:'uka4v.dat', ra:TEN(22,57,39), dec:-TEN(29,37,20), dist:7.7, mass:2.0, temp:8760., vsini:93., ldm:2.214, lde:0.023, uld:0.0, pa:156.0, calfor:'none'}   ; geom mean UD from Di Folco 04 based on long baselines
;alf_psa = {name:'alf_PsA', spectrum:'uka4v.dat', ra:TEN(22,57,39), dec:-TEN(29,37,20), dist:7.7, mass:2.0, temp:8760., vsini:93., ldm:2.245, lde:0.023, uld:0.189, pa:156.0, calfor:'none'} ; geom mean LD from Di Folco 04 based on long baselines
alf_psa = {name:'alf_PsA', spectrum:'uka4v.dat', ra:TEN(22,57,39), dec:-TEN(29,37,20), dist:7.7, mass:2.0, temp:8760., vsini:93., ldm:2.223, lde:0.022, uld:0.189, pa:156.0, calfor:'none'} ; based on the 2008 VINCI study

; SURVEY 2008
gam_ser = {name:'gam_Ser', spectrum:'ukf6v.dat', ra:TEN(15,56,27), dec:TEN(15,39,42), dist:11.14, mass:1.16, temp:6262, vsini:10.9, ldm:1.157, lde:0.059, uld:0.26, pa:0.0, calfor:'none'}
; Mass, Teff and vsini from Valenti05 -- uld TBC
O70_oph = {name:'70_Oph', spectrum:'ukk0v.dat', ra:TEN(18,05,27), dec:TEN(02,30,00), dist:5.10, mass:0.97, temp:5370, vsini:1, ldm:1.986, lde:0.101, uld:0.33, pa:0.0, calfor:'none'}
; Mass and Teff from AllendePrieto99 -- uld TBC
lam_ser = {name:'lam_Ser', spectrum:'ukg0v.dat', ra:TEN(15,46,27), dec:TEN(07,21,11), dist:11.78, mass:1.01, temp:5801, vsini:1, ldm:1.067, lde:0.055, uld:0.28, pa:0.0, calfor:'none'}
; Mass and Teff from Lambert04 -- uld TBC
alf_cep = {name:'alf_Cep', spectrum:'uka7v.dat', ra:TEN(21,18,35), dec:TEN(62,35,08), dist:14.99, mass:2.0, temp:7773., vsini:196., ldm:1.557, lde:0.080, uld:0.22, pa:0.0, calfor:'none'}
; alf Cep: diameter ok, vsini ok, temp from Gray, mass from vanBelle06
eps_cep = {name:'eps_Cep', spectrum:'ukf0v.dat', ra:TEN(22,15,02), dec:TEN(57,02,37), dist:25.78, mass:1.76, temp:7586., vsini:91., ldm:0.722, lde:0.037, uld:0.22, pa:0.0, calfor:'none'}
; eps Cep: main parameters from AllendePrieto99
O54_psc = {name:'54_Psc', spectrum:'ukk0v.dat', ra:TEN(00,39,22), dec:TEN(21,15,01), dist:11.12, mass:1, temp:5000, vsini:1, ldm:0.721, lde:0.015, uld:0.33, pa:0.0, calfor:'none'}
; WARNING - NEEDS TO BE UPDATED!!!!
eta_cas = {name:'eta_Cas', spectrum:'ukg3v.dat', ra:TEN(00,49,06), dec:TEN(57,48,55), dist:5.96, mass:0.90, temp:5920., vsini:0., ldm:1.633, lde:0.083, uld:0.28, pa:0.0, calfor:'none'}
; eta Cas: mass from AllendePrieto99, Teff from Gray01, vsini not found, log g = 4.2
ups_and = {name:'ups_And', spectrum:'ukf9v.dat', ra:TEN(01,36,48), dec:TEN(41,24,20), dist:13.49, mass:1.37, temp:6210., vsini:0., ldm:1.106, lde:0.057, uld:0.27, pa:0.0, calfor:'none'}
; ups And: mass from AllendePrieto99, Teff from Gray01, vsini not found, log g = 4.1
tet_per = {name:'tet_Per', spectrum:'ukf7v.dat', ra:TEN(02,44,12), dec:TEN(49,13,42), dist:11.25, mass:1.15, temp:6330, vsini:0., ldm:1.017, lde:0.052, uld:0.27, pa:0.0, calfor:'none'}
; tet Per: mass from AllendePrieto99, Teff from Gray01, vsini not found, log g = 4.2
O10_tau = {name:'10_Tau', spectrum:'ukf8v.dat', ra:TEN(03,36,52), dec:TEN(00,24,06), dist:13.74, mass:1.25, temp:5960., vsini:0., ldm:1.100, lde:0.056, uld:0.28, pa:0.0, calfor:'none'}
; 10 Tau: mass from AllendePrieto99, Teff from Gray01, vsini not found, log g = 3.9
OO1_ori = {name:'1_Ori', spectrum:'ukf6v.dat', ra:TEN(04,49,50), dec:TEN(06,57,41), dist:8.04, mass:1.24, temp:6530., vsini:0., ldm:1.522, lde:0.078, uld:0.26, pa:0.0, calfor:'none'}
; 1 Ori: mass from AllendePrieto99, Teff from Gray01, vsini not found, log g = 4.3
zet_lep = {name:'zet_Lep', spectrum:'uka2v.dat', ra:TEN(05,46,57), dec:-TEN(14,49,19), dist:21.56, mass:1.85, temp:8511., vsini:229., ldm:0.759, lde:0.017, uld:0.20, pa:0.0, calfor:'none'}
;zet_lep = {name:'zet_Lep', spectrum:'uka2v.dat', ra:TEN(05,46,57), dec:-TEN(14,49,19), dist:21.56, mass:1.85, temp:8511., vsini:229., ldm:0.670, lde:0.140, uld:0.20, pa:0.0, calfor:'none'}
; zet Lep: mass and Teff from AllendePrieto99, log g = 4.3
eta_lep = {name:'eta_Lep', spectrum:'ukf2v.dat', ra:TEN(05,56,24), dec:-TEN(14,10,04), dist:15.07, mass:1.51, temp:7079., vsini:26., ldm:0.992, lde:0.022, uld:0.24, pa:0.0, calfor:'none'}
; eta Lep: mass and Teff from AllendePrieto99, log g = 4.2
lam_gem = {name:'lam_Gem', spectrum:'uka3v.dat', ra:TEN(07,18,06), dec:TEN(16,32,25), dist:28.96, mass:2.17, temp:8511., vsini:154., ldm:0.719, lde:0.037, uld:0.20, pa:0.0, calfor:'none'}
; lam Gem: mass and Teff from AllendePrieto99, log g = 4.1
ksi_boo = {name:'ksi_Boo', spectrum:'ukg8v.dat', ra:TEN(14,51,23), dec:TEN(19,06,02), dist:6.71, mass:0.98, temp:5570., vsini:4.6, ldm:1.431, lde:0.073, uld:0.33, pa:0.0, calfor:'none'}
; Mass, Teff and vsini from Valenti05 -- uld TBC
sig_dra = {name:'sig_Dra', spectrum:'ukg9v.dat', ra:TEN(19,32,22), dec:TEN(69,39,40), dist:5.77, mass:0.758, temp:5246, vsini:1.4, ldm:1.183, lde:0.060, uld:0.33, pa:0.0, calfor:'none'}
; Mass, Teff and vsini from Valenti05 -- uld TBC
O70_vir = {name:'70_Vir', spectrum:'ukg5v.dat', ra:TEN(13,28,25), dec:TEN(13,46,44), dist:18.14, mass:1.48, temp:5545, vsini:2.7, ldm:0.852, lde:0.044, uld:0.33, pa:0.0, calfor:'none'}
; Mass, Teff and vsini from Valenti05 -- uld TBC
iot_vir = {name:'iot_Vir', spectrum:'ukf7v.dat', ra:TEN(14,16,01), dec:-TEN(06,00,02), dist:21.43, mass:1.53, temp:6180, vsini:1, ldm:1.136, lde:0.058, uld:0.28, pa:0.0, calfor:'none'}
; Teff from Gray, Mass from AllendePrieto99 -- uld TBC
del_uma = {name:'del_UMa', spectrum:'uka3v.dat', ra:TEN(12,15,26), dec:TEN(57,01,57), dist:25.01, mass:2.1, temp:8709, vsini:233, ldm:0.820, lde:0.042, uld:0.22, pa:0.0, calfor:'none'}
; Mass and Teff from AllendePrieto99 -- uld TBC
HD_69830 = {name:'HD_69830', spectrum:'ukg8v.dat', ra:TEN(08,18,23), dec:-TEN(12,37,55), dist:12.6, mass:0.8, temp:5385, vsini:0, ldm:0.63, lde:0.1, uld:0.3, pa:0.0, calfor:'none'}
; EVERYTHING NEEDS TO BE CHECKED (only guessing here)

; A-M Lagrange data
;bet_pic = {name:'beta_pic', spectrum:'uka5v.dat', ra:TEN(05,47,17), dec:TEN(-51,03,59), dist:19.3, mass:1.75, temp:8200., vsini:130.0, ldm:0.712, lde:0.010, uld:0.24, pa:119.5, calfor:'none'}   ; voir infos dans mail Oli (04.11.11), mass from  Crifo et al. 1997, v.sini from Royer et al. 2007
bet_pic = {name:'beta_pic', spectrum:'uka7v.dat', ra:TEN(05,47,17), dec:TEN(-51,03,59), dist:19.3, mass:1.75, temp:8200., vsini:130.0, ldm:0.736, lde:0.019, uld:0.24, pa:119.5, calfor:'none'}   ; Diameter from Defrere et al. 2012
;bet_pic = {name:'beta_pic', spectrum:'uka7v.dat', ra:TEN(05,47,17), dec:TEN(-51,03,59), dist:19.3, mass:1.75, temp:8200., vsini:130.0, ldm:0.875, lde:0.019, uld:0.24, pa:119.5, calfor:'none'}    ; Best-fit without disk
bet_picK = {name:'bet_pic', spectrum:'uka7v.dat', ra:TEN(05,47,17), dec:TEN(-51,03,59), dist:19.3, mass:1.75, temp:8200., vsini:130.0, ldm:0.736, lde:0.019, uld:0.20, pa:119.5, calfor:'none'}   ; K band uld, diameter from Defrere et al. 2012
au_mic  = {name:'au_mic', spectrum:'ukm1v.dat', ra:TEN(20,45,09), dec:TEN(-31,20,27), dist:9.9, mass:0.4, temp:3750., vsini:0.0, ldm:0.785, lde:0.040, uld:0.50, pa:0.0, calfor:'none'}

; SURVEY PIONIER (WARNING, mass and temp taken from stars above with the same spectral type. uld not updated yet and LD diameter computed with 2mass K magnitudes)
hd142    = {name:'hd_142',    spectrum:'ukf6v.dat',  ra:TEN(00,06,19), dec:TEN(-49,04,31), dist:25.71, mass:1.16, temp:6213.,  vsini:11.0,  ldm:0.519, lde:0.008, uld:0.272, pa:0.0, calfor:'none'}
hd7570   = {name:'hd_7570',   spectrum:'ukf9v.dat',  ra:TEN(01,15,11), dec:TEN(-45,31,54), dist:15.11, mass:1.37, temp:6018.,  vsini:4.3,   ldm:0.762, lde:0.010, uld:0.293, pa:0.0, calfor:'none'}
hd7788   = {name:'hd_7788',   spectrum:'ukf6v.dat',  ra:TEN(01,15,46), dec:TEN(-68,52,33), dist:20.96, mass:1.24, temp:6366.,  vsini:61.0,  ldm:0.739, lde:0.011, uld:0.257, pa:0.0, calfor:'none'}
hd15008  = {name:'hd_15008',  spectrum:'uka3v.dat',  ra:TEN(02,21,45), dec:TEN(-68,39,34), dist:42.83, mass:2.17, temp:9001.,  vsini:190.0, ldm:0.541, lde:0.007, uld:0.191, pa:0.0, calfor:'none'}
hd90132  = {name:'hd_90132',  spectrum:'uka7v.dat',  ra:TEN(10,23,29), dec:TEN(-38,00,35), dist:41.46, mass:2.0,  temp:7408.,  vsini:270.0, ldm:0.426, lde:0.006, uld:0.246, pa:0.0, calfor:'none'}  ; no logg measurement found. guess: 4.0
hd91324  = {name:'hd_91324',  spectrum:'ukf8v.dat',  ra:TEN(10,31,22), dec:TEN(-53,42,56), dist:21.81, mass:1.25, temp:6096.,  vsini:8.3,   ldm:0.793, lde:0.080, uld:0.288, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd99211  = {name:'hd_99211',  spectrum:'uka7v.dat',  ra:TEN(11,24,53), dec:TEN(-17,41,02), dist:25.24, mass:2.00, temp:7620.,  vsini:7.3,   ldm:0.705, lde:0.188, uld:0.246, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd102365 = {name:'hd_102365', spectrum:'ukg2v.dat',  ra:TEN(11,46,31), dec:TEN(-40,30,01), dist:9.22,  mass:1.00, temp:5575.,  vsini:0.7,   ldm:0.943, lde:0.014, uld:0.327, pa:0.0, calfor:'none'}
hd104731 = {name:'hd_104731', spectrum:'ukf5v.dat',  ra:TEN(12,03,40), dec:TEN(-42,26,03), dist:25.32, mass:1.50, temp:6385.,  vsini:15.9,  ldm:0.604, lde:0.008, uld:0.257, pa:0.0, calfor:'none'}
hd108767 = {name:'hd_108767', spectrum:'uka0iv.dat', ra:TEN(12,29,52), dec:TEN(-16,30,56), dist:26.63, mass:2.9,  temp:10371., vsini:236.0, ldm:0.792, lde:0.011, uld:0.177, pa:0.0, calfor:'none'}
hd109787 = {name:'hd_109787', spectrum:'uka2v.dat',  ra:TEN(12,37,42), dec:TEN(-48,32,29), dist:40.2,  mass:0.70, temp:8762.,  vsini:296.8, ldm:0.577, lde:0.008, uld:0.196, pa:0.0, calfor:'none'}  ; no logg measurement found. guess: 4.0
hd115617 = {name:'hd115617',  spectrum:'ukg8v.dat',  ra:TEN(13,18,24), dec:TEN(-18,18,40), dist:8.58,  mass:0.83, temp:5557.,  vsini:3.9,   ldm:1.147, lde:0.129, uld:0.327, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd120136 = {name:'hd_120136', spectrum:'ukf5iv.dat', ra:TEN(12,47,16), dec:TEN(+17,27,25), dist:15.6,  mass:1.50, temp:6449.,  vsini:15.6,  ldm:0.856, lde:0.116, uld:0.263, pa:0.0, calfor:'none'}
hd128898 = {name:'hd_128898', spectrum:'uka7v.dat',  ra:TEN(14,42,30), dec:TEN(-64,58,30), dist:16.57, mass:2.0,  temp:7426.,  vsini:13.2,  ldm:1.005, lde:0.025, uld:0.246, pa:0.0, calfor:'none'}
hd129502 = {name:'hd_129502', spectrum:'ukf2v.dat',  ra:TEN(14,43,04), dec:TEN(-05,39,30), dist:18.3,  mass:1.50, temp:6640.,  vsini:47.0,  ldm:1.027, lde:0.014, uld:0.245, pa:0.0, calfor:'none'}
hd130109 = {name:'hd_130109', spectrum:'uka0v.dat',  ra:TEN(14,46,15), dec:TEN(+01,53,34), dist:41.24, mass:2.30, temp:9683.,  vsini:258.0, ldm:0.613, lde:0.021, uld:0.178, pa:0.0, calfor:'none'}
hd134083 = {name:'hd_134083', spectrum:'ukf5v.dat',  ra:TEN(15,07,19), dec:TEN(+24,52,09), dist:19.55, mass:1.50, temp:6567.,  vsini:44.4,  ldm:0.661, lde:0.007, uld:0.263, pa:0.0, calfor:'none'}
hd135379 = {name:'hd_135379', spectrum:'uka3v.dat',  ra:TEN(15,17,30), dec:TEN(-58,48,04), dist:30.55, mass:2.0,  temp:8695.,  vsini:68.5,  ldm:0.570, lde:0.008, uld:0.196, pa:0.0, calfor:'none'}
hd136202 = {name:'hd_136202', spectrum:'ukf8v.dat',  ra:TEN(15,19,19), dec:TEN(+01,45,55), dist:25.38, mass:1.15, temp:6129.,  vsini:4.8,   ldm:0.623, lde:0.085, uld:0.272, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd139664 = {name:'hd_139664', spectrum:'ukf5v.dat',  ra:TEN(15,41,11), dec:TEN(-44,39,40), dist:17.4,  mass:1.50, temp:6681.,  vsini:1.8,   ldm:0.723, lde:0.010, uld:0.250, pa:0.0, calfor:'none'}
hd141891 = {name:'hd_141891', spectrum:'ukf1v.dat',  ra:TEN(15,55,09), dec:TEN(-63,25,51), dist:12.4,  mass:1.50, temp:7182.,  vsini:92.0,  ldm:1.433, lde:0.020, uld:0.225, pa:0.0, calfor:'none'}
hd149661 = {name:'hd_149661', spectrum:'ukk2v.dat',  ra:TEN(16,36,21), dec:TEN(-02,19,29), dist:9.75,  mass:0.84, temp:5254.,  vsini:2.2,   ldm:0.776, lde:0.011, uld:0.300, pa:0.0, calfor:'none'}
hd152391 = {name:'hd_152391', spectrum:'ukg8v.dat',  ra:TEN(16,52,58), dec:TEN(-00,01,35), dist:17.25, mass:0.83, temp:5464.,  vsini:3.0,   ldm:0.487, lde:0.008, uld:0.327, pa:0.0, calfor:'none'}
hd160032 = {name:'hd_160032', spectrum:'ukf4v.dat',  ra:TEN(17,40,24), dec:TEN(-49,24,56), dist:21.5,  mass:1.50, temp:6579.,  vsini:16.3,  ldm:0.662, lde:0.094, uld:0.257, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd160915 = {name:'hd_160915', spectrum:'ukf5v.dat',  ra:TEN(17,43,26), dec:TEN(-21,41,00), dist:17.65, mass:1.50, temp:6222.,  vsini:12.4,  ldm:0.655, lde:0.088, uld:0.272, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd164259 = {name:'hd_164259', spectrum:'ukf2v.dat',  ra:TEN(18,00,29), dec:TEN(-03,41,25), dist:23.55, mass:1.44, temp:6709.,  vsini:69.3,  ldm:0.716, lde:0.010, uld:0.245, pa:0.0, calfor:'none'}
hd165777 = {name:'hd_165777', spectrum:'uka3v.dat',  ra:TEN(18,07,21), dec:TEN(+09,34,00), dist:26.63, mass:2.0,  temp:8260.,  vsini:65.0,  ldm:0.717, lde:0.012, uld:0.214, pa:0.0, calfor:'none'}
hd172555 = {name:'hd_172555', spectrum:'uka7v.dat',  ra:TEN(18,45,27), dec:TEN(-64,52,17), dist:28.55, mass:2.00, temp:7816.,  vsini:116.4, ldm:0.494, lde:0.009, uld:0.240, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd178253 = {name:'hd_178253', spectrum:'uka2v.dat',  ra:TEN(19,09,28), dec:TEN(-37,54,16), dist:38.43, mass:2.28, temp:8927.,  vsini:20.5,  ldm:0.514, lde:0.067, uld:0.191, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd182572 = {name:'hd_182572', spectrum:'ukg8v.dat',  ra:TEN(19,24,58), dec:TEN(+11,56,40), dist:15.18, mass:0.83, temp:5594.,  vsini:2.2,   ldm:0.858, lde:0.018, uld:0.327, pa:0.0, calfor:'none'}
hd188228 = {name:'hd_188228', spectrum:'uka0v.dat',  ra:TEN(20,00,36), dec:TEN(-72,54,37), dist:32.22, mass:2.9,  temp:10418., vsini:89.7,  ldm:0.587, lde:0.072, uld:0.172, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd192425 = {name:'hd_192425', spectrum:'uka2v.dat',  ra:TEN(20,14,17), dec:TEN(+15,11,51), dist:45.98, mass:2.28, temp:8996.,  vsini:180.0, ldm:0.378, lde:0.005, uld:0.191, pa:0.0, calfor:'none'}
hd195627 = {name:'hd_195627', spectrum:'ukf0v.dat',  ra:TEN(20,35,35), dec:TEN(-60,34,54), dist:27.79, mass:1.76, temp:7211.,  vsini:122.0, ldm:0.578, lde:0.065, uld:0.225, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd197157 = {name:'hd_197157', spectrum:'ukf0v.dat',  ra:TEN(20,44,02), dec:TEN(-51,55,15), dist:24.17, mass:1.76, temp:7448.,  vsini:150.0, ldm:0.639, lde:0.082, uld:0.246, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd197692 = {name:'hd_197692', spectrum:'ukf5v.dat',  ra:TEN(20,46,06), dec:TEN(-25,16,15), dist:14.7,  mass:1.50, temp:6605.,  vsini:41.7,  ldm:0.949, lde:0.119, uld:0.263, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd202730 = {name:'hd_202730', spectrum:'uka5v.dat',  ra:TEN(21,19,52), dec:TEN(-53,26,58), dist:30.28, mass:1.75, temp:7870.,  vsini:135.6, ldm:0.518, lde:0.008, uld:0.234, pa:0.0, calfor:'none'}
hd203608 = {name:'hd_203608', spectrum:'ukf9v.dat',  ra:TEN(21,26,27), dec:TEN(-65,21,58), dist:9.26,  mass:1.37, temp:6026.,  vsini:3.7,   ldm:1.078, lde:0.015, uld:0.293, pa:0.0, calfor:'none'}
hd206860 = {name:'hd_206860', spectrum:'ukg0v.dat',  ra:TEN(21,44,31), dec:TEN(+14,46,19), dist:17.89, mass:1.01, temp:5929.,  vsini:10.4,  ldm:0.514, lde:0.010, uld:0.293, pa:0.0, calfor:'none'}
hd207129 = {name:'hd_207129', spectrum:'ukg2v.dat',  ra:TEN(21,48,16), dec:TEN(-47,18,13), dist:15.99, mass:1.00, temp:5920.,  vsini:3.5,   ldm:0.627, lde:0.009, uld:0.293, pa:0.0, calfor:'none'}
hd210049 = {name:'hd_210049', spectrum:'uka2v.dat',  ra:TEN(22,08,23), dec:TEN(-32,59,18), dist:41.65, mass:2.28, temp:8937.,  vsini:307.7, ldm:0.455, lde:0.006, uld:0.195, pa:0.0, calfor:'none'}
hd210302 = {name:'hd_210302', spectrum:'ukf6v.dat',  ra:TEN(22,10,09), dec:TEN(-32,32,54), dist:18.28, mass:1.24, temp:6348.,  vsini:14.0,  ldm:0.705, lde:0.010, uld:0.257, pa:0.0, calfor:'none'}
hd210418 = {name:'hd_210418', spectrum:'uka0v.dat',  ra:TEN(22,10,12), dec:TEN(+06,11,52), dist:28.30, mass:2.9,  temp:8705.,  vsini:144.0, ldm:0.734, lde:0.010, uld:0.196, pa:0.0, calfor:'none'}
hd213845 = {name:'hd_213845', spectrum:'ukf7v.dat',  ra:TEN(22,34,42), dec:TEN(-20,42,30), dist:22.68, mass:1.53, temp:6551.,  vsini:35.7,  ldm:0.523, lde:0.082, uld:0.257, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd215789 = {name:'hd_215789', spectrum:'uka2v.dat',  ra:TEN(22,48,33), dec:TEN(-51,19,01), dist:39.53, mass:2.28, temp:8625.,  vsini:253.1, ldm:0.899, lde:0.042, uld:0.205, pa:0.0, calfor:'none'}
hd219482 = {name:'hd_219482', spectrum:'ukf6v.dat',  ra:TEN(23,16,58), dec:TEN(-62,00,04), dist:20.54, mass:1.24, temp:6318.,  vsini:7.5,   ldm:0.527, lde:0.007, uld:0.257, pa:0.0, calfor:'none'}  ; K mag from 2mass with rel. large uncertainty
hd219571 = {name:'hd_219571', spectrum:'ukf5v.dat',  ra:TEN(23,17,26), dec:TEN(-58,14,09), dist:23.06, mass:1.50, temp:6534.,  vsini:79.4,  ldm:1.003, lde:0.014, uld:0.257, pa:0.0, calfor:'none'}
hd224392 = {name:'hd_224392', spectrum:'uka0v.dat',  ra:TEN(23,57,35), dec:TEN(-64,17,54), dist:47.44, mass:2.9,  temp:8799.,  vsini:20.8,  ldm:0.368, lde:0.005, uld:0.196, pa:0.0, calfor:'none'}

; MAIN SURVEY PIONIER
;hd142    = {name:'hd142',    spectrum:'ukf6v.dat',  ra:TEN(00,06,19), dec:TEN(-49,04,31), dist:25.71, mass:1.29, temp:6213.,  vsini:11.0,  ldm:0.519, lde:0.008, uld:0.272, pa:0.0, calfor:'none'}
;hd1581   = {name:'hd1581',   spectrum:'ukf8v.dat',  ra:TEN(00,20,04), dec:TEN(-64,52,29), dist:8.59,  mass:1.00, temp:5925.,  vsini:2.3,   ldm:1.149, lde:0.016, uld:0.293, pa:0.0, calfor:'none'}
;hd2262   = {name:'hd2262',   spectrum:'uka5v.dat',  ra:TEN(00,26,12), dec:TEN(-43,40,47), dist:23.81, mass:1.69, temp:7935.,  vsini:225.0, ldm:0.698, lde:0.010, uld:0.224, pa:0.0, calfor:'none'}
;hd3302   = {name:'hd3302',   spectrum:'ukf5v.dat',  ra:TEN(00,35,41), dec:TEN(-48,00,03), dist:34.82, mass:1.48, temp:6696.,  vsini:17.8,  ldm:0.493, lde:0.008, uld:0.245, pa:0.0, calfor:'none'}
;hd3823   = {name:'hd3823',   spectrum:'ukg0v.dat',  ra:TEN(00,40,26), dec:TEN(-59,27,17), dist:24.96, mass:1.01, temp:6015.,  vsini:2.3,   ldm:0.532, lde:0.007, uld:0.293, pa:0.0, calfor:'none'}
;hd4150   = {name:'hd4150',   spectrum:'uka0v.dat',  ra:TEN(00,43,21), dec:TEN(-57,27,47), dist:75.53, mass:1.93, temp:9810.,  vsini:133.1, ldm:0.457, lde:0.007, uld:0.189, pa:0.0, calfor:'none'}
;hd7570   = {name:'hd7570',   spectrum:'ukf8v.dat',  ra:TEN(01,15,11), dec:TEN(-45,31,54), dist:15.11, mass:1.23, temp:6018.,  vsini:4.3,   ldm:0.762, lde:0.010, uld:0.293, pa:0.0, calfor:'none'}
;hd7788   = {name:'hd7788',   spectrum:'ukf6v.dat',  ra:TEN(01,15,46), dec:TEN(-68,52,33), dist:20.96, mass:1.37, temp:6366.,  vsini:61.0,  ldm:0.739, lde:0.011, uld:0.257, pa:0.0, calfor:'none'}
;hd10647  = {name:'hd10647',  spectrum:'ukf8v.dat',  ra:TEN(01,42,29), dec:TEN(-53,44,27), dist:17.43, mass:1.20, temp:6150.,  vsini:5.5,   ldm:0.547, lde:0.072, uld:0.277, pa:0.0, calfor:'none'}
;hd11171  = {name:'hd11171',  spectrum:'ukf2iii.dat', ra:TEN(01,49,35), dec:TEN(-10,41,11), dist:23.19, mass:1.52, temp:5833., vsini:58.5,  ldm:0.626, lde:0.120, uld:0.310, pa:0.0, calfor:'none'}
;hd14412  = {name:'hd14412',  spectrum:'ukg8v.dat',  ra:TEN(02,18,59), dec:TEN(-25,56,44), dist:12.67, mass:1.03, temp:5388.,  vsini:0.0,   ldm:0.552, lde:0.007, uld:0.327, pa:0.0, calfor:'none'}
;hd15008  = {name:'hd15008',  spectrum:'uka3v.dat',  ra:TEN(02,21,45), dec:TEN(-68,39,34), dist:42.83, mass:2.25,  temp:9001.,  vsini:190.0, ldm:0.541, lde:0.007, uld:0.191, pa:0.0, calfor:'none'}
;hd15798  = {name:'hd15798',  spectrum:'ukf5v.dat',  ra:TEN(02,32,05), dec:TEN(-15,14,41), dist:26.70, mass:1.38, temp:6349.,  vsini:4.6,   ldm:0.770, lde:0.011, uld:0.272, pa:0.0, calfor:'none'}
;hd16555  = {name:'hd16555',  spectrum:'uka7v.dat',  ra:TEN(02,37,24), dec:TEN(-52,32,35), dist:45.56, mass:1.75, temp:7139.,  vsini:6.60,  ldm:0.468, lde:0.007, uld:0.225, pa:0.0, calfor:'none'}
;hd17051  = {name:'hd17051',  spectrum:'ukf8v.dat',  ra:TEN(02,42,33), dec:TEN(-50,48,01), dist:17.17, mass:1.27, temp:6142.,  vsini:6.5,   ldm:0.632, lde:0.009, uld:0.277, pa:0.0, calfor:'none'}
;hd17925  = {name:'hd17925',  spectrum:'ukk2v.dat',  ra:TEN(02,52,32), dec:TEN(-12,46,11), dist:10.35, mass:0.84, temp:5197.,  vsini:4.9,   ldm:0.725, lde:0.010, uld:0.345, pa:0.0, calfor:'none'}
;hd19107  = {name:'hd19107',  spectrum:'uka7v.dat',  ra:TEN(03,04,17), dec:TEN(-07,36,03), dist:41.56, mass:1.77, temp:7714.,  vsini:150.0, ldm:0.405, lde:0.006, uld:0.240, pa:0.0, calfor:'none'}
;hd20766  = {name:'hd20766',  spectrum:'ukg5v.dat',  ra:TEN(03,17,46), dec:TEN(-62,34,31), dist:12.01, mass:0.96, temp:5684.,  vsini:2.7,   ldm:0.676, lde:0.010, uld:0.310, pa:0.0, calfor:'none'}
;hd20794  = {name:'hd20794',  spectrum:'ukg8v.dat',  ra:TEN(03,19,56), dec:TEN(-43,04,11), dist:6.04,  mass:0.70, temp:5431.,  vsini:2.0,   ldm:1.295, lde:0.172, uld:0.327, pa:0.0, calfor:'none'}
;hd20807  = {name:'hd20807',  spectrum:'ukg0v.dat',  ra:TEN(03,18,13), dec:TEN(-62,30,23), dist:12.03, mass:0.95, temp:5821.,  vsini:2.7,   ldm:0.747, lde:0.011, uld:0.310, pa:0.0, calfor:'none'}
;hd22001  = {name:'hd22001',  spectrum:'ukf5v.dat',  ra:TEN(03,29,23), dec:TEN(-62,56,15), dist:21.68, mass:1.39, temp:6732.,  vsini:13.2,  ldm:0.704, lde:0.010, uld:0.245, pa:0.0, calfor:'none'}
;hd23249  = {name:'hd23249',  spectrum:'ukk1iv.dat', ra:TEN(03,43,15), dec:TEN(-09,45,48), dist:9.04,  mass:0.83, temp:5025.,  vsini:2.3,   ldm:1.810, lde:0.025, uld:0.361, pa:0.0, calfor:'none'}
;hd25457  = {name:'hd25457',  spectrum:'ukf6v.dat',  ra:TEN(04,02,37), dec:TEN(-00,16,08), dist:18.83, mass:1.21, temp:6270.,  vsini:18.0,  ldm:0.591, lde:0.011, uld:0.277, pa:0.0, calfor:'none'}
;hd28355  = {name:'hd28355',  spectrum:'uka7v.dat',  ra:TEN(04,28,50), dec:TEN(+13,02,51), dist:48.85, mass:1.98, temp:7835.,  vsini:105.0, ldm:0.439, lde:0.005, uld:0.234, pa:0.0, calfor:'none'}
;hd29388  = {name:'hd29388',  spectrum:'uka7v.dat',  ra:TEN(04,38,09), dec:TEN(+12,30,39), dist:47.08, mass:2.28, temp:8338.,  vsini:89.0,  ldm:0.560, lde:0.007, uld:0.214, pa:0.0, calfor:'none'}
;hd30495  = {name:'hd30495',  spectrum:'ukg2v.dat',  ra:TEN(04,47,36), dec:TEN(-16,54,04), dist:13.28, mass:1.04, temp:5836.,  vsini:3.6,   ldm:0.675, lde:0.013, uld:0.310, pa:0.0, calfor:'none'}
;hd31295  = {name:'hd31295',  spectrum:'uka0v.dat',  ra:TEN(04,54,54), dec:TEN(+10,09,03), dist:35.66, mass:2.08, temp:8896.,  vsini:11.7,  ldm:0.448, lde:0.011, uld:0.191, pa:0.0, calfor:'none'}
;hd31925  = {name:'hd31925',  spectrum:'ukf5v.dat',  ra:TEN(04,59,01), dec:TEN(-16,22,34), dist:40.47, mass:1.44, temp:6496.,  vsini:7.2,   ldm:0.515, lde:0.007, uld:0.257, pa:0.0, calfor:'none'}
;hd33111  = {name:'hd33111',  spectrum:'uka3iii.dat', ra:TEN(05,07,51), dec:TEN(-05,05,11), dist:27.40, mass:2.19, temp:8002., vsini:180.0, ldm:1.180, lde:0.020, uld:0.224, pa:0.0, calfor:'none'}
;hd33262  = {name:'hd33262',  spectrum:'ukf8v.dat',  ra:TEN(05,05,31), dec:TEN(-57,28,22), dist:11.65, mass:1.12, temp:6146.,  vsini:15.4,  ldm:0.879, lde:0.098, uld:0.277, pa:0.0, calfor:'none'}
;hd34721  = {name:'hd34721',  spectrum:'ukg0v.dat',  ra:TEN(05,18,50), dec:TEN(-18,07,50), dist:25.03, mass:1.13, temp:6000.,  vsini:4.4,   ldm:0.494, lde:0.007, uld:0.288, pa:0.0, calfor:'none'}
;hd38858  = {name:'hd38858',  spectrum:'ukg5v.dat',  ra:TEN(05,48,35), dec:TEN(-04,05,41), dist:15.18, mass:0.90, temp:5732.,  vsini:0.3,   ldm:0.576, lde:0.008, uld:0.310, pa:0.0, calfor:'none'}
;hd39060  = {name:'hd39060',  spectrum:'uka5v.dat',  ra:TEN(05,47,17), dec:TEN(-51,03,59), dist:19.44, mass:1.76, temp:8045.,  vsini:13.3,  ldm:0.766, lde:0.010, uld:0.224, pa:0.0, calfor:'none'}
;hd40307  = {name:'hd40307',  spectrum:'ukk2v.dat',  ra:TEN(05,54,04), dec:TEN(-60,01,24), dist:13.00, mass:0.77, temp:4881.,  vsini:1.6,   ldm:0.546, lde:0.007, uld:0.382, pa:0.0, calfor:'none'}
;hd43162  = {name:'hd43162',  spectrum:'ukg5v.dat',  ra:TEN(06,13,45), dec:TEN(-23,51,43), dist:16.72, mass:0.98, temp:5579.,  vsini:5.5,   ldm:0.496, lde:0.006, uld:0.327, pa:0.0, calfor:'none'}
;hd45184  = {name:'hd45184',  spectrum:'ukg2v.dat',  ra:TEN(06,24,44), dec:TEN(-28,46,48), dist:21.88, mass:1.03, temp:5820.,  vsini:2.5,   ldm:0.453, lde:0.006, uld:0.310, pa:0.0, calfor:'none'}
;hd53705  = {name:'hd53705',  spectrum:'ukg0v.dat',  ra:TEN(07,03,57), dec:TEN(-43,36,29), dist:16.52, mass:0.96, temp:5528.,  vsini:1.60,  ldm:0.623, lde:0.011, uld:0.327, pa:0.0, calfor:'none'}
;; hd56537  = {name:'hd56537',  spectrum:'uka3v.dat',  ra:TEN(07,18,06), dec:TEN(+16,32,25), dist:30.93, mass:2.19, temp:8395.,  vsini:154.0, ldm:0.651, lde:0.081, uld:0.204, pa:0.0, calfor:'none'}
;hd56537  = {name:'hd56537',  spectrum:'uka3v.dat',  ra:TEN(07,18,06), dec:TEN(+16,32,25), dist:30.93, mass:2.19, temp:8395.,  vsini:154.0, ldm:0.835, lde:0.013, uld:0.204, pa:0.0, calfor:'none'}
;hd69830  = {name:'hd69830',  spectrum:'ukg8v.dat',  ra:TEN(08,18,24), dec:TEN(-12,37,56), dist:12.49, mass:0.82, temp:5420.,  vsini:1.60,  ldm:0.656, lde:0.009, uld:0.327, pa:0.0, calfor:'none'}
;hd71155  = {name:'hd71155',  spectrum:'uka0v.dat',  ra:TEN(08,25,40), dec:TEN(-03,54,23), dist:37.51, mass:2.40, temp:9843.,  vsini:14.30, ldm:0.534, lde:0.007, uld:0.182, pa:0.0, calfor:'none'}
;hd72673  = {name:'hd72673',  spectrum:'ukg8v.dat',  ra:TEN(08,32,51), dec:TEN(-31,30,03), dist:12.21, mass:0.70, temp:5247.,  vsini:0.0,   ldm:0.597, lde:0.012, uld:0.345, pa:0.0, calfor:'none'}
;hd76151  = {name:'hd76151',  spectrum:'ukg2v.dat',  ra:TEN(08,54,18), dec:TEN(-05,26,04), dist:17.39, mass:1.04, temp:5763.,  vsini:4.0,   ldm:0.537, lde:0.009, uld:0.310, pa:0.0, calfor:'none'}
;hd76932  = {name:'hd76932',  spectrum:'ukg2v.dat',  ra:TEN(08,58,44), dec:TEN(-16,07,58), dist:21.03, mass:0.86, temp:5850.,  vsini:2.6,   ldm:0.561, lde:0.012, uld:0.305, pa:0.0, calfor:'none'}
;hd82434  = {name:'hd82434',  spectrum:'ukf5v.dat',  ra:TEN(09,30,42), dec:TEN(-40,28,00), dist:18.81, mass:1.63, temp:6867.,  vsini:156.0, ldm:1.065, lde:0.018, uld:0.245, pa:0.0, calfor:'none'}
;hd88955  = {name:'hd88955',  spectrum:'uka2v.dat',  ra:TEN(10,14,44), dec:TEN(-42,07,19), dist:31.08, mass:2.22, temp:8953.,  vsini:10.5,  ldm:0.597, lde:0.010, uld:0.191, pa:0.0, calfor:'none'}
;hd90132  = {name:'hd90132',  spectrum:'uka7v.dat',  ra:TEN(10,23,29), dec:TEN(-38,00,35), dist:41.46, mass:1.74, temp:7408.,  vsini:270.0, ldm:0.426, lde:0.006, uld:0.222, pa:0.0, calfor:'none'}
;hd91324  = {name:'hd91324',  spectrum:'ukf8v.dat',  ra:TEN(10,31,22), dec:TEN(-53,42,56), dist:21.81, mass:1.30, temp:6096.,  vsini:8.3,   ldm:0.793, lde:0.080, uld:0.288, pa:0.0, calfor:'none'}
;hd99211  = {name:'hd99211',  spectrum:'uka7v.dat',  ra:TEN(11,24,53), dec:TEN(-17,41,02), dist:25.24, mass:1.80, temp:7620.,  vsini:7.3,   ldm:0.705, lde:0.188, uld:0.246, pa:0.0, calfor:'none'}
;hd102365 = {name:'hd102365', spectrum:'ukg2v.dat',  ra:TEN(11,46,31), dec:TEN(-40,30,01), dist:9.22,  mass:0.82, temp:5575.,  vsini:0.7,   ldm:0.943, lde:0.014, uld:0.327, pa:0.0, calfor:'none'}
;hd104731 = {name:'hd104731', spectrum:'ukf5v.dat',  ra:TEN(12,03,40), dec:TEN(-42,26,03), dist:25.32, mass:1.31, temp:6385.,  vsini:15.9,  ldm:0.604, lde:0.008, uld:0.257, pa:0.0, calfor:'none'}
;hd108767 = {name:'hd108767', spectrum:'uka0iv.dat', ra:TEN(12,29,52), dec:TEN(-16,30,56), dist:26.63, mass:1.59, temp:10371., vsini:236.0, ldm:0.792, lde:0.011, uld:0.177, pa:0.0, calfor:'none'}
;hd109787 = {name:'hd109787', spectrum:'uka2v.dat',  ra:TEN(12,37,42), dec:TEN(-48,32,29), dist:40.2,  mass:2.44, temp:8762.,  vsini:296.8, ldm:0.577, lde:0.008, uld:0.196, pa:0.0, calfor:'none'}
;hd115617 = {name:'hd115617', spectrum:'ukg8v.dat',  ra:TEN(13,18,24), dec:TEN(-18,18,40), dist:8.58,  mass:0.88, temp:5557.,  vsini:3.9,   ldm:1.147, lde:0.129, uld:0.327, pa:0.0, calfor:'none'}
;hd120136 = {name:'hd120136', spectrum:'ukf5iv.dat', ra:TEN(13,47,16), dec:TEN(+17,27,25), dist:15.6,  mass:1.30, temp:6449.,  vsini:15.6,  ldm:0.856, lde:0.116, uld:0.263, pa:0.0, calfor:'none'}
;hd128898 = {name:'hd128898', spectrum:'uka7v.dat',  ra:TEN(14,42,30), dec:TEN(-64,58,30), dist:16.57, mass:1.81, temp:7426.,  vsini:13.2,  ldm:1.005, lde:0.025, uld:0.246, pa:0.0, calfor:'none'}
;hd129502 = {name:'hd129502', spectrum:'ukf2v.dat',  ra:TEN(14,43,04), dec:TEN(-05,39,30), dist:18.3,  mass:1.58, temp:6640.,  vsini:47.0,  ldm:1.027, lde:0.014, uld:0.245, pa:0.0, calfor:'none'}
;hd130109 = {name:'hd130109', spectrum:'uka0v.dat',  ra:TEN(14,46,15), dec:TEN(+01,53,34), dist:41.24, mass:2.58, temp:9683.,  vsini:258.0, ldm:0.613, lde:0.021, uld:0.178, pa:0.0, calfor:'none'}
;hd134083 = {name:'hd134083', spectrum:'ukf5v.dat',  ra:TEN(15,07,19), dec:TEN(+24,52,09), dist:19.55, mass:1.34, temp:6567.,  vsini:44.4,  ldm:0.661, lde:0.007, uld:0.263, pa:0.0, calfor:'none'}
;hd135379 = {name:'hd135379', spectrum:'uka3v.dat',  ra:TEN(15,17,30), dec:TEN(-58,48,04), dist:30.55, mass:1.89, temp:8695.,  vsini:68.5,  ldm:0.570, lde:0.008, uld:0.196, pa:0.0, calfor:'none'}
;hd136202 = {name:'hd136202', spectrum:'ukf8v.dat',  ra:TEN(15,19,19), dec:TEN(+01,45,55), dist:25.38, mass:1.40, temp:6129.,  vsini:4.8,   ldm:0.623, lde:0.085, uld:0.272, pa:0.0, calfor:'none'}
;hd139664 = {name:'hd139664', spectrum:'ukf5v.dat',  ra:TEN(15,41,11), dec:TEN(-44,39,40), dist:17.4,  mass:1.33, temp:6681.,  vsini:1.8,   ldm:0.723, lde:0.010, uld:0.250, pa:0.0, calfor:'none'}
;hd141891 = {name:'hd141891', spectrum:'ukf0v.dat',  ra:TEN(15,55,09), dec:TEN(-63,25,51), dist:12.4,  mass:1.71, temp:7182.,  vsini:92.0,  ldm:1.433, lde:0.020, uld:0.225, pa:0.0, calfor:'none'}
;hd149661 = {name:'hd149661', spectrum:'ukk2v.dat',  ra:TEN(16,36,21), dec:TEN(-02,19,29), dist:9.75,  mass:0.86, temp:5254.,  vsini:2.2,   ldm:0.776, lde:0.011, uld:0.300, pa:0.0, calfor:'none'}
;hd152391 = {name:'hd152391', spectrum:'ukg8v.dat',  ra:TEN(16,52,58), dec:TEN(-00,01,35), dist:17.25, mass:0.89, temp:5464.,  vsini:3.0,   ldm:0.487, lde:0.008, uld:0.327, pa:0.0, calfor:'none'}
;hd160032 = {name:'hd160032', spectrum:'ukf5v.dat',  ra:TEN(17,40,24), dec:TEN(-49,24,56), dist:21.5,  mass:1.32, temp:6579.,  vsini:16.3,  ldm:0.662, lde:0.094, uld:0.257, pa:0.0, calfor:'none'}
;hd160915 = {name:'hd160915', spectrum:'ukf5v.dat',  ra:TEN(17,43,26), dec:TEN(-21,41,00), dist:17.65, mass:1.18, temp:6222.,  vsini:12.4,  ldm:0.655, lde:0.088, uld:0.272, pa:0.0, calfor:'none'}
;hd164259 = {name:'hd164259', spectrum:'ukf2v.dat',  ra:TEN(18,00,29), dec:TEN(-03,41,25), dist:23.55, mass:1.50, temp:6709.,  vsini:69.3,  ldm:0.716, lde:0.010, uld:0.245, pa:0.0, calfor:'none'}
;hd165777 = {name:'hd165777', spectrum:'uka3v.dat',  ra:TEN(18,07,21), dec:TEN(+09,34,00), dist:26.63, mass:1.99, temp:8260.,  vsini:65.0,  ldm:0.717, lde:0.012, uld:0.214, pa:0.0, calfor:'none'}
;hd172555 = {name:'hd172555', spectrum:'uka7v.dat',  ra:TEN(18,45,27), dec:TEN(-64,52,17), dist:28.55, mass:1.67, temp:7816.,  vsini:116.4, ldm:0.494, lde:0.009, uld:0.240, pa:0.0, calfor:'none'}
;hd178253 = {name:'hd178253', spectrum:'uka2v.dat',  ra:TEN(19,09,28), dec:TEN(-37,54,16), dist:38.43, mass:2.22, temp:8927.,  vsini:20.5,  ldm:0.514, lde:0.067, uld:0.191, pa:0.0, calfor:'none'}
;hd182572 = {name:'hd182572', spectrum:'ukg8v.dat',  ra:TEN(19,24,58), dec:TEN(+11,56,40), dist:15.18, mass:1.06, temp:5594.,  vsini:2.2,   ldm:0.858, lde:0.018, uld:0.327, pa:0.0, calfor:'none'}
;hd188228 = {name:'hd188228', spectrum:'uka0v.dat',  ra:TEN(20,00,36), dec:TEN(-72,54,37), dist:32.22, mass:2.34, temp:10418., vsini:89.7,  ldm:0.587, lde:0.072, uld:0.172, pa:0.0, calfor:'none'}
;hd192425 = {name:'hd192425', spectrum:'uka2v.dat',  ra:TEN(20,14,17), dec:TEN(+15,11,51), dist:45.98, mass:2.06, temp:8996.,  vsini:180.0, ldm:0.378, lde:0.005, uld:0.191, pa:0.0, calfor:'none'}
;hd195627 = {name:'hd195627', spectrum:'ukf0v.dat',  ra:TEN(20,35,35), dec:TEN(-60,34,54), dist:27.79, mass:1.59, temp:7211.,  vsini:122.0, ldm:0.578, lde:0.065, uld:0.225, pa:0.0, calfor:'none'}
;hd197157 = {name:'hd197157', spectrum:'ukf0v.dat',  ra:TEN(20,44,02), dec:TEN(-51,55,15), dist:24.17, mass:1.58, temp:7448.,  vsini:150.0, ldm:0.639, lde:0.082, uld:0.246, pa:0.0, calfor:'none'}
;hd197692 = {name:'hd197692', spectrum:'ukf5v.dat',  ra:TEN(20,46,06), dec:TEN(-25,16,15), dist:14.7,  mass:1.35, temp:6605.,  vsini:41.7,  ldm:0.949, lde:0.119, uld:0.263, pa:0.0, calfor:'none'}
;hd202730 = {name:'hd202730', spectrum:'uka5v.dat',  ra:TEN(21,19,52), dec:TEN(-53,26,58), dist:30.28, mass:1.88, temp:7870.,  vsini:135.6, ldm:0.518, lde:0.008, uld:0.234, pa:0.0, calfor:'none'}
;hd203608 = {name:'hd203608', spectrum:'ukf8v.dat',  ra:TEN(21,26,27), dec:TEN(-65,21,58), dist:9.26,  mass:0.85, temp:6026.,  vsini:3.7,   ldm:1.078, lde:0.015, uld:0.293, pa:0.0, calfor:'none'}
;hd206860 = {name:'hd206860', spectrum:'ukg0v.dat',  ra:TEN(21,44,31), dec:TEN(+14,46,19), dist:17.89, mass:1.04, temp:5929.,  vsini:10.4,  ldm:0.514, lde:0.010, uld:0.293, pa:0.0, calfor:'none'}
;hd207129 = {name:'hd207129', spectrum:'ukg2v.dat',  ra:TEN(21,48,16), dec:TEN(-47,18,13), dist:15.99, mass:1.05, temp:5920.,  vsini:3.5,   ldm:0.627, lde:0.009, uld:0.293, pa:0.0, calfor:'none'}
;hd210049 = {name:'hd210049', spectrum:'uka2v.dat',  ra:TEN(22,08,23), dec:TEN(-32,59,18), dist:41.65, mass:2.15, temp:8937.,  vsini:307.7, ldm:0.455, lde:0.006, uld:0.195, pa:0.0, calfor:'none'}
;hd210277 = {name:'hd210277', spectrum:'ukg0v.dat',  ra:TEN(22,09,30), dec:TEN(-07,32,55), dist:21.56, mass:0.91, temp:5544.,  vsini:1.8,   ldm:0.488, lde:0.007, uld:0.327, pa:0.0, calfor:'none'}
;hd210302 = {name:'hd210302', spectrum:'ukf6v.dat',  ra:TEN(22,10,09), dec:TEN(-32,32,54), dist:18.28, mass:1.30, temp:6348.,  vsini:14.0,  ldm:0.705, lde:0.010, uld:0.257, pa:0.0, calfor:'none'}
;hd210418 = {name:'hd210418', spectrum:'uka0v.dat',  ra:TEN(22,10,12), dec:TEN(+06,11,52), dist:28.30, mass:2.29, temp:8705.,  vsini:144.0, ldm:0.734, lde:0.010, uld:0.196, pa:0.0, calfor:'none'}
;hd213845 = {name:'hd213845', spectrum:'ukf8v.dat',  ra:TEN(22,34,42), dec:TEN(-20,42,30), dist:22.68, mass:1.37, temp:6551.,  vsini:35.7,  ldm:0.523, lde:0.082, uld:0.257, pa:0.0, calfor:'none'}
;hd214953 = {name:'hd214953', spectrum:'ukg0v.dat',  ra:TEN(22,42,37), dec:TEN(-47,12,39), dist:23.64, mass:1.13, temp:6021.,  vsini:4.5,   ldm:0.526, lde:0.007, uld:0.293, pa:0.0, calfor:'none'}
;hd215648 = {name:'hd215648', spectrum:'ukf8v.dat',  ra:TEN(22,46,42), dec:TEN(+12,10,22), dist:16.30, mass:1.33, temp:6076.,  vsini:6.7,   ldm:1.089, lde:0.041, uld:0.288, pa:0.0, calfor:'none'}
;hd215789 = {name:'hd215789', spectrum:'uka2v.dat',  ra:TEN(22,48,33), dec:TEN(-51,19,01), dist:39.53, mass:2.37, temp:8625.,  vsini:253.1, ldm:0.899, lde:0.042, uld:0.205, pa:0.0, calfor:'none'}
;hd216435 = {name:'hd216435', spectrum:'ukg0v.dat',  ra:TEN(22,53,38), dec:TEN(-48,35,54), dist:32.62, mass:1.33, temp:5918.,  vsini:5.7,   ldm:0.472, lde:0.008, uld:0.288, pa:0.0, calfor:'none'}
;hd219482 = {name:'hd219482', spectrum:'ukf6v.dat',  ra:TEN(23,16,58), dec:TEN(-62,00,04), dist:20.54, mass:1.05, temp:6318.,  vsini:7.5,   ldm:0.527, lde:0.007, uld:0.257, pa:0.0, calfor:'none'}
;hd219571 = {name:'hd219571', spectrum:'ukf5v.dat',  ra:TEN(23,17,26), dec:TEN(-58,14,09), dist:23.06, mass:1.63, temp:6534.,  vsini:79.4,  ldm:1.003, lde:0.014, uld:0.257, pa:0.0, calfor:'none'}
;hd224392 = {name:'hd224392', spectrum:'uka0v.dat',  ra:TEN(23,57,35), dec:TEN(-64,17,54), dist:47.44, mass:2.10, temp:8799.,  vsini:20.8,  ldm:0.368, lde:0.005, uld:0.196, pa:0.0, calfor:'none'}

;hd174429 = {name:'hip92680', spectrum:'ukg8v.dat',  ra:TEN(18,53,06), dec:TEN(-50,10,50), dist:51.49, mass:0.97, temp:5150.,  vsini:66.20, ldm:0.250, lde:0.005, uld:0.345, pa:0.0, calfor:'none'}

; PIONIER H+K survey (H data)
;hd2262   = {name:'hd2262',   spectrum:'uka5v.dat',  ra:TEN(00,26,12), dec:TEN(-43,40,47), dist:23.81, mass:1.69, temp:7935.,  vsini:225.0, ldm:0.698, lde:0.010, uld:0.224, pa:0.0, calfor:'none'}
;hd7788   = {name:'hd7788',   spectrum:'ukf6v.dat',  ra:TEN(01,15,46), dec:TEN(-68,52,33), dist:20.96, mass:1.37, temp:6366.,  vsini:61.0,  ldm:0.739, lde:0.011, uld:0.257, pa:0.0, calfor:'none'}
;hd10700  = {name:'hd10700',  spectrum:'ukg8v.dat',  ra:TEN(01,44,04), dec:TEN(-15,56,15), dist:3.65,  mass:0.76, temp:5390.,  vsini:1.8,   ldm:2.074, lde:0.030, uld:0.327, pa:0.0, calfor:'none'}
;hd20794  = {name:'hd20794',  spectrum:'ukg8v.dat',  ra:TEN(03,19,56), dec:TEN(-43,04,11), dist:6.04,  mass:0.70, temp:5431.,  vsini:2.0,   ldm:1.295, lde:0.172, uld:0.327, pa:0.0, calfor:'none'}
;hd22484  = {name:'hd22484',  spectrum:'ukf8v.dat',  ra:TEN(03,36,52), dec:TEN(+00,24,05), dist:13.96, mass:1.14, temp:5971.,  vsini:4.0,   ldm:1.102, lde:0.018, uld:0.288, pa:0.0, calfor:'none'}
;hd39060  = {name:'hd39060',  spectrum:'uka5v.dat',  ra:TEN(05,47,17), dec:TEN(-51,03,59), dist:19.44, mass:1.76, temp:8045.,  vsini:13.3,  ldm:0.766, lde:0.010, uld:0.224, pa:0.0, calfor:'none'}
;hd158643 = {name:'hd158643', spectrum:'uka0v.dat',  ra:TEN(17,31,25), dec:TEN(-23,57,45), dist:124.4, mass:3.27, temp:9840.,  vsini:228.0, ldm:0.491, lde:0.006, uld:0.199, pa:0.0, calfor:'none'}
;hd161868 = {name:'hd161868', spectrum:'uka0v.dat',  ra:TEN(17,47,54), dec:TEN(+02,42,26), dist:31.52, mass:2.28, temp:9311.,  vsini:210.0, ldm:0.623, lde:0.013, uld:0.187, pa:0.0, calfor:'none'}
;hd172555 = {name:'hd172555', spectrum:'uka7v.dat',  ra:TEN(18,45,27), dec:TEN(-64,52,17), dist:28.55, mass:1.67, temp:7816.,  vsini:116.4, ldm:0.494, lde:0.009, uld:0.240, pa:0.0, calfor:'none'}
;hd173667 = {name:'hd173667', spectrum:'ukf6v.dat',  ra:TEN(18,45,40), dec:TEN(+20,32,47), dist:19.21, mass:1.46, temp:6467.,  vsini:18.0,  ldm:0.987, lde:0.013, uld:0.257, pa:0.0, calfor:'none'}
;hd177724 = {name:'hd177724', spectrum:'uka0iv.dat', ra:TEN(19,05,25), dec:TEN(+13,51,49), dist:25.46, mass:2.43, temp:9528.,  vsini:317.0, ldm:0.878, lde:0.012, uld:0.184, pa:0.0, calfor:'none'}
;hd197481 = {name:'hd197481', spectrum:'ukm1v.dat',  ra:TEN(20,45,10), dec:TEN(-31,20,27), dist:3.57,  mass:0.47, temp:3493.,  vsini:8.2,   ldm:0.826, lde:0.015, uld:0.247, pa:0.0, calfor:'none'}
;hd210302 = {name:'hd210302', spectrum:'ukf6v.dat',  ra:TEN(22,10,09), dec:TEN(-32,32,54), dist:18.28, mass:1.30, temp:6348.,  vsini:14.0,  ldm:0.705, lde:0.010, uld:0.257, pa:0.0, calfor:'none'}
;hd224392 = {name:'hd224392', spectrum:'uka0v.dat',  ra:TEN(23,57,35), dec:TEN(-64,17,54), dist:47.44, mass:2.10, temp:8799.,  vsini:20.8,  ldm:0.368, lde:0.005, uld:0.196, pa:0.0, calfor:'none'}

; PIONIER H+K survey (K data)
hd2262   = {name:'hd2262',   spectrum:'uka5v.dat',  ra:TEN(00,26,12), dec:TEN(-43,40,47), dist:23.81, mass:1.69, temp:7935.,  vsini:225.0, ldm:0.698, lde:0.010, uld:0.195, pa:0.0, calfor:'none'}
hd4150   = {name:'hd4150',   spectrum:'uka0v.dat',  ra:TEN(00,43,21), dec:TEN(-57,27,47), dist:75.53, mass:1.93, temp:9810.,  vsini:133.1, ldm:0.457, lde:0.007, uld:0.160, pa:0.0, calfor:'none'}
hd7788   = {name:'hd7788',   spectrum:'ukf6v.dat',  ra:TEN(01,15,46), dec:TEN(-68,52,33), dist:20.96, mass:1.37, temp:6366.,  vsini:61.0,  ldm:0.739, lde:0.011, uld:0.226, pa:0.0, calfor:'none'}
hd15798  = {name:'hd15798',  spectrum:'ukf5v.dat',  ra:TEN(02,32,05), dec:TEN(-15,14,41), dist:26.70, mass:1.38, temp:6349.,  vsini:4.6,   ldm:0.770, lde:0.011, uld:0.239, pa:0.0, calfor:'none'}
hd14412  = {name:'hd14412',  spectrum:'ukg8v.dat',  ra:TEN(02,18,59), dec:TEN(-25,56,44), dist:12.67, mass:1.03, temp:5388.,  vsini:0.0,   ldm:0.552, lde:0.007, uld:0.284, pa:0.0, calfor:'none'}
hd20794  = {name:'hd20794',  spectrum:'ukg8v.dat',  ra:TEN(03,19,56), dec:TEN(-43,04,11), dist:6.04,  mass:0.70, temp:5431.,  vsini:2.0,   ldm:1.295, lde:0.172, uld:0.284, pa:0.0, calfor:'none'}
hd22484  = {name:'hd22484',  spectrum:'ukf8v.dat',  ra:TEN(03,36,52), dec:TEN(+00,24,05), dist:13.96, mass:1.14, temp:5971.,  vsini:4.0,   ldm:1.102, lde:0.018, uld:0.252, pa:0.0, calfor:'none'}
hd39060  = {name:'hd39060',  spectrum:'uka5v.dat',  ra:TEN(05,47,17), dec:TEN(-51,03,59), dist:19.44, mass:1.76, temp:8045.,  vsini:13.3,  ldm:0.766, lde:0.010, uld:0.195, pa:0.0, calfor:'none'}
hd158643 = {name:'hd158643', spectrum:'uka0v.dat',  ra:TEN(17,31,25), dec:TEN(-23,57,45), dist:124.4, mass:3.27, temp:9840.,  vsini:228.0, ldm:0.491, lde:0.006, uld:0.167, pa:0.0, calfor:'none'}
hd161868 = {name:'hd161868', spectrum:'uka0v.dat',  ra:TEN(17,47,54), dec:TEN(+02,42,26), dist:31.52, mass:2.28, temp:9311.,  vsini:210.0, ldm:0.623, lde:0.013, uld:0.162, pa:0.0, calfor:'none'}
hd172555 = {name:'hd172555', spectrum:'uka7v.dat',  ra:TEN(18,45,27), dec:TEN(-64,52,17), dist:28.55, mass:1.67, temp:7816.,  vsini:116.4, ldm:0.494, lde:0.009, uld:0.207, pa:0.0, calfor:'none'}
hd173667 = {name:'hd173667', spectrum:'ukf6v.dat',  ra:TEN(18,45,40), dec:TEN(+20,32,47), dist:19.21, mass:1.46, temp:6467.,  vsini:18.0,  ldm:0.987, lde:0.013, uld:0.226, pa:0.0, calfor:'none'}
hd177724 = {name:'hd177724', spectrum:'uka0iv.dat', ra:TEN(19,05,25), dec:TEN(+13,51,49), dist:25.46, mass:2.43, temp:9528.,  vsini:317.0, ldm:0.878, lde:0.012, uld:0.159, pa:0.0, calfor:'none'}
hd197481 = {name:'hd197481', spectrum:'ukm1v.dat',  ra:TEN(20,45,10), dec:TEN(-31,20,27), dist:3.57,  mass:0.47, temp:3493.,  vsini:8.2,   ldm:0.826, lde:0.015, uld:0.207, pa:0.0, calfor:'none'}
hd210302 = {name:'hd210302', spectrum:'ukf6v.dat',  ra:TEN(22,10,09), dec:TEN(-32,32,54), dist:18.28, mass:1.30, temp:6348.,  vsini:14.0,  ldm:0.705, lde:0.010, uld:0.226, pa:0.0, calfor:'none'}
hd216956 = {name:'hd216956', spectrum:'uka3v.dat',  ra:TEN(22,57,39), dec:TEN(-29,37,20), dist:7.70,  mass:2.0,  temp:8610.,  vsini:93.0,  ldm:2.069, lde:0.030, uld:0.178, pa:0.0, calfor:'none'}
hd224392 = {name:'hd224392', spectrum:'uka0v.dat',  ra:TEN(23,57,35), dec:TEN(-64,17,54), dist:47.44, mass:2.10, temp:8799.,  vsini:20.8,  ldm:0.368, lde:0.005, uld:0.171, pa:0.0, calfor:'none'}

; LBTI HOSTS survey 
; *****************

; Note: Effective temperature and distance are important to compute the EEID, which determines the right photometric aperture to use.
; RA and DEC are important to compute the hour angle. 
; Limb-darkening diameter and error are important to compute the geometric leakage.

; Engineering targets
alp_aql = {name:'alf_aql' , spectrum:'uka7v.dat', ra:TEN(19,50,47.00), dec:TEN(08,52,05.96), dist:5.13, mass:1.68, temp:7650., vsini:217.0, ldm:3.35, lde:0.13, uld: 0.00, pa:0.0, calfor:'none'} ; ldm from surface brightness relations (Absil et al. 2013)
alp_boo = {name:'alp_boo' , spectrum:'ukk2iii.dat', ra:TEN(14,15,39.67), dec:TEN(19,10,56.67), dist:11.3, mass:0., temp:4380., vsini:0.0, ldm:20.2, lde:0.08, uld: 0.00, pa:0.0, calfor:'none'}   ; ldm from CHARM
alp_ori = {name:'alp_ori' , spectrum:'ukm2i.dat', ra:TEN(0,0,0), dec:TEN(0,0,0), dist:0.0, mass:0.0, temp:3520. , vsini:0.0, ldm:45.2, lde:0.2, uld: 0.00, pa:0.0, calfor:'none'}
alp_ser = {name:'alp_ser' , spectrum:'ukk2iii.dat', ra:TEN(15,44,16.1), dec:TEN(+06,25,32.3), dist:0.0, mass:0.0, temp:3520. , vsini:0.0, ldm:4.79, lde:0.05, uld: 0.00, pa:0.0, calfor:'none'}   ; ldm from CHARM
alp_tau = {name:'alp_tau' , spectrum:'ukk5iii.dat', ra:TEN(0,0,0), dec:TEN(0,0,0), dist:0.0, mass:0.0, temp:3950., vsini:0.0, ldm:20.58, lde:0.03, uld: 0.00, pa:0.0, calfor:'none'}
bet_and2 = {name:'bet_And', spectrum:'ukm0iii.dat', ra:TEN(01,09,43.92388), dec:TEN(35,37,14.0075), dist:60.53, mass:0., temp:3650., vsini:0.0, ldm:13.749, lde:0.137, uld: 0.00, pa:0.0, calfor:'none'} ; LDD from Morzenwich 2003
bet_and = {name:'betaAnd',  spectrum:'ukm0iii.dat', ra:TEN(01,09,43.92388), dec:TEN(35,37,14.0075), dist:60.53, mass:0., temp:3650., vsini:0.0, ldm:13.749, lde:0.137, uld: 0.00, pa:0.0, calfor:'gamPer'} ; LDD from Morzenwich 2003
bet_boo = {name:'betaBoo' , spectrum:'ukg0iv.dat', ra:TEN(13,54,41.07892), dec:TEN(18,23,51.7946), dist:11.3, mass:0.0, temp:6100., vsini:0.0, ldm:2.21, lde:0.09, uld: 0.00, pa:0.0, calfor:'none'} ; All junk
bet_gem = {name:'bet_gem' , spectrum:'ukk0iii.dat', ra:TEN(07,45,18.94987), dec:TEN(28,01,34.3160), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:7.9, lde:0.31, uld: 0.00, pa:0.0, calfor:'none'}
bet_peg = {name:'bet_peg' , spectrum:'ukm2iii.dat', ra:TEN(0,0,0), dec:TEN(0,0,0), dist:0.0, mass:0.0, temp:3689. , vsini:0.0, ldm:14.6, lde:0.7, uld: 0.00, pa:0.0, calfor:'none'}
eps_vir = {name:'eps_vir' , spectrum:'ukg8iii.dat', ra:TEN(13,02,10.6), dec:TEN(10,57,32.9), dist:0.0, mass:0.0, temp:3689. , vsini:0.0, ldm:3.26, lde:0.7, uld: 0.00, pa:0.0, calfor:'none'}      ; ldm from SBR dirty
eps_leo = {name:'eps_leo' , spectrum:'ukg2i.dat', ra:TEN(09,45,51.07330), dec:TEN(23,46,27.3208), dist:76.9, mass:0.0, temp:5248., vsini:0.0, ldm:2.56, lde:0.03, uld: 0.00, pa:0.0, calfor:'none'} ; lde and temp garbage
eta_boo = {name:'eta_boo' , spectrum:'ukg0iv.dat', ra:TEN(13,54,41.07892), dec:TEN(18,23,51.7946), dist:11.3, mass:0.0, temp:6100., vsini:0.0, ldm:2.21, lde:0.09, uld: 0.00, pa:0.0, calfor:'none'}
gam_per = {name:'gam Per' , spectrum:'ukg8iii.dat', ra:TEN(03,04,47.8), dec:TEN(53,30,23.2), dist:0.0, mass:0.0, temp:8000.0, vsini:0., ldm:3.56, lde:0.013, uld:0., pa:0., calfor:'none'}        ; Bogus!!
hp_boo  = {name:'hp_boo'  , spectrum:'ukg0v.dat',   ra:TEN(14,50,15.81), dec:TEN(23,54,42.64), dist:18.2, mass:0., temp:5900., vsini:0.0, ldm:0.548, lde:0.010, uld: 0.00, pa:0.0, calfor:'none'} ; ldm from PTI (van Belle et al. 2008)
lam_aur = {name:'lamaur'  , spectrum:'ukg0v.dat',   ra:TEN(14,50,15.81), dec:TEN(23,54,42.64), dist:18.2, mass:0., temp:5900., vsini:0.0, ldm:0.548, lde:0.010, uld: 0.00, pa:0.0, calfor:'none'} ; ldm from PTI (van Belle et al. 2008)
lam_aur2 = {name:'lam_Aur'  , spectrum:'ukg0v.dat',   ra:TEN(14,50,15.81), dec:TEN(23,54,42.64), dist:18.2, mass:0., temp:5900., vsini:0.0, ldm:0.548, lde:0.010, uld: 0.00, pa:0.0, calfor:'none'} ; ldm from PTI (van Belle et al. 2008)
mu_gem  = {name:'mu_gem'  , spectrum:'ukm3iii.dat', ra:TEN(06,22,57.0), dec:TEN(22,30,48.9), dist:71.1, mass:0., temp:3650., vsini:0.0, ldm:14.64, lde:0.96, uld: 0.00, pa:0.0, calfor:'none'}    ; ldm to be checked
mu_uma  = {name:'mu_uma'  , spectrum:'ukm0iii.dat', ra:TEN(10,22,19.73), dec:TEN(41,29,58.26), dist:71.0, mass:0., temp:3700., vsini:0.0, ldm:8.54, lde:0.09, uld: 0.00, pa:0.0, calfor:'none'}   ; ldm from CHARM
procyon = {name:'procyon' , spectrum:'ukf5v.dat',   ra:TEN(0,0,0), dec:TEN(0,0,0), dist:0.0, mass:0.0, temp:6500. , vsini:0.0, ldm:5.5, lde:0.17, uld: 0.00, pa:0.0, calfor:'none'}
rigel   = {name:'rigel'   , spectrum:'ukb8i.dat',   ra:TEN(0,0,0), dec:TEN(0,0,0), dist:0.0, mass:0.0, temp:12130., vsini:0.0, ldm:2.55, lde:0.05, uld: 0.00, pa:0.0, calfor:'none'}
sirius  = {name:'sirius'  , spectrum:'uka0v.dat',   ra:TEN(06,45,08.92), dec:TEN(-16,42,58.01), dist:0.0, mass:0.0, temp:9940. , vsini:0.0, ldm:5.89, lde:0.16, uld: 0.00, pa:0.0, calfor:'none'}
ups_boo = {name:'ups_boo' , spectrum:'ukk5iii.dat', ra:TEN(13,47,28.54), dec:TEN(15,47,52.4), dist:81.0, mass:0., temp:3950.0, vsini:0.0, ldm:4.62, lde:0.05, uld: 0.00, pa:0.0, calfor:'none'}   ; ldm from CHARM

; HOSTS stars
alf_lyr = {name:'alf_Lyr', spectrum:'uka0v.dat', ra:TEN(18,36,56), dec:TEN(38,47,01), dist:7.76, mass:2.3, temp:10000., vsini:21.9, ldm:3.305, lde:0.01, uld:0.361, pa:8.6, calfor:'none'} ; best fit to K-band limb profile
bet_Eri = {name:'bet_Eri', spectrum:'uka3iii.dat', ra:TEN(05,07,51), dec:TEN(-05,05,11), dist:27.4, mass:0.0, temp:8137., vsini:0., ldm:1.16, lde:0.12, uld:0., pa:0., calfor:'none'} ; best fit to K-band limb profile
beta_leo = {name:'beta Leo', spectrum:'uka3v.dat', ra:TEN(11,49,03), dec:TEN(14,34,19), dist:11.1, mass:2.0, temp:8580., vsini:120., ldm:1.341, lde:0.013, uld:0.192, pa:118.0, calfor:'none'} ; Akeson09
bet_leo  = {name:'bet_Leo',  spectrum:'uka3v.dat', ra:TEN(11,49,03), dec:TEN(14,34,19), dist:11.1, mass:2.0, temp:8580., vsini:120., ldm:1.341, lde:0.013, uld:0.192, pa:118.0, calfor:'none'} ; Akeson09
beta_uma = {name:'beta_uma', spectrum:'uka2v.dat', ra:TEN(11,01,50), dec:TEN(56,22,56), dist:24.3, mass:2.28, temp:9600., vsini:46., ldm:1.095, lde:0.022, uld:0.175, pa:113.8, calfor:'none'} ; PA from Booth et al. 2013, No A1V spectrum in Pickles (aprroximated by an A2V here) 
bet_uma = {name:'bet_uma', spectrum:'uka2v.dat', ra:TEN(11,01,50), dec:TEN(56,22,56), dist:24.3, mass:2.28, temp:9600., vsini:46., ldm:1.095, lde:0.022, uld:0.175, pa:113.8, calfor:'none'} ; PA from Booth et al. 2013, No A1V spectrum in Pickles (aprroximated by an A2V here)
del_Crv = {name:'del_Crv', spectrum:'uka0iv.dat', ra:TEN(11,01,50), dec:TEN(56,22,56), dist:18.2, mass:0., temp:6875., vsini:0., ldm:0.82, lde:0.12, uld:0., pa:0., calfor:'none'} ;
del_leo = {name:'del_Leo', spectrum:'uka3v.dat', ra:TEN(11,14,07), dec:TEN(20,31,25), dist:17.9, mass:2.1, temp:8120., vsini:180., ldm:1.165, lde:0.022, uld:0.192, pa:0.0, calfor:'none'}   ; diameter from Akeson09
del_Uma = {name:'del_UMa', spectrum:'ukf0iii.dat', ra:TEN(12,32,04), dec:-TEN(16,11,46), dist:24.7, mass:1.44, temp:8725., vsini:92., ldm:0.81, lde:0.13, uld:0.00, pa:0.0, calfor:'none'} ; Dummy
eps_eri  = {name:'eps_eri', spectrum:'ukk2v.dat', ra:TEN(03,32,56), dec:-TEN(09,27,30), dist:3.22, mass:0.84, temp:5122., vsini:2.0, ldm:2.126, lde:0.014, uld:0.333, pa:0.0, calfor:'none'} ; DiFolco07
ksi_gem  = {name:'ksi_gem',  spectrum:'ukf5iv.dat', ra:TEN(06,45,17.4), dec:TEN(12,53,44), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.96, lde:0.032, uld: 0.00, pa:0.0, calfor:'none'}  ; ldm to be checked
tau_cet  = {name:'tau_Cet', spectrum:'ukg8v.dat', ra:TEN(01,44,04), dec:-TEN(15,56,15), dist:3.65, mass:0.83, temp:5377., vsini:0.5, ldm:2.015, lde:0.011, uld:0.333, pa:0.0, calfor:'none'} ; DiFolco07
gam_uma = {name:'gam_UMa', spectrum:'uka0v.dat', ra:TEN(11,53,50), dec:TEN(53,41,41), dist:25.6, mass:2.9, temp:9440., vsini:178., ldm:1.078, lde:0.055, uld:0.0, pa:0.0, calfor:'none'}
eta_crv = {name:'eta_crv', spectrum:'ukf2v.dat', ra:TEN(12,32,04), dec:-TEN(16,11,46), dist:18.2, mass:1.44, temp:7000., vsini:92., ldm:0.819, lde:0.119, uld:0.236, pa:116.3, calfor:'none'} ; PA from Duchene et al. 2014, diameter from Absil 2013
psc_107 = {name:'107_psc', spectrum:'ukk1iv.dat', ra:TEN(12,32,04), dec:-TEN(16,11,46), dist:7.5, mass:1.44, temp:5181., vsini:92., ldm:1.02, lde:0.13, uld:0.236, pa:116.3, calfor:'none'} ; Dummy
ksi_peg = {name:'ksi_peg', spectrum:'ukf6v.dat', ra:TEN(12,32,04), dec:-TEN(16,11,46), dist:16.3, mass:1.44, temp:6204., vsini:92., ldm:1.04, lde:0.14, uld:0.236, pa:116.3, calfor:'none'} ; Dummy
ksi_peg_A = {name:'ksi_peg_A', spectrum:'ukf6v.dat', ra:TEN(12,32,04), dec:-TEN(16,11,46), dist:16.3, mass:1.44, temp:6204., vsini:92., ldm:1.04, lde:0.14, uld:0.236, pa:116.3, calfor:'none'} ; Dummy
gj_105A = {name:'gj_105A', spectrum:'ukk3v.dat', ra:TEN(12,32,04), dec:-TEN(16,11,46), dist:7.18, mass:1.44, temp:4866., vsini:92., ldm:1.000, lde:0.10, uld:0.236, pa:116.3, calfor:'none'} ; Dummy
mu_vir  = {name:'mu_vir'  , spectrum:'ukf2v.dat', ra:TEN(10,22,19.73), dec:TEN(41,29,58.26), dist:18.3, mass:0., temp:7079., vsini:0.0, ldm:0.94, lde:0.12, uld: 0.00, pa:0.0, calfor:'none'}   ; 
sig_boo = {name:'sig_Boo', spectrum:'ukf2v.dat', ra:TEN(14,34,41), dec:TEN(29,44,42), dist:15.5, mass:1.28, temp:6400., vsini:15., ldm:0.783, lde:0.016, uld:0.259, pa:0.0, calfor:'none'}
iot_psc = {name:'iot_psc', spectrum:'ukf6v.dat', ra:TEN(12,32,04), dec:-TEN(16,11,46), dist:13.7, mass:1.44, temp:6397., vsini:92., ldm:1.01, lde:0.14, uld:0.236, pa:116.3, calfor:'none'} ; Dummy
tet_Boo = {name:'tet_Boo ', spectrum:'ukf6v.dat', ra:TEN(12,32,04), dec:-TEN(16,11,46), dist:14.5, mass:1.44, temp:6247., vsini:92., ldm:1.17, lde:0.19, uld:0.236, pa:116.3, calfor:'none'} ; Dummy
uma_23  = {name:'23_uma', spectrum:'ukf0iii.dat', ra:TEN(12,32,04), dec:-TEN(16,11,46), dist:23.8, mass:1.44, temp:7186., vsini:92., ldm:1.01, lde:0.14, uld:0.236, pa:116.3, calfor:'none'} ; Dummy
ori_1   = {name:'1_Ori', spectrum:'ukf0iii.dat', ra:TEN(12,32,04), dec:-TEN(16,11,46), dist:8.07, mass:1.44, temp:6392., vsini:92., ldm:1.56, lde:0.05, uld:0.00, pa:0.0, calfor:'none'} ; Dummy
Leo_40  = {name:'40_Leo', spectrum:'ukf0iii.dat', ra:TEN(12,32,04), dec:-TEN(16,11,46), dist:21.36, mass:1.44, temp:7193., vsini:92., ldm:0.59, lde:0.09, uld:0.00, pa:0.0, calfor:'none'} ; Dummy

; Calibrators
hd99196  = {name:'HD99196', spectrum:'ukk4iii.dat', ra:TEN(11,24,58.9), dec:TEN(11,25,49.0), dist:150.4, mass:0., temp:3882., vsini:0.0, ldm:1.55, lde:0.06, uld: 0.00, pa:0.0, calfor:'bet_Leo'}  ; from Colavita et al. 2009
hd104979 = {name:'HD104979', spectrum:'ukg8iii.dat', ra:TEN(12,05,12.5), dec:TEN(08,43,58.7), dist:50.05, mass:0., temp:3882., vsini:0.0, ldm:1.64, lde:0.40, uld: 0.00, pa:0.0, calfor:'bet_Leo'} ; no info on this star
hd97603  = {name:'HD97603', spectrum:'uka5v.dat', ra:TEN(11,14,07), dec:TEN(20,31,25), dist:17.7, mass:2.1, temp:8450., vsini:180., ldm:1.165, lde:0.022, uld:0.192, pa:0.0, calfor:'bet_Leo'}   ; This is det_leo -- Akeson09
hd14872  = {name:'HD14872', spectrum:'ukk4iii.dat', ra:TEN(02,25,37.42351), dec:TEN(50,16,43.0645), dist:113.25, mass:0., temp:4097., vsini:0., ldm:3.28, lde:0.056, uld:0.0, pa:0.0, calfor:'gamPer'} ; ldm and lde from Borde  
hd108522  = {name:'HD108522', spectrum:'ukk4iii.dat', ra:TEN(12,28,01.8), dec:TEN(-14,27,03.4), dist:313.48, mass:0., temp:4200., vsini:0.0, ldm:1.204, lde:0.016, uld: 0.00, pa:0.0, calfor:'none'}  ; Teff and distance from McDonald 2012
hd107418  = {name:'HD107418', spectrum:'ukk0iii.dat', ra:TEN(12,20,55.71287), dec:TEN(-13,33,56.6100), dist:61.69, mass:0., temp:4948., vsini:0.0, ldm:1.335, lde:0.092, uld: 0.00, pa:0.0, calfor:'none'}  ;  Teff and distance from McDonald 2012
hd109272  = {name:'HD109272', spectrum:'ukg8iii.dat', ra:TEN(12,33,34.25847), dec:TEN(-12,49,48.7268), dist:49.36, mass:0., temp:5300., vsini:0.0, ldm:0.907, lde:0.063, uld: 0.00, pa:0.0, calfor:'none'}  ;  Teff and distance from McDonald 2012
hd92095   = {name:'HD92095', spectrum:'ukk3iii.dat', ra:TEN(11,24,58.9), dec:TEN(11,25,49.0), dist:150.4, mass:0., temp:3882., vsini:0.0, ldm:1.57, lde:0.02, uld: 0.00, pa:0.0, calfor:'none'}  ; Bogus (ldm and lde from Mérand)
hd102328  = {name:'HD102328', spectrum:'ukk3iii.dat', ra:TEN(11,24,58.9), dec:TEN(11,25,49.0), dist:150.4, mass:0., temp:3882., vsini:0.0, ldm:1.483, lde:0.11, uld: 0.00, pa:0.0, calfor:'none'}  ; Bogus (but ldm and lde from Colavita et al. 2009)
hd163770  = {name:'HD163770', spectrum:'ukk1iii.dat', ra:TEN(17,56,16), dec:TEN(37,15,02), dist:150.4, mass:0., temp:4467., vsini:0., ldm:3.15, lde:0.034, uld:0.403, pa: 0.0, calfor:'none'}  ; Bogus (but ldm and lde from Borde)
hd168775  = {name:'HD168775', spectrum:'ukk2iii.dat', ra:TEN(18,19,52), dec:TEN(36,03,52), dist:150.4, mass:0., temp:4380., vsini:0., ldm:2.28, lde:0.025, uld:0.410, pa: 0.0, calfor:'none'}  ; Bogus (but ldm and lde from Borde)
hd69830   = {name:'HD69830', spectrum:'ukg8v.dat', ra:TEN(08,18,23), dec:-TEN(12,37,55), dist:12.6, mass:0.8, temp:5385., vsini:0., ldm:0.63, lde:0.1, uld:0.3, pa:0.0, calfor:'none'} ; Mass and Teff from AllendePrieto99 -- uld TBC
hd65098   = {name:'HD65098', spectrum:'ukg8v.dat', ra:TEN(08,18,23), dec:-TEN(12,37,55), dist:12.6, mass:0.8, temp:5385., vsini:0., ldm:0.63, lde:0.1, uld:0.3, pa:0.0, calfor:'none'} ; Bogus
hd70409   = {name:'HD70409', spectrum:'ukg8v.dat', ra:TEN(08,18,23), dec:-TEN(12,37,55), dist:12.6, mass:0.8, temp:5385., vsini:0., ldm:0.63, lde:0.1, uld:0.3, pa:0.0, calfor:'none'} ; Bogus
NAC       = {name:'NAC', spectrum:'ukg8v.dat', ra:TEN(08,18,23), dec:-TEN(12,37,55), dist:0., mass:0., temp:5500., vsini:0., ldm:0.1, lde:0.0001, uld:0., pa:0.0, calfor:'none'}
hd18322   = {name:'HD18322', spectrum:'ukk1iii.dat', ra:TEN(02,56,25), dec:TEN(08,53,53), dist:0.0, mass:0., temp:4500., vsini:0.0, ldm:2.580, lde:0.06, uld: 0.00, pa:0.0, calfor:'none'} ; ldm and dle from Colavita et al 2009
hd23249   = {name:'HD23249', spectrum:'ukk1iii.dat', ra:TEN(03,43,15), dec:TEN(-09,45,48), dist:0.0, mass:0., temp:4500., vsini:0.0, ldm:1.810, lde:0.025, uld: 0.00, pa:0.0, calfor:'none'} 
hd29065   = {name:'HD29065', spectrum:'ukk1iii.dat', ra:TEN(03,43,15), dec:TEN(-09,45,48), dist:0.0, mass:0., temp:4500., vsini:0.0, ldm:2.318, lde:0.05, uld: 0.00, pa:0.0, calfor:'none'}  ; Garbage (but ldm and lde from Colavita et al. 2009)
hd92424   = {name:'HD92424', spectrum:'ukk1iii.dat', ra:TEN(03,43,15), dec:TEN(-09,45,48), dist:0.0, mass:0., temp:4500., vsini:0.0, ldm:1.635, lde:0.02, uld: 0.00, pa:0.0, calfor:'none'}  ; Garbage (but ldm/lde from Mérand)
hd112769  = {name:'HD112769', spectrum:'ukk5iii.dat', ra:TEN(12,58,55.4), dec:TEN(17,24,34.0), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:4.6, lde:0.05, uld: 0.00, pa:0.0, calfor:'eps_vir'}  ; ldm dirty from SBR
hd109742  = {name:'HD109742', spectrum:'ukk4iii.dat', ra:TEN(03,43,15), dec:TEN(-09,45,48), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.8, lde:0.021, uld: 0.00, pa:0.0, calfor:'bet_Leo'}
hd108381  = {name:'HD108381', spectrum:'ukk1iii.dat', ra:TEN(03,43,15), dec:TEN(-09,45,48), dist:0.0, mass:0., temp:4508., vsini:0.0, ldm:2.15, lde:0.024, uld: 0.00, pa:0.0, calfor:'bet_Leo'} ; Borde
hd48433   = {name:'HD48433',  spectrum:'ukk1iii.dat', ra:TEN(06,43,59.3), dec:TEN(13,13,40.9), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.07, lde:0.027, uld: 0.00, pa:0.0, calfor:'ksi_gem'} ; Borde (from Bordé)
hd52976   = {name:'HD52976',  spectrum:'ukk5iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.96, lde:0.032, uld: 0.00, pa:0.0, calfor:'ksi_gem'} ; ldm and sp from Borde
HD198149  = {name:'HD198149',  spectrum:'ukk5iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.68, lde:0.029, uld: 0.00, pa:0.0, calfor:'alf_cep'} ; ldm/e from Mérand
HD209960  = {name:'HD209960',  spectrum:'ukk5iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'alf_cep'} ; ldm dirty from SBR 
HD49968   = {name:'HD49968',  spectrum:'ukk5iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'107_psc'} ; dummy
HD218792  = {name:'HD218792',  spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'ksi_peg'} ; dummy
HD7087    = {name:'HD7087',  spectrum:'ukg8iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'ksi_peg'} ; dummy
HD209167  = {name:'HD209167',  spectrum:'ukk5iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'ksi_peg'} ; dummy
HD220009  = {name:'HD220009',  spectrum:'ukk1iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'iot_Psc'} ; dummy
HD21051   = {name:'HD21051',  spectrum:'ukk0iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'GJ_105A'} ; dummy
HD13596   = {name:'HD13596',  spectrum:'ukm0iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'GJ_105A'} ; dummy
HD52960   = {name:'HD52960',  spectrum:'ukk3iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'ksi_gem'} ; dummy
HD86378   = {name:'HD86378',  spectrum:'ukk5iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'bet_Uma'} ; dummy !Also for 23_Uma!!
HD7318    = {name:'HD7318',  spectrum:'ukk5iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'107_psc'} ; dummy
HD6953    = {name:'HD6953',  spectrum:'ukk5iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:2.44, lde:0.032, uld: 0.00, pa:0.0, calfor:'107_psc'} ; dummy
HD38656   = {name:'HD38656',  spectrum:'ukg8iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'lam_Aur'} ; dummy
HD40441   = {name:'HD40441',  spectrum:'ukk3iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'lam_Aur'} ; dummy
HD31421   = {name:'HD31421',  spectrum:'ukk2iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'1_Ori'}   ; dummy
HD31767   = {name:'HD31767',  spectrum:'ukk0iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'1_Ori'}   ; or bet_Eri!!! dummy
HD89024   = {name:'HD89024',  spectrum:'ukk2iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'40_Leo'}  ; dummy
HD93257   = {name:'HD93257',  spectrum:'ukk2iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'40_Leo'}  ; dummy
HD107465  = {name:'HD107465', spectrum:'ukk2iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'del_uma'}  ; dummy
HD102328  = {name:'HD102328', spectrum:'ukk2iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'del_uma'}  ; d
HD128902  = {name:'HD128902', spectrum:'ukk2iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'tet_Boo'}  ; 
HD128000  = {name:'HD128000', spectrum:'ukk2iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'tet_Boo'}  ;
HD36780   = {name:'HD36780', spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'bet_Eri'}  ; dummy 
HD31767   = {name:'HD31767',  spectrum:'ukk0iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'bet_Eri'}   ; or 1_Ori!!! dummy (comment out the one yo d'ont need)
HD94336   = {name:'HD94336', spectrum:'ukm1iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'del_leo'}  ; dummy
HD99902   = {name:'HD99902', spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'del_leo'}  ; dummy
HD114113  = {name:'HD114113', spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'del_Crv'}  ; dummy
HD111500  = {name:'HD111500', spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'del_Crv'}  ; dummy
HD131477  = {name:'HD131477', spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'mu_Vir'} 
HD133165  = {name:'HD133165', spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'mu_Vir'} 
HD130952  = {name:'HD130952', spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'mu_Vir'} 
HD95212   = {name:'HD95212',  spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'bet_uma'}
HD94247   = {name:'HD94247',  spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'bet_uma'}
HD164646  = {name:'HD164646',  spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'sig_boo'}
HD133392  = {name:'HD133392',  spectrum:'ukk4iii.dat', ra:TEN(07,03,51.6), dec:TEN(12,35,39.3), dist:0.0, mass:0., temp:3650., vsini:0.0, ldm:1.74, lde:0.12, uld: 0.00, pa:0.0, calfor:'sig_boo'}

; IOTA observations 
; *****************

pi_her   = {name:'pi_Her', spectrum:'ukk3ii.dat', ra:TEN(17,15,03), dec:TEN(36,48,33), temp:4246, vsini:0, ldm:5.29, lde:0.055, uld:0.420, calfor:'none'}   ; Borde
lam_lyr  = {name:'lam_Lyr', spectrum:'ukm2iii.dat', ra:TEN(19,00,01), dec:TEN(32,08,44), temp:3690, vsini:0, ldm:2.41, lde:0.026, uld:0.415, calfor:'none'}  ; Borde
kap_lyr  = {name:'kap_Lyr', spectrum:'ukk2iii.dat', ra:TEN(18,19,52), dec:TEN(36,03,52), temp:4380, vsini:0, ldm:2.28, lde:0.025, uld:0.410, calfor:'none'}   ; Borde
thet_her = {name:'thet_Her', spectrum:'ukk1ii.dat', ra:TEN(17,56,16), dec:TEN(37,15,02), temp:4467, vsini:0, ldm:3.15, lde:0.034, uld:0.403, calfor:'none'} ; Borde
; Data sylvestre
HD120477 = {name:'HD120477', ra:TEN(13,49,28), dec:TEN(15,47,52), vsini:0.0, temp:4002, ldm:4.62, lde:0.050, uld:0.436, calfor:'none'} ; Borde&Merand
HD129972 = {name:'HD129972', ra:TEN(14,45,14), dec:TEN(16,57,51), vsini:0.0, temp:5535, ldm:1.571, lde:0.020, uld:0.4, calfor:'none'}   ; Borde&Merand for ldm&lde, approximation for temp&uld
HD125560 = {name:'HD125560', ra:TEN(14,19,45), dec:TEN(16,18,25), vsini:0.0, temp:4256, ldm:1.96, lde:0.021, uld:0.420, calfor:'none'}   ; Borde&Merand


; Define array of targets
; ***********************

IF NOT KEYWORD_SET(CALIB) THEN BEGIN
  tgt_star = [gam_uma, eta_crv, sig_boo, alf_crb, gam_oph, zet_aql, bet_leo, del_leo, alf_lyr, tau_cet, eps_eri, alf_psa, $
              bet_pic, bet_picK, au_mic,$                                                                                                       ; AM Lagrange data
              hd142, hd7570  , hd7788  , hd15008 , hd90132 , hd91324 , hd99211 , hd102365, hd104731, hd108767, hd109787, hd115617, hd120136, $
              hd128898, hd129502, hd130109, hd134083, hd135379, hd136202, hd139664, hd141891, hd149661, hd152391, hd160032, hd160915, hd164259, $
              hd165777, hd172555, hd178253, hd182572, hd188228, hd192425, hd195627, hd197157, hd197692, hd202730, hd203608, hd206860, hd207129, $
              hd210049, hd210302, hd210418, hd213845, hd215789, hd219482, hd219571, hd224392, $                                                 ; PIONIER survey 
              alp_aql, alp_boo, alf_cep, alp_ser, alp_tau, alp_ori, bet_and, bet_and2, bet_boo, bet_Eri, bet_gem, bet_peg, bet_uma, del_Crv, del_uma, eps_vir, eps_leo, $
              eta_boo, gj_105A, iot_psc, ksi_boo, ksi_gem, ksi_peg, ksi_peg_A, lam_aur, lam_aur2, leo_40, mu_uma, mu_gem, mu_vir, NAC, ori_1, psc_107, procyon, $
              rigel, sirius, tet_Boo, uma_23, ups_boo, $
              hd99196, hd104979, hd97603,hd14872,gam_per,beta_leo,hd108522,hd107418,hd109272,hd92095,beta_uma,hd102328,hd163770,hd168775,hd69830,$
              hd65098,hd70409,hd18322,hd23249,hd29065,hd92424,hd112769,hd109742,hd108381,hd48433,hd52976,HD198149,HD209960,HD49968,HD218792,HD7087,$
              HD220009,HD21051,HD13596,HD52960,HD86378,HD209167,HD7318,HD6953,HD38656,HD40441,HD31421,HD31767, HD89024, HD93257, HD107465, HD102328, $
              HD128902, HD128000, HD36780, HD94336, HD99902, HD114113, HD111500, HD131477, HD133165, HD130952,HD95212,HD94247,HD133392,HD164646]                               ; LBTI survey                                                                                                 
ENDIF ELSE BEGIN
  tgt_star = [pi_her, lam_lyr, kap_lyr, thet_her]
ENDELSE      

; Find the input target
; *********************

idx = WHERE(STRCOMPRESS(STRLOWCASE(tgt_star.name), /REMOVE_ALL) EQ STRCOMPRESS(STRLOWCASE(input_name), /REMOVE_ALL), n_tgt)
IF n_tgt LE 0 THEN PRINT, 'Input target does not exist in the data base : ' + input_name
IF n_tgt GT 1 THEN PRINT, 'Input target exists with various definitions. Pick one above.'
IF n_tgt EQ 1 THEN tgt_out = tgt_star[idx] ELSE tgt_out = {name:'not found', calfor:'none'}
END