;+
; NAME: GET_PATH
; 
; PURPOSE:
;   Returns the path to specified procedure.
;   Used to find local path, independent of overall directory layout

; KEYWORDS
;   N_DIR_UP: go up by this number of directories
;   FUNC:     must be set to use GET_PATH with a function instead of a procedure
;
; MODIFICATION HISTORY:
;   Version 1.0,  09-DEC-2003, R. den Hartog, ESA/ESTEC//SCI-A/Genie Team, rdhartog@rssd.esa.int
;   Version 1.1,  12-MAR-2004, OA: implemented keyword FUNC
;   Version 1.2,  02-DEC-2015, DD: cleanup comments

FUNCTION GET_PATH, routine, N_DIR_UP=n_dir_up, FUNC=func

ON_ERROR,2

; Convert to uppercase, because that's how it is stored in the structure
name=STRLOWCASE(routine)

; Remove the '.pro' when asking for routine info
p=STRPOS(name,'.pro')
IF p GT 0 THEN name=STRMID(name, 0, p)
r = ROUTINE_INFO(name, /SOURCE, FUNC=func)
path=r.path
IF NOT KEYWORD_SET(n_dir_up) THEN n_dir_up=0

; Get the separator between directories right
IF !VERSION.OS_FAMILY EQ 'unix' THEN sep='/' ELSE sep='\'

; Go up number of specified directories
FOR i=0, n_dir_up DO BEGIN
  p=STRPOS(path, sep, /REVERSE_SEARCH)
  path=STRMID(path, 0, p)
ENDFOR

; Put one separator back to make it a useful path
RETURN, path+sep
END

;--------------------------------------------

PRO TEST_PATH
PRINT, GET_PATH('get_path.pro', N_DIR_UP=3)
 print, routine_info('generate_cloud', /source)
END