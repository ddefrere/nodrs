;+
; NAME: DEFINE_KEYMAP
;
; PURPOSE:
;   Simple procedure to define the keyword number mapping for faster reading of the L0 header
; 
; INPUT/OUTPUT
;   key_map : Input (resp, output) data structure to save (resp. to read)
;   p(n)    : Position number of the nth keyword (to save or to read)
;      
; KEYWORD
;   SAVE    : If set, key_map is the output and the p(n) the inputs. It is the other way around if not set.
;
; MODIFICATION HISTORY:
;   Version 1.0, 02-DEC-2015, by Denis Defrère, University of Arizona, ddefrere@email.arizona.edu (adapted from former routine 'nomic_null.pro')
;   Version 1.1, 16-DEC-2015, DD: added 5 more inputs
;   Version 1.2, 06-JUL-2016, DD: added 1 more parameter
;   Version 1.3, 19-OCT-2016, DD: added 2 more parameters
;   Version 1.4, 09-NOV-2016, DD: added 6 more parameters
;   Version 1.5, 30-JAN-2017, DD: added p90 to p99

PRO DEFINE_KEYMAP, key_map, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, $
                            p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54,  $
                            p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, $
                            p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96, p97, p98, p99, SAVE=save
                             
IF KEYWORD_SET(SAVE) THEN BEGIN
  key_map = {p1:p1, p2:p2, p3:p3, p4:p4, p5:p5, p6:p6, p7:p7, p8:p8, p9:p9, $
             p10:p10, p11:p11, p12:p12, p13:p13, p14:p14, p15:p15, p16:p16, p17:p17, p18:p18, p19:p19, $
             p20:p20, p21:p21, p22:p22, p23:p23, p24:p24, p25:p25, p26:p26, p27:p27, p28:p28, p29:p29, $
             p30:p30, p31:p31, p32:p32, p33:p33, p34:p34, p35:p35, p36:p36, p37:p37, p38:p38, p39:p39, $
             p40:p40, p41:p41, p42:p42, p43:p43, p44:p44, p45:p45, p46:p46, p47:p47, p48:p48, p49:p49, $
             p50:p50, p51:p51, p52:p52, p53:p53, p54:p54, p55:p55, p56:p56, p57:p57, p58:p58, p59:p59, $
             p60:p60, p61:p61, p62:p62, p63:p63, p64:p64, p65:p65, p66:p66, p67:p67, p68:p68, p69:p69, $
             p70:p70, p71:p71, p72:p72, p73:p73, p74:p74, p75:p75, p76:p76, p77:p77, p78:p78, p79:p79, $
             p80:p80, p81:p81, p82:p82, p83:p83, p84:p84, p85:p85, p86:p86, p87:p87, p88:p88, p89:p89, $
             p90:p90, p91:p91, p92:p92, p93:p93, p94:p94, p95:p95, p96:p96, p97:p97, p98:p98, p99:p99}
ENDIF ELSE BEGIN
  p1  = key_map.p1  & p2 = key_map.p2   & p3  = key_map.p3  & p4  = key_map.p4  & p6  = key_map.p6  & p6  = key_map.p6  & p7  = key_map.p7  & p8  = key_map.p8  & p9  = key_map.p9
  p10 = key_map.p10 & p11 = key_map.p11 & p12 = key_map.p12 & p13 = key_map.p13 & p14 = key_map.p14 & p15 = key_map.p15 & p16 = key_map.p16 & p17 = key_map.p17 & p18 = key_map.p18 & p19 = key_map.p19
  p20 = key_map.p20 & p21 = key_map.p21 & p22 = key_map.p22 & p23 = key_map.p23 & p24 = key_map.p24 & p25 = key_map.p25 & p26 = key_map.p26 & p27 = key_map.p27 & p28 = key_map.p28 & p29 = key_map.p29
  p30 = key_map.p30 & p31 = key_map.p31 & p32 = key_map.p32 & p33 = key_map.p33 & p34 = key_map.p34 & p35 = key_map.p35 & p36 = key_map.p36 & p37 = key_map.p37 & p38 = key_map.p38 & p39 = key_map.p39
  p40 = key_map.p40 & p41 = key_map.p41 & p42 = key_map.p42 & p43 = key_map.p43 & p44 = key_map.p44 & p45 = key_map.p45 & p46 = key_map.p46 & p47 = key_map.p47 & p48 = key_map.p48 & p49 = key_map.p49
  p50 = key_map.p50 & p51 = key_map.p51 & p52 = key_map.p52 & p53 = key_map.p53 & p54 = key_map.p54 & p55 = key_map.p55 & p56 = key_map.p56 & p57 = key_map.p57 & p58 = key_map.p58 & p59 = key_map.p59
  p60 = key_map.p60 & p61 = key_map.p61 & p62 = key_map.p62 & p63 = key_map.p63 & p64 = key_map.p64 & p65 = key_map.p65 & p66 = key_map.p66 & p67 = key_map.p67 & p68 = key_map.p68 & p69 = key_map.p69
  p70 = key_map.p70 & p71 = key_map.p71 & p72 = key_map.p72 & p73 = key_map.p73 & p74 = key_map.p74 & p75 = key_map.p75 & p76 = key_map.p76 & p77 = key_map.p77 & p78 = key_map.p78 & p79 = key_map.p79
  p80 = key_map.p80 & p81 = key_map.p81 & p82 = key_map.p82 & p83 = key_map.p83 & p84 = key_map.p84 & p85 = key_map.p85 & p86 = key_map.p86 & p87 = key_map.p87 & p88 = key_map.p88 & p89 = key_map.p89
  p90 = key_map.p90 & p91 = key_map.p91 & p92 = key_map.p92 & p93 = key_map.p93 & p94 = key_map.p94 & p95 = key_map.p95 & p96 = key_map.p96 & p97 = key_map.p97 & p98 = key_map.p98 & p99 = key_map.p99
ENDELSE
                             
                                                                                                               
END
