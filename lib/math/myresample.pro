;+
;  NAME: MYRESAMPLE
;  
;  DESCRITPION
;  My own version to resample and replacement a vector with N points. So you start with a vector with N indices, and from that vector you randomly build a new
;  vector of indices where all the the indices are one by one randomly drawn from the original ones (you can randomly draw the same indice several times).
;
;  MODIFICATION HISTORY
;   Version 1.0, 27-OCT-2012, by Olivier Absil, ULg, absil@astro.ulg.ac.be
;   Version 1.1, 27-NOV-2015, DD: cleaned up code

FUNCTION MYRESAMPLE, n, seed

;ind=uintarr(N) ; an array of N integers
ind=ulon64arr(N)
;for i=ulong(0),N-1 do begin
;ind(i)=fix(0.9999*randomu(seed+i,1)*(N))
;endfor

ind=ulong(randomu(seed,N)*float(N-1))

RETURN, ind

;print,ind
end