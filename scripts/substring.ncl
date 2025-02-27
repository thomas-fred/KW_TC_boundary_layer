;
; substring.ncl -- Extract a substring from another string.
;
; 2007-mar-05	Original version.  By Dave Allured.
;
; Usage:  outstr = substring (instr, first, last)
;
; Input:  instr = input string
;	  first = starting index of substring (zero based)
;	  last = ending index of substring (zero based); if last is
;	     less than first, the entire right substring is returned
;
; Output: outstr = selected substring
;
; General notes:
;
; Range checking is not done in this version.  Specifying "first" or
; "last" out of range will result in a non-elegant subscript error.
;
; Organization:
; University of Colorado, CIRES Climate Diagnostics Center (CDC)
; NOAA/ESRL/PSD, Climate Analysis Branch (CAB)
;
;----------------------------------------------------------------------------

function substring (instr[1]:string, first[1]:numeric, last[1]:numeric)

local instr, first, last, main, p2

begin
   main = stringtochar (instr)		; convert input string to char array
					; terminal null character is included   
   
   if (last .ge. first) then		; check requested end position
      p2 = last				; go to last position specified
   else					; but if less than first:
      p2 = dimsizes (main) - 2		; go to last avail char in main string
   end if
   
   return (chartostring (main(first:p2)))	; extract substring
end
