function Halphaline
; calculation taken from the Paper "An entropy based classification scheme for land applications of Polarimetric SAR", by Shane Cloude, 1997
  N = 100
  HmI       = fltarr(N)
  HmII1     = fltarr(N)
  HmII2     = fltarr(N)
  alphamI   = fltarr(N)
  alphamII1 = fltarr(N)
  alphamII2 = fltarr(N)
  
  for j=1,N-1 do begin
    m = 1.0*j/(N-1)
    m2 = (m/2+0.5)
    mlambdaI   = [1,m,m]
    mlambdaII1 = [0,1,2*m]
    mlambdaII2 = [2*m2-1,1,1]
    HmI[j]       = -total(mlambdaI/(2*m+1)*alog(mlambdaI/(2*m+1))/alog(3))
;    HmII1[j]     = -total(mlambdaII1/(2*m+1)*alog(mlambdaII1/(2*m+1))/alog(3))
    HmII2[j]     = -total(mlambdaII2/(2*m2+1)*alog(mlambdaII2/(2*m2+1))/alog(3))
    alphamI[j]   = 180*m/(2*m+1)
;    alphamII1[j] = (2*m+1)/(2*m+1)*90
    alphamII2[j] = 180/(2*m2+1)    
  endfor
  alphamI[0]   = 0
  alphamII2[0] = 90
  HmI[0]       = 0
  HmII2[0]     = -total([1,1]/(2.0)*alog([1,1]/(2.0))/alog(3)) 
  return, [[HmI,reverse(HmII2)],[alphamI,reverse(alphamII2)]]
end