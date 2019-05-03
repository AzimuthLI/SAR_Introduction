function Haline
; calculation taken from the Paper "An entropy based classification scheme for land applications of Polarimetric SAR", by Shane Cloude, 1997
  N = 1000
  HmI       = fltarr(N)
  HmII1     = fltarr(N)
  HmII2     = fltarr(N)
  AI   = fltarr(N)
  AII1 = fltarr(N)
  AII2 = fltarr(N)
  j = 0
  for j=1,N-1 do begin
    m = 1.0*j/(N-1)
    m1 = m/2 
    m2 = m/2+0.5
    mlambdaI   = [1,m,m]
    mlambdaII1 = [0,1,2*m1] ; 
    mlambdaII2 = [2*m2-1,1,1]
;    HmI[j]       = -total(mlambdaI/(2*m+1)*alog(mlambdaI/(2*m+1))/alog(3))
;    HmII1[j]     = -total(mlambdaII1[1:2]/(2*m1+1)*alog(mlambdaII1[1:2]/(2*m1+1))/alog(3))
    HmII2[j]     = -total(mlambdaII2/(2*m2+1)*alog(mlambdaII2/(2*m2+1))/alog(3))
;    AI[j]   = 0
;    AII1[j]   = (2*m1)/(2*m1)
    AII2[j] = (2-2*m2)/(2*m2)    
  endfor
  AII2[0] = 1
  HmII2[0] = -total([1,1]/(2.0)*alog([1,1]/2.0)/alog(3))  
  return, [[HmII2],[AII2]]
  end