function enhancepx, epsilon
ddim = size(epsilon,/dimensions)
eps = fltarr([ddim,9])
eps[*,*,0] = epsilon
eps[*,*,1] = shift(epsilon,0,1)
eps[*,*,2] = shift(epsilon,0,-1)
eps[*,*,3] = shift(epsilon,1,1)
eps[*,*,4] = shift(epsilon,1,0)
eps[*,*,5] = shift(epsilon,1,-1)
eps[*,*,6] = shift(epsilon,-1,1)
eps[*,*,7] = shift(epsilon,-1,0)
eps[*,*,8] = shift(epsilon,-1,-1)

epsA = fltarr(ddim)
epsN = intarr(ddim)
for j=0,8 do begin
  epsA += eps[*,*,j]
  epsN += eps[*,*,j] gt 0
endfor
epsI = where(epsN gt 0)
epsA[epsI] = epsA[epsI]/epsN[epsI]
return, epsA 
end
